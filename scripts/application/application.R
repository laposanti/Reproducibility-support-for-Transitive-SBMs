# ======================================================================
# Applications driver: LOOIC + unified order-conformity diagnostics
# Mirrors simulation study (same functions, names, logic)
# ======================================================================
cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]),
                               winslash = "/", mustWork = TRUE)
  setwd(normalizePath(file.path(dirname(script_path), "../.."),
                      winslash = "/", mustWork = TRUE))
}

suppressPackageStartupMessages({
  library(Matrix)
  library(mcclust)
  library(mcclust.ext)   # minVI
  library(coda)
  library(dplyr)
  library(tidyr)
  library(loo)
  library(here)
  library(fs)
  library(readr)
  library(salso)
  library(glue)
  library(label.switching)  
  library(future.apply)
  library(purrr)
})

# ----------------------------------------------------------------------
# Empty diagnostics placeholder
# ----------------------------------------------------------------------

.empty_diag_vec <- function() {
  base_stats <- c("mean", "lo", "hi", "ess", "mcse")
  diag_prefixes <- c(
    "thetaW_block_model", "thetaS_block_model",
    "thetaW_block_emp", "thetaS_block_emp",
    "thetaW_item_emp", "thetaS_item_emp",
    "thetaW_block_model_all", "thetaS_block_model_all",
    "thetaW_block_emp_all", "thetaS_block_emp_all",
    "thetaW_item_emp_all", "thetaS_item_emp_all"
  )
  diag_names <- unlist(lapply(diag_prefixes, function(prefix) {
    paste0(prefix, "_", base_stats)
  }))
  tail_names <- c(
    "prem_block_model_avg", "prem_block_emp_avg", "prem_item_emp_avg",
    "coverage_block_model_avg", "coverage_block_emp_avg", "coverage_item_emp_avg",
    "transitive_triads_mean", "transitive_triads_lo", "transitive_triads_hi",
    "cycle_mass_weighted_mean",
    "min_backward_weight_mean", "min_backward_weight_norm_mean",
    "hierarchy_energy_mean", "hierarchy_energy_lo", "hierarchy_energy_hi",
    "curl_energy_mean",
    "p_post_wst", "p_prior_wst", "bf_wst_0", "log_bf_wst_0",
    "p_post_sst", "p_prior_sst", "bf_sst_0", "log_bf_sst_0",
    "bf_sst_wst", "log_bf_sst_wst", "bf_prior_mc",
    "violation_rate_mean", "violation_rate_lo", "violation_rate_hi",
    "violation_count_mean", "cross_mass_mean",
    "violation_rate_zhat", "violation_count_zhat"
  )
  all_names <- c(diag_names, tail_names)
  out <- rep(NA_real_, length(all_names))
  names(out) <- all_names
  out
}

.normalize_app_ordered_prior <- function(partition_prior) {
  prior <- toupper(trimws(partition_prior))
  if (!nzchar(prior)) prior <- "OCRP"
  if (prior %in% c("R-OCRP", "REVERSED_OCRP", "REVERSED OCRP",
                   "MIRRORED_OCRP", "MIRRORED OCRP")) {
    prior <- "OCRP"
  }
  if (!prior %in% c("OCRP", "ROCRP")) {
    stop("APP_OSBM_PARTITION_PRIOR must be one of OCRP or ROCRP.", call. = FALSE)
  }
  prior
}

# ---------------- Parallel plan ----------------------------------------

RNGkind("L'Ecuyer-CMRG")    # parallel-safe RNG
set.seed(42)

# ======================================================================
# Data ingest (unchanged)
# ======================================================================

choose_dataset <- function(dataset = c("mountain_goats","citations_data","macaques_data","high_school")) {
  dataset <- match.arg(dataset)
  if (dataset == "mountain_goats") {
    matrix_files <- list.files("./data/ShizukaMcDonald_Data",
                               full.names = TRUE, pattern = "\\.csv$")
    n_each <- vapply(matrix_files, function(f) nrow(read.csv(f, row.names = 1)),
                     FUN.VALUE = integer(1))
    A <- as.matrix(read.csv(matrix_files[which.max(n_each)],
                            row.names = 1, check.names = FALSE))
  } else if (dataset == "citations_data") {
    A <- read.csv("./data/Citations_application/cross-citation-matrix.csv",
                  row.names = 1, header = TRUE, check.names = FALSE)
    diag(A) <- 0
  } else if (dataset == "macaques_data") {
    edge_list <- read.table("./data/macaques/out.moreno.txt")
    nodes <- sort(unique(c(edge_list[[1]], edge_list[[2]])))
    A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (i in seq_len(nrow(edge_list))) {
      A[edge_list[i,1], edge_list[i,2]] <- edge_list[i,"V3"]
    }
  } else if (dataset == "high_school") {
    edges <- read.csv("./data/high-school/edges.csv",
                      header = FALSE,
                      comment.char = "#",
                      strip.white = TRUE)
    if (ncol(edges) < 3L) stop("high_school edges.csv must have source,target,weight columns")
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)

    if (anyNA(edges$source) || anyNA(edges$target) || anyNA(edges$weight)) {
      stop("high_school edges.csv contains non-integer source/target/weight values")
    }

    # KONECT indices are 0-based in this file; convert to 1-based for R matrices.
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }

    n_nodes <- max(c(edges$source, edges$target))
    node_ids <- as.character(seq_len(n_nodes))
    A <- matrix(0L, n_nodes, n_nodes, dimnames = list(node_ids, node_ids))

    # Directed weighted nominations: A_ij = number of surveys where i named j.
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]
      j <- edges$target[r]
      w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w
    }
    diag(A) <- 0L
  } else stop("Unknown dataset.")
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  colnames(A) <- rownames(A)
  A
}

# ======================================================================
# Fit + score (per dataset, per {model,K}), unified with sim pipeline
# ======================================================================

fit_and_score_one <- function(A, dataset, fit_model, K,
                              n_iter, burn, thin, verbose, seed,
                              osbm_partition_prior = "OCRP",
                              theta_ocrp_val = 0.5,
                              alpha_crp_val = 0.5,
                              m_item_triplets = 2000L,
                              alpha_empirical = 0.5,
                              ppc_run_dir = here::here("output", "application", "ppc"),
                              do_ppc = FALSE,
                              ppc_draws = 200L,
                              ppc_triples = 2000L) {
  
  Sys.setenv("OMP_NUM_THREADS"="1","MKL_NUM_THREADS"="1",
             "OPENBLAS_NUM_THREADS"="1","VECLIB_MAXIMUM_THREADS"="1")
  options(mc.cores = 1)
  source('./core/ppc_checks.R')
  source("./helper_folder/sim_study_helper.R")
  source("./helper_folder/transitivity_check_helper.R")
  source("./core/DCSBM_varK.R")
  source("./helper_folder/helper.R")
  source("./helper_folder/SST_helpers.R")
  source("./helper_folder/WST_helpers.R")
  source("./helper_folder/Hyper_setup.R")
  source("./core/my_best_try_so_far.R")

  if (!exists(".empty_diag_vec", mode = "function")) {
    .empty_diag_vec <- function() {
      base_stats <- c("mean", "lo", "hi", "ess", "mcse")
      diag_prefixes <- c(
        "thetaW_block_model", "thetaS_block_model",
        "thetaW_block_emp", "thetaS_block_emp",
        "thetaW_item_emp", "thetaS_item_emp",
        "thetaW_block_model_all", "thetaS_block_model_all",
        "thetaW_block_emp_all", "thetaS_block_emp_all",
        "thetaW_item_emp_all", "thetaS_item_emp_all"
      )
      diag_names <- unlist(lapply(diag_prefixes, function(prefix) {
        paste0(prefix, "_", base_stats)
      }))
      tail_names <- c(
        "prem_block_model_avg", "prem_block_emp_avg", "prem_item_emp_avg",
        "coverage_block_model_avg", "coverage_block_emp_avg", "coverage_item_emp_avg",
        "transitive_triads_mean", "transitive_triads_lo", "transitive_triads_hi",
        "cycle_mass_weighted_mean",
        "min_backward_weight_mean", "min_backward_weight_norm_mean",
        "hierarchy_energy_mean", "hierarchy_energy_lo", "hierarchy_energy_hi",
        "curl_energy_mean",
        "p_post_wst", "p_prior_wst", "bf_wst_0", "log_bf_wst_0",
        "p_post_sst", "p_prior_sst", "bf_sst_0", "log_bf_sst_0",
        "bf_sst_wst", "log_bf_sst_wst", "bf_prior_mc",
        "violation_rate_mean", "violation_rate_lo", "violation_rate_hi",
        "violation_count_mean", "cross_mass_mean",
        "violation_rate_zhat", "violation_count_zhat"
      )
      all_names <- c(diag_names, tail_names)
      out <- rep(NA_real_, length(all_names))
      names(out) <- all_names
      out
    }
  }
  
  t0 <- proc.time()
  ok <- TRUE; msg <- NA_character_
  out <- NULL
  K_init_fit <- nrow(A)
  K_expected_fit <- max(2L, min(as.integer(ceiling(sqrt(nrow(A)))), K_init_fit))

  if (isTRUE(verbose)) {
    cat(glue::glue("[FIT] dataset={dataset} model={fit_model} seed={seed} n_iter={n_iter} burn={burn} thin={thin} K_init={K_init_fit} K_expected={K_expected_fit}\n"))
  }
  
  ## ---- make hypers available for ALL models ----
  ordering_prior_mode <- "equivalence_class"
  hypers <- get_principled_hypers_v2(
    A          = A,
    K_expected = K_expected_fit,
    c_kappa    = 3
  )
  hypers$tau0 <- max(hypers$tau0, 0.2)
  gamma_gn_tuned <- hypers$gamma
  
  tryCatch({
    
    if (fit_model %in% c("WST","SST")) {
      
      hypers$b_kappa <- 1
      hypers$b_eta <- 1
      
      out <- modular_osbm_sampler(
        A = A,
        K = nrow(A),
        n_iter = n_iter,
        burn   = burn,
        thin   = thin,
        verbose = verbose,
        a_kappa = 1,
        b_kappa = hypers$b_kappa,
        a_eta   = 1,
        b_eta   = hypers$b_eta,
        psi_constraint = fit_model ,
        partition_prior = osbm_partition_prior,
        theta_ocrp      = theta_ocrp_val,
        eta_identifiability = "block_sum",
        mu0     = hypers$mu0,
        sigma0  = hypers$sigma0,
        tau0    = hypers$tau0,
        use_mixing_moves = TRUE,
        seed = seed
      )
      
      if (is.null(out)) {
        stop(sprintf("[%s|%s K=%d] modular_osbm_sampler failed",
                     dataset, fit_model, K))
      }
      
      
    } else if (fit_model %in% c("DC-SBM","DCSBM")) {
      
      fit_model <- "DCSBM"
      b_eta = 1
      b_lambda = 1

      out <- fit_dcsbm_gibbs_gnedin(
        as.matrix(A),
        K_init     = nrow(A),
        priors = list(
          a_eta       = 1,  b_eta = b_eta,   # <- fixed
          a_lambda    = 1, b_lambda = b_lambda,
          alpha_crp = alpha_crp_val,
          partition_prior = "CRP"
        ),
        iters   = n_iter,
        burn_in = burn,
        thin    = thin,
        verbose = if (isTRUE(verbose)) 50 else 0,
        seed    = seed
      )
      
    } else stop("Unknown fit_model: ", fit_model)
    
  }, error = function(e) { ok <<- FALSE; msg <<- conditionMessage(e) })
  
  elapsed <- (proc.time() - t0)[["elapsed"]]

  if (isTRUE(verbose) && ok) {
    cat(glue::glue("[FIT] completed dataset={dataset} model={fit_model} elapsed_sec={round(elapsed, 2)}\n"))
  }
  
  if (!ok) {
    return(tibble::tibble(
      dataset = dataset, gen_model = NA_character_, fit_model = fit_model, K = K,
      gamma_gn = gamma_gn_tuned,
      ordering_prior_mode = ordering_prior_mode,
      difficulty = NA_character_, rep = NA_integer_,
      ari = NA_real_, vi = NA_real_,
      elpd_loo = NA_real_, p_loo = NA_real_, looic = NA_real_,
      pareto_k_max = NA_real_, pareto_k_share = NA_real_,
      ess_eta = NA_real_, ac_eta = NA_real_,
      ess_psi = NA_real_, ac_psi = NA_real_,
      ess_kappa = NA_real_, ac_kappa = NA_real_,
      !!!as.list(.empty_diag_vec()),
      time_sec = elapsed,
      note = msg
    ))
  }
  
  # --- Relabel BEFORE getting z_hat
  out_relab <- if(fit_model=='SST'){
    out
  }else if(fit_model == "WST") {
    out
  } else {
    relabel_by_ECR(mcmc_out = out,model = 'DCSBM')
  }
  
  # ------------------------------------------------------------------
  # SAVE FIT + PPC (only when requested)
  # ------------------------------------------------------------------
  if (do_ppc && fit_model %in% c("WST", "SST", "DCSBM")) {

    # One clean output folder per dataset/model under posterior_predictive_checks/
    ppc_dir <- file.path(ppc_run_dir, dataset, fit_model)
    fs::dir_create(ppc_dir)

    # Save observed adjacency once per dataset
    A_path <- file.path(ppc_run_dir, dataset, paste0(dataset, "_A_obs.rds"))
    if (!file.exists(A_path)) {
      fs::dir_create(dirname(A_path))
      saveRDS(A, A_path)
    }

    # Save fitted object
    fit_path <- file.path(ppc_dir,
                          paste0(dataset, "_", fit_model, "_fit.rds"))
    saveRDS(out_relab, fit_path)

    # Run PPC and save outputs
    ppc_tag <- paste0("_", dataset)
    ppc_res <- .run_ppc_model(
      A_obs     = A,
      out_relab = out_relab,
      regime    = fit_model,
      n_draws   = ppc_draws,
      m_triples = ppc_triples,
      seed      = seed + 999,
      out_dir   = ppc_dir,
      tag       = ppc_tag
    )

    saveRDS(ppc_res,
            file.path(ppc_dir, paste0(dataset, "_", fit_model, "_ppc.rds")))
    readr::write_csv(ppc_res$pvals,
                     file.path(ppc_dir, paste0(dataset, "_", fit_model, "_ppc_pvals.csv")))
    readr::write_csv(ppc_res$draw_stats,
                     file.path(ppc_dir, paste0(dataset, "_", fit_model, "_ppc_drawstats.csv")))

    cat(glue::glue("[PPC] Saved fit:  {fit_path}\n"))
    cat(glue::glue("[PPC] Saved plots: {ppc_dir}\n"))
  }
  if (is.matrix(out_relab$z)) {
    z_chain <- out_relab$z              # S x n
  } else {
    z_chain <- do.call(rbind, out_relab$z)  # list of length S -> S x n
  }
  S <- nrow(z_chain)
  K_occ_chain <- apply(z_chain, 1L, function(z) length(unique(z[!is.na(z)])))
  salso_k_cap <- as.integer(max(
    3L,
    min(
      ncol(z_chain),
      stats::quantile(K_occ_chain, probs = 0.99, na.rm = TRUE, type = 1) + 1L
    )
  ))
  
  # --- Partition estimate on relabelled chains
  z_hat <- salso::salso(
    z_chain,
    loss = salso::VI(),
    nRuns = 1L,
    nCores = 1L,
    maxNClusters = salso_k_cap,
    maxZealousAttempts = 100000L
  )
  K = length(unique(z_hat))
  # --- Pointwise log-lik (edgewise), unified naming
  n <- nrow(A)
  DI <- build_dyad_index(n)
  
  
  LL <- if (fit_model == "WST") {
    loglik_matrix_modular(A, out, regime = fit_model, dyad_index = DI)
  } else if (fit_model == "SST") {
    loglik_matrix_modular(A, out, regime = fit_model, dyad_index = DI)
  } else {
    loglik_matrix_dcsbm(A = as.matrix(A), dcsbm_out = out, dyad_index = DI)
  }
  
  # --- LOO (PSIS)
  loo_fit <- loo::loo(LL)
  elpd_loo <- loo_fit$estimates["elpd_loo","Estimate"]
  p_loo    <- loo_fit$estimates["p_loo","Estimate"]
  looic    <- loo_fit$estimates["looic","Estimate"]
  pareto_k_max   <- max(loo_fit$diagnostics$pareto_k)
  pareto_k_share <- mean(loo_fit$diagnostics$pareto_k > 0.7)
  
  # # --- ESS / AC (parameters) — only for OSBM where shapes are standard
  # ess_eta <- ac_eta <- ess_psi <- ac_psi <- ess_kappa <- ac_kappa <- NA_real_
  # if (fit_model %in% c("WST","SST")) {
  #   ess_eta <- mean(coda::effectiveSize(coda::mcmc(out_relab$eta)))
  #   ac_eta  <- mean(coda::autocorr.diag(coda::mcmc(out_relab$eta))[,2])
  #   kappamat <- matrix(out_relab$kappa, ncol = K*K)
  #   ess_kappa <- mean(coda::effectiveSize(coda::mcmc(kappamat)))
  #   ac_kappa  <- mean(coda::autocorr.diag(coda::mcmc(kappamat))[,2])
  #   if (!is.null(out_relab$psi)) {
  #     if (is.matrix(out_relab$psi)) {
  #       ess_psi <- mean(coda::effectiveSize(coda::mcmc(out_relab$psi)))
  #       ac_psi  <- mean(coda::autocorr.diag(coda::mcmc(out_relab$psi))[,2])
  #     } else {
  #       psimat <- if (length(dim(out_relab$psi)) == 3)
  #         matrix(out_relab$psi, ncol = K*K) else as.matrix(out_relab$psi)
  #       ess_psi <- mean(coda::effectiveSize(coda::mcmc(psimat)))
  #       ac_psi  <- mean(coda::autocorr.diag(coda::mcmc(psimat))[,2])
  #     }
  #   }
  # }
  
  # --- Order-conformity summaries (model-based + empirical)
  diag_note <- NA_character_
  diag_summ <- tryCatch({
    if (!any(is.na(z_hat))) {
      if (fit_model %in% c("WST","SST")) {
        summarise_osbm_diagnostics(out_relab, regime = fit_model, K_max_hint = K,
                                   z_hat = z_hat, n = n, m_items = m_item_triplets,
                                   alpha = alpha_empirical, A = A,seed_block = 123,
                                   method_order = 'bt',T_block = NULL)
      } else {
        summarise_dcsbm_diagnostics(fit = out_relab, z_hat = z_hat, K = K,
                                    n = n, m_items = m_item_triplets,
                                    alpha = alpha_empirical, A = A,seed_block = 123,method_order = 'bt',
                                    T_block = NULL)
      }
    } else {
      .empty_diag_vec()
    }
  }, error = function(e) {
    diag_note <<- paste0("diag_error: ", conditionMessage(e))
    if (isTRUE(verbose)) {
      cat(glue::glue("[WARN] diagnostics failed dataset={dataset} model={fit_model}: {conditionMessage(e)}\n"))
    }
    .empty_diag_vec()
  })

  if (fit_model == "DCSBM") {
    key_diag_cols <- c("thetaW_block_model_mean", "violation_rate_mean", "transitive_triads_mean")
    if (all(is.na(diag_summ[key_diag_cols]))) {
      na_msg <- "diag_all_na_dcsbm"
      diag_note <- if (is.na(diag_note)) na_msg else paste(diag_note, na_msg, sep = "; ")
      if (isTRUE(verbose)) {
        cat(glue::glue("[WARN] diagnostics all-NA dataset={dataset} model={fit_model}\n"))
      }
    }
  }
  
  K_hat_VI     <- length(unique(z_hat))
  z_hat_alt    <- salso::salso(
    z_chain,
    loss = salso::VI(a = 0.5),
    nRuns = 1L,
    nCores = 1L,
    maxNClusters = salso_k_cap,
    maxZealousAttempts = 100000L
  )
  K_hat_VI_a05 <- length(unique(z_hat_alt))
  
  # --- Assemble a single row
  row <- tibble::tibble(
    dataset = dataset,
    gen_model = NA_character_,            # not applicable in applications
    fit_model = fit_model,
    K = K,
    gamma_gn = gamma_gn_tuned,
    ordering_prior_mode = ordering_prior_mode,
    K_hat_VI = K_hat_VI,
    K_hat_VI_a05 = K_hat_VI_a05,
    difficulty = NA_character_,
    rep = NA_integer_,
    ari = NA_real_, vi = NA_real_,        # no truth in applications
    elpd_loo = elpd_loo, p_loo = p_loo, looic = looic,
    pareto_k_max = pareto_k_max, pareto_k_share = pareto_k_share,
    !!!as.list(diag_summ),
    time_sec = elapsed,
    note = diag_note
  )
  
  row
}

# ======================================================================
# Main sweep over datasets (sequential) and over {model,K} (parallel)
# ======================================================================

run_start_time <- Sys.time()
run_stamp <- format(run_start_time, "%Y%m%d_%H%M%S")
run_id <- paste0("application_run_", run_stamp)
run_dir <- here::here("output", "application", "raw", run_id)
fs::dir_create(run_dir)

# ---- Log to file + console ----
log_path <- file.path(run_dir, "run_log.txt")
log_con  <- file(log_path, open = "wt")
sink(log_con, type = "output",  split = TRUE)
sink(log_con, type = "message", append = TRUE)

# PPC goes into a single stable folder (no timestamp) so repeated
# runs overwrite rather than spawning duplicate directories.
ppc_run_dir <- here::here("output", "application", "ppc")
fs::dir_create(ppc_run_dir)

cat(glue("Run directory: {run_dir}\n"))
cat(glue("Posterior predictive directory: {ppc_run_dir}\n"))

datasets_default <- c('citations_data','mountain_goats','macaques_data','high_school')
datasets_env <- trimws(Sys.getenv("APP_DATASETS", unset = ""))
datasets <- if (nzchar(datasets_env)) strsplit(datasets_env, ",")[[1]] else datasets_default

models   <- c("SST",'WST','DC-SBM')
Ks       <- c(3:10)

# ---- Parallel plan: pick number of workers based on grid size ----


grid <- expand.grid(fit_model = models, stringsAsFactors = FALSE)
n_tasks <- nrow(grid)

slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")))
slurm_ntasks <- suppressWarnings(as.integer(Sys.getenv("SLURM_NTASKS", unset = "")))
env_workers <- suppressWarnings(as.integer(Sys.getenv("APP_WORKERS", unset = "")))

if (!is.na(env_workers) && env_workers >= 1L) {
  max_cores <- env_workers
} else if (!is.na(slurm_cpus) && slurm_cpus >= 1L) {
  max_cores <- slurm_cpus
} else if (!is.na(slurm_ntasks) && slurm_ntasks >= 1L) {
  max_cores <- slurm_ntasks
} else {
  max_cores <- parallelly::availableCores()
}

workers   <- max(1L, min(max_cores, n_tasks))
if (workers <= 1L) {
  plan(sequential)
} else {
  plan(multisession, workers = workers)
}

cat(sprintf("Using %d/%d cores (n_tasks=%d)\n", workers, max_cores, n_tasks))

# Hyperparams (shared with sim study style)
n_iter <- as.integer(Sys.getenv("APP_N_ITER", unset = "10000"))
burn <- as.integer(Sys.getenv("APP_BURN", unset = "3000"))
thin <- as.integer(Sys.getenv("APP_THIN", unset = "2"))
theta_ocrp_val <- as.numeric(Sys.getenv("APP_THETA_OCRP", unset = "0.5"))
alpha_crp_val  <- as.numeric(Sys.getenv("APP_ALPHA_CRP",  unset = "0.5"))
app_prior_env <- Sys.getenv("APP_OSBM_PARTITION_PRIOR",
                            unset = Sys.getenv("APP_PARTITION_PRIOR", unset = "OCRP"))
osbm_partition_prior <- .normalize_app_ordered_prior(app_prior_env)
verbose <- tolower(trimws(Sys.getenv("APP_VERBOSE", unset = "true"))) %in% c("1", "true", "yes")

cat(sprintf("OSBM partition prior: %s (theta=%.2f); DCSBM partition prior: CRP (alpha=%.2f)\n",
            osbm_partition_prior, theta_ocrp_val, alpha_crp_val))
cat(sprintf("Verbose logging: %s\n", ifelse(verbose, "ON", "OFF")))

if (!is.finite(n_iter) || n_iter < 1L) n_iter <- 10000L
if (!is.finite(burn) || burn < 0L || burn >= n_iter) burn <- max(0L, floor(n_iter / 5))
if (!is.finite(thin) || thin < 1L) thin <- 2L

big_res <- list()
do_ppc_global <- tolower(trimws(Sys.getenv("APP_DO_PPC", unset = "true"))) %in% c("1", "true", "yes")
# "all" = run PPC for every dataset; or set to a single name e.g. "mountain_goats"
ppc_dataset <- trimws(Sys.getenv("APP_PPC_DATASET", unset = "all"))

for (dataset in datasets) {
  cat("\n================ Dataset:", dataset, "================\n")
  A <- choose_dataset(dataset)
  chosen_dataset = dataset
  
  # Parallel over rows of grid within this dataset (or sequential when workers=1)
  run_one <- function(i) {
    set.seed(10000 + i)
    fm <- grid$fit_model[i]
    K <- max(2L, ceiling(sqrt(nrow(A))))
    
    do_ppc_here <- (do_ppc_global &&
                     (ppc_dataset == "all" || dataset == ppc_dataset) &&
                     fm %in% c("SST", "WST", "DC-SBM", "DCSBM"))
    
    fit_and_score_one(
      A, dataset, fm, K,
      n_iter, burn, thin, verbose,
      seed = 100 + K,
      osbm_partition_prior = osbm_partition_prior,
      theta_ocrp_val = theta_ocrp_val,
      alpha_crp_val = alpha_crp_val,
      m_item_triplets = 2000L,
      alpha_empirical = 0.5,
      ppc_run_dir = ppc_run_dir,
      do_ppc = do_ppc_here,
      ppc_draws = 200L,
      ppc_triples = 2000L
    )
  }

  res_list <- if (workers <= 1L) {
    lapply(seq_len(nrow(grid)), run_one)
  } else {
    future.apply::future_lapply(
      seq_len(nrow(grid)),
      run_one,
      future.seed = TRUE,
      future.stdout = TRUE
    )
  }
  
  dataset_res <- bind_rows(res_list)
  
  # Save per-dataset table (optional)
  out_csv <- file.path(run_dir, glue("{dataset}_applications_summary_{run_id}.csv"))
  readr::write_csv(dataset_res, out_csv)
  cat(glue("Saved: {out_csv}\n"))
  
  big_res[[dataset]] <- dataset_res
}

# ----- One big CSV (same schema as sim study + dataset) -----------------
app_res_summary <- bind_rows(big_res) %>%
  arrange(dataset, fit_model, K)

summary_csv <- file.path(run_dir, glue("applications_results_summary_{run_id}.csv"))
summary_rds <- file.path(run_dir, glue("applications_results_summary_{run_id}.rds"))
readr::write_csv(app_res_summary, summary_csv)
saveRDS(app_res_summary, summary_rds)

git_hash <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1],
  error = function(e) NA_character_
)
if (!is.character(git_hash) || !length(git_hash) || is.na(git_hash) || !nzchar(git_hash)) {
  git_hash <- "NA"
}

run_end_time <- Sys.time()
model_hyper_lines <- c(
  sprintf("hyperparams.SST: a_kappa=1,b_kappa=1,a_eta=1,b_eta=1,psi_constraint=SST,partition_prior=%s,theta_ocrp=%.2f", tolower(osbm_partition_prior), theta_ocrp_val),
  sprintf("hyperparams.WST: a_kappa=1,b_kappa=1,a_eta=1,b_eta=1,psi_constraint=WST,partition_prior=%s,theta_ocrp=%.2f", tolower(osbm_partition_prior), theta_ocrp_val),
  sprintf("hyperparams.DCSBM: a_eta=1,b_eta=1,a_lambda=1,b_lambda=1,partition_prior=CRP,alpha_crp=%.2f", alpha_crp_val),
  "hyperparams.shared: ordering_prior_mode=equivalence_class,get_principled_hypers_v2(K_expected=ceiling(sqrt(n_items)),c_kappa=3),tau0_floor=0.2,K_init=n_items"
)
manifest_lines <- c(
  paste0("run_id: ", run_id),
  paste0("script: application.R"),
  paste0("run_dir: ", run_dir),
  paste0("git_commit_short: ", git_hash),
  paste0("start_time: ", format(run_start_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("end_time: ", format(run_end_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("elapsed_seconds: ", round(as.numeric(difftime(run_end_time, run_start_time, units = "secs")), 3)),
  paste0("datasets: ", paste(datasets, collapse = ",")),
  paste0("fit_models: ", paste(models, collapse = ",")),
  model_hyper_lines,
  paste0("n_iter: ", n_iter),
  paste0("burn: ", burn),
  paste0("thin: ", thin),
  paste0("K_init_policy: n_items"),
  paste0("K_expected_policy: ceiling(sqrt(n_items))"),
  paste0("tau0_floor: 0.2"),
  paste0("ordering_prior_mode: equivalence_class"),
  paste0("summary_csv: ", summary_csv),
  paste0("summary_rds: ", summary_rds),
  paste0("ppc_run_dir: ", ppc_run_dir),
  paste0("ppc_enabled: ", do_ppc_global),
  paste0("ppc_dataset: ", ppc_dataset),
  paste0("n_rows_summary: ", nrow(app_res_summary)),
  paste0("n_error_notes: ", sum(!is.na(app_res_summary$note) & nzchar(app_res_summary$note))),
  paste0("R_version: ", R.version.string),
  paste0("platform: ", R.version$platform),
  paste0("packages: ", paste(sapply(c("coda","mcclust","mcclust.ext","salso","loo","dplyr","future.apply"), function(p) paste0(p, "=", packageVersion(p))), collapse = ", "))
)
writeLines(manifest_lines, con = file.path(run_dir, "run_manifest.txt"))

cat(glue("\nSaved: {summary_csv}\n"))

cat("\n Applications run complete\n")

best_csv = read.csv(summary_csv)

app_res_summary%>%
  group_by(dataset)%>%
  slice_min(looic)%>%
  select(looic, fit_model,K)

# ---- close log ----
sink(type = "message")
sink(type = "output")
close(log_con)
