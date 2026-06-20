############################################################
## 0. Libraries & helpers
############################################################

# Canonical fixed-K simulation driver for the bundled paper results.
# scripts/simulation/sim_study.R is the older legacy driver.

get_repo_root <- function(relative_to_script = "../..") {
  env_root <- Sys.getenv("TRANSITIVE_SBM_BUNDLE_ROOT", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
  }
  wd <- getwd()
  if (dir.exists(file.path(wd, "scripts")) && dir.exists(file.path(wd, "helper_folder"))) {
    return(normalizePath(wd, winslash = "/", mustWork = TRUE))
  }
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grepl("^--file=", cmd)]
  if (!length(file_arg)) {
    return(normalizePath(wd, winslash = "/", mustWork = TRUE))
  }
  script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg[1L])),
                               winslash = "/", mustWork = TRUE)
  normalizePath(file.path(dirname(script_path), relative_to_script),
                winslash = "/", mustWork = TRUE)
}

setwd(get_repo_root())

library(stats)
library(salso)
library(BayesLogit)
library(truncnorm)
library(fossil)
library(loo)
library(mcclust)
library(dplyr)
library(tidyr)
library(ggplot2)
library(kableExtra)
library(Matrix)
library(parallel)
library(Rcpp)

# --- Load your existing helpers (unchanged) ---

source("./helper_folder/simulation/simulation_study_helpers.R")
source("./helper_folder/models/ordered_sbm/shared_sampler_helpers.R")
source("./helper_folder/diagnostics/transitivity_diagnostics.R")
source("./core/DCSBM_varK.R")        # fit_dcsbm, dcsbm_relabel, etc.
source("./helper_folder/models/ordered_sbm/sst_helpers.R")
source("./helper_folder/models/ordered_sbm/wst_helpers.R")
source("./helper_folder/config/hyperparameter_setup.R")
source("./core/transitive_sbm_sampler.R")

############################################################
## 1. Generative priors: kappa, psi_WST, psi_SST (mean/var)
############################################################

# Symmetric kappa: Gamma with given mean and variance
sample_kappa_prior <- function(K, mean, var) {
  if (mean <= 0) stop("'mean' must be > 0 for Gamma.")
  if (var  <= 0) stop("'var'  must be > 0 for Gamma.")
  
  shape <- mean^2 / var    # alpha
  rate  <- mean / var      # beta
  
  KAP <- matrix(0, K, K)
  for (k in 1:K) {
    for (l in k:K) {
      KAP[k, l] <- rgamma(1L, shape = shape, rate = rate)
    }
  }
  KAP[lower.tri(KAP)] <- t(KAP)[lower.tri(KAP)]
  KAP
}

# WST psi: block-pair specific, skew-symmetric
# psi_{kℓ} (k<ℓ) ~ N^+(mu_psi, sigma_psi^2)
sample_psi_wst_prior <- function(K, mu_psi, var_psi) {
  if (var_psi <= 0) stop("'var_psi' must be > 0")
  sigma_psi <- sqrt(var_psi)
  
  PSI <- matrix(0, K, K)
  PSI[upper.tri(PSI)] <- truncnorm::rtruncnorm(
    n    = K * (K - 1) / 2,
    a    = 0,
    mean = mu_psi,
    sd   = sigma_psi
  )
  PSI[lower.tri(PSI)] <- -t(PSI)[lower.tri(PSI)]
  PSI
}

# SST psi: distance-based, psi_d = sum_{i<=d} delta_i, delta_i ~ N^+(mu_psi, tau^2)
# Both mu and tau are scaled by 1/(K-1) so that E[psi_{K-1}] is constant across K.
sample_psi_sst_prior <- function(K, mu_psi, var_psi) {
  if (K <= 1) return(numeric(0))
  if (var_psi <= 0) stop("'var_psi' must be > 0")
  
  mu_delta  <- mu_psi  / (K - 1)
  tau_delta <- sqrt(var_psi) / (K - 1)
  
  deltas <- truncnorm::rtruncnorm(
    n    = K - 1,
    a    = 0,
    mean = mu_delta,
    sd   = tau_delta
  )
  cumsum(deltas)
}

demo_osbm_partition_prior <- toupper(trimws(Sys.getenv("DEMOKVAR_OSBM_PARTITION_PRIOR", unset = "OCRP")))
if (!demo_osbm_partition_prior %in% c("OCRP", "ROCRP")) {
  stop("DEMOKVAR_OSBM_PARTITION_PRIOR must be OCRP or ROCRP.", call. = FALSE)
}
demo_theta_ocrp <- as.numeric(Sys.getenv("DEMOKVAR_THETA_OCRP", unset = "0.5"))
demo_alpha_crp  <- as.numeric(Sys.getenv("DEMOKVAR_ALPHA_CRP",  unset = "0.5"))

############################################################
## 3. One replicate: generate + fit all three models
############################################################

run_replicate_crossfit <- function(
    gen_model,
    K_true,
    n,
    kappa_mean,
    kappa_var,
    psi_mean,
    psi_var,
    hierch_label,
    seed,
    n_iter = 5000,
    burn   = 1000,
    thin   = 2
) {
  set.seed(seed)
  
  # --- 1. Generate data from OSBM prior ---
  z_true     <- rep(seq_len(K_true), length.out = n)
  kappa_true <- sample_kappa_prior(K_true, mean = kappa_mean, var = kappa_var)
  eta_true   <- runif(n, 0.9, 1.1)
  
  if (gen_model == "WST") {
    psi_true <- sample_psi_wst_prior(K_true, mu_psi = psi_mean, var_psi = psi_var)
    A        <- simulate_osbm(n, K_true, z_true, eta_true, kappa_true, psi_true, "WST")
  } else {
    psi_true <- sample_psi_sst_prior(K_true, mu_psi = psi_mean, var_psi = psi_var)
    A        <- simulate_osbm(n, K_true, z_true, eta_true, kappa_true, psi_true, "SST")
  }
  A <- as.matrix(A)
  
  # --- 2. Model hypers (principled Gibbs calibration) ---
  K_init_fit <- nrow(A)
  K_expected_fit <- K_true
  
  hypers <- get_principled_hypers_v2(
    A = A,
    K_expected = K_expected_fit,
    c_kappa = 3,
    ordering_prior_mode = "equivalence_class"
  )
  hypers$tau0 <- max(hypers$tau0, 0.2)
  gamma_gn_tuned <- hypers$gamma

  # K is learned from an overfitted start; keep volume/degree priors shared.
  hypers_wst <- hypers
  hypers_wst$b_kappa <- 1
  hypers_wst$b_eta <- 1

  hypers_sst <- hypers
  hypers_sst$b_kappa <- 1
  hypers_sst$b_eta <- 1
  
  # --- 3. Fit models: WST, SST, DC-SBM ---
  models_to_fit <- c("WST", "SST", "DC-SBM")
  results_list  <- list()
  DI            <- build_dyad_index(nrow(A))
  
  for (fit_model in models_to_fit) {
    loop_seed <- seed + which(models_to_fit == fit_model)
    set.seed(loop_seed)
    start_time <- proc.time()
    
    mcmc_out_for_diag <- NULL
    regime_for_diag   <- if (fit_model == "SST") "SST" else "WST"
    
    # --- 3a. DC-SBM: initialize from K = n and learn the occupied K ---
    if (fit_model == "DC-SBM") {
      
      fit <- fit_dcsbm_gibbs_gnedin(
        as.matrix(A),
        K_init     = K_init_fit,
        priors = list(
          a_eta         = 1,
          b_eta         = 1,
          a_lambda      = 1,
          b_lambda      = 1,
          alpha_crp     = demo_alpha_crp,
          partition_prior = "CRP"
        ),
        iters   = n_iter,
        burn_in = burn,
        thin    = thin,
        verbose = 200,
        seed    = loop_seed
      )
      
      relab_out <- relabel_chain_by_block_score(fit$z, fit$lambda, A, A + t(A))
      Z_chain_relab    <- relab_out$Z      # matrix S x n or list of length S
      Lambda_chain_rel <- relab_out$Lambda # list length S, each K_s x K_s
      
      
      # Coerce Z_chain_relab to an S x n matrix for SALSO and diagnostics
      if (is.matrix(Z_chain_relab)) {
        z_chain <- Z_chain_relab              # S x n
      } else {
        z_chain <- do.call(rbind, Z_chain_relab)  # list of length S -> S x n
      }
      S <- nrow(z_chain)
      
      
      LL <- loglik_matrix_dcsbm(A, fit, dyad_index = DI)
      
      ## --- Pseudo-psi + ordered diagnostics (treat DC-SBM as WST on blocks) ---
      psi_list <- vector("list", S)
      z_list   <- vector("list", S)
      
      for (s in seq_len(S)) {
        kap   <- Lambda_chain_rel[[s]]          # K_s x K_s
        K_curr <- nrow(kap)
        vol   <- kap + t(kap)
        
        rho <- matrix(0.5, nrow = K_curr, ncol = K_curr)
        mask <- vol > 0
        rho[mask] <- kap[mask] / vol[mask]
        rho_clamped <- pmin(pmax(rho, 1e-6), 1 - 1e-6)
        
        psi_list[[s]] <- qlogis(rho_clamped)   # WST-style: full K_s x K_s logits
        z_list[[s]]   <- z_chain[s, ]
      }
      
      mcmc_out_for_diag <- list(psi = psi_list, z = z_list)
      
      
    } else {
      # --- 3b. WST / SST OSBM (unchanged sampler) ---

      # #check SST
      # res <- run_fast_bruteforce_tests(A = A, N = A+t(A), z = z_true, psi_vec = psi_true, 
      #                                  i_set = sample(seq_along(z_true), 5),)
      # 
      # #check WST
      # run_fast_bruteforce_tests_wst(
      #   A = A, N = A+t(A), z = z_true, psi_mat = psi_true,
      #   mu0 = 0, sigma0 = 1,
      #   tol = 1e-8, verbose = TRUE
      # )
      
      fail_on_osbm_error <- tolower(Sys.getenv("DEMOKVAR_FAIL_ON_OSBM_ERROR", unset = "false")) %in% c("1", "true", "yes")

      fit_hypers <- if (fit_model == "SST") hypers_sst else hypers_wst
      out <- tryCatch({
        modular_osbm_sampler(
          A = A,
          K = K_init_fit,
          n_iter = n_iter,
          burn   = burn,
          thin   = thin,
          verbose = TRUE,
          a_kappa = 1,
          b_kappa = fit_hypers$b_kappa,
          a_eta   = 1,
          b_eta   = fit_hypers$b_eta,
          psi_constraint = fit_model,     # "WST" or "SST"
          partition_prior = demo_osbm_partition_prior,
          theta_ocrp      = demo_theta_ocrp,
          eta_identifiability = "block_sum",
          use_mixing_moves = TRUE,
          gamma_gn        = gamma_gn_tuned,
          mu0   = fit_hypers$mu0,
          sigma0 = fit_hypers$sigma0,
          tau0   = fit_hypers$tau0,
          seed = loop_seed
        )
      }, error = function(e) {
        message(sprintf(
          "OSBM fit failed | gen=%s fit=%s K_true=%d hierch=%s kappa_mean=%.3f psi_mean=%.3f seed=%d | %s",
          gen_model, fit_model, K_true, hierch_label, kappa_mean, psi_mean, loop_seed, conditionMessage(e)
        ))
        if (isTRUE(fail_on_osbm_error)) {
          stop(e)
        }
        NULL
      })
      
      if (is.null(out)) next
      
      z_chain <- do.call(rbind, lapply(out$z, as.integer))
      LL      <- loglik_matrix_modular(A, out, regime = fit_model, dyad_index = DI)
      mcmc_out_for_diag <- out
    }
    
    time_taken <- (proc.time() - start_time)["elapsed"]
    
    # --- 4. Clustering metrics ---
    z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L)
    ari   <- fossil::adj.rand.index(z_hat,z_true )
    vi    <- mcclust::vi.dist(z_hat, z_true)
    K_hat <- length(unique(z_hat))
    K_count = function(x)length(unique(x))
    K_p.mean = mean(apply(z_chain, MARGIN = 1, K_count))
    
    
    loo_res <- tryCatch(loo::loo(LL),
                        error = function(e) NULL)
    
    if (is.null(loo_res)) {
      elpd   <- NA_real_
      looic  <- NA_real_
      pk_bad <- NA_real_
    } else {
      elpd   <- loo_res$estimates["elpd_loo", "Estimate"]
      looic  <- loo_res$estimates["looic",   "Estimate"]
      pk_bad <- mean(loo_res$diagnostics$pareto_k > 0.7)
    }
    
    
    diag_res <- summarise_osbm_diagnostics(
      out_relab = mcmc_out_for_diag,
      regime    = regime_for_diag,
      K_max_hint = K_true,
      z_hat     = z_hat,
      n         = n,
      A         = A,
      m_items   = 500
    )
    
    res_row <- data.frame(
      gen_model  = gen_model,
      fit_model  = fit_model,
      K_true     = K_true,
      K_hat      = K_hat,
      K_p.mean   = K_p.mean,
      hierch     = hierch_label,           # weak / strong (via psi mean)
      kappa_mean = kappa_mean,
      kappa_var  = kappa_var,
      psi_mean   = psi_mean,
      psi_var    = psi_var,
      rate_kap   = kappa_mean,            # for backwards-compatible plots
      ari        = ari,
      vi         = vi,
      elpd       = elpd,
      looic      = looic,
      pk_bad     = pk_bad,
      time       = time_taken,
      seed       = loop_seed,
      stringsAsFactors = FALSE
    )
    print(res_row)
    
    results_list[[length(results_list) + 1L]] <-
      cbind(res_row, as.data.frame(as.list(diag_res)))
  }
  
  do.call(rbind, results_list)
}

############################################################
## 4. Simple file lock for progressive writing
############################################################

# Very simple directory-based lock: no extra packages needed.
with_ilock <- function(lock_base, code) {
  lock_dir <- paste0(lock_base, ".lock")
  repeat {
    ok <- dir.create(lock_dir, showWarnings = FALSE)
    if (ok) break
    Sys.sleep(runif(1, 0.01, 0.05))
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  code()
}

append_results_with_lock <- function(res_df, out_file) {
  with_ilock(out_file, function() {
    if (file.exists(out_file)) {
      old <- read.csv(out_file, stringsAsFactors = FALSE)
      new <- dplyr::bind_rows(old, res_df)
    } else {
      new <- res_df
    }
    write.csv(new, out_file, row.names = FALSE)
  })
}

############################################################
## 5. Simulation grid (mean/var scenarios)
############################################################

# Better-balanced default grid:
# - 3 K levels (small/medium/larger)
# - 4 density levels (sparse -> dense)
# - 3 hierarchy levels (weak/medium/strong)
# - both generative mechanisms (SST/WST)
K_grid <- c(3, 5, 8)
kappa_mean_grid <- c(0.75, 1.5, 3.0, 6.0)
kappa_cv <- 0.6

hierarchy_specs <- data.frame(
  hierch = c("weak", "medium", "strong"),
  psi_mean = c(0.20, 0.70, 1.30),
  psi_sd = c(0.30, 0.30, 0.30),
  stringsAsFactors = FALSE
)

gen_models <- c("SST", "WST")
# Keep the default paper reproduction lightweight enough for a fresh rerun.
n_rep      <- 3
n_items    <- 60

run_start_time <- Sys.time()
run_stamp <- format(run_start_time, "%Y%m%d_%H%M%S")
run_id <- paste0("DemoKvar_run_", run_stamp)
run_dir <- file.path("output", "simulation", "raw", "DemoKvar_runs", run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
message("Run directory: ", run_dir)

if (tolower(Sys.getenv("DEMOKVAR_SMOKE", unset = "false")) %in% c("1", "true", "yes")) {
  K_grid <- c(3)
  kappa_mean_grid <- c(1)
  hierarchy_specs <- data.frame(
    hierch = "weak",
    psi_mean = 0.25,
    psi_sd = 0.30,
    stringsAsFactors = FALSE
  )
  gen_models <- c("WST")
  n_rep <- 1
  n_items <- 24
}

n_rep_env <- suppressWarnings(as.integer(Sys.getenv("DEMOKVAR_N_REP", unset = "")))
if (!is.na(n_rep_env) && n_rep_env >= 1L) n_rep <- n_rep_env

n_items_env <- suppressWarnings(as.integer(Sys.getenv("DEMOKVAR_N_ITEMS", unset = "")))
if (!is.na(n_items_env) && n_items_env >= 4L) n_items <- n_items_env

demo_n_iter <- suppressWarnings(as.integer(Sys.getenv("DEMOKVAR_N_ITER", unset = "10000")))
demo_burn   <- suppressWarnings(as.integer(Sys.getenv("DEMOKVAR_BURN", unset = "1000")))
demo_thin   <- suppressWarnings(as.integer(Sys.getenv("DEMOKVAR_THIN", unset = "2")))

if (is.na(demo_n_iter) || demo_n_iter < 1L) demo_n_iter <- 10000L
if (is.na(demo_burn) || demo_burn < 0L || demo_burn >= demo_n_iter) demo_burn <- max(0L, floor(demo_n_iter / 5))
if (is.na(demo_thin) || demo_thin < 1L) demo_thin <- 2L

grid_base <- expand.grid(
  K_true     = K_grid,
  kappa_mean = kappa_mean_grid,
  gen_model  = gen_models,
  stringsAsFactors = FALSE
) %>%
  tidyr::crossing(hierarchy_specs) %>%
  arrange(gen_model, K_true, kappa_mean, factor(hierch, levels = c("weak", "medium", "strong"))) %>%
  mutate(
    psi_var   = psi_sd^2,
    kappa_var = (kappa_cv * kappa_mean)^2,
    rate_kap  = kappa_mean,         # keep old name for compatibility
    scenario_id = sprintf("scn_%03d", dplyr::row_number())
  )

grid_val <- grid_base %>%
  select(K_true, kappa_mean, hierch, gen_model, psi_mean, psi_var, kappa_var, rate_kap, scenario_id)

message(
  sprintf(
    "Grid: %d scenarios x %d reps = %d tasks",
    nrow(grid_val), n_rep, nrow(grid_val) * n_rep
  )
)

# Task grid: replicate x scenario
task_grid <- expand.grid(
  rep_id = seq_len(n_rep),
  idx    = seq_len(nrow(grid_val)),
  stringsAsFactors = FALSE
) %>%
  arrange(rep_id, idx)

out_file <- file.path(run_dir, paste0("full_simulation_crossfit_progress_", run_id, ".csv"))

############################################################
## 6. Parallel execution with progressive writing
############################################################

slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")))
slurm_ntasks <- suppressWarnings(as.integer(Sys.getenv("SLURM_NTASKS", unset = "")))
env_cores <- suppressWarnings(as.integer(Sys.getenv("DEMOKVAR_CORES", unset = "")))

if (!is.na(env_cores) && env_cores >= 1L) {
  n_cores <- env_cores
} else if (!is.na(slurm_cpus) && slurm_cpus >= 1L) {
  n_cores <- slurm_cpus
} else if (!is.na(slurm_ntasks) && slurm_ntasks >= 1L) {
  n_cores <- slurm_ntasks
} else {
  n_cores <- parallel::detectCores(logical = TRUE)
}

if (is.na(n_cores) || n_cores < 1L) {
  n_cores <- 1L
}
n_cores <- max(1L, min(n_cores, nrow(task_grid)))

message("Using ", n_cores, " cores.")

res_list <- mclapply(
  X = seq_len(nrow(task_grid)),
  FUN = function(tid) {
    rep_id <- task_grid$rep_id[tid]
    row_id <- task_grid$idx[tid]
    row    <- grid_val[row_id, ]
    
    seed <- 1000 + tid
    
    cat(sprintf(
      "Rep %d | Gen %s | K %d | kappa_mean %.2f | hierch %s | psi_mean %.2f\n",
      rep_id, row$gen_model, row$K_true, row$kappa_mean, row$hierch, row$psi_mean
    ))
    
    res_df <- tryCatch(
      run_replicate_crossfit(
        gen_model   = row$gen_model,
        K_true      = row$K_true,
        n           = n_items,
        kappa_mean  = row$kappa_mean,
        kappa_var   = row$kappa_var,
        psi_mean    = row$psi_mean,
        psi_var     = row$psi_var,
        hierch_label = row$hierch,
        seed        = seed,
        n_iter      = demo_n_iter,
        burn        = demo_burn,
        thin        = demo_thin
      ),
      error = function(e) {
        message(
          sprintf(
            "Error in tid=%d, rep=%d, gen=%s, K=%d, kappa_mean=%.2f, hierch=%s:\n  %s",
            tid, rep_id, row$gen_model, row$K_true,
            row$kappa_mean, row$hierch,
            conditionMessage(e)
          )
        )
        NULL
      }
    )
    
    if (!is.null(res_df)) {
      res_df$rep_id <- rep_id
      # redundant but explicit
      res_df$kappa_mean <- row$kappa_mean
      res_df$kappa_var  <- row$kappa_var
      res_df$psi_mean   <- row$psi_mean
      res_df$psi_var    <- row$psi_var
      res_df$hierch     <- row$hierch
      res_df$rate_kap   <- row$kappa_mean
      
      append_results_with_lock(res_df, out_file)
    }
    
    res_df
  },
  mc.cores      = n_cores,
  mc.preschedule = FALSE
)

final_table <- dplyr::bind_rows(res_list)

# Also write final combined table (in case you prefer a clean version)
final_file <- file.path(run_dir, paste0("full_simulation_crossfit_final_", run_id, ".csv"))
write.csv(final_table, final_file, row.names = FALSE)

git_hash <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1],
  error = function(e) NA_character_
)
if (!is.character(git_hash) || !length(git_hash) || is.na(git_hash) || !nzchar(git_hash)) {
  git_hash <- "NA"
}

run_end_time <- Sys.time()
manifest_lines <- c(
  paste0("run_id: ", run_id),
  paste0("script: DemoKvar.R"),
  paste0("sim_study_obsolete: true"),
  paste0("run_dir: ", run_dir),
  paste0("git_commit_short: ", git_hash),
  paste0("start_time: ", format(run_start_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("end_time: ", format(run_end_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("elapsed_seconds: ", round(as.numeric(difftime(run_end_time, run_start_time, units = "secs")), 3)),
  paste0("smoke_mode: ", tolower(Sys.getenv("DEMOKVAR_SMOKE", unset = "false"))),
  paste0("gen_models: ", paste(gen_models, collapse = ",")),
  paste0("K_grid: ", paste(K_grid, collapse = ",")),
  paste0("n_rep: ", n_rep),
  paste0("n_items: ", n_items),
  paste0("n_iter: ", demo_n_iter),
  paste0("burn: ", demo_burn),
  paste0("thin: ", demo_thin),
  paste0("K_init_ordered_models: n_items"),
  paste0("K_init_dcsbm: n_items"),
  paste0("K_expected_policy: K_true per replicate"),
  paste0("tau0_floor: 0.2"),
  paste0("b_kappa_ordered: 1"),
  paste0("b_eta_ordered: 1"),
  paste0("ordering_prior_mode: equivalence_class"),
  paste0("osbm_partition_prior: ", demo_osbm_partition_prior),
  paste0("theta_ocrp: ", demo_theta_ocrp),
  paste0("alpha_crp: ", demo_alpha_crp),
  paste0("progress_csv: ", out_file),
  paste0("final_csv: ", final_file),
  paste0("n_rows_final: ", nrow(final_table))
)
writeLines(manifest_lines, con = file.path(run_dir, "run_manifest.txt"))

############################################################
## 7. Post-processing (same spirit as before)
############################################################
# 
# # Clean types
# final_table <- final_table %>%
#   select(-starts_with("Unnamed")) %>%
#   mutate(
#     gen_model  = factor(gen_model,  levels = c("SST", "WST")),
#     fit_model  = factor(fit_model,  levels = c("SST", "WST", "DC-SBM")),
#     hierch     = factor(hierch,     levels = c("weak", "strong")),
#     K_true     = as.integer(K_true),
#     K_hat      = as.integer(K_hat),
#     rate_kap   = as.numeric(rate_kap),
#     kappa_mean = as.numeric(kappa_mean),
#     kappa_var  = as.numeric(kappa_var),
#     psi_mean   = as.numeric(psi_mean),
#     psi_var    = as.numeric(psi_var),
#     seed       = as.integer(seed)
#   )
# 
# # Main performance summary
# summary_main <- final_table %>%
#   group_by(gen_model, hierch, rate_kap, K_true, fit_model) %>%
#   summarise(
#     n_runs       = n(),
#     mean_K_hat   = mean(K_hat, na.rm = TRUE),
#     sd_K_hat     = sd(K_hat, na.rm = TRUE),
#     prop_K_true  = mean(K_hat == K_true, na.rm = TRUE),
#     prop_K_under = mean(K_hat < K_true, na.rm = TRUE),
#     prop_K_over  = mean(K_hat > K_true, na.rm = TRUE),
#     
#     mean_ari     = mean(ari, na.rm = TRUE),
#     sd_ari       = sd(ari, na.rm = TRUE),
#     mean_vi      = mean(vi, na.rm = TRUE),
#     
#     mean_elpd    = mean(elpd, na.rm = TRUE),
#     sd_elpd      = sd(elpd, na.rm = TRUE),
#     mean_looic   = mean(looic, na.rm = TRUE),
#     mean_pk_bad  = mean(pk_bad, na.rm = TRUE),
#     
#     .groups = "drop"
#   )
# 
# summary_main
# 
# # Diagnostics summary
# summary_diag <- final_table %>%
#   group_by(gen_model, hierch, rate_kap, K_true, fit_model) %>%
#   summarise(
#     mean_violation_rate = mean(violation_rate_mean, na.rm = TRUE),
#     sd_violation_rate   = sd(violation_rate_mean, na.rm = TRUE),
#     mean_viol_count     = mean(violation_count_mean, na.rm = TRUE),
#     mean_cross_mass     = mean(cross_mass_mean, na.rm = TRUE),
#     mean_cov_block_model = mean(coverage_block_model_avg, na.rm = TRUE),
#     mean_cov_block_emp   = mean(coverage_block_emp_avg,   na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# summary_diag
# 
# ############################################################
# ## 8. Example tables (same idea as your old ones)
# ############################################################
# 
# make_latex_table <- function(gen, K, hier) {
#   summary_main %>%
#     filter(gen_model == gen,
#            K_true    == K,
#            hierch    == hier) %>%
#     arrange(rate_kap, fit_model) %>%
#     select(rate_kap, fit_model,
#            mean_K_hat, prop_K_true,
#            mean_ari, mean_elpd, mean_pk_bad) %>%
#     kable(
#       format  = "latex",
#       booktabs = TRUE,
#       digits  = c(1, 1, 2, 2, 2, 1, 2),
#       caption = paste0(
#         gen, " generativo, $K^\\star=", K, "$, gerarchia ", hier,
#         ": recupero di $K$, ARI ed ELPD per modello di fit (prior mean/var)."
#       )
#     )
# }
# 
# # Examples:
# make_latex_table("SST", 3, "strong")
# make_latex_table("SST", 5, "strong")
# make_latex_table("WST", 3, "strong")
# make_latex_table("WST", 5, "strong")
# 
# tab_elpd_wide <- summary_main %>%
#   select(gen_model, hierch, K_true, rate_kap, fit_model, mean_elpd) %>%
#   mutate(
#     fit_model = as.character(fit_model),
#     rate_kap  = as.factor(rate_kap)
#   ) %>%
#   pivot_wider(
#     names_from  = fit_model,
#     values_from = mean_elpd
#   ) %>%
#   arrange(gen_model, hierch, K_true, rate_kap)
# 
# kable(
#   tab_elpd_wide,
#   format  = "latex",
#   booktabs = TRUE,
#   digits  = 1,
#   caption = "ELPD medio (prior mean/var) per combinazione di modello generativo, gerarchia, $K^\\star$ e densità (kappa_mean), a confronto tra modelli di fit."
# )
# 
# tab_violation <- summary_diag %>%
#   arrange(gen_model, hierch, K_true, rate_kap, fit_model) %>%
#   select(gen_model, hierch, K_true, rate_kap, fit_model,
#          mean_violation_rate, mean_viol_count,
#          mean_cross_mass, mean_cov_block_model)
# 
# kable(
#   tab_violation,
#   format  = "latex",
#   booktabs = TRUE,
#   digits  = c(0, 0, 0, 1, 0, 2, 1, 1, 2),
#   caption = "Statistiche di violazione dell'ordine e coverage (prior mean/var) per combinazione di scenario e modello di fit."
# )
# 
# ############################################################
# ## 9. Plots (same structure as before)
# ############################################################
# 
# ggplot(summary_main,
#        aes(x = rate_kap, y = mean_K_hat,
#            color = fit_model, group = interaction(fit_model, K_true))) +
#   geom_line() +
#   geom_point(aes(shape = factor(K_true)), size = 2) +
#   facet_grid(gen_model ~ hierch) +
#   labs(
#     x = "kappa_mean (densità blocchi)",
#     y = "K_hat medio",
#     color = "Modello di fit",
#     shape = "K_true",
#     title = "Recupero di K per scenario generativo (prior mean/var)"
#   ) +
#   theme_bw()
# 
# ggplot(summary_main,
#        aes(x = rate_kap, y = prop_K_true,
#            color = fit_model, group = interaction(fit_model, K_true))) +
#   geom_line() +
#   geom_point(aes(shape = factor(K_true)), size = 2) +
#   facet_grid(gen_model ~ hierch) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(
#     x = "kappa_mean",
#     y = "Pr(K_hat = K_true)",
#     color = "Modello di fit",
#     shape = "K_true",
#     title = "Probabilità di recupero esatto di K (prior mean/var)"
#   ) +
#   theme_bw()
# 
# ggplot(summary_main,
#        aes(x = rate_kap, y = mean_elpd,
#            color = fit_model, group = interaction(fit_model, K_true))) +
#   geom_line() +
#   geom_point(aes(shape = factor(K_true)), size = 2) +
#   facet_grid(gen_model ~ hierch) +
#   labs(
#     x = "kappa_mean",
#     y = "ELPD medio (loo)",
#     color = "Modello di fit",
#     shape = "K_true",
#     title = "Confronto predittivo tra modelli di fit (prior mean/var)"
#   ) +
#   theme_bw()
# 
# ggplot(summary_diag,
#        aes(x = rate_kap, y = mean_violation_rate,
#            color = fit_model, group = fit_model)) +
#   geom_line() +
#   geom_point(size = 2) +
#   facet_grid(gen_model ~ hierch) +
#   labs(
#     x = "kappa_mean",
#     y = "Violazione media (psi)",
#     color = "Modello di fit",
#     title = "Grado di violazione dell'ordine tra blocchi (prior mean/var)"
#   ) +
#   theme_bw()
