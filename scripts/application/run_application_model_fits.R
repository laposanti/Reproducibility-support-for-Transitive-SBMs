#!/usr/bin/env Rscript
# ======================================================================
# Application: CRP / ordered-CRP priors
#
# This driver extends the main application script to a larger dataset set.
# For the overlapping datasets, keep the same defaults as application.R so
# the model comparisons are directly comparable.
# ======================================================================
resolve_bundle_root <- function(relative_to_script = "../..") {
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

setwd(resolve_bundle_root())


run_start_time <- Sys.time()
run_stamp <- format(run_start_time, "%Y%m%d_%H%M%S")
run_id   <- paste0("application_run_", run_stamp)
run_dir  <- file.path("output", "application", "raw", run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Log to file + console ----
log_path <- file.path(run_dir, "run_log.txt")
log_con  <- file(log_path, open = "wt")
sink(log_con, type = "output",  split = TRUE)
sink(log_con, type = "message", append = TRUE)

cat("=== CRP / OCRP Application Run ===\n")
cat("Start:", format(run_start_time), "\n")
cat("Run directory:", run_dir, "\n\n")

suppressPackageStartupMessages({
  library(coda)
  library(mcclust)
  library(mcclust.ext)
  library(salso)
  library(dplyr)
  library(loo)
})

source("./helper_folder/models/ordered_sbm/shared_sampler_helpers.R")
source("./helper_folder/models/ordered_sbm/sst_helpers.R")
source("./helper_folder/models/ordered_sbm/wst_helpers.R")
source("./helper_folder/config/hyperparameter_setup.R")
source("./helper_folder/io/application_data_loader.R")
source("./helper_folder/simulation/simulation_study_helpers.R")
source("./helper_folder/diagnostics/transitivity_diagnostics.R")
source("./core/ppc_checks.R")
source("./core/transitive_sbm_sampler.R")
source("./core/DCSBM_varK.R")
if (!exists("log1pexp")) log1pexp <- function(x) ifelse(x < 35, log1p(exp(x)), x)

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

# ======================================================================
# Settings
# ======================================================================
C_KAPPA       <- 3
A_KAPPA       <- 1
THETA_OCRP    <- as.numeric(Sys.getenv("APP_THETA_OCRP", unset = "0.5"))
ALPHA_CRP     <- as.numeric(Sys.getenv("APP_ALPHA_CRP", unset = "0.5"))
N_ITER        <- as.integer(Sys.getenv("APP_N_ITER", unset = "10000"))
BURN          <- as.integer(Sys.getenv("APP_BURN", unset = "3000"))
THIN          <- as.integer(Sys.getenv("APP_THIN", unset = "2"))
SEED          <- as.integer(Sys.getenv("APP_SEED", unset = "42"))
MAX_LOO_CELLS <- as.numeric(Sys.getenv("APP_MAX_LOO_CELLS", unset = "2e8"))
OSBM_PARTITION_PRIOR <- .normalize_app_ordered_prior(
  Sys.getenv("APP_OSBM_PARTITION_PRIOR",
             unset = Sys.getenv("APP_PARTITION_PRIOR", unset = "OCRP"))
)

DATASET_CHOICES <- c(
   "moreno_sheep", "strauss_2019b", "mountain_goats", "citations_data", "macaques_data", "high_school"
)
datasets_env <- trimws(Sys.getenv("APP_DATASETS", unset = ""))
datasets <- if (nzchar(datasets_env)) strsplit(datasets_env, ",")[[1]] else DATASET_CHOICES
models <- c("WST", "SST", "DCSBM")

cat(sprintf("Config: %s theta=%.1f  |  CRP alpha=%.1f\n",
            OSBM_PARTITION_PRIOR, THETA_OCRP, ALPHA_CRP))
cat(sprintf("c_kappa=%.1f  a_kappa=%d  n_iter=%d  burn=%d  thin=%d  seed=%d\n",
            C_KAPPA, A_KAPPA, N_ITER, BURN, THIN, SEED))
cat(sprintf("max LOO log-lik cells: %.0f\n", MAX_LOO_CELLS))
cat(sprintf("DCSBM: CRP alpha=%.1f, baseline b_lambda=1  (citations: b_eta=11, b_lambda=7)\n",
            ALPHA_CRP))
cat(sprintf("Datasets: %s\n\n", paste(datasets, collapse = ", ")))

# ======================================================================
# Data loader
# ======================================================================
choose_dataset <- function(dataset = DATASET_CHOICES[1]) {
  load_application_adjacency(match.arg(dataset, choices = DATASET_CHOICES))
}

# ======================================================================
# Hyperparameters: OSBM aligned with application.R
# ======================================================================
get_aligned_osbm_hypers <- function(A, dataset, fit_model,
                                    theta_ocrp = 0.5,
                                    partition_prior = "OCRP",
                                    K_expected = NULL) {
  if (is.null(K_expected)) {
    K_expected <- max(2L, ceiling(sqrt(nrow(A))))
  }
  K_expected <- max(1L, min(as.integer(K_expected), nrow(A)))
  base <- get_principled_hypers_v2(
    A = A,
    K_expected = K_expected,
    c_kappa = C_KAPPA
  )
  base$tau0 <- max(base$tau0, 0.2)

  list(
    a_kappa = 1,
    b_kappa = 1,
    a_eta = 1,
    b_eta = 1,
    mu0 = base$mu0,
    sigma0 = base$sigma0,
    tau0 = base$tau0,
    theta_ocrp = theta_ocrp,
    partition_prior = partition_prior,
    meta = base$meta
  )
}

# ======================================================================
# DCSBM config (CRP prior)
# ======================================================================
build_dcsbm_config <- function(A, dataset, alpha_crp = 1.0) {
  b_eta    <- 1
  b_lambda <- 1
  list(
    a_eta = 1, b_eta = b_eta,
    a_lambda = 1, b_lambda = b_lambda,
    alpha_crp = alpha_crp,
    partition_prior = "CRP"
  )
}

# ======================================================================
# LOOIC
# ======================================================================
compute_looic <- function(A, out, fit_model, max_cells = Inf) {
  n <- nrow(A)
  DI <- build_dyad_index(n)
  n_saved <- if (is.matrix(out$z)) nrow(out$z) else length(out$z)
  n_dyads <- if (is.list(DI) && !is.null(DI$i)) {
    length(DI$i)
  } else if (is.matrix(DI)) {
    nrow(DI)
  } else {
    length(DI)
  }
  ll_cells <- as.double(n_saved) * as.double(n_dyads)
  if (is.finite(max_cells) && is.finite(ll_cells) && ll_cells > max_cells) {
    stop(sprintf(
      "Skipping LOOIC: %d saved draws x %d dyads = %.0f log-likelihood cells.",
      as.integer(n_saved), as.integer(n_dyads), ll_cells
    ), call. = FALSE)
  }
  LL <- if (fit_model %in% c("WST", "SST")) {
    loglik_matrix_modular(A, out, regime = fit_model, dyad_index = DI)
  } else {
    loglik_matrix_dcsbm(A = as.matrix(A), dcsbm_out = out, dyad_index = DI)
  }
  loo_fit <- loo::loo(LL)
  list(
    elpd_loo = loo_fit$estimates["elpd_loo", "Estimate"],
    p_loo    = loo_fit$estimates["p_loo",    "Estimate"],
    looic    = loo_fit$estimates["looic",    "Estimate"],
    pk_max   = max(loo_fit$diagnostics$pareto_k),
    pk_bad   = mean(loo_fit$diagnostics$pareto_k > 0.7)
  )
}

# ======================================================================
# K summary
# ======================================================================
extract_K_summary <- function(out) {
  z_chain <- if (is.matrix(out$z)) out$z else do.call(rbind, out$z)
  S <- nrow(z_chain)
  K_occ <- apply(z_chain, 1L, function(z) length(unique(z[!is.na(z)])))

  salso_cap <- as.integer(max(3L, min(ncol(z_chain),
    quantile(K_occ, 0.99, na.rm = TRUE, type = 1) + 1L)))
  z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L, nCores = 1L,
                        maxNClusters = salso_cap, maxZealousAttempts = 100000L)

  list(
    K_hat  = length(unique(z_hat)),
    K_mean = round(mean(K_occ), 2),
    K_lo   = quantile(K_occ, 0.025),
    K_hi   = quantile(K_occ, 0.975),
    z_hat  = z_hat
  )
}

# ======================================================================
# Fit one model
# ======================================================================
fit_one <- function(A, dataset, fit_model, hypers, n_iter, burn, thin, seed) {
  n <- nrow(A)
  t0 <- proc.time()
  out <- NULL; ok <- TRUE; msg <- NA_character_

  tryCatch({
    if (fit_model %in% c("WST", "SST")) {
      out <- modular_osbm_sampler(
        A = A, K = n,
        n_iter = n_iter, burn = burn, thin = thin,
        verbose = TRUE,
        psi_constraint = fit_model ,
        partition_prior = hypers$partition_prior,
        theta_ocrp = hypers$theta_ocrp,
        a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
        a_eta = hypers$a_eta, b_eta = hypers$b_eta,
        mu0 = hypers$mu0, sigma0 = hypers$sigma0, tau0 = hypers$tau0,
        eta_identifiability = "block_sum",
        use_mixing_moves = TRUE,
        seed = seed
      )
    } else {
      out <- fit_dcsbm_gibbs_gnedin(
        A = as.matrix(A),
        K_init = n,
        priors = list(
          a_eta = hypers$a_eta, b_eta = hypers$b_eta,
          a_lambda = hypers$a_lambda, b_lambda = hypers$b_lambda,
          alpha_crp = hypers$alpha_crp,
          partition_prior = hypers$partition_prior
        ),
        iters = n_iter, burn_in = burn, thin = thin,
        verbose = 50, seed = seed
      )
    }
  }, error = function(e) { ok <<- FALSE; msg <<- conditionMessage(e) })
  elapsed <- (proc.time() - t0)[["elapsed"]]

  if (!ok || is.null(out)) {
    cat(sprintf("    *** FAILED: %s ***\n", msg))
    return(data.frame(
      dataset = dataset, fit_model = fit_model,
      K_hat = NA, K_mean = NA, K_lo = NA, K_hi = NA,
      elpd_loo = NA, p_loo = NA, looic = NA,
      pk_max = NA, pk_bad = NA,
      b_kappa = if (is.null(hypers$b_kappa)) NA else hypers$b_kappa,
      time_sec = round(elapsed, 1), note = msg,
      stringsAsFactors = FALSE
    ))
  }

  # K summary
  ks <- tryCatch(extract_K_summary(out), error = function(e) {
    list(K_hat = NA, K_mean = NA, K_lo = NA, K_hi = NA, z_hat = NULL)
  })

  # LOOIC
  loo_note <- NA_character_
  loo_res <- tryCatch(compute_looic(A, out, fit_model, max_cells = MAX_LOO_CELLS), error = function(e) {
    loo_note <<- conditionMessage(e)
    cat(sprintf("    LOOIC error: %s\n", e$message))
    list(elpd_loo = NA, p_loo = NA, looic = NA, pk_max = NA, pk_bad = NA)
  })

  # Save fit
  fit_path <- file.path(run_dir, sprintf("%s_%s_fit.rds", dataset, fit_model))
  saveRDS(out, fit_path)

  b_val <- if (is.null(hypers$b_kappa)) NA_real_ else hypers$b_kappa

  data.frame(
    dataset = dataset, fit_model = fit_model,
    K_hat = ks$K_hat, K_mean = ks$K_mean, K_lo = ks$K_lo, K_hi = ks$K_hi,
    elpd_loo = round(loo_res$elpd_loo, 1),
    p_loo    = round(loo_res$p_loo, 1),
    looic    = round(loo_res$looic, 1),
    pk_max   = round(loo_res$pk_max, 3),
    pk_bad   = round(loo_res$pk_bad, 3),
    b_kappa  = round(b_val, 2),
    time_sec = round(elapsed, 1),
    note     = loo_note,
    stringsAsFactors = FALSE
  )
}

# ======================================================================
# Main loop
# ======================================================================
all_results <- list()

for (ds in datasets) {
  cat("\n================================================================\n")
  cat(sprintf("Dataset: %s\n", ds))
  cat("================================================================\n")

  A <- tryCatch(choose_dataset(ds), error = function(e) {
    cat(sprintf("  *** Cannot load %s: %s ***\n", ds, e$message)); NULL
  })
  if (is.null(A)) next

  n <- nrow(A)
  mean_deg <- sum(A) / n
  K_sp <- estimate_K_spectral(A)
  density_val <- mean_deg / (n - 1)
  deg <- rowSums(A) + colSums(A)
  cv_deg <- round(sd(deg) / mean(deg), 3)
  cat(sprintf("  n=%d  sum(A)=%d  mean_deg=%.1f  density=%.3f  CV(deg)=%.3f  K_sp=%d\n\n",
              n, sum(A), mean_deg, density_val, cv_deg, K_sp))

  ds_results <- list()

  for (regime in c("WST", "SST")) {
    hypers <- get_aligned_osbm_hypers(
      A,
      dataset = ds,
      fit_model = regime,
      theta_ocrp = THETA_OCRP,
      partition_prior = OSBM_PARTITION_PRIOR,
      K_expected = K_sp
    )
    cat(sprintf("  --- %s ---\n", regime))
    cat(sprintf(
      "    %s theta=%.1f  a_kappa=%.1f  b_kappa=%.2f  a_eta=%.1f  b_eta=%.1f\n",
      hypers$partition_prior, hypers$theta_ocrp,
      hypers$a_kappa, hypers$b_kappa, hypers$a_eta, hypers$b_eta
    ))
    res <- fit_one(A, ds, regime, hypers, N_ITER, BURN, THIN, SEED)
    cat(sprintf("    => K_hat=%s  K_mean=%s  LOOIC=%s  p_loo=%s  pk>0.7=%s  (%ss)\n",
                res$K_hat, res$K_mean, res$looic, res$p_loo,
                round(100*res$pk_bad, 1), res$time_sec))
    ds_results[[length(ds_results) + 1]] <- res
  }

  # --- DCSBM with CRP ---
  cat("  --- DCSBM ---\n")
  dcsbm_h <- build_dcsbm_config(A, ds, alpha_crp = ALPHA_CRP)
  cat(sprintf("    CRP alpha=%.1f  b_eta=%d  b_lambda=%d\n",
              dcsbm_h$alpha_crp, dcsbm_h$b_eta, dcsbm_h$b_lambda))
  res <- fit_one(A, ds, "DCSBM", dcsbm_h, N_ITER, BURN, THIN, SEED)
  cat(sprintf("    => K_hat=%s  K_mean=%s  LOOIC=%s  p_loo=%s  pk>0.7=%s  (%ss)\n",
              res$K_hat, res$K_mean, res$looic, res$p_loo,
              round(100*res$pk_bad, 1), res$time_sec))
  ds_results[[length(ds_results) + 1]] <- res

  ds_df <- do.call(rbind, ds_results)
  all_results[[ds]] <- ds_df
  write.csv(ds_df, file.path(run_dir, sprintf("%s_results.csv", ds)), row.names = FALSE)
  saveRDS(ds_df, file.path(run_dir, sprintf("%s_results.rds", ds)))
  cat(sprintf("  Saved: %s_results.csv and %s_results.rds\n", ds, ds))
}

# ======================================================================
# Combined results + leaderboard
# ======================================================================
final <- if (length(all_results) > 0L) do.call(rbind, all_results) else data.frame()
if (nrow(final) > 0L) rownames(final) <- NULL
summary_csv <- file.path(run_dir, sprintf("applications_results_summary_%s.csv", run_id))
summary_rds <- file.path(run_dir, sprintf("applications_results_summary_%s.rds", run_id))
write.csv(final, summary_csv, row.names = FALSE)
saveRDS(final, summary_rds)

# Backward-compatible filenames used by some downstream scripts.
write.csv(final, file.path(run_dir, "all_results.csv"), row.names = FALSE)
saveRDS(final, file.path(run_dir, "all_results.rds"))

cat("\n\n================================================================\n")
cat("FULL RESULTS\n")
cat("================================================================\n\n")
if (nrow(final) > 0L) {
  print(final[, c("dataset", "fit_model", "K_hat", "K_mean", "p_loo", "looic",
                  "pk_bad", "b_kappa", "time_sec")], row.names = FALSE)
} else {
  cat("No successful dataset fits; results tables are empty.\n")
}

cat("\n\n================================================================\n")
cat("LEADERBOARD (best LOOIC per dataset)\n")
cat("================================================================\n\n")
if (nrow(final) > 0L) {
  for (ds in datasets) {
    sub <- final[final$dataset == ds & !is.na(final$looic), ]
    if (nrow(sub) == 0) { cat(ds, ": no valid results\n"); next }
    best <- sub[which.min(sub$looic), ]
    cat(sprintf("  %s: %s  LOOIC=%.1f  K=%s  p_loo=%.1f\n",
                ds, best$fit_model, best$looic, best$K_hat, best$p_loo))
  }
} else {
  cat("No leaderboard available because final results are empty.\n")
}

# ======================================================================
# Pairwise LOO comparison (delta elpd with SE via loo::loo_compare)
# ======================================================================
if (nrow(final) > 0L) {
  cat("\n\n================================================================\n")
  cat("MODEL COMPARISON (best vs second-best, loo::loo_compare)\n")
  cat("================================================================\n\n")

  comp_rows <- list()

  for (ds in datasets) {
    A_ds <- tryCatch(choose_dataset(ds), error = function(e) NULL)
    if (is.null(A_ds)) { cat(sprintf("  %s: cannot reload data\n", ds)); next }

    DI_ds <- build_dyad_index(nrow(A_ds))
    n_obs <- if (is.list(DI_ds) && !is.null(DI_ds$i)) {
      length(DI_ds$i)
    } else if (is.matrix(DI_ds)) {
      nrow(DI_ds)
    } else {
      length(DI_ds)
    }

    loo_list <- list()
    for (mod in models) {
      fp <- file.path(run_dir, sprintf("%s_%s_fit.rds", ds, mod))
      if (!file.exists(fp)) next
      fit <- readRDS(fp)
      loo_obj <- tryCatch({
        LL <- if (mod %in% c("WST", "SST")) {
          loglik_matrix_modular(A_ds, fit, regime = mod, dyad_index = DI_ds)
        } else {
          loglik_matrix_dcsbm(A = as.matrix(A_ds), dcsbm_out = fit, dyad_index = DI_ds)
        }
        loo::loo(LL)
      }, error = function(e) NULL)
      if (!is.null(loo_obj)) loo_list[[mod]] <- loo_obj
      rm(fit); gc(verbose = FALSE)
    }

    if (length(loo_list) < 2L) next

    comp <- loo::loo_compare(loo_list)
    best_mod   <- rownames(comp)[1]
    second_mod <- rownames(comp)[2]
    d_elpd     <- comp[2, "elpd_diff"]
    se_d       <- comp[2, "se_diff"]
    d_looic    <- -2 * d_elpd
    se_looic   <- 2 * se_d
    z_ratio    <- abs(d_elpd / se_d)
    d_per_obs  <- d_elpd / n_obs
    verdict    <- if (z_ratio > 4) "clear winner" else if (z_ratio > 2) "weak edge" else "indistinguishable"

    third_mod    <- if (nrow(comp) >= 3) rownames(comp)[3] else NA_character_
    d_elpd_3     <- if (nrow(comp) >= 3) comp[3, "elpd_diff"] else NA_real_
    se_d_3       <- if (nrow(comp) >= 3) comp[3, "se_diff"]   else NA_real_

    cat(sprintf("  %s:\n", ds))
    cat(sprintf("    1st: %-6s   2nd: %-6s   delta_elpd = %.1f (SE = %.1f)  |z| = %.1f  =>  %s\n",
                best_mod, second_mod, d_elpd, se_d, z_ratio, verdict))
    cat(sprintf("    delta_LOOIC = %.1f (SE = %.1f)   delta_elpd/obs = %.4f   n_obs = %d\n",
                d_looic, se_looic, d_per_obs, n_obs))
    if (!is.na(third_mod)) {
      cat(sprintf("    3rd: %-6s   delta_elpd = %.1f (SE = %.1f)\n",
                  third_mod, d_elpd_3, se_d_3))
    }
    cat("\n")

    comp_rows[[ds]] <- data.frame(
      dataset = ds, best_model = best_mod, second_model = second_mod,
      delta_elpd = round(d_elpd, 2), se_delta_elpd = round(se_d, 2),
      delta_looic = round(d_looic, 2), se_delta_looic = round(se_looic, 2),
      z_ratio = round(z_ratio, 2), delta_elpd_per_obs = round(d_per_obs, 4),
      n_obs = n_obs, verdict = verdict,
      third_model = third_mod,
      delta_elpd_3rd = round(d_elpd_3, 2), se_delta_3rd = round(se_d_3, 2),
      stringsAsFactors = FALSE
    )
  }

  if (length(comp_rows) > 0L) {
    comp_table <- do.call(rbind, comp_rows)
    rownames(comp_table) <- NULL
    comp_csv <- file.path(run_dir, "model_comparison_loo.csv")
    write.csv(comp_table, comp_csv, row.names = FALSE)
    saveRDS(comp_table, file.path(run_dir, "model_comparison_loo.rds"))
    cat(sprintf("Saved: %s\n", comp_csv))
  }
}

# Manifest
run_end_time <- Sys.time()
git_hash <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1],
  error = function(e) NA_character_
)
if (!is.character(git_hash) || !length(git_hash) || is.na(git_hash) || !nzchar(git_hash)) {
  git_hash <- "NA"
}

manifest_lines <- c(
  paste0("run_id: ", run_id),
  paste0("script: scripts/application/application_crp_ocrp.R"),
  paste0("run_dir: ", run_dir),
  paste0("git_commit_short: ", git_hash),
  paste0("start_time: ", format(run_start_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("end_time: ", format(run_end_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("elapsed_seconds: ", round(as.numeric(difftime(run_end_time, run_start_time, units = "secs")), 3)),
  paste0("datasets: ", paste(datasets, collapse = ",")),
  paste0("fit_models: ", paste(models, collapse = ",")),
  paste0("hyperparams.WST: partition_prior=", tolower(OSBM_PARTITION_PRIOR), ",theta_ocrp=", THETA_OCRP,
         ",a_kappa=1,b_kappa=1,a_eta=1,b_eta=1,mu0/sigma0/tau0=get_principled_hypers_v2(K_expected=K_sp,c_kappa=3),tau0_floor=0.2,K_init=n_items"),
  paste0("hyperparams.SST: partition_prior=", tolower(OSBM_PARTITION_PRIOR), ",theta_ocrp=", THETA_OCRP,
         ",a_kappa=1,b_kappa=1,a_eta=1,b_eta=1,mu0/sigma0/tau0=get_principled_hypers_v2(K_expected=K_sp,c_kappa=3),tau0_floor=0.2,K_init=n_items"),
  paste0("hyperparams.DCSBM: partition_prior=CRP,alpha_crp=", ALPHA_CRP,
         ",a_eta=1,b_eta=1,a_lambda=1,b_lambda=1"),
  paste0("max_loo_cells: ", format(MAX_LOO_CELLS, scientific = FALSE)),
  paste0("n_iter: ", N_ITER),
  paste0("burn: ", BURN),
  paste0("thin: ", THIN),
  paste0("seed: ", SEED),
  paste0("summary_csv: ", summary_csv),
  paste0("summary_rds: ", summary_rds),
  paste0("n_rows_summary: ", nrow(final)),
  paste0("n_error_notes: ", if ("note" %in% names(final)) sum(!is.na(final$note) & nzchar(final$note)) else 0),
  paste0("R_version: ", R.version.string),
  paste0("platform: ", R.version$platform),
  paste0("packages: ", paste(sapply(c("coda","mcclust","mcclust.ext","salso","loo","dplyr"), function(p) paste0(p, "=", packageVersion(p))), collapse = ", "))
)
writeLines(manifest_lines, con = file.path(run_dir, "run_manifest.txt"))

sink(type = "message")
sink(type = "output")
close(log_con)
