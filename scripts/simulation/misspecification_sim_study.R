#!/usr/bin/env Rscript
# =============================================================================
# Misspecification Simulation Study
#
# PURPOSE:
#   Generate networks from a NON-HIERARCHICAL model (psi = 0, rho = 0.5 for
#   all off-diagonal pairs) and fit WST-OSBM, SST-OSBM, and DC-SBM.
#   Diagnostics reveal whether the models correctly detect the absence of
#   ordering structure.
#
# METRICS:
#   K_hat, ARI, NVI, ELPD_loo, LOOIC, violation_rate, EBF (WST/SST),
#   thetaW/thetaS transitivity rates
#
# USAGE:
#   Rscript misspecification_sim_study.R                     # full run
#   SMOKE=1 Rscript misspecification_sim_study.R             # quick check
# =============================================================================

## ---- Environment knobs ---------------------------------------------------
SMOKE       <- as.integer(Sys.getenv("SMOKE", unset = "0")) == 1L
N_REP       <- as.integer(Sys.getenv("N_REP",   unset = if (SMOKE) "1" else "5"))
N_ITER      <- as.integer(Sys.getenv("N_ITER",  unset = if (SMOKE) "600" else "8000"))
BURN        <- as.integer(Sys.getenv("BURN",    unset = if (SMOKE) "200" else "2000"))
THIN        <- as.integer(Sys.getenv("THIN",    unset = if (SMOKE) "1"  else "2"))
N_SIM       <- as.integer(Sys.getenv("N_SIM",   unset = if (SMOKE) "60" else "60"))
BASE_SEED   <- as.integer(Sys.getenv("BASE_SEED", unset = "7777"))

## ---- Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(BayesLogit)
  library(truncnorm)
  library(salso)
  library(fossil)
  library(mcclust)
  library(coda)
  library(loo)
  library(dplyr)
  library(tidyr)
})

## ---- Source project code -------------------------------------------------
local_dir <- "/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/"
if (dir.exists(local_dir)) {
  setwd(local_dir)
} else {
  # On the cluster, SLURM sets SLURM_SUBMIT_DIR to the submission directory
  submit_dir <- Sys.getenv("SLURM_SUBMIT_DIR", unset = "")
  if (nzchar(submit_dir)) setwd(submit_dir)
}
source("helper_folder/sim_study_helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("helper_folder/transitivity_check_helper.R")
source("helper_folder/helper.R")
source("core/my_best_try_so_far.R")
source("core/DCSBM_varK.R")

## ---- Output setup --------------------------------------------------------
out_dir <- file.path("output", "simulation", "raw", paste0("misspec_sim_",
                     format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
log_con  <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)

cat(sprintf("=== Misspecification simulation study started: %s ===\n", Sys.time()))
cat(sprintf("SMOKE=%s  N_REP=%d  N_ITER=%d  BURN=%d  THIN=%d  N_SIM=%d  SEED=%d\n",
            SMOKE, N_REP, N_ITER, BURN, THIN, N_SIM, BASE_SEED))

sgn      <- function(x) ifelse(x > 0, 1L, ifelse(x < 0, -1L, 0L))
logistic <- function(x) 1 / (1 + exp(-x))

# =============================================================================
# Reuse helper functions from extended_sim_study_v2.R
# =============================================================================

setup_hypers <- function(A, method, gen_model, K_true) {
  n <- nrow(A)
  a_kappa <- 1;  b_kappa <- 1
  a_eta   <- 1;  b_eta   <- 1
  mu0     <- 0;  sigma0  <- 2
  tau0    <- 0.15
  gamma_gn <- 0.5
  sample_b   <- FALSE
  alpha0_bk  <- 1.0
  beta0_bk   <- 0.01
  use_moves  <- TRUE

  if (method == "spectral_ruleA") {
    K_spec   <- estimate_K_spectral(A)
    b_kappa  <- 5 * n / max(K_spec, 2)
    gamma_gn <- choose_gamma_from_K_expected(n, K_expected = max(K_spec, 2),
                  ordering_prior_mode = "equivalence_class")

  } else if (method == "calibrate_full") {
    cal <- calibrate_osbm_hypers(A, regime = gen_model)
    a_kappa  <- cal$a_kappa;  b_kappa <- cal$b_kappa
    a_eta    <- cal$a_eta;    b_eta   <- cal$b_eta
    mu0      <- cal$mu0;      sigma0  <- cal$sigma0
    tau0     <- cal$tau0
    gamma_gn <- cal$gamma

  } else {
    stop("Unknown method: ", method)
  }

  list(
    a_kappa = a_kappa, b_kappa = b_kappa,
    a_eta = a_eta, b_eta = b_eta,
    mu0 = mu0, sigma0 = sigma0, tau0 = tau0,
    gamma_gn = gamma_gn,
    sample_b_kappa = sample_b,
    alpha0_bkappa  = alpha0_bk,
    beta0_bkappa   = beta0_bk,
    use_mixing_moves = use_moves
  )
}

# =============================================================================
# OSBM fit + diagnostics (K, ARI, VI, LOOIC, transitivity, EBF, violation)
# =============================================================================
run_osbm_fit_diag <- function(A, K_true, z_true, fit_model, hypers, seed,
                              n_iter = N_ITER, burn = BURN, thin_ = THIN) {
  n <- nrow(A)
  K_init <- max(K_true + 3L, 6L)

  t0 <- proc.time()[3]
  out <- tryCatch(
    modular_osbm_sampler(
      A = A, K = K_init,
      n_iter = n_iter, burn = burn, thin = thin_,
      psi_constraint = fit_model ,
      partition_prior = "GN",
      gamma_gn       = hypers$gamma_gn,
      a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
      a_eta   = hypers$a_eta,   b_eta   = hypers$b_eta,
      mu0     = hypers$mu0,     sigma0  = hypers$sigma0,
      tau0    = hypers$tau0,
      use_mixing_moves = hypers$use_mixing_moves,
      sample_b_kappa   = hypers$sample_b_kappa,
      alpha0_bkappa    = hypers$alpha0_bkappa,
      beta0_bkappa     = hypers$beta0_bkappa,
      seed    = seed,
      verbose = FALSE
    ),
    error = function(e) { message("  OSBM ERROR: ", e$message); NULL }
  )
  elapsed <- proc.time()[3] - t0

  empty_row <- data.frame(
    K_hat = NA, K_post_mean = NA, K_post_mode = NA,
    ari = NA, vi = NA, elpd_loo = NA, looic = NA,
    ess_K = NA, bk_final = NA, time_sec = elapsed,
    stringsAsFactors = FALSE
  )

  if (is.null(out)) {
    diag_vec <- .empty_diag_vec()
    return(cbind(empty_row, as.data.frame(t(diag_vec))))
  }

  # K trace and z chain
  z_chain <- do.call(rbind, lapply(out$z, as.integer))
  K_trace <- out$K_trace
  S <- length(K_trace)

  # Partition estimate
  K_occ <- apply(z_chain, 1, function(z) length(unique(z)))
  salso_cap <- as.integer(max(3L, min(ncol(z_chain),
                  quantile(K_occ, 0.99, type = 1) + 1L)))
  z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L,
                         maxNClusters = salso_cap)
  K_hat <- length(unique(z_hat))
  K_post_mean <- mean(K_trace)
  K_post_mode <- as.integer(names(which.max(table(K_trace))))

  ari <- if (!is.null(z_true)) fossil::adj.rand.index(z_hat, z_true) else NA_real_
  vi  <- if (!is.null(z_true)) mcclust::vi.dist(z_hat, z_true) else NA_real_

  # LOO-IC
  DI <- build_dyad_index(n)
  LL <- tryCatch(
    loglik_matrix_modular(A, out, regime = fit_model, dyad_index = DI),
    error = function(e) NULL
  )
  elpd_loo <- looic_val <- NA_real_
  if (!is.null(LL)) {
    loo_fit <- tryCatch(loo::loo(LL), error = function(e) NULL)
    if (!is.null(loo_fit)) {
      elpd_loo  <- loo_fit$estimates["elpd_loo", "Estimate"]
      looic_val <- loo_fit$estimates["looic", "Estimate"]
    }
  }

  ess_K <- tryCatch(
    as.numeric(coda::effectiveSize(coda::mcmc(K_trace))),
    error = function(e) NA_real_
  )
  bk_final <- hypers$b_kappa

  base_row <- data.frame(
    K_hat = K_hat, K_post_mean = round(K_post_mean, 2),
    K_post_mode = K_post_mode,
    ari = round(ari, 4), vi = round(vi, 4),
    elpd_loo = round(elpd_loo, 2), looic = round(looic_val, 2),
    ess_K = round(ess_K, 1), bk_final = round(bk_final, 2),
    time_sec = round(elapsed, 1),
    stringsAsFactors = FALSE
  )

  # Transitivity / EBF / violation diagnostics
  K_max_hint <- max(K_trace, na.rm = TRUE)
  diag_vec <- tryCatch(
    summarise_osbm_diagnostics(
      out_relab  = out,
      regime     = fit_model,
      K_max_hint = K_max_hint,
      z_hat      = z_hat,
      n          = n,
      A          = A,
      m_items    = 500L,
      T_block    = 500L
    ),
    error = function(e) {
      message("  DIAG ERROR (OSBM): ", e$message)
      .empty_diag_vec()
    }
  )

  cbind(base_row, as.data.frame(t(diag_vec)))
}

# =============================================================================
# DCSBM fit + diagnostics
# =============================================================================
run_dcsbm_fit_diag <- function(A, K_true, z_true, gamma_gn, seed,
                               n_iter = N_ITER, burn = BURN, thin_ = THIN) {
  n <- nrow(A)
  t0 <- proc.time()[3]
  out <- tryCatch(
    fit_dcsbm_gibbs_gnedin(
      A = as.matrix(A), K_init = 15L,
      priors = list(a_eta = 1, b_eta = 1, a_lambda = 1, b_lambda = 1,
                    gamma_gnedin = gamma_gn),
      iters = n_iter, burn_in = burn, thin = thin_,
      verbose = 0, seed = seed
    ),
    error = function(e) { message("  DCSBM ERROR: ", e$message); NULL }
  )
  elapsed <- proc.time()[3] - t0

  empty_row <- data.frame(
    K_hat = NA, K_post_mean = NA, K_post_mode = NA,
    ari = NA, vi = NA, elpd_loo = NA, looic = NA,
    ess_K = NA, bk_final = NA, time_sec = elapsed,
    stringsAsFactors = FALSE
  )

  if (is.null(out)) {
    diag_vec <- .empty_diag_vec()
    return(cbind(empty_row, as.data.frame(t(diag_vec))))
  }

  # z chain
  if (is.matrix(out$z)) {
    z_chain <- out$z
  } else {
    z_chain <- do.call(rbind, lapply(out$z, as.integer))
  }
  K_occ <- apply(z_chain, 1, function(z) length(unique(z)))
  K_trace <- K_occ
  S <- nrow(z_chain)

  salso_cap <- as.integer(max(3L, min(ncol(z_chain),
                  quantile(K_occ, 0.99, type = 1) + 1L)))
  z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L,
                         maxNClusters = salso_cap)
  K_hat <- length(unique(z_hat))
  K_post_mean <- mean(K_trace)
  K_post_mode <- as.integer(names(which.max(table(K_trace))))

  ari <- if (!is.null(z_true)) fossil::adj.rand.index(z_hat, z_true) else NA_real_
  vi  <- if (!is.null(z_true)) mcclust::vi.dist(z_hat, z_true) else NA_real_

  DI <- build_dyad_index(n)
  LL <- tryCatch(
    loglik_matrix_dcsbm(A = as.matrix(A), dcsbm_out = out, dyad_index = DI),
    error = function(e) NULL
  )
  elpd_loo <- looic_val <- NA_real_
  if (!is.null(LL)) {
    loo_fit <- tryCatch(loo::loo(LL), error = function(e) NULL)
    if (!is.null(loo_fit)) {
      elpd_loo  <- loo_fit$estimates["elpd_loo", "Estimate"]
      looic_val <- loo_fit$estimates["looic", "Estimate"]
    }
  }

  ess_K <- tryCatch(
    as.numeric(coda::effectiveSize(coda::mcmc(K_trace))),
    error = function(e) NA_real_
  )

  base_row <- data.frame(
    K_hat = K_hat, K_post_mean = round(K_post_mean, 2),
    K_post_mode = K_post_mode,
    ari = round(ari, 4), vi = round(vi, 4),
    elpd_loo = round(elpd_loo, 2), looic = round(looic_val, 2),
    ess_K = round(ess_K, 1), bk_final = NA_real_,
    time_sec = round(elapsed, 1),
    stringsAsFactors = FALSE
  )

  # Transitivity / EBF / violation diagnostics
  diag_vec <- tryCatch(
    summarise_dcsbm_diagnostics(
      fit     = out,
      z_hat   = z_hat,
      K       = K_hat,
      n       = n,
      A       = A,
      m_items = 500L,
      T_block = 500L
    ),
    error = function(e) {
      message("  DIAG ERROR (DCSBM): ", e$message)
      .empty_diag_vec()
    }
  )

  cbind(base_row, as.data.frame(t(diag_vec)))
}


# #########################################################################
# SIMULATION
# #########################################################################

cat("\n", strrep("=", 70), "\n")
cat("MISSPECIFICATION STUDY: psi=0  (no hierarchy)\n")
cat(strrep("=", 70), "\n")

## ---- Scenario grid -------------------------------------------------------
K_vals     <- if (SMOKE) c(3) else c(3, 5, 7)
osbm_methods <- c("spectral_ruleA", "calibrate_full")

sim_grid <- expand.grid(
  K_true = K_vals,
  rep    = seq_len(N_REP),
  stringsAsFactors = FALSE
)

cat(sprintf("Scenario grid: K ∈ {%s} × %d reps = %d datasets\n",
            paste(K_vals, collapse = ","), N_REP, nrow(sim_grid)))
cat(sprintf("Fits per dataset: %d OSBM (WST+SST × %d methods) + 1 DCSBM\n",
            2 * length(osbm_methods), length(osbm_methods)))

set.seed(BASE_SEED)
sim_seeds <- sample.int(1e6, nrow(sim_grid))

all_rows <- list()

for (i in seq_len(nrow(sim_grid))) {
  Kt  <- sim_grid$K_true[i]
  rep <- sim_grid$rep[i]
  seed_data <- sim_seeds[i]

  cat(sprintf("\n--- [%d/%d] K=%d rep=%d seed=%d ---\n",
              i, nrow(sim_grid), Kt, rep, seed_data))

  ## ---- Generate NON-HIERARCHICAL data (psi = 0) --------------------------
  set.seed(seed_data)
  z_true  <- rep(1:Kt, length.out = N_SIM)
  kappa_t <- make_kappa(Kt, mean = 2.0, var = 1.0)
  eta_t   <- runif(N_SIM, 0.8, 1.2)

  # psi = 0 => rho = 0.5 everywhere (no dominance)
  psi_zero_wst <- matrix(0, Kt, Kt)
  A <- as.matrix(simulate_osbm(N_SIM, Kt, z_true, eta_t, kappa_t,
                                psi_zero_wst, "WST"))

  cat(sprintf("  sum(A)=%d  density=%.3f  [psi=0, symmetric directions]\n",
              sum(A), sum(A > 0) / (N_SIM * (N_SIM - 1))))

  ## ---- OSBM fits: WST and SST with each method --------------------------
  fit_models_osbm <- c("WST", "SST")

  for (method in osbm_methods) {
    for (fm in fit_models_osbm) {
      hypers <- setup_hypers(A, method, fm, Kt)
      seed_fit <- seed_data + match(method, osbm_methods) * 100L +
                  match(fm, fit_models_osbm)

      cat(sprintf("  %s / fit=%s  b_kappa=%.1f gamma=%.3f ... ",
                  method, fm, hypers$b_kappa, hypers$gamma_gn))
      flush.console()

      res <- run_osbm_fit_diag(A, Kt, z_true, fm, hypers, seed_fit)

      cat(sprintf("K=%s ARI=%.3f LOOIC=%s viol=%.3f (%.1fs)\n",
                  res$K_hat, res$ari,
                  ifelse(is.na(res$looic), "NA", sprintf("%.0f", res$looic)),
                  ifelse(is.na(res$violation_rate_mean), NA,
                         res$violation_rate_mean),
                  res$time_sec))

      row <- cbind(
        data.frame(
          gen_model = "NONE", fit_model = fm, K_true = Kt,
          strength = "none", rep = rep,
          method = method, seed = seed_data,
          stringsAsFactors = FALSE
        ),
        res
      )
      all_rows[[length(all_rows) + 1L]] <- row

      # Checkpoint
      write.csv(dplyr::bind_rows(all_rows),
                file.path(out_dir, "misspec_results_streaming.csv"),
                row.names = FALSE)
    }
  }

  ## ---- DCSBM fit ---------------------------------------------------------
  gamma_dc <- choose_gamma_from_K_expected(N_SIM, K_expected = Kt,
                ordering_prior_mode = "equivalence_class")
  seed_dc <- seed_data + 999L

  cat(sprintf("  DCSBM  gamma=%.3f ... ", gamma_dc))
  flush.console()

  res_dc <- run_dcsbm_fit_diag(A, Kt, z_true, gamma_dc, seed_dc)

  cat(sprintf("K=%s ARI=%.3f LOOIC=%s (%.1fs)\n",
              res_dc$K_hat, res_dc$ari,
              ifelse(is.na(res_dc$looic), "NA", sprintf("%.0f", res_dc$looic)),
              res_dc$time_sec))

  row_dc <- cbind(
    data.frame(
      gen_model = "NONE", fit_model = "DCSBM", K_true = Kt,
      strength = "none", rep = rep,
      method = "dcsbm_varK", seed = seed_data,
      stringsAsFactors = FALSE
    ),
    res_dc
  )
  all_rows[[length(all_rows) + 1L]] <- row_dc

  write.csv(dplyr::bind_rows(all_rows),
            file.path(out_dir, "misspec_results_streaming.csv"),
            row.names = FALSE)
}

## ---- Final save and summary ----------------------------------------------
results_df <- dplyr::bind_rows(all_rows)
write.csv(results_df, file.path(out_dir, "misspec_results.csv"),
          row.names = FALSE)

cat(sprintf("\n%s\nResults: %d rows saved to %s\n",
            strrep("=", 70), nrow(results_df), out_dir))

## ---- Quick summary table -------------------------------------------------
key_cols <- c("fit_model", "method", "K_true", "K_hat", "ari", "vi",
              "elpd_loo", "looic",
              "violation_rate_mean", "bf_wst_0", "bf_sst_0",
              "thetaW_block_model_mean", "thetaS_block_model_mean")
avail <- intersect(key_cols, names(results_df))
summary_df <- results_df %>%
  group_by(fit_model, method, K_true) %>%
  summarise(
    n_rep             = n(),
    K_hat_mean        = mean(K_hat, na.rm = TRUE),
    ari_mean          = mean(ari, na.rm = TRUE),
    vi_mean           = mean(vi, na.rm = TRUE),
    elpd_mean         = mean(elpd_loo, na.rm = TRUE),
    looic_mean        = mean(looic, na.rm = TRUE),
    viol_rate_mean    = mean(violation_rate_mean, na.rm = TRUE),
    bf_wst_mean       = mean(bf_wst_0, na.rm = TRUE),
    bf_sst_mean       = mean(bf_sst_0, na.rm = TRUE),
    thetaW_model_mean = mean(thetaW_block_model_mean, na.rm = TRUE),
    thetaS_model_mean = mean(thetaS_block_model_mean, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Summary (averaged over reps) ---\n")
print(as.data.frame(summary_df), digits = 3)

write.csv(summary_df, file.path(out_dir, "misspec_summary.csv"),
          row.names = FALSE)

## ---- Cleanup -------------------------------------------------------------
sink(type = "message")
sink(type = "output")
close(log_con)

cat(sprintf("Done. Results in %s/\n", out_dir))
