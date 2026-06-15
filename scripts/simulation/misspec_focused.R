#!/usr/bin/env Rscript
# =============================================================================
# Focused misspecification experiment
#
# Generates a directed DC-SBM with:
#   K* = 3, n = 60, theta_mean = 1, theta_sd = 0.12,
#   and a mixed lambda matrix with no global ordering.
# Fits WST-OSBM, SST-OSBM, and DC-SBM.
# 1 replicate, 3000 iterations.
# =============================================================================
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

N_ITER <- 3000L
BURN   <- 1000L
THIN   <- 2L
N_SIM  <- 60L
SEED   <- 7777L

suppressPackageStartupMessages({
  library(Matrix); library(BayesLogit); library(truncnorm)
  library(salso); library(fossil); library(mcclust); library(coda); library(loo)
  library(dplyr); library(tidyr); library(mcclust.ext)
})

source("helper_folder/sim_study_helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("helper_folder/transitivity_check_helper.R")
source("helper_folder/helper.R")
source("core/my_best_try_so_far.R")
source("core/DCSBM_varK.R")

sgn      <- function(x) ifelse(x > 0, 1L, ifelse(x < 0, -1L, 0L))
logistic <- function(x) 1 / (1 + exp(-x))

out_dir <- file.path("output", "simulation", "raw",
                     paste0("misspec_focused_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
log_con  <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)

cat(sprintf("=== Focused misspecification study: %s ===\n", Sys.time()))
cat(sprintf("N_ITER=%d  BURN=%d  THIN=%d  N_SIM=%d  SEED=%d\n",
            N_ITER, BURN, THIN, N_SIM, SEED))

# ---------- Generate data ----------
set.seed(SEED)
Kt <- 3L
z_true <- rep(1:Kt, length.out = N_SIM)

theta_mean <- 1.0
theta_sd   <- 0.12
theta_t <- rgamma(N_SIM, shape = theta_mean^2 / theta_sd^2, rate = theta_mean / theta_sd^2)
theta_t <- normalize_block_theta(z_true, theta_t)

lambda_t <- matrix(0.05, Kt, Kt)
diag(lambda_t) <- 0.15
lambda_t[1, 2] <- 0.10; lambda_t[2, 1] <- 0.50
lambda_t[2, 3] <- 0.23; lambda_t[3, 2] <- 0.34
lambda_t[1, 3] <- 0.90; lambda_t[3, 1] <- 0.22

sim <- simulate_dcsbm(N_SIM, Kt, z_true, theta_t, lambda_t, seed = SEED)
A <- as.matrix(sim$A)

cat(sprintf("K*=%d  theta_mean=%.2f  theta_sd=%.2f  n=%d\n", Kt, theta_mean, theta_sd, N_SIM))
cat(sprintf("sum(A)=%d  density=%.3f\n", sum(A), sum(A > 0) / (N_SIM * (N_SIM - 1))))

# ---------- Build dyad index once ----------
DI <- build_dyad_index(N_SIM)

# ---------- Helper to fit OSBM ----------
fit_osbm <- function(A, fit_model, z_true, Kt, seed) {
  n <- nrow(A)
  K_spec <- estimate_K_spectral(A)
  b_kappa <- 5 * n / max(K_spec, 2)
  gamma_gn <- choose_gamma_from_K_expected(n, K_expected = max(K_spec, 2),
                ordering_prior_mode = "equivalence_class")

  t0 <- proc.time()[3]
  out <- tryCatch(
    modular_osbm_sampler(
      A = A, K = max(Kt + 3L, 6L),
      n_iter = N_ITER, burn = BURN, thin = THIN,
      psi_constraint = fit_model , partition_prior = "GN",
      gamma_gn = gamma_gn,
      a_kappa = 1, b_kappa = b_kappa,
      a_eta = 1, b_eta = 1,
      mu0 = 0, sigma0 = 2, tau0 = 0.15,
      use_mixing_moves = TRUE,
      sample_b_kappa = FALSE,
      seed = seed, verbose = FALSE
    ),
    error = function(e) { message(sprintf("  %s ERROR: %s", fit_model, e$message)); NULL }
  )
  elapsed <- proc.time()[3] - t0

  if (is.null(out)) {
    return(data.frame(fit_model = fit_model, K_hat = NA, K_post_mean = NA,
                      ari = NA, vi = NA, looic = NA,
                      violation_rate = NA, time_sec = round(elapsed, 1),
                      stringsAsFactors = FALSE))
  }

  z_chain <- do.call(rbind, lapply(out$z, as.integer))
  K_occ <- apply(z_chain, 1, function(z) length(unique(z)))
  salso_cap <- as.integer(max(3L, min(ncol(z_chain), quantile(K_occ, 0.99, type = 1) + 1L)))
  z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L, maxNClusters = salso_cap)
  K_hat <- length(unique(z_hat))
  K_post_mean <- mean(K_occ)

  ari <- fossil::adj.rand.index(z_hat, z_true)
  vi  <- mcclust::vi.dist(z_hat, z_true)

  LL <- tryCatch(
    loglik_matrix_modular(A, out, regime = fit_model, dyad_index = DI),
    error = function(e) NULL
  )
  looic_val <- NA_real_
  if (!is.null(LL)) {
    loo_fit <- tryCatch(loo::loo(LL), error = function(e) NULL)
    if (!is.null(loo_fit)) looic_val <- loo_fit$estimates["looic", "Estimate"]
  }

  # Violation rate
  viol <- tryCatch({
    diag_vec <- summarise_osbm_diagnostics(
      out_relab = out, regime = fit_model,
      K_max_hint = max(K_occ, na.rm = TRUE),
      z_hat = z_hat, n = N_SIM, A = A,
      m_items = 500L, T_block = 500L
    )
    diag_vec["violation_rate_mean"]
  }, error = function(e) NA_real_)

  data.frame(fit_model = fit_model, K_hat = K_hat,
             K_post_mean = round(K_post_mean, 2),
             ari = round(ari, 4), vi = round(vi, 4),
             looic = round(looic_val, 2),
             violation_rate = round(as.numeric(viol), 3),
             time_sec = round(elapsed, 1),
             stringsAsFactors = FALSE)
}

# ---------- Helper to fit DC-SBM ----------
fit_dcsbm <- function(A, z_true, Kt, seed) {
  n <- nrow(A)
  gamma_dc <- choose_gamma_from_K_expected(n, K_expected = Kt,
                ordering_prior_mode = "equivalence_class")

  t0 <- proc.time()[3]
  out <- tryCatch(
    fit_dcsbm_gibbs_gnedin(
      A = as.matrix(A), K_init = 15L,
      priors = list(a_eta = 1, b_eta = 1, a_lambda = 1, b_lambda = 1,
                    gamma_gnedin = gamma_dc),
      iters = N_ITER, burn_in = BURN, thin = THIN,
      verbose = FALSE, seed = seed
    ),
    error = function(e) { message(sprintf("  DCSBM ERROR: %s", e$message)); NULL }
  )
  elapsed <- proc.time()[3] - t0

  if (is.null(out)) {
    return(data.frame(fit_model = "DCSBM", K_hat = NA, K_post_mean = NA,
                      ari = NA, vi = NA, looic = NA,
                      violation_rate = NA, time_sec = round(elapsed, 1),
                      stringsAsFactors = FALSE))
  }

  z_chain <- if (is.matrix(out$z)) out$z else do.call(rbind, lapply(out$z, as.integer))
  K_occ <- apply(z_chain, 1, function(z) length(unique(z)))
  salso_cap <- as.integer(max(3L, min(ncol(z_chain), quantile(K_occ, 0.99, type = 1) + 1L)))
  z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L, maxNClusters = salso_cap)
  K_hat <- length(unique(z_hat))
  K_post_mean <- mean(K_occ)

  ari <- fossil::adj.rand.index(z_hat, z_true)
  vi  <- mcclust::vi.dist(z_hat, z_true)

  LL <- tryCatch(
    loglik_matrix_dcsbm(A = as.matrix(A), dcsbm_out = out, dyad_index = DI),
    error = function(e) NULL
  )
  looic_val <- NA_real_
  if (!is.null(LL)) {
    loo_fit <- tryCatch(loo::loo(LL), error = function(e) NULL)
    if (!is.null(loo_fit)) looic_val <- loo_fit$estimates["looic", "Estimate"]
  }

  # Violation rate from DC-SBM diagnostics
  viol <- tryCatch({
    diag_vec <- summarise_dcsbm_diagnostics(
      fit = out, z_hat = z_hat, K = K_hat,
      n = N_SIM, A = A, m_items = 500L, T_block = 500L
    )
    diag_vec["violation_rate_mean"]
  }, error = function(e) NA_real_)

  data.frame(fit_model = "DCSBM", K_hat = K_hat,
             K_post_mean = round(K_post_mean, 2),
             ari = round(ari, 4), vi = round(vi, 4),
             looic = round(looic_val, 2),
             violation_rate = round(as.numeric(viol), 3),
             time_sec = round(elapsed, 1),
             stringsAsFactors = FALSE)
}

# ---------- Run all three models ----------
cat("\n--- WST-OSBM ---\n")
res_wst <- fit_osbm(A, "WST", z_true, Kt, SEED + 1L)
cat(sprintf("  K=%s  ARI=%.3f  VI=%.3f  LOOIC=%s  viol=%.3f  (%.1fs)\n",
            res_wst$K_hat, res_wst$ari, res_wst$vi, res_wst$looic,
            res_wst$violation_rate, res_wst$time_sec))

cat("\n--- SST-OSBM ---\n")
res_sst <- fit_osbm(A, "SST", z_true, Kt, SEED + 2L)
cat(sprintf("  K=%s  ARI=%.3f  VI=%.3f  LOOIC=%s  viol=%.3f  (%.1fs)\n",
            res_sst$K_hat, res_sst$ari, res_sst$vi, res_sst$looic,
            res_sst$violation_rate, res_sst$time_sec))

cat("\n--- DC-SBM ---\n")
res_dc <- fit_dcsbm(A, z_true, Kt, SEED + 3L)
cat(sprintf("  K=%s  ARI=%.3f  VI=%.3f  LOOIC=%s  viol=%.3f  (%.1fs)\n",
            res_dc$K_hat, res_dc$ari, res_dc$vi, res_dc$looic,
            res_dc$violation_rate, res_dc$time_sec))

# ---------- Summary ----------
results <- rbind(res_wst, res_sst, res_dc)
cat("\n=== Summary ===\n")
print(results, row.names = FALSE)

write.csv(results, file.path(out_dir, "misspec_results.csv"), row.names = FALSE)

cat(sprintf("\nResults saved to %s\n", out_dir))

sink(type = "message")
sink(type = "output")
close(log_con)
cat("Done.\n")
