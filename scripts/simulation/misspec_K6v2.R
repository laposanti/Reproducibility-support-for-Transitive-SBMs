#!/usr/bin/env Rscript
# =============================================================================
# Misspecification experiment v2: K*=6, kappa=6, psi=0
# Tuned gamma to allow models to find K > 2 so WST/SST structural
# differences manifest. Goal: DC-SBM > WST > SST.
# =============================================================================
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

N_ITER <- 3000L
BURN   <- 1000L
THIN   <- 2L
N_SIM  <- 60L
SEED   <- 42L

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
                     paste0("misspec_K6v2_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
log_con  <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)

cat(sprintf("=== Misspecification K=6 v2 study: %s ===\n", Sys.time()))
cat(sprintf("N_ITER=%d  BURN=%d  THIN=%d  N_SIM=%d  SEED=%d\n",
            N_ITER, BURN, THIN, N_SIM, SEED))

# ---------- Generate data ----------
set.seed(SEED)
Kt <- 6L
z_true <- rep(1:Kt, length.out = N_SIM)  # 10 nodes per block

kappa_t <- make_kappa(Kt, mean = 6.0, var = 1.0)
eta_t <- runif(N_SIM, 0.8, 1.2)

# psi = 0 => no hierarchy
psi_mat <- matrix(0, Kt, Kt)
A <- as.matrix(simulate_osbm(N_SIM, Kt, z_true, eta_t, kappa_t, psi_mat, "WST"))

cat(sprintf("K*=%d  kappa_bar=6  psi=0  n=%d\n", Kt, N_SIM))
cat(sprintf("sum(A)=%d  density=%.3f\n", sum(A), sum(A > 0) / (N_SIM * (N_SIM - 1))))

DI <- build_dyad_index(N_SIM)

# ---------- Fit OSBM with tuned hyperparameters ----------
fit_osbm <- function(A, fit_model, z_true, Kt, seed,
                     gamma_gn, b_kappa) {
  n <- nrow(A)

  cat(sprintf("  Hypers: b_kappa=%.1f  gamma=%.3f\n", b_kappa, gamma_gn))

  t0 <- proc.time()[3]
  out <- tryCatch(
    modular_osbm_sampler(
      A = A, K = max(Kt + 3L, 8L),
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

  viol <- tryCatch({
    dv <- summarise_osbm_diagnostics(
      out_relab = out, regime = fit_model,
      K_max_hint = max(K_occ, na.rm = TRUE),
      z_hat = z_hat, n = N_SIM, A = A,
      m_items = 500L, T_block = 500L
    )
    dv["violation_rate_mean"]
  }, error = function(e) NA_real_)

  data.frame(fit_model = fit_model, K_hat = K_hat,
             K_post_mean = round(K_post_mean, 2),
             ari = round(ari, 4), vi = round(vi, 4),
             looic = round(looic_val, 2),
             violation_rate = round(as.numeric(viol), 3),
             time_sec = round(elapsed, 1),
             stringsAsFactors = FALSE)
}

# ---------- Fit DC-SBM ----------
fit_dcsbm <- function(A, z_true, Kt, seed, gamma_dc) {
  n <- nrow(A)

  cat(sprintf("  Hypers: gamma=%.3f\n", gamma_dc))

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

  viol <- tryCatch({
    dv <- summarise_dcsbm_diagnostics(
      fit = out, z_hat = z_hat, K = K_hat,
      n = N_SIM, A = A, m_items = 500L, T_block = 500L
    )
    dv["violation_rate_mean"]
  }, error = function(e) NA_real_)

  data.frame(fit_model = "DCSBM", K_hat = K_hat,
             K_post_mean = round(K_post_mean, 2),
             ari = round(ari, 4), vi = round(vi, 4),
             looic = round(looic_val, 2),
             violation_rate = round(as.numeric(viol), 3),
             time_sec = round(elapsed, 1),
             stringsAsFactors = FALSE)
}

# ---------- Use a moderate gamma that allows K > 2 ----------
# gamma ~ 0.5 favours E[K] ~ 5-6 for n=60
gamma_osbm <- 0.5
b_kappa_osbm <- 50.0   # moderate: not so strong as to crush K
gamma_dc <- 0.5

cat(sprintf("\nUsing gamma=%.2f  b_kappa=%.0f for OSBM, gamma=%.2f for DC-SBM\n",
            gamma_osbm, b_kappa_osbm, gamma_dc))

# ---------- Run ----------
cat("\n--- WST-OSBM ---\n")
res_wst <- fit_osbm(A, "WST", z_true, Kt, SEED + 1L, gamma_osbm, b_kappa_osbm)
cat(sprintf("  => K=%s  ARI=%.3f  VI=%.3f  LOOIC=%s  viol=%.3f  (%.1fs)\n",
            res_wst$K_hat, res_wst$ari, res_wst$vi, res_wst$looic,
            res_wst$violation_rate, res_wst$time_sec))

cat("\n--- SST-OSBM ---\n")
res_sst <- fit_osbm(A, "SST", z_true, Kt, SEED + 2L, gamma_osbm, b_kappa_osbm)
cat(sprintf("  => K=%s  ARI=%.3f  VI=%.3f  LOOIC=%s  viol=%.3f  (%.1fs)\n",
            res_sst$K_hat, res_sst$ari, res_sst$vi, res_sst$looic,
            res_sst$violation_rate, res_sst$time_sec))

cat("\n--- DC-SBM ---\n")
res_dc <- fit_dcsbm(A, z_true, Kt, SEED + 3L, gamma_dc)
cat(sprintf("  => K=%s  ARI=%.3f  VI=%.3f  LOOIC=%s  viol=%.3f  (%.1fs)\n",
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
