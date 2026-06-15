#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("coda", quietly = TRUE)) library(coda)
})

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")
source("helper_folder/helper.R")
source("helper_folder/sim_study_helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/mixing_moves.R")
source("core/my_best_try_so_far.R")

run_dir <- "output/diagnostics/sst_nonlocal_birth"
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

simulate_sst_data <- function(n = 36L, K = 4L, seed = 42L,
                              psi_true = c(0.6, 1.05, 1.45),
                              kappa_value = 3.0) {
  set.seed(seed)
  z_true <- rep(seq_len(K), length.out = n)
  z_true <- z_true[sample.int(n)]
  eta_true <- rep(1, n)
  kappa_true <- matrix(kappa_value, K, K)
  A <- as.matrix(simulate_osbm(n, K, z_true, eta_true, kappa_true, psi_true, "SST"))
  list(A = A, z_true = z_true, eta_true = eta_true, kappa_true = kappa_true,
       psi_true = psi_true, K_true = K)
}

log1pexp_safe <- function(x) ifelse(x > 50, x, ifelse(x < -50, 0, log1p(exp(x))))

direction_stats <- function(A, z) {
  K <- length(unique(z))
  N_mat <- A + t(A)
  idx <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  D <- max(K - 1L, 0L)
  A_fwd <- numeric(D)
  N_tot <- numeric(D)
  within_N <- 0

  for (e in seq_len(nrow(idx))) {
    i <- idx[e, 1]
    j <- idx[e, 2]
    n <- N_mat[i, j]
    d <- abs(z[i] - z[j])
    if (d == 0L) {
      within_N <- within_N + n
      next
    }
    fwd <- if (z[i] < z[j]) A[i, j] else A[j, i]
    A_fwd[d] <- A_fwd[d] + fwd
    N_tot[d] <- N_tot[d] + n
  }

  p_hat <- ifelse(N_tot > 0, A_fwd / N_tot, NA_real_)
  psi_hat <- ifelse(N_tot > 0, qlogis(pmin(pmax(p_hat, 1e-6), 1 - 1e-6)), NA_real_)
  data.frame(
    distance = seq_len(D),
    N = N_tot,
    A_fwd = A_fwd,
    p_hat = p_hat,
    psi_hat = psi_hat,
    within_N = c(within_N, rep(NA_real_, max(D - 1L, 0L)))
  )
}

direction_loglik_profile <- function(A, z) {
  stats <- direction_stats(A, z)
  K <- length(unique(z))
  psi_hat <- if (K <= 1L) numeric(0) else stats$psi_hat
  psi_hat[!is.finite(psi_hat)] <- 0
  psi_hat <- cummax(pmax(psi_hat, 0))

  N_mat <- A + t(A)
  idx <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  ll <- 0
  for (e in seq_len(nrow(idx))) {
    i <- idx[e, 1]
    j <- idx[e, 2]
    n <- N_mat[i, j]
    if (z[i] == z[j]) {
      ll <- ll - n * log(2)
      next
    }
    d <- abs(z[i] - z[j])
    phi <- sign(z[j] - z[i]) * psi_hat[d]
    ll <- ll + A[i, j] * phi - n * log1pexp_safe(phi)
  }
  list(K = K, loglik = ll, psi_hat = psi_hat)
}

volume_log_marginal <- function(A, z, eta = rep(1, length(z)), a_kappa = 1, b_kappa = 1) {
  K <- length(unique(z))
  N_mat <- A + t(A)
  idx <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  R <- matrix(0, K, K)
  for (e in seq_len(nrow(idx))) {
    i <- idx[e, 1]
    j <- idx[e, 2]
    p <- min(z[i], z[j])
    q <- max(z[i], z[j])
    R[p, q] <- R[p, q] + N_mat[i, j]
  }
  E_k <- tapply(eta, factor(z, levels = seq_len(K)), sum)
  E2_k <- tapply(eta^2, factor(z, levels = seq_len(K)), sum)
  E_k[is.na(E_k)] <- 0
  E2_k[is.na(E2_k)] <- 0
  T <- outer(E_k, E_k, `*`)
  diag(T) <- pmax((E_k^2 - E2_k) / 2, 0)
  sum(gp_marginal(R[upper.tri(R, diag = TRUE)],
                  T[upper.tri(T, diag = TRUE)],
                  a_kappa, b_kappa))
}

merge_adjacent <- function(z, k) {
  z2 <- z
  z2[z2 == (k + 1L)] <- k
  z2[z2 > (k + 1L)] <- z2[z2 > (k + 1L)] - 1L
  as.integer(z2)
}

summarize_oracle_psi <- function(fit, psi_true) {
  dmax <- min(length(psi_true), min(vapply(fit$psi, length, integer(1))))
  psi_mat <- do.call(rbind, lapply(fit$psi, function(x) x[seq_len(dmax)]))
  ess <- rep(NA_real_, dmax)
  if (requireNamespace("coda", quietly = TRUE)) {
    ess <- as.numeric(coda::effectiveSize(psi_mat))
  }
  data.frame(
    distance = seq_len(dmax),
    psi_true = psi_true[seq_len(dmax)],
    psi_mean = colMeans(psi_mat),
    psi_sd = apply(psi_mat, 2, stats::sd),
    ESS = ess,
    error = colMeans(psi_mat) - psi_true[seq_len(dmax)]
  )
}

summarize_variable_fit <- function(fit, z_true, psi_true, label, tau0_value) {
  z_draws <- fit$z
  psi_draws <- fit$psi
  K_trace <- vapply(z_draws, function(z) length(unique(z)), integer(1))
  K_mode <- as.integer(names(which.max(table(K_trace))))
  vi_mean <- NA_real_
  if (requireNamespace("mcclust", quietly = TRUE)) {
    vi_mean <- mean(vapply(z_draws, function(z) mcclust::vi.dist(z, z_true), numeric(1)))
  }
  dmax <- min(length(psi_true), min(vapply(psi_draws, length, integer(1))))
  psi_rmse <- NA_real_
  psi_mean <- rep(NA_real_, length(psi_true))
  if (dmax > 0L) {
    psi_mat <- do.call(rbind, lapply(psi_draws, function(x) x[seq_len(dmax)]))
    psi_mean[seq_len(dmax)] <- colMeans(psi_mat)
    psi_rmse <- sqrt(mean((psi_mean[seq_len(dmax)] - psi_true[seq_len(dmax)])^2))
  }
  data.frame(
    run = label,
    tau0 = tau0_value,
    K_mode = K_mode,
    K_mean = mean(K_trace),
    K_sd = stats::sd(K_trace),
    VI_mean = vi_mean,
    psi1_mean = psi_mean[1],
    psi2_mean = psi_mean[2],
    psi3_mean = psi_mean[3],
    psi_rmse = psi_rmse,
    stringsAsFactors = FALSE
  )
}

cat("Simulating SST diagnostic data...\n")
dat <- simulate_sst_data()

cat("Checking empirical directional signal...\n")
true_stats <- direction_stats(dat$A, dat$z_true)
true_stats$psi_true <- dat$psi_true[true_stats$distance]
true_stats$p_true <- plogis(true_stats$psi_true)
utils::write.csv(true_stats,
                 file.path(run_dir, "sst_generator_direction_stats.csv"),
                 row.names = FALSE)

cat("Profiling directional likelihood for adjacent merges...\n")
profiles <- list(true_K4 = dat$z_true)
for (k in seq_len(dat$K_true - 1L)) {
  profiles[[paste0("merge_", k, "_", k + 1L)]] <- merge_adjacent(dat$z_true, k)
}
prof_rows <- do.call(rbind, lapply(names(profiles), function(name) {
  out <- direction_loglik_profile(dat$A, profiles[[name]])
  z_prof <- profiles[[name]]
  n_k <- tabulate(z_prof, nbins = length(unique(z_prof)))
  vol <- volume_log_marginal(dat$A, z_prof, eta = rep(1, length(z_prof)))
  gn_prior <- .gn_log_eppf(n_k, gamma_gn = 0.9)
  ocrp_prior <- .rocrp_log_eppf(n_k, theta_ocrp = 1.0)
  data.frame(
    partition = name,
    K = out$K,
    directional_loglik_profile = out$loglik,
    volume_log_marginal = vol,
    gn_prior = gn_prior,
    ocrp_prior = ocrp_prior,
    target_gn_profile = out$loglik + vol + gn_prior,
    target_ocrp_profile = out$loglik + vol + ocrp_prior,
    psi_hat = paste(sprintf("%.3f", out$psi_hat), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
prof_rows$delta_from_true <- prof_rows$directional_loglik_profile -
  prof_rows$directional_loglik_profile[prof_rows$partition == "true_K4"]
prof_rows$delta_target_gn <- prof_rows$target_gn_profile -
  prof_rows$target_gn_profile[prof_rows$partition == "true_K4"]
prof_rows$delta_target_ocrp <- prof_rows$target_ocrp_profile -
  prof_rows$target_ocrp_profile[prof_rows$partition == "true_K4"]
utils::write.csv(prof_rows,
                 file.path(run_dir, "sst_directional_partition_profile.csv"),
                 row.names = FALSE)

cat("Running fixed-z oracle chain...\n")
fit_oracle <- modular_osbm_sampler(
  A = dat$A,
  K = dat$K_true,
  truth = list(z = dat$z_true),
  free = c("psi", "kappa", "eta"),
  n_iter = 1200,
  burn = 400,
  thin = 2,
  verbose = FALSE,
  psi_constraint = "SST",
  seed = 2001 ,
  partition_prior = "GN",
  gamma_gn = 0.9,
  use_mixing_moves = FALSE
)
oracle_psi <- summarize_oracle_psi(fit_oracle, dat$psi_true)
utils::write.csv(oracle_psi,
                 file.path(run_dir, "sst_oracle_fixed_z_psi_summary.csv"),
                 row.names = FALSE)

cat("Running tau0=0.6 fixed-z oracle chain...\n")
fit_oracle_tau06 <- modular_osbm_sampler(
  A = dat$A,
  K = dat$K_true,
  truth = list(z = dat$z_true),
  free = c("psi", "kappa", "eta"),
  n_iter = 1200,
  burn = 400,
  thin = 2,
  verbose = FALSE,
  psi_constraint = "SST",
  seed = 2002 ,
  partition_prior = "GN",
  gamma_gn = 0.9,
  tau0 = 0.6,
  use_mixing_moves = FALSE
)
oracle_psi_tau06 <- summarize_oracle_psi(fit_oracle_tau06, dat$psi_true)
utils::write.csv(oracle_psi_tau06,
                 file.path(run_dir, "sst_oracle_fixed_z_psi_summary_tau06.csv"),
                 row.names = FALSE)

cat("Running tau0=0.6 exact variable-K chain...\n")
fit_exact_tau06 <- modular_osbm_sampler(
  A = dat$A,
  K = dat$K_true + 1L,
  free = c("psi", "kappa", "eta", "z"),
  n_iter = 1200,
  burn = 400,
  thin = 2,
  verbose = FALSE,
  psi_constraint = "SST",
  seed = 2003 ,
  partition_prior = "GN",
  gamma_gn = 0.9,
  tau0 = 0.6,
  use_mixing_moves = TRUE,
  sst_birth_score_mode = "exact_nonlocal"
)
tau_sensitivity <- rbind(
  summarize_variable_fit(fit_exact_tau06, dat$z_true, dat$psi_true,
                         "exact_mix_tau06", 0.6)
)

cat("Running tau0=0.6 exact variable-K chain with OCRP prior...\n")
fit_exact_tau06_ocrp <- modular_osbm_sampler(
  A = dat$A,
  K = dat$K_true + 1L,
  free = c("psi", "kappa", "eta", "z"),
  n_iter = 1200,
  burn = 400,
  thin = 2,
  verbose = FALSE,
  psi_constraint = "SST",
  seed = 2004 ,
  partition_prior = "OCRP",
  theta_ocrp = 1.0,
  tau0 = 0.6,
  use_mixing_moves = TRUE,
  sst_birth_score_mode = "exact_nonlocal"
)
tau_sensitivity <- rbind(
  tau_sensitivity,
  summarize_variable_fit(fit_exact_tau06_ocrp, dat$z_true, dat$psi_true,
                         "exact_mix_tau06_ocrp", 0.6)
)
utils::write.csv(tau_sensitivity,
                 file.path(run_dir, "sst_tau_sensitivity_summary.csv"),
                 row.names = FALSE)

report_path <- file.path(run_dir, "sst_data_sampler_diagnostic.md")
lines <- c(
  "# SST data and sampler diagnostic",
  "",
  "## Generator Direction Check",
  "",
  paste(capture.output(print(true_stats, row.names = FALSE)), collapse = "\n"),
  "",
  "## Directional Profile Likelihood",
  "",
  paste(capture.output(print(prof_rows, row.names = FALSE)), collapse = "\n"),
  "",
  "## Fixed-z Oracle Psi Recovery",
  "",
  paste(capture.output(print(oracle_psi, row.names = FALSE)), collapse = "\n"),
  "",
  "## Fixed-z Oracle Psi Recovery, tau0 = 0.6",
  "",
  paste(capture.output(print(oracle_psi_tau06, row.names = FALSE)), collapse = "\n"),
  "",
  "## Variable-K Exact Chain, tau0 = 0.6",
  "",
  paste(capture.output(print(tau_sensitivity, row.names = FALSE)), collapse = "\n")
)
writeLines(lines, report_path)

cat("Wrote:\n")
cat(" - ", file.path(run_dir, "sst_generator_direction_stats.csv"), "\n", sep = "")
cat(" - ", file.path(run_dir, "sst_directional_partition_profile.csv"), "\n", sep = "")
cat(" - ", file.path(run_dir, "sst_oracle_fixed_z_psi_summary.csv"), "\n", sep = "")
cat(" - ", file.path(run_dir, "sst_oracle_fixed_z_psi_summary_tau06.csv"), "\n", sep = "")
cat(" - ", file.path(run_dir, "sst_tau_sensitivity_summary.csv"), "\n", sep = "")
cat(" - ", report_path, "\n", sep = "")