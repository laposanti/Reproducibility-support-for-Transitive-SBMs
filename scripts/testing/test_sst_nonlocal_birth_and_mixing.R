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

check_nonlocal_helpers <- function() {
  psi <- c(0.4, 0.8, 1.1)
  K <- 4L
  S <- matrix(0, K, K)
  O <- matrix(0, K, K)

  # d = K-1 pair contributes to old-old distance-K totals for middle insertions.
  S[1, 4] <- 2.0
  O[1, 4] <- 3.0

  # d = 2 pair contributes fixed correction over r in [a+1, b].
  S[1, 3] <- 1.5
  O[1, 3] <- 2.0

  out <- .sst_oldold_nonlocal_birth_terms(S, O, psi)

  # For (1,4), insertion ranks 2:4 should collect distance-K old-old mass.
  stopifnot(all(abs(out$YK_old[2:4] - 2.0) < 1e-12))
  stopifnot(all(abs(out$OK_old[2:4] - 3.0) < 1e-12))

  # For (1,3), correction applied to r=2 and r=3 only.
  delta_13 <- S[1, 3] * (psi[3] - psi[2]) - 0.5 * O[1, 3] * (psi[3]^2 - psi[2]^2)
  stopifnot(abs(out$corr_fixed[2] - delta_13) < 1e-12)
  stopifnot(abs(out$corr_fixed[3] - delta_13) < 1e-12)
  stopifnot(abs(out$corr_fixed[1]) < 1e-12)
  stopifnot(abs(out$corr_fixed[5]) < 1e-12)

  TRUE
}

simulate_sst_data <- function(n = 36L, K = 4L, seed = 42L) {
  set.seed(seed)
  z_true <- rep(seq_len(K), length.out = n)
  z_true <- z_true[sample.int(n)]
  eta_true <- rep(1, n)
  kappa_true <- matrix(3.0, K, K)
  psi_true <- c(0.6, 1.05, 1.45)
  A <- as.matrix(simulate_osbm(n, K, z_true, eta_true, kappa_true, psi_true, "SST"))
  list(A = A, z_true = z_true, psi_true = psi_true, K_true = K)
}

summarize_fit <- function(fit, z_true, psi_true, label, mixing_flag, score_mode) {
  z_draws <- fit$z
  psi_draws <- fit$psi

  K_trace <- vapply(z_draws, function(z) length(unique(z)), integer(1))
  K_mode <- as.integer(names(which.max(table(K_trace))))

  ess_K <- NA_real_
  if (requireNamespace("coda", quietly = TRUE)) {
    ess_K <- as.numeric(coda::effectiveSize(K_trace))
  }
  acf1_K <- if (length(K_trace) > 2L) as.numeric(stats::acf(K_trace, plot = FALSE, lag.max = 1)$acf[2]) else NA_real_

  vi_mean <- NA_real_
  if (requireNamespace("mcclust", quietly = TRUE)) {
    vi_mean <- mean(vapply(z_draws, function(z) mcclust::vi.dist(z, z_true), numeric(1)))
  }

  dmax <- min(length(psi_true), min(vapply(psi_draws, length, integer(1))))
  psi_mean <- rep(NA_real_, length(psi_true))
  ess_psi1 <- NA_real_
  if (dmax > 0L) {
    psi_mat <- do.call(rbind, lapply(psi_draws, function(v) v[seq_len(dmax)]))
    psi_mean[seq_len(dmax)] <- colMeans(psi_mat)
    if (requireNamespace("coda", quietly = TRUE)) {
      ess_psi1 <- as.numeric(coda::effectiveSize(psi_mat[, 1]))
    }
  }

  psi_rmse <- if (dmax > 0L) sqrt(mean((psi_mean[seq_len(dmax)] - psi_true[seq_len(dmax)])^2)) else NA_real_

  data.frame(
    run = label,
    score_mode = score_mode,
    use_mixing_moves = mixing_flag,
    K_mode = K_mode,
    K_mean = mean(K_trace),
    K_sd = stats::sd(K_trace),
    ESS_K = ess_K,
    ACF1_K = acf1_K,
    VI_mean = vi_mean,
    psi1_mean = psi_mean[1],
    psi_rmse = psi_rmse,
    ESS_psi1 = ess_psi1,
    stringsAsFactors = FALSE
  )
}

run_chain <- function(dat, seed, score_mode = c("exact_nonlocal", "local_approx"), use_mixing_moves = TRUE) {
  score_mode <- match.arg(score_mode)
  modular_osbm_sampler(
    A = dat$A,
    K = dat$K_true + 1L,
    free = c("psi", "kappa", "eta", "z"),
    n_iter = 1200,
    burn = 400,
    thin = 2,
    verbose = FALSE,
    psi_constraint = "SST",
    seed = seed ,
    partition_prior = "GN",
    gamma_gn = 0.9,
    use_mixing_moves = use_mixing_moves,
    sst_birth_score_mode = score_mode
  )
}

cat("Running helper checks...\n")
stopifnot(check_nonlocal_helpers())

cat("Simulating SST benchmark data...\n")
dat <- simulate_sst_data()

cat("Running local-approx + mixing chain...\n")
fit_local <- run_chain(dat, seed = 1001, score_mode = "local_approx", use_mixing_moves = TRUE)

cat("Running exact-nonlocal + mixing chain...\n")
fit_exact_mix <- run_chain(dat, seed = 1002, score_mode = "exact_nonlocal", use_mixing_moves = TRUE)

cat("Running exact-nonlocal + no-mixing chain...\n")
fit_exact_nomix <- run_chain(dat, seed = 1003, score_mode = "exact_nonlocal", use_mixing_moves = FALSE)

res <- rbind(
  summarize_fit(fit_local, dat$z_true, dat$psi_true, "local_mix", TRUE, "local_approx"),
  summarize_fit(fit_exact_mix, dat$z_true, dat$psi_true, "exact_mix", TRUE, "exact_nonlocal"),
  summarize_fit(fit_exact_nomix, dat$z_true, dat$psi_true, "exact_nomix", FALSE, "exact_nonlocal")
)

csv_path <- file.path(run_dir, "sst_nonlocal_birth_mixing_summary.csv")
utils::write.csv(res, csv_path, row.names = FALSE)

mk_path <- file.path(run_dir, "sst_nonlocal_birth_mixing_report.md")
lines <- c(
  "# SST non-local birth: recovery and mixing check",
  "",
  "## Setup",
  "- Data: synthetic SST network with n=36, K*=4, psi_true=(0.6, 1.05, 1.45)",
  "- Sampler: `modular_osbm_sampler`, SST variable-K, 1200 iter, burn 400, thin 2",
  "- Chains compared:",
  "  - `local_approx` + mixing moves",
  "  - `exact_nonlocal` + mixing moves",
  "  - `exact_nonlocal` + no mixing moves",
  "",
  "## Results",
  paste0("- CSV summary: `", csv_path, "`"),
  "",
  "| run | score_mode | mixing | K_mode | K_mean | ESS_K | ACF1_K | VI_mean | psi1_mean | psi_rmse | ESS_psi1 |",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  lines <- c(lines, sprintf(
    "| %s | %s | %d | %d | %.3f | %.2f | %.3f | %.3f | %.3f | %.3f | %.2f |",
    r$run, r$score_mode, as.integer(r$use_mixing_moves), r$K_mode, r$K_mean,
    r$ESS_K, r$ACF1_K, r$VI_mean, r$psi1_mean, r$psi_rmse, r$ESS_psi1
  ))
}

lines <- c(lines, "", "## Interpretation",
  "- If `exact_nonlocal` improves or matches VI and psi RMSE vs `local_approx`, the non-local correction is behaving as intended.",
  "- Compare `exact_mix` vs `exact_nomix` on ESS/ACF: if mixing moves improve ESS and reduce lag-1 autocorrelation, no immediate split-merge redesign is needed for correctness.",
  "- Any need for split-merge adjustment would then be efficiency-oriented (proposal quality), not MH-validity-oriented."
)

writeLines(lines, mk_path)
cat("Wrote:\n")
cat(" - ", csv_path, "\n", sep = "")
cat(" - ", mk_path, "\n", sep = "")
