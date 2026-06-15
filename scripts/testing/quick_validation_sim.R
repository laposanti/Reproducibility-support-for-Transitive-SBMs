#!/usr/bin/env Rscript
############################################################
## Quick validation sim: WST slot fix + mixing moves
## Compares K recovery for WST with old (restricted) vs new (full) slot_rad
## Uses 2 reps, 2000 iterations, WST K=3 scenario
############################################################
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

source("helper_folder/sim_study_helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("core/my_best_try_so_far.R")

library(ggplot2)
library(dplyr)
library(truncnorm)

sample_kappa_prior <- function(K, mean, var) {
  shape <- mean^2 / var; rate <- mean / var
  KAP <- matrix(0, K, K)
  for (k in 1:K) for (l in k:K) KAP[k, l] <- rgamma(1, shape = shape, rate = rate)
  KAP[lower.tri(KAP)] <- t(KAP)[lower.tri(KAP)]
  KAP
}

sample_psi_wst_prior <- function(K, mu_psi, var_psi) {
  sigma_psi <- sqrt(var_psi)
  PSI <- matrix(0, K, K)
  PSI[upper.tri(PSI)] <- truncnorm::rtruncnorm(
    n = K * (K - 1) / 2, a = 0, mean = mu_psi, sd = sigma_psi)
  PSI[lower.tri(PSI)] <- -t(PSI)[lower.tri(PSI)]
  PSI
}

set.seed(2025)
n <- 30; K_true <- 4
n_iter <- 2000L; burn <- 500L; thin <- 2L
n_rep <- 2L

results <- list()
K_traces_all <- list()

for (rep_id in seq_len(n_rep)) {
  seed <- 100 + rep_id
  set.seed(seed)

  z_true <- rep(1:K_true, length.out = n)
  kappa_true <- sample_kappa_prior(K_true, mean = 1.5, var = 0.3)
  eta_true <- runif(n, 0.8, 1.2)
  psi_true <- sample_psi_wst_prior(K_true, mu_psi = 1.2, var_psi = 0.3)
  A <- as.matrix(simulate_osbm(n, K_true, z_true, eta_true, kappa_true, psi_true, "WST"))

  hypers <- get_corollary_calibrated_hypers(A, K_expected = 1,
    ordering_prior_mode = "equivalence_class", a_kappa = 2, a_eta = 2,
    mu0 = 1.0, gamma_bounds = c(0.3, 0.7))

  cat(sprintf("Rep %d: Running WST with mixing moves (K_init=5)...\n", rep_id))
  set.seed(seed + 1000)
  t0 <- proc.time()[3]
  out <- tryCatch(
    modular_osbm_sampler(
      A = A, K = 5, n_iter = n_iter, burn = burn, thin = thin,
      verbose = FALSE, psi_constraint = "WST" , partition_prior = "GN",
      gamma_gn = hypers$gamma, a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
      a_eta = hypers$a_eta, b_eta = hypers$b_eta,
      mu0 = hypers$mu0, sigma0 = hypers$sigma0, tau0 = hypers$tau0,
      use_mixing_moves = TRUE, seed = seed + 1000
    ),
    error = function(e) { message("ERROR: ", e$message); NULL }
  )
  elapsed <- proc.time()[3] - t0

  if (is.null(out)) { cat("  FAILED\n"); next }

  z_chain <- do.call(rbind, lapply(out$z, as.integer))
  K_trace <- apply(z_chain, 1, function(x) length(unique(x)))
  z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L)
  ari <- fossil::adj.rand.index(z_hat, z_true)
  vi <- mcclust::vi.dist(z_hat, z_true)
  K_hat <- length(unique(z_hat))
  K_mode <- as.integer(names(which.max(table(K_trace))))

  cat(sprintf("  ARI=%.3f VI=%.3f K_hat=%d K_mode=%d (%.1fs)\n",
              ari, vi, K_hat, K_mode, elapsed))

  results[[length(results) + 1]] <- data.frame(
    rep = rep_id, ari = ari, vi = vi, K_hat = K_hat,
    K_mode = K_mode, K_mean = mean(K_trace), time = elapsed)

  K_traces_all[[length(K_traces_all) + 1]] <- data.frame(
    rep = rep_id, iter = seq_along(K_trace), K = K_trace)
}

results_df <- bind_rows(results)
K_traces_df <- bind_rows(K_traces_all)

cat("\n========== QUICK VALIDATION RESULTS ==========\n")
cat(sprintf("K_true = %d, n = %d, WST with mixing moves\n", K_true, n))
print(results_df)
cat(sprintf("\nMean ARI: %.3f\n", mean(results_df$ari)))
cat(sprintf("K_mode correct: %d/%d\n", sum(results_df$K_mode == K_true), nrow(results_df)))

# Plot K traces
dir.create("plots", showWarnings = FALSE)
p <- ggplot(K_traces_df, aes(x = iter, y = K, color = factor(rep))) +
  geom_line(alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = K_true, linetype = "dashed", color = "black") +
  labs(x = "Iteration (post burn-in)", y = "K",
       title = sprintf("WST K trace (K_true=%d, n=%d, mixing moves ON)", K_true, n),
       color = "Rep") +
  theme_minimal()
ggsave("output/simulation/plots/quick_validation_K_trace.png", p, width = 8, height = 4, dpi = 150)
cat("Saved: output/simulation/plots/quick_validation_K_trace.png\n")
cat("Done.\n")
