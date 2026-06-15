#!/usr/bin/env Rscript
# Diagnostic script: investigate K-scaling anomaly in simulation study
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

d <- read.csv("output/simulation/raw/full_simulation_crossfit_final_DemoKvar_run_20260302_153429.csv",
              stringsAsFactors = FALSE)

cat("=== Dataset overview ===\n")
cat("Total rows:", nrow(d), "\n")
cat("K values:", sort(unique(d$K_true)), "\n")
cat("hierch:", unique(d$hierch), "\n")
cat("gen_model:", unique(d$gen_model), "\n")
cat("fit_model:", unique(d$fit_model), "\n\n")

# ============================================================
# 1. ARI and VI by K, for each (gen_model, fit_model, hierch)
#    Focus on matching gen/fit models + DC-SBM
# ============================================================
cat("=== Table 1: ARI / VI / VI_norm by K, hierch (matched models only) ===\n")
s1 <- d %>%
  filter(gen_model == fit_model | fit_model == "DC-SBM") %>%
  group_by(gen_model, fit_model, hierch, K_true) %>%
  summarise(
    n         = n(),
    mean_ari  = mean(ari, na.rm = TRUE),
    sd_ari    = sd(ari, na.rm = TRUE),
    mean_vi   = mean(vi, na.rm = TRUE),
    vi_norm   = mean(vi / log2(K_true), na.rm = TRUE),
    mean_K_hat = mean(K_hat, na.rm = TRUE),
    prop_K_ok  = mean(K_hat == K_true, na.rm = TRUE),
    mean_viol  = mean(violation_rate_zhat, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(gen_model, fit_model, hierch, K_true)

options(width = 180)
print(as.data.frame(s1), digits = 3)

# ============================================================
# 2. Also look at cross-fit: WST data fit with SST and vice versa
# ============================================================
cat("\n=== Table 2: Cross-fit ARI / VI by K (all combos) ===\n")
s2 <- d %>%
  group_by(gen_model, fit_model, hierch, K_true) %>%
  summarise(
    n         = n(),
    mean_ari  = mean(ari, na.rm = TRUE),
    mean_vi   = mean(vi, na.rm = TRUE),
    vi_norm   = mean(vi / log2(K_true), na.rm = TRUE),
    mean_K_hat = mean(K_hat, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(gen_model, fit_model, hierch, K_true)

print(as.data.frame(s2), digits = 3)

# ============================================================
# 3. Now the key: simulate the DGP statistics to understand 
#    what changes with K. No MCMC needed - just the data itself.
# ============================================================
cat("\n\n=== DATA GENERATING PROCESS DIAGNOSTICS ===\n")
cat("Simulating networks at K=3,5,8 to measure signal strength\n\n")

source("./helper_folder/sim_study_helper.R")
# Inline the DGP functions instead of sourcing DemoKvar.R (which runs the full sim)
library(truncnorm)
library(Matrix)

sample_kappa_prior <- function(K, mean, var) {
  shape <- mean^2 / var; rate <- mean / var
  KAP <- matrix(0, K, K)
  for (k in 1:K) for (l in k:K) KAP[k, l] <- rgamma(1L, shape = shape, rate = rate)
  KAP[lower.tri(KAP)] <- t(KAP)[lower.tri(KAP)]; KAP
}
sample_psi_sst_prior <- function(K, mu_psi, var_psi) {
  if (K <= 1) return(numeric(0))
  tau <- sqrt(var_psi)
  cumsum(truncnorm::rtruncnorm(K - 1, a = 0, mean = mu_psi, sd = tau))
}
sample_psi_wst_prior <- function(K, mu_psi, var_psi) {
  sigma_psi <- sqrt(var_psi)
  PSI <- matrix(0, K, K)
  PSI[upper.tri(PSI)] <- truncnorm::rtruncnorm(K*(K-1)/2, a = 0, mean = mu_psi, sd = sigma_psi)
  PSI[lower.tri(PSI)] <- -t(PSI)[lower.tri(PSI)]; PSI
}

set.seed(42)
n <- 60
kappa_mean <- 3.0
kappa_cv <- 0.6
kappa_var <- (kappa_cv * kappa_mean)^2

for (hierch_label in c("weak", "medium", "strong")) {
  psi_mean_val <- switch(hierch_label, weak = 0.2, medium = 0.7, strong = 1.3)
  psi_var_val <- 0.09

  cat(sprintf("\n--- Hierarchy: %s (psi_mean = %.1f) ---\n", hierch_label, psi_mean_val))
  
  for (K in c(3, 5, 8)) {
    z_true <- rep(seq_len(K), length.out = n)
    n_per_block <- table(z_true)
    
    # Generate kappa and psi
    kappa_true <- sample_kappa_prior(K, mean = kappa_mean, var = kappa_var)
    
    # SST psi
    psi_sst <- sample_psi_sst_prior(K, mu_psi = psi_mean_val, var_psi = psi_var_val)
    
    # WST psi
    psi_wst <- sample_psi_wst_prior(K, mu_psi = psi_mean_val, var_psi = psi_var_val)
    
    # Simulate networks
    A_sst <- as.matrix(simulate_osbm(n, K, z_true, rep(1, n), kappa_true, psi_sst, "SST"))
    A_wst <- as.matrix(simulate_osbm(n, K, z_true, rep(1, n), kappa_true, psi_wst, "WST"))
    
    # Compute statistics
    N_sst <- A_sst + t(A_sst)  # total interactions
    N_wst <- A_wst + t(A_wst)
    
    # Per block-pair statistics
    n_block_pairs <- K * (K - 1) / 2
    
    # For SST: compute per-distance statistics
    cat(sprintf("\n  K=%d: blocks of size %s, n_block_pairs=%d\n", 
                K, paste(n_per_block, collapse="/"), n_block_pairs))
    
    # SST psi values
    cat(sprintf("    SST psi values: %s\n", paste(round(psi_sst, 3), collapse=", ")))
    cat(sprintf("    SST max psi = %.3f, implied p(fwd) = %.3f\n", 
                max(psi_sst), plogis(max(psi_sst))))
    cat(sprintf("    SST min psi = %.3f, implied p(fwd) = %.3f\n", 
                min(psi_sst), plogis(min(psi_sst))))
    
    # WST psi: upper triangle values
    wst_vals <- psi_wst[upper.tri(psi_wst)]
    cat(sprintf("    WST psi range: [%.3f, %.3f], mean=%.3f\n", 
                min(wst_vals), max(wst_vals), mean(wst_vals)))
    
    # Total edges and edges per block-pair
    total_edges_sst <- sum(A_sst)
    total_edges_wst <- sum(A_wst)
    
    # Between-block edges per pair
    between_N_sst <- 0
    between_pairs <- 0
    within_N_sst <- 0
    for (k in 1:K) {
      idx_k <- which(z_true == k)
      within_N_sst <- within_N_sst + sum(N_sst[idx_k, idx_k]) / 2
      if (k < K) {
        for (l in (k+1):K) {
          idx_l <- which(z_true == l)
          between_N_sst <- between_N_sst + sum(N_sst[idx_k, idx_l])
          between_pairs <- between_pairs + 1
        }
      }
    }
    avg_N_per_pair <- between_N_sst / between_pairs
    
    # Expected dyads per between-block pair
    nk <- as.numeric(n_per_block[1])
    nl <- as.numeric(n_per_block[2])  # approximately equal
    expected_dyads <- nk * nl
    expected_N <- expected_dyads * kappa_mean
    
    cat(sprintf("    Expected dyads per block-pair: %d (n_k * n_l = %d * %d)\n",
                expected_dyads, nk, nl))
    cat(sprintf("    Expected total N per pair: %.1f (dyads * kappa_mean)\n", expected_N))
    cat(sprintf("    Observed avg N per between-pair: %.1f\n", avg_N_per_pair))
    cat(sprintf("    Total between-block edges (SST): %d\n", as.integer(between_N_sst)))
    
    # Key ratio: signal per dyad vs noise
    # For between-block pair at distance d, the "signal" is plogis(psi_d) - 0.5
    # The "noise" in estimating the proportion is ~1/sqrt(N_pair)
    for (dist in c(1, min(K-1, 4))) {
      p_fwd <- plogis(psi_sst[dist])
      signal <- p_fwd - 0.5
      noise <- sqrt(0.25 / (expected_N + 1))  # SE of binomial proportion
      snr <- signal / noise
      cat(sprintf("    SST dist=%d: p(fwd)=%.3f, signal=%.3f, noise=%.3f, SNR=%.1f\n",
                  dist, p_fwd, signal, noise, snr))
    }
    
    # Number of "directionality profiles" per node
    # Each node interacts with K-1 other blocks. Its profile is the vector of
    # forward-proportions to each other block.
    cat(sprintf("    Directional profiles dimension: %d (K-1)\n", K-1))
    
    # For WST: all psi are drawn independently, so more pairs = more info
    # Count how many block-pairs have "strong enough" signal
    n_strong_wst <- sum(abs(wst_vals) > 0.5)
    cat(sprintf("    WST: %d/%d pairs with |psi| > 0.5\n", n_strong_wst, length(wst_vals)))
  }
}

# ============================================================
# 4. The crucial analysis: compute VI upper bound and ARI lower bound
#    See if the *normalized* metrics fix the pattern
# ============================================================
cat("\n\n=== Table 3: Normalized VI comparison ===\n")
s3 <- d %>%
  filter(gen_model == fit_model | fit_model == "DC-SBM") %>%
  mutate(vi_norm = vi / log2(K_true)) %>%
  group_by(gen_model, fit_model, hierch, K_true) %>%
  summarise(
    mean_vi      = mean(vi, na.rm = TRUE),
    mean_vi_norm = mean(vi_norm, na.rm = TRUE),
    mean_ari     = mean(ari, na.rm = TRUE),
    .groups      = "drop"
  ) %>%
  arrange(gen_model, fit_model, hierch, K_true)

print(as.data.frame(s3), digits = 3)

# ============================================================
# 5. Check: is it a K-recovery issue or a clustering issue?
#    Look at K_hat vs K_true patterns
# ============================================================
cat("\n\n=== Table 4: K_hat distribution by scenario ===\n")
s4 <- d %>%
  filter(gen_model == fit_model | fit_model == "DC-SBM") %>%
  group_by(gen_model, fit_model, hierch, K_true) %>%
  summarise(
    mean_K_hat  = mean(K_hat, na.rm = TRUE),
    sd_K_hat    = sd(K_hat, na.rm = TRUE),
    prop_under  = mean(K_hat < K_true),
    prop_exact  = mean(K_hat == K_true),
    prop_over   = mean(K_hat > K_true),
    .groups     = "drop"
  ) %>%
  arrange(gen_model, fit_model, hierch, K_true)

print(as.data.frame(s4), digits = 3)

# ============================================================
# 6. Edge density analysis across K
# ============================================================
cat("\n\n=== Table 5: Network density statistics from results ===\n")
# We can compute expected density: n=60, kappa_mean varies
# E[total_edges] = sum_{i<j} E[N_ij] = sum_{i<j} eta_i*eta_j*kappa[z_i,z_j]
# With eta≈1, kappa≈kappa_mean: 
# Between: n_block_pairs * (n/K)^2 * kappa_mean
# Within:  K * C(n/K, 2) * kappa_mean
cat("Expected total interactions (approx, eta=1, all kappa=kappa_mean):\n")
for (K in c(3, 5, 8)) {
  nk <- 60 / K
  n_between_pairs <- K * (K - 1) / 2
  n_within_pairs <- K
  between_dyads <- n_between_pairs * nk^2
  within_dyads <- n_within_pairs * nk * (nk - 1) / 2
  total_dyads <- between_dyads + within_dyads
  for (km in c(0.75, 1.5, 3, 6)) {
    total_N <- total_dyads * km
    avg_N_between <- nk^2 * km
    cat(sprintf("  K=%d, kappa=%.2f: total_N=%.0f, avg_N_per_between_pair=%.1f, n_between_pairs=%d\n",
                K, km, total_N, avg_N_between, n_between_pairs))
  }
}
