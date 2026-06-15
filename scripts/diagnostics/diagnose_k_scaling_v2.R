#!/usr/bin/env Rscript
# Diagnostic: investigate K-scaling anomaly — CSV analysis + standalone DGP sim
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(truncnorm); library(Matrix) })

d <- read.csv("output/simulation/raw/full_simulation_crossfit_final_DemoKvar_run_20260302_153429.csv",
              stringsAsFactors = FALSE)

cat("=== Table 1: RAW ARI and VI (matched gen==fit, plus DC-SBM) ===\n")
s1 <- d %>%
  filter(gen_model == fit_model | fit_model == "DC-SBM") %>%
  mutate(vi_norm = vi / log2(K_true)) %>%
  group_by(gen_model, fit_model, hierch, K_true) %>%
  summarise(n=n(), mean_ari=mean(ari), mean_vi=mean(vi), mean_vi_norm=mean(vi_norm),
            mean_K_hat=mean(K_hat), prop_K_ok=mean(K_hat==K_true), .groups="drop") %>%
  arrange(gen_model, fit_model, factor(hierch, levels=c("weak","medium","strong")), K_true)
options(width=180)
print(as.data.frame(s1), digits=3)

# Key patterns to check:
# 1. Does VI_norm fix the apparent anomaly?
# 2. Does strong hierarchy degrade with K while weak improves?
cat("\n=== Summary: ARI trend with K (3→5→8) ===\n")
trends <- s1 %>%
  group_by(gen_model, fit_model, hierch) %>%
  summarise(
    ari_K3 = mean_ari[K_true==3],
    ari_K5 = mean_ari[K_true==5],
    ari_K8 = mean_ari[K_true==8],
    vi_K3  = mean_vi[K_true==3],
    vi_K5  = mean_vi[K_true==5],
    vi_K8  = mean_vi[K_true==8],
    vinorm_K3 = mean_vi_norm[K_true==3],
    vinorm_K5 = mean_vi_norm[K_true==5],
    vinorm_K8 = mean_vi_norm[K_true==8],
    ari_improves = ari_K8 > ari_K3,
    vi_improves = vi_K8 < vi_K3,
    vinorm_improves = vinorm_K8 < vinorm_K3,
    .groups="drop"
  )
print(as.data.frame(trends), digits=3)

cat("\n=== Key question: does VI_norm resolve the anomaly? ===\n")
cat("Cases where raw VI drops (improves) with K: ", sum(trends$vi_improves), "/", nrow(trends), "\n")
cat("Cases where VI_norm drops (improves) with K: ", sum(trends$vinorm_improves), "/", nrow(trends), "\n")
cat("Cases where ARI rises (improves) with K:     ", sum(trends$ari_improves), "/", nrow(trends), "\n")

# ============================================================
# DGP ANALYSIS: understand the signal structure
# ============================================================
cat("\n\n=== DGP SIGNAL ANALYSIS ===\n")

set.seed(42)
n <- 60

sample_kappa_prior <- function(K, mean, var) {
  shape <- mean^2 / var; rate <- mean / var
  KAP <- matrix(0, K, K)
  for (k in 1:K) for (l in k:K) KAP[k, l] <- rgamma(1L, shape = shape, rate = rate)
  KAP[lower.tri(KAP)] <- t(KAP)[lower.tri(KAP)]; KAP
}
sample_psi_sst_prior <- function(K, mu_psi, var_psi) {
  if (K <= 1) return(numeric(0))
  cumsum(truncnorm::rtruncnorm(K - 1, a = 0, mean = mu_psi, sd = sqrt(var_psi)))
}
sample_psi_wst_prior <- function(K, mu_psi, var_psi) {
  PSI <- matrix(0, K, K)
  PSI[upper.tri(PSI)] <- truncnorm::rtruncnorm(K*(K-1)/2, a = 0, mean = mu_psi, sd = sqrt(var_psi))
  PSI[lower.tri(PSI)] <- -t(PSI)[lower.tri(PSI)]; PSI
}

kappa_mean <- 3.0; kappa_cv <- 0.6; kappa_var <- (kappa_cv * kappa_mean)^2

# Average over many draws to get expected psi values
n_draws <- 500
cat("\nExpected SST psi values (averaged over", n_draws, "draws):\n")
for (hierch_label in c("weak", "medium", "strong")) {
  psi_mean_val <- switch(hierch_label, weak = 0.2, medium = 0.7, strong = 1.3)
  cat(sprintf("\n  %s (mu_psi = %.1f):\n", hierch_label, psi_mean_val))
  
  for (K in c(3, 5, 8)) {
    psi_draws <- replicate(n_draws, sample_psi_sst_prior(K, psi_mean_val, 0.09))
    # psi_draws is (K-1) x n_draws matrix
    mean_psi <- rowMeans(psi_draws)
    cat(sprintf("    K=%d: E[psi] = %s\n", K, paste(round(mean_psi, 3), collapse=", ")))
    cat(sprintf("           max E[psi] = %.3f, p(fwd|max_dist) = %.3f\n", 
                max(mean_psi), plogis(max(mean_psi))))
    
    # Number of block-pair N: with n=60, block_size = 60/K, n_dyads = (60/K)^2
    nk <- 60 / K
    n_dyads_per_pair <- nk^2
    n_between_pairs <- K * (K - 1) / 2
    expected_N <- n_dyads_per_pair * kappa_mean
    
    # SNR for distance-1 pair
    p_d1 <- plogis(mean_psi[1])
    signal_d1 <- p_d1 - 0.5
    se_d1 <- sqrt(p_d1 * (1-p_d1) / expected_N)
    snr_d1 <- signal_d1 / se_d1
    
    # SNR for max-distance pair
    p_max <- plogis(max(mean_psi))
    signal_max <- p_max - 0.5
    se_max <- sqrt(p_max * (1-p_max) / expected_N)
    snr_max <- signal_max / se_max
    
    cat(sprintf("           n_dyads/pair = %.0f, E[N]/pair = %.0f\n", n_dyads_per_pair, expected_N))
    cat(sprintf("           d=1: p=%.3f, signal=%.3f, SE=%.4f, SNR=%.1f\n", p_d1, signal_d1, se_d1, snr_d1))
    cat(sprintf("           d=max: p=%.3f, signal=%.3f, SE=%.4f, SNR=%.1f\n", p_max, signal_max, se_max, snr_max))
  }
}

# WST analysis
cat("\n\nExpected WST psi values:\n")
for (hierch_label in c("weak", "medium", "strong")) {
  psi_mean_val <- switch(hierch_label, weak = 0.2, medium = 0.7, strong = 1.3)
  cat(sprintf("\n  %s (mu_psi = %.1f):\n", hierch_label, psi_mean_val))
  
  for (K in c(3, 5, 8)) {
    psi_vals_all <- replicate(n_draws, {
      P <- sample_psi_wst_prior(K, psi_mean_val, 0.09)
      P[upper.tri(P)]
    })
    mean_psi_wst <- mean(psi_vals_all)
    min_psi_wst <- mean(apply(psi_vals_all, 2, min))
    max_psi_wst <- mean(apply(psi_vals_all, 2, max))
    
    nk <- 60 / K
    n_between_pairs <- K * (K - 1) / 2
    cat(sprintf("    K=%d: n_between_pairs=%d, mean_psi=%.3f, min=%.3f, max=%.3f, p(fwd|mean)=%.3f\n",
                K, n_between_pairs, mean_psi_wst, min_psi_wst, max_psi_wst, plogis(mean_psi_wst)))
    
    # Key insight for WST: all psi values are drawn iid, so distribution doesn't change with K
    # But the number of "discriminating features" per node grows as K-1
    # And the N per pair shrinks as (n/K)^2
    expected_N_pair <- nk^2 * kappa_mean
    p_fwd <- plogis(mean_psi_wst)
    signal <- p_fwd - 0.5
    se <- sqrt(p_fwd*(1-p_fwd) / expected_N_pair)
    snr <- signal / se
    cat(sprintf("           N/pair=%.0f, SNR=%.1f, K-1 profiles=%d\n", expected_N_pair, snr, K-1))
    
    # Total information: SNR^2 * n_pairs
    total_info <- snr^2 * n_between_pairs
    cat(sprintf("           Total Fisher info proxy (SNR^2 * n_pairs) = %.1f\n", total_info))
  }
}

# ============================================================
# KEY INSIGHT: Compute total Fisher information about partition
# For each node i, it has K-1 between-block interaction partners.
# The info about z_i from the directionality of edges to block l is:
#   I_{i,l} = N_{i,l} * p(1-p)  (Fisher info for binomial proportion)
# Where N_{i,l} = sum_{j in block l} N_{ij} ≈ n_l * kappa
# And p = logistic(psi)
#
# Total info per node: sum_{l != z_i} I_{i,l}
# With K blocks: (K-1) terms, each with N ≈ (n/K) * kappa
# ============================================================
cat("\n\n=== TOTAL FISHER INFORMATION PER NODE ===\n")
cat("(Measures how identifiable each node's block membership is)\n\n")

for (hierch_label in c("weak", "medium", "strong")) {
  psi_mean_val <- switch(hierch_label, weak = 0.2, medium = 0.7, strong = 1.3)
  cat(sprintf("  %s (mu_psi=%.1f):\n", hierch_label, psi_mean_val))
  
  for (K in c(3, 5, 8)) {
    nk <- 60 / K  # block size
    
    # SST: each node interacts with K-1 other blocks at various distances
    # Node in block k interacts with blocks at distances 1,2,...,K-1
    # N per partner block ≈ nk * kappa
    N_per_partner <- nk * kappa_mean
    
    # Expected psi at each distance
    psi_d <- cumsum(rep(truncnorm::etruncnorm(a=0, mean=psi_mean_val, sd=sqrt(0.09)), K-1))
    p_d <- plogis(psi_d)
    
    # Fisher info per distance
    fisher_per_d <- N_per_partner * p_d * (1 - p_d)
    total_fisher_sst <- sum(fisher_per_d)
    
    # The "signal" component: deviation from 0.5
    signal_per_d <- N_per_partner * (p_d - 0.5)^2 / (p_d * (1-p_d))
    total_signal_sst <- sum(signal_per_d)
    
    cat(sprintf("    SST K=%d: N/partner=%.0f, total_Fisher=%.1f, total_signal=%.1f\n",
                K, N_per_partner, total_fisher_sst, total_signal_sst))
    
    # WST: same idea but all distances have ~same psi
    e_psi <- truncnorm::etruncnorm(a=0, mean=psi_mean_val, sd=sqrt(0.09))
    p_wst <- plogis(e_psi)
    fisher_per_pair_wst <- N_per_partner * p_wst * (1-p_wst)
    total_fisher_wst <- (K-1) * fisher_per_pair_wst
    signal_per_pair_wst <- N_per_partner * (p_wst-0.5)^2 / (p_wst*(1-p_wst))
    total_signal_wst <- (K-1) * signal_per_pair_wst
    
    cat(sprintf("    WST K=%d: N/partner=%.0f, total_Fisher=%.1f, total_signal=%.1f\n",
                K, N_per_partner, total_fisher_wst, total_signal_wst))
  }
  cat("\n")
}

cat("\n=== DECOMPOSITION: N/partner scales as n/K, but terms scale as K-1 ===\n")
cat("SST: N/partner * sum_{d=1}^{K-1} f(d) where f grows with d\n")
cat("  So both MORE terms and STRONGER terms as K grows => total info INCREASES\n")
cat("WST: N/partner * (K-1) * f(constant) \n")
cat("  N/partner = (n/K)*kappa ∝ 1/K\n")
cat("  # terms = K-1\n")
cat("  Product ∝ (K-1)/K ≈ 1 - 1/K => nearly constant, slight increase\n")
cat("  BUT: with larger K, more block-pairs create richer structure\n\n")

cat("=== STRONG HIERARCHY REVERSAL ===\n")
cat("At strong psi=1.3: p(fwd) for d=1 already ≈ 0.79\n")
cat("For K=3: max d=2, p ≈ 0.93 — very identifiable\n")
cat("For K=8: max d=7, p ≈ 0.9999 — saturated, no more info!\n")
cat("The bottleneck shifts to d=1 pairs, which have LESS data at large K\n")
cat("=> Strong hierarchy degrades because the signal saturates while data thins\n")
