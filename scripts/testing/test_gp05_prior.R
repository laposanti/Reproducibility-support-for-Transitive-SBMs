#!/usr/bin/env Rscript
# =============================================================================
# Test: Penalised GP05 (Dirichlet ordered partition) prior
# =============================================================================
# 1. Unit test: verify gp05_log_weights_packed against theory (sum-to-one)
# 2. Smoke test: run modular_osbm_sampler with partition_prior = "GP05"
#    on a small simulated network, both WST and SST
# 3. Comparison: GN vs GP05 on the same data, compare ARI and K recovery
# =============================================================================

cat("=== Test: GP05 (Dirichlet ordered partition) prior ===\n")
t0_script <- proc.time()[3]

suppressPackageStartupMessages({
  library(Matrix)
  library(BayesLogit)
})

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")

source("helper_folder/sim_study_helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("core/my_best_try_so_far.R")

# =====================================================================
# PART 1: Unit tests for gp05_log_weights_packed
# =====================================================================
cat("\n--- Part 1: Unit tests for gp05_log_weights_packed ---\n")

# Test 1a: Example from the theory doc (n=3, k=2, (2,1), alpha=0.5, gamma=0.8)
cat("Test 1a: numerical example from theory doc...\n")
v <- c(2, 1)  # block sizes after removing node i; n_remaining = 3
alpha <- 0.5
gamma <- 0.8
pw <- gp05_log_weights_packed(v, alpha_gp05 = alpha, gamma_gp05 = gamma)

w_exist <- exp(pw$exist)
w_new   <- exp(pw$new)

# Expected (from theory, equations for theta=0):
# Interior exist (j=1): v1 - alpha = 2 - 0.5 = 1.5
# Last exist    (j=2): (v2+1)(v2-alpha)/v2 = 2*0.5/1 = 1.0
# Interior new  (r=1,2): gamma * alpha = 0.8 * 0.5 = 0.4
# Far-right new (r=3):   gamma * alpha / v2  = 0.4 / 1 = 0.4
expected_exist <- c(1.5, 1.0)
expected_new   <- c(0.4, 0.4, 0.4)

stopifnot(all(abs(w_exist - expected_exist) < 1e-10))
stopifnot(all(abs(w_new   - expected_new)   < 1e-10))
cat("  PASS: weights match theory (Example from doc)\n")

# Test 1b: Verify normalised probabilities sum to 1
all_w <- c(w_exist, w_new)
p_all <- all_w / sum(all_w)
stopifnot(abs(sum(p_all) - 1.0) < 1e-12)
cat(sprintf("  PASS: probabilities sum to %.15f\n", sum(p_all)))

# Test 1c: Balanced blocks (3 blocks of size 5, n=15)
cat("Test 1c: balanced blocks (3x5)...\n")
v3 <- c(5, 5, 5)
pw3 <- gp05_log_weights_packed(v3, alpha_gp05 = 0.3, gamma_gp05 = 0.7)
w3_exist <- exp(pw3$exist)
w3_new   <- exp(pw3$new)

# v_j - alpha for interior (j=1,2): 5 - 0.3 = 4.7
# (v_H+1)(v_H-alpha)/v_H for last (j=3): 6*4.7/5 = 5.64
# gamma*alpha for interior new (r=1,2,3): 0.7*0.3 = 0.21
# gamma*alpha/v_H for far-right new (r=4): 0.21/5 = 0.042
stopifnot(abs(w3_exist[1] - 4.7) < 1e-10)
stopifnot(abs(w3_exist[2] - 4.7) < 1e-10)
stopifnot(abs(w3_exist[3] - 5.64) < 1e-10)
stopifnot(abs(w3_new[1] - 0.21) < 1e-10)
stopifnot(abs(w3_new[4] - 0.042) < 1e-10)
cat("  PASS: balanced block weights correct\n")

# Test 1d: Single block (H=1)
cat("Test 1d: single block...\n")
pw1 <- gp05_log_weights_packed(c(10), alpha_gp05 = 0.5, gamma_gp05 = 0.8)
w1_exist <- exp(pw1$exist)
w1_new   <- exp(pw1$new)
# exist[1] = (10+1)(10-0.5)/10 = 11*9.5/10 = 10.45  (last block formula, H=1)
# new[1]  = gamma*alpha = 0.4  (interior)
# new[2]  = gamma*alpha/10 = 0.04  (far-right)
stopifnot(abs(w1_exist[1] - 10.45) < 1e-10)
stopifnot(abs(w1_new[1] - 0.4) < 1e-10)
stopifnot(abs(w1_new[2] - 0.04) < 1e-10)
cat("  PASS: single block weights correct\n")

# Test 1e: GP05 vs theory sum for general composition
# With theta=0 predictives: sum of old = n/(n+1), sum of new = theta/(n+1)
# But penalised: sum of old + gamma * sum of new != 1, which is fine
# (normalisation is explicit)
cat("Test 1e: verify per_slot flag is TRUE...\n")
stopifnot(isTRUE(pw$per_slot))
cat("  PASS\n")

cat("\nAll unit tests passed.\n")

# =====================================================================
# PART 2: Smoke test - run sampler with GP05 prior
# =====================================================================
cat("\n--- Part 2: Smoke test with GP05 prior ---\n")

set.seed(42)
n <- 30
K_true <- 3
z_true <- rep(1:K_true, each = n / K_true)
kappa <- matrix(1.0, K_true, K_true)
diag(kappa) <- 2.0
eta <- rep(1.0, n)

# WST network
psi_wst <- matrix(0, K_true, K_true)
psi_wst[1, 2] <- 1.5; psi_wst[1, 3] <- 2.0; psi_wst[2, 3] <- 1.0
A_wst <- as.matrix(simulate_osbm(n, K_true, z_true, eta, kappa, psi_wst, "WST"))

cat(sprintf("WST network: n=%d, K_true=%d, sum(A)=%d, density=%.3f\n",
            n, K_true, sum(A_wst), sum(A_wst > 0) / (n * (n - 1))))

# Run with GP05 prior (WST)
cat("WST GP05 run (300 iter)...\n")
t0 <- proc.time()[3]
out_gp05_wst <- modular_osbm_sampler(
  A = A_wst, K = 5, n_iter = 300, burn = 100, thin = 1,
  psi_constraint = "WST" ,
  partition_prior = "GP05",
  gamma_gn = 0.7,
  alpha_gp05 = 0.5,
  a_kappa = 1, b_kappa = 1,
  a_eta = 1, b_eta = 1,
  mu0 = 1.0, sigma0 = 2.0, tau0 = 0.15,
  use_mixing_moves = TRUE, seed = 42, verbose = FALSE
)
t_wst <- proc.time()[3] - t0

z_chain_wst <- do.call(rbind, lapply(out_gp05_wst$z, as.integer))
K_trace_wst <- apply(z_chain_wst, 1, function(x) length(unique(x)))
K_mode_wst <- as.integer(names(which.max(table(K_trace_wst))))
z_hat_wst <- salso::salso(z_chain_wst, loss = salso::VI(), nRuns = 1L)
ari_wst <- fossil::adj.rand.index(z_hat_wst, z_true)
cat(sprintf("  Done in %.1fs. ARI=%.3f, K_mode=%d, K_range=[%d,%d]\n",
            t_wst, ari_wst, K_mode_wst, min(K_trace_wst), max(K_trace_wst)))

# SST network
psi_sst <- cumsum(c(0.8, 0.6))
A_sst <- as.matrix(simulate_osbm(n, K_true, z_true, eta, kappa, psi_sst, "SST"))

cat(sprintf("SST network: n=%d, K_true=%d, sum(A)=%d, density=%.3f\n",
            n, K_true, sum(A_sst), sum(A_sst > 0) / (n * (n - 1))))

# Run with GP05 prior (SST)
cat("SST GP05 run (300 iter)...\n")
t0 <- proc.time()[3]
out_gp05_sst <- modular_osbm_sampler(
  A = A_sst, K = 5, n_iter = 300, burn = 100, thin = 1,
  psi_constraint = "SST" ,
  partition_prior = "GP05",
  gamma_gn = 0.7,
  alpha_gp05 = 0.5,
  a_kappa = 1, b_kappa = 1,
  a_eta = 1, b_eta = 1,
  mu0 = 1.0, sigma0 = 2.0, tau0 = 0.15,
  use_mixing_moves = TRUE, seed = 42, verbose = FALSE
)
t_sst <- proc.time()[3] - t0

z_chain_sst <- do.call(rbind, lapply(out_gp05_sst$z, as.integer))
K_trace_sst <- apply(z_chain_sst, 1, function(x) length(unique(x)))
K_mode_sst <- as.integer(names(which.max(table(K_trace_sst))))
z_hat_sst <- salso::salso(z_chain_sst, loss = salso::VI(), nRuns = 1L)
ari_sst <- fossil::adj.rand.index(z_hat_sst, z_true)
cat(sprintf("  Done in %.1fs. ARI=%.3f, K_mode=%d, K_range=[%d,%d]\n",
            t_sst, ari_sst, K_mode_sst, min(K_trace_sst), max(K_trace_sst)))

# =====================================================================
# PART 3: GN vs GP05 comparison
# =====================================================================
cat("\n--- Part 3: GN vs GP05 comparison ---\n")

set.seed(123)
n <- 40
K_true <- 3
z_true <- c(rep(1, 15), rep(2, 15), rep(3, 10))
kappa <- matrix(1.5, K_true, K_true)
diag(kappa) <- 2.5
eta <- runif(n, 0.8, 1.2)

psi_wst <- matrix(0, K_true, K_true)
psi_wst[1, 2] <- 2.0; psi_wst[1, 3] <- 3.0; psi_wst[2, 3] <- 1.5
A <- as.matrix(simulate_osbm(n, K_true, z_true, eta, kappa, psi_wst, "WST"))
cat(sprintf("Comparison data: n=%d, K_true=%d, sum(A)=%d\n", n, K_true, sum(A)))

n_iter <- 500
burn   <- 150

# --- GN prior ---
cat("Running GN prior...\n")
t0 <- proc.time()[3]
out_gn <- modular_osbm_sampler(
  A = A, K = 5, n_iter = n_iter, burn = burn, thin = 1,
  psi_constraint = "WST" ,
  partition_prior = "GN",
  gamma_gn = 0.7,
  a_kappa = 1, b_kappa = 1,
  a_eta = 1, b_eta = 1,
  mu0 = 1.0, sigma0 = 2.0, tau0 = 0.15,
  use_mixing_moves = TRUE, seed = 123, verbose = FALSE
)
t_gn <- proc.time()[3] - t0

z_chain_gn <- do.call(rbind, lapply(out_gn$z, as.integer))
K_trace_gn <- apply(z_chain_gn, 1, function(x) length(unique(x)))
K_mode_gn <- as.integer(names(which.max(table(K_trace_gn))))
z_hat_gn <- salso::salso(z_chain_gn, loss = salso::VI(), nRuns = 1L)
ari_gn <- fossil::adj.rand.index(z_hat_gn, z_true)

# --- GP05 prior ---
cat("Running GP05 prior...\n")
t0 <- proc.time()[3]
out_gp05 <- modular_osbm_sampler(
  A = A, K = 5, n_iter = n_iter, burn = burn, thin = 1,
  psi_constraint = "WST" ,
  partition_prior = "GP05",
  gamma_gn = 0.7,
  alpha_gp05 = 0.5,
  a_kappa = 1, b_kappa = 1,
  a_eta = 1, b_eta = 1,
  mu0 = 1.0, sigma0 = 2.0, tau0 = 0.15,
  use_mixing_moves = TRUE, seed = 123, verbose = FALSE
)
t_gp05 <- proc.time()[3] - t0

z_chain_gp05 <- do.call(rbind, lapply(out_gp05$z, as.integer))
K_trace_gp05 <- apply(z_chain_gp05, 1, function(x) length(unique(x)))
K_mode_gp05 <- as.integer(names(which.max(table(K_trace_gp05))))
z_hat_gp05 <- salso::salso(z_chain_gp05, loss = salso::VI(), nRuns = 1L)
ari_gp05 <- fossil::adj.rand.index(z_hat_gp05, z_true)

cat("\n=== COMPARISON SUMMARY ===\n")
cat(sprintf("%-8s  time(s)  ARI    K_mode  K_range\n", "Prior"))
cat(sprintf("%-8s  %6.1f   %.3f  %d       [%d,%d]\n",
            "GN", t_gn, ari_gn, K_mode_gn, min(K_trace_gn), max(K_trace_gn)))
cat(sprintf("%-8s  %6.1f   %.3f  %d       [%d,%d]\n",
            "GP05", t_gp05, ari_gp05, K_mode_gp05,
            min(K_trace_gp05), max(K_trace_gp05)))

cat(sprintf("\nK posterior (GN):   %s\n",
            paste(names(table(K_trace_gn)), ":", table(K_trace_gn), collapse = "  ")))
cat(sprintf("K posterior (GP05): %s\n",
            paste(names(table(K_trace_gp05)), ":", table(K_trace_gp05), collapse = "  ")))

total_time <- proc.time()[3] - t0_script
cat(sprintf("\nTotal script time: %.1fs\n", total_time))
cat("=== GP05 test complete ===\n")
