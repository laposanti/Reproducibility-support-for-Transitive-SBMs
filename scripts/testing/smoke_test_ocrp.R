#!/usr/bin/env Rscript
# Quick smoke test: source the sampler, run OCRP prior for a few iterations
cat("=== OCRP Smoke Test ===\n")

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")
source("helper_folder/helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("core/my_best_try_so_far.R")

# --- Test 1: ocrp_log_weights_packed basic correctness ---
cat("\n--- Test 1: ocrp_log_weights_packed ---\n")
v <- c(5, 3, 2)
theta <- 1.0
pw <- ocrp_log_weights_packed(v, theta_ocrp = theta)
cat("  exist (log):", pw$exist, "\n")
cat("  new   (log):", pw$new, "\n")
cat("  per_slot:", pw$per_slot, "\n")

# Check sum-to-1
p_exist <- exp(pw$exist)
p_new   <- exp(pw$new)
total <- sum(p_exist) + sum(p_new)
cat("  sum of probs:", total, "\n")
stopifnot(abs(total - 1.0) < 1e-10)
cat("  PASS: probabilities sum to 1\n")

# --- Test 2: run OCRP sampler (SST) for a few iterations ---
cat("\n--- Test 2: SST sampler with OCRP prior (tiny network) ---\n")
set.seed(42)
n <- 20
K_true <- 3
z_true <- rep(1:K_true, length.out = n)
A <- matrix(0, n, n)
for (i in 1:(n-1)) for (j in (i+1):n) {
  rate <- if (z_true[i] == z_true[j]) 5 else 1
  A[i, j] <- rpois(1, rate)
}
A[lower.tri(A)] <- 0

tryCatch({
  fit <- modular_osbm_sampler(
    A = A, K = 5,
    psi_constraint = "SST",
    partition_prior = "OCRP",
    theta_ocrp = 1.0,
    n_iter = 50, burn = 10, thin = 1,
    verbose = FALSE,
    use_mixing_moves = FALSE,
    seed = 123
  )
  cat("  Sampler ran successfully!\n")
  cat("  Final K values:", tail(fit$K_trace, 10), "\n")
  cat("  PASS\n")
}, error = function(e) {
  cat("  FAIL:", conditionMessage(e), "\n")
})

# --- Test 3: run OCRP sampler (WST) for a few iterations ---
cat("\n--- Test 3: WST sampler with OCRP prior (tiny network) ---\n")
tryCatch({
  fit_wst <- modular_osbm_sampler(
    A = A, K = 5,
    psi_constraint = "WST",
    partition_prior = "OCRP",
    theta_ocrp = 1.0,
    n_iter = 50, burn = 10, thin = 1,
    verbose = FALSE,
    use_mixing_moves = FALSE,
    seed = 456
  )
  cat("  Sampler ran successfully!\n")
  cat("  Final K values:", tail(fit_wst$K_trace, 10), "\n")
  cat("  PASS\n")
}, error = function(e) {
  cat("  FAIL:", conditionMessage(e), "\n")
})

# --- Test 4: SST with mixing moves ---
cat("\n--- Test 4: SST sampler with OCRP + mixing moves ---\n")
tryCatch({
  fit_mm <- modular_osbm_sampler(
    A = A, K = 5,
    psi_constraint = "SST",
    partition_prior = "OCRP",
    theta_ocrp = 1.0,
    n_iter = 50, burn = 10, thin = 1,
    verbose = FALSE,
    use_mixing_moves = TRUE,
    seed = 789
  )
  cat("  Sampler ran successfully!\n")
  cat("  Final K values:", tail(fit_mm$K_trace, 10), "\n")
  cat("  PASS\n")
}, error = function(e) {
  cat("  FAIL:", conditionMessage(e), "\n")
})

cat("\n=== Smoke test complete ===\n")
