#!/usr/bin/env Rscript
# =============================================================================
# Smoke-test: OSBM Hyperparameter Sensitivity
# Plan: "Resolve OSBM Hyperparameter Sensitivity" (March 2026)
#
# PURPOSE: Validate that (A) hierarchical b_kappa and (B) Rule A + spectral K
# both solve cross-domain hyperparameter sensitivity without manual tuning.
#
# SCENARIOS (mimicking real datasets):
#   A. Mountain-like : n=45,  K=4,  sparse   (~ShizukaMcDonald data)
#   B. Citation-like : n=47,  K=7,  dense    (~Citations data)
#   C. Medium        : n=70,  K=5,  moderate (~High School data)
#   D. Many-group    : n=50,  K=10, sparse   (stress test)
#
# METHODS:
#   1. Default b=1 (baseline — known to be bad)
#   2. Oracle Rule A: b = 5 * n / K_true (cheats with true K)
#   3. Rule A + spectral: b = 5 * n / K_spectral (no K knowledge)
#   4. Hierarchical b_kappa: alpha0=1, beta0=0.01 (learns b during MCMC)
#   5. Alt hyperprior:       alpha0=0.5, beta0=0.001 (more diffuse)
#   6. Rule A + fixed K=5:  b = 5 * n / 5 (wrong K, tests robustness)
#
# SUCCESS CRITERIA:
#   - K_mode == K_true ± 1 for methods 2-5 in all scenarios
#   - Hierarchical b_kappa mixes within burn-in (trace stabilises quickly)
#   - Rule A + spectral K_spectral is within ± 2 of K_true
# =============================================================================

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

suppressPackageStartupMessages({
  source("helper_folder/sim_study_helper.R")
  source("helper_folder/SST_helpers.R")
  source("helper_folder/WST_helpers.R")
  source("helper_folder/Hyper_setup.R")
  source("core/my_best_try_so_far.R")
})

# =============================================================================
# UNIT TEST: update_b_kappa posterior correctness
# Verify: E[b_kappa | kappa] = (alpha0 + P*a) / (beta0 + sum(kappa))
# =============================================================================
cat("=== Unit test: update_b_kappa ===\n")
set.seed(99)
K_test <- 4; a_test <- 1
alpha0_test <- 1; beta0_test <- 0.01
kappa_test <- matrix(rgamma(K_test^2, shape = 2, rate = 1), K_test, K_test)
kappa_test <- (kappa_test + t(kappa_test)) / 2
ut_idx <- upper.tri(kappa_test, diag = TRUE)
P_test   <- sum(ut_idx)
sum_kappa <- sum(kappa_test[ut_idx])
expected_mean <- (alpha0_test + P_test * a_test) / (beta0_test + sum_kappa)
samples_b <- replicate(5000, update_b_kappa(kappa_test, a_kappa = a_test,
                                            alpha_b = alpha0_test,
                                            beta_b  = beta0_test))
cat(sprintf("  E[b_kappa] analytic = %.4f, Monte Carlo = %.4f  (should match)\n",
            expected_mean, mean(samples_b)))
stopifnot(abs(mean(samples_b) - expected_mean) / expected_mean < 0.05)
cat("  PASSED\n\n")


# =============================================================================
# UNIT TEST: estimate_K_spectral sanity check
# =============================================================================
cat("=== Unit test: estimate_K_spectral ===\n")
set.seed(42)
for (K_chk in c(3, 5, 7)) {
  n_chk <- 60
  z_chk <- rep(1:K_chk, length.out = n_chk)
  kap_chk <- diag(3, K_chk) + 0.3   # strong assortative structure
  eta_chk <- rep(1, n_chk)
  psi_chk <- make_psi_wst(K_chk, "strong")
  A_chk   <- as.matrix(simulate_osbm(n_chk, K_chk, z_chk, eta_chk,
                                     kap_chk, psi_chk, "WST"))
  K_hat <- estimate_K_spectral(A_chk)
  cat(sprintf("  K_true=%d -> K_spectral=%d  (expect within 2)\n", K_chk, K_hat))
}
cat("  DONE\n\n")


# =============================================================================
# Helper: run one sampler with given b_kappa strategy
# =============================================================================
run_method <- function(A, K_true, method_label,
                       b_kappa_val    = NULL,   # fixed b_kappa (methods 1,2,3,6)
                       sample_b       = FALSE,  # TRUE for hierarchical methods (4,5)
                       alpha0_bk      = 1.0,
                       beta0_bk       = 0.01,
                       n_iter         = 1500,
                       burn           = 500,
                       thin           = 1,
                       seed           = 1L) {

  n <- nrow(A)
  a_kappa <- 1

  # Eta hyperparameters: match degree CV (robust across scenarios)
  deg <- rowSums(A) + colSums(A)
  cv  <- sd(deg) / max(mean(deg), 1e-9)
  cv  <- max(cv, 0.1)
  a_eta <- min(max(1 / cv^2, 0.5), 10)
  b_eta <- a_eta

  # Gnedin gamma: use K_true as reference (methods using spectral will differ in b only)
  gamma_gn <- choose_gamma_from_K_expected(n, K_expected = K_true,
                                           ordering_prior_mode = "equivalence_class")

  if (is.null(b_kappa_val)) {
    if (!sample_b) stop("Either b_kappa_val or sample_b=TRUE must be supplied.")
    # Hierarchical: start b_kappa at oracle (Rule A, c=5) for warm start
    b_kappa_val <- 5 * n / K_true
  }

  t0 <- proc.time()[3]
  out <- modular_osbm_sampler(
    A  = A,
    K  = max(K_true + 3L, 6L),   # generous initial K
    n_iter = n_iter, burn = burn, thin = thin,
    psi_constraint = "WST" ,
    partition_prior = "GN",
    gamma_gn       = gamma_gn,
    a_kappa        = a_kappa,
    b_kappa        = b_kappa_val,
    a_eta          = a_eta,
    b_eta          = b_eta,
    mu0            = 0,
    sigma0         = 2,
    tau0           = 0.15,
    use_mixing_moves = TRUE,
    sample_b_kappa = sample_b,
    alpha0_bkappa  = alpha0_bk,
    beta0_bkappa   = beta0_bk,
    seed           = seed,
    verbose        = FALSE
  )
  elapsed <- proc.time()[3] - t0

  K_trace <- out$K_trace
  K_mode  <- as.integer(names(which.max(table(K_trace))))
  K_q     <- quantile(K_trace, c(0.025, 0.975))

  # b_kappa summary
  bk_trace <- out$b_kappa_trace
  bk_mean  <- if (sample_b) mean(bk_trace) else b_kappa_val
  bk_sd    <- if (sample_b) sd(bk_trace)   else 0

  list(
    method   = method_label,
    K_mode   = K_mode,
    K_lo     = as.integer(K_q[1]),
    K_hi     = as.integer(K_q[2]),
    bk_mean  = bk_mean,
    bk_sd    = bk_sd,
    elapsed  = round(elapsed, 1),
    K_trace  = K_trace,
    bk_trace = bk_trace
  )
}


# =============================================================================
# Scenario definitions
# =============================================================================
scenarios <- list(
  A = list(label = "Mountain-like", n = 45, K_true = 4,  kappa_mean = 1.2, kappa_var = 0.5),
  B = list(label = "Citation-like", n = 47, K_true = 7,  kappa_mean = 4.0, kappa_var = 1.0),
  C = list(label = "Medium",        n = 70, K_true = 5,  kappa_mean = 2.0, kappa_var = 0.8),
  D = list(label = "Many-group",    n = 50, K_true = 10, kappa_mean = 0.8, kappa_var = 0.3)
)

N_REP   <- as.integer(Sys.getenv("N_REP", unset = "3"))
N_ITER  <- as.integer(Sys.getenv("N_ITER", unset = "1500"))
BURN    <- as.integer(Sys.getenv("BURN", unset = "500"))

# Collect all results
all_results <- list()

set.seed(42)
rep_seeds <- sample.int(1e5, N_REP)

for (sc_key in names(scenarios)) {
  sc <- scenarios[[sc_key]]
  cat(sprintf("\n========== Scenario %s: %s (n=%d, K=%d) ==========\n",
              sc_key, sc$label, sc$n, sc$K_true))

  sc_results <- list()

  for (rep_idx in seq_len(N_REP)) {
    seed_rep <- rep_seeds[rep_idx]
    set.seed(seed_rep)

    # Generate data
    z_true  <- rep(1:sc$K_true, length.out = sc$n)
    kappa_t <- make_kappa(sc$K_true, mean = sc$kappa_mean, var = sc$kappa_var)
    eta_t   <- runif(sc$n, 0.8, 1.2)
    psi_t   <- make_psi_wst(sc$K_true, "weak")   # weak signal — hard case
    A       <- as.matrix(simulate_osbm(sc$n, sc$K_true, z_true, eta_t, kappa_t, psi_t, "WST"))

    # Spectral K estimate for this dataset
    K_spec <- estimate_K_spectral(A)

    cat(sprintf("  Rep %d: K_spectral=%d ", rep_idx, K_spec))
    flush.console()

    rep_res <- list()

    # Method 1: Default b=1
    rep_res[["default_b1"]] <- run_method(A, sc$K_true,
      "1_default_b1", b_kappa_val = 1, seed = seed_rep,
      n_iter = N_ITER, burn = BURN)

    # Method 2: Oracle Rule A (knows K_true)
    b_oracle <- 5 * sc$n / sc$K_true
    rep_res[["oracle"]] <- run_method(A, sc$K_true,
      "2_oracle_ruleA", b_kappa_val = b_oracle, seed = seed_rep,
      n_iter = N_ITER, burn = BURN)

    # Method 3: Rule A + spectral K
    b_spectral <- 5 * sc$n / K_spec
    rep_res[["spectral"]] <- run_method(A, sc$K_true,
      "3_ruleA_spectral", b_kappa_val = b_spectral, seed = seed_rep,
      n_iter = N_ITER, burn = BURN)

    # Method 4: Hierarchical b_kappa (alpha0=1, beta0=0.01)
    rep_res[["hier_default"]] <- run_method(A, sc$K_true,
      "4_hier_a1_b001", sample_b = TRUE,
      alpha0_bk = 1.0, beta0_bk = 0.01,
      seed = seed_rep, n_iter = N_ITER, burn = BURN)

    # Method 5: Alt hyperprior (alpha0=0.5, beta0=0.001)
    rep_res[["hier_alt"]] <- run_method(A, sc$K_true,
      "5_hier_a05_b0001", sample_b = TRUE,
      alpha0_bk = 0.5, beta0_bk = 0.001,
      seed = seed_rep, n_iter = N_ITER, burn = BURN)

    # Method 6: Rule A with wrong fixed K=5 (robustness baseline)
    b_wrong <- 5 * sc$n / 5
    rep_res[["wrongK"]] <- run_method(A, sc$K_true,
      "6_ruleA_K5fixed", b_kappa_val = b_wrong, seed = seed_rep,
      n_iter = N_ITER, burn = BURN)

    cat(sprintf("| K_modes: %s\n",
      paste(sapply(rep_res, `[[`, "K_mode"), collapse = " ")))

    sc_results[[rep_idx]] <- list(rep = rep_idx, seed = seed_rep,
                                  K_spec = K_spec, methods = rep_res)
  }
  all_results[[sc_key]] <- list(scenario = sc, reps = sc_results)
}


# =============================================================================
# Summary table
# =============================================================================
cat("\n\n")
cat(strrep("=", 90), "\n")
cat("SUMMARY TABLE\n")
cat(strrep("=", 90), "\n")

method_labels <- c("default_b1", "oracle", "spectral",
                   "hier_default", "hier_alt", "wrongK")

header <- sprintf("%-20s  %-18s  %s", "Scenario", "Method", "K_mode (mean ± sd)  |  b_kappa_mean  |  time(s)")
cat(header, "\n")
cat(strrep("-", 90), "\n")

for (sc_key in names(all_results)) {
  sc    <- all_results[[sc_key]]$scenario
  reps  <- all_results[[sc_key]]$reps
  K_true <- sc$K_true

  first <- TRUE
  for (mname in method_labels) {
    K_modes <- sapply(reps, function(r) r$methods[[mname]]$K_mode)
    bk_mus  <- sapply(reps, function(r) r$methods[[mname]]$bk_mean)
    times   <- sapply(reps, function(r) r$methods[[mname]]$elapsed)
    method_lab <- reps[[1]]$methods[[mname]]$method

    K_ok <- mean(abs(K_modes - K_true) <= 1)  # fraction within ±1

    line <- sprintf("%-20s  %-18s  %4.1f ± %4.1f  (±1 hit: %d/%d)  |  b=%8.1f  |  %.1fs",
      if (first) sprintf("%s %s (%d)", sc_key, sc$label, K_true) else "",
      method_lab,
      mean(K_modes), sd(K_modes),
      sum(abs(K_modes - K_true) <= 1), N_REP,
      mean(bk_mus),
      mean(times))
    cat(line, "\n")
    first <- FALSE
  }
  cat(strrep("-", 90), "\n")
}


# =============================================================================
# Hyperprior sensitivity check (method 4 vs 5 across all scenarios)
# =============================================================================
cat("\n=== HYPERPRIOR SENSITIVITY (method 4 vs 5) ===\n")
cat("Posterior K_mode should vary <10% between alpha0=1,beta0=0.01 and alpha0=0.5,beta0=0.001\n\n")
for (sc_key in names(all_results)) {
  sc   <- all_results[[sc_key]]$scenario
  reps <- all_results[[sc_key]]$reps
  K4   <- sapply(reps, function(r) r$methods$hier_default$K_mode)
  K5   <- sapply(reps, function(r) r$methods$hier_alt$K_mode)
  pct  <- if (mean(K4) > 0) 100 * abs(mean(K4) - mean(K5)) / mean(K4) else NA
  cat(sprintf("  Scenario %s: hier_default K_mode=%.1f, hier_alt K_mode=%.1f  (diff=%.1f%%)\n",
              sc_key, mean(K4), mean(K5), pct))
}


# =============================================================================
# Backward-compatibility check: sample_b_kappa=FALSE must give same K_trace
# (seed-matched single run)
# =============================================================================
cat("\n=== BACKWARD COMPATIBILITY CHECK ===\n")
set.seed(1); n_bc <- 30; K_bc <- 3
z_bc    <- rep(1:K_bc, length.out = n_bc)
kap_bc  <- make_kappa(K_bc, 2, 1)
eta_bc  <- rep(1, n_bc)
psi_bc  <- make_psi_wst(K_bc, "strong")
A_bc    <- as.matrix(simulate_osbm(n_bc, K_bc, z_bc, eta_bc, kap_bc, psi_bc, "WST"))

gamma_bc <- choose_gamma_from_K_expected(n_bc, K_expected = K_bc,
                                         ordering_prior_mode = "equivalence_class")
b_bc <- 5 * n_bc / K_bc

out_fixed <- modular_osbm_sampler(
  A = A_bc, K = K_bc + 2L, n_iter = 300, burn = 100, thin = 1,
  psi_constraint = "WST" , partition_prior = "GN",
  gamma_gn = gamma_bc, a_kappa = 1, b_kappa = b_bc,
  a_eta = 1, b_eta = 1, mu0 = 0, sigma0 = 2, tau0 = 0.15,
  use_mixing_moves = TRUE, sample_b_kappa = FALSE, seed = 7L, verbose = FALSE)

out_hier <- modular_osbm_sampler(
  A = A_bc, K = K_bc + 2L, n_iter = 300, burn = 100, thin = 1,
  psi_constraint = "WST" , partition_prior = "GN",
  gamma_gn = gamma_bc, a_kappa = 1, b_kappa = b_bc,
  a_eta = 1, b_eta = 1, mu0 = 0, sigma0 = 2, tau0 = 0.15,
  use_mixing_moves = TRUE, sample_b_kappa = TRUE,
  alpha0_bkappa = 1.0, beta0_bkappa = 0.01,
  seed = 7L, verbose = FALSE)

# b_kappa_trace should be constant in fixed mode
stopifnot(all(out_fixed$b_kappa_trace == b_bc))
cat("  PASSED: b_kappa_trace is constant when sample_b_kappa=FALSE\n")

# Hierarchical should have varying b_kappa trace
stopifnot(sd(out_hier$b_kappa_trace) > 0)
cat("  PASSED: b_kappa_trace varies when sample_b_kappa=TRUE\n")

cat(sprintf("  Fixed K_mode=%d | Hierarchical K_mode=%d  (both should ≈ %d)\n",
  as.integer(names(which.max(table(out_fixed$K_trace)))),
  as.integer(names(which.max(table(out_hier$K_trace)))),
  K_bc))

cat("\n=== ALL CHECKS COMPLETE ===\n")
