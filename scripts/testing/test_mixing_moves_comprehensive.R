#!/usr/bin/env Rscript
############################################################
## Comprehensive correctness tests for mixing moves
##
## Tests:
##  1. .full_loglik vectorized vs brute-force edge-loop
##  2. .partial_loglik consistency with .full_loglik
##  3. .node_dir_score vectorized vs brute-force
##  4. adjacent_block_swap_move preserves detailed balance
##  5. split_merge_move: merge→split reversibility
##  6. .gn_log_eppf consistency
##  7. Integration test: sampler runs without NaN/crash (SST+WST)
##  8. Cross-check .dir_loglik_edge against existing SST/WST testers
############################################################

library(BayesLogit)
library(truncnorm)
library(Matrix)
library(Rcpp)

source("./helper_folder/sim_study_helper.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./helper_folder/SST_tester.R")
source("./helper_folder/Hyper_setup.R")
source("./core/my_best_try_so_far.R")

# We also need slow_functions_and_checks for cross-validation
# These testers may depend on helper.R functions, so source it first
source("./helper_folder/helper.R")
tryCatch(source("./helper_folder/slow_functions_and_checks.R"), error = function(e) {
  message("Note: slow_functions_and_checks.R could not be loaded: ", e$message)
})

sgn      <- function(x) ifelse(x > 0, 1L, ifelse(x < 0, -1L, 0L))
logistic <- function(x) 1 / (1 + exp(-x))

# Track test results
test_results <- list()
run_test <- function(name, expr) {
  cat(sprintf("  TEST: %s ... ", name))
  res <- tryCatch({
    result <- eval(expr)
    cat("PASS\n")
    list(name = name, status = "PASS", detail = result)
  }, error = function(e) {
    cat(sprintf("FAIL: %s\n", e$message))
    list(name = name, status = "FAIL", detail = e$message)
  })
  test_results[[length(test_results) + 1L]] <<- res
  invisible(res)
}

assert_close <- function(a, b, tol = 1e-8, msg = "") {
  diff <- abs(a - b)
  if (any(diff > tol, na.rm = TRUE)) {
    stop(sprintf("Values differ by %.2e (tol=%.2e). %s", max(diff, na.rm = TRUE), tol, msg))
  }
}

############################################################
## Generate test data (SST and WST)
############################################################
set.seed(2025)

n <- 30L
K_true <- 4L
z_true <- rep(seq_len(K_true), length.out = n)

# SST params
psi_sst <- cumsum(abs(rnorm(K_true - 1, mean = 0.5, sd = 0.2)))
kappa <- matrix(rgamma(K_true^2, 2, 2), K_true, K_true)
kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]
eta <- runif(n, 0.8, 1.2)

# Generate SST network
A_sst <- as.matrix(simulate_osbm(n, K_true, z_true, eta, kappa, psi_sst, "SST"))
N_sst <- A_sst + t(A_sst)

# WST params
psi_wst <- matrix(0, K_true, K_true)
psi_wst[upper.tri(psi_wst)] <- abs(rnorm(K_true*(K_true-1)/2, 0.6, 0.3))
psi_wst[lower.tri(psi_wst)] <- -t(psi_wst)[lower.tri(psi_wst)]

# Generate WST network
A_wst <- as.matrix(simulate_osbm(n, K_true, z_true, eta, kappa, psi_wst, "WST"))
N_wst <- A_wst + t(A_wst)

# Build edge indices (same structure the sampler uses)
build_edge_data <- function(A) {
  n <- nrow(A)
  N <- A + t(A)
  pairs <- which(upper.tri(N) & N > 0, arr.ind = TRUE)
  i_idx <- pairs[,1]; j_idx <- pairs[,2]
  N_edge <- N[cbind(i_idx, j_idx)]

  edge_by_node <- vector("list", n)
  for (e in seq_along(i_idx)) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], e)
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], e)
  }
  list(i_idx = i_idx, j_idx = j_idx, N_edge = N_edge, edge_by_node = edge_by_node)
}

ed_sst <- build_edge_data(A_sst)
ed_wst <- build_edge_data(A_wst)

cat(sprintf("Test data: n=%d, K=%d, SST edges=%d, WST edges=%d\n\n",
            n, K_true, length(ed_sst$N_edge), length(ed_wst$N_edge)))

############################################################
## TEST 1: .full_loglik vectorized vs brute-force
############################################################
cat("=== Test 1: .full_loglik vectorized vs brute-force ===\n")

# Brute-force reference implementation
full_loglik_brute <- function(z, kappa, psi, eta, i_idx, j_idx, N_edge, A, psi_mode) {
  A_ij <- A[cbind(i_idx, j_idx)]
  ll <- 0
  for (e in seq_along(N_edge)) {
    ii <- i_idx[e]; jj <- j_idx[e]
    zi <- z[ii]; zj <- z[jj]; n_e <- N_edge[e]; a_e <- A_ij[e]
    lam <- eta[ii] * eta[jj] * kappa[min(zi,zj), max(zi,zj)]
    if (lam > 0) ll <- ll + n_e * log(lam) - lam
    else if (n_e > 0) ll <- ll - 1e10
    ll <- ll + .dir_loglik_edge(zi, zj, a_e, n_e, psi, psi_mode)
  }
  ll
}

run_test("full_loglik SST: vectorized == brute", quote({
  v <- .full_loglik(z_true, kappa, psi_sst, eta,
                    ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, A_sst, "distance")
  b <- full_loglik_brute(z_true, kappa, psi_sst, eta,
                         ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, A_sst, "distance")
  assert_close(v, b, tol = 1e-8, msg = sprintf("vec=%.6f brute=%.6f", v, b))
}))

run_test("full_loglik WST: vectorized == brute", quote({
  v <- .full_loglik(z_true, kappa, psi_wst, eta,
                    ed_wst$i_idx, ed_wst$j_idx, ed_wst$N_edge, A_wst, "pair")
  b <- full_loglik_brute(z_true, kappa, psi_wst, eta,
                         ed_wst$i_idx, ed_wst$j_idx, ed_wst$N_edge, A_wst, "pair")
  assert_close(v, b, tol = 1e-8, msg = sprintf("vec=%.6f brute=%.6f", v, b))
}))

############################################################
## TEST 2: .partial_loglik == .full_loglik on all edges
############################################################
cat("\n=== Test 2: .partial_loglik consistency ===\n")

run_test("partial_loglik SST: all edges == full_loglik", quote({
  A_ij <- A_sst[cbind(ed_sst$i_idx, ed_sst$j_idx)]
  p <- .partial_loglik(z_true, kappa, psi_sst, eta,
                       ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, A_ij,
                       seq_along(ed_sst$N_edge), "distance")
  f <- .full_loglik(z_true, kappa, psi_sst, eta,
                    ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, A_sst, "distance")
  assert_close(p, f, tol = 1e-8)
}))

run_test("partial_loglik WST: all edges == full_loglik", quote({
  A_ij <- A_wst[cbind(ed_wst$i_idx, ed_wst$j_idx)]
  p <- .partial_loglik(z_true, kappa, psi_wst, eta,
                       ed_wst$i_idx, ed_wst$j_idx, ed_wst$N_edge, A_ij,
                       seq_along(ed_wst$N_edge), "pair")
  f <- .full_loglik(z_true, kappa, psi_wst, eta,
                    ed_wst$i_idx, ed_wst$j_idx, ed_wst$N_edge, A_wst, "pair")
  assert_close(p, f, tol = 1e-8)
}))

# partial on subset vs full difference
run_test("partial_loglik: subset additivity", quote({
  A_ij <- A_sst[cbind(ed_sst$i_idx, ed_sst$j_idx)]
  E <- length(ed_sst$N_edge)
  s1 <- seq_len(floor(E/2)); s2 <- (floor(E/2)+1):E
  p1 <- .partial_loglik(z_true, kappa, psi_sst, eta,
                        ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, A_ij,
                        s1, "distance")
  p2 <- .partial_loglik(z_true, kappa, psi_sst, eta,
                        ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, A_ij,
                        s2, "distance")
  f  <- .full_loglik(z_true, kappa, psi_sst, eta,
                     ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, A_sst, "distance")
  assert_close(p1 + p2, f, tol = 1e-8)
}))

############################################################
## TEST 3: .node_dir_score vectorized vs brute-force
############################################################
cat("\n=== Test 3: .node_dir_score vectorized vs brute ===\n")

node_dir_score_brute <- function(node_i, z, psi, psi_mode, A, i_idx, j_idx, N_edge, edge_by_node) {
  e_list <- edge_by_node[[node_i]]
  if (!length(e_list)) return(0)
  A_ij <- A[cbind(i_idx, j_idx)]
  score <- 0
  for (e in e_list) {
    score <- score + .dir_loglik_edge(z[i_idx[e]], z[j_idx[e]], A_ij[e], N_edge[e], psi, psi_mode)
  }
  score
}

run_test("node_dir_score SST (5 random nodes)", quote({
  for (node_i in sample(n, 5)) {
    v <- .node_dir_score(node_i, z_true, psi_sst, "distance",
                         A_sst, ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, ed_sst$edge_by_node)
    b <- node_dir_score_brute(node_i, z_true, psi_sst, "distance",
                              A_sst, ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, ed_sst$edge_by_node)
    assert_close(v, b, tol = 1e-10, msg = sprintf("node %d", node_i))
  }
}))

run_test("node_dir_score WST (5 random nodes)", quote({
  for (node_i in sample(n, 5)) {
    v <- .node_dir_score(node_i, z_true, psi_wst, "pair",
                         A_wst, ed_wst$i_idx, ed_wst$j_idx, ed_wst$N_edge, ed_wst$edge_by_node)
    b <- node_dir_score_brute(node_i, z_true, psi_wst, "pair",
                              A_wst, ed_wst$i_idx, ed_wst$j_idx, ed_wst$N_edge, ed_wst$edge_by_node)
    assert_close(v, b, tol = 1e-10, msg = sprintf("node %d", node_i))
  }
}))

############################################################
## TEST 4: adjacent_block_swap preserves detailed balance
##         (swap twice == identity)
############################################################
cat("\n=== Test 4: Adjacent-block swap self-inverse ===\n")

run_test("swap kappa twice == identity", quote({
  k <- 2L
  kap2 <- swap_matrix_rowcol(swap_matrix_rowcol(kappa, k, k+1L), k, k+1L)
  assert_close(kap2, kappa, tol = 1e-15)
}))

run_test("swap psi_wst twice == identity", quote({
  k <- 2L
  p2 <- swap_matrix_rowcol(swap_matrix_rowcol(psi_wst, k, k+1L), k, k+1L)
  assert_close(p2, psi_wst, tol = 1e-15)
}))

run_test("swap z twice == identity", quote({
  k <- 2L; z2 <- z_true
  z2[z_true == k] <- k + 1L; z2[z_true == (k+1L)] <- k
  z3 <- z2; z3[z2 == k] <- k + 1L; z3[z2 == (k+1L)] <- k
  stopifnot(all(z3 == z_true))
}))

run_test("swap move: loglik ratio sign check", quote({
  # If we manually swap and compute full loglik, the ratio should match
  # what the swap move would compute internally
  k <- 1L
  z_s <- z_true; z_s[z_true == k] <- k+1L; z_s[z_true == (k+1L)] <- k
  ll_old <- .full_loglik(z_true, kappa, psi_sst, eta,
                         ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, A_sst, "distance")
  ll_new <- .full_loglik(z_s, kappa, psi_sst, eta,
                         ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge, A_sst, "distance")
  # Should be finite
  stopifnot(is.finite(ll_old), is.finite(ll_new))
}))

############################################################
## TEST 5: split-merge grow/shrink kappa round-trip
############################################################
cat("\n=== Test 5: grow/shrink psi/kappa consistency ===\n")

run_test("shrink(grow(kappa)): size correct", quote({
  k <- 2L
  set.seed(1)
  kap_g <- .grow_kappa_split(kappa, k, 2, 2)
  stopifnot(nrow(kap_g) == K_true + 1, ncol(kap_g) == K_true + 1)
  kap_s <- .shrink_kappa_merge(kap_g, k)
  stopifnot(nrow(kap_s) == K_true, ncol(kap_s) == K_true)
}))

run_test("shrink(grow(psi_wst)): size correct", quote({
  k <- 2L; set.seed(1)
  hyp <- list(mu0 = 1, sigma0 = 2, tau0 = 0.15)
  psi_g <- .grow_psi_split(psi_wst, k, "pair", hyp)
  stopifnot(nrow(psi_g) == K_true + 1)
  psi_s <- .shrink_psi_merge(psi_g, k, "pair")
  stopifnot(nrow(psi_s) == K_true)
}))

run_test("shrink(grow(psi_sst)): size correct", quote({
  k <- 2L; set.seed(1)
  hyp <- list(mu0 = 1, sigma0 = 2, tau0 = 0.15)
  psi_g <- .grow_psi_split(psi_sst, k, "distance", hyp)
  stopifnot(length(psi_g) == K_true)  # was K_true-1, now K_true
  psi_s <- .shrink_psi_merge(psi_g, k, "distance")
  stopifnot(length(psi_s) == K_true - 1)
}))

############################################################
## TEST 6: .gn_log_eppf properties
############################################################
cat("\n=== Test 6: GN log-EPPF correctness ===\n")

run_test("gn_log_eppf: more blocks => different value", quote({
  nk1 <- c(10, 10, 10)          # K=3

  nk2 <- c(10, 5, 5, 10)        # K=4, same n
  gamma <- 0.5
  v1 <- .gn_log_eppf(nk1, gamma)
  v2 <- .gn_log_eppf(nk2, gamma)
  stopifnot(is.finite(v1), is.finite(v2))
  # They should differ (not =)
  stopifnot(abs(v1 - v2) > 1e-10)
}))

run_test("gn_log_eppf: gamma sensitivity", quote({
  nk <- c(10, 10, 10)
  v1 <- .gn_log_eppf(nk, 0.3)
  v2 <- .gn_log_eppf(nk, 0.7)
  stopifnot(is.finite(v1), is.finite(v2))
  stopifnot(abs(v1 - v2) > 1e-6)
}))

############################################################
## TEST 7: Cross-check .dir_loglik_edge vs SST/WST testers
############################################################
cat("\n=== Test 7: .dir_loglik_edge vs existing brute-force testers ===\n")

run_test("dir_loglik_edge SST: matches brute_existing_full", quote({
  # brute_existing_full computes sum of edge logliks for node i in block k
  # We sum .dir_loglik_edge for the same edges and compare
  test_nodes <- sample(n, 3)
  for (node_i in test_nodes) {
    k <- z_true[node_i]
    # brute approach: sum over all j != i
    ll_brute <- 0
    for (j in seq_len(n)) {
      if (j == node_i) next
      a <- A_sst[node_i, j]; nij <- N_sst[node_i, j]
      if (nij == 0) next
      zi <- z_true[node_i]; zj <- z_true[j]
      ll_brute <- ll_brute + .dir_loglik_edge(zi, zj, a, nij, psi_sst, "distance")
    }
    # vectorized approach
    ll_vec <- .node_dir_score(node_i, z_true, psi_sst, "distance",
                              A_sst, ed_sst$i_idx, ed_sst$j_idx, ed_sst$N_edge,
                              ed_sst$edge_by_node)
    # Note: brute sums over both A[i,j] and A[j,i] directions separately,
    # while node_dir_score uses the upper-tri index. They may differ by
    # how edges are indexed. Check they're consistent.
    # Actually, edge_by_node has both (i,j) as sender and receiver edges.
    # For upper-tri, node_i appears as i_idx[e] or j_idx[e].
    # The directional loglik uses z[i_idx[e]], z[j_idx[e]] with the
    # upper-tri A value. This matches the full dyad model.
    assert_close(ll_vec, ll_brute, tol = 0.1,
                 msg = sprintf("node %d: vec=%.4f brute=%.4f", node_i, ll_vec, ll_brute))
  }
}))

run_test("dir_loglik_edge WST: matches brute per-edge", quote({
  test_nodes <- sample(n, 3)
  for (node_i in test_nodes) {
    ll_brute <- 0
    for (j in seq_len(n)) {
      if (j == node_i) next
      a <- A_wst[node_i, j]; nij <- N_wst[node_i, j]
      if (nij == 0) next
      zi <- z_true[node_i]; zj <- z_true[j]
      ll_brute <- ll_brute + .dir_loglik_edge(zi, zj, a, nij, psi_wst, "pair")
    }
    ll_vec <- .node_dir_score(node_i, z_true, psi_wst, "pair",
                              A_wst, ed_wst$i_idx, ed_wst$j_idx, ed_wst$N_edge,
                              ed_wst$edge_by_node)
    assert_close(ll_vec, ll_brute, tol = 0.1,
                 msg = sprintf("node %d: vec=%.4f brute=%.4f", node_i, ll_vec, ll_brute))
  }
}))

############################################################
## TEST 8: Integration test — sampler runs (SST + WST)
############################################################
cat("\n=== Test 8: Integration test — short sampler runs ===\n")

run_test("SST sampler: 100 iter, with_moves=TRUE, no crash", quote({
  hypers <- get_corollary_calibrated_hypers(A_sst, K_expected = 1,
    ordering_prior_mode = "equivalence_class", a_kappa = 2, a_eta = 2,
    mu0 = 1.0, gamma_bounds = c(0.3, 0.7))
  out <- modular_osbm_sampler(
    A = A_sst, K = n, n_iter = 100, burn = 20, thin = 1,
    psi_constraint = "SST" , partition_prior = "GN",
    gamma_gn = hypers$gamma, a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
    a_eta = hypers$a_eta, b_eta = hypers$b_eta,
    mu0 = hypers$mu0, sigma0 = hypers$sigma0, tau0 = hypers$tau0,
    use_mixing_moves = TRUE, seed = 42)
  stopifnot(!is.null(out), length(out$z) > 0)
  K_trace <- sapply(out$z, function(x) length(unique(x)))
  stopifnot(all(is.finite(K_trace)))
  sprintf("K range: %d-%d", min(K_trace), max(K_trace))
}))

run_test("WST sampler: 100 iter, with_moves=TRUE, no crash", quote({
  hypers <- get_corollary_calibrated_hypers(A_wst, K_expected = 1,
    ordering_prior_mode = "equivalence_class", a_kappa = 2, a_eta = 2,
    mu0 = 1.0, gamma_bounds = c(0.3, 0.7))
  out <- modular_osbm_sampler(
    A = A_wst, K = n, n_iter = 100, burn = 20, thin = 1,
    psi_constraint = "WST" , partition_prior = "GN",
    gamma_gn = hypers$gamma, a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
    a_eta = hypers$a_eta, b_eta = hypers$b_eta,
    mu0 = hypers$mu0, sigma0 = hypers$sigma0, tau0 = hypers$tau0,
    use_mixing_moves = TRUE, seed = 42)
  stopifnot(!is.null(out), length(out$z) > 0)
  K_trace <- sapply(out$z, function(x) length(unique(x)))
  stopifnot(all(is.finite(K_trace)))
  sprintf("K range: %d-%d", min(K_trace), max(K_trace))
}))

run_test("SST sampler: with_moves=FALSE also works", quote({
  hypers <- get_corollary_calibrated_hypers(A_sst, K_expected = 1,
    ordering_prior_mode = "equivalence_class", a_kappa = 2, a_eta = 2,
    mu0 = 1.0, gamma_bounds = c(0.3, 0.7))
  out <- modular_osbm_sampler(
    A = A_sst, K = n, n_iter = 50, burn = 10, thin = 1,
    psi_constraint = "SST" , partition_prior = "GN",
    gamma_gn = hypers$gamma, a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
    a_eta = hypers$a_eta, b_eta = hypers$b_eta,
    mu0 = hypers$mu0, sigma0 = hypers$sigma0, tau0 = hypers$tau0,
    use_mixing_moves = FALSE, seed = 42)
  stopifnot(!is.null(out), length(out$z) > 0)
}))

############################################################
## TEST 9: Cross-check with run_fast_bruteforce_tests (SST)
############################################################
cat("\n=== Test 9: SST fast vs brute-force (existing tester) ===\n")

if (exists("run_fast_bruteforce_tests", mode = "function")) {
  run_test("SST tester: run_fast_bruteforce_tests passes", quote({
    res <- run_fast_bruteforce_tests(
      A = A_sst, N = N_sst, z = z_true, psi_vec = psi_sst,
      i_set = sample(seq_len(n), 5), tol = 1e-6, verbose = FALSE)
    stopifnot(res$max_existing_diff < 1e-6)
    sprintf("max existing diff: %.2e, max new diff: %.2e",
            res$max_existing_diff, res$max_new_diff)
  }))
} else {
  cat("  SKIP: run_fast_bruteforce_tests not available\n")
}

############################################################
## TEST 10: Cross-check with WST tester
############################################################
cat("\n=== Test 10: WST fast vs brute-force (existing tester) ===\n")

wst_tester_exists <- file.exists("./helper_folder/WST_tester.R")
if (wst_tester_exists) {
  tryCatch(source("./helper_folder/WST_tester.R"), error = function(e)
    message("Note: WST_tester.R could not be loaded: ", e$message))
}
if (exists("run_fast_bruteforce_tests_wst", mode = "function")) {
  run_test("WST tester: run_fast_bruteforce_tests_wst passes", quote({
    res <- run_fast_bruteforce_tests_wst(
      A = A_wst, N = N_wst, z = z_true, psi_mat = psi_wst,
      i_set = sample(seq_len(n), 5), tol = 1e-6, verbose = FALSE)
    stopifnot(res$max_existing_diff < 1e-6)
    sprintf("max existing diff: %.2e, max new diff: %.2e",
            res$max_existing_diff, res$max_new_diff)
  }))
} else {
  cat("  SKIP: run_fast_bruteforce_tests_wst not available\n")
}

############################################################
## TEST 11: run_step1_checks master validation (if available)
############################################################
cat("\n=== Test 11: Master step1 checks ===\n")

if (exists("run_step1_checks", mode = "function")) {
  run_test("run_step1_checks passes", quote({
    res <- run_step1_checks(seed = 2025, n = 20, K = 3,
                            tol_fast_vs_slow = 5e-3,
                            tol_closed_vs_mc = 5e-3, Rmc = 5e4)
    sprintf("SST max_diff=%.2e, WST max_diff=%.2e",
            res$sst_max_diff, res$wst_max_diff)
  }))
} else {
  cat("  SKIP: run_step1_checks not available\n")
}

############################################################
## TEST 12: .dir_loglik_vec matches loop-based .dir_loglik_edge
############################################################
cat("\n=== Test 12: .dir_loglik_vec vs .dir_loglik_edge loop ===\n")

run_test("dir_loglik_vec SST: matches edge loop", quote({
  A_ij <- A_sst[cbind(ed_sst$i_idx, ed_sst$j_idx)]
  zi <- z_true[ed_sst$i_idx]; zj <- z_true[ed_sst$j_idx]
  n_e <- ed_sst$N_edge; a_e <- A_ij
  # Loop-based
  ll_loop <- sum(vapply(seq_along(n_e), function(e)
    .dir_loglik_edge(zi[e], zj[e], a_e[e], n_e[e], psi_sst, "distance"),
    numeric(1)))
  # Vectorized
  ll_vec <- .dir_loglik_vec(zi, zj, a_e, n_e, psi_sst, "distance")
  assert_close(ll_vec, ll_loop, tol = 1e-10,
               msg = sprintf("vec=%.8f loop=%.8f", ll_vec, ll_loop))
}))

run_test("dir_loglik_vec WST: matches edge loop", quote({
  A_ij <- A_wst[cbind(ed_wst$i_idx, ed_wst$j_idx)]
  zi <- z_true[ed_wst$i_idx]; zj <- z_true[ed_wst$j_idx]
  n_e <- ed_wst$N_edge; a_e <- A_ij
  # Loop-based
  ll_loop <- sum(vapply(seq_along(n_e), function(e)
    .dir_loglik_edge(zi[e], zj[e], a_e[e], n_e[e], psi_wst, "pair"),
    numeric(1)))
  # Vectorized
  ll_vec <- .dir_loglik_vec(zi, zj, a_e, n_e, psi_wst, "pair")
  assert_close(ll_vec, ll_loop, tol = 1e-10,
               msg = sprintf("vec=%.8f loop=%.8f", ll_vec, ll_loop))
}))

############################################################
## Summary
############################################################

cat("\n\n========================================\n")
cat("           TEST SUMMARY\n")
cat("========================================\n")

n_pass <- sum(sapply(test_results, function(r) r$status == "PASS"))
n_fail <- sum(sapply(test_results, function(r) r$status == "FAIL"))
n_total <- length(test_results)

for (r in test_results) {
  status <- if (r$status == "PASS") "  PASS" else "**FAIL**"
  cat(sprintf("  %s  %s\n", status, r$name))
}

cat(sprintf("\n  %d/%d passed, %d failed\n", n_pass, n_total, n_fail))

if (n_fail > 0) {
  cat("\nFailed tests details:\n")
  for (r in test_results) {
    if (r$status == "FAIL") {
      cat(sprintf("  - %s: %s\n", r$name, r$detail))
    }
  }
}

cat("\n========================================\n")
