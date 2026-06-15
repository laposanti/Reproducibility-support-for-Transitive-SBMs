###########################################################################
# Smoke test for mixing improvements
# Tests: adjacent-block swap, split-merge, n_inner_sweeps, random node
# order, wider slot_radius — on a small synthetic network.
###########################################################################

cat("=== Loading libraries ===\n")
suppressPackageStartupMessages({
  library(BayesLogit)
  library(truncnorm)
  library(Matrix)
  library(Rcpp)
})

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

cat("=== Sourcing helpers ===\n")
source("./helper_folder/helper.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./core/my_best_try_so_far.R")

# ---- Generate a small directed network with known block structure ----
set.seed(42)
K_true <- 3
n <- 30
z_true <- rep(1:K_true, each = n / K_true)

# Interaction rates: blocks interact more within
kappa_true <- matrix(0.5, K_true, K_true)
diag(kappa_true) <- 3
eta_true <- rep(1, n)

# Directional bias: block 1 beats 2, 2 beats 3
psi_true_sst <- c(1.0, 1.5)  # ψ_1 = 1.0, ψ_2 = 1.5 (SST monotone)

A <- matrix(0L, n, n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    ki <- z_true[i]; kj <- z_true[j]
    lam <- eta_true[i] * eta_true[j] * kappa_true[ki, kj]
    N_ij <- rpois(1, lam)
    if (N_ij > 0) {
      d <- abs(ki - kj)
      if (d == 0) {
        rho <- 0.5
      } else {
        psi_d <- psi_true_sst[d]
        sgn <- sign(kj - ki)
        rho <- 1 / (1 + exp(-sgn * psi_d))
      }
      a_ij <- rbinom(1, N_ij, rho)
      A[i, j] <- a_ij
      A[j, i] <- N_ij - a_ij
    }
  }
}

cat("Network: n =", n, ", edges =", sum(A > 0), "\n")

# =====================================================================
# Test 1: Unit test for adjacent_block_swap_move
# =====================================================================
cat("\n=== Test 1: Adjacent-block swap move ===\n")

z_test <- z_true
K_test <- K_true
kappa_test <- kappa_true
psi_test_sst <- psi_true_sst

# Build edge list
N_mat <- A + t(A); diag(N_mat) <- 0
idx <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
i_idx <- idx[, 1]; j_idx <- idx[, 2]; N_edge <- as.numeric(N_mat[idx])
edge_by_node <- replicate(n, integer(0), simplify = FALSE)
for (e in seq_len(nrow(idx))) {
  edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], as.integer(e))
  edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], as.integer(e))
}

# SST swap
n_accepted <- 0
for (rep in 1:50) {
  res <- adjacent_block_swap_move(
    z = z_test, kappa = kappA_test <- kappa_true,
    psi = psi_test_sst, eta = eta_true,
    A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
    edge_by_node = edge_by_node, psi_mode = "distance"
  )
  if (isTRUE(res$accepted)) n_accepted <- n_accepted + 1
}
cat("  SST: accepted", n_accepted, "/ 50 swaps\n")

# WST swap
psi_test_wst <- matrix(0, K_test, K_test)
psi_test_wst[1, 2] <- 1.0; psi_test_wst[1, 3] <- 1.5; psi_test_wst[2, 3] <- 1.0
psi_test_wst <- psi_test_wst + t(psi_test_wst)
diag(psi_test_wst) <- 0
# WST stores upper-tri only (symmetric)
psi_test_wst[lower.tri(psi_test_wst)] <- t(psi_test_wst)[lower.tri(psi_test_wst)]

n_accepted_wst <- 0
for (rep in 1:50) {
  res <- adjacent_block_swap_move(
    z = z_test, kappa = kappa_test,
    psi = psi_test_wst, eta = eta_true,
    A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
    edge_by_node = edge_by_node, psi_mode = "pair"
  )
  if (isTRUE(res$accepted)) n_accepted_wst <- n_accepted_wst + 1
}
cat("  WST: accepted", n_accepted_wst, "/ 50 swaps\n")

# =====================================================================
# Test 2: Unit test for split_merge_move
# =====================================================================
cat("\n=== Test 2: Split-merge move ===\n")

n_accepted_sm <- 0
n_splits <- 0; n_merges <- 0
for (rep in 1:20) {
  res <- tryCatch(
    split_merge_move(
      z = z_test, kappa = kappa_test,
      psi = psi_test_sst, eta = eta_true,
      A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = 1, b_kappa = 1, gamma_gn = 0.8,
      psi_mode = "distance",
      hyper_psi = list(mu0 = 1, sigma0 = 2, tau0 = 0.15),
      n_restricted_scans = 3L
    ),
    error = function(e) {
      cat("  Error:", conditionMessage(e), "\n")
      list(accepted = FALSE, move = "error")
    }
  )
  if (isTRUE(res$accepted)) n_accepted_sm <- n_accepted_sm + 1
  if (res$move == "split") n_splits <- n_splits + 1
  if (res$move == "merge") n_merges <- n_merges + 1
}
cat("  SST split-merge: accepted", n_accepted_sm, "/ 20 (",
    n_splits, "splits,", n_merges, "merges)\n")

# =====================================================================
# Test 3: Full sampler run (SST)
# =====================================================================
cat("\n=== Test 3: Full SST sampler (200 iter) ===\n")

t0 <- proc.time()[3]
res_sst <- tryCatch(
  modular_osbm_sampler(
    A = A, K = 5,
    free = c("psi", "kappa", "eta", "z"),
    n_iter = 200, burn = 50, thin = 1,
    verbose = TRUE,
    psi_constraint = "SST",
    seed = 123,
    gamma_gn = 0.8,
    mu0 = 1, sigma0 = 2, tau0 = 0.15
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
t1 <- proc.time()[3]

if (!is.null(res_sst)) {
  cat("\n  SST sampler completed in", round(t1 - t0, 1), "s\n")
  cat("  K trace:", head(res_sst$K_trace, 20), "...\n")
  cat("  K range:", range(res_sst$K_trace), "\n")
} else {
  cat("  SST sampler FAILED\n")
}

# =====================================================================
# Test 4: Full sampler run (WST)
# =====================================================================
cat("\n=== Test 4: Full WST sampler (200 iter) ===\n")

t0 <- proc.time()[3]
res_wst <- tryCatch(
  modular_osbm_sampler(
    A = A, K = 5,
    free = c("psi", "kappa", "eta", "z"),
    n_iter = 200, burn = 50, thin = 1,
    verbose = TRUE,
    psi_constraint = "WST",
    seed = 456,
    gamma_gn = 0.8,
    mu0 = 1, sigma0 = 2, tau0 = 0.15
  ),
  error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  }
)
t1 <- proc.time()[3]

if (!is.null(res_wst)) {
  cat("\n  WST sampler completed in", round(t1 - t0, 1), "s\n")
  cat("  K trace:", head(res_wst$K_trace, 20), "...\n")
  cat("  K range:", range(res_wst$K_trace), "\n")
} else {
  cat("  WST sampler FAILED\n")
}

cat("\n=== All tests finished ===\n")
