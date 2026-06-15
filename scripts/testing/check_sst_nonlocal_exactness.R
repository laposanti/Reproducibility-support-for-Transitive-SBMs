#!/usr/bin/env Rscript

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")
source("helper_folder/helper.R")
source("helper_folder/sim_study_helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")

build_edge_data <- function(A) {
  n <- nrow(A)
  N_mat <- A + t(A)
  diag(N_mat) <- 0
  idx <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  i_idx <- idx[, 1]
  j_idx <- idx[, 2]
  N_edge <- as.numeric(N_mat[idx])
  edge_by_node <- replicate(n, integer(0), simplify = FALSE)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], as.integer(e))
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], as.integer(e))
  }
  list(i_idx = i_idx, j_idx = j_idx, N_edge = N_edge, edge_by_node = edge_by_node)
}

score_error_for_node <- function(A, z, psi, i, tau0) {
  K <- max(z)
  edge_data <- build_edge_data(A)
  omega_edge <- draw_edge_omega(
    z, psi, edge_data$i_idx, edge_data$j_idx, edge_data$N_edge,
    mode = "distance"
  )

  n_minus_full <- tabulate(z[-i], nbins = K)
  keep_full <- sort(which(n_minus_full > 0L))
  K_minus <- length(keep_full)
  if (K_minus < 1L) return(NULL)

  map <- integer(K)
  map[keep_full] <- seq_len(K_minus)
  z_packed <- map[z]
  psi_minus <- reindex_psi_sst_keep(psi_old = psi, K_old = K, keep = keep_full)

  c_plus <- numeric(K_minus)
  N_tot <- numeric(K_minus)
  for (pos in seq_len(K_minus)) {
    lab <- keep_full[pos]
    nodes <- which(z == lab)
    if (!length(nodes)) next
    c_plus[pos] <- sum(A[i, nodes])
    N_tot[pos] <- sum(A[i, nodes] + A[nodes, i])
  }
  Omega_blk <- aggregate_omega_by_block(
    i, z, omega_edge, edge_data$i_idx, edge_data$j_idx,
    edge_data$edge_by_node, K
  )[keep_full]
  oldold_stats <- .sst_oldold_pair_stats_excluding_i(
    i = i,
    z_packed = z_packed,
    A = A,
    i_idx = edge_data$i_idx,
    j_idx = edge_data$j_idx,
    N_edge = edge_data$N_edge,
    omega_edge = omega_edge,
    K_minus = K_minus
  )
  r_set <- seq_len(K_minus + 1L)

  lp_fast <- dir_pg_SST_new_exact_nonlocal(
    c_plus = c_plus,
    N_tot = N_tot,
    Omega_blk = Omega_blk,
    psi_vec = psi_minus,
    S_old = oldold_stats$S_old,
    O_old = oldold_stats$O_old,
    r_set = r_set,
    tau0 = tau0
  )
  lp_brute <- dir_pg_SST_new_exact_nonlocal_bruteforce(
    i = i,
    z_packed = z_packed,
    A = A,
    i_idx = edge_data$i_idx,
    j_idx = edge_data$j_idx,
    N_edge = edge_data$N_edge,
    omega_edge = omega_edge,
    psi_vec = psi_minus,
    r_set = r_set,
    tau0 = tau0
  )

  c(
    max_abs_error = max(abs(lp_fast - lp_brute)),
    max_abs_error_centered = max(abs((lp_fast - lp_fast[1]) - (lp_brute - lp_brute[1])))
  )
}

make_dense_stress_state <- function(n = 36L, K = 12L) {
  z <- rep(seq_len(K), each = n / K)
  A <- matrix(0, n, n)
  for (u in seq_len(n - 1L)) {
    for (v in seq.int(u + 1L, n)) {
      N_uv <- 3L + ((u + v) %% 3L)
      fwd <- 1L + ((2L * u + v) %% max(1L, N_uv - 1L))
      A[u, v] <- fwd
      A[v, u] <- N_uv - fwd
    }
  }
  list(A = A, z = as.integer(z), psi = cumsum(rep(0.2, K - 1L)))
}

set.seed(123)
tau0 <- 0.15
errors <- list()

boundary <- make_dense_stress_state(n = 12L, K = 2L)
boundary$z <- as.integer(c(1L, rep(2L, 11L)))
boundary$psi <- c(0.9)
errors[[length(errors) + 1L]] <- score_error_for_node(
  boundary$A, boundary$z, boundary$psi, i = 1L, tau0 = tau0
)

for (rep in 1:20) {
  n <- 14L
  K <- 4L
  z <- sample(rep(seq_len(K), length.out = n))
  eta <- rep(1, n)
  kappa <- matrix(2.5, K, K)
  psi <- c(0.5, 0.9, 1.2)
  A <- as.matrix(simulate_osbm(n, K, z, eta, kappa, psi, "SST"))
  errors[[length(errors) + 1L]] <- score_error_for_node(A, z, psi, sample.int(n, 1), tau0)
}

stress <- make_dense_stress_state()
for (i in seq_len(nrow(stress$A))) {
  errors[[length(errors) + 1L]] <- score_error_for_node(stress$A, stress$z, stress$psi, i, tau0)
}

err_mat <- do.call(rbind, Filter(Negate(is.null), errors))
max_err <- max(err_mat[, "max_abs_error"])
max_err_centered <- max(err_mat[, "max_abs_error_centered"])

cat(sprintf("max_abs_error=%.12f\n", max_err))
cat(sprintf("max_abs_error_centered=%.12f\n", max_err_centered))

if (max_err > 1e-8 || max_err_centered > 1e-8) {
  stop("Fast exact SST birth scorer disagrees with brute-force oracle.", call. = FALSE)
}
