# Metropolis-Hastings mixing moves used by the ordered SBM samplers.
#
###########################################################################
#  Mixing moves for the Ordered SBM sampler
#  - Adjacent-block swap (MH)  [WST + SST]
#  - Split-merge (MH)          [WST + SST]
#
#  These are called once per iteration from the main Gibbs loop,
#  AFTER the single-site z-sweep and BEFORE the kappa/eta updates.
###########################################################################

# =====================================================================
# (1)  ADJACENT-BLOCK SWAP  (Metropolis–Hastings)
# =====================================================================
#
# Propose swapping labels k <-> k+1 for a uniformly chosen k in 1:(K-1).
# The proposal is deterministic & symmetric, so acceptance = min(1, lik ratio).
#
# Volume part (kappa): we also swap rows/cols of kappa, so the volume
# likelihood is *unchanged* — every dyad (i,j) still sees the same
# kappa_{z_i, z_j} because both labels and matrix are permuted.
#
# Directional part: this is where the change happens.
# We only need to evaluate the directional log-lik difference for
# dyads that involve at least one node in block k or k+1.
#
# Under the uniform-over-orderings prior, the prior ratio is 1.
# Under GN the partition prior depends only on block sizes (exchangeable),
# so the prior ratio is also 1.

adjacent_block_swap_move <- function(
    z, kappa, psi, eta,
    A, i_idx, j_idx, N_edge, edge_by_node,
    psi_mode = c("pair", "distance")
) {
  psi_mode <- match.arg(psi_mode)
  K <- if (psi_mode == "pair") nrow(psi) else length(psi) + 1L
  if (K < 2L) return(list(z = z, kappa = kappa, psi = psi, accepted = FALSE))

  # Pick a random adjacent pair

  k <- sample.int(K - 1L, 1L)

  # Identify edges affected: at least one endpoint in block k or k+1
  zl <- z[i_idx]; zr <- z[j_idx]
  affected <- (zl == k | zl == (k + 1L) | zr == k | zr == (k + 1L))

  if (!any(affected)) {
    # No edges between/within these blocks => swap is free, always accept
    z_new <- z
    z_new[z == k]       <- k + 1L
    z_new[z == (k + 1L)] <- k
    kappa_new <- swap_matrix_rowcol(kappa, k, k + 1L)
    psi_new <- if (psi_mode == "pair") swap_matrix_rowcol(psi, k, k + 1L) else psi
    return(list(z = z_new, kappa = kappa_new, psi = psi_new, accepted = TRUE))
  }

  # --- Compute directional log-lik ratio for affected edges (VECTORIZED) ---
  A_ij <- A[cbind(i_idx, j_idx)]
  idx_aff <- which(affected)

  ii_a <- i_idx[idx_aff]; jj_a <- j_idx[idx_aff]
  zi_a <- z[ii_a]; zj_a <- z[jj_a]
  n_a  <- N_edge[idx_aff]; a_a <- A_ij[idx_aff]

  # Old directional log-lik (fully vectorized)
  loglik_old <- .dir_loglik_vec(zi_a, zj_a, a_a, n_a, psi, psi_mode)

  # Proposed labels
  zi_new <- ifelse(zi_a == k, k + 1L, ifelse(zi_a == (k + 1L), k, zi_a))
  zj_new <- ifelse(zj_a == k, k + 1L, ifelse(zj_a == (k + 1L), k, zj_a))

  if (psi_mode == "pair") {
    # For WST: swap psi rows/cols once, then use vectorized computation
    psi_swapped <- swap_matrix_rowcol(psi, k, k + 1L)
    loglik_new <- .dir_loglik_vec(zi_new, zj_new, a_a, n_a, psi_swapped, psi_mode)
  } else {
    loglik_new <- .dir_loglik_vec(zi_new, zj_new, a_a, n_a, psi, psi_mode)
  }

  log_alpha <- loglik_new - loglik_old

  if (log(runif(1)) < log_alpha) {
    z_new <- z
    z_new[z == k]       <- k + 1L
    z_new[z == (k + 1L)] <- k
    kappa_new <- swap_matrix_rowcol(kappa, k, k + 1L)
    psi_new <- if (psi_mode == "pair") swap_matrix_rowcol(psi, k, k + 1L) else psi
    return(list(z = z_new, kappa = kappa_new, psi = psi_new, accepted = TRUE))
  }

  list(z = z, kappa = kappa, psi = psi, accepted = FALSE)
}


# --- helpers for adjacent swap ---

.swap_label <- function(zi, k) {
  if (zi == k) k + 1L
  else if (zi == (k + 1L)) k
  else zi
}

swap_matrix_rowcol <- function(M, k, l) {
  M_new <- M
  M_new[k, ] <- M[l, ]
  M_new[l, ] <- M[k, ]
  tmp <- M_new[, k]
  M_new[, k] <- M_new[, l]
  M_new[, l] <- tmp
  M_new
}

# Directional log-likelihood for a single edge (i,j) with labels (zi, zj)
# WST: psi is KxK matrix. SST: psi is length K-1 vector.
# psi_swap_k: if not NULL, use psi with rows/cols k,k+1 swapped (for WST proposal)
.dir_loglik_edge <- function(zi, zj, a_ij, n_ij, psi, psi_mode,
                             psi_swap_k = NULL) {
  if (n_ij == 0) return(0)
  if (zi == zj) return(-n_ij * log(2))  # within-block: rho = 0.5

  if (psi_mode == "pair") {
    k_lo <- min(zi, zj); k_hi <- max(zi, zj)
    sgn <- sign(zj - zi)
    if (!is.null(psi_swap_k)) {
      # Look up psi value after swapping rows/cols psi_swap_k and psi_swap_k+1
      # This is equivalent to remapping the indices
      k_lo_orig <- .swap_label(k_lo, psi_swap_k)
      k_hi_orig <- .swap_label(k_hi, psi_swap_k)
      psi_val <- psi[min(k_lo_orig, k_hi_orig), max(k_lo_orig, k_hi_orig)]
    } else {
      psi_val <- psi[k_lo, k_hi]
    }
  } else {
    d <- abs(zi - zj)
    if (d < 1 || d > length(psi)) return(0)
    psi_val <- psi[d]
    sgn <- sign(zj - zi)
  }

  phi <- sgn * psi_val
  # log p(A_ij | N_ij, phi) = A_ij * phi - N_ij * log(1 + exp(phi))
  a_ij * phi - n_ij * .log1pexp_safe(phi)
}

.log1pexp_safe <- function(x) {
  ifelse(x > 50, x, ifelse(x < -50, 0, log1p(exp(x))))
}

# Vectorized directional log-likelihood for a batch of edges
# Returns scalar sum. zi, zj, a, n are equal-length vectors.
.dir_loglik_vec <- function(zi, zj, a, n, psi, psi_mode) {
  same <- (zi == zj)
  ll <- ifelse(same, -n * log(2), 0)

  diff_mask <- !same & (n > 0)
  if (any(diff_mask)) {
    zi_d <- zi[diff_mask]; zj_d <- zj[diff_mask]
    a_d  <- a[diff_mask];  n_d  <- n[diff_mask]

    if (psi_mode == "pair") {
      klo <- pmin(zi_d, zj_d); khi <- pmax(zi_d, zj_d)
      psi_val <- psi[cbind(klo, khi)]
    } else {
      d <- abs(zi_d - zj_d)
      d <- pmin(d, length(psi)); d <- pmax(d, 1L)
      psi_val <- psi[d]
    }
    sgn <- sign(zj_d - zi_d)
    phi <- sgn * psi_val
    ll[diff_mask] <- a_d * phi - n_d * .log1pexp_safe(phi)
  }
  sum(ll)
}


# =====================================================================
# (2)  SPLIT-MERGE MOVE  (Metropolis–Hastings)
# =====================================================================
#
# Two sub-moves, chosen with equal probability:
#   MERGE: pick two adjacent blocks k, k+1 and merge them.
#   SPLIT: pick a block k with ≥ 2 nodes and split it into k, k+1.
#
# For dimension-matching, split proposes node assignments via a
# restricted Gibbs scan, and merge is the reverse.
#
# The volume part (kappa) is handled by collapsing kappa out
# (using the Gamma-Poisson marginal), so we don't need to match
# kappa proposals. For psi, merge removes one parameter and split
# adds one — we integrate/propose appropriately.

split_merge_move <- function(
    z, kappa, psi, eta,
    A, i_idx, j_idx, N_edge, edge_by_node,
    a_kappa, b_kappa,
    gamma_gn,
    psi_mode = c("pair", "distance"),
    hyper_psi = list(mu0 = 1, sigma0 = 2, tau0 = 0.15),
    n_restricted_scans = 3L,
    partition_prior = "OCRP",
    theta_ocrp = 1.0
) {
  psi_mode <- match.arg(psi_mode)
  K <- if (psi_mode == "pair") nrow(psi) else length(psi) + 1L
  n <- length(z)

  # Choose split or merge with equal probability
  # Merge requires K >= 2; Split requires at least one block with >= 2 nodes
  n_k <- tabulate(z, nbins = K)
  can_merge <- (K >= 2L)
  can_split <- any(n_k >= 2L)

  if (!can_merge && !can_split) {
    return(list(z = z, kappa = kappa, psi = psi, accepted = FALSE, move = "none"))
  }

  if (can_merge && can_split) {
    do_split <- (runif(1) < 0.5)
  } else if (can_split) {
    do_split <- TRUE
  } else {
    do_split <- FALSE
  }

  if (do_split) {
    .split_move(z, kappa, psi, eta, A, i_idx, j_idx, N_edge, edge_by_node,
                a_kappa, b_kappa, gamma_gn, psi_mode, hyper_psi, n_restricted_scans,
                partition_prior = partition_prior, theta_ocrp = theta_ocrp)
  } else {
    .merge_move(z, kappa, psi, eta, A, i_idx, j_idx, N_edge, edge_by_node,
                a_kappa, b_kappa, gamma_gn, psi_mode, hyper_psi, n_restricted_scans,
                partition_prior = partition_prior, theta_ocrp = theta_ocrp)
  }
}


# --- MERGE: merge adjacent blocks k and k+1 into block k ---
.merge_move <- function(z, kappa, psi, eta, A, i_idx, j_idx, N_edge, edge_by_node,
                        a_kappa, b_kappa, gamma_gn, psi_mode, hyper_psi,
                        n_restricted_scans,
                        partition_prior = "OCRP", theta_ocrp = 1.0) {
  K <- if (psi_mode == "pair") nrow(psi) else length(psi) + 1L
  if (K < 2L) return(list(z = z, kappa = kappa, psi = psi, accepted = FALSE, move = "merge"))
  n <- length(z)

  # Pick random adjacent pair to merge
  k <- sample.int(K - 1L, 1L)
  nodes_k  <- which(z == k)
  nodes_k1 <- which(z == (k + 1L))
  nodes_merged <- c(nodes_k, nodes_k1)
  n_merged <- length(nodes_merged)

  if (n_merged < 2L) {
    return(list(z = z, kappa = kappa, psi = psi, accepted = FALSE, move = "merge"))
  }

  # --- Proposed state: merge k+1 into k, then shift labels > k+1 down ---
  z_prop <- z
  z_prop[z == (k + 1L)] <- k
  z_prop[z_prop > (k + 1L)] <- z_prop[z_prop > (k + 1L)] - 1L
  K_prop <- K - 1L

  # Shrink kappa: merge rows/cols k and k+1
  kappa_prop <- .shrink_kappa_merge(kappa, k)

  # Shrink psi
  psi_prop <- .shrink_psi_merge(psi, k, psi_mode)

  # --- Log-likelihood ratio (DIFFERENTIAL — only affected edges) ---
  # Merge only changes edges involving blocks k or k+1
  zl <- z[i_idx]; zr <- z[j_idx]
  affected <- (zl == k | zl == (k + 1L) | zr == k | zr == (k + 1L))
  idx_aff <- which(affected)

  if (length(idx_aff) > 0) {
    A_ij <- A[cbind(i_idx, j_idx)]
    ll_old <- .partial_loglik(z, kappa, psi, eta, i_idx, j_idx, N_edge, A_ij,
                              idx_aff, psi_mode)
    ll_new <- .partial_loglik(z_prop, kappa_prop, psi_prop, eta, i_idx, j_idx, N_edge, A_ij,
                              idx_aff, psi_mode)
  } else {
    ll_old <- 0; ll_new <- 0
  }

  # --- Prior ratio ---
  n_k_old <- tabulate(z, nbins = K)
  n_k_new <- tabulate(z_prop, nbins = K_prop)
  prior_type <- if (exists(".normalize_partition_prior", mode = "function", inherits = TRUE)) {
    .normalize_partition_prior(partition_prior)
  } else {
    toupper(trimws(partition_prior))
  }
  if (identical(prior_type, "OCRP")) {
    log_prior_ratio <- .ocrp_log_eppf(n_k_new, theta_ocrp) - .ocrp_log_eppf(n_k_old, theta_ocrp)
  } else if (identical(prior_type, "ROCRP")) {
    log_prior_ratio <- .rocrp_log_eppf(n_k_new, theta_ocrp) - .rocrp_log_eppf(n_k_old, theta_ocrp)
  } else {
    log_prior_ratio <- .gn_log_eppf(n_k_new, gamma_gn) - .gn_log_eppf(n_k_old, gamma_gn)
  }

  # --- Proposal ratio: merge is deterministic, split requires calculating
  #     the probability that the reverse split would regenerate this exact partition.
  #     For simplicity we use a symmetric approximation: the proposal ratio
  #     accounts for the split probability of choosing block k and reproducing
  #     the original assignment via restricted Gibbs.
  log_split_prob <- .restricted_gibbs_logprob(
    nodes_merged, z, k, kappa, psi, eta,
    A, i_idx, j_idx, N_edge, edge_by_node,
    psi_mode, n_restricted_scans
  )

  # P(choose merge | K) / P(choose split | K-1) adjustment
  # From K: merge probability = 0.5 (if both available), pick 1/(K-1) pairs
  # From K-1: split probability = 0.5, pick 1/(# blocks with >=2 nodes)
  n_splittable_new <- sum(n_k_new >= 2L)
  if (n_splittable_new == 0) n_splittable_new <- 1  # guard
  log_proposal_ratio <- log(K - 1L) - log(n_splittable_new) + log_split_prob

  log_alpha <- ll_new - ll_old + log_prior_ratio + log_proposal_ratio

  if (log(runif(1)) < log_alpha) {
    return(list(z = z_prop, kappa = kappa_prop, psi = psi_prop,
                accepted = TRUE, move = "merge"))
  }
  list(z = z, kappa = kappa, psi = psi, accepted = FALSE, move = "merge")
}


# --- SPLIT: split block k into k and k+1, shifting labels ---
.split_move <- function(z, kappa, psi, eta, A, i_idx, j_idx, N_edge, edge_by_node,
                        a_kappa, b_kappa, gamma_gn, psi_mode, hyper_psi,
                        n_restricted_scans,
                        partition_prior = "OCRP", theta_ocrp = 1.0) {
  K <- if (psi_mode == "pair") nrow(psi) else length(psi) + 1L
  n <- length(z)
  n_k <- tabulate(z, nbins = K)

  # Pick a block with >= 2 nodes
  splittable <- which(n_k >= 2L)
  if (!length(splittable)) {
    return(list(z = z, kappa = kappa, psi = psi, accepted = FALSE, move = "split"))
  }
  k <- splittable[sample.int(length(splittable), 1L)]
  nodes_in_k <- which(z == k)

  # --- Grow kappa and psi FIRST so restricted Gibbs uses correct dimensions ---
  K_prop <- K + 1L
  kappa_prop <- .grow_kappa_split(kappa, k, a_kappa, b_kappa)
  psi_prop <- .grow_psi_split(psi, k, psi_mode, hyper_psi)

  # --- Propose a binary split via restricted Gibbs ---
  # Randomly initialize: half to k, half to k+1
  n_in_k <- length(nodes_in_k)
  init_assign <- sample(c(0L, 1L), n_in_k, replace = TRUE)
  # Ensure both sides are non-empty
  if (all(init_assign == 0L)) init_assign[sample.int(n_in_k, 1L)] <- 1L
  if (all(init_assign == 1L)) init_assign[sample.int(n_in_k, 1L)] <- 0L

  # Build proposed z: nodes with assign=0 stay in k, assign=1 go to k+1
  z_prop <- z
  # First shift all labels > k up by 1
  z_prop[z > k] <- z_prop[z > k] + 1L
  # Assign the split
  z_prop[nodes_in_k[init_assign == 1L]] <- k + 1L

  # Do restricted Gibbs sweeps to refine the split (using GROWN psi_prop)
  for (sweep in seq_len(n_restricted_scans)) {
    for (idx in sample(seq_along(nodes_in_k))) {
      node_i <- nodes_in_k[idx]
      cur_label <- z_prop[node_i]

      # Score assigning node_i to k vs k+1
      lp <- numeric(2)
      for (opt in 0:1) {
        z_prop[node_i] <- k + opt

        # Quick score: directional contribution only for edges of this node
        lp[opt + 1L] <- .node_dir_score(node_i, z_prop, psi_prop, psi_mode,
                                         A, i_idx, j_idx, N_edge, edge_by_node)
        # Add CRP-like size preference (avoid copy with z[-node_i])
        n_opt <- sum(z_prop == (k + opt)) - 1L  # subtract self
        lp[opt + 1L] <- lp[opt + 1L] + log(max(n_opt, 0.5))
      }

      p <- exp(lp - max(lp))
      chosen <- sample(0:1, 1L, prob = p)
      z_prop[node_i] <- k + chosen
    }
  }

  # Ensure both sub-blocks are non-empty after restricted Gibbs
  n_stay <- sum(z_prop[nodes_in_k] == k)
  n_move <- sum(z_prop[nodes_in_k] == (k + 1L))
  if (n_stay == 0 || n_move == 0) {
    return(list(z = z, kappa = kappa, psi = psi, accepted = FALSE, move = "split"))
  }

  # --- Acceptance ratio (DIFFERENTIAL — only edges involving block k or k+1) ---
  zl_old <- z[i_idx]; zr_old <- z[j_idx]
  # In old z, only block k matters (nodes_in_k were all in block k)
  # In new z_prop, blocks k and k+1 are affected, plus shifted labels
  zl_new <- z_prop[i_idx]; zr_new <- z_prop[j_idx]
  affected_old <- (zl_old == k | zr_old == k)
  affected_new <- (zl_new == k | zl_new == (k + 1L) | zr_new == k | zr_new == (k + 1L))
  affected <- affected_old | affected_new
  # Also include edges where labels shifted (> k in old z)
  shifted <- (zl_old > k | zr_old > k)
  affected <- affected | shifted
  idx_aff <- which(affected)

  A_ij <- A[cbind(i_idx, j_idx)]
  ll_old <- .partial_loglik(z, kappa, psi, eta, i_idx, j_idx, N_edge, A_ij,
                            idx_aff, psi_mode)
  ll_new <- .partial_loglik(z_prop, kappa_prop, psi_prop, eta, i_idx, j_idx, N_edge, A_ij,
                            idx_aff, psi_mode)

  n_k_old <- tabulate(z, nbins = K)
  n_k_new <- tabulate(z_prop, nbins = K_prop)
  prior_type <- if (exists(".normalize_partition_prior", mode = "function", inherits = TRUE)) {
    .normalize_partition_prior(partition_prior)
  } else {
    toupper(trimws(partition_prior))
  }
  if (identical(prior_type, "OCRP")) {
    log_prior_ratio <- .ocrp_log_eppf(n_k_new, theta_ocrp) - .ocrp_log_eppf(n_k_old, theta_ocrp)
  } else if (identical(prior_type, "ROCRP")) {
    log_prior_ratio <- .rocrp_log_eppf(n_k_new, theta_ocrp) - .rocrp_log_eppf(n_k_old, theta_ocrp)
  } else {
    log_prior_ratio <- .gn_log_eppf(n_k_new, gamma_gn) - .gn_log_eppf(n_k_old, gamma_gn)
  }

  # Proposal ratio (reverse of split is a merge)
  log_split_prob <- .restricted_gibbs_logprob(
    nodes_in_k, z_prop, k, kappa_prop, psi_prop, eta,
    A, i_idx, j_idx, N_edge, edge_by_node,
    psi_mode, n_restricted_scans
  )

  n_mergeable_new <- K_prop - 1L  # number of adjacent pairs
  n_splittable_old <- length(splittable)
  log_proposal_ratio <- log(n_splittable_old) - log(n_mergeable_new) - log_split_prob

  log_alpha <- ll_new - ll_old + log_prior_ratio + log_proposal_ratio

  if (log(runif(1)) < log_alpha) {
    return(list(z = z_prop, kappa = kappa_prop, psi = psi_prop,
                accepted = TRUE, move = "split"))
  }
  list(z = z, kappa = kappa, psi = psi, accepted = FALSE, move = "split")
}


# =====================================================================
# Helper functions for split-merge
# =====================================================================

# Partial log-likelihood over a subset of edges (vectorized)
# Used for differential merge/split acceptance ratios
.partial_loglik <- function(z, kappa, psi, eta, i_idx, j_idx, N_edge, A_ij,
                            edge_subset, psi_mode) {
  ii <- i_idx[edge_subset]; jj <- j_idx[edge_subset]
  zi <- z[ii]; zj <- z[jj]
  n_e <- N_edge[edge_subset]; a_e <- A_ij[edge_subset]

  k_lo <- pmin(zi, zj); k_hi <- pmax(zi, zj)
  kap_vals <- kappa[cbind(k_lo, k_hi)]
  lam <- eta[ii] * eta[jj] * kap_vals

  ll_vol <- ifelse(lam > 0, n_e * log(lam) - lam,
                   ifelse(n_e > 0, -1e10, 0))

  same <- (zi == zj)
  ll_dir <- ifelse(same, -n_e * log(2), 0)

  diff_mask <- !same & (n_e > 0)
  if (any(diff_mask)) {
    zi_d <- zi[diff_mask]; zj_d <- zj[diff_mask]
    a_d  <- a_e[diff_mask]; n_d  <- n_e[diff_mask]
    if (psi_mode == "pair") {
      klo <- pmin(zi_d, zj_d); khi <- pmax(zi_d, zj_d)
      psi_val <- psi[cbind(klo, khi)]
    } else {
      d <- abs(zi_d - zj_d)
      d <- pmin(d, length(psi)); d <- pmax(d, 1L)
      psi_val <- psi[d]
    }
    sgn_d <- sign(zj_d - zi_d)
    phi <- sgn_d * psi_val
    ll_dir[diff_mask] <- a_d * phi - n_d * .log1pexp_safe(phi)
  }
  sum(ll_vol) + sum(ll_dir)
}

# Full directional + volume log-likelihood (up to constants that cancel)
# VECTORIZED for speed — no R-level for-loop over edges
.full_loglik <- function(z, kappa, psi, eta, i_idx, j_idx, N_edge, A, psi_mode) {
  A_ij <- A[cbind(i_idx, j_idx)]
  zi <- z[i_idx]; zj <- z[j_idx]
  n_e <- N_edge; a_e <- A_ij

  # Volume (Poisson): lam = eta_i * eta_j * kappa[min(zi,zj), max(zi,zj)]
  k_lo <- pmin(zi, zj); k_hi <- pmax(zi, zj)
  kap_vals <- kappa[cbind(k_lo, k_hi)]
  lam <- eta[i_idx] * eta[j_idx] * kap_vals

  # log-lik for volume part
  ll_vol <- ifelse(lam > 0, n_e * log(lam) - lam,
                   ifelse(n_e > 0, -1e10, 0))

  # Direction part, vectorized
  same <- (zi == zj)
  # within-block: -N * log(2)
  ll_dir <- ifelse(same, -n_e * log(2), 0)

  diff_mask <- !same & (n_e > 0)
  if (any(diff_mask)) {
    zi_d <- zi[diff_mask]; zj_d <- zj[diff_mask]
    a_d  <- a_e[diff_mask]; n_d <- n_e[diff_mask]

    if (psi_mode == "pair") {
      klo <- pmin(zi_d, zj_d); khi <- pmax(zi_d, zj_d)
      psi_val <- psi[cbind(klo, khi)]
    } else {
      d <- abs(zi_d - zj_d)
      d <- pmin(d, length(psi))  # safety
      d <- pmax(d, 1L)
      psi_val <- psi[d]
    }
    sgn_d <- sign(zj_d - zi_d)
    phi <- sgn_d * psi_val
    ll_dir[diff_mask] <- a_d * phi - n_d * .log1pexp_safe(phi)
  }

  sum(ll_vol) + sum(ll_dir)
}

# Node-level directional score (sum over edges of node i)
# VECTORIZED — no R-level for-loop over incident edges
.node_dir_score <- function(node_i, z, psi, psi_mode,
                            A, i_idx, j_idx, N_edge, edge_by_node) {
  e_list <- edge_by_node[[node_i]]
  if (!length(e_list)) return(0)

  ii <- i_idx[e_list]; jj <- j_idx[e_list]
  zi <- z[ii]; zj <- z[jj]
  A_vals <- A[cbind(ii, jj)]
  n_e <- N_edge[e_list]

  same <- (zi == zj)
  score <- sum(ifelse(same, -n_e * log(2), 0))

  diff_mask <- !same & (n_e > 0)
  if (any(diff_mask)) {
    zi_d <- zi[diff_mask]; zj_d <- zj[diff_mask]
    a_d  <- A_vals[diff_mask]; n_d <- n_e[diff_mask]

    if (psi_mode == "pair") {
      klo <- pmin(zi_d, zj_d); khi <- pmax(zi_d, zj_d)
      psi_val <- psi[cbind(klo, khi)]
    } else {
      d <- abs(zi_d - zj_d)
      d <- pmin(d, length(psi))
      d <- pmax(d, 1L)
      psi_val <- psi[d]
    }
    sgn_d <- sign(zj_d - zi_d)
    phi <- sgn_d * psi_val
    score <- score + sum(a_d * phi - n_d * .log1pexp_safe(phi))
  }
  score
}

# GN log-EPPF (unnormalized, up to constants that depend on n but not on the partition)
.gn_log_eppf <- function(n_k, gamma_gn) {
  n_k <- n_k[n_k > 0]
  K <- length(n_k)
  if (K == 0) return(-Inf)
  # log EPPF ∝ log[(1-γ)_{K-1}] + Σ log(n_k!) + K² terms
  # Using the form from the OSBM paper Eq (gn_eppf):
  # p(n_1,...,n_K) ∝ γ(1-γ)_{K-1} * ∏(n_k!) / [(n-1)!(1+γ)_{n-1}]
  # The denominator depends only on n, so it cancels in ratios.
  # Numerator (log): log γ + Σ_{j=0}^{K-2} log(1-γ+j) + Σ_k log(n_k!)
  val <- log(gamma_gn)
  if (K >= 2L) {
    val <- val + sum(log((1 - gamma_gn) + 0:(K - 2L)))
  }
  val <- val + sum(lfactorial(n_k))
  val
}

# OCRP log-EPPF (unnormalized, up to constants that depend on n but not on the partition)
# p(n_1,...,n_K) = theta^K / (theta)_{(n)} * n! / (N_1 ... N_K)
# where N_j = n_j + ... + n_K (right cumulative sums).
# The denominator (theta)_{(n)} depends only on n, so cancels in ratios.
.rocrp_log_eppf <- function(n_k, theta_ocrp) {
  n_k <- n_k[n_k > 0]
  K <- length(n_k)
  if (K == 0) return(-Inf)
  n <- sum(n_k)
  # Right cumulative sums
  N_j <- rev(cumsum(rev(n_k)))
  # log p = K*log(theta) + log(n!) - sum(log(N_j)) - log((theta)_{(n)})
  # Since (theta)_{(n)} = theta*(theta+1)*...*(theta+n-1) depends only on n,
  # it cancels in ratios. Keep it for absolute value:
  val <- K * log(theta_ocrp) + lfactorial(n) - sum(log(N_j))
  val
}

.ocrp_log_eppf <- function(n_k, theta_ocrp) {
  n_k <- n_k[n_k > 0]
  K <- length(n_k)
  if (K == 0) return(-Inf)
  n <- sum(n_k)
  S_j <- cumsum(n_k)
  val <- K * log(theta_ocrp) + lfactorial(n) - sum(log(S_j))
  val
}

# Restricted Gibbs log-probability of producing a specific binary split
.restricted_gibbs_logprob <- function(nodes, z_target, k, kappa, psi, eta,
                                      A, i_idx, j_idx, N_edge, edge_by_node,
                                      psi_mode, n_sweeps) {
  # Approximate: compute log probability of the last sweep reproducing z_target
  lp_total <- 0

  for (idx in seq_along(nodes)) {
    node_i <- nodes[idx]
    target_label <- z_target[node_i]

    lp <- numeric(2)
    for (opt in 0:1) {
      z_target[node_i] <- k + opt

      lp[opt + 1L] <- .node_dir_score(node_i, z_target, psi, psi_mode,
                                        A, i_idx, j_idx, N_edge, edge_by_node)
      n_opt <- sum(z_target == (k + opt)) - 1L
      lp[opt + 1L] <- lp[opt + 1L] + log(max(n_opt, 0.5))
    }

    # Restore the target label
    z_target[node_i] <- target_label

    log_norm <- .log_sum_exp(lp)
    chosen_idx <- if (target_label == k) 1L else 2L
    lp_total <- lp_total + (lp[chosen_idx] - log_norm)
  }
  lp_total
}

.log_sum_exp <- function(x) {
  m <- max(x[is.finite(x)])
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(x - m)))
}

# Shrink kappa by merging rows/cols k and k+1 (average the rates)
.shrink_kappa_merge <- function(kappa, k) {
  K <- nrow(kappa)
  K_new <- K - 1L
  if (K_new < 1L) return(matrix(0, 1, 1))

  # Map: old labels {1,..,k, k+1,..,K} -> {1,..,k,..,K-1}
  # Blocks k and k+1 merge into k; blocks > k+1 shift down
  kap_new <- matrix(0, K_new, K_new)

  map_old2new <- function(j) {
    if (j <= k) j
    else if (j == k + 1L) k  # merges into k
    else j - 1L
  }

  for (a in seq_len(K)) {
    for (b in a:K) {
      a_new <- map_old2new(a)
      b_new <- map_old2new(b)
      p <- min(a_new, b_new); q <- max(a_new, b_new)
      # For merged pairs, average the rates (crude but simple)
      kap_new[p, q] <- kap_new[p, q] + kappa[a, b]
    }
  }
  # Normalize: count how many old pairs mapped to each new pair
  counts <- matrix(0, K_new, K_new)
  for (a in seq_len(K)) {
    for (b in a:K) {
      a_new <- map_old2new(a)
      b_new <- map_old2new(b)
      p <- min(a_new, b_new); q <- max(a_new, b_new)
      counts[p, q] <- counts[p, q] + 1
    }
  }
  counts[counts == 0] <- 1
  kap_new <- kap_new / counts

  # Symmetrize
  kap_new[lower.tri(kap_new)] <- t(kap_new)[lower.tri(kap_new)]
  kap_new
}

# Shrink psi for merge (WST: remove row/col k+1; SST: average adjacent distances)
.shrink_psi_merge <- function(psi, k, psi_mode) {
  if (psi_mode == "pair") {
    K <- nrow(psi)
    keep <- setdiff(seq_len(K), k + 1L)
    psi_new <- psi[keep, keep, drop = FALSE]
    diag(psi_new) <- 0
    return(psi_new)
  }

  # SST: psi is length K-1 distance vector. After merging blocks k and k+1,
  # the max distance decreases by 1.
  D_old <- length(psi)
  if (D_old <= 1L) return(numeric(0))
  # Remove one element. Since pairs at old distance d map to:
  # - same distance d if both endpoints are outside {k, k+1}
  # - potentially d-1 if one endpoint was in {k+1} and shifted down
  # Simplest correct approach: just drop the last distance
  psi_new <- psi[seq_len(D_old - 1L)]
  psi_new
}

# Grow kappa for split (duplicate row/col k, draw new cross-block rate)
.grow_kappa_split <- function(kappa, k, a_kappa, b_kappa) {
  K <- nrow(kappa)
  K_new <- K + 1L

  kap_new <- matrix(0, K_new, K_new)

  # Shift labels > k up by 1
  map_old2new <- function(j) if (j <= k) j else j + 1L

  for (a in seq_len(K)) {
    for (b in seq_len(K)) {
      kap_new[map_old2new(a), map_old2new(b)] <- kappa[a, b]
    }
  }

  # New block k+1 inherits from old block k, with noise
  for (j in seq_len(K_new)) {
    if (j == k + 1L) next
    # Cross-rate between new block k+1 and block j: start from kappa[k, map_back(j)]
    val <- kap_new[k, j]  # inherited from old block k
    kap_new[k + 1L, j] <- kap_new[j, k + 1L] <- val
  }
  # Within new block k+1: draw from prior
  kap_new[k + 1L, k + 1L] <- rgamma(1, shape = a_kappa, rate = b_kappa)
  # Cross between k and k+1: draw from prior
  kap_new[k, k + 1L] <- kap_new[k + 1L, k] <- rgamma(1, shape = a_kappa, rate = b_kappa)

  kap_new
}

# Grow psi for split (WST: duplicate row/col; SST: extend distance vector)
.grow_psi_split <- function(psi, k, psi_mode, hyper_psi) {
  if (psi_mode == "pair") {
    K <- nrow(psi)
    K_new <- K + 1L

    psi_new <- matrix(0, K_new, K_new)
    map_old2new <- function(j) if (j <= k) j else j + 1L

    for (a in seq_len(K)) {
      for (b in seq_len(K)) {
        psi_new[map_old2new(a), map_old2new(b)] <- psi[a, b]
      }
    }

    # New block k+1 inherits psi from old k
    for (j in seq_len(K_new)) {
      if (j == k + 1L || j == k) next
      psi_new[min(k + 1L, j), max(k + 1L, j)] <-
        psi_new[min(k, j), max(k, j)]
    }
    # psi between k and k+1: small value (close blocks)
    mu0 <- hyper_psi$mu0
    sig0 <- hyper_psi$sigma0
    psi_new[k, k + 1L] <- abs(rnorm(1, mean = mu0 * 0.5, sd = sig0 * 0.5))
    psi_new[k + 1L, k] <- psi_new[k, k + 1L]
    diag(psi_new) <- 0

    return(psi_new)
  }

  # SST: extend distance vector by 1
  D_old <- length(psi)
  tau0 <- hyper_psi$tau0
  # Add a new distance: last_psi + delta, delta ~ TN+(0, tau0^2)
  last_val <- if (D_old > 0) psi[D_old] else 0
  delta <- abs(rnorm(1, mean = 0, sd = tau0))
  psi_new <- c(psi, last_val + delta)
  psi_new
}
