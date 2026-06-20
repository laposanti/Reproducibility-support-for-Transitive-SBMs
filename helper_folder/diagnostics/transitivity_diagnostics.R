# Order-conformity and transitivity diagnostics shared by analysis and post-processing scripts.
#
# ---------------- New helpers: diagnostics (robust, extended) -----------------
order_from_rho <- function(rho, method = c("mean","bt","identity"),
                           order_direction = c("strong_to_weak","weak_to_strong"),
                           eps = 1e-8) {
  method <- match.arg(method)
  order_direction <- match.arg(order_direction)
  K <- nrow(rho)
  stopifnot(ncol(rho) == K)
  if (K <= 1) return(list(order = seq_len(K), scores = rep(0, K), direction = order_direction))

  if (method == "identity") {
    return(list(order = seq_len(K), scores = seq_len(K), direction = order_direction))
  }
  
  # clamp to (eps, 1-eps) & set diag to 0.5
  R <- pmin(pmax(rho, eps), 1 - eps)
  diag(R) <- 0.5
  
  if (method == "mean") {
    r <- rowMeans(R)
  } else {
    # Bradley–Terry via least squares with regularization fallback
    logit <- function(p) log(p/(1 - p))
    L <- matrix(0, K, K); b <- numeric(K)
    for (k in 1:(K - 1)) for (l in (k + 1):K) {
      y <- logit(R[k, l])                    # ≈ β_k - β_l
      L[k, k] <- L[k, k] + 1
      L[l, l] <- L[l, l] + 1
      L[k, l] <- L[k, l] - 1
      L[l, k] <- L[l, k] - 1
      b[k]    <- b[k] + y
      b[l]    <- b[l] - y
    }
    A <- rbind(L, rep(1, K)); y <- c(b, 0)
    beta <- tryCatch(
      as.numeric(qr.solve(A, y)),
      error = function(e) {
        ridge <- 1e-6
        as.numeric(qr.solve(rbind(L + ridge*diag(K), rep(1, K)), y))
      }
    )
    r <- beta
  }
  
  # ---- enforce direction: block 1 = strongest ----
  if (order_direction == "strong_to_weak") {
    ord <- order(-r, seq_len(K))  # highest score first → index 1 is strongest
  } else {
    ord <- order( r, seq_len(K))  # lowest first → index 1 is weakest
  }
  list(order = ord, scores = r, direction = order_direction)
}

.ordered_rho <- function(rho,
                         method_order = c("mean","bt","identity"),
                         order_direction = c("strong_to_weak","weak_to_strong"),
                         eps = 1e-8) {
  method_order <- match.arg(method_order)
  order_direction <- match.arg(order_direction)
  if (method_order == "identity") {
    R <- as.matrix(rho)
    diag(R) <- 0.5
    return(R)
  }
  ord <- order_from_rho(
    rho, method = method_order,
    order_direction = order_direction, eps = eps
  )$order
  R <- rho[ord, ord, drop = FALSE]
  diag(R) <- 0.5
  R
}

.mean_forward_rho <- function(rho,
                              method_order = c("mean","bt","identity"),
                              order_direction = c("strong_to_weak","weak_to_strong"),
                              eps = 1e-8) {
  K <- nrow(rho)
  if (K < 2L) return(NA_real_)
  R <- .ordered_rho(rho, method_order = method_order,
                    order_direction = order_direction, eps = eps)
  mean(R[upper.tri(R)], na.rm = TRUE)
}

sigm <- function(x) 1/(1+exp(-x))

# --- Safe compaction: ensure labels are 1..Kc and handle NA gracefully ---
.compact_labels <- function(z) {
  if (anyNA(z)) {
    keep <- which(!is.na(z))
    z <- z[keep]
    if (length(z) == 0L) {
      return(list(zc = integer(0), occ = integer(0), Kc = 0L, keep = integer(0)))
    }
    occ <- sort(unique(z))
    return(list(zc = match(z, occ), occ = occ, Kc = length(occ), keep = keep))
  } else {
    occ <- sort(unique(z))
    return(list(zc = match(z, occ), occ = occ, Kc = length(occ), keep = seq_along(z)))
  }
}

.logit_safe <- function(p, eps = 1e-8) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

.diag_from_rho_invariant <- function(rho, eps = 1e-8) {
  K <- nrow(rho)
  if (K < 2L) {
    return(list(
      transitive_triads = NA_real_,
      cycle_mass_weighted = NA_real_,
      min_backward_weight = NA_real_,
      min_backward_weight_norm = NA_real_,
      hierarchy_energy = NA_real_,
      curl_energy = NA_real_
    ))
  }

  R <- pmin(pmax(rho, eps), 1 - eps)
  diag(R) <- 0.5
  B <- (R > 0.5)
  diag(B) <- FALSE

  w <- abs(.logit_safe(R, eps = eps))
  diag(w) <- 0
  w_total <- sum(w[upper.tri(w)])

  # Approximate minimum feedback arc set using score sort + local adjacent swaps.
  score <- rowSums(B) - colSums(B)
  ord <- order(-score, seq_len(K))
  perm_sign <- 1L
  improved <- TRUE
  passes <- 0L
  while (improved && passes < 5L) {
    improved <- FALSE
    passes <- passes + 1L
    for (i in seq_len(K - 1L)) {
      o1 <- ord
      o2 <- ord
      o2[c(i, i + 1L)] <- o2[c(i + 1L, i)]

      backward_mass <- function(o) {
        v <- 0
        for (a in seq_len(K - 1L)) for (b in (a + 1L):K) {
          if (B[o[b], o[a]]) v <- v + w[o[b], o[a]]
        }
        v
      }

      m1 <- backward_mass(o1)
      m2 <- backward_mass(o2)
      if (m2 + 1e-12 < m1) {
        ord <- o2
        perm_sign <- -perm_sign
        improved <- TRUE
      }
    }
  }

  min_backward <- 0
  for (a in seq_len(K - 1L)) for (b in (a + 1L):K) {
    if (B[ord[b], ord[a]]) min_backward <- min_backward + w[ord[b], ord[a]]
  }
  min_backward_norm <- if (w_total > 0) min_backward / w_total else NA_real_

  # Bradley-Terry / Hodge hierarchy energy on antisymmetric log-odds.
  L <- .logit_safe(R, eps = eps)
  diag(L) <- 0
  L <- 0.5 * (L - t(L))

  Lap <- matrix(0, K, K)
  b <- numeric(K)
  for (i in seq_len(K - 1L)) for (j in (i + 1L):K) {
    wij <- w[i, j]
    Lap[i, i] <- Lap[i, i] + wij
    Lap[j, j] <- Lap[j, j] + wij
    Lap[i, j] <- Lap[i, j] - wij
    Lap[j, i] <- Lap[j, i] - wij
    b[i] <- b[i] + wij * L[i, j]
    b[j] <- b[j] - wij * L[i, j]
  }
  A_sys <- rbind(cbind(Lap, rep(1, K)), c(rep(1, K), 0))
  rhs <- c(b, 0)
  sol <- tryCatch(qr.solve(A_sys, rhs), error = function(e) rep(0, K + 1L))
  s_hat <- sol[seq_len(K)]

  num <- 0
  den <- 0
  for (i in seq_len(K - 1L)) for (j in (i + 1L):K) {
    wij <- w[i, j]
    r_ij <- L[i, j] - (s_hat[i] - s_hat[j])
    num <- num + wij * r_ij^2
    den <- den + wij * L[i, j]^2
  }
  hierarchy_energy <- if (den > 0) 1 - num / den else NA_real_

  # Triad diagnostics
  tri_total <- choose(K, 3)
  cyc_count <- 0L
  cyc_mass <- 0
  curl_acc <- 0
  if (tri_total > 0) {
    for (a in seq_len(K - 2L)) for (b in (a + 1L):(K - 1L)) for (c in (b + 1L):K) {
      cyc_abc <- as.integer(B[a, b] && B[b, c] && B[c, a])
      cyc_acb <- as.integer(B[a, c] && B[c, b] && B[b, a])
      if ((cyc_abc + cyc_acb) > 0L) {
        cyc_count <- cyc_count + 1L
        tri_w <- w[a, b] + w[b, c] + w[c, a]
        cyc_mass <- cyc_mass + tri_w
      }
      curl <- L[a, b] + L[b, c] + L[c, a]
      curl_acc <- curl_acc + curl^2
    }
  }

  transitive_triads <- if (tri_total > 0) 1 - cyc_count / tri_total else NA_real_
  cycle_mass_weighted <- if (tri_total > 0) cyc_mass / tri_total else NA_real_
  curl_energy <- if (tri_total > 0) curl_acc / tri_total else NA_real_

  list(
    transitive_triads = transitive_triads,
    cycle_mass_weighted = cycle_mass_weighted,
    min_backward_weight = min_backward,
    min_backward_weight_norm = min_backward_norm,
    hierarchy_energy = hierarchy_energy,
    curl_energy = curl_energy
  )
}

.check_wst_sst_sorted <- function(rho, tol = 1e-12) {
  K <- nrow(rho)
  if (K < 2L) return(list(wst = NA, sst = NA))

  R <- as.matrix(rho)
  diag(R) <- 0.5
  ord <- order(-rowMeans(R), seq_len(K))
  Rs <- R[ord, ord, drop = FALSE]

  wst_ok <- all(Rs[upper.tri(Rs)] + tol >= 0.5)
  if (K < 3L) return(list(wst = wst_ok, sst = NA))

  sst_ok <- wst_ok
  for (a in seq_len(K - 2L)) {
    for (b in (a + 1L):(K - 1L)) {
      for (c in (b + 1L):K) {
        rac <- Rs[a, c]
        rab <- Rs[a, b]
        rbc <- Rs[b, c]
        if (!(rac + tol >= max(rab, rbc))) sst_ok <- FALSE
        if (!sst_ok) return(list(wst = wst_ok, sst = FALSE))
      }
    }
  }
  list(wst = wst_ok, sst = sst_ok)
}

.prior_constraint_prob_mc <- function(K, n_mc = 2000L, a_lambda = 1, b_lambda = 1,
                                      seed = 12345L) {
  if (K < 3L || n_mc < 1L) {
    return(list(p_wst = NA_real_, p_sst = NA_real_, n_mc = as.integer(n_mc)))
  }
  if (!is.null(seed)) set.seed(seed)

  wst_hits <- logical(n_mc)
  sst_hits <- logical(n_mc)
  for (t in seq_len(n_mc)) {
    lam <- matrix(stats::rgamma(K * K, shape = a_lambda, rate = b_lambda), nrow = K, ncol = K)
    rho <- lam / (lam + t(lam) + 1e-12)
    diag(rho) <- 0.5
    chk <- .check_wst_sst_sorted(rho)
    wst_hits[t] <- isTRUE(chk$wst)
    sst_hits[t] <- isTRUE(chk$sst)
  }

  list(
    p_wst = mean(wst_hits),
    p_sst = mean(sst_hits),
    n_mc = as.integer(n_mc)
  )
}

# Build a draw-specific strong->weak label ranking from empirical block rates,
# then compute hierarchy-violation mass in a label-invariant way.
.rank_partition_by_strength <- function(A, z_vec, alpha = 0.5,
                                        order_direction = c("strong_to_weak","weak_to_strong")) {
  order_direction <- match.arg(order_direction)
  valid <- !is.na(z_vec)
  if (sum(valid) < 2L) return(rep(NA_integer_, length(z_vec)))

  A_use <- A[valid, valid, drop = FALSE]
  z_use <- z_vec[valid]
  comp <- .compact_labels(z_use)
  zc <- comp$zc
  Kc <- comp$Kc
  if (Kc < 1L) return(rep(NA_integer_, length(z_vec)))

  C_blk <- matrix(0, Kc, Kc)
  M_blk <- matrix(0, Kc, Kc)
  n_use <- nrow(A_use)
  for (i in seq_len(n_use)) {
    for (j in seq_len(n_use)) if (i != j) {
      ki <- zc[i]; kj <- zc[j]
      C_blk[ki, kj] <- C_blk[ki, kj] + A_use[i, j]
      M_blk[ki, kj] <- M_blk[ki, kj] + A_use[i, j] + A_use[j, i]
    }
  }

  P_blk <- (C_blk + alpha) / (pmax(M_blk, 0) + 2 * alpha)
  diag(P_blk) <- 0.5
  blk_scores <- rowMeans(P_blk)
  if (order_direction == "strong_to_weak") {
    blk_ord <- order(-blk_scores, seq_len(Kc))
  } else {
    blk_ord <- order(blk_scores, seq_len(Kc))
  }

  # rank map on compact labels: 1 = strongest (or weakest if configured)
  rank_map <- integer(Kc)
  rank_map[blk_ord] <- seq_len(Kc)
  z_rank <- rep(NA_integer_, length(z_vec))
  z_rank[valid] <- rank_map[zc]
  z_rank
}

.violation_stats_from_ranked <- function(A, z_ranked) {
  rows <- row(A); cols <- col(A); vals <- A
  z_i <- z_ranked[rows]; z_j <- z_ranked[cols]
  cross_mask <- (z_i != z_j) & (vals > 0) & !is.na(z_i) & !is.na(z_j)
  if (!any(cross_mask, na.rm = TRUE)) {
    return(list(rate = NA_real_, count = NA_real_, cross_mass = NA_real_))
  }
  # ranked labels: smaller index means stronger
  viol_mask <- cross_mask & (z_i > z_j)
  v_mass <- sum(vals[viol_mask], na.rm = TRUE)
  t_mass <- sum(vals[cross_mask], na.rm = TRUE)
  list(
    rate = if (t_mass > 0) v_mass / t_mass else NA_real_,
    count = v_mass,
    cross_mass = t_mass
  )
}

.empirical_block_rates_one <- function(A, z_vec,
                                       alpha = 0.5,
                                       T_block = 1000L,
                                       seed = 123,
                                       method_order = c("mean","bt","identity"),
                                       order_direction = c("strong_to_weak","weak_to_strong")) {
  method_order    <- match.arg(method_order)
  order_direction <- match.arg(order_direction)
  
  ## keep only non-NA nodes
  valid <- !is.na(z_vec)
  if (sum(valid) < 3L) {
    return(list(
      thetaW     = NA_real_,
      thetaS     = NA_real_,
      prem       = 0L,
      thetaW_all = NA_real_,
      thetaS_all = NA_real_,
      coverage   = 0,
      total      = 0L
    ))
  }
  
  A_use <- A[valid, valid, drop = FALSE]
  z_use <- z_vec[valid]
  n_use <- nrow(A_use)
  if (n_use < 3L) {
    return(list(
      thetaW     = NA_real_,
      thetaS     = NA_real_,
      prem       = 0L,
      thetaW_all = NA_real_,
      thetaS_all = NA_real_,
      coverage   = 0,
      total      = 0L
    ))
  }
  
  ## compact labels {1,...,Kc}
  comp <- .compact_labels(z_use)
  zc   <- comp$zc
  Kc   <- comp$Kc
  
  if (Kc < 3L) {
    return(list(
      thetaW     = NA_real_,
      thetaS     = NA_real_,
      prem       = 0L,
      thetaW_all = NA_real_,
      thetaS_all = NA_real_,
      coverage   = 0,
      total      = 0L
    ))
  }
  
  ## block-level counts C_blk(k,l) = sum_{i in k, j in l} A_ij
  ## M_blk(k,l) = total volume between blocks k,l (A_ij + A_ji)
  C_blk <- matrix(0, Kc, Kc)
  M_blk <- matrix(0, Kc, Kc)
  
  for (i in seq_len(n_use)) {
    for (j in seq_len(n_use)) if (i != j) {
      ki <- zc[i]; kj <- zc[j]
      C_blk[ki, kj] <- C_blk[ki, kj] + A_use[i, j]
      M_blk[ki, kj] <- M_blk[ki, kj] + A_use[i, j] + A_use[j, i]
    }
  }
  
  rhohat <- (C_blk + alpha) / pmax(M_blk + 2 * alpha, .Machine$double.eps)
  diag(rhohat) <- 0.5
  
  block_diag_rates_ext(
    rhohat,
    T_block        = T_block,
    seed           = seed,
    method_order   = method_order,
    order_direction = order_direction
  )
}

# -- small utility: coerce possibly transposed traces --------------------
.coerce_Sxn <- function(M, n_expected, name = "trace") {
  if (!is.matrix(M)) stop(name, " must be a matrix.")
  if (nrow(M) == n_expected) {
    return(t(M))   # n x S  -> S x n
  } else if (ncol(M) == n_expected) {
    return(M)      # S x n  -> ok
  } else {
    stop(sprintf("%s has incompatible dimensions: got %d x %d, expected one dim == %d",
                 name, nrow(M), ncol(M), n_expected))
  }
}


#' Convert OSBM parameters (psi) to Success Probability Matrix (rho)
#' Robust to K=1 and Ragged Inputs
rho_from_osbm_draw <- function(psi_draw, regime, K) {
  
  # Case 1: K=1 (Degenerate)
  if (is.null(K) || K < 2) {
    return(matrix(0.5, 1, 1))
  }
  
  rho <- matrix(0.5, K, K)
  
  if (regime == "WST") {
    # psi_draw should be K x K matrix
    if (!is.matrix(psi_draw)) {
      # Fallback for scalar/vector edge cases
      if (length(psi_draw) == K*K) psi_draw <- matrix(psi_draw, K, K)
      else return(rho) # Safety return
    }
    
    # Fill upper triangle
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        # Check bounds safety
        if (i <= nrow(psi_draw) && j <= ncol(psi_draw)) {
          val <- psi_draw[i, j]
          rho[i, j] <- 1 / (1 + exp(-val)) # Logit to Prob (i > j)
          rho[j, i] <- 1 - rho[i, j]
        }
      }
    }
    
  } else { # SST
    # psi_draw is vector of length K-1
    if (length(psi_draw) < (K-1)) return(rho) # Safety
    
    for (i in 1:(K-1)) {
      for (j in (i+1):K) {
        d <- abs(i - j)
        if (d <= length(psi_draw)) {
          val <- psi_draw[d] # Distance parameter
          # SST direction: i < j implies i is stronger (distance d)
          # sgn(j-i) = 1. 
          rho[i, j] <- 1 / (1 + exp(-val))
          rho[j, i] <- 1 - rho[i, j]
        }
      }
    }
  }
  
  return(rho)
}


.empirical_item_rates_one <- function(A, z_vec,
                                      T = 2000L,
                                      alpha = 0.5,
                                      seed = NULL,
                                      method_order = c("mean","bt","identity"),
                                      order_direction = c("strong_to_weak","weak_to_strong")) {
  method_order <- match.arg(method_order)
  order_direction <- match.arg(order_direction)
  
  valid <- !is.na(z_vec)
  if (sum(valid) < 3L) {
    return(list(thetaW = NA_real_, thetaS = NA_real_,
                prem = 0L, coverage = 0,
                thetaW_all = NA_real_, thetaS_all = NA_real_))
  }
  
  A_use <- A[valid, valid, drop = FALSE]
  z_use <- z_vec[valid]
  n     <- nrow(A_use)
  if (n < 3L) {
    return(list(thetaW = NA_real_, thetaS = NA_real_,
                prem = 0L, coverage = 0,
                thetaW_all = NA_real_, thetaS_all = NA_real_))
  }
  
  ## --- 1. Compact labels & block-level empirical order ---
  comp  <- .compact_labels(z_use)
  zc    <- comp$zc        # labels in {1,...,Kc}
  occ   <- comp$occ
  Kc    <- comp$Kc
  if (Kc < 3L) {
    return(list(thetaW = NA_real_, thetaS = NA_real_,
                prem = 0L, coverage = 0,
                thetaW_all = NA_real_, thetaS_all = NA_real_))
  }
  
  # Block-level aggregated counts
  C_blk <- matrix(0, Kc, Kc)
  M_blk <- matrix(0, Kc, Kc)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) if (i != j) {
      ki <- zc[i]; kj <- zc[j]
      C_blk[ki, kj] <- C_blk[ki, kj] + A_use[i, j]
      M_blk[ki, kj] <- M_blk[ki, kj] + A_use[i, j] + A_use[j, i]
    }
  }
  # Smoothed block-level probabilities
  P_blk <- (C_blk + alpha) / (pmax(M_blk, 0) + 2 * alpha)
  diag(P_blk) <- 0.5
  
  # Determine block order. WST/SST use intrinsic labels; DC-SBM can be
  # reordered by an empirical/model-derived strength score at the caller.
  blk_ord <- order_from_rho(
    P_blk, method = method_order,
    order_direction = order_direction
  )$order
  
  # Map each block to its rank
  blk_rank <- integer(Kc)
  blk_rank[blk_ord] <- seq_len(Kc)  # strongest block has rank 1
  
  ## --- 2. Build per-block item lists (by ordered rank) ---
  block_items <- lapply(seq_len(Kc), function(rk) {
    # original block index corresponding to rank rk
    b <- which(blk_rank == rk)
    which(zc == b)
  })
  
  # Keep only ranks with at least 1 node
  nonempty <- which(vapply(block_items, length, integer(1)) > 0L)
  if (length(nonempty) < 3L) {
    return(list(thetaW = NA_real_, thetaS = NA_real_,
                prem = 0L, coverage = 0,
                thetaW_all = NA_real_, thetaS_all = NA_real_))
  }
  
  ## --- 3. Smoothed item-level probabilities ---
  matches <- A_use + t(A_use)
  P_hat   <- (A_use + alpha) / (matches + 2 * alpha)
  diag(P_hat) <- 0.5
  
  ## --- 4. Sample item triples consistent with block order ---
  if (!is.null(seed)) set.seed(seed)
  
  # all block triples in rank space
  blk_triples <- utils::combn(nonempty, 3L)
  n_blk_triples <- ncol(blk_triples)
  if (n_blk_triples == 0L) {
    return(list(thetaW = NA_real_, thetaS = NA_real_,
                prem = 0L, coverage = 0,
                thetaW_all = NA_real_, thetaS_all = NA_real_))
  }
  
  T_use <- min(T, n_blk_triples)
  pick  <- if (n_blk_triples > T_use) {
    sort(sample.int(n_blk_triples, T_use, replace = FALSE))
  } else {
    seq_len(n_blk_triples)
  }
  blk_triples <- blk_triples[, pick, drop = FALSE]
  
  triples_items <- matrix(NA_integer_, nrow = 3L, ncol = ncol(blk_triples))
  for (t in seq_len(ncol(blk_triples))) {
    rk_a <- blk_triples[1, t]
    rk_b <- blk_triples[2, t]
    rk_c <- blk_triples[3, t]
    
    ia <- sample(block_items[[rk_a]], 1L)
    ib <- sample(block_items[[rk_b]], 1L)
    ic <- sample(block_items[[rk_c]], 1L)
    
    triples_items[, t] <- c(ia, ib, ic)
  }
  
  ## --- 5. Score triples via WST/SST implication ---
  sc <- .score_triples(P_hat, triples_items, eps = 1e-8, inclusive = TRUE)
  
  thetaW_cond <- if (sc$prem > 0L) sc$okW_prem / sc$prem else NA_real_
  thetaS_cond <- if (sc$prem > 0L) sc$okS_prem / sc$prem else NA_real_
  thetaW_all  <- if (sc$total > 0L) sc$okW_prem / sc$total else NA_real_
  thetaS_all  <- if (sc$total > 0L) sc$okS_prem / sc$total else NA_real_
  coverage    <- if (sc$total > 0L) sc$prem / sc$total else 0
  
  list(
    thetaW     = thetaW_cond,
    thetaS     = thetaS_cond,
    prem       = sc$prem,
    coverage   = coverage,
    thetaW_all = thetaW_all,
    thetaS_all = thetaS_all
  )
}



# --- Utility for summarizing chains ---
summarise_chain <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(c(mean=NA, lo=NA, hi=NA, ess=NA, mcse=NA))
  
  m <- mean(x)
  qs <- quantile(x, c(0.025, 0.975))
  
  # Fast ESS approximation if coda fails
  ess <- tryCatch(coda::effectiveSize(x), error=function(e) length(x))
  if(length(ess) > 1) ess <- mean(ess)
  
  mcse <- sd(x) / sqrt(ess)
  
  c(mean=m, lo=qs[1], hi=qs[2], ess=ess, mcse=mcse)
}


rho_from_dcsbm_draw <- function(theta_out_draw, theta_in_draw, lambda_draw,
                                z_hat, K, use_propensity = F, eps = 1e-12) {
  if (is.null(dim(lambda_draw))) {
    if (length(lambda_draw) == K * K) {
      lambda_draw <- matrix(lambda_draw, nrow = K, ncol = K)
    } else if (length(lambda_draw) == 1L) {
      lambda_draw <- matrix(lambda_draw, nrow = K, ncol = K)
    } else {
      stop("lambda_draw must be a K x K matrix or length K^2 vector.")
    }
  }
  
  if (length(dim(lambda_draw)) > 2L) {
    if (dim(lambda_draw)[3] == 1L) {
      lambda_draw <- lambda_draw[, , 1L, drop = TRUE]
    } else {
      stop("lambda_draw must be a K x K matrix or a 3D array with a single slice.")
    }
  }
  
  lambda_draw <- as.matrix(lambda_draw)
  stopifnot(all(dim(lambda_draw) == c(K, K)))
  rho <- matrix(0.5, K, K)
  if (use_propensity) {
    labs  <- factor(z_hat, levels = seq_len(K))
    out_k <- as.numeric(tapply(theta_out_draw, labs, mean, na.rm = TRUE))
    in_k  <- as.numeric(tapply(theta_in_draw,  labs, mean, na.rm = TRUE))
    g_out <- mean(theta_out_draw, na.rm = TRUE); g_in <- mean(theta_in_draw, na.rm = TRUE)
    out_k[!is.finite(out_k)] <- g_out; in_k[!is.finite(in_k)] <- g_in
  }
  for (k in 1:K) for (l in 1:K) if (k != l) {
    if (use_propensity) {
      mu_kl <- lambda_draw[k, l] * out_k[k] * in_k[l]
      mu_lk <- lambda_draw[l, k] * out_k[l] * in_k[k]
    } else {
      mu_kl <- lambda_draw[k, l]; mu_lk <- lambda_draw[l, k]
    }
    den <- mu_kl + mu_lk + eps
    rho[k, l] <- mu_kl / den
  }
  diag(rho) <- 0.5
  rho
}

# ---------------- Block-level scoring: conditional + unconditional -----------
# ---- core triple scorer (preserves SST ⇒ WST under the SST premise) ----
.score_triples <- function(R, triples, eps = 1e-8, inclusive = T) {
  R <- pmin(pmax(R, eps), 1 - eps); diag(R) <- 0.5
  gt <- if (inclusive) `>=` else `>`
  
  total <- ncol(triples)
  prem <- okW_prem <- okS_prem <- 0L
  okW_uncond <- okS_uncond <- 0L
  
  for (c in seq_len(total)) {
    a <- triples[1, c]; b <- triples[2, c]; d <- triples[3, c]
    rab <- R[a,b]; rbd <- R[b,d]; rad <- R[a,d]
    # unconditional "success" counters (diagnostic; does NOT preserve implication)
    if (gt(rad, 0.5))             okW_uncond <- okW_uncond + 1L
    if (gt(rad, max(rab, rbd)))   okS_uncond <- okS_uncond + 1L
    
    # conditional (premise) counters (the ones we use for reporting)
    if (gt(rab, 0.5) && gt(rbd, 0.5)) {
      prem <- prem + 1L
      if (gt(rad, 0.5))           okW_prem <- okW_prem + 1L
      if (gt(rad, max(rab, rbd))) okS_prem <- okS_prem + 1L
    }
  }
  list(
    total = total,
    prem = prem,
    okW_prem = okW_prem, okS_prem = okS_prem,
    okW_uncond = okW_uncond, okS_uncond = okS_uncond
  )
}

# DC-SBM relabeling is defined once in
# helper_folder/models/ordered_sbm/shared_sampler_helpers.R.
# This diagnostics file intentionally avoids carrying a second copy.

summarise_dcsbm_diagnostics <- function(fit, z_hat, K, n,
                                        m_items = 2000L, alpha = 0.5, A,
                                        T_block = 1000L, seed_block = 123, seed_items = 123,
                                        method_order = c("mean","bt","identity"),
                                        order_direction = c("strong_to_weak","weak_to_strong"),
                                        use_propensity = FALSE,
                                        bf_prior_mc = 2000L,
                                        a_lambda_prior = 1,
                                        b_lambda_prior = 1,
                                        bf_seed = 12345L) {
  method_order    <- match.arg(method_order)
  order_direction <- match.arg(order_direction)
  
  K_hat <- if (!is.null(z_hat)) {
    length(unique(z_hat[!is.na(z_hat)]))
  } else {
    NA_integer_
  }
  
  Z <- fit$z
  THO <- fit$theta_out
  THI <- fit$theta_in
  if (is.null(THO) || is.null(THI)) {
    TH <- fit$theta
    if (is.null(TH)) {
      stop("fit must provide theta_out/theta_in or a single theta matrix.")
    }
    THO <- TH
    THI <- TH
  }
  LAM <- fit$lambda
  
  if (is.matrix(Z)) {
    S <- nrow(Z)
  } else if (is.list(Z)) {
    S <- length(Z)
  } else {
    stop("fit$z must be either a matrix (S x n) or a list of length S.")
  }
  
  if (is.list(LAM)) {
    if (length(LAM) != S) {
      stop("When lambda is a list, length(lambda) must equal the number of iterations inferred from z.")
    }
  } else if (length(dim(LAM)) == 3L) {
    if (dim(LAM)[1] != S) {
      stop("First dimension of lambda array must match the number of iterations inferred from z.")
    }
    LAM <- lapply(seq_len(S), function(s) LAM[s, , ])
  } else {
    stop("fit$lambda must be either a list of matrices or a 3D array (S x K x K).")
  }
  
  get_theta_row <- function(TH, s, n) {
    if (is.matrix(TH)) {
      if (ncol(TH) != n) {
        stop("theta_out/theta_in matrix must have n columns equal to the number of nodes.")
      }
      TH[s, ]
    } else if (is.list(TH)) {
      TH[[s]]
    } else {
      stop("theta_out/theta_in must be either a matrix (S x n) or a list of length S.")
    }
  }
  
  get_z <- function(s) {
    if (is.list(Z)) Z[[s]] else Z[s, ]
  }
  
  thetaW_b_mod <- thetaS_b_mod <- prem_b_mod <- rep(NA_real_, S)
  thetaW_b_emp <- thetaS_b_emp <- prem_b_emp <- rep(NA_real_, S)
  thetaW_i_emp <- thetaS_i_emp <- prem_i_emp <- rep(NA_real_, S)
  
  thetaW_b_mod_all <- thetaS_b_mod_all <- cov_b_mod <- rep(NA_real_, S)
  thetaW_b_emp_all <- thetaS_b_emp_all <- cov_b_emp <- rep(NA_real_, S)
  thetaW_i_emp_all <- thetaS_i_emp_all <- cov_i_emp <- rep(NA_real_, S)
  
  violation_rate_est  <- rep(NA_real_, S)
  violation_count_est <- rep(NA_real_, S)
  cross_mass_est      <- rep(NA_real_, S)
  K_occ_draw          <- rep(NA_integer_, S)
  post_wst_pass       <- rep(NA_real_, S)
  post_sst_pass       <- rep(NA_real_, S)
  bar_rho_model       <- rep(NA_real_, S)

  trans_tri <- rep(NA_real_, S)
  cycle_mass_w <- rep(NA_real_, S)
  min_back_w <- rep(NA_real_, S)
  min_back_w_norm <- rep(NA_real_, S)
  hier_energy <- rep(NA_real_, S)
  curl_energy <- rep(NA_real_, S)
  
  for (s in seq_len(S)) {
    z_s <- get_z(s)
    lam_s <- as.matrix(LAM[[s]])
    occ_s <- sort(unique(z_s[!is.na(z_s)]))
    K_occ <- length(occ_s)
    K_occ_draw[s] <- K_occ

    if (K_occ < 1L) {
      next
    }

    z_s_compact <- match(z_s, occ_s)

    if (all(occ_s <= nrow(lam_s)) && all(occ_s <= ncol(lam_s))) {
      lam_occ <- lam_s[occ_s, occ_s, drop = FALSE]
    } else if (nrow(lam_s) == K_occ && ncol(lam_s) == K_occ) {
      lam_occ <- lam_s
    } else {
      stop(
        "Cannot align lambda dimensions with occupied labels at draw ", s,
        ": dim(lambda)=(", nrow(lam_s), ",", ncol(lam_s), "), K_occ=", K_occ,
        ", max_label=", max(occ_s)
      )
    }
    
    tho <- get_theta_row(THO, s, n)
    thi <- get_theta_row(THI, s, n)
    rho_occ <- rho_from_dcsbm_draw(tho, thi, lam_occ, z_hat = z_s_compact, K = K_occ,
                                    use_propensity = use_propensity)
    
    if (K_occ >= 3L && (is.na(K_hat) || K_occ == K_hat)) {
      blk_mod <- block_diag_rates_ext(
        rho_occ, T_block = T_block, seed = seed_block,
        method_order = method_order, order_direction = order_direction
      )
      bar_rho_model[s] <- .mean_forward_rho(
        rho_occ, method_order = method_order,
        order_direction = order_direction
      )
      chk <- .check_wst_sst_sorted(rho_occ)
      post_wst_pass[s] <- as.numeric(isTRUE(chk$wst))
      post_sst_pass[s] <- as.numeric(isTRUE(chk$sst))
      inv_diag <- .diag_from_rho_invariant(rho_occ)
      trans_tri[s]       <- inv_diag$transitive_triads
      cycle_mass_w[s]    <- inv_diag$cycle_mass_weighted
      min_back_w[s]      <- inv_diag$min_backward_weight
      min_back_w_norm[s] <- inv_diag$min_backward_weight_norm
      hier_energy[s]     <- inv_diag$hierarchy_energy
      curl_energy[s]     <- inv_diag$curl_energy
      thetaW_b_mod[s]     <- blk_mod$thetaW
      thetaS_b_mod[s]     <- blk_mod$thetaS
      prem_b_mod[s]       <- blk_mod$prem
      thetaW_b_mod_all[s] <- blk_mod$thetaW_all
      thetaS_b_mod_all[s] <- blk_mod$thetaS_all
      cov_b_mod[s]        <- blk_mod$coverage
    } else {
      thetaW_b_mod[s]     <- NA_real_
      thetaS_b_mod[s]     <- NA_real_
      prem_b_mod[s]       <- NA_real_
      thetaW_b_mod_all[s] <- NA_real_
      thetaS_b_mod_all[s] <- NA_real_
      cov_b_mod[s]        <- NA_real_
      post_wst_pass[s]    <- NA_real_
      post_sst_pass[s]    <- NA_real_
      bar_rho_model[s]    <- NA_real_
    }
    
    if (K_occ >= 3L && (is.na(K_hat) || K_occ == K_hat)) {
      blk_emp <- .empirical_block_rates_one(
        A = A, z_vec = z_s, alpha = alpha,
        T_block = T_block, seed = seed_block,
        method_order = method_order, order_direction = order_direction
      )
      itm_emp <- .empirical_item_rates_one(
        A = A, z_vec = z_s, T = m_items,
        alpha = alpha, seed = seed_items,
        method_order = method_order, order_direction = order_direction
      )
      
      thetaW_b_emp[s]     <- blk_emp$thetaW
      thetaS_b_emp[s]     <- blk_emp$thetaS
      prem_b_emp[s]       <- blk_emp$prem
      thetaW_b_emp_all[s] <- blk_emp$thetaW_all
      thetaS_b_emp_all[s] <- blk_emp$thetaS_all
      cov_b_emp[s]        <- blk_emp$coverage
      
      thetaW_i_emp[s]     <- itm_emp$thetaW
      thetaS_i_emp[s]     <- itm_emp$thetaS
      prem_i_emp[s]       <- itm_emp$prem
      thetaW_i_emp_all[s] <- itm_emp$thetaW_all
      thetaS_i_emp_all[s] <- itm_emp$thetaS_all
      cov_i_emp[s]        <- itm_emp$coverage
    } else {
      thetaW_b_emp[s]     <- NA_real_
      thetaS_b_emp[s]     <- NA_real_
      prem_b_emp[s]       <- NA_real_
      thetaW_b_emp_all[s] <- NA_real_
      thetaS_b_emp_all[s] <- NA_real_
      cov_b_emp[s]        <- NA_real_
      thetaW_i_emp[s]     <- NA_real_
      thetaS_i_emp[s]     <- NA_real_
      prem_i_emp[s]       <- NA_real_
      thetaW_i_emp_all[s] <- NA_real_
      thetaS_i_emp_all[s] <- NA_real_
      cov_i_emp[s]        <- NA_real_
    }
    
    z_ranked <- .rank_partition_by_strength(
      A = A, z_vec = z_s, alpha = alpha, order_direction = order_direction
    )
    rows <- row(A); cols <- col(A); vals <- A
    z_i  <- z_ranked[rows]; z_j <- z_ranked[cols]
    
    cross_mask <- (z_i != z_j) & (vals > 0)
    if (any(cross_mask, na.rm = TRUE)) {
      viol_mask <- cross_mask & (z_i > z_j)
      v_mass <- sum(vals[viol_mask], na.rm = TRUE)
      t_mass <- sum(vals[cross_mask], na.rm = TRUE)
      violation_count_est[s] <- v_mass
      cross_mass_est[s]      <- t_mass
      violation_rate_est[s]  <- if (t_mass > 0) v_mass / t_mass else NA_real_
    } else {
      violation_rate_est[s]  <- NA_real_
      violation_count_est[s] <- NA_real_
      cross_mass_est[s]      <- NA_real_
    }
  }

  p_post_wst <- mean(post_wst_pass, na.rm = TRUE)
  p_post_sst <- mean(post_sst_pass, na.rm = TRUE)

  K_bf <- suppressWarnings(as.integer(K_hat))
  if (!is.finite(K_bf) || K_bf < 3L) {
    K_cand <- K_occ_draw[is.finite(K_occ_draw) & K_occ_draw >= 3L]
    if (length(K_cand)) K_bf <- as.integer(round(stats::median(K_cand)))
  }

  prior_probs <- if (is.finite(K_bf) && K_bf >= 3L) {
    .prior_constraint_prob_mc(
      K = K_bf,
      n_mc = as.integer(bf_prior_mc),
      a_lambda = a_lambda_prior,
      b_lambda = b_lambda_prior,
      seed = bf_seed
    )
  } else {
    list(p_wst = NA_real_, p_sst = NA_real_, n_mc = as.integer(bf_prior_mc))
  }

  p_prior_wst <- prior_probs$p_wst
  p_prior_sst <- prior_probs$p_sst

  bf_ratio <- function(num, den) {
    if (is.na(num) || is.na(den)) return(NA_real_)
    if (den == 0 && num > 0) return(Inf)
    if (den == 0 && num == 0) return(NA_real_)
    num / den
  }

  bf_wst_0 <- bf_ratio(p_post_wst, p_prior_wst)
  bf_sst_0 <- bf_ratio(p_post_sst, p_prior_sst)
  bf_sst_wst <- bf_ratio(bf_sst_0, bf_wst_0)

  safe_log <- function(x) ifelse(is.finite(x) & x > 0, log(x), NA_real_)
  
  res <- c(
    setNames(summarise_chain(thetaW_b_mod),     paste0("thetaW_block_model_", c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_b_mod),     paste0("thetaS_block_model_", c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaW_b_emp),     paste0("thetaW_block_emp_",   c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_b_emp),     paste0("thetaS_block_emp_",   c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaW_i_emp),     paste0("thetaW_item_emp_",    c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_i_emp),     paste0("thetaS_item_emp_",    c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(bar_rho_model),    paste0("bar_rho_model_",      c("mean","lo","hi","ess","mcse"))),
    
    setNames(summarise_chain(thetaW_b_mod_all), paste0("thetaW_block_model_all_", c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_b_mod_all), paste0("thetaS_block_model_all_", c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaW_b_emp_all), paste0("thetaW_block_emp_all_",   c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_b_emp_all), paste0("thetaS_block_emp_all_",   c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaW_i_emp_all), paste0("thetaW_item_emp_all_",    c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_i_emp_all), paste0("thetaS_item_emp_all_",    c("mean","lo","hi","ess","mcse"))),
    
    prem_block_model_avg      = mean(prem_b_mod, na.rm = TRUE),
    prem_block_emp_avg        = mean(prem_b_emp,  na.rm = TRUE),
    prem_item_emp_avg         = mean(prem_i_emp,  na.rm = TRUE),
    coverage_block_model_avg  = mean(cov_b_mod,   na.rm = TRUE),
    coverage_block_emp_avg    = mean(cov_b_emp,   na.rm = TRUE),
    coverage_item_emp_avg     = mean(cov_i_emp,   na.rm = TRUE),

    transitive_triads_mean         = mean(trans_tri, na.rm = TRUE),
    transitive_triads_lo           = quantile(trans_tri, 0.025, na.rm = TRUE),
    transitive_triads_hi           = quantile(trans_tri, 0.975, na.rm = TRUE),
    cycle_mass_weighted_mean       = mean(cycle_mass_w, na.rm = TRUE),
    min_backward_weight_mean       = mean(min_back_w, na.rm = TRUE),
    min_backward_weight_norm_mean  = mean(min_back_w_norm, na.rm = TRUE),
    hierarchy_energy_mean          = mean(hier_energy, na.rm = TRUE),
    hierarchy_energy_lo            = quantile(hier_energy, 0.025, na.rm = TRUE),
    hierarchy_energy_hi            = quantile(hier_energy, 0.975, na.rm = TRUE),
    curl_energy_mean               = mean(curl_energy, na.rm = TRUE),

    p_post_wst                     = p_post_wst,
    p_prior_wst                    = p_prior_wst,
    bf_wst_0                       = bf_wst_0,
    log_bf_wst_0                   = safe_log(bf_wst_0),
    p_post_sst                     = p_post_sst,
    p_prior_sst                    = p_prior_sst,
    bf_sst_0                       = bf_sst_0,
    log_bf_sst_0                   = safe_log(bf_sst_0),
    bf_sst_wst                     = bf_sst_wst,
    log_bf_sst_wst                 = safe_log(bf_sst_wst),
    bf_prior_mc                    = as.numeric(prior_probs$n_mc),
    
    violation_rate_mean       = mean(violation_rate_est,  na.rm = TRUE),
    violation_rate_lo         = quantile(violation_rate_est, 0.025, na.rm = TRUE),
    violation_rate_hi         = quantile(violation_rate_est, 0.975, na.rm = TRUE),
    violation_count_mean      = mean(violation_count_est, na.rm = TRUE),
    cross_mass_mean           = mean(cross_mass_est,      na.rm = TRUE),
    n_draws_K_match           = sum(!is.na(K_occ_draw) & (is.na(K_hat) | K_occ_draw == K_hat))
  )
  
  if (!is.null(z_hat)) {
    rows <- row(A); cols <- col(A); vals <- A
    z_i  <- z_hat[rows]; z_j <- z_hat[cols]
    cross_mask <- (z_i != z_j) & (vals > 0)
    if (any(cross_mask)) {
      viol_mask <- cross_mask & (z_i > z_j)
      v_mass <- sum(vals[viol_mask])
      t_mass <- sum(vals[cross_mask])
      res["violation_rate_zhat"]   <- v_mass / t_mass
      res["violation_count_zhat"]  <- v_mass
    } else {
      res["violation_rate_zhat"]   <- 0
      res["violation_count_zhat"]  <- 0
    }
  }
  
  res
}
summarise_osbm_diagnostics <- function(out_relab, regime, K_max_hint, z_hat, n,
                                       m_items = 2000L, alpha = 0.5, A,
                                       T_block = 1000L, seed_block = 123, seed_items = 123,
                                       method_order = c("mean","bt","identity"),
                                       order_direction = c("strong_to_weak","weak_to_strong")) {
  
  method_order    <- match.arg(method_order)
  order_direction <- match.arg(order_direction)
  
  ## --- K_hat from point estimate -------------------------
  K_hat <- if (!is.null(z_hat)) {
    length(unique(z_hat[!is.na(z_hat)]))
  } else {
    NA_integer_
  }
  
  ## --- 1. Handle Ragged Input (Variable K) ---------------
  is_ragged <- is.list(out_relab$psi) && !is.array(out_relab$psi)
  S <- if (is_ragged) length(out_relab$psi) else dim(out_relab$psi)[1]
  
  get_rho <- function(s) {
    if (is_ragged) {
      psi_s <- out_relab$psi[[s]]
      if (regime == "WST") {
        K_s <- nrow(psi_s)
        rho_from_osbm_draw(psi_s, regime, K_s)
      } else {
        K_s <- length(psi_s) + 1L
        rho_from_osbm_draw(psi_s, regime, K_s)
      }
    } else {
      if (regime == "WST") {
        rho_from_osbm_draw(out_relab$psi[s, , ], regime, K_max_hint)
      } else {
        rho_from_osbm_draw(out_relab$psi[s, ], regime, K_max_hint)
      }
    }
  }
  
  get_z <- function(s) {
    if (is_ragged) out_relab$z[[s]] else out_relab$z[s, ]
  }
  
  ## --- 2. Holders ---------------------------------------
  thetaW_b_mod <- thetaS_b_mod <- prem_b_mod <- rep(NA_real_, S)
  thetaW_b_emp <- thetaS_b_emp <- prem_b_emp <- rep(NA_real_, S)
  thetaW_i_emp <- thetaS_i_emp <- prem_i_emp <- rep(NA_real_, S)
  
  thetaW_b_mod_all <- thetaS_b_mod_all <- cov_b_mod <- rep(NA_real_, S)
  thetaW_b_emp_all <- thetaS_b_emp_all <- cov_b_emp <- rep(NA_real_, S)
  thetaW_i_emp_all <- thetaS_i_emp_all <- cov_i_emp <- rep(NA_real_, S)
  
  violation_rate_est  <- rep(NA_real_, S)
  violation_count_est <- rep(NA_real_, S)
  cross_mass_est      <- rep(NA_real_, S)
  bar_rho_model       <- rep(NA_real_, S)
  K_occ_draw          <- rep(NA_integer_, S)

  trans_tri <- rep(NA_real_, S)
  cycle_mass_w <- rep(NA_real_, S)
  min_back_w <- rep(NA_real_, S)
  min_back_w_norm <- rep(NA_real_, S)
  hier_energy <- rep(NA_real_, S)
  curl_energy <- rep(NA_real_, S)
  
  ## --- 3. Loop over draws -------------------------------
  for (s in seq_len(S)) {
    rho_full <- get_rho(s)
    z_s      <- get_z(s)
    
    occ_s <- sort(unique(z_s[!is.na(z_s)]))
    K_s   <- length(occ_s)
    K_occ_draw[s] <- K_s
    
    ## --- A. Model-based transitivity diagnostics (only if K_s >= 3 and (optionally) matches K_hat) ---
    if (K_s >= 3L && (is.na(K_hat) || K_s == K_hat)) {
      rho_occ <- rho_full[occ_s, occ_s, drop = FALSE]
      blk_mod <- block_diag_rates_ext(
        rho_occ, T_block = T_block, seed = seed_block,
        method_order = method_order, order_direction = order_direction
      )
      bar_rho_model[s] <- .mean_forward_rho(
        rho_occ, method_order = method_order,
        order_direction = order_direction
      )
      inv_diag <- .diag_from_rho_invariant(rho_occ)
      trans_tri[s]       <- inv_diag$transitive_triads
      cycle_mass_w[s]    <- inv_diag$cycle_mass_weighted
      min_back_w[s]      <- inv_diag$min_backward_weight
      min_back_w_norm[s] <- inv_diag$min_backward_weight_norm
      hier_energy[s]     <- inv_diag$hierarchy_energy
      curl_energy[s]     <- inv_diag$curl_energy
      thetaW_b_mod[s]     <- blk_mod$thetaW
      thetaS_b_mod[s]     <- blk_mod$thetaS
      prem_b_mod[s]       <- blk_mod$prem
      thetaW_b_mod_all[s] <- blk_mod$thetaW_all
      thetaS_b_mod_all[s] <- blk_mod$thetaS_all
      cov_b_mod[s]        <- blk_mod$coverage
    } else {
      thetaW_b_mod[s]     <- NA_real_
      thetaS_b_mod[s]     <- NA_real_
      prem_b_mod[s]       <- NA_real_
      thetaW_b_mod_all[s] <- NA_real_
      thetaS_b_mod_all[s] <- NA_real_
      cov_b_mod[s]        <- NA_real_
      bar_rho_model[s]    <- NA_real_
    }
    
    ## --- B. Empirical (block + item) diagnostics, restricted to same-K draws ---
    if (K_s >= 3L && (is.na(K_hat) || K_s == K_hat)) {
      blk_emp <- .empirical_block_rates_one(
        A = A, z_vec = z_s, alpha = alpha,
        T_block = T_block, seed = seed_block,
        method_order = method_order, order_direction = order_direction
      )
      itm_emp <- .empirical_item_rates_one(
        A = A, z_vec = z_s, T = m_items,
        alpha = alpha, seed = seed_items,
        method_order = method_order, order_direction = order_direction
      )
      
      thetaW_b_emp[s]     <- blk_emp$thetaW
      thetaS_b_emp[s]     <- blk_emp$thetaS
      prem_b_emp[s]       <- blk_emp$prem
      thetaW_b_emp_all[s] <- blk_emp$thetaW_all
      thetaS_b_emp_all[s] <- blk_emp$thetaS_all
      cov_b_emp[s]        <- blk_emp$coverage
      
      thetaW_i_emp[s]     <- itm_emp$thetaW
      thetaS_i_emp[s]     <- itm_emp$thetaS
      prem_i_emp[s]       <- itm_emp$prem
      thetaW_i_emp_all[s] <- itm_emp$thetaW_all
      thetaS_i_emp_all[s] <- itm_emp$thetaS_all
      cov_i_emp[s]        <- itm_emp$coverage
    }
    
    ## --- C. Hierarchy-violating mass under the intrinsic OSBM order ---
    vstats <- .violation_stats_from_ranked(A = A, z_ranked = z_s)
    violation_rate_est[s]  <- vstats$rate
    violation_count_est[s] <- vstats$count
    cross_mass_est[s]      <- vstats$cross_mass
  }
  
  
  ## --- 4. Summarise over s (dropping NAs) ----------------
  res <- c(
    setNames(summarise_chain(thetaW_b_mod),     paste0("thetaW_block_model_", c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_b_mod),     paste0("thetaS_block_model_", c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaW_b_emp),     paste0("thetaW_block_emp_",   c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_b_emp),     paste0("thetaS_block_emp_",   c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaW_i_emp),     paste0("thetaW_item_emp_",    c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_i_emp),     paste0("thetaS_item_emp_",    c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(bar_rho_model),    paste0("bar_rho_model_",      c("mean","lo","hi","ess","mcse"))),
    
    setNames(summarise_chain(thetaW_b_mod_all), paste0("thetaW_block_model_all_", c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_b_mod_all), paste0("thetaS_block_model_all_", c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaW_b_emp_all), paste0("thetaW_block_emp_all_",   c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_b_emp_all), paste0("thetaS_block_emp_all_",   c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaW_i_emp_all), paste0("thetaW_item_emp_all_",    c("mean","lo","hi","ess","mcse"))),
    setNames(summarise_chain(thetaS_i_emp_all), paste0("thetaS_item_emp_all_",    c("mean","lo","hi","ess","mcse"))),
    
    prem_block_model_avg      = mean(prem_b_mod, na.rm = TRUE),
    prem_block_emp_avg        = mean(prem_b_emp,  na.rm = TRUE),
    prem_item_emp_avg         = mean(prem_i_emp,  na.rm = TRUE),
    coverage_block_model_avg  = mean(cov_b_mod,   na.rm = TRUE),
    coverage_block_emp_avg    = mean(cov_b_emp,   na.rm = TRUE),
    coverage_item_emp_avg     = mean(cov_i_emp,   na.rm = TRUE),

    transitive_triads_mean         = mean(trans_tri, na.rm = TRUE),
    transitive_triads_lo           = quantile(trans_tri, 0.025, na.rm = TRUE),
    transitive_triads_hi           = quantile(trans_tri, 0.975, na.rm = TRUE),
    cycle_mass_weighted_mean       = mean(cycle_mass_w, na.rm = TRUE),
    min_backward_weight_mean       = mean(min_back_w, na.rm = TRUE),
    min_backward_weight_norm_mean  = mean(min_back_w_norm, na.rm = TRUE),
    hierarchy_energy_mean          = mean(hier_energy, na.rm = TRUE),
    hierarchy_energy_lo            = quantile(hier_energy, 0.025, na.rm = TRUE),
    hierarchy_energy_hi            = quantile(hier_energy, 0.975, na.rm = TRUE),
    curl_energy_mean               = mean(curl_energy, na.rm = TRUE),

    p_post_wst                     = NA_real_,
    p_prior_wst                    = NA_real_,
    bf_wst_0                       = NA_real_,
    log_bf_wst_0                   = NA_real_,
    p_post_sst                     = NA_real_,
    p_prior_sst                    = NA_real_,
    bf_sst_0                       = NA_real_,
    log_bf_sst_0                   = NA_real_,
    bf_sst_wst                     = NA_real_,
    log_bf_sst_wst                 = NA_real_,
    bf_prior_mc                    = NA_real_,
    
    violation_rate_mean       = mean(violation_rate_est,  na.rm = TRUE),
    violation_rate_lo         = quantile(violation_rate_est, 0.025, na.rm = TRUE),
    violation_rate_hi         = quantile(violation_rate_est, 0.975, na.rm = TRUE),
    violation_count_mean      = mean(violation_count_est, na.rm = TRUE),
    cross_mass_mean           = mean(cross_mass_est,      na.rm = TRUE),
    n_draws_K_match           = sum(!is.na(K_occ_draw) & (is.na(K_hat) | K_occ_draw == K_hat))
  )
  
  ## --- 5. Point-estimate violation rate for z_hat --------
  if (!is.null(z_hat)) {
    vstats <- .violation_stats_from_ranked(A = A, z_ranked = z_hat)
    res["violation_rate_zhat"]   <- vstats$rate
    res["violation_count_zhat"]  <- vstats$count
  }
  
  res
}
block_diag_rates_ext <- function(rho,
                                 T_block = 1000L,
                                 seed = NULL,
                                 method_order = c("mean","bt","identity"),
                                 order_direction = c("strong_to_weak","weak_to_strong"),
                                 eps = 1e-8) {
  method_order    <- match.arg(method_order)
  order_direction <- match.arg(order_direction)
  if (!is.null(seed)) set.seed(seed)
  
  K <- nrow(rho)
  if (K < 3L) {
    return(list(
      thetaW    = NA_real_,
      thetaS    = NA_real_,
      prem      = 0L,
      coverage  = 0,
      thetaW_all = NA_real_,
      thetaS_all = NA_real_
    ))
  }
  
  ## 1. Order blocks. WST/SST diagnostics pass method_order="identity";
  ## DC-SBM diagnostics pass a strength-based ordering.
  R <- .ordered_rho(
    rho, method_order = method_order,
    order_direction = order_direction, eps = eps
  )
  
  ## 2. Triple set (possibly subsampled)
  all_triples <- utils::combn(K, 3L)
  total <- ncol(all_triples)
  if (!is.null(T_block) && T_block < total) {
    keep <- sort(sample.int(total, T_block, replace = FALSE))
    triples <- all_triples[, keep, drop = FALSE]
  } else {
    triples <- all_triples
  }
  
  ## 3. Score triples: WST/SST under the premise
  sc <- .score_triples(R, triples, eps = eps, inclusive = TRUE)
  
  thetaW_cond <- if (sc$prem > 0L) sc$okW_prem / sc$prem else NA_real_
  thetaS_cond <- if (sc$prem > 0L) sc$okS_prem / sc$prem else NA_real_
  thetaW_all  <- if (sc$total > 0L) sc$okW_prem / sc$total else NA_real_
  thetaS_all  <- if (sc$total > 0L) sc$okS_prem / sc$total else NA_real_
  coverage    <- if (sc$total > 0L) sc$prem / sc$total else 0
  
  list(
    thetaW    = thetaW_cond,
    thetaS    = thetaS_cond,
    prem      = sc$prem,
    coverage  = coverage,
    thetaW_all = thetaW_all,
    thetaS_all = thetaS_all
  )
}

.empty_diag_vec <- function() {
  base_stats <- c("mean", "lo", "hi", "ess", "mcse")
  diag_prefixes <- c(
    "thetaW_block_model", "thetaS_block_model",
    "thetaW_block_emp", "thetaS_block_emp",
    "thetaW_item_emp", "thetaS_item_emp",
    "bar_rho_model",
    "thetaW_block_model_all", "thetaS_block_model_all",
    "thetaW_block_emp_all", "thetaS_block_emp_all",
    "thetaW_item_emp_all", "thetaS_item_emp_all"
  )
  diag_names <- unlist(lapply(diag_prefixes, function(prefix) {
    paste0(prefix, "_", base_stats)
  }))
  tail_names <- c(
    "prem_block_model_avg", "prem_block_emp_avg", "prem_item_emp_avg",
    "coverage_block_model_avg", "coverage_block_emp_avg", "coverage_item_emp_avg",
    "transitive_triads_mean", "transitive_triads_lo", "transitive_triads_hi",
    "cycle_mass_weighted_mean",
    "min_backward_weight_mean", "min_backward_weight_norm_mean",
    "hierarchy_energy_mean", "hierarchy_energy_lo", "hierarchy_energy_hi",
    "curl_energy_mean",
    "p_post_wst", "p_prior_wst", "bf_wst_0", "log_bf_wst_0",
    "p_post_sst", "p_prior_sst", "bf_sst_0", "log_bf_sst_0",
    "bf_sst_wst", "log_bf_sst_wst", "bf_prior_mc",
    "violation_rate_mean", "violation_rate_lo", "violation_rate_hi",
    "violation_count_mean", "cross_mass_mean", "n_draws_K_match",
    "violation_rate_zhat", "violation_count_zhat"
  )
  all_names <- c(diag_names, tail_names)
  out <- rep(NA_real_, length(all_names))
  names(out) <- all_names
  out
}

# DC-SBM dyad log-likelihood helpers are defined once in
# helper_folder/models/ordered_sbm/shared_sampler_helpers.R.




loglik_matrix_modular <- function(A, mcmc_out, regime = c("WST", "SST"), dyad_index = NULL) {
  regime <- match.arg(regime)
  
  # 1. Setup Data Structures
  A_mat <- as.matrix(A)
  n <- nrow(A_mat)
  
  # Create dyad indices (upper triangle only for counting pairs)
  # Note: We compute LL for the *pair* (A_ij, A_ji) to be consistent with dyadic independence
  if (is.null(dyad_index)) {
    dyad_index <- which(upper.tri(A_mat), arr.ind = TRUE)
  }
  n_dyads <- nrow(dyad_index)
  
  # Extract chains
  z_chain     <- mcmc_out$z
  eta_chain   <- mcmc_out$eta
  kappa_chain <- mcmc_out$kappa
  psi_chain   <- mcmc_out$psi
  
  n_samples <- length(z_chain)
  
  # Pre-allocate result matrix (Samples x Dyads)
  # Initialize with -Inf to catch errors
  ll_matrix <- matrix(-Inf, nrow = n_samples, ncol = n_dyads)
  
  # Pre-fetch data vectors
  row_idx <- dyad_index[, 1]
  col_idx <- dyad_index[, 2]
  
  y_ij <- A_mat[cbind(row_idx, col_idx)]
  y_ji <- A_mat[cbind(col_idx, row_idx)]
  N_dyad <- y_ij + y_ji
  
  # 2. Loop over MCMC Samples
  for (t in seq_len(n_samples)) {
    
    # Get params for this iteration
    z_t     <- z_chain[[t]] # Vector length n
    eta_t   <- eta_chain[[t]] # Vector length n
    kappa_t <- kappa_chain[[t]] # Matrix K_t x K_t
    psi_t   <- psi_chain[[t]]   # WST: Matrix, SST: Vector
    
    # Current K
    K_t <- nrow(kappa_t)
    
    # --- A. Compute Volume Intensities (Lambda) ---
    # lambda_ij = eta_i * eta_j * kappa_{z_i, z_j}
    # Map z_t to indices
    k_i <- z_t[row_idx]
    k_j <- z_t[col_idx]
    
    # Vectorized Kappa lookup
    # Note: kappa is symmetric
    kappa_vals <- kappa_t[cbind(k_i, k_j)]
    
    lambda_vals <- eta_t[row_idx] * eta_t[col_idx] * kappa_vals
    
    # --- B. Compute Directional Probabilities (Rho) ---
    psi_vals <- numeric(n_dyads)
    
    # Direction indicators: s = 1 if z_i < z_j (i is stronger), -1 if z_i > z_j
    # Note: In your model, lower index = stronger.
    # s_{k,l} = sgn(l - k). 
    # If k=1 (top), l=2 (bottom), l-k = 1. s=1. 
    # Prob(i->j) = logit^-1(s * psi)
    
    s_vals <- sign(k_j - k_i) 
    
    if (regime == "WST") {
      # WST: psi is a KxK matrix. We need psi_{min, max}
      k_min <- pmin(k_i, k_j)
      k_max <- pmax(k_i, k_j)
      
      # Handle singleton K (if K=1, psi might be empty/null or 1x1 zero)
      if (K_t > 1) {
        # We only need lookups where blocks differ
        diff_mask <- k_i != k_j
        if (any(diff_mask)) {
          # Extract upper tri values
          psi_vals[diff_mask] <- psi_t[cbind(k_min[diff_mask], k_max[diff_mask])]
        }
      }
      
    } else {
      # SST: psi is a vector of distances
      # distance d = |k_i - k_j|
      d_vals <- abs(k_i - k_j)
      
      if (K_t > 1) {
        diff_mask <- d_vals > 0
        if (any(diff_mask)) {
          # psi_t has length K-1. psi_t[d] gives log-odds for distance d
          # Safety check for index bounds
          valid_d <- d_vals[diff_mask]
          psi_vals[diff_mask] <- psi_t[valid_d] 
        }
      }
    }
    
    # Compute Rho (Prob i -> j)
    # logit(rho) = s * psi
    # If same block (s=0), logit=0 -> rho=0.5. Correct.
    logits <- s_vals * psi_vals
    
    # Stable Log-Sum-Exp for log probabilities to avoid underflow
    # log_rho = -log(1 + exp(-logits))
    # log_1m_rho = -log(1 + exp(logits))
    log_rho     <- -log1pexp(-logits)
    log_1m_rho  <- -log1pexp(logits)
    
    # --- C. Compute Total Log-Likelihood per Dyad ---
    # LL_Poisson = dpois(N, lambda, log=T)
    # LL_Binom   = dbinom(A_ij, N, rho, log=T)
    
    ll_pois <- dpois(N_dyad, lambda_vals, log = TRUE)
    
    # Manual Binomial LogLik to use stable log_rho
    # log(choose(N, k)) + k*log(p) + (N-k)*log(1-p)
    # We ignore log(choose) because it cancels in LOO comparison (constant data), 
    # BUT for absolute LOO-IC it's needed. 'dbinom' handles it.
    # We calculate A_ij successes with prob rho
    ll_binom <- lchoose(N_dyad, y_ij) + y_ij * log_rho + y_ji * log_1m_rho
    
    # Combine
    ll_matrix[t, ] <- ll_pois + ll_binom
  }
  
  return(ll_matrix)
}


  
compute_block_scores <- function(z, A, N) {
  n <- length(z)
  K <- length(unique(z))
  r_bar <- rep(NA_real_, K)
  
  for (k in seq_len(K)) {
    idx_k <- which(z == k)
    nk <- length(idx_k)
    if (nk == 0L) {
      next
    }
    
    per_node_scores <- rep(NA_real_, nk)
    
    for (h in seq_along(idx_k)) {
      i <- idx_k[h]
      j_all <- setdiff(seq_len(n), i)
      
      Nij <- N[i, j_all]
      Aij <- A[i, j_all]
      
      valid <- !is.na(Nij) & Nij > 0
      if (any(valid)) {
        per_node_scores[h] <- sum(Aij[valid] / Nij[valid])
      } else {
        per_node_scores[h] <- NA_real_
      }
    }
    
    r_bar[k] <- mean(per_node_scores, na.rm = TRUE)
  }
  
  r_bar
}

relabel_by_block_score_single <- function(z, lambda, A, N) {
  # z: length n, labels in {1,...,K_t}
  # lambda: K_t x K_t matrix
  # A, N: n x n matrices
  
  r_bar <- compute_block_scores(z, A, N)
  K <- length(r_bar)
  
  # permutation: highest score -> label 1, etc.
  perm <- order(r_bar, decreasing = TRUE, na.last = TRUE)
  
  # inverse permutation: new_label = inv_perm[old_label]
  inv_perm <- integer(K)
  inv_perm[perm] <- seq_len(K)
  
  z_new <- inv_perm[z]
  lambda_new <- lambda[perm, perm, drop = FALSE]
  
  list(
    z      = z_new,
    lambda = lambda_new,
    r_bar  = r_bar,
    perm   = perm
  )
}

relabel_chain_by_block_score <- function(Z_chain, Lambda_chain, A, N) {
  # Z_chain:
  #   - matrix: n_iter x n
  #   - OR list: length n_iter, each an integer vector of length n
  # Lambda_chain:
  #   - list of length n_iter, each element a K_t x K_t matrix
  # A, N: n x n matrices
  
  ## --- handle Z_chain type ---
  is_mat_z <- is.matrix(Z_chain)
  if (is_mat_z) {
    n_iter <- nrow(Z_chain)
    n      <- ncol(Z_chain)
  } else if (is.list(Z_chain)) {
    n_iter <- length(Z_chain)
    n      <- length(Z_chain[[1]])
  } else {
    stop("Z_chain must be either a matrix (n_iter x n) or a list of length n_iter.")
  }
  
  ## --- Lambda as list of matrices with varying K_t ---
  if (!is.list(Lambda_chain)) {
    stop("Lambda_chain must be a list of K_t x K_t matrices (one per iteration).")
  }
  if (length(Lambda_chain) != n_iter) {
    stop("Length of Lambda_chain list must match the number of iterations in Z_chain.")
  }
  
  ## --- outputs ---
  if (is_mat_z) {
    Z_out <- matrix(NA_integer_, nrow = n_iter, ncol = n)
  } else {
    Z_out <- vector("list", n_iter)
  }
  
  Lambda_out <- vector("list", n_iter)
  scores     <- vector("list", n_iter)  # each element: numeric vector length K_t
  perms      <- vector("list", n_iter)  # each element: integer vector length K_t
  
  ## --- loop over iterations ---
  for (t in seq_len(n_iter)) {
    z_t <- if (is_mat_z) Z_chain[t, ] else Z_chain[[t]]
    lambda_t <- Lambda_chain[[t]]
    
    rel <- relabel_by_block_score_single(
      z      = z_t,
      lambda = lambda_t,
      A      = A,
      N      = N
    )
    
    if (is_mat_z) {
      Z_out[t, ] <- rel$z
    } else {
      Z_out[[t]] <- rel$z
    }
    
    Lambda_out[[t]] <- rel$lambda
    scores[[t]]     <- rel$r_bar
    perms[[t]]      <- rel$perm
  }
  
  list(
    Z      = Z_out,
    Lambda = Lambda_out,
    scores = scores,
    perms  = perms
  )
}
