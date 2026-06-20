# Simulation-study helpers for generating, relabelling, and summarising ordered SBM runs.
#
# relabel_osbm <- function(mcmc_out,
#                          A,
#                          ordering = c("WST","SST","NONE"),
#                          score = c("success","outdeg","indeg","netdeg")) {
#   ordering <- match.arg(ordering)
#   score    <- match.arg(score)
#   
#   # SST or NONE: do nothing (we don't need perm in this script)
#   if (ordering != "WST") {
#     attr(mcmc_out, "relabel_ordering") <- ordering
#     attr(mcmc_out, "relabel_score")    <- NA_character_
#     return(mcmc_out)
#   }
#   
#   # --- WST relabelling by block score (A-driven) -----------------------
#   if (inherits(A, "Matrix")) A <- as.matrix(A)
#   stopifnot(is.matrix(A), nrow(A) == ncol(A))
#   n <- nrow(A)
#   
#   # per-node metrics (fixed, from data)
#   outdeg   <- rowSums(A)
#   indeg    <- colSums(A)
#   netdeg   <- outdeg - indeg
#   matches  <- A + t(A)
#   success  <- outdeg / pmax(rowSums(matches), 1)  # row-wise success rate
#   
#   node_metric <- switch(score,
#                         outdeg  = outdeg,
#                         indeg   = indeg,
#                         netdeg  = netdeg,
#                         success = success)
#   
#   # z is assumed matrix: S x n
#   draws_z <- mcmc_out$z
#   if (!is.matrix(draws_z)) {
#     draws_z <- as.matrix(draws_z)
#   }
#   S <- nrow(draws_z)
#   
#   # kappa / psi may be arrays OR lists OR absent; we handle robustly
#   draws_kap <- mcmc_out$kappa
#   draws_psi <- mcmc_out$psi
#   
#   is_kap_list <- is.list(draws_kap)
#   is_psi_list <- is.list(draws_psi)
#   
#   relab_z   <- matrix(NA_integer_, nrow = S, ncol = n)
#   relab_kap <- if (is_kap_list) vector("list", S) else draws_kap
#   relab_psi <- if (is_psi_list) vector("list", S) else draws_psi
#   perms     <- vector("list", S)      # each draw can have its own K_s
#   
#   for (s in seq_len(S)) {
#     z_s <- as.integer(draws_z[s, ])
#     K_s <- max(z_s, na.rm = TRUE)
#     
#     # block means; empty blocks -> -Inf so they go to the tail
#     means <- rep(-Inf, K_s)
#     for (k in seq_len(K_s)) {
#       idx <- which(z_s == k)
#       if (length(idx) > 0) means[k] <- mean(node_metric[idx])
#     }
#     
#     # ord_old: old labels sorted by decreasing block score (strongest first)
#     ord_old <- order(means, decreasing = TRUE, na.last = TRUE)
#     # mapping old -> new (p[old] = new)
#     p_old2new <- integer(K_s)
#     p_old2new[ord_old] <- seq_len(K_s)
#     
#     # apply to z
#     relab_z[s, ] <- p_old2new[z_s]
#     
#     # relabel kappa (if available)
#     if (!is.null(draws_kap)) {
#       if (is_kap_list) {
#         kap_s <- draws_kap[[s]]
#         if (!is.null(kap_s) && length(kap_s) > 0) {
#           kap_s <- as.matrix(kap_s)
#           if (nrow(kap_s) >= K_s && ncol(kap_s) >= K_s) {
#             relab_kap[[s]] <- kap_s[ord_old, ord_old, drop = FALSE]
#           } else {
#             relab_kap[[s]] <- kap_s  # fallback: leave as is
#           }
#         }
#       } else if (!is.null(dim(draws_kap))) {
#         # array S x K_max x K_max; we only relabel the first K_s blocks
#         kap_s <- draws_kap[s, , , drop = FALSE]
#         K_max <- dim(draws_kap)[2]
#         if (K_s <= K_max) {
#           # reorder the active sub-block; keep rest untouched
#           kap_mat <- kap_s[1, , ]
#           idx_active <- ord_old
#           kap_active <- kap_mat[idx_active, idx_active, drop = FALSE]
#           # embed back in K_max x K_max (simple approach: overwrite top-left)
#           kap_new <- kap_mat
#           kap_new[seq_len(K_s), seq_len(K_s)] <- kap_active
#           relab_kap[s, , ] <- kap_new
#         }
#       }
#     }
#     
#     # relabel psi (if available) – same idea as kappa
#     if (!is.null(draws_psi)) {
#       if (is_psi_list) {
#         psi_s <- draws_psi[[s]]
#         if (!is.null(psi_s) && length(psi_s) > 0) {
#           psi_s <- as.matrix(psi_s)
#           if (nrow(psi_s) >= K_s && ncol(psi_s) >= K_s) {
#             relab_psi[[s]] <- psi_s[ord_old, ord_old, drop = FALSE]
#           } else {
#             relab_psi[[s]] <- psi_s
#           }
#         }
#       } else if (!is.null(dim(draws_psi))) {
#         psi_s <- draws_psi[s, , , drop = FALSE]
#         K_max <- dim(draws_psi)[2]
#         if (K_s <= K_max) {
#           psi_mat <- psi_s[1, , ]
#           idx_active <- ord_old
#           psi_active <- psi_mat[idx_active, idx_active, drop = FALSE]
#           psi_new <- psi_mat
#           psi_new[seq_len(K_s), seq_len(K_s)] <- psi_active
#           relab_psi[s, , ] <- psi_new
#         }
#       }
#     }
#     
#     perms[[s]] <- p_old2new
#   }
#   
#   mcmc_out$z <- relab_z
#   if (!is.null(draws_kap)) mcmc_out$kappa <- relab_kap
#   if (!is.null(draws_psi)) mcmc_out$psi   <- relab_psi
#   
#   # perm is a list because K_s can vary; downstream code here never uses it
#   mcmc_out$perm <- perms
#   
#   attr(mcmc_out, "relabel_ordering") <- "WST"
#   attr(mcmc_out, "relabel_score")    <- score
#   mcmc_out
# }


# ---------- Helpers: dyads, logistic, safe choose(K,3) sampler ----------
logistic <- function(x) 1/(1+exp(-x))

build_dyad_index <- function(n) {
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  cbind(i = idx[,1], j = idx[,2])
}

# ---------- Block-level P matrices ----------
# Posterior-mean block win probs for OSBM (from relabeled draws)
block_P_osbm <- function(out_relab, regime, K) {
  if (regime == "WST") {
    # out_relab$psi: [iter, K, K]
    psi_bar <- apply(out_relab$psi, c(2,3), mean)
    P <- matrix(0.5, K, K)
    off <- row(P) != col(P)
    P[off] <- logistic(psi_bar[off])
    diag(P) <- 0.5
    return(P)
  } else {
    # regime == "SST"; out_relab$psi: [iter, K-1]
    psi_bar <- colMeans(out_relab$psi)
    P <- matrix(0.5, K, K)
    for (k in 1:K) for (l in 1:K) if (k != l) {
      d <- abs(k - l); s <- sign(l - k)
      P[k, l] <- logistic(s * psi_bar[d])
    }
    diag(P) <- 0.5
    return(P)
  }
}

# Empirical P from data and a partition (for DC-SBM or PPC-style data-level check)
block_P_empirical <- function(A, z_hat, K) {
  A <- as.matrix(A)
  num <- den <- matrix(0, K, K)
  n <- nrow(A)
  for (i in 1:(n-1)) for (j in (i+1):n) {
    ki <- z_hat[i]; lj <- z_hat[j]; if (ki == lj) next
    Nij <- A[i,j] + A[j,i]; if (Nij == 0) next
    num[ki, lj] <- num[ki, lj] + A[i,j]
    num[lj, ki] <- num[lj, ki] + A[j,i]
    den[ki, lj] <- den[ki, lj] + Nij
    den[lj, ki] <- den[lj, ki] + Nij
  }
  P <- ifelse(den > 0, num/den, 0.5)
  diag(P) <- 0.5
  P
}

# ---------- Triplet conformity on a permutation pi ----------
.triplet_conf_WST <- function(P, a, b, c, pi) {
  # assume a<b<c in position order of pi
  pab <- P[pi[a], pi[b]]; pbc <- P[pi[b], pi[c]]; pac <- P[pi[a], pi[c]]
  !( (pab >= .5) && (pbc >= .5) && (pac < .5) )
}

.triplet_conf_SST <- function(P, a, b, c, pi) {
  pab <- P[pi[a], pi[b]]; pbc <- P[pi[b], pi[c]]; pac <- P[pi[a], pi[c]]
  pac >= max(pab, pbc)
}

# Exact or Monte Carlo estimator of TOC (and SE)
estimate_TOC <- function(P, type = c("WST","SST"), pi = NULL, m = NULL) {
  type <- match.arg(type); K <- nrow(P)
  if (is.null(pi)) pi <- 1:K
  if (K < 3) return(list(est = NA_real_, se = NA_real_, m = 0L, exact = TRUE))
  
  # all triples or sampled triples
  all_tr <- t(combn(K, 3))
  n_all <- nrow(all_tr)
  if (is.null(m) || m >= n_all) {
    m_idx <- seq_len(n_all)
    exact <- TRUE
  } else {
    m_idx <- sample.int(n_all, m)
    exact <- FALSE
  }
  
  ok <- logical(length(m_idx))
  if (type == "WST") {
    for (t in seq_along(m_idx)) {
      a <- all_tr[m_idx[t], 1]; b <- all_tr[m_idx[t], 2]; c <- all_tr[m_idx[t], 3]
      ok[t] <- .triplet_conf_WST(P, a, b, c, pi)
    }
  } else {
    for (t in seq_along(m_idx)) {
      a <- all_tr[m_idx[t], 1]; b <- all_tr[m_idx[t], 2]; c <- all_tr[m_idx[t], 3]
      ok[t] <- .triplet_conf_SST(P, a, b, c, pi)
    }
  }
  
  p_hat <- mean(ok)
  se <- sqrt(p_hat*(1 - p_hat) / length(ok))
  list(est = p_hat, se = se, m = length(ok), exact = exact)
}

# ---------- Directional scores at dyad level (log loss, Brier) ----------
direction_scores <- function(A, z_hat, Pblk) {
  A <- as.matrix(A); n <- nrow(A)
  y <- p <- NULL; idx <- 0L
  for (i in 1:(n-1)) for (j in (i+1):n) {
    Nij <- A[i,j] + A[j,i]
    if (Nij == 0L) next
    ki <- z_hat[i]; lj <- z_hat[j]
    pij <- if (ki == lj) 0.5 else Pblk[ki, lj]
    yij <- as.integer(A[i,j] > A[j,i])
    idx <- idx + 1L
    y[idx] <- yij
    p[idx] <- pij
  }
  if (idx == 0L) return(c(logloss = NA_real_, brier = NA_real_))
  eps <- 1e-12; p <- pmin(pmax(p, eps), 1 - eps)
  logloss <- -mean(y*log(p) + (1 - y)*log(1 - p))
  brier   <- mean((y - p)^2)
  c(logloss = logloss, brier = brier)
}

# ---------- WST: test if a permutation exists (DAG) and get topological order ----------
wst_topo_from_P <- function(P, tol = 0) {
  K <- nrow(P)
  # Edge k->l if P[k,l] > 0.5 + tol
  B <- (P > (0.5 + tol)) * 1L
  diag(B) <- 0L
  indeg <- colSums(B)
  Q <- which(indeg == 0L)
  order <- integer(0)
  while (length(Q) > 0) {
    v <- Q[1]; Q <- Q[-1]
    order <- c(order, v)
    outs <- which(B[v,] == 1L)
    for (w in outs) {
      indeg[w] <- indeg[w] - 1L
      if (indeg[w] == 0L) Q <- c(Q, w)
    }
  }
  is_dag <- length(order) == K
  list(is_dag = is_dag, pi = if (is_dag) match(seq_len(K), order) else 1:K)
}

# ---------- SST: violation score and greedy adjacent-swap improvement ----------
sst_violation_score <- function(P, pi) {
  K <- nrow(P); V <- 0L; Ttot <- 0L
  for (i in 1:(K-1)) for (j in (i+1):K) {
    for (k in 1:K) if (k != i && k != j) {
      # rows: P[pi[i], pi[k]] >= P[pi[j], pi[k]]
      V <- V + as.integer(P[pi[i], pi[k]] < P[pi[j], pi[k]])
      # cols: P[pi[k], pi[i]] <= P[pi[k], pi[j]]
      V <- V + as.integer(P[pi[k], pi[i]] > P[pi[k], pi[j]])
      Ttot <- Ttot + 2L
    }
  }
  list(viol = V, rate = if (Ttot>0) V/Ttot else NA_real_)
}

borda_init <- function(P) {
  # higher row-sum => stronger
  order(order(-rowSums(P)))
}

sst_greedy_order <- function(P, max_iter = NULL) {
  K <- nrow(P); if (is.null(max_iter)) max_iter <- 5L*K*K
  pi <- borda_init(P)
  best <- sst_violation_score(P, pi)$viol
  it <- 0L; improved <- TRUE
  while (improved && it < max_iter) {
    improved <- FALSE; it <- it + 1L
    for (pos in 1:(K-1)) {
      pi2 <- pi; pi2[c(pos, pos+1)] <- pi2[c(pos+1, pos)]
      v2 <- sst_violation_score(P, pi2)$viol
      if (v2 < best) { pi <- pi2; best <- v2; improved <- TRUE }
    }
  }
  out <- sst_violation_score(P, pi)
  list(pi = pi, viol = out$viol, rate = out$rate)
}

#simulating a network
simulate_osbm <- function(n, K, z, eta, kappa, psi, regime = c("WST","SST")) {
  regime <- match.arg(regime)
  A <- Matrix::Matrix(0, n, n, sparse = TRUE)
  
  # helper: signed logit and safe index for WST psi
  psi_pair <- function(k, l, psi_mat) psi_mat[pmin(k,l), pmax(k,l)]
  sgn <- function(x) if (x > 0) 1L else if (x < 0) -1L else 0L
  
  for (i in 1:(n-1)) for (j in (i+1):n) {
    lam <- eta[i] * eta[j] * kappa[z[i], z[j]]         # symmetric kappa assumed
    N   <- rpois(1L, lam); if (N == 0L) next
    
    if (z[i] == z[j]) {
      p <- 0.5
      fwd <- rbinom(1L, N, p)
      # within-block: direction meaningless but keep symmetry
      A[i,j] <- fwd
      A[j,i] <- N - fwd
    } else if (regime == "WST") {
      # correct WST branch
      th   <- psi[pmin(z[i], z[j]), pmax(z[i], z[j])]   # ψ_{lower,higher} ≥ 0
      pFwd <- plogis(th)
      fwd  <- rbinom(1L, N, pFwd)
      
      if (z[i] < z[j]) {
        A[i,j] <- fwd;     A[j,i] <- N - fwd
      } else {
        A[i,j] <- N - fwd; A[j,i] <- fwd
      }
      
    } else { # SST
      # correct SST branch
      d    <- abs(z[i] - z[j])
      pFwd <- plogis(psi[d])            # forward = lower block -> higher block
      fwd  <- rbinom(1L, N, pFwd)
      
      if (z[i] < z[j]) {                # i’s block is lower
        A[i,j] <- fwd;     A[j,i] <- N - fwd
      } else {                          # j’s block is lower
        A[i,j] <- N - fwd; A[j,i] <- fwd
      }
      
    }
  }
  A
}
#generating a KxK matrix psi ~SST
make_psi_sst <- function(K, strength = c("strong", "weak"), manual_scale = NULL) {
  
  mode <- match.arg(strength)
  D <- K - 1  # Number of distances
  
  if (D < 1) return(numeric(0))
  
  # 1. Determine Increment parameters (Delta)
  #    psi_d = sum_{i=1}^d delta_i
  if (!is.null(manual_scale)) {
    d_mean <- manual_scale[1]
    d_sd   <- if (length(manual_scale) > 1) manual_scale[2] else d_mean / 4
  } else {
    if (mode == "strong") {
      # Strong: Each step up the ladder is a decisive advantage.
      # e.g., d=1 => logit ~ 0.8 (p=0.69)
      #       d=4 => logit ~ 3.2 (p=0.96)
      d_mean <- 0.8
      d_sd   <- 0.2
    } else {
      # Weak: Steps are very subtle. Hierarchy is "flat".
      # e.g., d=1 => logit ~ 0.15 (p=0.54)
      #       d=4 => logit ~ 0.60 (p=0.65) -> HARD for models to detect!
      d_mean <- 0.15
      d_sd   <- 0.05
    }
  }
  
  # 2. Draw positive increments (Truncated Normal implicitly via abs)
  #    We use abs() to ensure strictly non-decreasing, though with these means
  #    negative draws are unlikely.
  deltas <- abs(rnorm(D, mean = d_mean, sd = d_sd))
  
  # 3. Cumulative Sum to get Distance Parameters (psi_1, psi_2, ...)
  psi_vec <- cumsum(deltas)
  
  return(psi_vec)
}
#generating a KxK matrix psi ~WST
make_psi_wst <- function(K, strength = c("strong", "weak"), manual_scale = NULL) {
  
  PSI <- matrix(0, K, K)
  n_upper <- K * (K - 1) / 2
  
  # Define Logit Ranges
  if (!is.null(manual_scale)) {
    lo <- manual_scale[1]
    hi <- manual_scale[2]
  } else {
    mode <- match.arg(strength)
    if (mode == "strong") {
      # Strong: Probabilities ~ 0.80 to 0.99
      # Logits: 1.4 to 4.5
      lo <- 1.4
      hi <- 4.5
    } else {
      # Weak: Probabilities ~ 0.52 to 0.65
      # Logits: 0.1 to 0.6
      # This is the "Hard" regime for WST vs SST
      lo <- 0.1
      hi <- 0.6
    }
  }
  
  # Draw random logits for upper triangle (i < j)
  # We use runif to simulate "unstructured" WST (no distance dependence)
  # This distinguishes it mathematically from SST (which is monotonic).
  psi_vals <- runif(n_upper, min = lo, max = hi)
  
  # Fill Upper Triangle
  PSI[upper.tri(PSI)] <- psi_vals
  
  # Force Skew-Symmetry (Lower triangle = -Upper)
  # PSI[j, i] = -PSI[i, j]
  PSI <- PSI - t(PSI)
  
  return(PSI)
}
#generating a KxK matrix kappa which is symmetric
make_kappa <- function(K, mean = 2, var = 1) {
  if (mean <= 0) stop("'mean' must be > 0 for a Gamma distribution.")
  if (var  <= 0) stop("'var' must be > 0 for a Gamma distribution.")
  
  shape <- mean^2 / var   # alpha
  rate  <- mean / var     # beta (R's rate parametrization)
  
  KAP <- matrix(0, K, K)
  for (k in 1:K) {
    for (l in k:K) {
      KAP[k, l] <- rgamma(1, shape = shape, rate = rate)
    }
  }
  KAP[lower.tri(KAP)] <- t(KAP)[lower.tri(KAP)]
  KAP
}

# ---- OSBM wrappers (assumes your modular_osbm_sampler follows the new WST/SST conventions) ----
# fit_osbm <- function(A, K, truth, regime = c("WST","SST"),
#                      n_iter=3500, burn=500, thin=2, verbose=TRUE) {
#   regime <- match.arg(regime)
#   free <- c("z","kappa","psi","eta")
#   psi_constraint <- if (regime=="WST") "WST" else "SST"
#   out <- modular_osbm_sampler(
#     A=A, K=K, truth=truth, free=free,
#     n_iter=n_iter, burn=burn, thin=thin,
#     verbose=verbose, psi_constraint=psi_constraint
#   )
#   out
# }

# ---- DC–SBM wrapper (calls your provided Gibbs) ----
# fit_dcsbm <- function(A, K,
#                       iters=3500, burn=500, thin=2, verbose=250, seed=NULL) {
#   fit <- fit_dirdcsbm_gibbs(
#     A, K,
#     iters = iters, burn_in = burn, thin = thin,
#     priors = list(a_out=1, b_out=1, a_in=1, b_in=1, a_lambda=1, b_lambda=1, alpha=rep(1,K)),
#     normalize_each_iter = TRUE,
#     verbose = verbose,
#     seed = seed
#   )
#   fit
# }
