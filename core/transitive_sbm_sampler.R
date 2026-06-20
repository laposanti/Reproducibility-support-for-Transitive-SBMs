

Rcpp::sourceCpp("core/block_totals_for_poisson_cpp.cpp")
source("helper_folder/mixing_moves.R")
#' Compute Pointwise Log-Likelihood Matrix for Variable-K OSBM
#' 
#' @param A Adjacency matrix (n x n)
#' @param mcmc_out List returned by modular_osbm_sampler
#' @param regime "WST" or "SST"
#' @param dyad_index Optional: matrix with columns (row, col) for upper tri dyads
#' 
#' @return S x N_dyads matrix suitable for loo::loo()
E_by_block_excluding_i <- function(i, z, eta, K) {
  E <- tapply(eta, factor(z, levels = seq_len(K)), sum)
  E[is.na(E)] <- 0
  ki <- z[i]
  if (!is.na(ki) && ki >= 1L && ki <= K) E[ki] <- E[ki] - eta[i]
  E
}

pack_state_sst <- function(z, kappa, psi) {
  z <- as.integer(z)
  K_old <- nrow(kappa)
  stopifnot(max(z) <= K_old, min(z) >= 1L)
  
  keep  <- sort(unique(z))          # occupied labels in increasing order
  K_new <- length(keep)
  
  # old -> new map
  map <- integer(K_old)
  map[keep] <- seq_len(K_new)
  z_new <- map[z]
  
  # shrink kappa
  kappa_new <- kappa[keep, keep, drop = FALSE]
  
  # Recompute psi via its increment (delta) representation so that when
  # interior labels are removed, distances are reindexed instead of
  # blindly truncating/padding (which can change distance semantics).
  if (K_new <= 1L) {
    psi_new <- numeric(0)
  } else {
    D_old <- max(K_old - 1L, 0L)
    if (D_old == 0L) {
      psi_old <- numeric(0)
    } else {
      psi_old <- psi
      if (length(psi_old) < D_old) {
        last <- if (length(psi_old)) tail(psi_old, 1L) else 0
        psi_old <- c(psi_old, rep(last, D_old - length(psi_old)))
      } else {
        psi_old <- psi_old[seq_len(D_old)]
      }
    }
    
    if (D_old > 0L) {
      delta_old <- numeric(D_old)
      delta_old[1] <- max(psi_old[1], 0)
      if (D_old >= 2L) {
        delta_old[2:D_old] <- pmax(psi_old[2:D_old] - psi_old[1:(D_old - 1L)], 0)
      }
      psi_old <- cumsum(delta_old)
    }
    
    D_new <- K_new - 1L
    psi_new <- numeric(D_new)
    if (D_old > 0L) {
      for (d in seq_len(D_new)) {
        old_dists <- keep[(1L + d):K_new] - keep[1L:(K_new - d)]
        old_dists <- pmax(pmin(old_dists, D_old), 1L)
        psi_new[d] <- mean(psi_old[old_dists])
      }
      psi_new <- cummax(psi_new)
    }
  }
  
  changed <- !(K_new == K_old && identical(keep, seq_len(K_old)))
  
  list(
    z = as.integer(z_new),
    kappa = kappa_new,
    psi = as.numeric(psi_new),
    K = K_new,
    changed = changed
  )
}



assert_sst_invariants <- function(z, kappa, psi, where = "") {
  K <- nrow(kappa)
  u <- sort(unique(z))
  
  missing <- setdiff(seq_len(K), u)
  extra   <- setdiff(u, seq_len(K))
  
  if (length(extra) || length(missing)) {
    stop(sprintf(
      "SST invariant failed%s\n  K=%d\n  unique(z)=%s\n  missing=%s\n  extra=%s\n  range(z)=[%s,%s]",
      if (nzchar(where)) paste0(" @ ", where) else "",
      K,
      paste(u, collapse=","),
      paste(missing, collapse=","),
      paste(extra, collapse=","),
      min(z), max(z)
    ), call. = FALSE)
  }
  
  if (max(z) != K) {
    stop(sprintf("SST invariant failed%s: max(z)=%d != K=%d",
                 if (nzchar(where)) paste0(" @ ", where) else "", max(z), K),
         call. = FALSE)
  }
  
  if (length(psi) != max(K - 1L, 0L)) {
    stop(sprintf("SST invariant failed%s: length(psi)=%d != K-1=%d",
                 if (nzchar(where)) paste0(" @ ", where) else "", length(psi), max(K-1L,0L)),
         call. = FALSE)
  }
  
  invisible(TRUE)
}


get_principled_hypers_v4 <- function(
    A,
    # Target prior weight for a "typical" block pair:
    # w_kappa_star = prior weight in posterior mean of kappa
    w_kappa_star   = 0.2,
    # Prior mean intensity: still shrunk version of global mean
    sparsity_factor  = 5.0,
    # Expected number of blocks under the partition prior.
    # If NULL, we use a simple heuristic sqrt(n).
    K0              = NULL,
    # Directional priors
    dominance_target = 0.75,  # for SST: target extreme dominance prob
    N0_psi           = 13     # for WST: prior effective comparisons per block pair
) {
  if (w_kappa_star <= 0 || w_kappa_star >= 1) {
    stop("w_kappa_star must be in (0,1).")
  }
  
  n <- nrow(A)
  if (is.null(n) || n != ncol(A)) {
    stop("A must be a square adjacency / count matrix.")
  }
  
  ## ---- Symmetrised counts for summaries ----
  N_mat <- A + t(A)
  diag(N_mat) <- 0
  mean_edges <- mean(N_mat[upper.tri(N_mat)])
  if (!is.finite(mean_edges) || mean_edges < 1e-8) {
    mean_edges <- 1e-8
  }
  
  ## ---- Degrees for eta prior ----
  degrees <- rowSums(A) + colSums(A)
  deg_mean <- mean(degrees)
  deg_sd   <- sd(degrees)
  cv_deg   <- deg_sd / deg_mean
  if (!is.finite(cv_deg) || cv_deg < 0.1) {
    cv_deg <- 0.1
  }
  
  ## ---- Typical block-pair exposure T_typ ----
  # Use prior expected K0 if given, otherwise a simple heuristic.
  if (is.null(K0)) {
    K0 <- sqrt(n)   # you can replace this with E[K_n] under Gnedin if you have it
  }
  # Typical block size ~ n / K0  ->  T_typ ~ (n / K0)^2
  T_typ <- (n / K0)^2
  
  ## ---- VOLUME: calibrate via shrinkage weight w_kappa_star ----
  # Prior mean kappa from global intensity and sparsity_factor
  mu_kappa <- mean_edges / sparsity_factor
  
  # Solve w_kappa_star = b / (b + T_typ)  for b:
  b_kappa <- (w_kappa_star / (1 - w_kappa_star)) * T_typ
  # Then a_kappa from mean = a / b = mu_kappa
  a_kappa <- mu_kappa * b_kappa
  
  # Derived CV_kappa (for diagnostics / curiosity)
  CV_kappa_derived <- 1 / sqrt(a_kappa)
  
  ## ---- DEGREES: match empirical CV of degrees ----
  a_eta <- 1 / (cv_deg^2)
  b_eta <- a_eta
  
  ## ---- DIRECTIONAL (WST): sigma0 via effective sample size N0_psi ----
  # N0_psi = 4 / sigma0^2  ->  sigma0 = sqrt(4 / N0_psi)
  sigma0 <- sqrt(4 / N0_psi)
  mu0    <- 0
  
  ## ---- DIRECTIONAL (SST): tau0 via dominance_target and K0 ----
  # Use distance ~ (K0 - 1) instead of fixed 10
  logit_target <- qlogis(dominance_target)
  D_dist <- max(1, K0 - 1)  # typical extreme distance in the order
  tau0 <- logit_target / (D_dist * sqrt(2 / pi))
  
  list(
    # volume
    a_kappa = a_kappa,
    b_kappa = b_kappa,
    mu_kappa = mu_kappa,
    CV_kappa = CV_kappa_derived,
    w_kappa_star = w_kappa_star,
    T_typ = T_typ,
    
    # degrees
    a_eta = a_eta,
    b_eta = b_eta,
    cv_deg = cv_deg,
    
    # directional
    tau0 = tau0,
    sigma0 = sigma0,
    mu0 = mu0,
    
    # for reference
    K0_used = K0,
    sparsity_factor = sparsity_factor,
    dominance_target = dominance_target,
    N0_psi = N0_psi
  )
}

# slots adjacent to current block k_cur; radius=1 -> {k_cur, k_cur+1}
.adjacent_slots <- function(k_cur, K, radius = 1L) {
  stopifnot(radius >= 1L)
  lo <- k_cur - (radius - 1L)
  hi <- k_cur + radius
  seq.int(max(1L, lo), min(K + 1L, hi))
}

# --- UPDATED: Principled Hyperparameters with Sparsity Boost ---
# get_principled_hypers <- function(A, K_expected = 5, dominance_target = 0.75) {
#   N_mat <- A + t(A); diag(N_mat) <- 0
#   mean_edges <- mean(N_mat[upper.tri(N_mat)])
#   if (mean_edges < 1e-6) mean_edges <- 1e-6
#   
#   degrees <- rowSums(A) + colSums(A)
#   cv_deg  <- sd(degrees) / mean(degrees)
#   if (is.na(cv_deg) || cv_deg < 0.1) cv_deg <- 0.1
#   
#   # VOLUME: Stricter Regularization against Shattering
#   # a_kappa = 1.5 (Flexible shape)
#   # b_kappa: Originally a/mean. 
#   # NEW: We multiply b_kappa by 5.0. 
#   # This pushes the prior mean towards 0.2 * mean_edges.
#   # This "Under-fitting" prior forces the model to only create blocks 
#   # if the data STRONGLY contradicts the sparsity assumption.
#   a_kappa <- 1.5
#   b_kappa <- (a_kappa / mean_edges) * 5.0 
#   
#   # DEGREES
#   a_eta <- 1 / (cv_deg^2)
#   b_eta <- a_eta 
#   
#   # SST tau0: Shallower slope to prevent rapid separation
#   logit_target <- qlogis(dominance_target)
#   tau0 <- logit_target / (10 * sqrt(2/pi))
#   
#   # WST sigma0: Tight regularization
#   sigma0 <- 0.45
#   
#   list(a_kappa=a_kappa, b_kappa=b_kappa, a_eta=a_eta, b_eta=b_eta, 
#        tau0=tau0, sigma0=sigma0, mu0=0)
# }

# get_principled_hypers_v2 -- REMOVED: canonical version lives in Hyper_setup.R
# (kept as comment to avoid breaking git blame)
#
# get_principled_hypers_v2 <- function(
#     A, 
#     K_expected      = 5,
#     dominance_target = 0.75,
#     kappa_frac      = 0.003,
#     sparsity_factor = 5.0
# ) { ... }



# Aggregate exact directed counts between node i and blocks (length-K guaranteed)
# c_plus[ℓ] = total "i -> j" wins over j in block ℓ
# N_tot [ℓ] = total matches vs block ℓ
# counts_by_block_exact <- function(i, A, z, i_idx, j_idx, N_edge, edge_by_node, K = max(z)) {
#   c_plus <- numeric(K)
#   N_tot  <- numeric(K)
#   e_list <- edge_by_node[[i]]
#   if (!length(e_list)) return(list(c_plus = c_plus, N_tot = N_tot))
#   for (e in e_list) {
#     u <- i_idx[e]; v <- j_idx[e]
#     j <- if (i == u) v else u
#     ell <- z[j]
#     if (ell < 1L || ell > K) next  # very defensive
#     n_e <- N_edge[e]
#     a_uv <- A[u, v]
#     a_i_to_j <- if (i == u) a_uv else (n_e - a_uv)  # orient as i -> j
#     c_plus[ell] <- c_plus[ell] + a_i_to_j
#     N_tot [ell] <- N_tot [ell] + n_e
#   }
#   list(c_plus = c_plus, N_tot = N_tot)
# }
Rcpp::sourceCpp("core/counts_by_block_exact_cpp.cpp")

# Sum of neighbour etas by block for node i, using observed dyads only.
# Returns numeric vector length K: S_eta_blk[ell] = sum_{j in block ell, (i,j) observed} eta[j]
neighbor_eta_by_block <- function(i, z, eta, i_idx, j_idx, edge_by_node, K) {
  out <- numeric(K)
  e_list <- edge_by_node[[i]]
  if (!length(e_list)) return(out)
  
  for (e in e_list) {
    u <- i_idx[e]; v <- j_idx[e]
    j <- if (i == u) v else u
    ell <- z[j]
    if (ell >= 1L && ell <= K) out[ell] <- out[ell] + eta[j]
  }
  out
}


# 1) PG draws -- check again!
draw_pg_latents <- function(i_idx, j_idx, N_edge, z, psi, mode = c("pair","distance")) {
  mode <- match.arg(mode)
  m <- length(N_edge)
  if (m == 0L) return(numeric(0))
  c_e <- numeric(m)
  
  if (mode == "pair") {
    k   <- pmin(z[i_idx], z[j_idx])
    ell <- pmax(z[i_idx], z[j_idx])
    c_e <- abs(psi[cbind(k, ell)])
  } else {
    d <- abs(z[i_idx] - z[j_idx])
    maskd <- (d > 0L)
    if (any(maskd)) {
      md <- max(d[maskd])
      if (md > length(psi)) {
        stop(sprintf("[draw_pg_latents] max(d)=%d > length(psi)=%d. (This is the REAL bug.)",
                     md, length(psi)), call. = FALSE)
      }
      c_e[maskd] <- abs(psi[d[maskd]])
    }
    c_e[!maskd] <- 0
  }
  
  omega <- numeric(m)
  maskN <- (N_edge > 0)
  if (any(maskN)) {
    nn <- as.integer(sum(maskN))
    z_sub <- c_e[maskN]
    h_sub <- N_edge[maskN]
    omega[maskN] <- BayesLogit::rpg(num = nn, h = h_sub, z = z_sub)
  }
  omega
}


# --- insert a new block label at slot r; shift existing labels >= r by +1 ----
insert_block_labels <- function(z, K, r) {
  stopifnot(is.integer(K) || is.numeric(K), K >= 1L,
            is.integer(r) || is.numeric(r), r >= 1L, r <= (K + 1L))
  z_new <- ifelse(z < r, z, z + 1L)
  list(z_new = as.integer(z_new))
}



# Pool N, Y, and Ω for node i against each block ℓ using realized edge-level ω
pool_node_block_stats <- function(i, z, i_idx, j_idx, N_edge, omega_e, A, K, edge_by_node) {
  N_blk     <- numeric(K)
  Y_blk     <- numeric(K)   # unsigned, i->j orientation; sign is applied later by slot
  Omega_blk <- numeric(K)   # sum of realized PG ω over edges from i to block ℓ
  
  e_list <- edge_by_node[[i]]
  if (!length(e_list)) return(list(N_blk = N_blk, Y_blk = Y_blk, Omega_blk = Omega_blk))
  
  for (e in e_list) {
    u <- i_idx[e]; v <- j_idx[e]
    j <- if (i == u) v else u
    ell <- z[j]
    if (ell < 1L || ell > K) next
    
    n_e  <- N_edge[e]
    a_uv <- A[u, v]
    a_i_to_j <- if (i == u) a_uv else (n_e - a_uv)
    
    N_blk    [ell] <- N_blk    [ell] + n_e
    Y_blk    [ell] <- Y_blk    [ell] + (a_i_to_j - 0.5 * n_e)
    om <- omega_e[e]
    if (!is.finite(om)) om <- 0
    Omega_blk[ell] <- Omega_blk[ell] + om
    
  }
  list(N_blk = N_blk, Y_blk = Y_blk, Omega_blk = Omega_blk)
}



# ====== NEW-BLOCK: collapsed κ piece (independent of slot r) ===============
log_kappa_new_collapsed <- function(i, z, eta, N_blk, a_kappa, b_kappa) {
  K <- max(z)
  E_l <- tapply(eta[-i], factor(z[-i], levels = seq_len(K)), sum); E_l[is.na(E_l)] <- 0
  eta_i <- eta[i]
  lp <- 0
  for (l in seq_len(K)) {
    dR <- N_blk[l]; dT <- eta_i * E_l[l]
    if (dR > 0 || dT > 0) {
      lp <- lp +
        lgamma(a_kappa + dR) - lgamma(a_kappa) +   # shape ratio
        a_kappa * log(b_kappa) -                   # include constant for stability
        (a_kappa + dR) * log(b_kappa + dT)
    }
  }
  lp
}

.log_crp_weights <- function(n_minus, K, alpha0, d = 0) {
  w_exist <- log(pmax(n_minus - d, 0))
  w_exist[!is.finite(w_exist)] <- -Inf
  w_new   <- log(alpha0 + d*K)
  list(exist = w_exist, new = w_new)
}

# ----------------------------------------------------------
# Helper: safe θ(s, d) even if psi_vec has length 0
#   - returns 0 when psi_vec is empty (neutral direction)
# ----------------------------------------------------------
theta_at <- function(s, d, psi_vec) {
  if (!length(psi_vec)) return(0.0)
  d_use <- min(d, length(psi_vec))
  s * psi_vec[d_use]
}

update_eta_observed <- function(eta, z, kappa,
                                a_eta, b_eta,
                                i_idx, j_idx, N_edge, edge_by_node) {
  n <- length(eta)
  K <- nrow(kappa)
  eta_new <- eta
  
  # shape increments: total observed volume incident to each node
  # (each dyad contributes N_edge to both endpoints)
  G_i <- numeric(n)
  for (e in seq_along(N_edge)) {
    u <- i_idx[e]; v <- j_idx[e]
    G_i[u] <- G_i[u] + N_edge[e]
    G_i[v] <- G_i[v] + N_edge[e]
  }
  
  for (i in seq_len(n)) {
    k_i <- z[i]
    if (is.na(k_i) || k_i < 1L || k_i > K) next
    
    e_list <- edge_by_node[[i]]
    rate_i <- b_eta
    if (length(e_list)) {
      for (e in e_list) {
        u <- i_idx[e]; v <- j_idx[e]
        j <- if (i == u) v else u
        k_j <- z[j]
        if (is.na(k_j) || k_j < 1L || k_j > K) next
        rate_i <- rate_i + kappa[k_i, k_j] * eta_new[j]
      }
    }
    rate_i <- max(rate_i, 1e-12)
    
    eta_new[i] <- rgamma(1, shape = a_eta + G_i[i], rate = rate_i)
  }
  
  eta_new
}


# --------------------------------------------------------------
# Gamma–Poisson marginal log-likelihood piece, per (R, T)
# f(R,T) = lgamma(a+R) - lgamma(a) + a*log b - (a+R)*log(b+T)
# --------------------------------------------------------------
# log marginal for Gamma-Poisson:
# kappa ~ Gamma(a_kappa, b_kappa) with "rate" parameterisation
# y | kappa ~ Poisson(kappa * T)  (T >= 0)
# returns log ∫ p(y | kappa,T) p(kappa) dkappa  up to exact constants (this is exact).
gp_marginal <- function(R, T, a_kappa, b_kappa) {
  if (any(T < 0)) stop("gp_marginal: T must be non-negative.")
  if (any(R < 0)) stop("gp_marginal: R must be non-negative.")
  # handle T=0 safely: if T=0, likelihood is Poisson with mean 0, so only R=0 allowed
  out <- rep(-Inf, length(R))
  ok0 <- (T == 0)
  out[ok0 & (R == 0)] <- 0.0  # log(1)
  ok <- !ok0
  if (any(ok)) {
    out[ok] <- lgamma(a_kappa + R[ok]) - lgamma(a_kappa) +
      a_kappa * log(b_kappa) - (a_kappa + R[ok]) * log(b_kappa + T[ok])
  }
  out
}

# AFTER (correct)
gp_delta <- function(R0, T0, r_add, t_add, a, b) {
  gp_marginal(R0 + r_add, T0 + t_add, a, b) - gp_marginal(R0, T0, a, b)
}

# Constant directional baseline for node i (independent of r)
dir_const_node <- function(i, i_idx, j_idx, N_edge, A, edge_by_node) {
  e_list <- edge_by_node[[i]]
  if (!length(e_list)) return(0)
  u <- i_idx[e_list]; v <- j_idx[e_list]
  n <- N_edge[e_list]
  a_uv <- A[cbind(u, v)]
  # orient "i -> j" consistently with counts_by_block_exact
  a_i_to_j <- ifelse(i == u, a_uv, n - a_uv)
  sum(lchoose(n, a_i_to_j) - n * log(2))
}




# ====== shrink: remove empty blocks and re-pack parameters ==================
remove_empty_blocks <- function(z, kappa, psi, regime = c("WST","SST")) {
  regime <- match.arg(regime)
  K <- nrow(kappa)
  
  if (regime == "SST") {
    # SST: use your packer (ordering-preserving) as the only authority
    if (length(psi) != max(K - 1L, 0L)) {
      stop(sprintf("[remove_empty_blocks] SST: length(psi)=%d but K-1=%d (K=%d)",
                   length(psi), max(K-1L,0L), K), call. = FALSE)
    }
    packed <- pack_state_sst(z, kappa, psi)
    return(list(z = packed$z, kappa = packed$kappa, psi = packed$psi,
                K = packed$K, changed = packed$changed))
  }
  
  # WST case (matrix psi)
  nz <- tabulate(z, nbins = K)
  keep <- which(nz > 0L)
  if (length(keep) == K) {
    return(list(z = z, kappa = kappa, psi = psi, K = K, changed = FALSE))
  }
  
  p <- integer(K); p[keep] <- seq_along(keep)
  z_new <- p[z]
  kap2  <- kappa[keep, keep, drop = FALSE]
  psi2  <- psi[keep, keep, drop = FALSE]
  
  list(z = z_new, kappa = kap2, psi = psi2, K = nrow(kap2), changed = TRUE)
}

# ====== block-wise edge      summaries for node i (same as in your z-updaters) ===
summarize_node_vs_blocks <- function(i, A, z) {
  e_out <- A[i, ];  e_in <- A[, i];  Nij <- e_out + e_in
  K <- max(z)
  N_blk <- A_out_blk <- numeric(K)
  if (any(Nij > 0)) {
    ell <- z[which(Nij > 0)]
    N_blk_tmp     <- tapply(Nij[Nij > 0], ell, sum)
    A_out_blk_tmp <- tapply(e_out[Nij > 0], ell, sum)
    N_blk[as.integer(names(N_blk_tmp))]         <- N_blk_tmp
    A_out_blk[as.integer(names(A_out_blk_tmp))] <- A_out_blk_tmp
  }
  list(N_blk = N_blk, A_out_blk = A_out_blk)
}


.softsample <- function(logw) {
  if (!any(is.finite(logw))) return(sample.int(length(logw), 1L))
  p <- exp(logw - lse(logw))
  p[!is.finite(p)] <- 0
  sample.int(length(logw), 1L, prob = p)
}

# block_totals_for_poisson_legacy <- function(z, eta, i_idx, j_idx, N_edge, K) {
#   p <- pmin(z[i_idx], z[j_idx]); q <- pmax(z[i_idx], z[j_idx])
#   p <- factor(p, levels = 1:K);  q <- factor(q, levels = 1:K)
#   Rkl <- as.matrix(xtabs(N_edge ~ p + q, drop.unused.levels = FALSE))
#   E_k   <- tapply(eta,   factor(z, levels = 1:K), sum);   E_k[is.na(E_k)] <- 0
#   eta2k <- tapply(eta^2, factor(z, levels = 1:K), sum); eta2k[is.na(eta2k)] <- 0
#   Tkl <- outer(E_k, E_k, `*`); diag(Tkl) <- pmax((E_k^2 - eta2k)/2, 0)
#   list(Rkl = Rkl, Tkl = Tkl, E_k = E_k, eta2k = eta2k)
# }
# block_totals_for_poisson <- function(z, eta, i_idx, j_idx, N_edge, K) {
#   # assume cluster labels are in {1, ..., K} or NA
#   z <- as.integer(z)
#   
#   ## 1. E_k and eta2k via rowsum (NA groups ignored)
#   E_k   <- numeric(K)
#   eta2k <- numeric(K)
#   
#   valid_z <- !is.na(z)
#   if (any(valid_z)) {
#     # sums of eta by z
#     tmp  <- rowsum(eta[valid_z], z[valid_z], reorder = FALSE)
#     g    <- as.integer(rownames(tmp))           # group labels actually present
#     E_k[g] <- as.numeric(tmp)
#     
#     # sums of eta^2 by z
#     tmp2   <- rowsum((eta[valid_z])^2, z[valid_z], reorder = FALSE)
#     g2     <- as.integer(rownames(tmp2))
#     eta2k[g2] <- as.numeric(tmp2)
#   }
#   
#   ## 2. Rkl via integer indexing + rowsum (mirroring xtabs over p,q)
#   p <- pmin(z[i_idx], z[j_idx])
#   q <- pmax(z[i_idx], z[j_idx])
#   
#   valid_e <- !(is.na(p) | is.na(q))
#   Rkl <- matrix(0, nrow = K, ncol = K)
#   
#   if (any(valid_e)) {
#     p2   <- p[valid_e]
#     q2   <- q[valid_e]
#     # correct linear index for a KxK matrix (rows = p, cols = q, column-major)
#     idx  <- p2 + (q2 - 1L) * K
#     
#     R_vec <- numeric(K * K)
#     tmpR  <- rowsum(N_edge[valid_e], idx, reorder = FALSE)
#     gR    <- as.integer(rownames(tmpR))
#     R_vec[gR] <- as.numeric(tmpR)
#     dim(R_vec) <- c(K, K)
#     Rkl <- R_vec
#   }
#   
#   ## 3. Tkl as in original (but using tcrossprod)
#   Tkl <- tcrossprod(E_k)  # same as outer(E_k, E_k, `*`)
#   diag(Tkl) <- pmax((E_k^2 - eta2k) / 2, 0)
#   
#   list(Rkl = Rkl, Tkl = Tkl, E_k = E_k, eta2k = eta2k)
# }






# ---------- GN prior weights (log-scale), robust to empty blocks ----------
# n_minus: integer vector of sizes excluding i for labels 1..K (may include zeros)
# gamma_gn in (0,1)
urn_GN <- function(v_minus, gamma_GN){
  H <- length(v_minus)
  c(
    (v_minus + 1) * (sum(v_minus) - H + gamma_GN),  # existing blocks (length H)
    H^2 - H * gamma_GN                              # new block (scalar)
  )
}

.normalize_partition_prior <- function(partition_prior) {
  if (is.null(partition_prior) || !nzchar(trimws(partition_prior))) {
    return("GN")
  }

  prior <- toupper(trimws(partition_prior))
  if (prior %in% c("R-OCRP", "REVERSED_OCRP", "REVERSED OCRP",
                   "MIRRORED_OCRP", "MIRRORED OCRP")) {
    return("OCRP")
  }

  prior
}

# ---------- GN / CRP prior weights (log-scale) ----------
# n_minus: integer vector (length K) of block sizes excluding i, may contain zeros
.log_partition_weights <- function(n_minus, K, psi_hyper) {
  prior <- .normalize_partition_prior(psi_hyper$partition_prior)
  
  
  n_minus <- as.numeric(n_minus)
  
  gamma_gn <- psi_hyper$gamma_gn
  
  # work only with non-empty blocks
  pos <- which(n_minus > 0)
  exist_log <- rep(-Inf, K)
  
  if (length(pos) == 0L) {
    # no extant blocks -> only "new" event
    return(list(exist = exist_log, new = 0.0))  # log(1)
  }
  
  v_minus <- n_minus[pos]
  w_vec   <- urn_GN(v_minus, gamma_GN = gamma_gn)
  H       <- length(v_minus)
  
  w_exist <- w_vec[seq_len(H)]
  w_new   <- w_vec[H + 1L]
  
  exist_log[pos] <- ifelse(w_exist > 0, log(w_exist), -Inf)
  new_log        <- if (w_new > 0) log(w_new) else -Inf
  
  list(
    exist =  exist_log,
    new   =  new_log
  )
  
}

update_eta_all_dyads <- function(eta, z, kappa,
                                 a_eta, b_eta,
                                 i_idx, j_idx, N_edge) {
  n <- length(eta)
  K <- nrow(kappa)
  eta_new <- eta
  
  # shape increments: sum_j N_ij (zeros add nothing)
  G_i <- numeric(n)
  for (e in seq_along(N_edge)) {
    u <- i_idx[e]; v <- j_idx[e]
    G_i[u] <- G_i[u] + N_edge[e]
    G_i[v] <- G_i[v] + N_edge[e]
  }
  
  # maintain block sums sequentially for correct conditioning
  E_k <- tapply(eta_new, factor(z, levels = seq_len(K)), sum)
  E_k[is.na(E_k)] <- 0
  
  for (i in seq_len(n)) {
    k_i <- z[i]
    if (is.na(k_i) || k_i < 1L || k_i > K) next
    
    # exclude i from its own block sum in the rate
    E_k[k_i] <- E_k[k_i] - eta_new[i]
    
    rate_i <- b_eta + sum(kappa[k_i, ] * E_k)
    rate_i <- max(rate_i, 1e-12)
    
    eta_new[i] <- rgamma(1, shape = a_eta + G_i[i], rate = rate_i)
    
    # put it back (with the new value)
    E_k[k_i] <- E_k[k_i] + eta_new[i]
  }
  
  eta_new
}


update_b_kappa <- function(kappa, a_kappa, alpha_b, beta_b) {
  # kappa is a KxK matrix; only upper-triangular (incl diag) are unique
  if (!is.matrix(kappa)) stop("kappa must be a matrix")
  if (a_kappa <= 0) stop("a_kappa must be > 0")
  if (alpha_b <= 0 || beta_b <= 0) stop("alpha_b and beta_b must be > 0")
  
  idx <- upper.tri(kappa, diag = TRUE)
  k_vec <- kappa[idx]
  
  if (any(!is.finite(k_vec)) || any(k_vec < 0)) stop("kappa contains invalid entries")
  
  M <- length(k_vec)
  shape_post <- alpha_b + a_kappa * M
  rate_post  <- beta_b  + sum(k_vec)
  
  rgamma(1, shape = shape_post, rate = rate_post)
}

gn_log_weights_packed <- function(v_minus, gamma_gn) {
  v_minus <- as.numeric(v_minus)
  H <- length(v_minus)
  if (H < 1L) stop("gn_log_weights_packed: H<1 (degenerate).", call. = FALSE)
  if (!(gamma_gn > 0 && gamma_gn < 1)) stop("gamma_gn must be in (0,1).", call. = FALSE)
  
  w_exist <- (v_minus + 1) * (sum(v_minus) - H + gamma_gn)
  w_new   <- H^2 - H * gamma_gn
  
  if (any(w_exist <= 0) || w_new <= 0)
    stop("GN weights became non-positive; check gamma_gn and block sizes.", call. = FALSE)
  
  list(exist = log(w_exist), new = log(w_new))
}

# ---------------------------------------------------------------------------
# Penalised GP05 (Dirichlet ordered partition) prior log-weights
# ---------------------------------------------------------------------------
# Returns per-block and per-slot log-weights based on the GP05 composition
# structure (theta=0) with a geometric penalty gamma on new-block events.
#
# Existing-block weights (Theorem 8, theta=0):
#   Interior j < H:  w_j = v_j - alpha
#   Last     j = H:  w_H = (v_H + 1)(v_H - alpha) / v_H
#
# New-block slot weights (with gamma penalty):
#   Interior r <= H:   w_r = gamma * alpha
#   Far-right r = H+1: w_r = gamma * alpha / v_H
#
# These are unnormalised log-weights; normalisation happens at sampling time.
# Returns list(exist = numeric(H), new = numeric(H+1), per_slot = TRUE).
# The per_slot = TRUE flag tells the z-update that new-block weights are
# per-slot (no uniform slot averaging needed).
gp05_log_weights_packed <- function(v_minus, alpha_gp05, gamma_gp05) {
  v_minus <- as.numeric(v_minus)
  H <- length(v_minus)
  if (H < 1L) stop("gp05_log_weights_packed: H<1 (degenerate).", call. = FALSE)
  if (!(alpha_gp05 > 0 && alpha_gp05 < 1))
    stop("alpha_gp05 must be in (0,1).", call. = FALSE)
  if (!(gamma_gp05 > 0 && gamma_gp05 < 1))
    stop("gamma_gp05 must be in (0,1).", call. = FALSE)

  # --- existing block weights ---
  w_exist <- numeric(H)
  for (j in seq_len(H)) {
    if (j < H) {
      w_exist[j] <- v_minus[j] - alpha_gp05
    } else {
      w_exist[j] <- (v_minus[H] + 1) * (v_minus[H] - alpha_gp05) / v_minus[H]
    }
  }
  if (any(w_exist <= 0))
    stop("GP05 existing-block weights non-positive; alpha too large for block sizes.",
         call. = FALSE)

  # --- new block slot weights ---
  w_new <- numeric(H + 1L)
  for (r in seq_len(H + 1L)) {
    if (r <= H) {
      w_new[r] <- gamma_gp05 * alpha_gp05
    } else {
      w_new[r] <- gamma_gp05 * alpha_gp05 / v_minus[H]
    }
  }

  list(exist = log(w_exist), new = log(w_new), per_slot = TRUE)
}

# ---------------------------------------------------------------------------
# Ordered CRP (regenerative DP composition) prior log-weights
# ---------------------------------------------------------------------------
# Predictive probabilities for the ordered CRP associated to the
# regenerative composition p(n_1,...,n_k) = theta^k / (theta)_{(n)} * n! / (N_1...N_k)
# where N_j = n_j + ... + n_k (right cumulative sums).
#
# Existing block j:
#   p_j^old = (n_j + 1) / (theta + n) * prod_{i=1}^{j} N_i / (N_i + 1)
#
# New singleton at slot r (r = 1, ..., H+1):
#   q_r^new = theta / (theta + n) * prod_{i=1}^{r-1} N_i/(N_i+1) * 1/(N_r + 1)
#   with N_{H+1} = 0.
#
# Returns list(exist = numeric(H), new = numeric(H+1), per_slot = TRUE).
rocrp_log_weights_packed <- function(v_minus, theta_ocrp) {
  v_minus <- as.numeric(v_minus)
  H <- length(v_minus)
  if (H < 1L) stop("rocrp_log_weights_packed: H<1 (degenerate).", call. = FALSE)
  if (theta_ocrp <= 0) stop("theta_ocrp must be > 0.", call. = FALSE)

  n <- sum(v_minus)  # total count excluding node i

  # Right cumulative sums: V_j = v_j + ... + v_H
  V <- rev(cumsum(rev(v_minus)))  # V[1] = n, V[H] = v_H
  # Extend: V[H+1] = 0
  V_ext <- c(V, 0)

  # Precompute cumulative log-ratio: A_j = prod_{i=1}^{j-1} V[i]/(V[i]+1)
  # log_A[1] = 0, log_A[j+1] = log_A[j] + log(V[j]) - log(V[j]+1)
  log_A <- numeric(H + 1L)  # log_A[1] = 0
  for (j in seq_len(H)) {
    log_A[j + 1L] <- log_A[j] + log(V[j]) - log(V[j] + 1)
  }

  # Existing block weights (log scale)
  # p_j^old = (v_j + 1) / (theta + n) * A_{j+1}
  # (A_{j+1} = prod_{i=1}^{j} V_i/(V_i+1))
  lp_exist <- log(v_minus + 1) - log(theta_ocrp + n) + log_A[2:(H + 1L)]

  # New block slot weights (log scale)
  # q_r^new = theta / (theta + n) * A_r / (V_r + 1)
  # r = 1, ..., H+1
  lp_new <- numeric(H + 1L)
  for (r in seq_len(H + 1L)) {
    lp_new[r] <- log(theta_ocrp) - log(theta_ocrp + n) +
                 log_A[r] - log(V_ext[r] + 1)
  }

  list(exist = lp_exist, new = lp_new, per_slot = TRUE)
}

ocrp_log_weights_packed <- function(v_minus, theta_ocrp) {
  v_minus <- as.numeric(v_minus)
  K <- length(v_minus)
  if (K < 1L) stop("ocrp_log_weights_packed: K<1 (degenerate).", call. = FALSE)
  if (theta_ocrp <= 0) stop("theta_ocrp must be > 0.", call. = FALSE)

  n <- sum(v_minus)

  # Left cumulative sums: S_j = v_1 + ... + v_j
  S <- cumsum(v_minus)
  S_prev <- c(0, S)

  # B_j = prod_{i=j}^K S_i/(S_i+1), with B_{K+1} = 1
  log_B <- numeric(K + 1L)
  for (j in K:1L) {
    log_B[j] <- log_B[j + 1L] + log(S[j]) - log(S[j] + 1)
  }

  # Existing block weights:
  # p_j^old = (v_j + 1)/(theta+n) * prod_{i=j}^K S_i/(S_i+1)
  lp_exist <- log(v_minus + 1) - log(theta_ocrp + n) + log_B[1:K]

  # New block slot weights:
  # q_r^new = theta/(theta+n) * 1/(S_{r-1}+1) * prod_{i=r}^K S_i/(S_i+1)
  lp_new <- numeric(K + 1L)
  for (r in seq_len(K + 1L)) {
    lp_new[r] <- log(theta_ocrp) - log(theta_ocrp + n) -
      log(S_prev[r] + 1) + log_B[r]
  }

  list(exist = lp_exist, new = lp_new, per_slot = TRUE)
}

# ---------------------------------------------------------------------------
# Expected K under penalised GP05 prior (exact + heuristic)
# ---------------------------------------------------------------------------
# Exact computation via O(n^2) recursion for C_{n,m}(alpha).
# P(K=m) = gamma^{m-1} C_{n,m}(alpha) / Z_n(alpha,gamma)
# E[K] = sum_m m * P(K=m)
#
# log_q_gp05(n, m, alpha): log decrement matrix q_{alpha,0}(n:m)
.log_q_gp05 <- function(n, m, alpha) {
  if (n < 1L || m < 1L || m > n) return(-Inf)
  if (m == n) {
    # diagonal: (1-alpha)_{n-1} / (n-1)!
    if (n == 1L) return(0)  # q(1:1) = 1
    return(sum(log(seq(1 - alpha, length.out = n - 1L, by = 1))) - lgamma(n))
  }
  # off-diagonal (m < n):
  # C(n,m) * (1-alpha)_{m-1} / (n-m)_m * (n-m)*alpha / n
  lq <- lchoose(n, m)
  if (m >= 2L) lq <- lq + sum(log(seq(1 - alpha, length.out = m - 1L, by = 1)))
  # (n-m)_m = (n-m)(n-m+1)...(n-1) = Gamma(n)/Gamma(n-m)
  lq <- lq - (lgamma(n) - lgamma(n - m))
  lq <- lq + log(n - m) + log(alpha) - log(n)
  lq
}

gp05_expected_K <- function(n, alpha, gamma) {
  stopifnot(n >= 1L, alpha > 0, alpha < 1, gamma > 0, gamma < 1)
  if (n == 1L) return(list(EK = 1, pmf = 1, heuristic = 1))

  # Precompute all log_q values: log_q_mat[i, j] = log q_{alpha,0}(i:j)
  log_q_mat <- matrix(-Inf, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      log_q_mat[i, j] <- .log_q_gp05(i, j, alpha)
    }
  }

  # Recursion for log C_{n,m}(alpha)
  # C_{j,1} = q_{alpha,0}(j:j)  for j >= 1
  # C_{n,m} = sum_{j=1}^{n-m+1} q_{alpha,0}(n:j) * C_{n-j, m-1}
  log_C <- matrix(-Inf, n, n)

  # Base case: m=1
  for (nn in seq_len(n)) {
    log_C[nn, 1L] <- log_q_mat[nn, nn]
  }

  # Fill m = 2, ..., n
  for (mm in 2L:n) {
    for (nn in mm:n) {
      j_max <- nn - mm + 1L
      if (j_max < 1L) next
      log_terms <- numeric(j_max)
      for (j in seq_len(j_max)) {
        remainder <- nn - j
        if (remainder < mm - 1L) { log_terms[j] <- -Inf; next }
        log_terms[j] <- log_q_mat[nn, j] + log_C[remainder, mm - 1L]
      }
      mx <- max(log_terms)
      if (is.finite(mx)) {
        log_C[nn, mm] <- mx + log(sum(exp(log_terms - mx)))
      }
    }
  }

  # P(K=m) proportional to gamma^{m-1} * C_{n,m}
  log_unnorm <- numeric(n)
  for (mm in seq_len(n)) {
    log_unnorm[mm] <- (mm - 1L) * log(gamma) + log_C[n, mm]
  }

  mx <- max(log_unnorm[is.finite(log_unnorm)])
  probs <- exp(log_unnorm - mx)
  probs <- probs / sum(probs)

  EK <- sum(seq_len(n) * probs)
  list(
    EK    = EK,
    pmf   = probs,
    heuristic = 1 + alpha * gamma * sum(1 / seq_len(n - 1L))
  )
}


modular_osbm_sampler <- function(
    A, K, truth = NA,
    free   = c("psi","kappa","eta","z"),
    n_iter = 4000, burn = 500, thin = 1,
    verbose = FALSE,
    psi_constraint = c("WST","SST"),
    seed = NULL,
    eta_identifiability = c("none", "block_sum"),
    # NEW: controls
    shrink_when = c("after_z_sweep","end_of_iter","never"),
    refresh_pg_after_birth = TRUE,
    a_kappa = 1, b_kappa = 1,
    a_eta = 1, b_eta = 1,
    mu0 = 1, sigma0 = 2, tau0 = 0.15,
    alpha0 = 1.0, discount = 0.0,
    partition_prior = 'OCRP',
    gamma_gn = 0.8,
    alpha_gp05 = 0.5,
    theta_ocrp = 1.0,
    sst_birth_score_mode = c("exact_nonlocal", "local_approx"),
    slot_radius_after_burnin = NULL,
    sample_b_kappa = FALSE,
    alpha0_bkappa  = 1.0,    # hyperprior shape: b_kappa ~ Gamma(alpha0, beta0)
    beta0_bkappa   = 0.01,   # hyperprior rate:  prior mean = alpha0 / beta0 = 100
    # Mixing moves: set FALSE to disable swap/split-merge MH moves
    use_mixing_moves = TRUE,
    # DEBUG: dump sampler state to disk on NaN / invariant failure
    debug_dump_dir = NULL
){
  ## ---- hygiene -----------------------------------------------------------
  op <- options(
    warnPartialMatchArgs   = TRUE,
    warnPartialMatchDollar = TRUE,
    warnPartialMatchAttr   = TRUE
  )
  on.exit(options(op), add = TRUE)
  
  if (!is.null(seed)) set.seed(as.integer(seed))
  psi_constraint <- match.arg(psi_constraint)
  sst_birth_score_mode <- match.arg(sst_birth_score_mode)
  psi_mode <- if (psi_constraint == "WST") "pair" else "distance"
  eta_identifiability <- match.arg(eta_identifiability)
  shrink_when <- match.arg(shrink_when)
  partition_prior <- .normalize_partition_prior(partition_prior)
  
  if (!all(free %in% c("psi","kappa","eta","z")))
    stop("`free` must be subset of {'psi','kappa','eta','z'}.", call. = FALSE)
  
  needs <- setdiff(c("psi","kappa","eta","z"), free)
  if (length(needs) > 0) {
    if (exists(".assert_truth", mode = "function", inherits = TRUE)) {
      .assert_truth(truth, needs)
    } else {
      if (is.null(truth) || is.na(truth)) {
        stop("`truth` must be a named list when parameters are fixed.", call. = FALSE)
      }
      miss <- setdiff(needs, names(truth))
      if (length(miss)) {
        stop("`truth` is missing: ", paste(miss, collapse = ", "), call. = FALSE)
      }
    }
  }
  
  .assert_scalar_int(K, "K")
  .assert_scalar_int(n_iter, "n_iter")
  .assert_scalar_int(burn, "burn")
  .assert_scalar_int(thin, "thin")
  if (burn >= n_iter) stop("`burn` must be < `n_iter`.", call. = FALSE)
  if (thin < 1L)      stop("`thin` must be >= 1.", call. = FALSE)
  
  n <- nrow(A)
  .assert_scalar_int(n, "n (nrow(A))")
  .assert_matrix(A, n, "A")
  
  if (!requireNamespace("BayesLogit", quietly = TRUE))
    stop("Package 'BayesLogit' is required.", call. = FALSE)
  if (psi_constraint == "SST" &&
      !exists("update_psi_sst", mode = "function", inherits = TRUE))
    stop("`psi_constraint='SST'` requested but `update_psi_sst()` is not visible.", call. = FALSE)
  
  ## ---- hyper-parameters --------------------------------------------------
  hyper <- list(
    a_kappa = a_kappa, b_kappa = b_kappa,
    a_eta   = a_eta,   b_eta   = b_eta,
    mu0     = mu0,     sigma0  = sigma0,
    tau0    = tau0,
    partition_prior = partition_prior,
    gamma_gn = gamma_gn,
    alpha_gp05 = alpha_gp05,
    theta_ocrp = theta_ocrp
  )
  # 
  # hyper <- calibrate_osbm_hypers(A,
  #                              rho0 = 0.60,
  #                              r0 = 0.65,
  #                              a_eta = 2.0,
  #                              v_shape = 20,
  #                              gn_gamma = 0.95)
  # after hyper is created, once per run
  # hyper$alpha_b <- 10
  # hyper$beta_b  <- hyper$alpha_b / hyper$b_kappa  # prior mean = initial b_kappa
  # 
  psi_hyper <- list(
    mu0 = hyper$mu0, sigma0 = hyper$sigma0, tau0 = hyper$tau0,
    alpha0 = hyper$alpha0, discount = hyper$discount,
    partition_prior = hyper$partition_prior,
    gamma_gn = hyper$gamma_gn,
    alpha_gp05 = hyper$alpha_gp05,
    theta_ocrp = hyper$theta_ocrp
  )
  # right after you build hyper / psi_hyper
  
  
  
  prior_type <- partition_prior
  
  ## ---- flags -------------------------------------------------------------
  psi_est   <- "psi"   %in% free
  z_est     <- "z"     %in% free
  eta_est   <- "eta"   %in% free
  kappa_est <- "kappa" %in% free
  
  ## ---- edge list ---------------------------------------------------------
  A <- as.matrix(A)
  
  # total matches per dyad
  N_mat <- A + t(A)
  diag(N_mat) <- 0
  
  idx   <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  i_idx <- idx[,1]
  j_idx <- idx[,2]
  N_edge <- as.numeric(N_mat[idx])   # total matches
  
  # keep A as wins matrix for directional use
  
  edge_by_node <- replicate(n, integer(0), simplify = FALSE)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], as.integer(e))
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], as.integer(e))
  }
  
  
  ## ---- initial values ----------------------------------------------------
  # Save the requested initial K before any override
  K_init <- as.integer(K)
  # z: use the passed K_init to create a balanced starting partition
  z <- if (z_est) {
    as.integer(sample(seq_len(K_init), n, replace = TRUE))
  } else {
    .assert_len(truth$z, n, "truth$z"); as.integer(truth$z)
  }
  # K is now the number of *occupied* blocks (may be < K_init by chance)
  K <- max(tabulate(z, nbins = K_init) > 0L) * K_init  # at most K_init
  K <- length(unique(z))  # exact occupied count
  # Relabel z to 1..K (contiguous)
  uv <- sort(unique(z))
  if (!identical(uv, seq_len(K))) {
    remap <- integer(max(uv)); remap[uv] <- seq_len(K)
    z <- remap[z]
  }
  # ψ
  if (psi_est) {
    if (psi_mode == "pair") {
      psi <- matrix(0, K, K); psi[upper.tri(psi)] <- 0.25
    } else {
      psi <- cumsum(rep(0.25, max(K-1L,0)))
    }
  } else {
    if (psi_mode == "pair") {
      .assert_matrix(truth$psi, K, "truth$psi (WST pair)"); psi <- truth$psi
    } else {
      .assert_len(truth$psi, K-1, "truth$psi (SST distance)"); psi <- as.numeric(truth$psi)
    }
  }
  
  # η
  if (eta_est) {
    eta <- rgamma(n, hyper$a_eta, hyper$b_eta)
    if (identical(eta_identifiability, "block_sum")) {
      eta <- eta_rescale_by_block(eta, z, K)
    }
  } else {
    .assert_len(truth$eta, n, "truth$eta"); eta <- truth$eta
  }
  
  # κ
  if (kappa_est) {
    p <- pmin(z[i_idx], z[j_idx]); q <- pmax(z[i_idx], z[j_idx])
    p <- factor(p, levels = 1:K);  q <- factor(q, levels = 1:K)
    Rkl <- as.matrix(xtabs(N_edge ~ p + q, drop.unused.levels = FALSE))
    
    E_k   <- tapply(eta,   factor(z, levels = 1:K), sum);   E_k[is.na(E_k)] <- 0
    eta2k <- tapply(eta^2, factor(z, levels = 1:K), sum); eta2k[is.na(eta2k)] <- 0
    
    Tkl <- outer(E_k, E_k, `*`); diag(Tkl) <- pmax((E_k^2 - eta2k) / 2, 0)
    kappa <- (Rkl + hyper$a_kappa) / (Tkl + hyper$b_kappa)
    kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]
  } else {
    .assert_matrix(truth$kappa, K, "truth$kappa"); kappa <- truth$kappa
  }
  
  ## ---- storage (ragged -> lists) -----------------------------------------
  keep_seq <- seq(burn + 1L, n_iter, by = thin)
  n_keep   <- length(keep_seq)
  
  draws_z     <- vector("list", n_keep)
  draws_kappa <- vector("list", n_keep)
  draws_eta   <- vector("list", n_keep)
  draws_psi   <- vector("list", n_keep)
  K_trace     <- integer(n_keep)
  b_kappa_trace <- numeric(n_keep)   # b_kappa value at each saved iteration
  
  keep <- 0L
  t0 <- proc.time()[3]
  
  A_ij <- A[cbind(i_idx, j_idx)]
  
  ## ============================================================
  ## Main Gibbs loop (refactored SST variable-K using LOO updater)
  ## ============================================================
  
  # Optional: only keep omega_e if you truly need it elsewhere (e.g. diagnostics)
  need_omega_e <- FALSE  # set TRUE only if some later code uses omega_e
  
  .fast_sst_dimcheck <- function(z, kappa, psi) {
    K <- nrow(kappa)
    if (ncol(kappa) != K) stop("kappa not square.", call. = FALSE)
    if (length(psi) != max(K - 1L, 0L)) stop("psi length != K-1.", call. = FALSE)
    if (max(z) != K) stop("max(z) != K (labels not contiguous).", call. = FALSE)
    K
  }
  
for (it in seq_len(n_iter)) {

  ## --- (0) cheap coherence repair ONLY if needed --------------------------
  K <- nrow(kappa)

  if (psi_mode == "distance") {  # ===== SST =====
    # Goal: maintain the invariants required by the SST parameterisation:
    #   - labels contiguous 1..K
    #   - kappa is KxK
    #   - psi is length K-1 (distance vector)
    #
    # Why "only if needed"?
    # pack_state_sst() is not free (it relabels, shrinks, etc.).
    # Doing it every iter would be wasted work and can hide bugs.
    if (ncol(kappa) != K || length(psi) != max(K - 1L, 0L) || max(z) != K) {
      packed <- pack_state_sst(z, kappa, psi)
      z <- packed$z; kappa <- packed$kappa; psi <- packed$psi
      K <- packed$K
    }

    # One cheap check per iteration, not per node i.
    # If this fails, you want it to fail early rather than 200 lines later.
    K <- .fast_sst_dimcheck(z, kappa, psi)

  } else {
    # ===== WST =====
    # In WST, psi is a KxK matrix (pairwise).
    # There's no "psi length must be K-1" constraint.
    # You still want contiguous labels, but if you're using the state-returning
    # updater + pack_state_wst() after sweeps, you can often skip packing here.
    K <- nrow(kappa)
  }

  # omega_e is an OPTIONAL cache of per-edge Polya–Gamma variables.
  #
  # omega_dirty tracks whether the cached omega_e no longer matches (z, psi).
  # We reset to FALSE at the start of each iteration because:
  #   - omega_e, if it exists, is assumed correct at this moment
  #   - we will set omega_dirty = TRUE if we change z or psi later in the iter
  #
  # Important: omega_dirty is meaningless if we do not keep omega_e at all.
  omega_dirty <- FALSE


  ## ---- (A) psi block ------------------------------------------------------
  if (psi_est) {

    if (psi_mode == "pair") {  # ===== WST ψ update =====

      K <- nrow(kappa)

      # Packing here ensures the WST invariants before updating psi:
      #   - psi must be KxK
      #   - max(z) must equal K
      #
      # Why do it here rather than only at the top of the iter?
      # Because the previous iteration might have changed K via births/deaths,
      # and you want psi update to see a coherent (z, psi, kappa) state.
      if (nrow(psi) != K || ncol(psi) != K || max(z) != K) {
        packed <- pack_state_wst(z, kappa, psi)
        z <- packed$z; kappa <- packed$kappa; psi <- packed$psi; K <- packed$K
      }

      # This function changes psi (and only psi) by sampling ψ_{kℓ} | rest.
      psi <- update_psi_wst_pair(
        K = K, z = z,
        i_idx = i_idx, j_idx = j_idx,
        A_ij = A_ij, N_edge = N_edge,
        psi_curr = psi,
        mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2
      )

      # Why set omega_dirty = TRUE after psi changes?
      # Because ω_e ~ PG(N_e, θ_e) and θ_e depends on ψ (and z).
      # If you store ω_e and intend to reuse it later in this iteration,
      # it is now out-of-date.
      #
      # Why guard with if (need_omega_e)?
      # Because if you do NOT store ω_e, there is nothing to refresh.
      if (need_omega_e) omega_dirty <- TRUE

    } else {  # ===== SST ψ update =====

      # In SST you update the distance vector ψ_d using aggregated-by-distance
      # statistics (bar_y, bar_omega). This does NOT require ω_e per-edge;
      # it uses bar_omega ~ PG(B_d, ψ_d) (distance-pooled omega).
      #
      # Even though this update doesn't use ω_e, changing ψ still invalidates
      # any cached ω_e, if you happen to keep ω_e for some other reason.
      agg   <- aggregate_by_distance(
        K,
        z_i = z[i_idx], z_j = z[j_idx],
        A_ij = A_ij, N_edge = N_edge
      )
      bar_y <- agg$bar_y

      B_d <- distance_totals(K, z_i = z[i_idx], z_j = z[j_idx], N_edge = N_edge)
      bar_omega <- draw_omega_bar(B_d = B_d, psi = psi)

      psi <- update_psi_sst(
        K = K, bar_y = bar_y, bar_omega = bar_omega,
        psi_curr = psi,
        mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2,
        tau2_0 = hyper$tau0^2, n_inner_sweeps = 4L
      )

      # Same reasoning: psi changed => cached omega_e is now stale *if it exists*.
      if (need_omega_e) omega_dirty <- TRUE
    }
  }


  ## ---- (B) z block (variable-K) ------------------------------------------
  if (z_est) {

    if (psi_mode == "pair") {  # ===== WST z-update (variable K) =====

      # Ensure coherent state before doing many single-site updates.
      K <- nrow(kappa)
      if (nrow(psi) != K || ncol(psi) != K || max(z) != K) {
        packed <- pack_state_wst(z, kappa, psi)
        z <- packed$z; kappa <- packed$kappa; psi <- packed$psi; K <- packed$K
      }

      # Draw edge-level PG latents ONCE for the entire z-sweep.
      # This is the PG-augmented approach: omega ~ PG(N_edge, |phi|).
      omega_edge <- draw_edge_omega(z, psi, i_idx, j_idx, N_edge, mode = "pair")

      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)

      slot_rad <- K + 1L

      node_order <- sample.int(n)
      for (i in node_order) {
        z_old_i <- z[i]
        K_old   <- K

        res <- wst_update_i_with_birth_LOO_pg(
          i = i,
          A = A, z = z, eta = eta, kappa = kappa, psi = psi,
          Rkl = bt$Rkl, Tkl = bt$Tkl,
          i_idx = i_idx, j_idx = j_idx, N_edge = N_edge, edge_by_node = edge_by_node,
          a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
          gamma_gn = psi_hyper$gamma_gn,
          mu0 = hyper$mu0, sigma0 = hyper$sigma0, sig2_0 = hyper$sigma0^2,
          omega_edge = omega_edge,
          slot_radius = slot_rad,
          partition_prior = psi_hyper$partition_prior,
          alpha_gp05 = psi_hyper$alpha_gp05,
          theta_ocrp = psi_hyper$theta_ocrp
        )

        z     <- res$z
        kappa <- res$kappa
        psi   <- res$psi
        K     <- res$K

        if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
          bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
          if (need_omega_e) omega_dirty <- TRUE
        }
      }

      packed <- pack_state_wst(z, kappa, psi)
      z <- packed$z; kappa <- packed$kappa; psi <- packed$psi; K <- packed$K

    } else if (psi_mode == "distance") {  # ===== SST z-update (variable K) =====

      # Draw edge-level PG latents ONCE for the entire z-sweep.
      omega_edge <- draw_edge_omega(z, psi, i_idx, j_idx, N_edge, mode = "distance")

      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
      slot_rad <- K + 1

      node_order <- sample.int(n)
      for (i in node_order) {
        z_old_i <- z[i]
        K_old   <- K

        res <- sst_update_i_with_birth_LOO_pg(
          i = i,
          A = A, z = z, eta = eta, kappa = kappa, psi = psi,
          Rkl = bt$Rkl, Tkl = bt$Tkl,
          i_idx = i_idx, j_idx = j_idx, N_edge = N_edge, edge_by_node = edge_by_node,
          a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
          gamma_gn = psi_hyper$gamma_gn,
          omega_edge = omega_edge,
          slot_radius = slot_rad,
          tau0 = psi_hyper$tau0,
          mu0 = psi_hyper$mu0,
          sig2_0 = psi_hyper$sigma0^2,
          partition_prior = psi_hyper$partition_prior,
          alpha_gp05 = psi_hyper$alpha_gp05,
          theta_ocrp = psi_hyper$theta_ocrp,
          birth_score_mode = sst_birth_score_mode
        )

        z     <- res$z
        kappa <- res$kappa
        psi   <- res$psi
        K     <- res$K

        if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
          bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
          if (need_omega_e) omega_dirty <- TRUE
        }
      }

      # SST post-sweep invariant check
      K <- .fast_sst_dimcheck(z, kappa, psi)
    }
  }


  ## ---- (B2) Adjacent-block swap MH move -----------------------------------
  # Proposes swapping two adjacent blocks k <-> k+1.
  # Symmetric proposal; only directional likelihood changes.
  if (use_mixing_moves && z_est && K >= 2L) {
    swap_res <- adjacent_block_swap_move(
      z = z, kappa = kappa, psi = psi, eta = eta,
      A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      psi_mode = psi_mode
    )
    z <- swap_res$z; kappa <- swap_res$kappa; psi <- swap_res$psi
    if (isTRUE(swap_res$accepted) && need_omega_e) omega_dirty <- TRUE
  }


  ## ---- (B3) Split-merge MH move ------------------------------------------
  # Proposes splitting a block in two or merging two adjacent blocks.
  if (use_mixing_moves && z_est && K >= 1L) {
    sm_res <- split_merge_move(
      z = z, kappa = kappa, psi = psi, eta = eta,
      A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
      gamma_gn = psi_hyper$gamma_gn,
      psi_mode = psi_mode,
      hyper_psi = list(mu0 = hyper$mu0, sigma0 = hyper$sigma0, tau0 = hyper$tau0),
      n_restricted_scans = 3L,
      partition_prior = psi_hyper$partition_prior,
      theta_ocrp = psi_hyper$theta_ocrp
    )
    z <- sm_res$z; kappa <- sm_res$kappa; psi <- sm_res$psi
    if (isTRUE(sm_res$accepted)) {
      K <- if (psi_mode == "pair") nrow(kappa) else length(psi) + 1L
      if (need_omega_e) omega_dirty <- TRUE
    }
  }

  if (eta_est && identical(eta_identifiability, "block_sum")) {
    eta <- eta_rescale_by_block(eta, z, K)
  }


  ## ---- (C) omega_e (ONLY if you keep it) ---------------------------------
  # This is the *only* place where omega_e is actually refreshed.
  #
  # Why here?
  # Because you want to avoid expensive draws unless you have to,
  # and you want ω_e to correspond to the *final* (z, psi) state
  # after both psi-update and z-sweep in the current iteration.
  #
  # Why the AND condition?
  # - need_omega_e: only refresh if somebody will actually use omega_e later
  # - omega_dirty: only refresh if z/psi changed this iter
  if (need_omega_e && omega_dirty) {
    omega_e <- draw_pg_latents(i_idx, j_idx, N_edge, z, psi, mode = psi_mode)
    omega_dirty <- FALSE
  }


  ## ---- (D) kappa block ----------------------------------------------------
  if (kappa_est) {
    # κ depends on (z, eta) through Rkl, Tkl.
    # After the z-sweep (and possibly eta update later), you rebuild totals
    # and sample κ fresh.
    bt  <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
    Rkl <- bt$Rkl
    Tkl <- bt$Tkl

    kappa <- matrix(0, K, K)
    # Vectorized kappa sampling from Gamma posteriors
    shape_mat <- hyper$a_kappa + Rkl
    rate_mat  <- hyper$b_kappa + Tkl
    # Fill upper triangle (including diagonal)
    ut_idx <- which(upper.tri(kappa, diag = TRUE), arr.ind = TRUE)
    kappa[ut_idx] <- rgamma(nrow(ut_idx),
                            shape = shape_mat[ut_idx],
                            rate  = rate_mat[ut_idx])
    # Symmetrize
    kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]
  }

  ## ---- (D2) b_kappa update (hierarchical, optional) ----------------------
  # Partially collapsed Gibbs: collapse kappa for z-updates (step B), then
  # un-collapse here to draw b_kappa | kappa ~ Gamma(alpha0 + P*a_kappa,
  #                                                    beta0  + sum(kappa_{ut}))
  # where P = K*(K+1)/2 (number of unique block pairs).
  # Self-regulation: K↑ => more kappa params => sum(kappa)↑ => b_kappa↑ => births suppressed.
  if (sample_b_kappa && kappa_est) {
    hyper$b_kappa <- update_b_kappa(kappa,
                                    a_kappa = hyper$a_kappa,
                                    alpha_b = alpha0_bkappa,
                                    beta_b  = beta0_bkappa)
  }


  ## ---- (E) eta block ------------------------------------------------------
  if (eta_est) {
    # η depends on κ and z; it’s sampled after κ in this design.
    eta <- update_eta_all_dyads(
      eta = eta, z = z, kappa = kappa,
      a_eta = hyper$a_eta, b_eta = hyper$b_eta,
      i_idx = i_idx, j_idx = j_idx, N_edge = N_edge
    )
    if (identical(eta_identifiability, "block_sum")) {
      eta <- eta_rescale_by_block(eta, z, K)
    }
  }


  ## ---- save draws ---------------------------------------------------------
  if (it %in% keep_seq) {
    keep <- keep + 1L
    draws_z[[keep]]     <- as.integer(z)
    draws_kappa[[keep]] <- kappa
    draws_eta[[keep]]   <- as.numeric(eta)
    draws_psi[[keep]]   <- psi
    K_trace[keep]       <- nrow(kappa)
    b_kappa_trace[keep] <- hyper$b_kappa
  }

  # DEBUG: check for NaN/Inf in psi after each iteration
  if (!is.null(debug_dump_dir) || verbose) {
    bad_psi <- if (psi_mode == "distance") {
      anyNA(psi) || any(!is.finite(psi))
    } else {
      anyNA(psi) || any(!is.finite(psi[upper.tri(psi)]))
    }
    if (bad_psi) {
      msg <- sprintf("[DEBUG|it=%d] psi contains NA/Inf/NaN (K=%d, mode=%s)",
                     it, nrow(kappa), psi_mode)
      if (!is.null(debug_dump_dir)) {
        dir.create(debug_dump_dir, showWarnings = FALSE, recursive = TRUE)
        fn <- file.path(debug_dump_dir,
                        sprintf("bad_psi_it%05d.rds", it))
        saveRDS(list(it = it, K = nrow(kappa), z = z, psi = psi,
                     kappa = kappa, eta = eta, mode = psi_mode), fn)
        message(msg, "  -> state dumped to ", fn)
      } else {
        warning(msg)
      }
      # Hard repair: reset bad psi entries to small positive
      if (psi_mode == "distance") {
        psi[!is.finite(psi)] <- 0.01
        psi <- pmax(psi, 0)
      } else {
        psi[!is.finite(psi)] <- 0
        psi <- pmax(psi, 0)
      }
    }
  }

  if (verbose && it %% 20 == 0) {
    cat("iter", it, "/", n_iter, " saved:", keep, "  K=", nrow(kappa), "  ",
        round(proc.time()[3] - t0, 2), "s \r")
  }
}

  
  
  
  invisible(list(
    z        = draws_z,
    psi      = draws_psi,
    kappa    = draws_kappa,
    eta      = draws_eta,
    K_trace  = K_trace,
    b_kappa_trace = b_kappa_trace,
    keep     = keep_seq,
    meta     = list(psi_mode = psi_mode, psi_constraint = psi_constraint,
                    hyper = hyper, free = free,
                    eta_identifiability = eta_identifiability,
                    shrink_when = shrink_when,
                    refresh_pg_after_birth = refresh_pg_after_birth,
                    sample_b_kappa = sample_b_kappa,
                    alpha0_bkappa = alpha0_bkappa,
                    beta0_bkappa  = beta0_bkappa)
  ))
}




# ===== Minimal assertion helpers (no extra packages) =========================
.assert_scalar_int <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x != as.integer(x) || x < 1L)
    stop(sprintf("`%s` must be a positive integer scalar.", name), call. = FALSE)
  invisible(TRUE)
}
.assert_logical_scalar <- function(x, name){
  if (length(x) != 1L || is.na(x) || !is.logical(x))
    stop(sprintf("`%s` must be TRUE/FALSE.", name), call. = FALSE)
  invisible(TRUE)
}
.assert_len <- function(x, len, name){
  if (length(x) != len) stop(sprintf("`%s` must have length %d.", name, len), call. = FALSE)
  invisible(TRUE)
}
.assert_matrix <- function(A, n, name){
  if (!is.matrix(A) && !inherits(A,"Matrix"))
    stop(sprintf("`%s` must be a matrix or Matrix.", name), call. = FALSE)
  if (nrow(A) != n || ncol(A) != n)
    stop(sprintf("`%s` must be square %dx%d.", name, n, n), call. = FALSE)
  if (anyNA(A)) stop(sprintf("`%s` contains NA.", name), call. = FALSE)
  invisible(TRUE)
}
.assert_truth <- function(truth, needs){
  if (is.null(truth) || is.na(truth)) stop("`truth` must be a named list when parameters are fixed.", call. = FALSE)
  miss <- setdiff(needs, names(truth))
  if (length(miss)) stop("`truth` is missing: ", paste(miss, collapse=", "), call. = FALSE)
  invisible(TRUE)
}
.safe_log <- function(x) log(pmax(x, 1e-15))


update_psi_sst <- function(K, bar_y, bar_omega, psi_curr,
                           mu0, sig2_0, tau2_0,
                           n_inner_sweeps = 2)
{
  D <- K - 1
  if (D == 0) return(numeric(0))
  
  # Build T, Ω, Q, g, and prior precision/shift
  Tmat  <- make_T(K)
  Omega <- diag(bar_omega, nrow = D, ncol = D)
  Q     <- crossprod(Tmat, Omega %*% Tmat)   # T^T Ω T
  g     <- crossprod(Tmat, bar_y)            # T^T \bar y   (note: no Ω)
  
  Vinv  <- diag(c(1/sig2_0, rep(1/tau2_0, max(0, D-1))), D, D)  # diag precision
  muvec <- c(mu0, rep(0, max(0, D-1)))
  
  Qstar <- Q + Vinv
  gstar <- g + Vinv %*% muvec
  
  # Start at current δ (invert the cumulative sum): δ0 = ψ1, δr = ψ_{r+1} - ψ_r
  delta <- numeric(D)
  if (length(psi_curr) == D) {
    delta[1] <- max(psi_curr[1], 1e-10)
    if (D >= 2) {
      incs <- pmax(psi_curr[-1] - psi_curr[-D], 1e-10)
      delta[2:D] <- incs
    }
  } else {
    # fall back if no current ψ provided
    delta[] <- 0.1
  }
  
  # Coordinate-wise truncated-Normal updates on δ_r > 0
  for (s in seq_len(n_inner_sweeps)) {
    for (r in seq_len(D)) {
      a_rr <- Qstar[r, r]
      b_r  <- as.numeric(gstar[r] - sum(Qstar[r, -r, drop = FALSE] * delta[-r]))
      mean_r <- b_r / a_rr
      sd_r   <- 1 / sqrt(a_rr)
      delta[r] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = mean_r, sd = sd_r)
      if (!is.finite(delta[r]) || delta[r] <= 0) delta[r] <- 1e-10
    }
  }
  
  # Map back to ψ = T δ (monotone & positive by construction)
  as.numeric(Tmat %*% delta)
}

make_log1p_tables <- function(psi){ 
  list(pos = log1pexp( psi),
       neg = log1pexp(-psi))}


# Node-wise Gibbs update for z (distance ψ)
z_update_osbm_distance <- function(i, A, z, eta, kappa, psi, alpha_vec,
                                   log1p_tab, log_kappa,
                                   i_idx = NULL, j_idx = NULL, N_edge = NULL,
                                   edge_by_node = NULL)
{
  K <- length(alpha_vec)
  idxK <- seq_len(K)
  
  # dyad summaries wrt current z (same as your original code)
  e_out <- A[i, ];  e_in <- A[, i];  Nij <- e_out + e_in
  N_blk <- A_out_blk <- numeric(K)
  if (any(Nij > 0)) {
    ell <- z[which(Nij > 0)]
    N_blk_tmp     <- tapply(Nij[Nij > 0], ell, sum)
    A_out_blk_tmp <- tapply(e_out[Nij > 0], ell, sum)
    N_blk[as.integer(names(N_blk_tmp))]         <- N_blk_tmp
    A_out_blk[as.integer(names(A_out_blk_tmp))] <- A_out_blk_tmp
  }
  
  n_minus     <- tabulate(z[-i], nbins = K)
  eta_sum_blk <- tapply(eta[-i], factor(z[-i], levels = idxK), sum)
  eta_sum_blk[is.na(eta_sum_blk)] <- 0
  
  lp_vec <- sapply(
    idxK, log_post_fast_distance, i           = i,
    N_blk       = N_blk,
    A_out_blk   = A_out_blk,
    n_minus     = n_minus,
    eta_sum_blk = eta_sum_blk,
    eta_i       = eta[i],
    kappa       = kappa,
    psi         = psi,
    alpha_vec   = alpha_vec,
    log_kappa   = log_kappa,
    log1p_tab   = log1p_tab
  )
  
  p <- exp(lp_vec - max(lp_vec))
  sample.int(K, 1L, prob = p)
}


log_post_fast_distance <- function(i, k,
                                   N_blk, A_out_blk,
                                   n_minus,
                                   eta_sum_blk, eta_i,
                                   kappa, psi, alpha_vec,
                                   log_kappa, log1p_tab)
{
  log2 <- log(2)
  lp   <- log(alpha_vec[k] + n_minus[k])
  lambda_vec_k <- kappa %*% eta_sum_blk
  lp <- lp - eta_i * lambda_vec_k[k] + (log_kappa %*% N_blk)[k]
  
  idx        <- seq_along(alpha_vec)
  mask_cross <- idx != k
  d_vec      <- abs(idx - k)[mask_cross]
  
  sign_vec <- ifelse(idx[mask_cross] > k,  1, -1)
  
  A_dir      <- A_out_blk[mask_cross]
  lp <- lp +
    sum(sign_vec * A_dir * psi[d_vec]) -
    sum(N_blk[mask_cross] *
          ifelse(sign_vec == 1,
                 log1p_tab$pos[d_vec],
                 log1p_tab$neg[d_vec])) -
    N_blk[k] * log2
  lp
}
