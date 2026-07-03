# Ordered SBM sampler helpers: state packing, partition priors, and
# collapsed update utilities used by the variable-K transitive sampler.
# The functions in this file were extracted from core/transitive_sbm_sampler.R
# so the main sampler body can stay focused on the Gibbs loop itself.

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

  if (!prior %in% c("GN", "OCRP", "ROCRP")) {
    stop(
      "`partition_prior` must be one of 'GN', 'OCRP', or 'ROCRP' in the cleaned bundle.",
      call. = FALSE
    )
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

