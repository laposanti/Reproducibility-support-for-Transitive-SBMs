# SST-specific ordered SBM helpers for birth moves, directional scores, and psi updates.
#
#Returns a if not NULL, else b.
`%||%` <- function(a, b) if (!is.null(a)) a else b

#\(cosh(x)=\frac{e^{x}+e^{-x}}{2}\)#
# ---- stable log(2 cosh(x/2)) (used in SST) ---------------------------------
.log2cosh_half <- function(x) {
  u <- abs(x) / 2
  u + log1p(exp(-2*u))
}
make_T <- function(K){
  D <- K - 1
  if (D <= 0) return(matrix(0,0,0))
  T <- matrix(0, nrow = D, ncol = D)
  for (r in 1:D) T[r, 1:r] <- 1
  T
}
logsig     <- function(x) -log1pexp(-x)
log1m_sig  <- function(x) -log1pexp(x)
log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
# replace your lse with this:
lse <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(-Inf)
  m <- max(x[is.finite(x)])
  m + log(sum(exp(x[is.finite(x)] - m)))
}

# -------------------------------------------------------------------
# Directional: exact new-slot log-likelihood via aggregated counts
# ψ_vec is length K-1 with ψ_d for distance d
# Returns vector over r = 1..K+1
# -------------------------------------------------------------------

gamma_from_mean_var <- function(mean, var) {
  if (any(mean <= 0)) {
    stop("Mean must be strictly positive.")
  }
  if (any(var <= 0)) {
    stop("Variance must be strictly positive.")
  }
  
  shape <- mean^2 / var
  rate  <- mean / var
  
  list(shape = shape, rate = rate)
}

# =============================================================================
# PG-AUGMENTED SST DIRECTIONAL SCORES
# =============================================================================

# Existing-block score using quadratic PG kernel: Y*psi - 0.5*Omega*psi^2
dir_pg_SST_existing_counts <- function(k, c_plus, N_tot, Omega_blk, psi_vec) {
  K <- length(N_tot)
  if (K <= 1L) return(0)
  val <- 0
  for (ell in seq_len(K)) {
    if (ell == k) next
    n <- N_tot[ell]
    if (n == 0) next
    d  <- abs(ell - k)
    th <- psi_vec[d]
    A_fwd <- if (ell > k) c_plus[ell] else (n - c_plus[ell])
    Y <- A_fwd - 0.5 * n
    val <- val + Y * th - 0.5 * Omega_blk[ell] * th^2
  }
  val
}

# Closed-form integral for SST extreme distance (d = K_old) under PG kernel.
# Integrates delta ~ N+(m_inc, tau0^2) with PG kernel.
log_int_pg_sst_extreme <- function(Y_star, Omega_star, psiKm1, tau0,
                                   m_inc = 0) {
  tau2_0 <- tau0^2
  mu_u   <- psiKm1 + m_inc
  A_prec <- Omega_star + 1 / tau2_0
  B_loc  <- Y_star + mu_u / tau2_0
  m_post <- B_loc / A_prec
  s_post <- 1 / sqrt(A_prec)

  log_I <- -0.5 * log(tau2_0 * A_prec) +
    B_loc^2 / (2 * A_prec) - mu_u^2 / (2 * tau2_0) +
    pnorm((m_post - psiKm1) / s_post, log.p = TRUE) -
    pnorm(m_inc / tau0, log.p = TRUE)

  log_I
}

# New-block score for SST (local approximation): known distances use PG kernel,
# new farthest distance uses closed-form integral using only node-i terms.
dir_pg_SST_new_local_approx <- function(c_plus, N_tot, Omega_blk, psi_vec,
                                        r_set = NULL, tau0, m_inc = 0) {
  K_old <- length(N_tot)
  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)

  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(numeric(0))
  if (length(psi_vec) != (K_old - 1L))
    stop("SST-PG: psi_vec must have length K_old-1.", call. = FALSE)

  for (idx in seq_along(r_set)) {
    r <- r_set[idx]
    val <- 0
    A_K <- 0; n_K <- 0; Omega_K <- 0

    for (ell in seq_len(K_old)) {
      n <- N_tot[ell]
      if (n == 0) next
      pos_ell <- if (ell < r) ell else ell + 1L
      d <- abs(pos_ell - r)
      if (d == 0L) next
      A_fwd <- if (pos_ell > r) c_plus[ell] else (n - c_plus[ell])
      Y_d <- A_fwd - 0.5 * n

      if (d <= (K_old - 1L)) {
        th <- psi_vec[d]
        val <- val + Y_d * th - 0.5 * Omega_blk[ell] * th^2
      } else if (d == K_old) {
        A_K <- A_K + A_fwd
        n_K <- n_K + n
        Omega_K <- Omega_K + Omega_blk[ell]
      }
    }

    if (n_K > 0) {
      Y_K <- A_K - 0.5 * n_K
      psiKm1 <- if (K_old == 1L) 0 else psi_vec[K_old - 1L]
      val <- val + log_int_pg_sst_extreme(Y_K, Omega_K, psiKm1, tau0, m_inc)
    }

    lp[idx] <- val
  }
  lp
}

# Backward-compatible alias; prefer dir_pg_SST_new_exact_nonlocal in SST birth updates.
dir_pg_SST_new <- dir_pg_SST_new_local_approx

# Old-old directional sufficient statistics by packed block pair (a < b),
# excluding dyads that involve node i.
.sst_oldold_pair_stats_excluding_i <- function(i, z_packed, A,
                                               i_idx, j_idx, N_edge,
                                               omega_edge, K_minus) {
  S_old <- matrix(0, K_minus, K_minus)
  O_old <- matrix(0, K_minus, K_minus)
  E <- length(N_edge)
  for (e in seq_len(E)) {
    u <- i_idx[e]; v <- j_idx[e]
    if (u == i || v == i) next

    zu <- z_packed[u]; zv <- z_packed[v]
    if (zu == zv || zu <= 0L || zv <= 0L) next

    a <- min(zu, zv); b <- max(zu, zv)
    A_uv <- A[u, v]
    A_fwd <- if (zu < zv) A_uv else (N_edge[e] - A_uv)
    Y_e <- A_fwd - 0.5 * N_edge[e]

    S_old[a, b] <- S_old[a, b] + Y_e
    O_old[a, b] <- O_old[a, b] + omega_edge[e]
  }
  list(S_old = S_old, O_old = O_old)
}

# Non-local old-old birth terms over insertion ranks r = 1..K_minus+1.
# Returns:
# - corr_fixed[r]: old-old correction over fixed distances d <= K_minus-2.
# - YK_old[r], OK_old[r]: old-old contributions that become distance K_minus.
# - baseline_extreme[r]: old distance-(K_minus-1) baseline subtracted when those
#   pairs are folded into the new extreme-distance integral.
.sst_oldold_nonlocal_birth_terms <- function(S_old, O_old, psi_vec) {
  K_minus <- nrow(S_old)
  R <- K_minus + 1L
  corr_fixed <- rep(0, R)
  YK_old <- rep(0, R)
  OK_old <- rep(0, R)
  baseline_extreme <- rep(0, R)

  if (K_minus <= 1L) {
    return(list(corr_fixed = corr_fixed, YK_old = YK_old, OK_old = OK_old,
                baseline_extreme = baseline_extreme))
  }

  D <- rep(0, R + 1L)
  for (a in seq_len(K_minus - 1L)) {
    for (b in seq.int(a + 1L, K_minus)) {
      S_ab <- S_old[a, b]
      O_ab <- O_old[a, b]
      if (S_ab == 0 && O_ab == 0) next

      d <- b - a
      if (d <= (K_minus - 2L)) {
        th_d <- psi_vec[d]
        th_next <- psi_vec[d + 1L]
        delta_ab <- S_ab * (th_next - th_d) -
          0.5 * O_ab * (th_next^2 - th_d^2)

        lo <- a + 1L; hi <- b
        D[lo] <- D[lo] + delta_ab
        D[hi + 1L] <- D[hi + 1L] - delta_ab
      } else if (d == (K_minus - 1L)) {
        lo <- a + 1L; hi <- b
        YK_old[lo:hi] <- YK_old[lo:hi] + S_ab
        OK_old[lo:hi] <- OK_old[lo:hi] + O_ab

        th_old <- psi_vec[d]
        old_kernel <- S_ab * th_old - 0.5 * O_ab * th_old^2
        baseline_extreme[lo:hi] <- baseline_extreme[lo:hi] + old_kernel
      }
    }
  }

  corr_fixed <- cumsum(D[seq_len(R)])
  list(corr_fixed = corr_fixed, YK_old = YK_old, OK_old = OK_old,
       baseline_extreme = baseline_extreme)
}

# Exact non-local new-slot score for SST births.
# Adds old-old straddle correction and folds old-old distance-K terms into the
# new extreme-distance integral, relative to the leave-one-out old-old baseline.
dir_pg_SST_new_exact_nonlocal <- function(c_plus, N_tot, Omega_blk, psi_vec,
                                          S_old, O_old,
                                          r_set = NULL, tau0, m_inc = 0,
                                          first_mu0 = 0,
                                          first_sig2_0 = NULL) {
  K_old <- length(N_tot)
  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)

  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(numeric(0))
  if (length(psi_vec) != (K_old - 1L))
    stop("SST-PG: psi_vec must have length K_old-1.", call. = FALSE)
  if (!all(dim(S_old) == c(K_old, K_old)) || !all(dim(O_old) == c(K_old, K_old)))
    stop("SST-PG: old-old stats must be K_old x K_old matrices.", call. = FALSE)

  oldold_terms <- .sst_oldold_nonlocal_birth_terms(S_old, O_old, psi_vec)

  for (idx in seq_along(r_set)) {
    r <- r_set[idx]
    val <- oldold_terms$corr_fixed[r] - oldold_terms$baseline_extreme[r]
    YK_node <- 0
    OK_node <- 0

    for (ell in seq_len(K_old)) {
      n <- N_tot[ell]
      if (n == 0) next
      pos_ell <- if (ell < r) ell else ell + 1L
      d <- abs(pos_ell - r)
      if (d == 0L) next
      A_fwd <- if (pos_ell > r) c_plus[ell] else (n - c_plus[ell])
      Y_d <- A_fwd - 0.5 * n

      if (d <= (K_old - 1L)) {
        th <- psi_vec[d]
        val <- val + Y_d * th - 0.5 * Omega_blk[ell] * th^2
      } else if (d == K_old) {
        YK_node <- YK_node + Y_d
        OK_node <- OK_node + Omega_blk[ell]
      }
    }

    YK_tot <- YK_node + oldold_terms$YK_old[r]
    OK_tot <- OK_node + oldold_terms$OK_old[r]
    if (YK_tot != 0 || OK_tot != 0) {
      psiKm1 <- if (K_old == 1L) 0 else psi_vec[K_old - 1L]
      sd_ext <- if (K_old == 1L && !is.null(first_sig2_0)) sqrt(first_sig2_0) else tau0
      mean_ext <- if (K_old == 1L) first_mu0 else m_inc
      val <- val + log_int_pg_sst_extreme(YK_tot, OK_tot, psiKm1, sd_ext, mean_ext)
    }

    lp[idx] <- val
  }
  lp
}

# Brute-force exact non-local scorer over all dyads for candidate new slots.
# This is slower but directly evaluates the incremental target used in the
# single-site conditional: candidate birth log-kernel minus the LOO old-old
# baseline. Existing-block scores already omit that baseline.
dir_pg_SST_new_exact_nonlocal_bruteforce <- function(i, z_packed, A,
                                                      i_idx, j_idx, N_edge,
                                                      omega_edge, psi_vec,
                                                      r_set = NULL,
                                                      tau0, m_inc = 0,
                                                      first_mu0 = 0,
                                                      first_sig2_0 = NULL) {
  K_old <- max(z_packed)
  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)

  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(numeric(0))
  if (length(psi_vec) != (K_old - 1L)) {
    stop("SST-PG brute: psi_vec must have length K_old-1.", call. = FALSE)
  }

  E <- length(N_edge)
  baseline_oldold <- 0
  for (e in seq_len(E)) {
    u <- i_idx[e]
    v <- j_idx[e]
    if (u == i || v == i) next

    zu <- z_packed[u]
    zv <- z_packed[v]
    if (zu == zv) next

    d <- abs(zu - zv)
    th <- psi_vec[d]
    A_uv <- A[u, v]
    A_fwd <- if (zu < zv) A_uv else (N_edge[e] - A_uv)
    Y <- A_fwd - 0.5 * N_edge[e]
    baseline_oldold <- baseline_oldold + Y * th - 0.5 * omega_edge[e] * th^2
  }

  for (idx in seq_along(r_set)) {
    r <- r_set[idx]
    val <- 0
    YK <- 0
    OK <- 0

    for (e in seq_len(E)) {
      u <- i_idx[e]
      v <- j_idx[e]

      zu <- if (u == i) r else {
        z_u <- z_packed[u]
        if (z_u < r) z_u else z_u + 1L
      }
      zv <- if (v == i) r else {
        z_v <- z_packed[v]
        if (z_v < r) z_v else z_v + 1L
      }

      if (zu == zv) next

      a <- min(zu, zv)
      b <- max(zu, zv)
      d <- b - a

      A_uv <- A[u, v]
      A_fwd <- if (zu < zv) A_uv else (N_edge[e] - A_uv)
      Y <- A_fwd - 0.5 * N_edge[e]

      if (d <= (K_old - 1L)) {
        th <- psi_vec[d]
        val <- val + Y * th - 0.5 * omega_edge[e] * th^2
      } else if (d == K_old) {
        YK <- YK + Y
        OK <- OK + omega_edge[e]
      }
    }

    if (YK != 0 || OK != 0) {
      psiKm1 <- if (K_old == 1L) 0 else psi_vec[K_old - 1L]
      sd_ext <- if (K_old == 1L && !is.null(first_sig2_0)) sqrt(first_sig2_0) else tau0
      mean_ext <- if (K_old == 1L) first_mu0 else m_inc
      val <- val + log_int_pg_sst_extreme(YK, OK, psiKm1, sd_ext, mean_ext)
    }

    lp[idx] <- val - baseline_oldold
  }
  lp
}

# =============================================================================
# PG-AUGMENTED z-UPDATE WRAPPER — SST
# =============================================================================
sst_update_i_with_birth_LOO_pg <- function(
    i,
    A, z, eta, kappa, psi,
    Rkl, Tkl,
    i_idx, j_idx, N_edge, edge_by_node,
    a_kappa, b_kappa,
    gamma_gn,
    omega_edge,
    r_set = NULL,
    slot_radius = NULL,
    tau0,
    mu0 = 0,
    sig2_0 = 1,
    partition_prior = "OCRP",
    theta_ocrp = 1.0,
    birth_score_mode = c("exact_nonlocal", "local_approx")
) {
  birth_score_mode <- match.arg(birth_score_mode)
  K_full <- nrow(kappa)
  oldk  <- z[i]
  eta_i <- eta[i]

  # ---- LOO occupancy and packing ----
  n_minus_full <- tabulate(z[-i], nbins = K_full)
  keep_full <- sort(which(n_minus_full > 0L))
  K_minus <- length(keep_full)
  if (K_minus < 1L) stop("LOO state has K_minus<1.", call. = FALSE)

  map <- integer(K_full)
  map[keep_full] <- seq_len(K_minus)
  z_packed <- map[z]
  oldk_pos <- match(oldk, keep_full)

  # ---- pack parameters ----
  kappa_minus <- kappa[keep_full, keep_full, drop = FALSE]
  psi_minus   <- reindex_psi_sst_keep(psi_old = psi, K_old = K_full,
                                       keep = keep_full)
  if (length(psi_minus) != max(K_minus - 1L, 0L))
    stop("psi_minus length mismatch after packing.", call. = FALSE)

  # ---- node-vs-block counts ----
  C_i_full <- counts_by_block_exact_cpp(
    i = i, A = A, z = as.integer(z),
    i_idx = as.integer(i_idx), j_idx = as.integer(j_idx),
    N_edge = as.numeric(N_edge),
    edge_by_node = edge_by_node, K = K_full
  )
  c_plus <- as.numeric(C_i_full$c_plus)[keep_full]
  N_tot  <- as.numeric(C_i_full$N_tot)[keep_full]

  # ---- aggregate omega ----
  Omega_blk_full <- aggregate_omega_by_block(
    i, z, omega_edge, i_idx, j_idx, edge_by_node, K_full
  )
  Omega_blk <- Omega_blk_full[keep_full]

  # ---- exposure sums excluding i ----
  E_excl_full <- E_by_block_excluding_i(i, z, eta, K_full)
  E_excl <- as.numeric(E_excl_full[keep_full])
  r_add <- N_tot
  t_add <- eta_i * E_excl

  # ---- LOO base R/T ----
  R_minus <- Rkl[keep_full, keep_full, drop = FALSE]
  T_minus <- Tkl[keep_full, keep_full, drop = FALSE]
  if (!is.na(oldk_pos)) {
    for (ell_pos in seq_len(K_minus)) {
      subR <- N_tot[ell_pos]
      subT <- eta_i * E_excl[ell_pos]
      p <- min(oldk_pos, ell_pos); q <- max(oldk_pos, ell_pos)
      R_minus[p, q] <- R_minus[p, q] - subR
      T_minus[p, q] <- T_minus[p, q] - subT
      R_minus[q, p] <- R_minus[p, q]
      T_minus[q, p] <- T_minus[p, q]
      if (R_minus[p, q] < 0) R_minus[p, q] <- R_minus[q, p] <- 0
      if (T_minus[p, q] < 0) T_minus[p, q] <- T_minus[q, p] <- 0
    }
  }

  # ---- candidate slots ----
  if (is.null(r_set)) {
    if (!is.null(slot_radius) && !is.na(oldk_pos)) {
      lo <- max(1L, oldk_pos - (slot_radius - 1L))
      hi <- min(K_minus + 1L, oldk_pos + slot_radius)
      r_set <- seq.int(lo, hi)
    } else {
      r_set <- seq_len(K_minus + 1L)
    }
  }

  # ---- (A) PG-augmented directional scores ----
  lp_dir_exist <- vapply(seq_len(K_minus), function(kc) {
    dir_pg_SST_existing_counts(kc, c_plus, N_tot, Omega_blk, psi_minus)
  }, numeric(1))

  lp_dir_new <- rep(-Inf, K_minus + 1L)
  if (birth_score_mode == "exact_nonlocal") {
    oldold_stats <- .sst_oldold_pair_stats_excluding_i(
      i = i,
      z_packed = z_packed,
      A = A,
      i_idx = i_idx,
      j_idx = j_idx,
      N_edge = N_edge,
      omega_edge = omega_edge,
      K_minus = K_minus
    )
    lp_dir_new[r_set] <- dir_pg_SST_new_exact_nonlocal(
      c_plus = c_plus,
      N_tot = N_tot,
      Omega_blk = Omega_blk,
      psi_vec = psi_minus,
      S_old = oldold_stats$S_old,
      O_old = oldold_stats$O_old,
      r_set = r_set,
      tau0 = tau0,
      first_mu0 = mu0,
      first_sig2_0 = sig2_0
    )
  } else {
    lp_dir_new[r_set] <- dir_pg_SST_new_local_approx(
      c_plus, N_tot, Omega_blk, psi_minus, r_set, tau0
    )
  }

  # ---- (B) collapsed kappa ----
  lp_kappa_exist <- numeric(K_minus)
  for (kc in seq_len(K_minus)) {
    acc <- 0
    for (ell in seq_len(K_minus)) {
      if (r_add[ell] == 0 && t_add[ell] == 0) next
      p <- min(kc, ell); q <- max(kc, ell)
      acc <- acc + (gp_marginal(R_minus[p, q] + r_add[ell],
                                T_minus[p, q] + t_add[ell],
                                a_kappa, b_kappa) -
                      gp_marginal(R_minus[p, q], T_minus[p, q],
                                  a_kappa, b_kappa))
    }
    lp_kappa_exist[kc] <- acc
  }
  lp_kappa_new <- sum(gp_marginal(r_add, t_add, a_kappa, b_kappa))

  # ---- (C) partition prior ----
  v_minus <- n_minus_full[keep_full]
  prior_type <- if (exists(".normalize_partition_prior", mode = "function", inherits = TRUE)) {
    .normalize_partition_prior(partition_prior)
  } else {
    toupper(trimws(partition_prior))
  }
  use_ocrp <- identical(prior_type, "OCRP")
  use_rocrp <- identical(prior_type, "ROCRP")

  if (use_ocrp) {
    pw <- ocrp_log_weights_packed(v_minus, theta_ocrp = theta_ocrp)
  } else if (use_rocrp) {
    pw <- rocrp_log_weights_packed(v_minus, theta_ocrp = theta_ocrp)
  } else {
    pw <- gn_log_weights_packed(v_minus, gamma_gn = gamma_gn)
  }
  lp_prior_exist <- pw$exist
  lp_prior_new   <- pw$new
  per_slot <- isTRUE(pw$per_slot)

  total_exist <- lp_dir_exist + lp_kappa_exist + lp_prior_exist
  total_new   <- lp_dir_new   + lp_kappa_new   + lp_prior_new

  logW_join <- lse(total_exist)
  if (per_slot) {
    logW_new <- lse(total_new[r_set])
  } else {
    logW_new <- lse(total_new[r_set] - log(length(r_set)))
  }

  if (!is.finite(logW_new) && !is.finite(logW_join))
    stop("Both existing and new weights are -Inf/NaN.", call. = FALSE)

  choose_new <- (log(runif(1)) < (logW_new - lse(c(logW_new, logW_join))))

  z_new <- integer(length(z))
  z_new[-i] <- z_packed[-i]

  if (!choose_new) {
    k_star <- .softsample(total_exist)
    z_new[i] <- k_star
    return(list(z = as.integer(z_new), kappa = kappa_minus,
                psi = as.numeric(psi_minus), K = K_minus,
                changed = (K_minus != K_full),
                move = list(type = "existing", k = k_star)))
  }

  # ---- birth ----
  if (per_slot) {
    r_star <- r_set[.softsample(total_new[r_set])]
  } else {
    r_star <- r_set[.softsample(total_new[r_set] - log(length(r_set)))]
  }
  z_shift <- z_new
  z_shift[-i] <- ifelse(z_new[-i] < r_star, z_new[-i], z_new[-i] + 1L)
  z_shift[i] <- r_star

  init_kap <- rgamma(K_minus, shape = a_kappa + r_add, rate = b_kappa + t_add)

  grown <- grow_sst_params(
    kappa = kappa_minus, psi_vec = psi_minus, r = r_star,
    init_kappa_rowcol = init_kap,
    a_kappa = a_kappa, b_kappa = b_kappa, tau0 = tau0,
    mu0 = mu0, sigma0 = sqrt(sig2_0)
  )

  list(z = as.integer(z_shift), kappa = grown$kappa,
       psi = as.numeric(grown$psi), K = nrow(grown$kappa),
       changed = TRUE, move = list(type = "new", r = r_star))
}

# =============================================================================
# HYBRID SST: exact Binom‐logistic for existing blocks AND for known distances
# in new‐block proposals; PG closed‐form for the EXTREME distance only.
# =============================================================================

# --- Hybrid new-block scoring: exact binomial for d=1..K-1, PG for d=K ---
# Same as dir_exact_SST_new_vec_counts, except the extreme-distance integral
# uses the PG closed-form (log_int_pg_sst_extreme) with a plug-in Omega_star
# instead of numerical quadrature (log_int_extreme_binom_truncnorm).
#
# The plug-in Omega_star = n_K * E[PG(1, |theta_0|)] where theta_0 is the
# prior-predicted log-odds at distance K: theta_0 = psi_{K-1} + m_inc.
# When theta_0 ≈ 0, E[PG(1,0)] = 1/4, so Omega_star = n_K/4.
# Otherwise E[PG(1,c)] = tanh(c/2) / (2c).
dir_hybrid_SST_new_vec_counts <- function(c_plus, N_tot, psi_vec,
                                          r_set = NULL,
                                          tau0 = 0.15,
                                          m_inc = 0) {
  K_old <- length(N_tot)
  stopifnot(length(c_plus) == K_old)

  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)

  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(numeric(0))
  if (K_old == 1L) { lp[] <- 0.0; return(lp) }

  if (length(psi_vec) != (K_old - 1L))
    stop("SST-hybrid: psi_vec must have length K_old-1.", call. = FALSE)

  for (idx in seq_along(r_set)) {
    r <- r_set[idx]
    val <- 0.0

    A_K <- 0.0; n_K <- 0.0

    for (ell in seq_len(K_old)) {
      n <- N_tot[ell]
      if (n == 0) next

      pos_ell <- if (ell < r) ell else ell + 1L
      d <- abs(pos_ell - r)
      if (d == 0L) next

      A_fwd <- if (pos_ell > r) c_plus[ell] else (n - c_plus[ell])

      if (d <= (K_old - 1L)) {
        # --- exact binomial kernel for known distances ---
        th <- psi_vec[d]
        val <- val + (A_fwd - 0.5 * n) * th - n * .log2cosh_half(th)
      } else if (d == K_old) {
        A_K <- A_K + A_fwd
        n_K <- n_K + n
      } else {
        stop("Impossible distance in SST hybrid new-slot scoring.", call. = FALSE)
      }
    }

    # --- PG closed-form integral for extreme distance ---
    if (n_K > 0) {
      psiKm1 <- psi_vec[K_old - 1L]
      # Plug-in Omega_star: evaluate E[PG(1, |theta_0|)] at prior centre
      theta_0 <- abs(psiKm1 + m_inc)
      if (theta_0 < 1e-8) {
        Omega_star <- n_K / 4            # E[PG(1, 0)] = 1/4
      } else {
        Omega_star <- n_K * tanh(theta_0 / 2) / (2 * theta_0)
      }
      Y_K <- A_K - 0.5 * n_K
      val <- val + log_int_pg_sst_extreme(Y_K, Omega_star, psiKm1, tau0, m_inc)
    }

    lp[idx] <- val
  }
  lp
}

# --- Hybrid z-update wrapper ---
# Identical to sst_update_i_with_birth_LOO except the new-block scoring uses
# dir_hybrid_SST_new_vec_counts (PG closed-form for extreme distance only).
sst_update_i_with_birth_LOO_hybrid <- function(
    i,
    A, z, eta, kappa, psi,
    Rkl, Tkl,
    i_idx, j_idx, N_edge, edge_by_node,
    a_kappa, b_kappa,
    gamma_gn,
    r_set = NULL,
    slot_radius = NULL,
    tau0,
    mu0 = 0,
    sig2_0 = 1,
    refresh_last_delta_on_birth = TRUE,
    allow_birth = TRUE
) {
  if (!requireNamespace("truncnorm", quietly = TRUE))
    stop("Package 'truncnorm' is required.", call. = FALSE)

  if (is.function(psi)) stop("psi is a function.", call. = FALSE)
  if (!is.numeric(psi)) stop("psi must be numeric.", call. = FALSE)

  K_full <- nrow(kappa)
  if (ncol(kappa) != K_full) stop("kappa not square.", call. = FALSE)
  if (length(psi) != max(K_full - 1L, 0L))
    stop(sprintf("SST invariant: length(psi)=%d but K-1=%d",
                 length(psi), max(K_full - 1L, 0L)), call. = FALSE)

  u_full <- sort(unique(z))
  if (!identical(u_full, seq_len(K_full)))
    stop(sprintf("SST invariant: labels not contiguous: unique(z)=%s, expected 1..%d",
                 paste(u_full, collapse = ","), K_full), call. = FALSE)

  oldk  <- z[i]
  eta_i <- eta[i]

  # ---- LOO occupancy and packing ----
  n_minus_full <- tabulate(z[-i], nbins = K_full)
  keep_full <- sort(which(n_minus_full > 0L))
  K_minus <- length(keep_full)
  if (K_minus < 1L) stop("LOO state has K_minus<1.", call. = FALSE)

  map <- integer(K_full); map[keep_full] <- seq_len(K_minus)
  z_packed <- map[z]
  oldk_pos <- match(oldk, keep_full)

  kappa_minus <- kappa[keep_full, keep_full, drop = FALSE]
  psi_minus   <- reindex_psi_sst_keep(psi_old = psi, K_old = K_full, keep = keep_full)
  if (length(psi_minus) != max(K_minus - 1L, 0L))
    stop("psi_minus length mismatch after packing.", call. = FALSE)

  # ---- node-vs-block counts ----
  C_i_full <- counts_by_block_exact_cpp(
    i = i, A = A, z = as.integer(z),
    i_idx = as.integer(i_idx), j_idx = as.integer(j_idx),
    N_edge = as.numeric(N_edge),
    edge_by_node = edge_by_node, K = K_full
  )
  c_plus <- as.numeric(C_i_full$c_plus)[keep_full]
  N_tot  <- as.numeric(C_i_full$N_tot)[keep_full]

  E_excl_full <- E_by_block_excluding_i(i, z, eta, K_full)
  E_excl <- as.numeric(E_excl_full[keep_full])
  r_add <- N_tot; t_add <- eta_i * E_excl

  # ---- LOO base R/T ----
  R_minus <- Rkl[keep_full, keep_full, drop = FALSE]
  T_minus <- Tkl[keep_full, keep_full, drop = FALSE]
  if (!is.na(oldk_pos)) {
    for (ell_pos in seq_len(K_minus)) {
      subR <- N_tot[ell_pos]; subT <- eta_i * E_excl[ell_pos]
      p <- min(oldk_pos, ell_pos); q <- max(oldk_pos, ell_pos)
      R_minus[p, q] <- R_minus[p, q] - subR
      T_minus[p, q] <- T_minus[p, q] - subT
      R_minus[q, p] <- R_minus[p, q]; T_minus[q, p] <- T_minus[p, q]
      if (R_minus[p, q] < 0) R_minus[p, q] <- R_minus[q, p] <- 0
      if (T_minus[p, q] < 0) T_minus[p, q] <- T_minus[q, p] <- 0
    }
  }

  # ---- candidate slots ----
  if (is.null(r_set)) {
    if (!is.null(slot_radius) && !is.na(oldk_pos)) {
      lo <- max(1L, oldk_pos - (slot_radius - 1L))
      hi <- min(K_minus + 1L, oldk_pos + slot_radius)
      r_set <- seq.int(lo, hi)
    } else {
      r_set <- seq_len(K_minus + 1L)
    }
  }

  # ---- (A) directional scores: EXACT existing, HYBRID new ----
  lp_dir_exist <- vapply(seq_len(K_minus), function(kc) {
    dir_exact_SST_existing_counts(kc, c_plus, N_tot, psi_minus)
  }, numeric(1))

  lp_dir_new <- rep(-Inf, K_minus + 1L)
  lp_dir_new[r_set] <- dir_hybrid_SST_new_vec_counts(
    c_plus, N_tot, psi_minus, r_set, tau0
  )

  # ---- (B) collapsed kappa ----
  lp_kappa_exist <- numeric(K_minus)
  for (kc in seq_len(K_minus)) {
    acc <- 0
    for (ell in seq_len(K_minus)) {
      if (r_add[ell] == 0 && t_add[ell] == 0) next
      p <- min(kc, ell); q <- max(kc, ell)
      acc <- acc + (gp_marginal(R_minus[p, q] + r_add[ell],
                                T_minus[p, q] + t_add[ell],
                                a_kappa, b_kappa) -
                      gp_marginal(R_minus[p, q], T_minus[p, q],
                                  a_kappa, b_kappa))
    }
    lp_kappa_exist[kc] <- acc
  }
  lp_kappa_new <- sum(gp_marginal(r_add, t_add, a_kappa, b_kappa))

  # ---- (C) GN prior ----
  v_minus <- n_minus_full[keep_full]
  pw <- gn_log_weights_packed(v_minus, gamma_gn = gamma_gn)
  lp_prior_exist <- pw$exist
  lp_prior_new   <- pw$new

  total_exist <- lp_dir_exist + lp_kappa_exist + lp_prior_exist
  total_new   <- lp_dir_new   + lp_kappa_new   + lp_prior_new

  logW_join <- lse(total_exist)
  logW_new  <- lse(total_new[r_set] - log(length(r_set)))

  choose_new <- (log(runif(1)) < (logW_new - lse(c(logW_new, logW_join))))
  if (!isTRUE(allow_birth)) choose_new <- FALSE

  z_new <- integer(length(z)); z_new[-i] <- z_packed[-i]

  if (!choose_new) {
    k_star <- .softsample(total_exist)
    z_new[i] <- k_star
    return(list(z = as.integer(z_new), kappa = kappa_minus,
                psi = as.numeric(psi_minus), K = K_minus,
                changed = (K_minus != K_full),
                move = list(type = "existing", k = k_star)))
  }

  # ---- birth ----
  r_star <- r_set[.softsample(total_new[r_set] - log(length(r_set)))]
  z_shift <- z_new
  z_shift[-i] <- ifelse(z_new[-i] < r_star, z_new[-i], z_new[-i] + 1L)
  z_shift[i] <- r_star

  init_kap <- rgamma(K_minus, shape = a_kappa + r_add, rate = b_kappa + t_add)

  grown <- grow_sst_params(
    kappa = kappa_minus, psi_vec = psi_minus, r = r_star,
    init_kappa_rowcol = init_kap,
    a_kappa = a_kappa, b_kappa = b_kappa, tau0 = tau0
  )

  # ---- refresh last increment (same as collapsed) ----
  if (isTRUE(refresh_last_delta_on_birth)) {
    K_new <- nrow(grown$kappa); D_new <- K_new - 1L
    if (D_new >= 1L) {
      A_ij <- A[cbind(i_idx, j_idx)]
      zi <- z_shift[i_idx]; zj <- z_shift[j_idx]
      agg <- aggregate_by_distance(K_new, z_i = zi, z_j = zj,
                                   A_ij = A_ij, N_edge = N_edge)
      B_d <- distance_totals(K_new, z_i = zi, z_j = zj, N_edge = N_edge)
      bar_omega <- draw_omega_bar(B_d = B_d, psi = grown$psi)
      Tmat  <- make_T(K_new)
      Omega <- diag(bar_omega, nrow = D_new, ncol = D_new)
      Q     <- crossprod(Tmat, Omega %*% Tmat)
      g     <- crossprod(Tmat, agg$bar_y)
      Vinv  <- diag(c(1 / sig2_0, rep(1 / (tau0^2), max(0, D_new - 1))),
                    D_new, D_new)
      muvec <- c(mu0, rep(0, max(0, D_new - 1)))
      Qstar <- Q + Vinv; gstar <- g + Vinv %*% muvec
      delta <- numeric(D_new)
      delta[1] <- max(grown$psi[1], 1e-10)
      if (D_new >= 2L)
        delta[2:D_new] <- pmax(grown$psi[2:D_new] - grown$psi[1:(D_new - 1L)], 1e-10)
      r <- D_new
      a_rr <- Qstar[r, r]
      b_r  <- as.numeric(gstar[r] - sum(Qstar[r, -r, drop = FALSE] * delta[-r]))
      mean_r <- b_r / a_rr; sd_r <- 1 / sqrt(a_rr)
      delta[r] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = mean_r, sd = sd_r)
      if (!is.finite(delta[r]) || delta[r] <= 0) delta[r] <- 1e-10
      grown$psi <- cumsum(delta)
    }
  }

  list(z = as.integer(z_shift), kappa = grown$kappa,
       psi = as.numeric(grown$psi), K = nrow(grown$kappa),
       changed = TRUE, move = list(type = "new", r = r_star))
}

dir_exact_SST_existing_counts <- function(k, c_plus, N_tot, psi_vec) {
  K <- length(N_tot)
  stopifnot(length(c_plus) == K)
  if (K <= 1L) return(0)
  
  # within-block term (theta = 0 => p = 1/2), applied ONCE
  val <- - N_tot[k] * log(2)
  
  for (ell in seq_len(K)) {
    if (ell == k) next
    n <- N_tot[ell]
    if (n == 0) next
    
    d  <- abs(ell - k)
    th <- psi_vec[d]
    
    # forward = lower label -> higher label
    A_fwd <- if (ell > k) c_plus[ell] else (n - c_plus[ell])
    
    val <- val + (A_fwd - 0.5 * n) * th - n * .log2cosh_half(th)
  }
  val
}


shrink_psi_sst <- function(psi_old, K_new) {
  D_new <- max(K_new - 1L, 0L)
  if (D_new == 0L) return(numeric(0))
  psi_old[seq_len(D_new)]
}

dir_exact_SST_new_vec_counts <- function(c_plus, N_tot, psi_vec,
                                         r_set = NULL,
                                         tau0 = 0.15,
                                         m_inc = 0) {
  K_old <- length(N_tot)
  stopifnot(length(c_plus) == K_old)
  
  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)
  
  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(numeric(0))
  if (K_old == 1L) { lp[] <- 0.0; return(lp) }
  
  if (length(psi_vec) != (K_old - 1L)) {
    stop("SST: psi_vec must have length K_old-1.", call. = FALSE)
  }
  
  for (idx in seq_along(r_set)) {
    r <- r_set[idx]
    val <- 0.0
    
    # collect the farthest-distance (d = K_old) block contribution, if it exists
    A_K <- 0.0
    n_K <- 0.0
    
    for (ell in seq_len(K_old)) {
      n <- N_tot[ell]
      if (n == 0) next
      
      pos_ell <- if (ell < r) ell else ell + 1L
      d <- abs(pos_ell - r)
      if (d == 0L) next
      
      # forward count lower->higher
      A_fwd <- if (pos_ell > r) c_plus[ell] else (n - c_plus[ell])
      
      if (d <= (K_old - 1L)) {
        th <- psi_vec[d]
        val <- val + (A_fwd - 0.5 * n) * th - n * .log2cosh_half(th)
      } else if (d == K_old) {
        A_K <- A_K + A_fwd
        n_K <- n_K + n
      } else {
        stop("Impossible distance in SST new-slot scoring.", call. = FALSE)
      }
    }
    
    # exact integral for ψ_K = ψ_{K-1} + δ, δ ~ N^+(m_inc, tau0^2)
    if (n_K > 0) {
      psiKm1 <- psi_vec[K_old - 1L]
      val <- val + log_int_extreme_binom_truncnorm(
        A = A_K, n = n_K, psiKm1 = psiKm1, tau0 = tau0, m_inc = m_inc
      )
    }
    
    lp[idx] <- val
  }
  
  lp
}

# ------- SST: extreme-slot integral for new-block directional score -----------

log_int_extreme_binom_truncnorm <- function(A, n, psiKm1, tau0, m_inc = 0,
                                            rel.tol = 1e-10) {
  stopifnot(n >= 0, A >= 0, A <= n, tau0 > 0)
  if (n == 0) return(0)
  
  # normalising constant for truncation δ > 0
  logZ <- pnorm(m_inc / tau0, log.p = TRUE)  # log Φ(m/τ)
  
  logf <- function(delta) {
    th <- psiKm1 + delta
    A * logsig(th) + (n - A) * log1m_sig(th) +
      dnorm(delta, mean = m_inc, sd = tau0, log = TRUE) - logZ
  }
  
  # Find a decent stabilising offset (mode-ish) on a wide interval
  ub <- 40  # logistic saturates long before this; safe upper search
  opt <- optimize(function(d) -logf(d), interval = c(0, ub))
  m0 <- max(logf(0), logf(opt$minimum), logf(ub))
  
  g <- function(delta) exp(logf(delta) - m0)
  out <- integrate(g, lower = 0, upper = Inf, rel.tol = rel.tol)
  
  if (!is.finite(out$value) || out$value <= 0) {
    stop("Extreme-slot integral failed (non-finite / non-positive).", call. = FALSE)
  }
  log(out$value) + m0
}


# Exact SST directional log-lik for inserting a NEW block at slot r
# Returns length-(K+1) vector over r = 1..K+1 unless r_set is supplied
# --- NEW: exact SST directional log-lik for inserting a new block at slot r
#     with proper extreme-slot integral for distance K_old via sst_increment_fast()
#
# c_plus[ell] = i -> block ell wins (as oriented by counts_by_block_exact_cpp)
# N_tot[ell]  = total matches vs block ell
# psi_vec length = K_old-1, psi_vec[d] = ψ_d (nondecreasing, >=0)
# tau0 = SD of new increment δ_{K-1} ~ N^+(0, tau0^2)
#


#------------------------------------------------------
# Directional: existing-slot PG surrogate using pooled Y/Ω
# pools$Y_blk and pools$Omega_blk are per old block ℓ around i
# --------------------------------------------------------------
dir_pg_SST_existing <- function(k, Y_blk, Omega_blk, psi_vec) {
  K_old <- length(Y_blk)
  lp <- 0
  for (ell in seq_len(K_old)) {
    if (ell == k) next
    d  <- abs(ell - k)
    if (d == 0L) next
    sgn <- if (ell > k) +1L else -1L
    th  <- theta_at(+1L, d, psi_vec)      # ψ_d ≥ 0
    lp  <- lp + (sgn * Y_blk[ell]) * th - 0.5 * Omega_blk[ell] * th^2
  }
  lp
}

.check_sst_state <- function(stage, it, i = NA_integer_,
                             z, kappa, psi, i_idx, j_idx,
                             dump_dir = NULL) {
  K <- nrow(kappa)
  D <- max(K - 1L, 0L)
  
  fail <- function(msg) {
    if (!is.null(dump_dir)) {
      dir.create(dump_dir, showWarnings = FALSE, recursive = TRUE)
      fn <- file.path(dump_dir, sprintf("FAIL_%s_it%05d_i%05d.rds", stage, it, i %||% -1L))
      saveRDS(list(stage = stage, it = it, i = i, K = K, D = D,
                   z = z, psi = psi, kappa = kappa),
              file = fn)
      message("State dumped to: ", fn)
    }
    stop(msg, call. = FALSE)
  }
  
  if (ncol(kappa) != K) fail(sprintf("[%s|it=%d|i=%s] kappa not square", stage, it, i))
  
  if (length(psi) != D) {
    fail(sprintf("[%s|it=%d|i=%s] length(psi)=%d, expected %d (K=%d)",
                 stage, it, i, length(psi), D, K))
  }
  
  if (anyNA(psi) || any(!is.finite(psi)) || any(psi < 0)) {
    fail(sprintf("[%s|it=%d|i=%s] psi has NA/Inf/neg (NA=%d, min=%g, max=%g)",
                 stage, it, i, sum(is.na(psi)), min(psi, na.rm=TRUE), max(psi, na.rm=TRUE)))
  }
  
  if (anyNA(z) || any(z < 1L) || any(z > K)) {
    fail(sprintf("[%s|it=%d|i=%s] z outside 1..K or NA (NA=%d, min=%d, max=%d, K=%d)",
                 stage, it, i, sum(is.na(z)), min(z, na.rm=TRUE), max(z, na.rm=TRUE), K))
  }

  u <- sort(unique(z))
  if (!identical(u, seq_len(K))) {
    fail(sprintf("[%s|it=%d|i=%s] labels not contiguous: unique(z)=%s, expected 1..%d",
                 stage, it, i, paste(u, collapse=","), K))
  }
  
  if (length(i_idx)) {
    d <- abs(z[i_idx] - z[j_idx])
    md <- max(d, na.rm = TRUE)
    if (md > D) {
      fail(sprintf("[%s|it=%d|i=%s] max distance=%d exceeds D=%d (K=%d)",
                   stage, it, i, md, D, K))
    }
  }
  
  invisible(TRUE)
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


# ===================================================================
# UPDATED BIRTH UPDATER (SST / "distance") — with PARTIALLY COLLAPSED κ
# ===================================================================
# What changed vs the uncollapsed version:
# - The existing-block κ contribution in the z-update now integrates out κ_{kℓ}
#   using block-pair sufficient statistics with node i removed.
# - This removes dependence on the current sampled κ in the z step.
#
# NOTE (important for correctness of partially-collapsed Gibbs):
# - Do NOT resample κ in between the single-site z_i updates. Update all z's first,
#   then draw κ | z,rest as usual.

# ---- Gamma–Poisson marginal kernel for a single block-pair (κ-dependent part)
# Model: N_{uv} | κ ~ Poisson(η_u η_v κ), κ ~ Gamma(a, b)   (shape a, rate b)
# Let R = sum N_{uv} over dyads in the pair, T = sum η_u η_v over those dyads.
# Then (up to constants not depending on κ / labels):
#   log m(R,T) = log Γ(a+R) - log Γ(a) + a log b - (a+R) log(b+T).
# gp_marginal <- function(R, T, a, b) {
#   stopifnot(all(R >= 0), all(T >= 0), a > 0, b > 0)
#   lgamma(a + R) - lgamma(a) + a * log(b) - (a + R) * log(b + T)
# }

## Partially-collapsed kappa term for z-update, with correct upper-triangular indexing

# ---------- helpers ----------



# Access an upper-triangular (incl diagonal) matrix M storing stats for unordered block-pairs:
# M[p,q] with p = min(k,l), q = max(k,l)
tri_get <- function(M, k, l) {
  p <- if (k < l) k else l
  q <- if (k < l) l else k
  M[p, q]
}
tri_set <- function(M, k, l, value) {
  p <- if (k < l) k else l
  q <- if (k < l) l else k
  M[p, q] <- value
  M
}

# ---------- main updater (patched) ----------
# This function assumes the same arguments as your existing z_update_osbm_distance_with_birth().
# If your signature differs, tell me and I’ll adapt it exactly.
#
# Key required objects:
# - A: adjacency/interaction matrix (counts), n x n
# - z: current block assignments, length n, integer in 1:K
# - eta: node propensities, length n
# - Rkl: upper-triangular matrix of block-pair counts (sum of A_{uv} over unordered pairs)
# - Tkl: upper-triangular matrix of block-pair exposures (sum of eta_u * eta_v over unordered pairs)
# - a_kappa, b_kappa: Gamma prior hyperparameters for kappa
# - log_dist_terms: whatever you already compute for the distance/CRP part (left untouched here)
# - plus whatever you need for directional term (left untouched)
#
# IMPORTANT: This patch only replaces the "existing blocks: kappa contribution" part.
#


# ---- patched updater ------------------------------------------------------
sst_update_i_with_birth_LOO <- function(
    i,
    A, z, eta, kappa, psi,
    Rkl, Tkl,
    i_idx, j_idx, N_edge, edge_by_node,
    a_kappa, b_kappa,
    gamma_gn,
    r_set = NULL,
    slot_radius = NULL,
    tau0,
    mu0 = 0,
    sig2_0 = 1,
    refresh_last_delta_on_birth = TRUE,
    # When FALSE the function behaves as a fixed-K sampler: the birth weights
    # are computed (for diagnostics) but choose_new is forced to FALSE.
    allow_birth = TRUE
) {
  if (!requireNamespace("truncnorm", quietly = TRUE))
    stop("Package 'truncnorm' is required.", call. = FALSE)
  
  # ---- fail-fast invariants on entry (full state) ----
  if (is.function(psi)) stop("psi is a function (digamma). You lost your numeric psi vector.", call. = FALSE)
  if (!is.numeric(psi)) stop("psi must be numeric.", call. = FALSE)
  
  K_full <- nrow(kappa)
  if (ncol(kappa) != K_full) stop("kappa not square.", call. = FALSE)
  if (length(psi) != max(K_full - 1L, 0L))
    stop(sprintf("SST invariant: length(psi)=%d but K-1=%d", length(psi), max(K_full-1L,0L)),
         call. = FALSE)
  
  u_full <- sort(unique(z))
  if (!identical(u_full, seq_len(K_full)))
    stop(sprintf("SST invariant: labels not contiguous: unique(z)=%s, expected 1..%d",
                 paste(u_full, collapse=","), K_full), call. = FALSE)
  
  oldk  <- z[i]
  eta_i <- eta[i]
  
  # ---- LOO occupancy and packing map ----
  n_minus_full <- tabulate(z[-i], nbins = K_full)
  keep_full <- sort(which(n_minus_full > 0L))
  K_minus <- length(keep_full)
  if (K_minus < 1L) stop("LOO state has K_minus<1 (n too small?).", call. = FALSE)
  
  map <- integer(K_full)
  map[keep_full] <- seq_len(K_minus)
  
  z_packed <- map[z]
  if (any(z_packed[-i] < 1L | z_packed[-i] > K_minus))
    stop("Packing produced out-of-range labels.", call. = FALSE)
  
  oldk_pos <- match(oldk, keep_full)
  
  # ---- pack parameters to the LOO occupied blocks ----
  kappa_minus <- kappa[keep_full, keep_full, drop = FALSE]
  psi_minus   <- reindex_psi_sst_keep(psi_old = psi, K_old = K_full, keep = keep_full)
  
  if (length(psi_minus) != max(K_minus - 1L, 0L))
    stop("psi_minus length mismatch after packing.", call. = FALSE)
  
  # ---- node-vs-block directed counts (FULL indexing), then subset ----
  C_i_full <- counts_by_block_exact_cpp(
    i = i, A = A, z = as.integer(z),
    i_idx = as.integer(i_idx), j_idx = as.integer(j_idx),
    N_edge = as.numeric(N_edge),
    edge_by_node = edge_by_node,
    K = K_full
  )
  c_plus <- as.numeric(C_i_full$c_plus)[keep_full]
  N_tot  <- as.numeric(C_i_full$N_tot )[keep_full]
  
  # ---- exposure sums excluding i (FULL indexing), then subset ----
  E_excl_full <- E_by_block_excluding_i(i, z, eta, K_full)
  E_excl <- as.numeric(E_excl_full[keep_full])
  
  r_add <- N_tot
  t_add <- eta_i * E_excl
  
  # ---- LOO base R/T among other nodes in packed indexing ----
  R_minus <- Rkl[keep_full, keep_full, drop = FALSE]
  T_minus <- Tkl[keep_full, keep_full, drop = FALSE]
  
  # subtract i's contribution if its old block survives in LOO
  if (!is.na(oldk_pos)) {
    for (ell_pos in seq_len(K_minus)) {
      subR <- N_tot[ell_pos]
      subT <- eta_i * E_excl[ell_pos]
      p <- min(oldk_pos, ell_pos)
      q <- max(oldk_pos, ell_pos)
      
      R_minus[p, q] <- R_minus[p, q] - subR
      T_minus[p, q] <- T_minus[p, q] - subT
      R_minus[q, p] <- R_minus[p, q]
      T_minus[q, p] <- T_minus[p, q]
      
      if (R_minus[p, q] < -1e-8 || T_minus[p, q] < -1e-8)
        stop(sprintf("Negative LOO totals after subtracting i: R=%g T=%g (p=%d,q=%d)",
                     R_minus[p,q], T_minus[p,q], p, q), call. = FALSE)
      
      if (R_minus[p, q] < 0) R_minus[p, q] <- R_minus[q, p] <- 0
      if (T_minus[p, q] < 0) T_minus[p, q] <- T_minus[q, p] <- 0
    }
  }
  
  # ---- candidate slots in LOO packed indexing ----
  if (is.null(r_set)) {
    if (!is.null(slot_radius) && !is.na(oldk_pos)) {
      lo <- max(1L, oldk_pos - (slot_radius - 1L))
      hi <- min(K_minus + 1L, oldk_pos + slot_radius)
      r_set <- seq.int(lo, hi)
    } else {
      r_set <- seq_len(K_minus + 1L)
    }
  } else {
    r_set <- as.integer(r_set)
    if (min(r_set) < 1L || max(r_set) > (K_minus + 1L))
      stop("r_set must lie in 1:(K_minus+1) for the LOO state.", call. = FALSE)
  }
  
  # ---- (A) directional scores ----
  lp_dir_exist <- vapply(seq_len(K_minus), function(kc) {
    dir_exact_SST_existing_counts(k = kc, c_plus = c_plus, N_tot = N_tot, psi_vec = psi_minus)
  }, numeric(1))
  
  lp_dir_new <- rep(-Inf, K_minus + 1L)
  lp_dir_new[r_set] <- dir_exact_SST_new_vec_counts(
    c_plus = c_plus, N_tot = N_tot, psi_vec = psi_minus, r_set = r_set,
    tau0 = tau0
  )
  
  # ---- (B) collapsed kappa scores ----
  lp_kappa_exist <- numeric(K_minus)
  for (kc in seq_len(K_minus)) {
    acc <- 0
    for (ell in seq_len(K_minus)) {
      if (r_add[ell] == 0 && t_add[ell] == 0) next
      p <- min(kc, ell)
      q <- max(kc, ell)
      acc <- acc + (gp_marginal(R_minus[p,q] + r_add[ell], T_minus[p,q] + t_add[ell], a_kappa, b_kappa) -
                      gp_marginal(R_minus[p,q],         T_minus[p,q],         a_kappa, b_kappa))
    }
    lp_kappa_exist[kc] <- acc
  }
  lp_kappa_new <- sum(gp_marginal(r_add, t_add, a_kappa, b_kappa))
  
  # ---- (C) GN prior weights (packed) ----
  v_minus <- n_minus_full[keep_full]
  pw <- gn_log_weights_packed(v_minus, gamma_gn = gamma_gn)
  lp_prior_exist <- pw$exist
  lp_prior_new   <- pw$new
  
  # ---- (D) sample move ----
  total_exist <- lp_dir_exist + lp_kappa_exist + lp_prior_exist
  total_new   <- lp_dir_new   + lp_kappa_new   + lp_prior_new
  
  logW_join <- lse(total_exist)
  logW_new  <- lse(total_new[r_set] - log(length(r_set)))

  choose_new <- (log(runif(1)) < (logW_new - lse(c(logW_new, logW_join))))
  if (!isTRUE(allow_birth)) choose_new <- FALSE  # fixed-K warm-up

  # ---- (E) apply move and return UPDATED STATE ----
  z_new <- integer(length(z))
  z_new[-i] <- z_packed[-i]
  
  if (!choose_new) {
    k_star <- .softsample(total_exist)
    z_new[i] <- k_star
    return(list(
      z = as.integer(z_new),
      kappa = kappa_minus,
      psi = as.numeric(psi_minus),
      K = K_minus,
      changed = (K_minus != K_full),
      move = list(type = "existing", k = k_star)
    ))
  }
  
  # ---------- birth ----------
  r_star <- r_set[ .softsample(total_new[r_set] - log(length(r_set))) ]
  
  z_shift <- z_new
  z_shift[-i] <- ifelse(z_new[-i] < r_star, z_new[-i], z_new[-i] + 1L)
  z_shift[i] <- r_star
  
  init_kap <- rgamma(
    K_minus,
    shape = a_kappa + r_add,
    rate  = b_kappa + t_add
  )
  
  grown <- grow_sst_params(
    kappa = kappa_minus,
    psi_vec = psi_minus,
    r = r_star,
    init_kappa_rowcol = init_kap,
    a_kappa = a_kappa,
    b_kappa = b_kappa,
    tau0 = tau0
  )
  
  # ---- refresh ONLY the new last increment δ_D using its TN full conditional ----
  if (isTRUE(refresh_last_delta_on_birth)) {
    K_new <- nrow(grown$kappa)
    D_new <- K_new - 1L
    if (D_new >= 1L) {
      A_ij <- A[cbind(i_idx, j_idx)]
      zi <- z_shift[i_idx]
      zj <- z_shift[j_idx]
      
      # global aggregated stats by distance under the new ordering
      agg <- aggregate_by_distance(K_new, z_i = zi, z_j = zj, A_ij = A_ij, N_edge = N_edge)
      bar_y <- agg$bar_y
      B_d   <- distance_totals(K_new, z_i = zi, z_j = zj, N_edge = N_edge)
      bar_omega <- draw_omega_bar(B_d = B_d, psi = grown$psi)
      
      # build Q*, g* as in update_psi_sst
      Tmat  <- make_T(K_new)
      Omega <- diag(bar_omega, nrow = D_new, ncol = D_new)
      Q     <- crossprod(Tmat, Omega %*% Tmat)
      g     <- crossprod(Tmat, bar_y)
      
      Vinv  <- diag(c(1/sig2_0, rep(1/(tau0^2), max(0, D_new-1))), D_new, D_new)
      muvec <- c(mu0, rep(0, max(0, D_new-1)))
      
      Qstar <- Q + Vinv
      gstar <- g + Vinv %*% muvec
      
      # current delta from psi
      psi_curr <- grown$psi
      delta <- numeric(D_new)
      delta[1] <- max(psi_curr[1], 1e-10)
      if (D_new >= 2L) {
        delta[2:D_new] <- pmax(psi_curr[2:D_new] - psi_curr[1:(D_new-1L)], 1e-10)
      }
      
      r <- D_new
      a_rr <- Qstar[r, r]
      b_r  <- as.numeric(gstar[r] - sum(Qstar[r, -r, drop = FALSE] * delta[-r]))
      mean_r <- b_r / a_rr
      sd_r   <- 1 / sqrt(a_rr)
      
      delta[r] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = mean_r, sd = sd_r)
      if (!is.finite(delta[r]) || delta[r] <= 0) delta[r] <- 1e-10
      
      psi_ref <- as.numeric(Tmat %*% delta)
      grown$psi <- cummax(pmax(psi_ref, 0))
    }
  }
  
  list(
    z = as.integer(z_shift),
    kappa = grown$kappa,
    psi = as.numeric(grown$psi),
    K = nrow(grown$kappa),
    changed = TRUE,
    move = list(type = "new", r = r_star)
  )
}



reindex_psi_sst_keep <- function(psi_old, K_old, keep) {
  keep <- sort(as.integer(keep))
  stopifnot(K_old >= 1L, all(keep %in% seq_len(K_old)))
  
  K_new <- length(keep)
  if (K_new <= 1L) return(numeric(0))
  
  D_old <- K_old - 1L
  if (length(psi_old) != D_old)
    stop(sprintf("reindex_psi_sst_keep: length(psi_old)=%d but K_old-1=%d", length(psi_old), D_old),
         call. = FALSE)
  
  # enforce monotone via delta representation
  delta_old <- numeric(D_old)
  delta_old[1] <- max(psi_old[1], 0)
  if (D_old >= 2L) delta_old[2:D_old] <- pmax(psi_old[2:D_old] - psi_old[1:(D_old - 1L)], 0)
  psi_mon <- cumsum(delta_old)
  
  D_new <- K_new - 1L
  psi_new <- numeric(D_new)
  
  for (d in seq_len(D_new)) {
    old_dists <- keep[(1L + d):K_new] - keep[1L:(K_new - d)]
    old_dists <- pmax(pmin(old_dists, D_old), 1L)
    psi_new[d] <- mean(psi_mon[old_dists])
  }
  
  cummax(psi_new)
}






# ====== NEW-BLOCK: SST directional score (distance ψ) =======================
# collapsed directional score for NEW block under SST (slot r)
# - N_blk[ell]: dyad totals between node i and block ell
# - A_out_blk[ell]: forward counts i->block(ell)
# - psi_vec: length K-1 (ψ_1,...,ψ_{K-1})
# - psi_hyper$tau0 is the SD for the increment prior δ_{K-1} ~ N^+(0, τ0^2)
log_dirlik_new_sst <- function(r, Y_blk, Omega_blk, psi_vec, psi_hyper) {
  # K_old = current number of blocks; psi_vec has length K_old - 1
  K <- length(psi_vec) + 1L
  stopifnot(length(Y_blk) == K, length(Omega_blk) == K)
  
  log_sum <- 0
  
  # Prior for the new increment δ_{K-1} (distance K): N^+(m_inc, v_inc)
  # In the appendix we take m_inc = 0 and variance τ0^2.
  tau0 <- if (!is.null(psi_hyper$tau0)) psi_hyper$tau0 else psi_hyper$sigma0
  cache_inc <- make_sst_cache(m_inc = 0, tau0 = tau0)
  
  # Aggregate PG statistics by distance d after inserting at slot r
  # Y_by_d[d]    = sum s_{sd} * Y_{ij} for all edges from i to blocks at distance d
  # Omega_by_d[d]= sum ω_{ij} for same edges
  Y_by_d     <- numeric(K)   # indices 1..K; d = K only at extremes
  Omega_by_d <- numeric(K)
  
  for (ell in seq_len(K)) {
    if (Y_blk[ell] == 0 && Omega_blk[ell] == 0) next
    
    # position of block ℓ after inserting the new block at slot r
    pos_ell <- if (ell < r) ell else ell + 1L
    d <- abs(pos_ell - r)
    if (d == 0L) next
    
    sgn <- if (pos_ell > r) +1L else -1L
    
    if (d <= K - 1L) {
      Y_by_d[d]     <- Y_by_d[d]     + sgn * Y_blk[ell]
      Omega_by_d[d] <- Omega_by_d[d] + Omega_blk[ell]
    } else if (d == K) {
      Y_by_d[K]     <- Y_by_d[K]     + sgn * Y_blk[ell]
      Omega_by_d[K] <- Omega_by_d[K] + Omega_blk[ell]
    } else {
      stop("log_dirlik_new_sst: impossible distance in SST")
    }
  }
  
  # Reuse existing ψ_d for distances d = 1,…,K−1
  if (K > 1L) {
    for (d in seq_len(K - 1L)) {
      if (Omega_by_d[d] <= 0 && Y_by_d[d] == 0) next
      th <- psi_vec[d]
      log_sum <- log_sum + (Y_by_d[d] * th - 0.5 * Omega_by_d[d] * th^2)
    }
  }
  
  # Distance K only appears when r is an extreme slot; integrate δ_{K-1}
  if (Omega_by_d[K] > 0 || Y_by_d[K] != 0) {
    psiKm1 <- if (K >= 2L) psi_vec[K - 1L] else 0
    log_sum <- log_sum + sst_increment_fast(
      yK        = Y_by_d[K],
      omegaK    = Omega_by_d[K],
      psiKm1    = psiKm1,
      cache     = cache_inc,
      return_log = TRUE
    )
  }
  
  as.numeric(log_sum)
}



# ------- SST: totals and PG draws for the ψ-update --------------------------
distance_totals <- function(K, z_i, z_j, N_edge) {
  B_d <- numeric(max(K - 1L, 0))
  for (e in seq_along(N_edge)) {
    d <- abs(z_i[e] - z_j[e])
    if (d == 0L) next
    B_d[d] <- B_d[d] + N_edge[e]
  }
  B_d
}

# --- grow κ and ψ for SST (distance ψ_d) ------------------------------------

grow_sst_params <- function(kappa, psi_vec, r, init_kappa_rowcol,
                            a_kappa = 1, b_kappa = 1,
                            tau0 = 0.15,
                            mu0 = 0, sigma0 = NULL) {
  K_old <- nrow(kappa)
  stopifnot(ncol(kappa) == K_old,
            length(psi_vec) == (K_old - 1L),
            length(init_kappa_rowcol) == K_old)
  
  K_new <- K_old + 1L
  map_new <- function(o) if (o < r) o else (o + 1L)
  
  # κ expand: carry old
  kap2 <- matrix(0, K_new, K_new)
  for (a in seq_len(K_old)) for (b in seq_len(K_old)) {
    kap2[map_new(a), map_new(b)] <- kappa[a, b]
  }
  # newborn row/col from posterior draws passed in
  for (ell in seq_len(K_old)) {
    j <- map_new(ell)
    kap2[r, j] <- kap2[j, r] <- init_kappa_rowcol[ell]
  }
  # within-new (unused until >=2 nodes): prior draw
  kap2[r, r] <- rgamma(1, shape = a_kappa, rate = b_kappa)
  
  # ψ extend by one. The first SST distance uses the δ_1 prior; later births
  # add the new extreme-distance increment with the τ prior.
  last <- if (length(psi_vec)) psi_vec[length(psi_vec)] else 0
  inc_mean <- if (length(psi_vec)) 0 else mu0
  inc_sd <- if (length(psi_vec)) tau0 else (sigma0 %||% tau0)
  delta <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = inc_mean, sd = inc_sd)
  if (!is.finite(delta) || delta <= 0) delta <- 1e-10
  psi_new <- c(psi_vec, last + delta)
  
  
  list(kappa = kap2, psi = psi_new)
}


# helpers for candidate slot r (SST new)
pos_r <- function(ell, r) if (ell < r) ell else ell + 1L
sgn_r <- function(ell, r) if (ell >= r) +1L else -1L


# cache for SST increment prior
make_sst_cache <- function(m_inc = 0, v_inc = NULL, tau0 = NULL, mu0 = NULL, sigma0 = NULL) {
  if (!is.null(tau0)) v_inc <- tau0^2
  if (is.null(v_inc) && !is.null(sigma0)) v_inc <- sigma0^2
  if (is.null(v_inc)) v_inc <- 1
  stopifnot(is.finite(m_inc), is.finite(v_inc), v_inc > 0)
  list(
    m_inc = m_inc,
    v_inc = v_inc,
    inv_v = 1 / v_inc,
    logZ_prior = pnorm(m_inc / sqrt(v_inc), log.p = TRUE)  # log Φ(m/√v)
  )
}
sst_increment_fast <- function(yK, omegaK, psiKm1, cache, return_log = TRUE) {
  m <- cache$m_inc
  v <- cache$v_inc
  inv_v <- cache$inv_v
  logZ_prior <- cache$logZ_prior
  
  lambda <- omegaK + inv_v
  b <- yK - omegaK * psiKm1 + inv_v * m
  
  logPhi_post <- pnorm(b / sqrt(lambda), log.p = TRUE)
  
  logI <- -0.5 * log(lambda) - 0.5 * log(v) +
    (b * b) / (2 * lambda) - (m * m) / (2 * v) +
    logPhi_post - logZ_prior
  
  if (return_log) logI else exp(logI)
}

aggregate_by_distance <- function(K, z_i, z_j, A_ij, N_edge){
  d_vec <- abs(z_i - z_j)
  keep  <- d_vec > 0L
  bar_y <- numeric(K - 1)
  if (any(keep)) {
    # forward = lower block -> higher block
    forward <- z_i[keep] < z_j[keep]
    A_fwd   <- ifelse(forward, A_ij[keep], N_edge[keep] - A_ij[keep])
    S_vec   <- A_fwd - N_edge[keep]/2
    by_y    <- tapply(S_vec, d_vec[keep], sum)
    bar_y[as.integer(names(by_y))] <- by_y
  }
  list(bar_y = bar_y)
}


draw_omega_bar <- function(B_d, psi, zmax = 40) {
  if (length(B_d) != length(psi)) {
    stop(sprintf(paste0(
      "[draw_omega_bar] length(B_d)=%d implies K=%d, but length(psi)=%d. ",
      "Labels must be contiguous 1:K."
    ), length(B_d), length(B_d) + 1L, length(psi)), call. = FALSE)
  }
  len <- length(B_d)
  if (len == 0L) return(numeric(0))
  
  out <- numeric(len)
  keep <- (B_d > 0)
  
  if (any(keep)) {
    z <- psi[keep]
    if (anyNA(z) || any(!is.finite(z))) {
      stop("[draw_omega_bar] non-finite tilt in psi[keep]", call. = FALSE)
    }
    z <- pmin(pmax(z, -zmax), zmax)
    out[keep] <- BayesLogit::rpg(num = sum(keep), h = B_d[keep], z = z)
  }
  out
}
