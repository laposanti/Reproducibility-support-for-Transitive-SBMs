# WST-specific ordered SBM helpers for birth moves, directional scores, and psi updates.
#
rtruncnorm_pos <- function(n = 1L, mean = 0, sd = 1) {
  stopifnot(sd > 0)
  if (!is.finite(sd)) {
    # Degenerate case: infinite variance -> draw from prior mean
    return(rep(max(mean, 1e-6), n))
  }
  if (!is.finite(mean)) mean <- 0
  a <- (0 - mean) / sd
  # Guard inversion when pnorm(a) numerically hits 1 in the extreme right tail.
  upper <- 1 - 2 * .Machine$double.eps
  lower <- pnorm(a)
  if (!is.finite(lower)) lower <- 0
  if (lower >= upper) {
    return(rep(pmax(mean, 0) + 1e-6, n))
  }
  u <- runif(n, min = lower, max = upper)
  x <- mean + sd * qnorm(u)
  # Final guard: replace any non-finite or negative with a small positive fallback
  bad <- !is.finite(x) | x < 0
  if (any(bad)) x[bad] <- pmax(mean, 0) + 1e-6
  x
}

pack_state_wst <- function(z, kappa, psi) {
  u <- sort(unique(z))
  K_new <- length(u)
  
  # map old labels -> 1..K_new
  map <- integer(max(u))
  map[u] <- seq_len(K_new)
  z_new <- map[z]
  
  kappa_new <- kappa[u, u, drop = FALSE]
  psi_new   <- psi[u, u, drop = FALSE]
  diag(psi_new) <- 0
  psi_new[lower.tri(psi_new)] <- t(psi_new)[lower.tri(psi_new)]
  
  list(z = as.integer(z_new), kappa = kappa_new, psi = psi_new, K = K_new)
}


.log2cosh_half <- function(x) {
  u <- abs(x) / 2
  u + log1p(exp(-2*u))
}

# =============================================================================
# PG-AUGMENTED z-UPDATE HELPERS (shared by WST and SST)
# =============================================================================

# Draw edge-level PG latents omega ~ PG(N_edge, |phi|) for the full edge list.
# Called ONCE per z-sweep; the resulting vector is passed into node-level updates.
# mode = "pair" for WST (psi_obj is KxK matrix), "distance" for SST (psi_obj is vector).
draw_edge_omega <- function(z, psi_obj, i_idx, j_idx, N_edge, mode = "pair") {
  n_edges <- length(N_edge)
  phi_abs <- numeric(n_edges)
  zi <- z[i_idx]; zj <- z[j_idx]

  if (mode == "pair") {
    psi_mat <- psi_obj
    for (e in seq_len(n_edges)) {
      p <- min(zi[e], zj[e]); q <- max(zi[e], zj[e])
      if (p == q) { phi_abs[e] <- 0 } else { phi_abs[e] <- psi_mat[p, q] }
    }
  } else {
    psi_vec <- psi_obj
    d <- abs(zi - zj)
    for (e in seq_len(n_edges)) {
      if (d[e] == 0L) { phi_abs[e] <- 0 } else { phi_abs[e] <- psi_vec[d[e]] }
    }
  }

  omega <- numeric(n_edges)
  pos <- N_edge > 0
  if (any(pos)) {
    tilt <- pmin(pmax(phi_abs[pos], 0), 40)
    omega[pos] <- BayesLogit::rpg(sum(pos), h = N_edge[pos], z = tilt)
    bad <- !is.finite(omega[pos]) | omega[pos] <= 0
    if (any(bad)) {
      idx_bad <- which(pos)[bad]
      omega[idx_bad] <- N_edge[idx_bad] / 4
    }
  }
  omega
}

# Aggregate edge-level omega by block for node i.
aggregate_omega_by_block <- function(i, z, omega_edge, i_idx, j_idx,
                                     edge_by_node, K) {
  Omega_blk <- numeric(K)
  edges_i <- edge_by_node[[i]]
  if (length(edges_i) == 0L) return(Omega_blk)

  for (e in edges_i) {
    j <- if (i == i_idx[e]) j_idx[e] else i_idx[e]
    ell <- z[j]
    if (ell >= 1L && ell <= K) {
      Omega_blk[ell] <- Omega_blk[ell] + omega_edge[e]
    }
  }
  Omega_blk
}

# =============================================================================
# PG-AUGMENTED WST DIRECTIONAL SCORES
# =============================================================================

# Existing-block score using quadratic PG kernel: Y*psi - 0.5*Omega*psi^2
dir_pg_WST_existing <- function(k, c_plus, N_tot, Omega_blk, psi_mat) {
  K <- length(N_tot)
  if (K <= 1L) return(0)
  val <- 0
  for (ell in seq_len(K)) {
    if (ell == k) next
    n <- N_tot[ell]
    if (n == 0) next
    th <- psi_mat[min(k, ell), max(k, ell)]
    A_fwd <- if (ell > k) c_plus[ell] else (n - c_plus[ell])
    Y <- A_fwd - 0.5 * n
    val <- val + Y * th - 0.5 * Omega_blk[ell] * th^2
  }
  val
}

# Closed-form integral for a single new-block pair (WST) under PG kernel.
# From supplement eqs. (scalar-int-unc-log) + (Isl-WST-ratio).
log_int_pg_closed <- function(Y_star, Omega_star, mu0, sigma0) {
  sig2_0 <- sigma0^2
  A_prec <- Omega_star + 1 / sig2_0
  B_loc  <- Y_star + mu0 / sig2_0
  m <- B_loc / A_prec
  s <- 1 / sqrt(A_prec)

  log_I_unc <- -0.5 * log(1 + sig2_0 * Omega_star) -
    mu0^2 / (2 * sig2_0) + B_loc^2 / (2 * A_prec)

  log_trunc <- pnorm(m / s, log.p = TRUE) - pnorm(mu0 / sigma0, log.p = TRUE)

  log_I_unc + log_trunc
}

# New-block score: closed-form integral over each new psi_{new, ell}.
dir_pg_WST_new <- function(c_plus, N_tot, Omega_blk, mu0, sigma0,
                           r_set = NULL) {
  K_old <- length(N_tot)
  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)

  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(lp)

  for (idx in seq_along(r_set)) {
    r <- r_set[idx]
    val <- 0
    for (ell in seq_len(K_old)) {
      n <- N_tot[ell]
      if (n == 0) next
      pos_ell <- if (ell < r) ell else ell + 1L
      A_fwd <- if (pos_ell > r) c_plus[ell] else (n - c_plus[ell])
      Y_star <- A_fwd - 0.5 * n
      val <- val + log_int_pg_closed(Y_star, Omega_blk[ell], mu0, sigma0)
    }
    lp[idx] <- val
  }
  lp
}

# =============================================================================
# PG-AUGMENTED z-UPDATE WRAPPER — WST
# =============================================================================
wst_update_i_with_birth_LOO_pg <- function(
    i,
    A, z, eta, kappa, psi,
    Rkl, Tkl,
    i_idx, j_idx, N_edge, edge_by_node,
    a_kappa, b_kappa,
    gamma_gn,
    mu0, sigma0, sig2_0,
    omega_edge,
    r_set = NULL,
    slot_radius = NULL,
    partition_prior = "OCRP",
    theta_ocrp = 1.0
) {
  K_full <- nrow(kappa)
  oldk <- z[i]
  eta_i <- eta[i]

  # ---- LOO occupancy and packing map ----
  n_minus_full <- tabulate(z[-i], nbins = K_full)
  keep_full <- sort(which(n_minus_full > 0L))
  K_minus <- length(keep_full)
  if (K_minus < 1L) stop("LOO state has K_minus < 1.", call. = FALSE)

  map <- integer(K_full)
  map[keep_full] <- seq_len(K_minus)

  z_packed <- map[z]
  oldk_pos <- match(oldk, keep_full)

  # ---- pack kappa, psi ----
  kappa_minus <- kappa[keep_full, keep_full, drop = FALSE]
  psi_minus   <- psi[keep_full, keep_full, drop = FALSE]
  diag(psi_minus) <- 0
  psi_minus[lower.tri(psi_minus)] <- t(psi_minus)[lower.tri(psi_minus)]

  # ---- node vs block counts ----
  C_i_full <- counts_by_block_exact_cpp(
    i = i, A = A, z = as.integer(z),
    i_idx = as.integer(i_idx), j_idx = as.integer(j_idx),
    N_edge = as.numeric(N_edge),
    edge_by_node = edge_by_node, K = K_full
  )
  c_plus <- as.numeric(C_i_full$c_plus)[keep_full]
  N_tot  <- as.numeric(C_i_full$N_tot)[keep_full]

  # ---- aggregate omega by block ----
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
    dir_pg_WST_existing(kc, c_plus, N_tot, Omega_blk, psi_minus)
  }, numeric(1))

  lp_dir_new <- rep(-Inf, K_minus + 1L)
  lp_dir_new[r_set] <- dir_pg_WST_new(
    c_plus, N_tot, Omega_blk, mu0, sigma0, r_set
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

  if (!is.finite(logW_new) && !is.finite(logW_join)) {
    stop("Both existing and new assignment weights are -Inf/NaN.", call. = FALSE)
  }

  choose_new <- (log(runif(1)) < (logW_new - lse(c(logW_new, logW_join))))

  z_new <- integer(length(z))
  z_new[-i] <- z_packed[-i]

  if (!choose_new) {
    k_star <- .softsample(total_exist)
    z_new[i] <- k_star
    return(list(z = as.integer(z_new), kappa = kappa_minus,
                psi = psi_minus, K = K_minus,
                changed = TRUE, move = list(type = "existing", k = k_star)))
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

  init_psi <- numeric(K_minus)
  for (ell in seq_len(K_minus)) {
    n <- N_tot[ell]
    if (n == 0) {
      init_psi[ell] <- rtruncnorm_pos(1, mean = mu0, sd = sqrt(sig2_0))
      next
    }
    pos_ell <- if (ell < r_star) ell else ell + 1L
    A_fwd <- if (pos_ell > r_star) c_plus[ell] else (n - c_plus[ell])
    y <- A_fwd - 0.5 * n
    omega_init <- Omega_blk[ell]
    if (!is.finite(omega_init) || omega_init <= 0) omega_init <- n / 4
    v <- 1 / (omega_init + 1 / sig2_0)
    m <- v * (y + mu0 / sig2_0)
    raw <- rtruncnorm_pos(1, mean = m, sd = sqrt(v))
    init_psi[ell] <- pmin(pmax(raw, 0), 20)
  }

  grown <- grow_wst_params(
    kappa = kappa_minus, psi_mat = psi_minus, r = r_star,
    init_kappa_rowcol = init_kap, init_psi_rowcol = init_psi,
    a_kappa = a_kappa, b_kappa = b_kappa
  )

  list(z = as.integer(z_shift), kappa = grown$kappa, psi = grown$psi,
       K = nrow(grown$kappa),
       changed = TRUE, move = list(type = "new", r = r_star))
}

# ---------------------------
# FAST: existing-block WST
# ---------------------------
# psi_mat is KxK symmetric, diag 0, with psi_mat[k,ell] = psi_{min(k,ell),max(k,ell)} >= 0
# c_plus[ell] = wins of node i vs block ell (A[i,j] summed over j in ell)
# N_tot[ell]  = total matches node i vs block ell (N[i,j] summed over j in ell)
#
# Kernel version excludes within-block (ell==k). If include_within=TRUE, add -N_tot[k]*log(2).
dir_exact_WST_existing_counts <- function(k, c_plus, N_tot, psi_mat, include_within = FALSE) {
  K <- length(N_tot)
  stopifnot(length(c_plus) == K, is.matrix(psi_mat), nrow(psi_mat) == K, ncol(psi_mat) == K)
  
  if (K <= 1L) return(0)
  
  val <- 0
  for (ell in seq_len(K)) {
    if (ell == k) next
    n <- N_tot[ell]
    if (n == 0) next
    
    th <- psi_mat[min(k, ell), max(k, ell)]
    
    # forward = lower label -> higher label
    A_fwd <- if (ell > k) c_plus[ell] else (n - c_plus[ell])
    
    # exact binomial log-lik written in canonical "cosh" form
    val <- val + (A_fwd - 0.5 * n) * th - n * .log2cosh_half(th)
  }
  
  if (include_within) {
    val <- val - N_tot[k] * log(2)
  }
  
  val
}


# new-slot: integrates ψ_{new,ell} under N^+(mu0, sigma0^2) per ell
log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
logsig   <- function(x) -log1pexp(-x)
log1m_sig<- function(x) -log1pexp(x)

# --- generic integral: Binomial-logistic likelihood + truncated Normal prior ---
# Integrates ψ ~ N(mu0, sigma0^2) truncated to ψ >= 0.
# Likelihood is Binom(n, logistic(+ψ)) evaluated at forward-count A (lower->higher wins).
#
# Uses a Laplace approximation: the integrand is log-concave so the approximation
# is very accurate for block-level counts. This avoids the O(quadrature_points) R
# closure calls that made the original integrate()-based version the dominant cost
# (55 % of WST wall-time on dense datasets such as citations_data).
#
# Interior mode (ψ* > 0):  log∫ ≈ logf(ψ*) + ½ log(2π/H)
# Boundary mode (ψ* ≈ 0):  log∫ ≈ logf(0)  + ½ log(π/(2H))
# where H = n σ(ψ*)(1−σ(ψ*)) + 1/σ₀²  (exact Hessian magnitude of log-integrand).
log_int_binom_truncnorm <- function(A, n, mu0, sigma0, rel.tol = 1e-10) {
  stopifnot(n >= 0, A >= 0, A <= n, sigma0 > 0)
  if (n == 0) return(0)

  # normalising constant for truncation ψ >= 0
  logZ <- pnorm(mu0 / sigma0, log.p = TRUE)  # log Φ(mu/sigma)

  logf <- function(psi) {
    A * logsig(psi) + (n - A) * log1m_sig(psi) +
      dnorm(psi, mean = mu0, sd = sigma0, log = TRUE) - logZ
  }

  # Find mode on [0, 40] (same interval as before)
  ub <- 40
  opt      <- optimize(function(x) -logf(x), interval = c(0, ub))
  psi_star <- opt$minimum

  # Exact Hessian magnitude of log-integrand at psi_star
  p_star <- 1 / (1 + exp(-psi_star))
  H      <- n * p_star * (1 - p_star) + 1 / sigma0^2

  logf_star <- logf(psi_star)

  # Laplace approximation: half-Gaussian correction when mode hits boundary
  if (psi_star < 1e-6) {
    logf_star + 0.5 * log(pi / (2 * H))
  } else {
    logf_star + 0.5 * log(2 * pi / H)
  }
}
# ---------------------------
# FAST: new-slot WST (birth)
# ---------------------------
# Inserting a new block at slot r creates NEW psi_{new,ell} for each existing block ell,
# each integrated out under ψ ~ N^+(mu0, sigma0^2), independent across ell.
# Hence this term depends on (c_plus, N_tot) and the slot mapping, not on existing psi_mat.
#
# Returns a vector over r = 1..K+1 unless r_set is provided.

dir_exact_WST_new_vec_counts <- function(c_plus, N_tot, mu0 = 0, sigma0 = 1,
                                         r_set = NULL, rel.tol = 1e-10) {
  K_old <- length(N_tot)
  stopifnot(length(c_plus) == K_old)
  
  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)
  r_set <- as.integer(r_set)
  
  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(lp)
  
  for (idx in seq_along(r_set)) {
    r <- r_set[idx]
    val <- 0
    
    for (ell in seq_len(K_old)) {
      n <- N_tot[ell]
      if (n == 0) next
      
      # after insertion at r, old block ell shifts if ell >= r
      pos_ell <- if (ell < r) ell else ell + 1L
      
      # forward = lower label -> higher label (lower beats higher)
      A_fwd <- if (pos_ell > r) c_plus[ell] else (n - c_plus[ell])
      
      val <- val + log_int_binom_truncnorm(A = A_fwd, n = n, mu0 = mu0, sigma0 = sigma0,
                                           rel.tol = rel.tol)
    }
    
    lp[idx] <- val
  }
  
  lp
}



update_psi_wst_pair <- function(K, z, i_idx, j_idx, A_ij, N_edge,
                                psi_curr, mu0, sig2_0) {
  stopifnot(is.matrix(psi_curr), nrow(psi_curr) == K, ncol(psi_curr) == K)

  # --- guard: repair any NA/Inf/negative entries before using as rpg tilt ---
  bad_idx <- !is.finite(psi_curr) | psi_curr < 0
  if (any(bad_idx)) {
    warning(paste0("[update_psi_wst_pair] psi_curr has ", sum(bad_idx),
                   " non-finite/negative entries - resetting to 0."))
    psi_curr[bad_idx] <- 0
  }

  psi_new <- psi_curr
  diag(psi_new) <- 0
  if (K < 2L) return(psi_new)   # <<< critical guard
  
  B <- matrix(0, K, K)
  Y <- matrix(0, K, K)
  
  zi <- z[i_idx]; zj <- z[j_idx]
  p <- pmin(zi, zj)
  q <- pmax(zi, zj)
  
  # extra invariant check: labels must be within 1..K
  if (anyNA(p) || anyNA(q) || max(q) > K || min(p) < 1L) {
    stop("update_psi_wst_pair(): z labels out of [1..K]. Pack state before psi update.", call. = FALSE)
  }
  
  A_fwd <- ifelse(zi < zj, A_ij, N_edge - A_ij)
  
  for (e in seq_along(N_edge)) {
    if (p[e] == q[e]) next
    B[p[e], q[e]] <- B[p[e], q[e]] + N_edge[e]
    Y[p[e], q[e]] <- Y[p[e], q[e]] + (A_fwd[e] - 0.5 * N_edge[e])
  }
  
  for (k in seq_len(K - 1L)) {
    for (l in (k + 1L):K) {
      if (B[k, l] > 0) {
        z_tilt <- pmax(psi_curr[k, l], 0)
        if (!is.finite(z_tilt)) z_tilt <- 0  # extra safety
        omega <- BayesLogit::rpg(1, h = B[k, l], z = z_tilt)
        if (!is.finite(omega) || omega <= 0) omega <- B[k, l] / 4  # PG mean at z=0 is h/4
        v <- 1 / (omega + 1 / sig2_0)
        m <- v * (Y[k, l] + mu0 / sig2_0)
        psi_new[k, l] <- rtruncnorm_pos(1, mean = m, sd = sqrt(v))
      } else {
        psi_new[k, l] <- rtruncnorm_pos(1, mean = mu0, sd = sqrt(sig2_0))
      }
      psi_new[l, k] <-  psi_new[k, l]
    }
  }
  
  psi_new
}


grow_wst_params <- function(kappa, psi_mat, r,
                            init_kappa_rowcol,
                            init_psi_rowcol,
                            a_kappa = 1, b_kappa = 1) {
  K_old <- nrow(kappa)
  stopifnot(ncol(kappa) == K_old,
            is.matrix(psi_mat), nrow(psi_mat) == K_old, ncol(psi_mat) == K_old,
            length(init_kappa_rowcol) == K_old,
            length(init_psi_rowcol) == K_old)
  
  K_new <- K_old + 1L
  map_new <- function(o) if (o < r) o else (o + 1L)
  
  # κ expand
  kap2 <- matrix(0, K_new, K_new)
  for (a in seq_len(K_old)) for (b in seq_len(K_old)) {
    kap2[map_new(a), map_new(b)] <- kappa[a, b]
  }
  for (ell in seq_len(K_old)) {
    j <- map_new(ell)
    kap2[r, j] <- kap2[j, r] <- init_kappa_rowcol[ell]
  }
  kap2[r, r] <- rgamma(1, shape = a_kappa, rate = b_kappa)
  
  # ψ expand
  psi2 <- matrix(0, K_new, K_new)
  for (a in seq_len(K_old)) for (b in seq_len(K_old)) {
    psi2[map_new(a), map_new(b)] <- psi_mat[a, b]
  }
  for (ell in seq_len(K_old)) {
    j <- map_new(ell)
    # ψ stored as unordered pair min/max in matrix; symmetric anyway
    psi2[r, j] <- psi2[j, r] <- init_psi_rowcol[ell]
  }
  diag(psi2) <- 0
  
  list(kappa = kap2, psi = psi2)
}


wst_update_i_with_birth_LOO <- function(
    i,
    A, z, eta, kappa, psi,
    Rkl, Tkl,
    i_idx, j_idx, N_edge, edge_by_node,
    a_kappa, b_kappa,
    gamma_gn,
    mu0, sigma0, sig2_0,
    r_set = NULL,
    slot_radius = NULL
) {
  K_full <- nrow(kappa)
  oldk <- z[i]
  eta_i <- eta[i]
  
  # ---- LOO occupancy and packing map ----
  n_minus_full <- tabulate(z[-i], nbins = K_full)
  keep_full <- which(n_minus_full > 0L)
  keep_full <- sort(keep_full)
  K_minus <- length(keep_full)
  if (K_minus < 1L) stop("LOO state has K_minus < 1.", call. = FALSE)
  
  map <- integer(K_full)
  map[keep_full] <- seq_len(K_minus)
  
  z_packed <- map[z]
  oldk_pos <- match(oldk, keep_full)  # NA if oldk was singleton
  
  # ---- pack κ, ψ ----
  kappa_minus <- kappa[keep_full, keep_full, drop = FALSE]
  psi_minus   <- psi[keep_full, keep_full, drop = FALSE]
  diag(psi_minus) <- 0
  psi_minus[lower.tri(psi_minus)] <- t(psi_minus)[lower.tri(psi_minus)]
  
  # ---- exact node-vs-block directed counts, then subset ----
  C_i_full <- counts_by_block_exact_cpp(
    i = i, A = A, z = as.integer(z),
    i_idx = as.integer(i_idx), j_idx = as.integer(j_idx),
    N_edge = as.numeric(N_edge),
    edge_by_node = edge_by_node,
    K = K_full
  )
  c_plus_full <- as.numeric(C_i_full$c_plus)
  N_tot_full  <- as.numeric(C_i_full$N_tot)
  
  c_plus <- c_plus_full[keep_full]
  N_tot  <- N_tot_full [keep_full]
  
  # ---- exposure sums excluding i (FULL indexing), then subset ----
  E_excl_full <- E_by_block_excluding_i(i, z, eta, K_full)
  E_excl <- as.numeric(E_excl_full[keep_full])
  
  # increments when assigning i to some block
  r_add <- N_tot
  t_add <- eta_i * E_excl
  
  # ---- build LOO base R/T among others ----
  R_minus <- Rkl[keep_full, keep_full, drop = FALSE]
  T_minus <- Tkl[keep_full, keep_full, drop = FALSE]
  
  # subtract i’s contribution from its old block if it survives
  if (!is.na(oldk_pos)) {
    for (ell_pos in seq_len(K_minus)) {
      subR <- N_tot[ell_pos]
      subT <- eta_i * E_excl[ell_pos]
      p <- min(oldk_pos, ell_pos); q <- max(oldk_pos, ell_pos)
      R_minus[p, q] <- R_minus[p, q] - subR
      T_minus[p, q] <- T_minus[p, q] - subT
      R_minus[q, p] <- R_minus[p, q]
      T_minus[q, p] <- T_minus[p, q]
      if (R_minus[p, q] < -1e-8 || T_minus[p, q] < -1e-8) stop("Negative LOO totals.", call. = FALSE)
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
  } else {
    r_set <- as.integer(r_set)
  }
  
  # ---- (A) directional ----
  lp_dir_exist <- vapply(seq_len(K_minus), function(kc) {
    dir_exact_WST_existing_counts(k = kc, c_plus = c_plus, N_tot = N_tot, psi_mat = psi_minus, include_within = TRUE)
  }, numeric(1))
  
  lp_dir_new <- rep(-Inf, K_minus + 1L)
  lp_dir_new[r_set] <- dir_exact_WST_new_vec_counts(
    c_plus = c_plus, N_tot = N_tot,
    mu0 = mu0, sigma0 = sigma0,
    r_set = r_set
  )
  
  # ---- (B) collapsed κ ----
  lp_kappa_exist <- numeric(K_minus)
  for (kc in seq_len(K_minus)) {
    acc <- 0
    for (ell in seq_len(K_minus)) {
      if (r_add[ell] == 0 && t_add[ell] == 0) next
      p <- min(kc, ell); q <- max(kc, ell)
      acc <- acc + (gp_marginal(R_minus[p, q] + r_add[ell], T_minus[p, q] + t_add[ell], a_kappa, b_kappa) -
                      gp_marginal(R_minus[p, q],         T_minus[p, q],         a_kappa, b_kappa))
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
  
  if (!is.finite(logW_new) && !is.finite(logW_join)) {
    # nothing is possible, this indicates upstream bug or all weights nuked
    stop("Both existing and new assignment weights are -Inf/NaN. Inspect totals.", call. = FALSE)
  }
  
  
  choose_new <- (log(runif(1)) < (logW_new - lse(c(logW_new, logW_join))))
  
  
  # ---- apply move ----
  z_new <- integer(length(z))
  z_new[-i] <- z_packed[-i]
  
  if (!choose_new) {
    k_star <- .softsample(total_exist)
    z_new[i] <- k_star
    return(list(z = as.integer(z_new), kappa = kappa_minus, psi = psi_minus, K = K_minus,
                changed = TRUE, move = list(type = "existing", k = k_star)))
  }
  
  # new block
  r_star <- r_set[ .softsample(total_new[r_set] - log(length(r_set))) ]
  
  z_shift <- z_new
  z_shift[-i] <- ifelse(z_new[-i] < r_star, z_new[-i], z_new[-i] + 1L)
  z_shift[i] <- r_star
  
  # κ init (same as SST)
  init_kap <- rgamma(K_minus, shape = a_kappa + r_add, rate = b_kappa + t_add)
  
  # ψ init for new block vs each existing block using 1-step PG Gaussian approx on aggregated Binomial
  init_psi <- numeric(K_minus)
  for (ell in seq_len(K_minus)) {
    n <- N_tot[ell]
    if (n == 0) {
      init_psi[ell] <- rtruncnorm_pos(1, mean = mu0, sd = sqrt(sig2_0))
      next
    }
    pos_ell <- if (ell < r_star) ell else ell + 1L
    A_fwd <- if (pos_ell > r_star) c_plus[ell] else (n - c_plus[ell])
    y <- A_fwd - 0.5 * n
    
    psi0 <- max(0, mu0)
    omega <- BayesLogit::rpg(1, h = n, z = psi0)
    if (!is.finite(omega) || omega <= 0) omega <- n / 4  # PG mean at z=0 is h/4
    v <- 1 / (omega + 1 / sig2_0)
    m <- v * (y + mu0 / sig2_0)

    raw <- rtruncnorm_pos(1, mean = m, sd = sqrt(v))
    # Clamp to a safe range to prevent extreme psi values entering the chain
    init_psi[ell] <- pmin(pmax(raw, 0), 20)
  }
  
  grown <- grow_wst_params(
    kappa = kappa_minus, psi_mat = psi_minus, r = r_star,
    init_kappa_rowcol = init_kap,
    init_psi_rowcol   = init_psi,
    a_kappa = a_kappa, b_kappa = b_kappa
  )
  
  list(z = as.integer(z_shift), kappa = grown$kappa, psi = grown$psi, K = nrow(grown$kappa),
       changed = TRUE, move = list(type = "new", r = r_star))
}
