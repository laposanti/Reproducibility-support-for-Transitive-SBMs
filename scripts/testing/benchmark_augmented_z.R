#!/usr/bin/env Rscript
# =============================================================================
# Benchmark: Collapsed vs PG-Augmented z-update for WST and SST
# =============================================================================
# Compares:
#   (a) collapsed z-update  – exact Binom-logistic kernel + numerical integration
#   (b) PG-augmented z-update – quadratic PG kernel + closed-form new-block integral
#
# Metrics:
#   - Wall time per z-sweep
#   - Trace of K over iterations
#   - Effective sample size (ESS) of K and log-lik
#   - Partition accuracy (ARI vs truth)
# =============================================================================

cat("=== Benchmark: collapsed vs PG-augmented z-update ===\n")
t0_script <- proc.time()[3]

suppressPackageStartupMessages({
  library(Matrix)
  library(BayesLogit)
  library(mcclust)    # for arandi
})

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")

# Source helpers
source("helper_folder/helper.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("helper_folder/mixing_moves.R")
Rcpp::sourceCpp("core/counts_by_block_exact_cpp.cpp")
Rcpp::sourceCpp("core/block_totals_for_poisson_cpp.cpp")

# Source main sampler for variable-K helpers
source("core/my_best_try_so_far.R")

cat("  Sources loaded.\n")

# =============================================================================
# (1) PG-AUGMENTED DIRECTIONAL FUNCTIONS — WST
# =============================================================================

# --- Draw edge-level omega for the whole graph (once per sweep) ---------------
draw_edge_omega <- function(z, psi_obj, i_idx, j_idx, N_edge, mode = "pair") {
  n_edges <- length(N_edge)
  phi_abs <- numeric(n_edges)

  zi <- z[i_idx]; zj <- z[j_idx]

  if (mode == "pair") {
    # WST: phi = psi_mat[min(zi,zj), max(zi,zj)]
    psi_mat <- psi_obj
    for (e in seq_len(n_edges)) {
      p <- min(zi[e], zj[e]); q <- max(zi[e], zj[e])
      if (p == q) { phi_abs[e] <- 0 } else { phi_abs[e] <- psi_mat[p, q] }
    }
  } else {
    # SST: phi = psi_vec[|zi-zj|]
    psi_vec <- psi_obj
    d <- abs(zi - zj)
    for (e in seq_len(n_edges)) {
      if (d[e] == 0L) { phi_abs[e] <- 0 } else { phi_abs[e] <- psi_vec[d[e]] }
    }
  }

  # draw omega ~ PG(N_edge, |phi|)
  omega <- numeric(n_edges)
  pos <- N_edge > 0
  if (any(pos)) {
    tilt <- pmin(pmax(phi_abs[pos], 0), 40)
    omega[pos] <- BayesLogit::rpg(sum(pos), h = N_edge[pos], z = tilt)
    # safety: replace non-finite with PG mean at z=0
    bad <- !is.finite(omega[pos]) | omega[pos] <= 0
    if (any(bad)) {
      idx_bad <- which(pos)[bad]
      omega[idx_bad] <- N_edge[idx_bad] / 4
    }
  }
  omega
}

# --- Aggregate omega by block for a single node i ----------------------------
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

# --- PG-augmented WST: existing-block directional score -----------------------
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
  # No within-block correction needed (PG kernel at psi=0 gives 0)
  val
}

# --- Log closed-form integral for a single new-block pair (WST) ---------------
# From supplement eq. (scalar-int-unc-log) + (Isl-WST-ratio):
#   log I = -1/2 log(1 + sigma0^2 * Omega*) - mu0^2/(2 sigma0^2)
#           + B^2/(2A) + log Phi(m/s) - log Phi(mu0/sigma0)
# where A = Omega* + 1/sigma0^2, B = Y* + mu0/sigma0^2, m = B/A, s = 1/sqrt(A)
log_int_pg_closed <- function(Y_star, Omega_star, mu0, sigma0) {
  sig2_0 <- sigma0^2
  A_prec <- Omega_star + 1 / sig2_0
  B_loc  <- Y_star + mu0 / sig2_0
  m <- B_loc / A_prec
  s <- 1 / sqrt(A_prec)

  log_I_unc <- -0.5 * log(1 + sig2_0 * Omega_star) -
    mu0^2 / (2 * sig2_0) + B_loc^2 / (2 * A_prec)

  # truncation ratio: Phi(m/s) / Phi(mu0/sigma0)
  log_trunc <- pnorm(m / s, log.p = TRUE) - pnorm(mu0 / sigma0, log.p = TRUE)

  log_I_unc + log_trunc
}

# --- PG-augmented WST: new-block directional score (closed form) --------------
dir_pg_WST_new <- function(c_plus, N_tot, Omega_blk, mu0, sigma0,
                           r_set = NULL) {
  K_old <- length(N_tot)
  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)

  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(lp)
  if (K_old == 1L) { lp[] <- 0; return(lp) }

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
# (2) PG-AUGMENTED z-UPDATE WRAPPER — WST
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
    slot_radius = NULL
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

  # ---- node vs block counts (FULL indexing), then subset ----
  C_i_full <- counts_by_block_exact_cpp(
    i = i, A = A, z = as.integer(z),
    i_idx = as.integer(i_idx), j_idx = as.integer(j_idx),
    N_edge = as.numeric(N_edge),
    edge_by_node = edge_by_node, K = K_full
  )
  c_plus <- as.numeric(C_i_full$c_plus)[keep_full]
  N_tot  <- as.numeric(C_i_full$N_tot)[keep_full]

  # ---- aggregate omega by block (FULL indexing), then subset ----
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
  r_star <- r_set[.softsample(total_new[r_set] - log(length(r_set)))]

  z_shift <- z_new
  z_shift[-i] <- ifelse(z_new[-i] < r_star, z_new[-i], z_new[-i] + 1L)
  z_shift[i] <- r_star

  init_kap <- rgamma(K_minus, shape = a_kappa + r_add, rate = b_kappa + t_add)

  # init psi via 1-step PG conditional (same as collapsed version)
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


# =============================================================================
# (3) PG-AUGMENTED DIRECTIONAL FUNCTIONS — SST
# =============================================================================

# --- PG-augmented SST: existing-block directional score -----------------------
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

# --- Closed-form integral for SST extreme distance (d = K_old) ----------------
# Integrates delta ~ N+(m_inc, tau0^2) with PG kernel:
#   exp{Y*(psi_{K-1}+delta) - 0.5*Omega*(psi_{K-1}+delta)^2}
# Substituting u = psi_{K-1}+delta, lower limit = psi_{K-1}:
#   A' = Omega + 1/tau0^2
#   B' = Y + (psi_{K-1} + m_inc)/tau0^2
#   m' = B'/A', s' = 1/sqrt(A')
#   log I = -0.5*log(tau0^2*A') + B'^2/(2A') - (psi_{K-1}+m_inc)^2/(2tau0^2)
#           + log Phi((m'-psi_{K-1})/s') - log Phi(m_inc/tau0)
log_int_pg_sst_extreme <- function(Y_star, Omega_star, psiKm1, tau0,
                                   m_inc = 0) {
  tau2_0 <- tau0^2
  mu_u   <- psiKm1 + m_inc  # prior mean of u = psi_{K-1}+delta
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

# --- PG-augmented SST: new-block directional score ----------------------------
dir_pg_SST_new <- function(c_plus, N_tot, Omega_blk, psi_vec,
                           r_set = NULL, tau0, m_inc = 0) {
  K_old <- length(N_tot)
  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)

  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(numeric(0))
  if (K_old == 1L) { lp[] <- 0; return(lp) }
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
        # known distance: use PG kernel with existing psi_d
        th <- psi_vec[d]
        val <- val + Y_d * th - 0.5 * Omega_blk[ell] * th^2
      } else if (d == K_old) {
        # new farthest distance: accumulate for closed-form integral
        A_K <- A_K + A_fwd
        n_K <- n_K + n
        Omega_K <- Omega_K + Omega_blk[ell]
      }
    }

    # integrate out the new extreme distance
    if (n_K > 0) {
      Y_K <- A_K - 0.5 * n_K
      psiKm1 <- psi_vec[K_old - 1L]
      val <- val + log_int_pg_sst_extreme(Y_K, Omega_K, psiKm1, tau0, m_inc)
    }

    lp[idx] <- val
  }
  lp
}


# =============================================================================
# (4) PG-AUGMENTED z-UPDATE WRAPPER — SST
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
    sig2_0 = 1
) {
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
  lp_dir_new[r_set] <- dir_pg_SST_new(
    c_plus, N_tot, Omega_blk, psi_minus, r_set, tau0
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

  list(z = as.integer(z_shift), kappa = grown$kappa,
       psi = as.numeric(grown$psi), K = nrow(grown$kappa),
       changed = TRUE, move = list(type = "new", r = r_star))
}


# =============================================================================
# (5) FULL SWEEP FUNCTIONS
# =============================================================================

# --- WST z-sweep: collapsed (current code) ---
sweep_z_wst_collapsed <- function(z, A, eta, kappa, psi,
                                  i_idx, j_idx, N_edge, edge_by_node,
                                  hyper, psi_hyper) {
  n <- length(z); K <- nrow(kappa)
  bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  slot_rad <- K + 1L

  for (i in sample.int(n)) {
    z_old_i <- z[i]; K_old <- K
    res <- wst_update_i_with_birth_LOO(
      i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
      Rkl = bt$Rkl, Tkl = bt$Tkl,
      i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
      gamma_gn = psi_hyper$gamma_gn,
      mu0 = hyper$mu0, sigma0 = hyper$sigma0,
      sig2_0 = hyper$sigma0^2,
      slot_radius = slot_rad
    )
    z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K
    if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
    }
  }
  packed <- pack_state_wst(z, kappa, psi)
  list(z = packed$z, kappa = packed$kappa, psi = packed$psi, K = packed$K)
}

# --- WST z-sweep: PG-augmented ---
sweep_z_wst_augmented <- function(z, A, eta, kappa, psi,
                                  i_idx, j_idx, N_edge, edge_by_node,
                                  hyper, psi_hyper) {
  n <- length(z); K <- nrow(kappa)

  # draw edge-level omega ONCE for the sweep
  omega_edge <- draw_edge_omega(z, psi, i_idx, j_idx, N_edge, mode = "pair")

  bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  slot_rad <- K + 1L

  for (i in sample.int(n)) {
    z_old_i <- z[i]; K_old <- K
    res <- wst_update_i_with_birth_LOO_pg(
      i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
      Rkl = bt$Rkl, Tkl = bt$Tkl,
      i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
      gamma_gn = psi_hyper$gamma_gn,
      mu0 = hyper$mu0, sigma0 = hyper$sigma0,
      sig2_0 = hyper$sigma0^2,
      omega_edge = omega_edge,
      slot_radius = slot_rad
    )
    z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K
    if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
    }
  }
  packed <- pack_state_wst(z, kappa, psi)
  list(z = packed$z, kappa = packed$kappa, psi = packed$psi, K = packed$K)
}

# --- WST z-sweep + mixing moves: collapsed ---
sweep_z_wst_collapsed_mix <- function(z, A, eta, kappa, psi,
                                      i_idx, j_idx, N_edge, edge_by_node,
                                      hyper, psi_hyper) {
  res <- sweep_z_wst_collapsed(z, A, eta, kappa, psi,
                               i_idx, j_idx, N_edge, edge_by_node,
                               hyper, psi_hyper)
  z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K

  # adjacent-block swap
  if (K >= 2L) {
    swap_res <- adjacent_block_swap_move(
      z = z, kappa = kappa, psi = psi, eta = eta,
      A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node, psi_mode = "pair"
    )
    z <- swap_res$z; kappa <- swap_res$kappa; psi <- swap_res$psi
  }

  # split-merge
  sm_res <- split_merge_move(
    z = z, kappa = kappa, psi = psi, eta = eta,
    A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
    edge_by_node = edge_by_node,
    a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
    gamma_gn = psi_hyper$gamma_gn, psi_mode = "pair",
    hyper_psi = list(mu0 = hyper$mu0, sigma0 = hyper$sigma0, tau0 = hyper$tau0),
    n_restricted_scans = 3L
  )
  z <- sm_res$z; kappa <- sm_res$kappa; psi <- sm_res$psi
  K <- nrow(kappa)

  packed <- pack_state_wst(z, kappa, psi)
  list(z = packed$z, kappa = packed$kappa, psi = packed$psi, K = packed$K)
}

# --- WST z-sweep + mixing moves: PG-augmented ---
sweep_z_wst_augmented_mix <- function(z, A, eta, kappa, psi,
                                      i_idx, j_idx, N_edge, edge_by_node,
                                      hyper, psi_hyper) {
  res <- sweep_z_wst_augmented(z, A, eta, kappa, psi,
                               i_idx, j_idx, N_edge, edge_by_node,
                               hyper, psi_hyper)
  z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K

  # adjacent-block swap
  if (K >= 2L) {
    swap_res <- adjacent_block_swap_move(
      z = z, kappa = kappa, psi = psi, eta = eta,
      A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node, psi_mode = "pair"
    )
    z <- swap_res$z; kappa <- swap_res$kappa; psi <- swap_res$psi
  }

  # split-merge
  sm_res <- split_merge_move(
    z = z, kappa = kappa, psi = psi, eta = eta,
    A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
    edge_by_node = edge_by_node,
    a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
    gamma_gn = psi_hyper$gamma_gn, psi_mode = "pair",
    hyper_psi = list(mu0 = hyper$mu0, sigma0 = hyper$sigma0, tau0 = hyper$tau0),
    n_restricted_scans = 3L
  )
  z <- sm_res$z; kappa <- sm_res$kappa; psi <- sm_res$psi
  K <- nrow(kappa)

  packed <- pack_state_wst(z, kappa, psi)
  list(z = packed$z, kappa = packed$kappa, psi = packed$psi, K = packed$K)
}

# --- SST z-sweep: collapsed (current code) ---
sweep_z_sst_collapsed <- function(z, A, eta, kappa, psi,
                                  i_idx, j_idx, N_edge, edge_by_node,
                                  hyper, psi_hyper) {
  n <- length(z); K <- nrow(kappa)
  bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  slot_rad <- K + 1L

  for (i in sample.int(n)) {
    z_old_i <- z[i]; K_old <- K
    res <- sst_update_i_with_birth_LOO(
      i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
      Rkl = bt$Rkl, Tkl = bt$Tkl,
      i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
      gamma_gn = psi_hyper$gamma_gn,
      slot_radius = slot_rad,
      tau0 = psi_hyper$tau0,
      mu0  = psi_hyper$mu0,
      sig2_0 = psi_hyper$sigma0
    )
    z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K
    if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
    }
  }
  K_out <- nrow(kappa)
  list(z = z, kappa = kappa, psi = psi, K = K_out)
}

# --- SST z-sweep: PG-augmented ---
sweep_z_sst_augmented <- function(z, A, eta, kappa, psi,
                                  i_idx, j_idx, N_edge, edge_by_node,
                                  hyper, psi_hyper) {
  n <- length(z); K <- nrow(kappa)

  omega_edge <- draw_edge_omega(z, psi, i_idx, j_idx, N_edge, mode = "distance")

  bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  slot_rad <- K + 1L

  for (i in sample.int(n)) {
    z_old_i <- z[i]; K_old <- K
    res <- sst_update_i_with_birth_LOO_pg(
      i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
      Rkl = bt$Rkl, Tkl = bt$Tkl,
      i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
      gamma_gn = psi_hyper$gamma_gn,
      omega_edge = omega_edge,
      slot_radius = slot_rad,
      tau0 = psi_hyper$tau0,
      mu0 = psi_hyper$mu0,
      sig2_0 = psi_hyper$sigma0
    )
    z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K
    if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
    }
  }
  K_out <- nrow(kappa)
  list(z = z, kappa = kappa, psi = psi, K = K_out)
}


# =============================================================================
# (6) COMMON PSI / KAPPA / ETA UPDATES (shared by both chains)
# =============================================================================

update_kappa_block <- function(z, eta, kappa, i_idx, j_idx, N_edge, hyper) {
  K <- nrow(kappa)
  p <- pmin(z[i_idx], z[j_idx]); q <- pmax(z[i_idx], z[j_idx])
  pf <- factor(p, levels = 1:K); qf <- factor(q, levels = 1:K)
  Rkl <- as.matrix(xtabs(N_edge ~ pf + qf, drop.unused.levels = FALSE))

  E_k   <- tapply(eta, factor(z, levels = 1:K), sum);   E_k[is.na(E_k)] <- 0
  eta2k <- tapply(eta^2, factor(z, levels = 1:K), sum); eta2k[is.na(eta2k)] <- 0
  Tkl <- outer(E_k, E_k, `*`)
  diag(Tkl) <- pmax((E_k^2 - eta2k) / 2, 0)

  kappa_new <- matrix(0, K, K)
  for (k in seq_len(K)) for (l in k:K) {
    val <- rgamma(1, shape = hyper$a_kappa + Rkl[k, l],
                  rate = hyper$b_kappa + Tkl[k, l])
    kappa_new[k, l] <- kappa_new[l, k] <- val
  }
  kappa_new
}

update_eta_block <- function(z, eta, kappa, A, hyper) {
  n <- length(z); K <- nrow(kappa)
  G_i <- Matrix::rowSums(A) + Matrix::colSums(A)
  E_k <- tapply(eta, factor(z, levels = 1:K), sum); E_k[is.na(E_k)] <- 0

  for (i in seq_len(n)) {
    k <- z[i]; E_k[k] <- E_k[k] - eta[i]
    rate_i <- hyper$b_eta + sum(kappa[k, ] * E_k)
    eta[i] <- rgamma(1, shape = hyper$a_eta + G_i[i],
                     rate = max(rate_i, 1e-10))
    E_k[k] <- E_k[k] + eta[i]
  }
  n_k <- as.integer(tabulate(z, nbins = K))
  for (k in seq_len(K)) {
    idx_k <- which(z == k); if (!length(idx_k)) next
    s_k <- sum(eta[idx_k])
    if (s_k > 0) eta[idx_k] <- n_k[k] * eta[idx_k] / s_k
    else         eta[idx_k] <- n_k[k] / length(idx_k)
  }
  eta
}

update_psi_wst_block <- function(z, psi, A, i_idx, j_idx, N_edge, hyper) {
  K <- nrow(psi)
  if (K < 2L) return(psi)
  bar_y <- aggregate_by_pair(K, z_i = z[i_idx], z_j = z[j_idx],
                             A_ij = A[cbind(i_idx, j_idx)], N_edge = N_edge)
  B_mat <- pair_totals(K, z_i = z[i_idx], z_j = z[j_idx], N_edge = N_edge)
  bar_omega <- draw_omega_pair(B_mat, psi_mat = psi)
  psi_new <- update_psi_pair(bar_y = bar_y, bar_omega = bar_omega,
                             mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2,
                             trunc = TRUE)
  psi_new[lower.tri(psi_new)] <- t(psi_new)[lower.tri(psi_new)]
  psi_new
}

update_psi_sst_block <- function(z, psi_vec, A, i_idx, j_idx, N_edge, hyper) {
  K <- max(z)
  D <- K - 1L
  if (D < 1L) return(numeric(0))
  agg <- aggregate_by_distance(K, z_i = z[i_idx], z_j = z[j_idx],
                               A_ij = A[cbind(i_idx, j_idx)], N_edge = N_edge)
  B_d <- distance_totals(K, z_i = z[i_idx], z_j = z[j_idx], N_edge = N_edge)
  bar_omega <- draw_omega_bar(B_d = B_d, psi = psi_vec)
  psi_new <- update_psi_sst(K = K, bar_y = agg$bar_y, bar_omega = bar_omega,
                            psi_curr = psi_vec,
                            mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2,
                            tau2_0 = hyper$tau0^2, n_inner_sweeps = 4L)
  psi_new
}


# =============================================================================
# (7) MINI-SAMPLER: runs n_iter sweeps with a chosen z-update strategy
# =============================================================================

run_chain <- function(A, z_init, kappa_init, psi_init, eta_init,
                      i_idx, j_idx, N_edge, edge_by_node,
                      hyper, psi_hyper,
                      mode = c("WST", "SST"),
                      method = c("collapsed", "augmented",
                                "collapsed_mix", "augmented_mix"),
                      n_iter = 500, burn = 100, seed = 42,
                      truth_z = NULL) {
  mode   <- match.arg(mode)
  method <- match.arg(method)
  set.seed(seed)

  z     <- z_init
  kappa <- kappa_init
  psi   <- psi_init
  eta   <- eta_init
  n     <- length(z)

  # pick z-sweep function
  z_sweep_fn <- switch(
    paste0(mode, "_", method),
    WST_collapsed      = sweep_z_wst_collapsed,
    WST_augmented      = sweep_z_wst_augmented,
    WST_collapsed_mix  = sweep_z_wst_collapsed_mix,
    WST_augmented_mix  = sweep_z_wst_augmented_mix,
    SST_collapsed      = sweep_z_sst_collapsed,
    SST_augmented      = sweep_z_sst_augmented
  )

  # storage
  K_trace   <- integer(n_iter)
  ari_trace <- numeric(n_iter)
  z_time    <- numeric(n_iter)
  psi_time  <- numeric(n_iter)

  for (it in seq_len(n_iter)) {
    # --- z-sweep (timed) ---
    t1 <- proc.time()[3]
    res <- z_sweep_fn(z, A, eta, kappa, psi,
                      i_idx, j_idx, N_edge, edge_by_node,
                      hyper, psi_hyper)
    z_time[it] <- proc.time()[3] - t1

    z     <- res$z
    kappa <- res$kappa
    psi   <- res$psi
    K     <- res$K

    # --- psi update ---
    t2 <- proc.time()[3]
    if (mode == "WST") {
      psi <- update_psi_wst_block(z, psi, A, i_idx, j_idx, N_edge, hyper)
    } else {
      psi <- update_psi_sst_block(z, psi, A, i_idx, j_idx, N_edge, hyper)
    }
    psi_time[it] <- proc.time()[3] - t2

    # --- kappa update ---
    kappa <- update_kappa_block(z, eta, kappa, i_idx, j_idx, N_edge, hyper)

    # --- eta update ---
    eta <- update_eta_block(z, eta, kappa, A, hyper)

    K_trace[it] <- K
    ari_trace[it] <- if (!is.null(truth_z)) mcclust::arandi(z, truth_z) else NA

    if (it %% 100 == 0)
      cat(sprintf("  [%s-%s] it=%d/%d K=%d z_sweep=%.3fs\n",
                  mode, method, it, n_iter, K, z_time[it]))
  }

  post_idx <- (burn + 1):n_iter

  list(
    K_trace   = K_trace,
    ari_trace = ari_trace,
    z_time    = z_time,
    psi_time  = psi_time,
    # post-burn summaries
    K_mean    = mean(K_trace[post_idx]),
    ari_mean  = mean(ari_trace[post_idx], na.rm = TRUE),
    z_time_mean = mean(z_time[post_idx]),
    z_time_sd   = sd(z_time[post_idx]),
    z_final   = z,
    kappa_final = kappa,
    psi_final  = psi,
    eta_final  = eta
  )
}


# =============================================================================
# (8) DATA GENERATION
# =============================================================================

logistic <- function(x) 1 / (1 + exp(-x))

generate_network <- function(n, K, psi_true, kappa_true, z_true, eta_true,
                             mode = "WST") {
  A <- Matrix(0, n, n, sparse = TRUE)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    lam <- eta_true[i] * eta_true[j] * kappa_true[z_true[i], z_true[j]]
    N <- rpois(1, lam)
    if (N > 0) {
      if (mode == "WST") {
        k <- min(z_true[i], z_true[j])
        l <- max(z_true[i], z_true[j])
        th <- if (k == l) 0 else psi_true[k, l]
      } else {
        d <- abs(z_true[i] - z_true[j])
        th <- if (d == 0) 0 else psi_true[d]
      }
      pr <- if (z_true[i] <= z_true[j]) logistic(th) else logistic(-th)
      fwd <- rbinom(1, N, pr)
      A[i, j] <- fwd
      A[j, i] <- N - fwd
    }
  }
  A
}


# =============================================================================
# (9) RUN BENCHMARK
# =============================================================================

set.seed(2026)

# ---- Generate test data (HARDER: weak signal, cold start) ----
n <- 100; K_true <- 4
z_true <- rep(1:K_true, each = n / K_true)
eta_true <- rgamma(n, 10, 10)
n_k <- tabulate(z_true, nbins = K_true)
for (k in seq_len(K_true)) {
  idx_k <- which(z_true == k)
  eta_true[idx_k] <- n_k[k] * eta_true[idx_k] / sum(eta_true[idx_k])
}

# --- WST truth (weak psi -> harder to detect) ---
psi_true_wst <- matrix(0, K_true, K_true)
psi_true_wst[1, 2] <- 0.8; psi_true_wst[1, 3] <- 1.5; psi_true_wst[1, 4] <- 2.2
psi_true_wst[2, 3] <- 0.6; psi_true_wst[2, 4] <- 1.3; psi_true_wst[3, 4] <- 0.7
psi_true_wst[lower.tri(psi_true_wst)] <- t(psi_true_wst)[lower.tri(psi_true_wst)]

kappa_true <- matrix(2, K_true, K_true)
diag(kappa_true) <- 4

cat("\n--- Generating WST network (n=", n, ", K=", K_true, ") ---\n")
A_wst <- generate_network(n, K_true, psi_true_wst, kappa_true, z_true, eta_true,
                          mode = "WST")
A_wst <- as.matrix(A_wst)  # dense for Rcpp compatibility
cat("  sum(A)=", sum(A_wst), " density=", round(sum(A_wst > 0) / (n * (n - 1)), 3), "\n")

# --- SST truth ---
psi_true_sst <- c(0.8, 1.5, 2.2)
cat("--- Generating SST network (n=", n, ", K=", K_true, ") ---\n")
A_sst <- generate_network(n, K_true, psi_true_sst, kappa_true, z_true, eta_true,
                          mode = "SST")
A_sst <- as.matrix(A_sst)  # dense for Rcpp compatibility
cat("  sum(A)=", sum(A_sst), " density=", round(sum(A_sst > 0) / (n * (n - 1)), 3), "\n")

# ---- Edge list ----
build_edges <- function(A, n) {
  A_sym <- A + t(A)
  idx   <- which(A_sym > 0 & upper.tri(A_sym), arr.ind = TRUE)
  i_idx <- idx[, 1]; j_idx <- idx[, 2]
  N_edge <- A_sym[idx]
  edge_by_node <- vector("list", n)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], e)
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], e)
  }
  list(i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
       edge_by_node = edge_by_node)
}

edges_wst <- build_edges(A_wst, n)
edges_sst <- build_edges(A_sst, n)

# ---- Hyperparameters ----
hyper <- list(
  a_kappa = 1, b_kappa = 1,
  a_eta = 2, b_eta = 2,
  mu0 = 0.5, sigma0 = 1.5,
  tau0 = 0.5
)
psi_hyper <- list(
  gamma_gn = 0.5,
  mu0 = hyper$mu0,
  sigma0 = hyper$sigma0,
  tau0 = hyper$tau0
)

# ---- Initial state: COLD START (K=1, random) ----
N_ITER <- 800
BURN   <- 200

cat("\n====== WST BENCHMARK (4 variants) ======\n")
cat("Running", N_ITER, "iterations (burn=", BURN, "), COLD start from K=1\n")

# Cold start: everyone in one block
z_init <- rep(1L, n)
kappa_init_1 <- matrix(rgamma(1, 2, 1), 1, 1)
psi_init_wst_1 <- matrix(0, 1, 1)
eta_init <- rep(1, n)

res_wst_collapsed <- run_chain(
  A_wst, z_init, kappa_init_1, psi_init_wst_1, eta_init,
  edges_wst$i_idx, edges_wst$j_idx, edges_wst$N_edge, edges_wst$edge_by_node,
  hyper, psi_hyper, mode = "WST", method = "collapsed",
  n_iter = N_ITER, burn = BURN, seed = 42, truth_z = z_true
)

res_wst_augmented <- run_chain(
  A_wst, z_init, kappa_init_1, psi_init_wst_1, eta_init,
  edges_wst$i_idx, edges_wst$j_idx, edges_wst$N_edge, edges_wst$edge_by_node,
  hyper, psi_hyper, mode = "WST", method = "augmented",
  n_iter = N_ITER, burn = BURN, seed = 42, truth_z = z_true
)

res_wst_collapsed_mix <- run_chain(
  A_wst, z_init, kappa_init_1, psi_init_wst_1, eta_init,
  edges_wst$i_idx, edges_wst$j_idx, edges_wst$N_edge, edges_wst$edge_by_node,
  hyper, psi_hyper, mode = "WST", method = "collapsed_mix",
  n_iter = N_ITER, burn = BURN, seed = 42, truth_z = z_true
)

res_wst_augmented_mix <- run_chain(
  A_wst, z_init, kappa_init_1, psi_init_wst_1, eta_init,
  edges_wst$i_idx, edges_wst$j_idx, edges_wst$N_edge, edges_wst$edge_by_node,
  hyper, psi_hyper, mode = "WST", method = "augmented_mix",
  n_iter = N_ITER, burn = BURN, seed = 42, truth_z = z_true
)


# =============================================================================
# (10) RESULTS
# =============================================================================

# ESS via batch means (simple)
ess_bm <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 10) return(NA)
  b <- floor(sqrt(n))
  a <- floor(n / b)
  batch_means <- sapply(seq_len(a), function(i) mean(x[((i-1)*b+1):(i*b)]))
  var_grand <- var(x)
  var_batch <- var(batch_means)
  if (var_batch <= 0 || var_grand <= 0) return(n)
  pmin(n, n * var_grand / (b * var_batch))
}

cat("\n")
cat("=" , rep("=", 76), "\n", sep = "")
cat("                    BENCHMARK RESULTS\n")
cat("=", rep("=", 76), "\n", sep = "")

results <- data.frame(
  regime  = rep("WST", 4),
  method  = c("collapsed", "augmented", "collapsed+mix", "augmented+mix"),
  K_mean  = c(res_wst_collapsed$K_mean, res_wst_augmented$K_mean,
              res_wst_collapsed_mix$K_mean, res_wst_augmented_mix$K_mean),
  ARI     = c(res_wst_collapsed$ari_mean, res_wst_augmented$ari_mean,
              res_wst_collapsed_mix$ari_mean, res_wst_augmented_mix$ari_mean),
  z_time_ms = c(res_wst_collapsed$z_time_mean, res_wst_augmented$z_time_mean,
                res_wst_collapsed_mix$z_time_mean, res_wst_augmented_mix$z_time_mean) * 1000,
  z_time_sd = c(res_wst_collapsed$z_time_sd, res_wst_augmented$z_time_sd,
                res_wst_collapsed_mix$z_time_sd, res_wst_augmented_mix$z_time_sd) * 1000,
  ESS_K   = c(ess_bm(res_wst_collapsed$K_trace[(BURN+1):N_ITER]),
              ess_bm(res_wst_augmented$K_trace[(BURN+1):N_ITER]),
              ess_bm(res_wst_collapsed_mix$K_trace[(BURN+1):N_ITER]),
              ess_bm(res_wst_augmented_mix$K_trace[(BURN+1):N_ITER])),
  stringsAsFactors = FALSE
)
results$speedup <- NA
results$speedup[2] <- results$z_time_ms[1] / results$z_time_ms[2]
results$speedup[4] <- results$z_time_ms[3] / results$z_time_ms[4]
results$ESS_per_s <- results$ESS_K / (results$z_time_ms / 1000 * (N_ITER - BURN))

print(results, digits = 3)

cat("\n--- WST z-sweep timing ---\n")
cat(sprintf("  Collapsed: %.1f +/- %.1f ms/sweep\n",
            res_wst_collapsed$z_time_mean * 1000,
            res_wst_collapsed$z_time_sd * 1000))
cat(sprintf("  Augmented:      %.1f +/- %.1f ms/sweep\n",
            res_wst_augmented$z_time_mean * 1000,
            res_wst_augmented$z_time_sd * 1000))
cat(sprintf("  Collapsed+mix:  %.1f +/- %.1f ms/sweep\n",
            res_wst_collapsed_mix$z_time_mean * 1000,
            res_wst_collapsed_mix$z_time_sd * 1000))
cat(sprintf("  Augmented+mix:  %.1f +/- %.1f ms/sweep\n",
            res_wst_augmented_mix$z_time_mean * 1000,
            res_wst_augmented_mix$z_time_sd * 1000))
cat(sprintf("  Speedup (no mix):   %.2fx\n",
            res_wst_collapsed$z_time_mean / res_wst_augmented$z_time_mean))
cat(sprintf("  Speedup (with mix): %.2fx\n",
            res_wst_collapsed_mix$z_time_mean / res_wst_augmented_mix$z_time_mean))

cat("\n--- K trace (post-burn) ---\n")
all_res <- list(res_wst_collapsed, res_wst_augmented,
                res_wst_collapsed_mix, res_wst_augmented_mix)
all_labels <- c("collapsed", "augmented", "collapsed+mix", "augmented+mix")
for (j in seq_along(all_res)) {
  r <- all_res[[j]]
  cat(sprintf("  WST %-16s: K_mean=%.2f  ESS(K)=%.0f  ARI=%.3f\n",
              all_labels[j], r$K_mean,
              ess_bm(r$K_trace[(BURN+1):N_ITER]),
              r$ari_mean))
}

# ---- Save diagnostic plots ----
out_dir <- "output/simulation/benchmark_augmented_z"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(out_dir, "K_traces.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(res_wst_collapsed$K_trace, type = "l", col = "steelblue",
     main = "WST collapsed: K trace", xlab = "iteration", ylab = "K")
abline(h = K_true, lty = 2, col = "red")

plot(res_wst_augmented$K_trace, type = "l", col = "darkorange",
     main = "WST augmented: K trace", xlab = "iteration", ylab = "K")
abline(h = K_true, lty = 2, col = "red")

plot(res_wst_collapsed_mix$K_trace, type = "l", col = "steelblue",
     main = "WST collapsed+mix: K trace", xlab = "iteration", ylab = "K")
abline(h = K_true, lty = 2, col = "red")

plot(res_wst_augmented_mix$K_trace, type = "l", col = "darkorange",
     main = "WST augmented+mix: K trace", xlab = "iteration", ylab = "K")
abline(h = K_true, lty = 2, col = "red")

dev.off()

pdf(file.path(out_dir, "z_sweep_times.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

boxplot(list(collapsed = res_wst_collapsed$z_time * 1000,
             augmented = res_wst_augmented$z_time * 1000),
        main = "WST (no mix): z-sweep time (ms)", ylab = "ms",
        col = c("steelblue", "darkorange"))

boxplot(list("collapsed+mix" = res_wst_collapsed_mix$z_time * 1000,
             "augmented+mix" = res_wst_augmented_mix$z_time * 1000),
        main = "WST (with mix): z-sweep time (ms)", ylab = "ms",
        col = c("steelblue", "darkorange"))

dev.off()

cat("\nPlots saved to:", out_dir, "\n")
cat("Total script time:", round(proc.time()[3] - t0_script, 1), "s\n")
