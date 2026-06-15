# ---- Utilities: Normal pdf/cdf wrappers (numerically stable) ---------------
.Norm_cdf <- function(x) pnorm(x)
.Norm_pdf <- function(x) dnorm(x)

# ---- WST: closed-form log I for one block-pair -----------------------------
# Prior: ψ ~ N^+(μ0, σ0^2) on (0, ∞).
# Kernel: exp{ S ψ - 0.5 * ω ψ^2 }.
logI_wst_pair_closed <- function(S, omega, mu0, sig2) {
  # handle trivial case
  if (omega < 0) stop("omega must be >= 0")
  if (sig2  <= 0) stop("sig2 must be > 0")
  sig0  <- sqrt(sig2)
  prec0 <- 1 / sig2
  prec  <- omega + prec0
  mstar <- (S + prec0 * mu0) / prec
  # lower limit 0 mapped to z0
  z0    <- - mstar * sqrt(prec)
  # final closed form (expectation wrt truncated prior)
  val <- -0.5 * log(prec) - log(sig0) +
    0.5 * (S + prec0 * mu0)^2 / prec - 0.5 * (mu0^2) * prec0 -
    log(.Norm_cdf(mu0 / sig0)) +
    log(.Norm_cdf(-z0))
  as.numeric(val)
}

# ---- WST: Monte Carlo check for one block-pair ------------------------------
# Draw ψ ~ N^+(μ0, σ0^2) and compute mean( exp{S ψ - 0.5 ω ψ^2} ).
logI_wst_pair_mc <- function(S, omega, mu0, sig2, R = 20000L, seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  psi <- truncnorm::rtruncnorm(R, a = 0, b = Inf, mean = mu0, sd = sqrt(sig2))
  val <- mean( exp(S * psi - 0.5 * omega * psi^2) )
  log(pmax(val, 1e-300))
}
# Closed-form log ∫_{δ>0} exp{ (yK - ωK ψ_{K-1}) δ - 0.5 ωK δ^2 }  ×  N^+(δ | m0, v0) dδ
# Returns the FULL log-integral including prior normalizer pieces.
# Closed-form log ∫_{δ>0} exp{ (yK - ωK ψ_{K-1}) δ - 0.5 ωK δ^2 } × N^+(δ | m0, v0) dδ
# FULL log-integral (same target as numeric/MC), with the 0.5*log(2π) REMOVED.
logI_closed_sst <- function(yK, omegaK, psi_prev, m0, v0) {
  lam <- omegaK + 1/v0
  mu  <- (yK - omegaK*psi_prev + m0/v0) / lam
  # truncation masses
  num <- pnorm(mu * sqrt(lam), log.p = TRUE)
  den <- pnorm(m0 / sqrt(v0),  log.p = TRUE)
  # final
  -0.5*log(lam) - 0.5*log(v0) +
    ((yK - omegaK*psi_prev + m0/v0)^2)/(2*lam) - (m0^2)/(2*v0) +
    (num - den)
}

# Deterministic numerical check of the same integral using base::integrate.
# This integrates the FULL integrand (likelihood × truncated-normal prior) over δ>0.
logI_num_sst <- function(yK, omegaK, psi_prev, m0, v0) {
  a  <- yK - omegaK*psi_prev
  sd0 <- sqrt(v0)
  # truncation normalizer of the prior
  Z  <- 1 - pnorm(0, mean = m0, sd = sd0)   # = Φ(m0/sd0)
  if (Z <= 0) return(-Inf)
  
  # change of variable δ = t/(1-t), t in (0,1); jacobian = 1/(1-t)^2
  f_t <- function(t) {
    delta <- t/(1 - t)
    kernel <- exp(a*delta - 0.5*omegaK*delta^2) *
      dnorm(delta, mean = m0, sd = sd0) / Z
    kernel / (1 - t)^2
  }
  
  # integrate with high accuracy
  res <- integrate(f_t, lower = 0, upper = 1,
                   rel.tol = 1e-10, abs.tol = 1e-12, stop.on.error = FALSE)
  if (!is.null(res$message) && nzchar(res$message)) warning(res$message)
  log(res$value)
}

# ---- SST (increment form): closed-form log I for δ>0 -----------------------
# Let ψ_K = ψ_{K-1} + δ, with δ ~ N^+(m0, v0).
# Kernel in δ: exp{ L δ - 0.5 * ω δ^2 }, where L = \bar y_K - ω ψ_{K-1}.
logI_sst_delta_closed <- function(L, omega, m0, v0) {
  if (omega < 0) stop("omega must be >= 0")
  if (v0    <= 0) stop("v0 must be > 0")
  s0    <- sqrt(v0)
  prec0 <- 1 / v0
  prec  <- omega + prec0
  mstar <- (L + prec0 * m0) / prec
  z0    <- - mstar * sqrt(prec)
  val <- -0.5 * log(prec) - 0.5 * log(v0) +
    0.5 * (L + prec0 * m0)^2 / prec - 0.5 * (m0^2) * prec0 -
    log(.Norm_cdf(m0 / s0)) +
    log(.Norm_cdf(-z0))
  as.numeric(val)
}
# Sample from N^+(m0, v0) via inversion; then average L(δ) only.
# Use antithetic uniforms (u, 1-u) for variance reduction.
# MC for the SST integral with an estimated SE for log I
# returns list(log = logI_hat, se = se_logI)
logI_mc_sst <- function(yK, omegaK, psi_prev, m0, v0, nsim = 2e5, seed = 123) {
  stopifnot(nsim %% 2 == 0)
  set.seed(seed)
  a   <- yK - omegaK * psi_prev
  sd0 <- sqrt(v0)
  
  # sample δ ~ N(m0, v0) truncated to (0, ∞) via inverse-CDF
  c0  <- pnorm(0, mean = m0, sd = sd0)
  u   <- runif(nsim %/% 2); u <- c(u, 1 - u)         # antithetic
  q   <- c0 + u * (1 - c0)
  delta <- m0 + sd0 * qnorm(q)
  
  lvals <- a * delta - 0.5 * omegaK * delta^2               # log kernel ratio
  # self-normalized IS is not needed here (proposal = target base)
  # We estimate log(mean(exp(lvals))) with delta-method SE
  m <- max(lvals)
  w <- exp(lvals - m)                                        # stabilized weights
  mean_w <- mean(w)
  var_w  <- var(w)
  
  logI_hat <- log(mean_w) + m
  # delta method: Var(log mean) ≈ Var(mean_w) / mean_w^2 = var_w / (n * mean_w^2)
  se_logI  <- sqrt(var_w / (length(w) * mean_w^2))
  
  list(log = logI_hat, se = se_logI)
}


# ---- SST (increment): Monte Carlo check ------------------------------------
logI_sst_delta_mc <- function(L, omega, m0, v0, R = 20000L, seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  delta <- truncnorm::rtruncnorm(R, a = 0, b = Inf, mean = m0, sd = sqrt(v0))
  val   <- mean( exp(L * delta - 0.5 * omega * delta^2) )
  log(pmax(val, 1e-300))
}
# ---- WST: slot r directional collapsed factor (sum over existing ℓ) -------
# Inputs:
#   bar_y_vec    : length-K numeric, entries S_{rℓ} = sum s_{rℓ}(A_ij - N_ij/2)
#   bar_omega_vec: length-K numeric, entries ω_{rℓ} = sum ω_ij
#   mu0, sig2    : half-Normal prior for ψ_{rℓ}
# Returns log-product over ℓ != r.
wst_slot_logdir_closed <- function(r, bar_y_vec, bar_omega_vec, mu0, sig2) {
  K <- length(bar_y_vec)
  keep <- setdiff(seq_len(K), r)
  sapply(keep, function(l)
    logI_wst_pair_closed(S = bar_y_vec[l], omega = bar_omega_vec[l], mu0 = mu0, sig2 = sig2)
  ) |> sum()
}

wst_slot_logdir_mc <- function(r, bar_y_vec, bar_omega_vec, mu0, sig2, R = 20000L, seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  K <- length(bar_y_vec); keep <- setdiff(seq_len(K), r)
  sapply(keep, function(l)
    logI_wst_pair_mc(S = bar_y_vec[l], omega = bar_omega_vec[l], mu0 = mu0, sig2 = sig2, R = R)
  ) |> sum()
}

# ---- SST: contribution of new distance K (extreme slots only) --------------
# Inputs:
#   yK, omegaK : local sums at distance K for the *birth proposal* (0 in middle slots)
#   psiKm1     : current ψ_{K-1}
#   m0, v0     : prior for δ_{K-1} (increment)
sst_newdistance_logdir_closed <- function(yK, omegaK, psiKm1, m0, v0) {
  if (omegaK == 0 && yK == 0) return(0)      # middle slots: no K-distance dyads locally
  L <- yK - omegaK * psiKm1
  logI_sst_delta_closed(L = L, omega = omegaK, m0 = m0, v0 = v0)
}

sst_newdistance_logdir_mc <- function(yK, omegaK, psiKm1, m0, v0, R = 20000L, seed = NULL) {
  if (omegaK == 0 && yK == 0) return(0)
  L <- yK - omegaK * psiKm1
  logI_sst_delta_mc(L = L, omega = omegaK, m0 = m0, v0 = v0, R = R, seed = seed)
}
# ---- Helper: stable log(1+exp(sgn * x)) -----------------------------------
.log1pexp_signed <- function(x, sgn) if (sgn > 0) log1p(exp(x)) else log1p(exp(-x))

# ---- WST: slow log-score for assigning node i -> existing block k ----------
# Uses current (η, κ, ψ) and raw A, no precomputed tables.
# Slow, reference-accurate log-score for assigning node i -> block k (WST, pairwise ψ_{kℓ})
# Matches: log_post_fast_pair()
slow_logscore_existing_wst <- function(i, k, A, z, eta, kappa, psi_mat, alpha = 1) {
  K <- max(z)
  idxK <- seq_len(K)
  log2 <- log(2)
  
  # sums towards each block from node i
  e_out <- as.numeric(A[i, ])
  e_in  <- as.numeric(A[, i])
  Nij   <- e_out + e_in
  
  N_blk     <- numeric(K)
  A_out_blk <- numeric(K)
  if (any(Nij > 0)) {
    ell <- z[which(Nij > 0)]
    tmpN <- tapply(Nij[Nij > 0], ell, sum)
    tmpA <- tapply(e_out[Nij > 0], ell, sum)
    N_blk[as.integer(names(tmpN))]     <- tmpN
    A_out_blk[as.integer(names(tmpA))] <- tmpA
  }
  
  # prior on z
  n_minus <- tabulate(z[-i], nbins = K)
  lp <- log(alpha + n_minus[k])
  
  # Poisson/Gamma volume part (conditioned on kappa)
  eta_sum_blk <- tapply(eta[-i], factor(z[-i], levels = idxK), sum)
  eta_sum_blk[is.na(eta_sum_blk)] <- 0
  lambda_vec_k <- as.numeric(kappa %*% eta_sum_blk)   # length K
  lp <- lp - eta[i] * lambda_vec_k[k] + sum(N_blk * log(pmax(kappa[k, ], 1e-15)))
  
  # Directional part (pairwise ψ_{min,max} with sign)
  mask_cross <- idxK != k
  if (any(mask_cross)) {
    l_idx <- idxK[mask_cross]
    a <- pmin(k, l_idx); b <- pmax(k, l_idx)
    sgn <- ifelse(l_idx > k, +1, -1)
    psi_use <- psi_mat[cbind(a, b)]
    A_dir   <- A_out_blk[mask_cross]
    log1p_pos <- log1p(exp(psi_use))
    log1p_neg <- log1p(exp(-psi_use))
    lp <- lp +
      sum(sgn * A_dir * psi_use) -
      sum(N_blk[mask_cross] * ifelse(sgn == +1, log1p_pos, log1p_neg))
  }
  
  lp <- lp - N_blk[k] * log2
  lp
}


# ---- SST: slow log-score for assigning node i -> existing block k ----------
# psi is length K-1 by distance; same logic, no caches.
# Slow, reference-accurate log-score for assigning node i -> block k (SST, distance ψ)
# Matches: log_post_fast_distance()
slow_logscore_existing_sst <- function(i, k, A, z, eta, kappa, psi, alpha = 1) {
  K <- max(z)
  idxK <- seq_len(K)
  log2 <- log(2)
  
  # sums towards each block from node i
  e_out <- as.numeric(A[i, ])
  e_in  <- as.numeric(A[, i])
  Nij   <- e_out + e_in
  
  N_blk     <- numeric(K)
  A_out_blk <- numeric(K)
  if (any(Nij > 0)) {
    ell <- z[which(Nij > 0)]
    tmpN <- tapply(Nij[Nij > 0], ell, sum)
    tmpA <- tapply(e_out[Nij > 0], ell, sum)
    N_blk[as.integer(names(tmpN))]     <- tmpN
    A_out_blk[as.integer(names(tmpA))] <- tmpA
  }
  
  # prior on z
  n_minus <- tabulate(z[-i], nbins = K)
  lp <- log(alpha + n_minus[k])
  
  # Poisson/Gamma (volume) part, conditioned on kappa:
  # - eta_i * sum_l kappa_{k l} * sum_{j:z_j=l, j!=i} eta_j + sum_l N_blk[l] * log kappa_{k l}
  eta_sum_blk <- tapply(eta[-i], factor(z[-i], levels = idxK), sum)
  eta_sum_blk[is.na(eta_sum_blk)] <- 0
  lambda_vec_k <- as.numeric(kappa %*% eta_sum_blk)   # length K
  lp <- lp - eta[i] * lambda_vec_k[k] + sum(N_blk * log(pmax(kappa[k, ], 1e-15)))
  
  # Directional (Binomial) part:
  # for l != k: A_out contributes with sign * ψ_d, and N_blk uses log(1+exp(sign*ψ_d))
  mask_cross <- idxK != k
  if (any(mask_cross)) {
    l_idx   <- idxK[mask_cross]
    d_vec   <- abs(l_idx - k)
    signVec <- ifelse(l_idx > k, +1, -1)
    A_dir   <- A_out_blk[mask_cross]
    # log(1 + exp(±ψ_d)) via stable computation
    log1p_pos <- log1p(exp(psi[d_vec]))  # for +ψ
    log1p_neg <- log1p(exp(-psi[d_vec])) # for -ψ
    lp <- lp +
      sum(signVec * A_dir * psi[d_vec]) -
      sum(N_blk[mask_cross] * ifelse(signVec == +1, log1p_pos, log1p_neg))
  }
  
  # within-block Bernoulli(1/2): sum_j N_ij * log(1/2) = -N_blk[k] * log 2
  lp <- lp - N_blk[k] * log2
  
  lp
}
check_integral_sst <- function(yK, omegaK, psi_prev, m0, v0,
                               tol_num = 1e-6,
                               mc_sigma_mult = 4,  # pass if within 4*SE
                               nsim = 2e5, seed = 42) {
  lc <- logI_closed_sst(yK, omegaK, psi_prev, m0, v0)
  ln <- logI_num_sst   (yK, omegaK, psi_prev, m0, v0)
  mc <- logI_mc_sst    (yK, omegaK, psi_prev, m0, v0, nsim, seed)
  
  cat("SST log-integral (closed, numeric, MC):\n")
  print(c(closed   = lc,
          numeric  = ln,
          mc       = mc$log,
          diff_num = lc - ln,
          diff_mc  = lc - mc$log,
          se_mc    = mc$se))
  
  stopifnot(abs(lc - ln) < tol_num)
  stopifnot(abs(lc - mc$log) <= mc_sigma_mult * mc$se)
  invisible(list(closed = lc, numeric = ln, mc = mc))
}



# ---- One-shot validator (WST) ----------------------------------------------
check_wst_slot_directional <- function(r, bar_y_vec, bar_omega_vec, mu0, sig2,
                                       R = 50000L, seed = 1) {
  a <- wst_slot_logdir_closed(r, bar_y_vec, bar_omega_vec, mu0, sig2)
  m <- wst_slot_logdir_mc(r, bar_y_vec, bar_omega_vec, mu0, sig2, R = R, seed = seed)
  c(closed = a, mc = m, diff = a - m)
}

# ---- One-shot validator (SST) ----------------------------------------------
check_sst_newdistance <- function(yK, omegaK, psiKm1, m0, v0, R = 50000L, seed = 1) {
  a <- sst_newdistance_logdir_closed(yK, omegaK, psiKm1, m0, v0)
  m <- sst_newdistance_logdir_mc(yK, omegaK, psiKm1, m0, v0, R = R, seed = seed)
  c(closed = a, mc = m, diff = a - m)
}



# =========================
# Step 1 — Validation harness
# =========================
run_step1_checks <- function(
    seed = 2025,
    n = 40, K = 3,
    tol_fast_vs_slow = 5e-3,    # acceptable abs diff (log-scale) for fast vs slow
    tol_closed_vs_mc = 5e-3,    # acceptable abs diff (log-scale) for closed vs MC
    Rmc = 5e4                   # MC draws for integral checks
){
  set.seed(seed)
  # ---- deps ----
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Need Matrix")
  if (!requireNamespace("truncnorm", quietly = TRUE)) stop("Need truncnorm")
  if (!requireNamespace("BayesLogit", quietly = TRUE)) stop("Need BayesLogit")
  
  # ---- simulate a small SST-style OSBM (then reuse for WST by mapping ψ_{kl}=ψ_{|k-l|}) ----
  z_t <- rep(1:K, length.out = n)
  z_t <- sample(z_t)  # shuffle nodes across blocks
  
  # degree factors (blockwise normalised)
  eta_t <- runif(n, 0.8, 1.2)
  for (k in 1:K) {
    idx <- which(z_t == k)
    if (length(idx)) eta_t[idx] <- length(idx) * eta_t[idx] / sum(eta_t[idx])
  }
  
  # symmetric κ (volumes) with mild between-block heterogeneity
  kappa_t <- matrix(0, K, K)
  for (k in 1:K) for (l in k:K) {
    kappa_t[k, l] <- 3 + (k == l) * 2 + abs(k - l) * 0.5
    kappa_t[l, k] <- kappa_t[k, l]
  }
  
  # SST distance ψ (strictly increasing)
  psi_dist_t <- seq(0.8, 2.2, length.out = K - 1)
  
  # simulate A
  A <- Matrix::Matrix(0, n, n, sparse = TRUE)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    lam <- eta_t[i] * eta_t[j] * kappa_t[z_t[i], z_t[j]]
    N_ij <- rpois(1, lam)
    if (N_ij > 0) {
      d <- abs(z_t[i] - z_t[j])
      p <- if (d == 0) 0.5 else {
        sgn <- ifelse(z_t[i] < z_t[j], +1, -1)  # forward = lower index -> higher index
        plogis(sgn * psi_dist_t[d])
      }
      a_ij <- rbinom(1, N_ij, p)
      A[i, j] <- a_ij
      A[j, i] <- N_ij - a_ij
    }
  }
  
  # ==== 1) EXISTING-BLOCK SCORE: FAST vs SLOW ====
  message("1) Fast vs slow existing-block scores…")
  
  # ---- SST (distance) fast vs slow ----
  # caches for fast scorer
  log1p_tab_sst <- make_log1p_tables(psi_dist_t)
  log_kappa     <- log(pmax(kappa_t, 1e-15))
  alpha_vec     <- rep(1, K)
  
  # pick a random node
  i <- sample.int(n, 1)
  eta <- eta_t; z <- z_t; kappa <- kappa_t; psi <- psi_dist_t
  
  # compute n_minus, eta_sum_blk for slow (internally recomputed anyway)
  # compare all candidate k
  diffs_sst <- numeric(K)
  for (k in 1:K) {
    slow <- slow_logscore_existing_sst(i, k, as.matrix(A), z, eta, kappa, psi, alpha = 1)
    fast <- log_post_fast_distance(i, k,
                                   N_blk       = { # construct once:
                                     e_out <- as.numeric(A[i, ])
                                     e_in  <- as.numeric(A[, i])
                                     Nij   <- e_out + e_in
                                     v <- numeric(K)
                                     if (any(Nij > 0)) {
                                       ell <- z[which(Nij > 0)]
                                       tmp <- tapply(Nij[Nij > 0], ell, sum)
                                       v[as.integer(names(tmp))] <- tmp
                                     }
                                     v
                                   },
                                   A_out_blk   = {
                                     e_out <- as.numeric(A[i, ])
                                     e_in  <- as.numeric(A[, i])
                                     Nij   <- e_out + e_in
                                     v <- numeric(K)
                                     if (any(Nij > 0)) {
                                       ell <- z[which(Nij > 0)]
                                       tmp <- tapply(e_out[Nij > 0], ell, sum)
                                       v[as.integer(names(tmp))] <- tmp
                                     }
                                     v
                                   },
                                   n_minus     = tabulate(z[-i], nbins = K),
                                   eta_sum_blk = {
                                     v <- tapply(eta[-i], factor(z[-i], levels = 1:K), sum)
                                     v[is.na(v)] <- 0
                                     as.numeric(v)
                                   },
                                   eta_i       = eta[i],
                                   kappa       = kappa,
                                   psi         = psi,
                                   alpha_vec   = alpha_vec,
                                   log_kappa   = log_kappa,
                                   log1p_tab   = log1p_tab_sst)
    diffs_sst[k] <- abs(slow - fast)
  }
  print(diffs_sst)
  stopifnot(max(diffs_sst) < tol_fast_vs_slow)
  
  # ---- WST (pair) fast vs slow ----
  # map a pairwise ψ matrix from distance ψ (upper-tri)
  psi_pair <- matrix(0, K, K)
  ut <- which(upper.tri(psi_pair), arr.ind = TRUE)
  psi_pair[ut] <- psi_dist_t[ut[, 2] - ut[, 1]]
  log1p_tab_wst <- make_log1p_tables_pair(psi_pair)
  
  diffs_wst <- numeric(K)
  for (k in 1:K) {
    slow <- slow_logscore_existing_wst(i, k, as.matrix(A), z, eta, kappa, psi_pair, alpha = 1)
    fast <- log_post_fast_pair(i, k,
                               N_blk       = { # reuse build above for consistency
                                 e_out <- as.numeric(A[i, ])
                                 e_in  <- as.numeric(A[, i])
                                 Nij   <- e_out + e_in
                                 v <- numeric(K)
                                 if (any(Nij > 0)) {
                                   ell <- z[which(Nij > 0)]
                                   tmp <- tapply(Nij[Nij > 0], ell, sum)
                                   v[as.integer(names(tmp))] <- tmp
                                 }
                                 v
                               },
                               A_out_blk   = {
                                 e_out <- as.numeric(A[i, ])
                                 e_in  <- as.numeric(A[, i])
                                 Nij   <- e_out + e_in
                                 v <- numeric(K)
                                 if (any(Nij > 0)) {
                                   ell <- z[which(Nij > 0)]
                                   tmp <- tapply(e_out[Nij > 0], ell, sum)
                                   v[as.integer(names(tmp))] <- tmp
                                 }
                                 v
                               },
                               n_minus     = tabulate(z[-i], nbins = K),
                               eta_sum_blk = {
                                 v <- tapply(eta[-i], factor(z[-i], levels = 1:K), sum)
                                 v[is.na(v)] <- 0
                                 as.numeric(v)
                               },
                               eta_i       = eta[i],
                               kappa       = kappa,
                               psi_mat     = psi_pair,
                               alpha_vec   = alpha_vec,
                               log_kappa   = log_kappa,
                               log1p_tab_pair = log1p_tab_wst)
    diffs_wst[k] <- abs(slow - fast)
  }
  print(diffs_wst)
  stopifnot(max(diffs_wst) < tol_fast_vs_slow)
  
  message("✓ Existing-block scores (SST/WST) agree within tolerance.")
  
  # ==== 2) CLOSED-FORM vs MC for WST slot & SST new-distance ====
  message("2) Closed-form vs Monte Carlo integral checks…")
  
  # -- Build WST aggregates (bar_y, B) from data; then PG to get bar_omega
  idx <- which((A + Matrix::t(A)) > 0 & upper.tri(A + Matrix::t(A)), arr.ind = TRUE)
  i_idx <- idx[, 1]; j_idx <- idx[, 2]
  N_edge <- (A + Matrix::t(A))[idx]
  bar_y_mat <- aggregate_by_pair(K,
                                 z_i    = z[i_idx],
                                 z_j    = z[j_idx],
                                 A_ij   = as.matrix(A)[cbind(i_idx, j_idx)],
                                 N_edge = N_edge)
  B_mat <- pair_totals(K, z_i = z[i_idx], z_j = z[j_idx], N_edge = N_edge)
  bar_omega_mat <- draw_omega_pair(B_mat, psi_mat = psi_pair)   # uses BayesLogit
  
  # pick a random slot r
  r <- sample.int(K, 1)
  bar_y_vec    <- bar_y_mat[r, ]
  bar_omega_vec<- bar_omega_mat[r, ]
  # ignore self-index
  res_wst <- check_wst_slot_directional(r, bar_y_vec, bar_omega_vec,
                                        mu0 = 0.5, sig2 = 2^2, R = Rmc, seed = seed + 1)
  print(res_wst)
  stopifnot(abs(res_wst["diff"]) < tol_closed_vs_mc)
  
  # -- SST new-distance (extreme slot) — use synthetic but consistent summaries
  # Construct yK, omegaK from random values or from aggregates; here: synthetic.
  # Example test values (the ones your harness uses)
  # inside run_step1_checks(), replace the SST integral check section with:
  
  yK <- 1.7; omegaK <- 4.3; psi_prev <- 0.6; m0 <- 0.5; v0 <- 2^2
  res_sst <- check_integral_sst(yK, omegaK, psi_prev, m0, v0,
                            tol_num = 1e-6, mc_sigma_mult = 4,
                            nsim = 2e5, seed = 42)
  
  
  
  
  message("✓ Closed-form vs MC integrals (WST/SST) agree within tolerance.")
  
  invisible(list(
    diffs_sst = diffs_sst,
    diffs_wst = diffs_wst,
    wst_integral = res_wst,
    sst_integral = res_sst
  ))
}

# ---- run once ----
out_chk <- run_step1_checks()
str(out_chk)
