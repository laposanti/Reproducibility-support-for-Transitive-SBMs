.reconcile_after_z_change <- function(z, kappa, psi, psi_constraint){
  stopifnot(length(z) > 0L)
  present <- sort(unique(as.integer(z)))
  Kz <- length(present)
  
  ## ----- κ: permute to 'present' (not just delete), then size Kz x Kz
  if (!identical(present, seq_len(nrow(kappa)))) {
    if (any(present > nrow(kappa)))
      stop("Reconcile: some labels exceed nrow(kappa).")
    kappa <- kappa[present, present, drop = FALSE]
  }
  if (nrow(kappa) != Kz || ncol(kappa) != Kz)
    stop("Reconcile: κ and z disagree after permutation/shrink.")
  
  ## ----- ψ
  if (psi_constraint == "WST") {
    # κ/ψ share the same block order; permute ψ the same way and keep upper-tri semantics
    if (!identical(present, seq_len(nrow(psi)))) {
      if (any(present > nrow(psi)))
        stop("Reconcile: some labels exceed nrow(psi) in WST.")
      psi <- psi[present, present, drop = FALSE]
    }
    if (nrow(psi) != Kz || ncol(psi) != Kz)
      stop("Reconcile: ψ(WST) and z disagree after permutation/shrink.")
  } else {
    # SST: distance ψ has length K-1; trim/pad to Kz-1
    want <- max(0L, Kz - 1L)
    have <- length(psi)
    if (have > want)      psi <- psi[seq_len(want)]
    else if (have < want) psi <- c(psi, rep(0, want - have))
  }
  
  list(kappa = kappa, psi = psi, K = Kz)
}




## ---- NEW: EPPF predictive weights ----
.pred_weights <- function(n_k, eppf, hyp){
  K <- length(n_k); n <- sum(n_k)
  if (eppf == "DM") {
    # DM/CRP with alpha0: p(join k) ∝ n_k, p(new) ∝ alpha0
    list(w_k = n_k, w_new = hyp$alpha0)
  } else if (eppf == "PY") {
    sigma <- hyp$sigma; theta <- hyp$theta
    stopifnot(theta > -sigma, sigma < 1)
    list(w_k = pmax(n_k - sigma, 0), w_new = theta + K * sigma)
  } else { # "GN" (simple, commonly used form)
    # One convenient Gnedin parametrization: p(new) ∝ gamma, p(join k) ∝ n_k + gamma/K
    g <- hyp$gamma
    list(w_k = n_k + g / max(K,1), w_new = g)
  }
}

## ---- NEW: slot set per node i ----
.slot_set <- function(K, k_prime, it, burn, policy){
  if (policy == "all" || it <= burn) return(seq_len(K+1L))
  # adjacent-only after burn: before/after the block where i came from (k')
  c(max(1L, k_prime), min(K+1L, k_prime + 1L))
}
## ===== Collapsed Poisson/Gamma volume (slot-invariant) ======================
collapsed_pois_volume_i <- function(i, A, z, eta, K, a_kappa, b_kappa){
  # R_{iℓ} = sum_j N_ij I(z_j = ℓ);  T_{iℓ} = η_i * sum_j η_j I(z_j = ℓ)
  e_out <- A[i, ]; e_in <- A[, i]; Nij <- e_out + e_in
  N_blk <- numeric(K)
  if (any(Nij > 0)) {
    ell <- z[which(Nij > 0)]
    N_tmp <- tapply(Nij[Nij > 0], ell, sum)
    N_blk[as.integer(names(N_tmp))] <- N_tmp
  }
  eta_sum <- tapply(eta, factor(z, levels = seq_len(K)), sum); eta_sum[is.na(eta_sum)] <- 0
  T_blk   <- eta[i] * eta_sum
  # ∏_ℓ Γ(a+R)/Γ(a) * b^a / (b+T)^(a+R)  (return log for stability)
  logC <- lgamma(a_kappa + N_blk) - lgamma(a_kappa) +
    a_kappa * log(b_kappa) - (a_kappa + N_blk) * log(b_kappa + T_blk)
  sum(logC)
}

## ===== Stable expected Polya–Gamma for proposals ============================
# E[ω | h=N, z] = N/(2z) * tanh(z/2); limit z→0 gives N/4.
.pg_expectation <- function(N, z){
  z <- as.numeric(z)
  N <- as.numeric(N)
  out <- numeric(length(N))
  small <- abs(z) < 1e-8
  out[small]  <- 0.25 * N[small]
  out[!small] <- 0.5 * N[!small] / z[!small] * tanh(z[!small] / 2)
  out
}
collapsed_dir_WST_slot_r <- function(
    i, r, A, z, K, 
    mu0, sigma0,                           # prior N^+(μ0, σ0^2)
    use_pg_expectation = TRUE, psi_ref = 0 # ẑ for ω; default 0 ⇒ N/4
){
  # For each ℓ, collect signed bar_y and bar_omega over j in block ℓ
  e_out <- A[i, ]; e_in <- A[, i]; Nij <- e_out + e_in
  idxJ  <- which(Nij > 0)
  if (!length(idxJ)) return(0) # log(1)
  zj    <- z[idxJ]
  Aj    <- e_out[idxJ]; Nj <- Nij[idxJ]
  # sign pattern s_{rℓ}: +1 if ℓ > r (forward), -1 if ℓ < r
  # bar_y_{rℓ} = Σ_j s_{rℓ}(A_ij - N_ij/2);  bar_omega_{rℓ} = Σ_j ω_ij
  logC <- 0
  for (ell in seq_len(K)) {
    J <- idxJ[zj == ell]
    if (!length(J)) next
    sgn <- if (ell > r) +1 else -1
    ybar <- sum( sgn * (A[i, J] - Nij[J]/2) )
    # omega: either expectation at ẑ, or draw PG; for proposals expectation is robust.
    if (use_pg_expectation) {
      zhat  <- abs(psi_ref) # any fixed reference; ẑ=0 ⇒ N/4
      obar  <- sum(.pg_expectation(Nij[J], zhat))
    } else {
      obar <- sum(BayesLogit::rpg(length(J), h = Nij[J], z = abs(psi_ref)))
    }
    # scalar half-line Gaussian integral with truncation (your Step 5)
    prec  <- obar + 1/(sigma0^2)
    meanN <- (ybar + (mu0/(sigma0^2))) / prec
    zcdf  <- (ybar + (mu0/(sigma0^2))) / sqrt(prec)
    logI  <- -0.5*log(prec) - log(sigma0) +
      0.5 * (ybar + mu0/(sigma0^2))^2 / prec - 0.5 * (mu0^2)/(sigma0^2) +
      pnorm(zcdf, log.p = TRUE) - pnorm(mu0/sigma0, log.p = TRUE)
    logC <- logC + logI
  }
  logC
}


collapsed_dir_SST_slot_r <- function(
    i, r, A, z, K, psi_vec,           # psi_1..psi_{K-1}
    mu0, sigma0,                      # (unused for shared part; kept for uniform sig)
    mKm1, vKm1,                       # prior for increment δ_{K-1}
    use_pg_expectation = TRUE
){
  e_out <- A[i, ]; e_in <- A[, i]; Nij <- e_out + e_in
  idxJ  <- which(Nij > 0)
  if (!length(idxJ)) {
    # middle: no ψ_K term; shared part is zero; extreme: Ik reduces to prior ratio
    return(0)
  }
  zj <- z[idxJ]
  
  # distance-wise accumulators
  D <- K - 1L
  S_d <- numeric(D)
  O_d <- numeric(D)
  
  for (jj in seq_along(idxJ)) {
    j  <- idxJ[jj]
    d  <- abs(z[j] - r)
    if (d == 0 || d > D) next
    sgn <- ifelse(z[j] > r, +1, -1)
    S_d[d] <- S_d[d] + sgn * (A[i, j] - Nij[j]/2)
    if (use_pg_expectation) {
      O_d[d] <- O_d[d] + .pg_expectation(Nij[j], abs(psi_vec[d]))
    } else {
      O_d[d] <- O_d[d] + BayesLogit::rpg(1, h = Nij[j], z = abs(psi_vec[d]))
    }
  }
  
  # Shared log-kernel: Σ_d [ S_d ψ_d - ½ O_d ψ_d^2 ]
  log_shared <- sum(S_d * psi_vec) - 0.5 * sum(O_d * (psi_vec^2))
  
  # New maximal distance only at extremes
  if (r > 1 && r < K + 1) return(log_shared)
  
  # Extreme slot: build yK, oK at distance K
  yK <- 0; oK <- 0
  for (jj in seq_along(idxJ)) {
    j <- idxJ[jj]
    d <- abs(z[j] - r)
    if (d == K) {
      sgn <- ifelse(z[j] > r, +1, -1)
      yK  <- yK + sgn * (A[i, j] - Nij[j]/2)
      oK  <- oK + if (use_pg_expectation) .pg_expectation(Nij[j], 0) else BayesLogit::rpg(1, h = Nij[j], z = 0)
    }
  }
  
  # Increment half-line integral (closed form)
  prec <- oK + 1/vKm1
  num  <- (yK - oK * psi_vec[K-1] + mKm1 / vKm1)
  logIk <- -0.5*log(prec) - 0.5*log(vKm1) +
    0.5*(num^2)/prec - 0.5*(mKm1^2)/vKm1 +
    pnorm(num / sqrt(prec), log.p = TRUE) - pnorm(mKm1 / sqrt(vKm1), log.p = TRUE)
  
  as.numeric(log_shared + logIk)
}



z_update_with_birth <- function(
    i, A, z, eta, kappa, psi, alpha_vec, log_kappa,
    psi_constraint = c("WST","SST"),
    it, burn,
    eppf, eppf_hyp,
    mu0, sigma0, tau0,
    hyper,
    slot_policy = c("all","adjacent"),
    use_pg_expectation = TRUE
){
  psi_constraint <- match.arg(psi_constraint)
  slot_policy    <- match.arg(slot_policy)
  
  K <- length(alpha_vec)
  n <- length(z)
  
  # --- remove i and compress if needed
  zi <- z[i]
  dropped_singleton <- (sum(z == zi) == 1L)
  
  z_wo_i <- z; z_wo_i[i] <- NA_integer_
  keep   <- !is.na(z_wo_i)
  z_wo_i <- as.integer(factor(z_wo_i[keep], levels = sort(unique(z_wo_i[keep]))))
  z_curr <- integer(n); z_curr[keep] <- z_wo_i; z_curr[i] <- 0L
  Kc <- max(z_curr)
  
  # k' for adjacent-slot policy
  uniq <- sort(unique(z_curr[z_curr > 0L]))
  if (!dropped_singleton) {
    k_prime <- which(uniq == zi)
  } else {
    k_prime <- min(max(1L, zi), Kc + 1L)
  }
  
  # EPPF predictive weights
  n_k <- tabulate(z_curr[keep], nbins = Kc)
  pw  <- .pred_weights(n_k, eppf, eppf_hyp)
  w_k <- pw$w_k; w_new <- pw$w_new
  
  # per-block summaries for i
  e_out <- A[i, ]; e_in <- A[, i]; Nij <- e_out + e_in
  N_blk <- A_out_blk <- numeric(Kc)
  if (any(Nij > 0)) {
    ell <- z_curr[which(Nij > 0)]
    N_tmp  <- tapply(Nij[Nij > 0], ell, sum)
    A_tmp  <- tapply(e_out[Nij > 0], ell, sum)
    N_blk[as.integer(names(N_tmp))] <- N_tmp
    A_out_blk[as.integer(names(A_tmp))] <- A_tmp
  }
  n_minus     <- tabulate(z_curr[keep], nbins = Kc)
  eta_sum_blk <- tapply(eta[keep], factor(z_curr[keep], levels = seq_len(Kc)), sum)
  eta_sum_blk[is.na(eta_sum_blk)] <- 0
  
  # IMPORTANT: neutralise DM prior inside fast scorers to avoid double counting
  alpha_vec_eff <- rep(0, Kc)
  
  if (psi_constraint == "WST") {
    lp_join <- sapply(seq_len(Kc), log_post_fast_pair, i = i,
                      N_blk = N_blk, A_out_blk = A_out_blk,
                      n_minus = n_minus, eta_sum_blk = eta_sum_blk, eta_i = eta[i],
                      kappa = kappa[1:Kc, 1:Kc, drop = FALSE],
                      psi_mat = psi[1:Kc, 1:Kc, drop = FALSE],
                      alpha_vec = alpha_vec_eff,
                      log_kappa = log_kappa[1:Kc, 1:Kc, drop = FALSE],
                      log1p_tab_pair = make_log1p_tables_pair(psi[1:Kc, 1:Kc, drop = FALSE]))
  } else {
    lp_join <- sapply(seq_len(Kc), log_post_fast_distance, i = i,
                      N_blk = N_blk, A_out_blk = A_out_blk,
                      n_minus = n_minus, eta_sum_blk = eta_sum_blk, eta_i = eta[i],
                      kappa = kappa[1:Kc, 1:Kc, drop = FALSE],
                      psi = as.numeric(psi[1:max(0L, Kc-1)]),
                      alpha_vec = alpha_vec_eff,
                      log_kappa = log_kappa[1:Kc, 1:Kc, drop = FALSE],
                      log1p_tab = make_log1p_tables(as.numeric(psi[1:max(0L, Kc-1)])))
  }
  
  # Total join weight (EPPF × likelihood)
  logW_join <- matrixStats::logSumExp(log(w_k) + lp_join)
  
  # Birth-by-slot collapsed predictives
  slots <- .slot_set(Kc, k_prime, it, burn, slot_policy)
  logu  <- numeric(length(slots))
  logC_vol <- collapsed_pois_volume_i(
    i, A, z_curr, eta, Kc, a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa
  )
  for (t in seq_along(slots)) {
    r <- slots[t]
    if (psi_constraint == "WST") {
      log_dir <- collapsed_dir_WST_slot_r(
        i, r, A, z_curr, Kc, mu0 = mu0, sigma0 = sigma0,
        use_pg_expectation = use_pg_expectation, psi_ref = 0
      )
    } else {
      log_dir <- collapsed_dir_SST_slot_r(
        i, r, A, z_curr, Kc,
        psi_vec = as.numeric(psi[1:max(0L, Kc-1)]),
        mu0 = mu0, sigma0 = sigma0,
        mKm1 = mu0, vKm1 = tau0^2,
        use_pg_expectation = use_pg_expectation
      )
    }
    # uniform over eligible slots
    logu[t] <- log(1/length(slots)) + (logC_vol + log_dir)
  }
  logW_new <- log(w_new) + matrixStats::logSumExp(logu)
  
  # Stage (i): join vs new
  jj <- matrixStats::logSumExp(c(logW_join, logW_new))
  p_new <- exp(logW_new - jj)
  if (runif(1) > p_new) {
    # JOIN route
    k <- sample.int(Kc, 1L, prob = exp((log(w_k) + lp_join) -
                                         matrixStats::logSumExp(log(w_k) + lp_join)))
    z_after <- z_curr; z_after[i] <- k
    return(list(
      decision = "join", k = k, r = NA_integer_,
      z_after = z_after, dropped = dropped_singleton,
      K_after = Kc
    ))
  }
  
  # Stage (ii): choose slot r for NEW block and relabel
  r <- sample(slots, 1L, prob = exp(logu - matrixStats::logSumExp(logu)))
  z_new <- z_curr
  z_new[z_new >= r] <- z_new[z_new >= r] + 1L
  z_new[i] <- r
  
  list(
    decision = "new", k = NA_integer_, r = r,
    z_after = z_new, dropped = FALSE,
    K_after = Kc + 1L
  )
}



## Materialize κ draws to an S x Kmax x Kmax array (padded with NA)
pad_kappa_array <- function(kappa_list, fill = NA_real_) {
  stopifnot(is.list(kappa_list), length(kappa_list) > 0)
  S <- length(kappa_list)
  Kmax <- max(vapply(kappa_list, function(M) nrow(M), integer(1)))
  arr <- array(fill, dim = c(S, Kmax, Kmax))
  for (s in seq_len(S)) {
    M <- kappa_list[[s]]
    Ks <- nrow(M)
    arr[s, 1:Ks, 1:Ks] <- M
  }
  arr
}

## WST: materialize ψ to S x Kmax x Kmax array (upper-tri meaningful; pad NA)
pad_psi_array_wst <- function(psi_list, fill = NA_real_) {
  stopifnot(is.list(psi_list), length(psi_list) > 0)
  S <- length(psi_list)
  Kmax <- max(vapply(psi_list, function(M) nrow(M), integer(1)))
  arr <- array(fill, dim = c(S, Kmax, Kmax))
  for (s in seq_len(S)) {
    M <- psi_list[[s]]
    Ks <- nrow(M)
    arr[s, 1:Ks, 1:Ks] <- M
  }
  arr
}

## SST: materialize ψ (distance) to S x (Kmax-1) matrix (pad NA in tail)
pad_psi_matrix_sst <- function(psi_list, K_trace, fill = NA_real_) {
  stopifnot(is.list(psi_list), length(psi_list) == length(K_trace))
  S <- length(psi_list)
  Kmax <- max(K_trace)
  Dmax <- max(0L, Kmax - 1L)
  mat <- matrix(fill, nrow = S, ncol = Dmax)
  for (s in seq_len(S)) {
    psi_s <- psi_list[[s]]
    if (length(psi_s) == 0L) next
    Ds <- length(psi_s)
    mat[s, 1:Ds] <- psi_s
  }
  mat
}


# --- helpers ---
# replace the old insert_rowcol_sym() with this
insert_rowcol_sym <- function(M, r, offdiag, selfval) {
  # M: K x K symmetric; returns (K+1) x (K+1)
  # r: slot (1..K+1) where the new block is inserted
  # offdiag: numeric vector of length K (entries to other blocks, in order of the old labels)
  # selfval: scalar for the new diagonal element κ_rr
  K <- nrow(M)
  stopifnot(length(offdiag) == K, length(selfval) == 1)
  
  M2 <- matrix(0, K + 1L, K + 1L)
  
  # copy blocks before r
  if (r > 1L) M2[1:(r-1), 1:(r-1)] <- M[1:(r-1), 1:(r-1)]
  
  # copy blocks after r (shift by +1)
  if (r <= K) {
    M2[(r+1):(K+1), (r+1):(K+1)] <- M[r:K, r:K]
    M2[(r+1):(K+1), 1:(r-1)]     <- M[r:K, 1:(r-1)]
    M2[1:(r-1),     (r+1):(K+1)] <- M[1:(r-1), r:K]
  }
  
  # set new row/col r off-diagonals (order: columns except r)
  M2[r, -r] <- offdiag
  M2[-r, r] <- offdiag
  
  # set new diagonal
  M2[r, r] <- selfval
  M2
}
.sample_new_sst_edge_psiK <- function(A, z, r, psi_curr, mu0, tau2_0, use_pg_expectation = TRUE){
  # If r ∈ {1, K+1}, sample the new maximal-distance parameter ψ_K from its FC.
  K <- max(z)
  if (r == 1L || r == (max(z))) {
    
  
  
  # Build sufficient stats at distance d = K (only pairs that are K apart)
  yK <- 0; oK <- 0
  idx_r <- which(z == r)        # size 1 at birth
  idx_far <- which(abs(z - r) == K)  # blocks at maximum distance
  if (length(idx_r) && length(idx_far)) {
    for (i in idx_r) for (j in idx_far) {
      Nij <- A[i, j] + A[j, i]
      sgn <- ifelse(z[j] > r, +1, -1)
      yK  <- yK + sgn * (A[i, j] - Nij/2)
      oK  <- oK + if (use_pg_expectation) .pg_expectation(Nij, 0) else BayesLogit::rpg(1, h = Nij, z = 0)
    }
  }
  # Gaussian half-line posterior for the increment beyond ψ_{K-1}
  prec <- oK + 1/tau2_0
  mean_num <- (yK - oK * (if (K-1L >= 1L) psi_curr[K-1L] else 0)) + mu0/tau2_0
  m <- mean_num/prec; s <- sqrt(1/prec)
  truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = m, sd = s)
  }else NULL
}

.sample_new_wst_psi_rowcol <- function(A, z, r, mu0, sig2_0, use_pg_expectation = TRUE) {
  # Returns a vector of length K-1 with ψ_{min(r,ℓ),max(r,ℓ)} for all ℓ≠r
  K <- max(z)
  stopifnot(r >= 1L, r <= K)
  get_stats <- function(ell) {
    # Signed pseudo-responses and Ω across edges between block r and block ell
    idx_r <- which(z == r); idx_l <- which(z == ell)
    if (!length(idx_r) || !length(idx_l)) return(c(ybar = 0, obar = 0))
    # All pairs (i in r, j in ell)
    A_rl <- as.matrix(A[idx_r, idx_l, drop = FALSE])
    A_lr <- as.matrix(A[idx_l, idx_r, drop = FALSE])
    # Nij and signed (A_ij - Nij/2): sign is +1 if r<ell (forward), else -1
    Nij  <- A_rl + t(A_lr)
    sgn  <- if (r < ell) +1 else -1
    ybar <- sgn * sum(A_rl - Nij/2)
    if (use_pg_expectation) {
      obar <- sum(.pg_expectation(as.numeric(Nij), 0))  # z=0 ⇒ N/4
    } else {
      obar <- sum(BayesLogit::rpg(length(Nij), h = as.numeric(Nij), z = 0))
    }
    c(ybar = ybar, obar = obar)
  }
  
  out <- numeric(K - 1L)
  pos <- 0L
  for (ell in seq_len(K)) if (ell != r) {
    pos <- pos + 1L
    st  <- get_stats(ell)
    prec <- st["obar"] + 1/sig2_0
    mean_num <- st["ybar"] + mu0/sig2_0
    # half-N( mean = mean_num/prec, var = 1/prec ) truncated at 0
    m  <- as.numeric(mean_num/prec)
    s  <- sqrt(1/prec)
    # draw from N^+(m,s^2) (use truncnorm)
    out[pos] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = m, sd = s)
  }
  out
}

.insert_wst_psi_fullcond_at_slot <- function(psi_old, r, new_vals) {
  # Place the K-1 values into the upper-tri at slot r
  K  <- nrow(psi_old)
  K1 <- K + 1L
  stopifnot(length(new_vals) == K)
  psi_new <- matrix(0, K1, K1)
  if (K > 0) psi_new[1:K, 1:K] <- psi_old
  # map values to the new row/col
  pos <- 0L
  for (ell in seq_len(K1)) if (ell != r) {
    pos <- pos + 1L
    if (min(r, ell) < max(r, ell)) {
      if (r < ell) psi_new[r, ell] <- new_vals[pos] else psi_new[ell, r] <- new_vals[pos]
    }
  }
  psi_new
}
