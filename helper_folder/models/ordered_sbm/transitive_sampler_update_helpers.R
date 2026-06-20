# Ordered SBM sampler helpers: assertion and SST single-site update logic.
# These functions were extracted from core/transitive_sbm_sampler.R to keep
# the main sampler implementation focused on iteration control.

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
