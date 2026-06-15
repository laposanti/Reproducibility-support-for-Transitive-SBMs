###########################################################################
# Ordered SBM – internal consistency checks
# 
###########################################################################

poisson_term_brute <- function(k, N_blk, eta_blk, eta_i, kappa) {
  acc <- 0
  for (ell in seq_along(N_blk)) {
    lam <- kappa[k, ell]
    acc <- acc + N_blk[ell] * log(pmax(lam, 1e-15)) -
      eta_i * lam * eta_blk[ell]
  }
  acc
}

binomial_term_brute <- function(k, N_blk, A_out_blk, psi) {
  acc <- - N_blk[k] * log(2)
  for (ell in seq_along(N_blk)) {
    if (ell == k) next
    d <- abs(k - ell)
    if (k < ell) {
      acc <- acc + A_out_blk[ell] *  psi[d] -
        N_blk[ell]     * log1pexp( psi[d])
    } else {
      acc <- acc - A_out_blk[ell] *  psi[d] -
        N_blk[ell]     * log1pexp(-psi[d])
    }
  }
  acc
}

aggregate_stats_brute <- function(i, A, z, eta, K) {
  N_blk <- A_out_blk <- eta_blk <- numeric(K)
  for (j in seq_len(nrow(A))) {
    if (j == i) next
    ell <- z[j]
    Nij <- A[i, j] + A[j, i]
    N_blk[ell]     <- N_blk[ell] + Nij
    A_out_blk[ell] <- A_out_blk[ell] + A[i, j]
    eta_blk[ell]   <- eta_blk[ell] + eta[j]
  }
  list(N_blk = N_blk, A_out_blk = A_out_blk, eta_blk = eta_blk)
}

run_osbm_checks <- function(A, z, eta, kappa, psi, alpha_vec,
                            tol = 1e-10, verbose = FALSE) {
  K          <- length(alpha_vec)
  log1p_tab  <- make_log1p_tables(psi)
  log_kappa  <- log(pmax(kappa, 1e-15))
  
  for (i in seq_len(nrow(A))) {
    stats       <- aggregate_stats_brute(i, A, z, eta, K)
    n_minus     <- tabulate(z[-i], nbins = K)
    eta_sum_blk <- tapply(eta[-i], factor(z[-i], levels = 1:K), sum)
    eta_sum_blk[is.na(eta_sum_blk)] <- 0
    
    lp_fast <- sapply(1:K, log_post_fast, i = i,
                      N_blk       = stats$N_blk,
                      A_out_blk   = stats$A_out_blk,
                      n_minus     = n_minus,
                      eta_sum_blk = eta_sum_blk,
                      eta_i       = eta[i],
                      kappa       = kappa,
                      psi         = psi,
                      alpha_vec   = alpha_vec,
                      log_kappa   = log_kappa,
                      log1p_tab   = log1p_tab)
    
    lp_brute <- sapply(1:K, function(k) {
      lp_prior <- log(alpha_vec[k] + n_minus[k])
      lp_pois  <- poisson_term_brute(k, stats$N_blk,
                                     eta_sum_blk, eta[i], kappa)
      lp_binom <- binomial_term_brute(k, stats$N_blk,
                                      stats$A_out_blk, psi)
      lp_prior + lp_pois + lp_binom
    })
    
    if (max(abs(lp_fast - lp_brute)) > tol)
      stop(sprintf("Consistency check failed at node %d", i))
  }
  
  if (verbose) cat("All fast/brute checks passed (tol =", tol, ").\n")
  invisible(TRUE)
}
