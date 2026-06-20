# Variable-K directed DCSBM helper functions.
#
# This file contains prior, normalisation, sufficient-statistic, and
# single-site update helpers extracted from core/DCSBM_varK.R.


# ============================================================
# Partition prior weights
# ============================================================
gnedin_prior_weights <- function(n_minus, n_k_minus, gamma) {
  active <- which(n_k_minus > 0L)
  k_active <- length(active)
  if (k_active == 0L) {
    return(list(active = integer(0), w_active = numeric(0), w_new = 1))
  }
  w_active <- (n_k_minus[active] + 1) * (n_minus - k_active + gamma)
  w_new <- k_active * (k_active - gamma)
  list(active = active, w_active = w_active, w_new = w_new)
}

crp_prior_weights <- function(n_k_minus, alpha_crp) {
  active <- which(n_k_minus > 0L)
  k_active <- length(active)
  if (k_active == 0L) {
    return(list(active = integer(0), w_active = numeric(0), w_new = alpha_crp))
  }
  w_active <- n_k_minus[active]
  w_new    <- alpha_crp
  list(active = active, w_active = w_active, w_new = w_new)
}

# ============================================================
# Utilities
# ============================================================
safe_log <- function(x) log(pmax(x, .Machine$double.xmin))

.safe_rate <- function(x) {
  if (!is.finite(x) || x <= 0) return(1e-8)
  x
}
.safe_shape <- function(x) {
  if (!is.finite(x) || x <= 0) return(1e-8)
  x
}

# ============================================================
# Shrink blocks (now only one theta)
# ============================================================
shrink_blocks <- function(z, lambda, theta) {
  labs <- sort(unique(z))
  K_new <- length(labs)
  z_new <- match(z, labs)
  lambda_new <- lambda[labs, labs, drop = FALSE]
  list(
    z      = z_new,
    lambda = lambda_new,
    theta  = theta,
    K      = K_new
  )
}

# ============================================================
# Normalisation: within-block sum(theta_i)=n_k
# (matches the identifiability choice in the LaTeX)
# ============================================================
normalize_block_theta <- function(z, theta) {
  K <- max(z, na.rm = TRUE)
  for (k in seq_len(K)) {
    idx <- which(z == k)
    if (!length(idx)) next
    s <- sum(theta[idx])
    n_k <- length(idx)
    if (is.finite(s) && s > 0) theta[idx] <- theta[idx] * (n_k / s)
  }
  theta
}

# ============================================================
# Sufficient statistics for lambda update
# R_{k\ell} and T_{k\ell}(theta) = sum_{i in k} sum_{j in ell, j!=i} theta_i theta_j
# ============================================================
R_counts <- function(A, z, K) {
  n <- nrow(A)
  R <- matrix(0, K, K)
  for (k in seq_len(K)) {
    I <- which(z == k)
    if (!length(I)) next
    for (l in seq_len(K)) {
      J <- which(z == l)
      if (!length(J)) next
      if (k != l) {
        R[k, l] <- sum(A[I, J, drop = FALSE])
      } else {
        # exclude i=j
        R[k, l] <- sum(A[I, I, drop = FALSE]) - sum(diag(A[I, I, drop = FALSE]))
      }
    }
  }
  R
}

T_exposures <- function(z, theta, K) {
  S <- tapply(theta, z, sum)
  S <- as.numeric(S)[seq_len(K)]
  S[is.na(S)] <- 0
  
  diag_sum <- tapply(theta^2, z, sum)
  diag_sum <- as.numeric(diag_sum)[seq_len(K)]
  diag_sum[is.na(diag_sum)] <- 0
  
  T <- outer(S, S, "*")
  diag(T) <- diag(T) - diag_sum
  T[!is.finite(T)] <- 0
  T <- pmax(T, 0)
  T
}

# ============================================================
# Parameter updates
# Model: A_ij ~ Pois(theta_i * theta_j * lambda[z_i, z_j]), i!=j
# Priors: theta_i ~ Gamma(a_eta, b_eta), lambda_kl ~ Gamma(a_lambda, b_lambda)
# ============================================================
update_lambda <- function(A, z, theta, a_lambda, b_lambda, K) {
  R <- R_counts(A, z, K)
  T <- T_exposures(z, theta, K)
  shape <- a_lambda + R
  rate  <- b_lambda + T
  matrix(
    rgamma(K * K, shape = as.vector(shape), rate = as.vector(rate)),
    K, K
  )
}

update_theta <- function(A, z, theta, lambda, a_eta, b_eta) {
  # Full conditional:
  # theta_i | rest ~ Gamma(a_eta + G_i, b_eta + sum_ell (lambda_{k ell}+lambda_{ell k}) * S^{-i}_ell)
  # where k=z_i, G_i = sum_{j!=i}(A_ij + A_ji), S^{-i}_ell = sum_{j:z_j=ell, j!=i} theta_j
  n <- nrow(A)
  K <- nrow(lambda)
  
  G <- rowSums(A) + colSums(A)
  
  # block sums S_ell = sum_{j:z_j=ell} theta_j
  S <- tapply(theta, z, sum)
  S <- as.numeric(S)[seq_len(K)]
  S[is.na(S)] <- 0
  
  for (i in seq_len(n)) {
    k <- z[i]
    
    # S^{-i}_ell: subtract theta_i from its own block
    S_minus_i <- S
    S_minus_i[k] <- S_minus_i[k] - theta[i]
    S_minus_i <- pmax(S_minus_i, 0)
    
    rate_i <- b_eta + sum((lambda[k, ] + lambda[, k]) * S_minus_i, na.rm = TRUE)
    rate_i <- .safe_rate(rate_i)
    
    shape_i <- .safe_shape(a_eta + G[i])
    val <- rgamma(1L, shape = shape_i, rate = rate_i)
    theta[i] <- if (is.finite(val) && val > 0) val else 1e-6
  }
  
  theta
}

# ============================================================
# New-block lambda instantiation from full conditional (single theta)
# ============================================================
sample_lambda_new_block_from_fc <- function(
    i, A, z, theta,
    lambda_old,
    a_lambda, b_lambda
) {
  n      <- nrow(A)
  K_old  <- nrow(lambda_old)
  K_new  <- K_old + 1L
  
  lambda_new <- matrix(0, K_new, K_new)
  lambda_new[1:K_old, 1:K_old] <- lambda_old
  
  k_new <- K_new
  
  for (ell in seq_len(K_old)) {
    # nodes in block ell excluding i
    J <- which(z == ell & seq_len(n) != i)
    
    R_out <- if (length(J)) sum(A[i, J]) else 0
    R_in  <- if (length(J)) sum(A[J, i]) else 0
    
    # exposure is the same form in both directions: theta_i * sum_{j in ell} theta_j
    T_i_ell <- if (length(J)) theta[i] * sum(theta[J]) else 0
    T_i_ell <- .safe_rate(T_i_ell)
    
    lambda_new[k_new, ell] <- rgamma(1L, a_lambda + R_out, b_lambda + T_i_ell)
    lambda_new[ell, k_new] <- rgamma(1L, a_lambda + R_in,  b_lambda + T_i_ell)
  }
  
  # singleton self-interaction (no self-edges anyway) => prior
  lambda_new[k_new, k_new] <- rgamma(1L, a_lambda, b_lambda)
  
  lambda_new
}

# ============================================================
# z-update likelihood pieces (single theta)
# ============================================================
loglik_i_in_block <- function(i, k, A, z, theta, lambda) {
  n <- nrow(A)
  idx <- setdiff(seq_len(n), i)
  
  # outgoing edges i -> j
  A_out <- A[i, idx]
  z_out <- z[idx]
  theta_j <- theta[idx]
  lam_out <- lambda[k, z_out]
  
  # incoming edges j -> i
  A_in <- A[idx, i]
  z_in <- z[idx]
  lam_in <- lambda[z_in, k]
  
  ll_out <- sum(
    A_out * (safe_log(theta[i]) + safe_log(theta_j) + safe_log(lam_out)) -
      theta[i] * theta_j * lam_out,
    na.rm = TRUE
  )
  
  ll_in <- sum(
    A_in * (safe_log(theta_j) + safe_log(theta[i]) + safe_log(lam_in)) -
      theta_j * theta[i] * lam_in,
    na.rm = TRUE
  )
  
  ll_out + ll_in
}

loglik_i_new_block_collapsed <- function(
    i, A, z, theta,
    a_lambda, b_lambda,
    include_factorial = TRUE
) {
  n <- nrow(A)
  K_curr <- length(unique(z))  # assuming labels are 1:K
  
  loglik <- 0
  
  for (ell in seq_len(K_curr)) {
    J <- which(z == ell & seq_len(n) != i)
    if (length(J) == 0L) next
    
    R_out <- sum(A[i, J])
    R_in  <- sum(A[J, i])
    
    T_i_ell <- theta[i] * sum(theta[J])
    T_i_ell <- .safe_rate(T_i_ell)
    
    # outgoing marginal
    loglik <- loglik +
      (lgamma(a_lambda + R_out) - lgamma(a_lambda)) +
      a_lambda * log(b_lambda) -
      (a_lambda + R_out) * log(b_lambda + T_i_ell)
    if (include_factorial) loglik <- loglik - lgamma(R_out + 1)
    
    # incoming marginal (same exposure)
    loglik <- loglik +
      (lgamma(a_lambda + R_in) - lgamma(a_lambda)) +
      a_lambda * log(b_lambda) -
      (a_lambda + R_in) * log(b_lambda + T_i_ell)
    if (include_factorial) loglik <- loglik - lgamma(R_in + 1)
  }
  
  loglik
}

update_z_gnedin_varK <- function(
    A, z, theta, lambda,
    gamma_gnedin = 0.8,
    alpha_crp = 1.0,
    partition_prior = "GN",
    a_lambda, b_lambda
) {
  n <- nrow(A)
  K_full <- nrow(lambda)
  
  for (i in sample.int(n)) {
    n_minus   <- n - 1L
    n_k_minus <- tabulate(z[-i], nbins = K_full)
    
    if (partition_prior == "CRP") {
      pw <- crp_prior_weights(
        n_k_minus = n_k_minus,
        alpha_crp = alpha_crp
      )
    } else {
      pw <- gnedin_prior_weights(
        n_minus   = n_minus,
        n_k_minus = n_k_minus,
        gamma     = gamma_gnedin
      )
    }
    active   <- pw$active
    k_active <- length(active)
    
    logps  <- c()
    labels <- c()
    
    if (k_active > 0L) {
      log_prior_active <- safe_log(pw$w_active)
      
      for (idx in seq_along(active)) {
        k <- active[idx]
        ll_k <- loglik_i_in_block(i, k, A, z, theta, lambda)
        logps  <- c(logps, log_prior_active[idx] + ll_k)
        labels <- c(labels, k)
      }
      
      if (pw$w_new > 0) {
        log_prior_new <- safe_log(pw$w_new)
        ll_new <- loglik_i_new_block_collapsed(
          i, A, z, theta,
          a_lambda = a_lambda, b_lambda = b_lambda
        )
        logps  <- c(logps, log_prior_new + ll_new)
        labels <- c(labels, K_full + 1L)
      }
    } else {
      ll_new <- loglik_i_new_block_collapsed(
        i, A, z, theta,
        a_lambda = a_lambda, b_lambda = b_lambda
      )
      logps  <- ll_new
      labels <- K_full + 1L
    }
    
    m <- max(logps)
    p <- exp(logps - m)
    p <- p / sum(p)
    new_label <- sample(labels, size = 1L, prob = p)
    
    if (new_label > K_full) {
      lambda <- sample_lambda_new_block_from_fc(
        i, A, z, theta,
        lambda_old = lambda,
        a_lambda   = a_lambda,
        b_lambda   = b_lambda
      )
      K_full <- nrow(lambda)
    }
    
    z[i] <- new_label
  }
  
  list(z = z, lambda = lambda)
}

# ============================================================
# Main Gibbs sampler (single theta)
# ============================================================
