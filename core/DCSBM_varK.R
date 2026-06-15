library(dplyr)
library(mcclust)
library(mcclust.ext)

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
fit_dcsbm_gibbs_gnedin <- function(
    A, K_init,
    iters = 4000, burn_in = floor(iters/2), thin = 10,
    priors = list(
      a_eta = 1, b_eta = 1,
      a_lambda = 1, b_lambda = 1,
      gamma_gnedin = 0.8,
      alpha_crp = 1.0,
      partition_prior = "GN"
    ),
    z_init = NULL,
    normalize_each_iter = TRUE,
    save_z = TRUE,
    save_theta = FALSE,
    seed = NULL, verbose = 100
) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(A)
  stopifnot(ncol(A) == n)
  diag(A) <- 0
  
  a_eta <- if (is.null(priors$a_eta)) 1 else priors$a_eta
  b_eta <- if (is.null(priors$b_eta)) 1 else priors$b_eta
  a_lambda <- if (is.null(priors$a_lambda)) 1 else priors$a_lambda
  b_lambda <- if (is.null(priors$b_lambda)) 1 else priors$b_lambda
  gamma_gnedin <- if (is.null(priors$gamma_gnedin)) 0.8 else priors$gamma_gnedin
  alpha_crp <- if (is.null(priors$alpha_crp)) 1.0 else priors$alpha_crp
  partition_prior <- if (is.null(priors$partition_prior)) "GN" else priors$partition_prior
  
  # init labels
  if (is.null(z_init)) {
    z <- sample(seq_len(K_init), size = n, replace = TRUE)
    z <- match(z, sort(unique(z)))           # relabel to contiguous 1:K_actual
  } else {
    z <- as.integer(z_init)
  }
  
  # init theta
  theta <- rgamma(n, shape = a_eta, rate = b_eta)
  if (normalize_each_iter) theta <- normalize_block_theta(z, theta)
  
  # init lambda according to current K
  K_curr <- length(unique(z))
  lambda <- matrix(
    rgamma(K_curr * K_curr, shape = a_lambda, rate = b_lambda),
    K_curr, K_curr
  )
  
  # storage
  keep_idx <- seq.int(burn_in + 1, iters, by = thin)
  S        <- length(keep_idx)
  z_trace     <- matrix(NA_integer_, S, n)
  theta_trace <- matrix(NA_real_,    S, n)
  lambda_trace <- vector("list", S)
  K_trace <- integer(S)
  
  s <- 0L
  for (t in seq_len(iters)) {
    K_curr <- nrow(lambda)
    
    # 1) update lambda | rest
    lambda <- update_lambda(A, z, theta, a_lambda, b_lambda, K_curr)
    
    # 2) update theta | rest
    theta <- update_theta(A, z, theta, lambda, a_eta, b_eta)
    
    if (normalize_each_iter) theta <- normalize_block_theta(z, theta)
    
    # 3) update z (variable K)
    uz <- update_z_gnedin_varK(
      A, z, theta, lambda,
      gamma_gnedin = gamma_gnedin,
      alpha_crp = alpha_crp,
      partition_prior = partition_prior,
      a_lambda = a_lambda, b_lambda = b_lambda
    )
    z      <- uz$z
    lambda <- uz$lambda
    
    # shrink labels / lambda to used blocks
    sb <- shrink_blocks(z, lambda, theta)
    z      <- sb$z
    lambda <- sb$lambda
    theta  <- sb$theta
    K_curr <- sb$K
    
    if (normalize_each_iter) theta <- normalize_block_theta(z, theta)
    
    # save
    if (t %in% keep_idx) {
      s <- s + 1L
      z_trace[s, ] <- z
      theta_trace[s, ] <- theta
      lambda_trace[[s]] <- lambda
      K_trace[s] <- K_curr
    }
    
    if (is.numeric(verbose) && verbose > 0 && (t %% verbose == 0)) {
      cat(sprintf("iter %d | K=%d", t, K_curr), "\r")
    }
  }
  
  list(
    z         = z_trace,
    theta     = theta_trace,
    lambda    = lambda_trace,
    K         = K_trace,
    keep_idx  = keep_idx
  )
}

# ------------------------------------------------------------
# Compatibility wrapper for older directed call sites
# ------------------------------------------------------------
fit_dirdcsbm_gibbs_gnedin <- function(
    A, K_init,
    iters = 4000, burn_in = floor(iters/2), thin = 10,
    priors = list(
      a_out = 1, b_out = 1,
      a_in = 1, b_in = 1,
      a_lambda = 1, b_lambda = 1,
      gamma_gnedin = 0.8
    ),
    z_init = NULL,
    normalize_each_iter = TRUE,
    save_z = TRUE,
    save_theta = FALSE,
    seed = NULL, verbose = 100
) {
  if (!is.null(priors$a_eta) || !is.null(priors$b_eta)) {
    a_eta <- if (is.null(priors$a_eta)) 1 else priors$a_eta
    b_eta <- if (is.null(priors$b_eta)) 1 else priors$b_eta
  } else {
    a_eta <- if (is.null(priors$a_out)) 1 else priors$a_out
    b_eta <- if (is.null(priors$b_out)) 1 else priors$b_out
  }
  a_lambda <- if (is.null(priors$a_lambda)) 1 else priors$a_lambda
  b_lambda <- if (is.null(priors$b_lambda)) 1 else priors$b_lambda
  gamma_gnedin <- if (is.null(priors$gamma_gnedin)) 0.8 else priors$gamma_gnedin
  alpha_crp <- if (is.null(priors$alpha_crp)) 1.0 else priors$alpha_crp
  partition_prior <- if (is.null(priors$partition_prior)) "GN" else priors$partition_prior
  
  fit_dcsbm_gibbs_gnedin(
    A = A,
    K_init = K_init,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    priors = list(
      a_eta = a_eta, b_eta = b_eta,
      a_lambda = a_lambda, b_lambda = b_lambda,
      gamma_gnedin = gamma_gnedin,
      alpha_crp = alpha_crp,
      partition_prior = partition_prior
    ),
    z_init = z_init,
    normalize_each_iter = normalize_each_iter,
    save_z = save_z,
    save_theta = save_theta,
    seed = seed,
    verbose = verbose
  )
}

# ============================================================
# Simulation (single theta)
# ============================================================
simulate_dcsbm <- function(n, K, z, theta, lambda, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(z) == n, length(theta) == n)
  stopifnot(all(z %in% seq_len(K)))
  
  A <- matrix(0L, n, n)
  diag(A) <- 0L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      k <- z[i]; l <- z[j]
      mu <- theta[i] * theta[j] * lambda[k, l]
      A[i, j] <- rpois(1L, mu)
    }
  }
  list(A = A, z = z, theta = theta, lambda = lambda)
}

ari <- function(z1, z2) {
  z1 <- as.integer(factor(z1))
  z2 <- as.integer(factor(z2))
  n <- length(z1)
  tab <- table(z1, z2)
  a <- rowSums(tab)
  b <- colSums(tab)
  comb2 <- function(x) ifelse(x >= 2, x * (x - 1) / 2, 0)
  sum_comb <- sum(comb2(tab))
  sum_a <- sum(comb2(a))
  sum_b <- sum(comb2(b))
  expected <- sum_a * sum_b / comb2(n)
  max_index <- 0.5 * (sum_a + sum_b)
  denom <- max_index - expected
  if (denom == 0) return(0)
  (sum_comb - expected) / denom
}


# --------------------------------------------------------
# Model summary (updated)
#
# Model (i != j):
#   A_ij ~ Poisson(theta[i] * theta[j] * lambda[z[i], z[j]])
# Priors (shape-rate):
#   theta[i]  ~ Gamma(a_eta, b_eta)
#   lambda[k,l] ~ Gamma(a_lambda, b_lambda)
# z follows Gnedin prior (variable K)
#
# Full conditionals:
#   lambda_{kl} ~ Gamma(a_lambda + R_{kl}, b_lambda + T_{kl}(theta))
#   theta_i ~ Gamma(a_eta + G_i, b_eta + sum_ell (lambda_{k ell}+lambda_{ell k}) S^{-i}_ell)
# where k = z[i],
#   R_{kl} = sum_{i:z_i=k} sum_{j:z_j=l, j!=i} A_ij
#   T_{kl} = sum_{i:z_i=k} sum_{j:z_j=l, j!=i} theta_i theta_j
#          = (sum_{i:z_i=k} theta_i)(sum_{j:z_j=l} theta_j) - 1{k=l} sum_{i:z_i=k} theta_i^2
#   G_i = sum_{j!=i} (A_ij + A_ji)
#   S^{-i}_ell = sum_{j:z_j=ell, j!=i} theta_j
#
# Identifiability (used here): within each block k, enforce sum_{i:z_i=k} theta_i = n_k
# --------------------------------------------------------



# ---------- Sufficient statistics ----------

block_sums <- function(z, theta_out, theta_in, K) {
  S_out <- tapply(theta_out, z, sum); S_out <- as.numeric(S_out)[seq_len(K)]; S_out[is.na(S_out)] <- 0
  S_in  <- tapply(theta_in,  z, sum); S_in  <- as.numeric(S_in)[seq_len(K)]; S_in[is.na(S_in)]  <- 0
  diag_sum <- tapply(theta_out * theta_in, z, sum); diag_sum <- as.numeric(diag_sum)[seq_len(K)]; diag_sum[is.na(diag_sum)] <- 0
  n_k <- tabulate(z, nbins = K)
  list(S_out = S_out, S_in = S_in, diag_sum = diag_sum, n_k = n_k)
}


make_lambda_true <- function(K, diag_level = 2.0, off_range = c(0.03, 0.20)) {
  lam <- matrix(runif(K * K, off_range[1], off_range[2]), K, K)
  diag(lam) <- diag_level + runif(K, 0, 0.4)
  # add a bit of directed asymmetry (otherwise it’s too “nice”)
  if (K >= 2) {
    for (k in 1:(K - 1)) lam[k, k + 1] <- off_range[2] + runif(1, 0, 0.2)
    for (k in 2:K)       lam[k, k - 1] <- off_range[1] + runif(1, 0, 0.05)
  }
  lam
}

demo_dirdcsbm_gnedin_small_study <- function(
    seed = 123,
    K_trues = rep(c(3, 5, 7),5),  # 15 runs
    n = 80,          # keeps n small
    iters = 2500, burn_in = 1250, thin = 1   # minimal-ish but end-to-end
) {
  set.seed(seed)
  
  results <- vector("list", length(K_trues))
  fits    <- vector("list", length(K_trues))
  
  for (r in seq_along(K_trues)) {
    K_true <- K_trues[r]
    set.seed(seed)
    
    z_true <- rep(1:K_true, length = n)
    
    theta_true <- runif(n, 0.8, 1.2)
    theta_true <- normalize_block_theta(z_true, theta_true)
    
    lambda_true <- matrix(0.05, K_true, K_true)
    diag(lambda_true) <- 0.15
    lambda_true[1, 2] <- 0.1;  lambda_true[2, 1] <- 0.5
    lambda_true[2, 3] <- 0.23; lambda_true[3, 2] <- 0.34
    lambda_true[1, 3] <- 0.9;  lambda_true[3, 1] <- 0.22
    
    sim <- simulate_dcsbm(n, K_true, z_true, theta_true, lambda_true)
    A <- sim$A
    
    
    K_init <- nrow(A) # over-specified initial K
    z_init <- sample(rep(seq_len(K_init), length.out = n))
    z_init <- as.integer(factor(z_init))  # ensure labels are 1:K_used
    
    fit <- fit_dcsbm_gibbs_gnedin(
      A, K_init = nrow(A),
      iters = 1500, burn_in = 750, thin = 1,
      priors = list(
        a_eta = 1, b_eta = 1,
        a_lambda = 1, b_lambda = 1,
        gamma_gnedin = 0.8
      ),
      normalize_each_iter = TRUE,
      verbose = 250,
      seed = seed
    )
    
    # VI point estimate from the posterior similarity matrix
    psm   <- comp.psm(fit$z)
    z_hat <- minVI(psm)$cl
    
    results[[r]] <- data.frame(
      K_true = K_true,
      n = n,
      K_init_used = max(z_init),
      ARI = fossil::adj.rand.index(z_true, z_hat),
      vi_dist = vi.dist(z_true, z_hat),
      K_post_mean = mean(fit$K),
      K_post_median = median(fit$K),
      K_post_last = tail(fit$K, 1)
    )
    
    fits[[r]] <- list(
      A = A,
      truth = list(z = z_true, lambda = lambda_true),
      fit = fit,
      z_hat = z_hat
    )
    
    cat(sprintf(
      "Done: K_true=%d (n=%d) | ARI=%.3f | mean K=%.2f\n",
      K_true, n, results[[r]]$ARI, results[[r]]$K_post_mean
    ))
  }
  
  summary_df <- dplyr::bind_rows(results)
  list(summary = summary_df, runs = fits)
}


if(F){
  out <- demo_dirdcsbm_gnedin_small_study()
  out$summary
  View(out$summary)  # uncomment to see full summary table
  write.csv(out$summary,file = "./Desktop/out_summary.csv")
  library(kableExtra)
  
  out_df = read.csv(file = "./Desktop/out_summary.csv")
  out_df%>%
    group_by(K_true)%>%
    summarise(mean_ARI = mean(ARI),
              sd_ARI = sd(ARI),
              mean_K_post_mean = mean(K_post_mean),
              sd_K_post_mean = sd(K_post_mean))%>%
    kable(format = 'latex', booktabs = T, digits = 3)
}
