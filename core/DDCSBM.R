# Directed Degree-Corrected SBM (Poisson) — Gibbs Sampler
# --------------------------------------------------------
# Implements the baseline model in Section 1 of the note.
# 
#
# Model (for i != j):
#   A_ij ~ Poisson(theta_out[i] * theta_in[j] * lambda[z[i], z[j]])
# Priors:
#   theta_out[i] ~ Gamma(a_out, b_out)  (shape-rate)
#   theta_in[i]  ~ Gamma(a_in,  b_in)
#   lambda[k,l]   ~ Gamma(a_lambda, b_lambda)
#   pi           ~ Dirichlet(alpha)
# z[i] | pi ~ Categorical(pi)
#
# Full conditionals (shape-rate):
#   lambda_{kl} ~ Gamma(a_lambda + R_{kl}, b_lambda + T_{kl}(theta))
#   theta_out[i] ~ Gamma(a_out + G_out[i], b_out + sum_l lambda_{k l}(S_in[l] - 1{l=k} * theta_in[i]))
#   theta_in[i]  ~ Gamma(a_in  + G_in[i],  b_in  + sum_k lambda_{k l}(S_out[k] - 1{k=l} * theta_out[i]))
# where k = z[i], l = z[i],
#   R_{kl} = sum_{i:z[i]=k} sum_{j:z[j]=l, j!=i} A_ij
#   T_{kl}(theta) = sum_{i:z[i]=k} sum_{j:z[j]=l, j!=i} theta_out[i] theta_in[j]
#                 = (sum_{i:z[i]=k} theta_out[i]) (sum_{j:z[j]=l} theta_in[j]) - 1{k=l} * sum_{i:z[i]=k} theta_out[i] theta_in[i]
#
# Identifiability: for each block k set sum_{i:z[i]=k} theta_out[i] = n_k and sum_{i:z[i]=k} theta_in[i] = n_k
# (performed after theta updates; skip empty blocks.)
# --------------------------------------------------------

# ---------- Utilities ----------

sample_dirichlet <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}

safe_log <- function(x) {
  # avoid -Inf for zero
  log(pmax(x, .Machine$double.xmin))
}

row_mode <- function(x) {
  # mode per row for integer matrix
  apply(x, 1, function(v) {
    t <- table(v)
    as.integer(names(t)[which.max(t)])
  })
}

# Adjusted Rand Index (ARI) — self-contained implementation
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

# ---------- Simulation ----------

simulate_dirdcsbm <- function(n, K, z, theta_out, theta_in, lambda, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(z) == n, length(theta_out) == n, length(theta_in) == n)
  stopifnot(all(z %in% seq_len(K)))
  A <- matrix(0L, n, n)
  diag(A) <- 0L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      k <- z[i]; l <- z[j]
      lam <- theta_out[i] * theta_in[j] * lambda[k, l]
      A[i, j] <- rpois(1L, lam)
    }
  }
  list(A = A, z = z, theta_out = theta_out, theta_in = theta_in, lambda = lambda)
}

# ---------- Sufficient statistics ----------

block_sums <- function(z, theta_out, theta_in, K) {
  S_out <- tapply(theta_out, z, sum); S_out <- as.numeric(S_out)[seq_len(K)]; S_out[is.na(S_out)] <- 0
  S_in  <- tapply(theta_in,  z, sum); S_in  <- as.numeric(S_in)[seq_len(K)]; S_in[is.na(S_in)]  <- 0
  diag_sum <- tapply(theta_out * theta_in, z, sum); diag_sum <- as.numeric(diag_sum)[seq_len(K)]; diag_sum[is.na(diag_sum)] <- 0
  n_k <- tabulate(z, nbins = K)
  list(S_out = S_out, S_in = S_in, diag_sum = diag_sum, n_k = n_k)
}


R_counts <- function(A, z, K) {
  n <- nrow(A)
  R <- matrix(0, K, K)
  for (k in seq_len(K)) {
    I <- which(z == k)
    if (length(I) == 0) next
    for (l in seq_len(K)) {
      J <- which(z == l)
      if (length(J) == 0) next
      if (k != l) {
        R[k, l] <- sum(A[I, J, drop = FALSE])
      } else {
        # exclude self-edges i=j
        R[k, l] <- sum(A[I, I, drop = FALSE]) - sum(diag(A[I, I, drop = FALSE]))
      }
    }
  }
  R
}

T_exposures <- function(z, theta_out, theta_in, K) {
  bs <- block_sums(z, theta_out, theta_in, K)
  S_out <- bs$S_out; S_in <- bs$S_in; diag_sum <- bs$diag_sum
  T <- outer(S_out, S_in, "*")
  diag(T) <- diag(T) - diag_sum
  T[!is.finite(T)] <- 0
  T <- pmax(T, 0)              # numeric guard
  T
}

# ---------- Parameter updates ----------

update_lambda <- function(A, z, theta_out, theta_in, a_lambda, b_lambda, K) {
  R <- R_counts(A, z, K)
  T <- T_exposures(z, theta_out, theta_in, K)
  shape <- a_lambda + R
  rate  <- b_lambda + T
  matrix(rgamma(K * K, shape = as.vector(shape), rate = as.vector(rate)), K, K)
}

.safe_rate <- function(x) {
  if (!is.finite(x) || x <= 0) return(1e-8)
  x
}
.safe_shape <- function(x) {
  if (!is.finite(x) || x <= 0) return(1e-8)
  x
}

update_theta_out <- function(A, z, theta_out, theta_in, lambda, a_out, b_out) {
  n <- nrow(A); K <- nrow(lambda)  # <- FIX
  G_out <- rowSums(A)
  S_in <- tapply(theta_in, z, sum)
  S_in <- as.numeric(S_in)[seq_len(K)]; S_in[is.na(S_in)] <- 0
  
  for (i in seq_len(n)) {
    k <- z[i]
    vec <- S_in - (seq_len(K) == k) * ifelse(is.finite(theta_in[i]), theta_in[i], 0)
    rate_i <- b_out + sum(lambda[k, ] * vec, na.rm = TRUE)   # lengths now match
    rate_i <- .safe_rate(rate_i)
    shape_i <- .safe_shape(a_out + G_out[i])
    val <- rgamma(1L, shape = shape_i, rate = rate_i)
    theta_out[i] <- if (is.finite(val) && val > 0) val else 1e-6
  }
  theta_out
}

update_theta_in <- function(A, z, theta_out, theta_in, lambda, a_in, b_in) {
  n <- nrow(A); K <- nrow(lambda)  # <- FIX
  G_in <- colSums(A)
  S_out <- tapply(theta_out, z, sum)
  S_out <- as.numeric(S_out)[seq_len(K)]; S_out[is.na(S_out)] <- 0
  
  for (i in seq_len(n)) {
    l <- z[i]
    vec <- S_out - (seq_len(K) == l) * ifelse(is.finite(theta_out[i]), theta_out[i], 0)
    rate_i <- b_in + sum(lambda[, l] * vec, na.rm = TRUE)    # lengths now match
    rate_i <- .safe_rate(rate_i)
    shape_i <- .safe_shape(a_in + G_in[i])
    val <- rgamma(1L, shape = shape_i, rate = rate_i)
    theta_in[i] <- if (is.finite(val) && val > 0) val else 1e-6
  }
  theta_in
}



normalize_block_thetas <- function(z, theta_out, theta_in) {
  K <- max(z, na.rm = TRUE)
  for (k in seq_len(K)) {
    idx <- which(z == k)
    if (!length(idx)) next
    so <- sum(theta_out[idx]); si <- sum(theta_in[idx]); n_k <- length(idx)
    if (is.finite(so) && so > 0) theta_out[idx] <- theta_out[idx] * (n_k / so)
    if (is.finite(si) && si > 0) theta_in[idx]  <- theta_in[idx]  * (n_k / si)
  }
  list(theta_out = theta_out, theta_in = theta_in)
}

update_z <- function(A, z, theta_out, theta_in, lambda, log_pi) {
  n <- nrow(A); K <- nrow(lambda)
  for (i in sample.int(n)) {
    Ai_row <- A[i, ]
    Ai_col <- A[, i]
    
    logp <- rep(-Inf, K)
    for (k in seq_len(K)) {
      # guard: if any NA creeps in, skip contribution
      term1 <- suppressWarnings(
        sum(Ai_row * (safe_log(theta_out[i]) + safe_log(theta_in) + safe_log(lambda[k, z])) -
              theta_out[i] * theta_in * lambda[k, z], na.rm = TRUE)
      )
      term1 <- term1 - (Ai_row[i] * (safe_log(theta_out[i]) + safe_log(theta_in[i]) + safe_log(lambda[k, z[i]])) -
                          theta_out[i] * theta_in[i] * lambda[k, z[i]])
      
      term2 <- suppressWarnings(
        sum(Ai_col * (safe_log(theta_out) + safe_log(theta_in[i]) + safe_log(lambda[z, k])) -
              theta_out * theta_in[i] * lambda[z, k], na.rm = TRUE)
      )
      term2 <- term2 - (Ai_col[i] * (safe_log(theta_out[i]) + safe_log(theta_in[i]) + safe_log(lambda[z[i], k])) -
                          theta_out[i] * theta_in[i] * lambda[z[i], k])
      
      logp[k] <- log_pi[k] + term1 + term2
    }
    
    # stabilize
    if (!all(is.finite(logp))) {
      # replace non-finite with a very small number
      logp[!is.finite(logp)] <- -1e300
    }
    m <- max(logp)
    if (!is.finite(m)) {
      # degenerate: keep current label
      next
    }
    p <- exp(logp - m)
    s <- sum(p)
    if (s <= 0 || !is.finite(s)) {
      next  # keep current label
    }
    z[i] <- sample.int(K, 1L, prob = p / s)
  }
  z
}
# Global identifiability: sum_i theta_out[i] = target_sum_out, sum_i theta_in[i] = target_sum_in
normalize_global_thetas <- function(theta_out, theta_in,
                                    target_sum_out = length(theta_out),
                                    target_sum_in  = length(theta_in)) {
  s_out <- sum(theta_out)
  s_in  <- sum(theta_in)
  if (is.finite(s_out) && s_out > 0) theta_out <- theta_out * (target_sum_out / s_out)
  if (is.finite(s_in)  && s_in  > 0) theta_in  <- theta_in  * (target_sum_in  / s_in)
  list(theta_out = theta_out, theta_in = theta_in)
}

# ---------- Main Gibbs sampler ----------

fit_dirdcsbm_gibbs <- function(
    A, K,
    iters = 4000, burn_in = floodraws_kapr(iters/2), thin = 10,
    priors = list(a_out = 1, b_out = 1, a_in = 1, b_in = 1, a_lambda = 1, b_lambda = 1, alpha = NULL),
    z_init = NULL,
    normalize_each_iter = TRUE,
    save_z = TRUE, save_lambda = TRUE,
    seed = NULL, verbose = 100
) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(A)
  stopifnot(ncol(A) == n)
  diag(A) <- 0
  G_out <- rowSums(A); G_in <- colSums(A) # data-only
  
  if (is.null(priors$alpha)) priors$alpha <- rep(1, K)
  with(priors, {
    a_out <<- a_out; b_out <<- b_out; a_in <<- a_in; b_in <<- b_in; a_lambda <<- a_lambda; b_lambda <<- b_lambda
  })
  
  # init labels
  if (is.null(z_init)) {
    # quick heuristic: k-means on (out,in) degrees
    X <- cbind(G_out + 1e-6, G_in + 1e-6)
    km <- kmeans(log(X), centers = K, nstart = 10)
    z <- as.integer(km$cluster)
  } else {
    z <- as.integer(z_init)
  }
  
  # init thetas (Gamma prior means)
  theta_out <- rgamma(n, shape = a_out, rate = b_out)
  theta_in  <- rgamma(n, shape = a_in,  rate = b_in)
  tmp       <- normalize_block_thetas(z, theta_out, theta_in)
  theta_out <- tmp$theta_out; theta_in <- tmp$theta_in
  
  # init lambda
  lambda <- matrix(rgamma(K*K, shape = a_lambda, rate = b_lambda), K, K)
  
  # init pi
  pi <- sample_dirichlet(priors$alpha)
  
  # storage
  keep_idx         <- seq.int(burn_in + 1, iters, by = thin)
  S                <- length(keep_idx)
  z_trace          <- matrix(NA_integer_, S, n) 
  lambda_trace     <- array(NA_real_, dim = c(S,K, K)) 
  theta_out_trace  <-  matrix(NA_real_, S, n) 
  theta_in_trace   <- matrix(NA_real_, S, n) 
  
  s <- 0L
  for (t in seq_len(iters)) {
    # 1) update lambda | rest
    lambda <- update_lambda(A, z, theta_out, theta_in, a_lambda, b_lambda, K)
    
    # 2) update theta_out, theta_in | rest
    theta_out <- update_theta_out(A, z, theta_out, theta_in, lambda, a_out, b_out)
    theta_in  <- update_theta_in(A, z, theta_out, theta_in, lambda, a_in,  b_in)
    
    if (normalize_each_iter) {
      tmp <- normalize_block_thetas(z, theta_out, theta_in)
      theta_out <- tmp$theta_out; theta_in <- tmp$theta_in
    }
    
    # 3) update z | rest
    log_pi <- safe_log(pi)
    z <- update_z(A, z, theta_out, theta_in, lambda, log_pi)
    
    if (normalize_each_iter) {
      tmp <- normalize_block_thetas(z, theta_out, theta_in)
      theta_out <- tmp$theta_out; theta_in <- tmp$theta_in
    }
    
    # 4) update pi | z
    n_k <- tabulate(z, nbins = K)
    pi <- sample_dirichlet(priors$alpha + n_k)
    
    # save
    if (t %in% keep_idx) {
      s <- s + 1L
      z_trace[s,] <- z
      lambda_trace[s, , ] <- lambda
      theta_out_trace[s,] <- theta_out
      theta_in_trace[s,]  <- theta_in
    }
    
    if (!isFALSE(verbose) && (t %% verbose == 0)) {
      cat(sprintf("iter %d | min(n_k)=%d max(n_k)=%d\r", t, min(tabulate(z, K)), max(tabulate(z, K))))
    }
  }
  
  list(z_last = z,
       theta_out_last = theta_out,
       theta_in_last = theta_in,
       lambda_last = lambda,
       pi_last = pi,
       z = z_trace,
       lambda = lambda_trace,
       theta_out = theta_out_trace,
       theta_in  = theta_in_trace,
       keep_idx = keep_idx)
}

# ---------- Simple demo ----------
# A small self-contained demo that (i) simulates a network; (ii) fits the model; (iii) reports ARI.

demo_dirdcsbm <- function(seed = 123) {
  set.seed(seed)
  n <- 60; K <- 3
  # true labels (balanced)
  z_true <- rep(1:K, each = n/K)
  # node effects
  theta_out_true <- runif(n,0.8,1.2)
  theta_in_true  <- runif(n,0.8,1.2)
  # within-block normalization for identifiability in the DGP (optional)
  tmp <- normalize_block_thetas(z_true, theta_out_true, theta_in_true)
  theta_out_true <- tmp$theta_out; theta_in_true <- tmp$theta_in
  # block intensities: strong diagonal, weak off-diagonal, asymmetric allowed
  lambda_true <- matrix(0.05, K, K)
  diag(lambda_true) <- 2.5
  # add some asymmetry
  lambda_true[1, 2] <- 0.1; lambda_true[2, 1] <- 0.5
  lambda_true[2, 3] <- 0.23; lambda_true[3, 2] <- 0.34
  lambda_true[1, 3] <- 0.9; lambda_true[3, 1] <- 0.22
  
  sim <- simulate_dirdcsbm(n, K, z_true, theta_out_true, theta_in_true, lambda_true)
  A <- sim$A
  
  fit <- fit_dirdcsbm_gibbs(
    A, K,
    iters = 1500, burn_in = 750, thin = 1,
    priors = list(a_out = 1, b_out = 1, a_in = 1, b_in = 1, a_lambda = 1, b_lambda = 1, alpha = rep(1, K)),
    normalize_each_iter = TRUE,
    verbose = 250,
    seed = seed
  )
  
  psm = comp.psm(fit$z_trace)
  
  z_hat <- minVI(psm)$cl
  
  cat(sprintf("ARI (last draw vs truth): %.3f\n", ari(z_true, z_hat)))
  
  if (!is.null(fit$lambda_trace)) {
    lambda_mean <- apply(fit$lambda_trace, c(2, 3), mean, na.rm = TRUE)
    cat("\nPosterior mean lambda (rounded):\n")
    print(round(lambda_mean, 2))
  }
  
  invisible(list(A = A, truth = list(z = z_true, theta_out = theta_out_true, theta_in = theta_in_true, lambda = lambda_true), fit = fit))
}

# To run the demo in an R session:
#   source("directed_dcsbm_gibbs.R")


 
 
 
 
 