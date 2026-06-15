#!/usr/bin/env Rscript
# test_fast_bruteforce_WST.R
# ------------------------------------------------------------
# WST fast-vs-brute-force directional score checks
# ------------------------------------------------------------

# stable log(1+exp(x))
log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
logsig   <- function(x) -log1pexp(-x)
log1m_sig<- function(x) -log1pexp(x)

# stable log(2 cosh(x/2))
.log2cosh_half <- function(x) {
  u <- abs(x) / 2
  u + log1p(exp(-2*u))
}




# ---------------------------
# BRUTE: existing-block WST (dyad-level)
# ---------------------------
brute_existing_wst_kernel <- function(i, k, A, N, z, psi_mat) {
  n <- length(z)
  val <- 0
  
  for (j in seq_len(n)) {
    if (j == i) next
    ell <- z[j]
    if (ell == k) next
    
    Nij <- N[i, j]
    if (Nij == 0) next
    
    th <- psi_mat[min(k, ell), max(k, ell)]
    theta <- if (ell > k) +th else -th
    
    Aij <- A[i, j]
    val <- val + Aij * logsig(theta) + (Nij - Aij) * log1m_sig(theta)
  }
  
  val
}

brute_existing_wst_full <- function(i, k, A, N, z, psi_mat) {
  n <- length(z)
  val <- 0
  
  for (j in seq_len(n)) {
    if (j == i) next
    ell <- z[j]
    
    Nij <- N[i, j]
    if (Nij == 0) next
    
    if (ell == k) {
      theta <- 0
    } else {
      th <- psi_mat[min(k, ell), max(k, ell)]
      theta <- if (ell > k) +th else -th
    }
    
    Aij <- A[i, j]
    val <- val + Aij * logsig(theta) + (Nij - Aij) * log1m_sig(theta)
  }
  
  val
}

# ---------------------------
# BRUTE: new-slot WST
# ---------------------------
# Brute loops dyads, aggregates (A_fwd, n) per block ell, then integrates per ell
brute_newslot_wst <- function(i, r, A, N, z, mu0 = 0, sigma0 = 1, rel.tol = 1e-10) {
  n <- length(z)
  K_old <- max(z)
  if (K_old <= 1L) return(0)
  
  A_fwd_blk <- numeric(K_old)
  n_blk     <- numeric(K_old)
  
  for (j in seq_len(n)) {
    if (j == i) next
    ell <- z[j]
    Nij <- N[i, j]
    if (Nij == 0) next
    
    pos_ell <- if (ell < r) ell else ell + 1L
    if (pos_ell == r) next
    
    Aij <- A[i, j]
    A_fwd <- if (pos_ell > r) Aij else (Nij - Aij)
    
    A_fwd_blk[ell] <- A_fwd_blk[ell] + A_fwd
    n_blk[ell]     <- n_blk[ell] + Nij
  }
  
  val <- 0
  for (ell in seq_len(K_old)) {
    if (n_blk[ell] == 0) next
    val <- val + log_int_binom_truncnorm(A = A_fwd_blk[ell], n = n_blk[ell],
                                         mu0 = mu0, sigma0 = sigma0, rel.tol = rel.tol)
  }
  
  val
}


inv_logit <- function(x) 1 / (1 + exp(-x))

simulate_wst_toy <- function(n = 14L, K = 4L, lambda_N = 3,
                             psi_mat = NULL, seed = 42L) {
  set.seed(seed)
  stopifnot(K >= 1L, n >= K)
  
  # nonempty partition
  z <- rep(seq_len(K), length.out = n)
  z <- sample(z, n, replace = FALSE)
  
  # psi matrix (symmetric, diag 0), positive above diagonal
  if (is.null(psi_mat)) {
    psi_mat <- matrix(0, K, K)
    if (K >= 2L) {
      for (k in 1:(K-1L)) for (ell in (k+1L):K) {
        psi_mat[k, ell] <- rgamma(1, shape = 2, rate = 2)  # mild hierarchy
        psi_mat[ell, k] <- psi_mat[k, ell]
      }
    }
  }
  diag(psi_mat) <- 0
  
  N <- matrix(0L, n, n)
  A <- matrix(0L, n, n)
  
  for (i in 1:(n - 1L)) for (j in (i + 1L):n) {
    Nij <- rpois(1L, lambda = lambda_N)
    if (Nij == 0L) next
    
    if (z[i] == z[j]) {
      p_ij <- 0.5
    } else {
      k <- min(z[i], z[j]); ell <- max(z[i], z[j])
      th <- psi_mat[k, ell]
      # lower label beats higher label with prob logistic(+psi)
      if (z[i] < z[j]) p_ij <- inv_logit(+th) else p_ij <- inv_logit(-th)
    }
    
    Aij <- rbinom(1L, size = Nij, prob = p_ij)
    A[i, j] <- Aij
    A[j, i] <- Nij - Aij
    N[i, j] <- Nij
    N[j, i] <- Nij
  }
  
  diag(A) <- 0L
  diag(N) <- 0L
  
  list(A = A, N = N, z = z, psi_mat = psi_mat)
}

counts_i_by_block <- function(i, A, N, z) {
  n <- length(z)
  K <- max(z)
  c_plus <- numeric(K)
  N_tot  <- numeric(K)
  
  for (ell in seq_len(K)) {
    js <- which(z == ell & seq_len(n) != i)
    if (length(js)) {
      N_tot[ell]  <- sum(N[i, js])
      c_plus[ell] <- sum(A[i, js])
    }
  }
  list(c_plus = c_plus, N_tot = N_tot)
}

run_fast_bruteforce_tests_wst <- function(A, N, z, psi_mat,
                                          i_set = NULL,
                                          mu0 = 0, sigma0 = 1,
                                          tol = 1e-8, verbose = TRUE) {
  n <- length(z)
  K <- max(z)
  stopifnot(is.matrix(A), is.matrix(N), all(dim(A) == dim(N)))
  stopifnot(is.matrix(psi_mat), nrow(psi_mat) == K, ncol(psi_mat) == K)
  
  if (is.null(i_set)) {
    tab <- table(z)
    ok  <- which(tab[z] >= 2L)
    i_set <- if (length(ok)) sample(ok, min(3L, length(ok))) else sample.int(n, min(3L, n))
  }
  
  out <- list()
  
  for (i in i_set) {
    cnt <- counts_i_by_block(i, A, N, z)
    c_plus <- cnt$c_plus
    N_tot  <- cnt$N_tot
    
    # --- existing block: fast vs brute
    fast_exist_kernel <- vapply(seq_len(K),
                                function(k) dir_exact_WST_existing_counts(k, c_plus, N_tot, psi_mat, include_within = FALSE),
                                numeric(1)
    )
    brute_exist_kernel <- vapply(seq_len(K),
                                 function(k) brute_existing_wst_kernel(i, k, A, N, z, psi_mat),
                                 numeric(1)
    )
    
    fast_exist_full <- vapply(seq_len(K),
                              function(k) dir_exact_WST_existing_counts(k, c_plus, N_tot, psi_mat, include_within = TRUE),
                              numeric(1)
    )
    brute_exist_full <- vapply(seq_len(K),
                               function(k) brute_existing_wst_full(i, k, A, N, z, psi_mat),
                               numeric(1)
    )
    
    diff_exist_kernel <- fast_exist_kernel - brute_exist_kernel
    diff_exist_full   <- fast_exist_full   - brute_exist_full
    
    max_abs_exist_kernel <- max(abs(diff_exist_kernel))
    max_abs_exist_full   <- max(abs(diff_exist_full))
    
    # --- new slot: fast vs brute
    fast_new <- dir_exact_WST_new_vec_counts(c_plus, N_tot, mu0 = mu0, sigma0 = sigma0)
    brute_new <- vapply(seq_len(K + 1L),
                        function(r) brute_newslot_wst(i, r, A, N, z, mu0 = mu0, sigma0 = sigma0),
                        numeric(1)
    )
    diff_new <- fast_new - brute_new
    max_abs_new <- max(abs(diff_new))
    
    if (verbose) {
      cat("\n============================================================\n")
      cat(sprintf("Node i = %d (current block z[i]=%d), K = %d\n", i, z[i], K))
      
      cat("\nExisting-block check (KERNEL, FAST - BRUTE):\n")
      print(data.frame(k = seq_len(K), fast = fast_exist_kernel,
                       brute = brute_exist_kernel, diff = diff_exist_kernel),
            row.names = FALSE)
      cat(sprintf("Max |diff| existing (kernel) = %.3e  (tol=%.1e)\n", max_abs_exist_kernel, tol))
      
      cat("\nExisting-block check (FULL, FAST - BRUTE):\n")
      print(data.frame(k = seq_len(K), fast = fast_exist_full,
                       brute = brute_exist_full, diff = diff_exist_full),
            row.names = FALSE)
      cat(sprintf("Max |diff| existing (full)   = %.3e  (tol=%.1e)\n", max_abs_exist_full, tol))
      
      cat("\nNew-slot insertion check (FAST - BRUTE):\n")
      print(data.frame(r = seq_len(K + 1L), fast = fast_new, brute = brute_new, diff = diff_new),
            row.names = FALSE)
      cat(sprintf("Max |diff| new-slot          = %.3e  (tol=%.1e)\n", max_abs_new, tol))
      
      if (max_abs_exist_kernel > tol || max_abs_exist_full > tol || max_abs_new > tol) {
        cat("\n⚠️  Mismatch above tolerance. That’s indexing/sign/integral logic, not randomness.\n")
      } else {
        cat("\n✅  FAST matches BRUTE within tolerance.\n")
      }
    }
    
    out[[as.character(i)]] <- list(
      i = i,
      max_abs_exist_kernel = max_abs_exist_kernel,
      max_abs_exist_full   = max_abs_exist_full,
      max_abs_new          = max_abs_new
    )
  }
  
  invisible(out)
}

# ---- Example run on toy instance ----
if (sys.nframe() == 0L) {
  toy <- simulate_wst_toy(n = 14L, K = 4L, lambda_N = 3, seed = 42L)
  cat("Toy psi upper triangle:\n")
  print(round(toy$psi_mat, 3))
  
  run_fast_bruteforce_tests_wst(
    A = toy$A, N = toy$N, z = toy$z, psi_mat = toy$psi_mat,
    mu0 = 0, sigma0 = 1,
    tol = 1e-8, verbose = TRUE
  )
}
