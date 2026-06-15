#!/usr/bin/env Rscript
# ------------------------------------------------------------
# SST fast-vs-brute-force directional score checks
#   - existing-block scores: dir_exact_SST_existing_counts()
#   - new-slot insertion scores: dir_exact_SST_new_vec_counts()
#
# What this catches immediately:
#   - wrong c_plus / N_tot construction
#   - wrong sign conventions
#   - wrong distance indexing / slot insertion logic
#   - (bonus) whether skipping within-block terms is biasing comparisons
# ------------------------------------------------------------

# ---- 0) Source your fast implementations (adjust path as needed) ------------
# This file must define:
#   logsig(), log1m_sig(), .log2cosh_half()
#   dir_exact_SST_existing_counts()
#   dir_exact_SST_new_vec_counts()
#   log_int_extreme_binom_truncnorm()


inv_logit <- function(x) 1 / (1 + exp(-x))

# ---- 1) Tiny SST data generator (optional, but very handy) ------------------
simulate_sst_toy <- function(n = 12L, K = 4L, lambda_N = 3,
                             psi_vec = NULL, seed = 1L) {
  set.seed(seed)
  
  stopifnot(K >= 1L, n >= K)
  
  # make a nonempty partition
  z <- rep(seq_len(K), length.out = n)
  z <- sample(z, n, replace = FALSE)
  
  # monotone nonnegative psi
  if (is.null(psi_vec)) {
    if (K <= 1L) {
      psi_vec <- numeric(0)
    } else {
      inc <- abs(rnorm(K - 1L, mean = 0.6, sd = 0.2))
      psi_vec <- cummax(cumsum(inc))
    }
  }
  stopifnot(length(psi_vec) == max(K - 1L, 0L))
  
  N <- matrix(0L, n, n)
  A <- matrix(0L, n, n)
  
  for (i in 1:(n - 1L)) for (j in (i + 1L):n) {
    Nij <- rpois(1L, lambda = lambda_N)
    if (Nij == 0L) next
    
    if (z[i] == z[j]) {
      p_ij <- 0.5
    } else {
      d <- abs(z[i] - z[j])
      th <- psi_vec[d]
      # lower label beats higher label with prob logistic(+psi_d)
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
  
  list(A = A, N = N, z = z, psi = psi_vec)
}

# ---- 2) Block-aggregated counts used by your fast score functions -----------
# c_plus[ell] = total wins of i against opponents currently in block ell
# N_tot[ell]  = total matches of i against block ell
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

# ---- 3) Brute force: existing-block score (dyad-level) ----------------------
# Two versions:
#  - kernel: excludes within-block dyads (ell == k) so it should match FAST
#  - full: includes within-block dyads at theta = 0 (i.e. p=0.5 giving -n log 2)
brute_existing_kernel <- function(i, k, A, N, z, psi_vec) {
  n <- length(z)
  val <- 0
  
  for (j in seq_len(n)) {
    if (j == i) next
    ell <- z[j]
    if (ell == k) next               # match your FAST function behaviour
    Nij <- N[i, j]
    if (Nij == 0) next
    
    d  <- abs(ell - k)
    th <- psi_vec[d]
    sgn <- if (ell > k) +1 else -1   # if i lower -> +th, else -th
    theta <- sgn * th
    Aij <- A[i, j]
    
    val <- val + Aij * logsig(theta) + (Nij - Aij) * log1m_sig(theta)
  }
  val
}

brute_existing_full <- function(i, k, A, N, z, psi_vec) {
  n <- length(z)
  K = length(unique(z))

  val <-0
  
  for (j in seq_len(n)) {
    if (j == i) next
    ell <- z[j]
    Nij <- N[i, j]
    if (Nij == 0) next
    
    if (ell == k) {
      theta <- 0
    } else {
      d  <- abs(ell - k)
      th <- psi_vec[d]
      sgn <- if (ell > k) +1 else -1
      theta <- sgn * th
    }
    
    Aij <- A[i, j]
    val <- val + Aij * logsig(theta) + (Nij - Aij) * log1m_sig(theta)
  }
  val
}

# ---- 4) Brute force: new-slot insertion score (dyad-level) ------------------
# This matches your fast new-slot logic, including the extreme-distance integral.
brute_newslot <- function(i, r, A, N, z, psi_vec, tau0 = 0.15, m_inc = 0) {
  n <- length(z)
  K_old <- max(z)
  if (K_old <= 1L) return(0)
  
  val <- 0
  A_K <- 0
  n_K <- 0
  
  for (j in seq_len(n)) {
    if (j == i) next
    ell <- z[j]
    Nij <- N[i, j]
    if (Nij == 0) next
    
    # after inserting new singleton at slot r, old block ell shifts if ell >= r
    pos_ell <- if (ell < r) ell else ell + 1L
    d <- abs(pos_ell - r)
    if (d == 0L) next
    
    Aij <- A[i, j]
    
    if (d <= (K_old - 1L)) {
      th <- psi_vec[d]
      theta <- if (pos_ell > r) +th else -th
      val <- val + Aij * logsig(theta) + (Nij - Aij) * log1m_sig(theta)
      
    } else if (d == K_old) {
      # aggregate forward wins (lower label -> higher label)
      A_fwd <- if (pos_ell > r) Aij else (Nij - Aij)
      A_K <- A_K + A_fwd
      n_K <- n_K + Nij
      
    } else {
      stop("Impossible distance encountered in brute_newslot().")
    }
  }
  
  if (n_K > 0) {
    psiKm1 <- psi_vec[K_old - 1L]
    val <- val + log_int_extreme_binom_truncnorm(
      A = A_K, n = n_K, psiKm1 = psiKm1, tau0 = tau0, m_inc = m_inc
    )
  }
  
  val
}

# ---- 5) Runner: compares FAST vs BRUTE and prints diagnostics ---------------
run_fast_bruteforce_tests <- function(A, N, z, psi_vec,
                                      i_set = NULL,
                                      tau0 = 0.15, m_inc = 0,
                                      tol = 1e-8, verbose = TRUE) {
  n <- length(z)
  K <- max(z)
  
  stopifnot(is.matrix(A), is.matrix(N), nrow(A) == n, ncol(A) == n)
  stopifnot(all(dim(N) == dim(A)))
  stopifnot(length(psi_vec) == max(K - 1L, 0L))
  
  if (is.null(i_set)) {
    # pick 3-ish nodes, avoiding singleton blocks where possible
    tab <- table(z)
    ok  <- which(tab[z] >= 2L)
    i_set <- if (length(ok)) sample(ok, min(3L, length(ok))) else sample.int(n, min(3L, n))
  }
  
  out <- list()
  
  for (i in i_set) {
    cnt <- counts_i_by_block(i, A, N, z)
    c_plus <- cnt$c_plus
    N_tot  <- cnt$N_tot
    
    # --- existing-block scores
    # --- existing-block scores
    fast_exist <- vapply(seq_len(K),
                         function(k) dir_exact_SST_existing_counts(k, c_plus, N_tot, psi_vec),
                         numeric(1))
    
    brute_exist_kernel <- vapply(seq_len(K),
                                 function(k) brute_existing_kernel(i, k, A, N, z, psi_vec),
                                 numeric(1))
    
    brute_exist_full <- vapply(seq_len(K),
                               function(k) brute_existing_full(i, k, A, N, z, psi_vec),
                               numeric(1))
    
    # THIS is the correct equality check now:
    diff_exist <- fast_exist - brute_exist_full
    max_abs_exist <- max(abs(diff_exist))
    
    # optional diagnostic: what kernel misses
    withinN <- vapply(seq_len(K), function(k) {
      js <- which(z == k & seq_len(n) != i)
      sum(N[i, js])
    }, numeric(1))
    
    diff_kernel <- fast_exist - brute_exist_kernel
    expected_kernel <- -withinN * log(2)
    
    
    # --- new-slot scores
    fast_new <- dir_exact_SST_new_vec_counts(c_plus, N_tot, psi_vec, tau0 = tau0, m_inc = m_inc)
    brute_new <- vapply(seq_len(K + 1L),
                        function(r) brute_newslot(i, r, A, N, z, psi_vec, tau0 = tau0, m_inc = m_inc),
                        numeric(1))
    diff_new <- fast_new - brute_new
    max_abs_new <- max(abs(diff_new))
    
    if (verbose) {
      cat("\n============================================================\n")
      cat(sprintf("Node i = %d (current block z[i]=%d), K = %d\n", i, z[i], K))
      cat("Existing-block check (FAST - BRUTE_FULL):\n")
      print(data.frame(
        k = seq_len(K),
        fast = fast_exist,
        brute_full = brute_exist_full,
        diff = diff_exist
      ), row.names = FALSE)
      
      cat(sprintf("Max |diff| existing = %.3e  (tol = %.1e)\n", max_abs_exist, tol))
      
      cat("\nKernel diagnostic (FAST - BRUTE_KERNEL) should equal -withinN*log(2):\n")
      print(data.frame(
        k = seq_len(K),
        diff_kernel = diff_kernel,
        expected = expected_kernel,
        diff_minus_expected = diff_kernel - expected_kernel
      ), row.names = FALSE)
      
      
      cat("\nNew-slot insertion check (FAST - BRUTE):\n")
      print(data.frame(
        r = seq_len(K + 1L),
        fast = fast_new,
        brute = brute_new,
        diff = diff_new
      ), row.names = FALSE)
      
      cat(sprintf("Max |diff| new-slot  = %.3e  (tol = %.1e)\n", max_abs_new, tol))
      
      if (max_abs_exist > tol || max_abs_new > tol) {
        cat("\n⚠️  Mismatch above tolerance. That’s not ‘Monte Carlo noise’, it’s algebra or indexing.\n")
      } else {
        cat("\n✅  FAST matches BRUTE within tolerance.\n")
      }
    }
    
    out[[as.character(i)]] <- list(
      i = i,
      fast_exist = fast_exist,
      brute_exist_kernel = brute_exist_kernel,
      brute_exist_full = brute_exist_full,
      max_abs_exist = max_abs_exist,
      fast_new = fast_new,
      brute_new = brute_new,
      max_abs_new = max_abs_new
    )
  }
  
  invisible(out)
}

# ---- 6) Example run on a toy instance --------------------------------------
# Comment this block out if you only want the function definitions.
if (sys.nframe() == 0L) {
  toy <- simulate_sst_toy(n = 14L, K = 4L, lambda_N = 3, seed = 42L)
  cat("Toy psi:", paste(round(toy$psi, 3), collapse = ", "), "\n")
  run_fast_bruteforce_tests(
    A = toy$A, N = toy$N, z = toy$z, psi_vec = toy$psi,
    tau0 = 0.15, m_inc = 0,
    tol = 1e-8, verbose = TRUE
  )
}
