#!/usr/bin/env Rscript

# Slow-reference checks for point-wise log-likelihood matrices:
#   - OSBM (WST and SST)
#   - DC-SBM
#
# Purpose:
#   Validate fast vectorized implementations against transparent dyad-by-dyad loops.

suppressPackageStartupMessages({
  library(Matrix)
})

source("./helper_folder/helper.R")
source("./helper_folder/transitivity_check_helper.R")

build_dyad_index_local <- function(n) {
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  cbind(i = idx[, 1], j = idx[, 2])
}

simulate_osbm_toy <- function(n, z, eta, kappa, regime = c("WST", "SST"), psi) {
  regime <- match.arg(regime)
  A <- matrix(0L, n, n)
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    lam <- eta[i] * eta[j] * kappa[z[i], z[j]]
    Nij <- rpois(1L, lam)
    if (Nij == 0L) next
    if (z[i] == z[j]) {
      p <- 0.5
    } else if (regime == "WST") {
      a <- min(z[i], z[j]); b <- max(z[i], z[j])
      p <- plogis(sign(z[j] - z[i]) * psi[a, b])
    } else {
      d <- abs(z[i] - z[j])
      p <- plogis(sign(z[j] - z[i]) * psi[d])
    }
    Aij <- rbinom(1L, Nij, p)
    A[i, j] <- Aij
    A[j, i] <- Nij - Aij
  }
  A
}

simulate_dcsbm_toy <- function(n, z, theta_out, theta_in, lambda) {
  A <- matrix(0L, n, n)
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    mu_ij <- theta_out[i] * theta_in[j] * lambda[z[i], z[j]]
    mu_ji <- theta_out[j] * theta_in[i] * lambda[z[j], z[i]]
    A[i, j] <- rpois(1L, mu_ij)
    A[j, i] <- rpois(1L, mu_ji)
  }
  A
}

slow_loglik_matrix_osbm <- function(A, mcmc_out, regime = c("WST", "SST"), dyad_index) {
  regime <- match.arg(regime)
  A <- as.matrix(A)
  n <- nrow(A)
  stopifnot(ncol(A) == n)

  Z <- mcmc_out$z
  ETA <- mcmc_out$eta
  KAP <- mcmc_out$kappa
  PSI <- mcmc_out$psi

  S <- nrow(Z)
  M <- nrow(dyad_index)
  LL <- matrix(NA_real_, nrow = S, ncol = M)

  for (s in seq_len(S)) {
    z <- Z[s, ]
    eta <- ETA[s, ]
    kappa <- KAP[s, , ]
    for (m in seq_len(M)) {
      i <- dyad_index[m, 1]
      j <- dyad_index[m, 2]
      Aij <- A[i, j]
      Aji <- A[j, i]
      Nij <- Aij + Aji

      lam <- eta[i] * eta[j] * kappa[z[i], z[j]]

      if (z[i] == z[j]) {
        p <- 0.5
      } else if (regime == "WST") {
        a <- min(z[i], z[j]); b <- max(z[i], z[j])
        p <- plogis(sign(z[j] - z[i]) * PSI[s, a, b])
      } else {
        d <- abs(z[i] - z[j])
        p <- plogis(sign(z[j] - z[i]) * PSI[s, d])
      }

      LL[s, m] <- dpois(Nij, lam, log = TRUE) + dbinom(Aij, Nij, p, log = TRUE)
    }
  }

  LL
}

slow_loglik_matrix_osbm_ragged <- function(A, mcmc_out, regime = c("WST", "SST"), dyad_index) {
  regime <- match.arg(regime)
  A <- as.matrix(A)
  n <- nrow(A)
  stopifnot(ncol(A) == n)

  Z <- mcmc_out$z
  ETA <- mcmc_out$eta
  KAP <- mcmc_out$kappa
  PSI <- mcmc_out$psi

  S <- length(Z)
  M <- nrow(dyad_index)
  LL <- matrix(NA_real_, nrow = S, ncol = M)

  for (s in seq_len(S)) {
    z <- as.integer(Z[[s]])
    eta <- as.numeric(ETA[[s]])
    kappa <- as.matrix(KAP[[s]])
    psi <- PSI[[s]]

    for (m in seq_len(M)) {
      i <- dyad_index[m, 1]
      j <- dyad_index[m, 2]
      Aij <- A[i, j]
      Aji <- A[j, i]
      Nij <- Aij + Aji

      lam <- eta[i] * eta[j] * kappa[z[i], z[j]]

      if (z[i] == z[j]) {
        p <- 0.5
      } else if (regime == "WST") {
        a <- min(z[i], z[j]); b <- max(z[i], z[j])
        p <- plogis(sign(z[j] - z[i]) * psi[a, b])
      } else {
        d <- abs(z[i] - z[j])
        p <- plogis(sign(z[j] - z[i]) * psi[d])
      }

      LL[s, m] <- dpois(Nij, lam, log = TRUE) + dbinom(Aij, Nij, p, log = TRUE)
    }
  }

  LL
}

slow_loglik_matrix_dcsbm <- function(A, dcsbm_out, dyad_index) {
  A <- as.matrix(A)
  n <- nrow(A)
  stopifnot(ncol(A) == n)

  Z <- dcsbm_out$z
  THO <- dcsbm_out$theta_out
  THI <- dcsbm_out$theta_in
  if (is.null(THO) || is.null(THI)) {
    TH <- dcsbm_out$theta
    THO <- TH
    THI <- TH
  }
  LAM <- dcsbm_out$lambda
  S <- if (is.matrix(Z)) nrow(Z) else length(Z)
  M <- nrow(dyad_index)
  LL <- matrix(NA_real_, nrow = S, ncol = M)

  get_row <- function(x, s) if (is.matrix(x)) x[s, ] else x[[s]]
  get_z <- function(z, s) if (is.matrix(z)) z[s, ] else z[[s]]
  get_lam <- function(lam, s) {
    if (is.list(lam)) return(lam[[s]])
    lam[s, , ]
  }

  for (s in seq_len(S)) {
    z <- get_z(Z, s)
    tho <- get_row(THO, s)
    thi <- get_row(THI, s)
    lam <- get_lam(LAM, s)

    for (m in seq_len(M)) {
      i <- dyad_index[m, 1]
      j <- dyad_index[m, 2]
      mu_ij <- tho[i] * thi[j] * lam[z[i], z[j]]
      mu_ji <- tho[j] * thi[i] * lam[z[j], z[i]]
      LL[s, m] <- dpois(A[i, j], mu_ij, log = TRUE) + dpois(A[j, i], mu_ji, log = TRUE)
    }
  }

  LL
}

run_application_loglik_checks <- function(seed = 20260304L, n = 10L, K = 3L, S = 8L, tol = 1e-10) {
  set.seed(seed)

  z_true <- rep(seq_len(K), length.out = n)
  eta_true <- runif(n, 0.8, 1.2)
  eta_true <- eta_true * n / sum(eta_true)
  kappa_true <- matrix(rgamma(K * K, shape = 2, rate = 1.5), K, K)
  kappa_true[lower.tri(kappa_true)] <- t(kappa_true)[lower.tri(kappa_true)]

  psi_wst_true <- matrix(0, K, K)
  psi_wst_true[upper.tri(psi_wst_true)] <- runif(K * (K - 1) / 2, 0.2, 1.5)
  psi_wst_true[lower.tri(psi_wst_true)] <- t(psi_wst_true)[lower.tri(psi_wst_true)]

  psi_sst_true <- cumsum(runif(K - 1L, 0.2, 0.8))

  A_wst <- simulate_osbm_toy(n, z_true, eta_true, kappa_true, regime = "WST", psi = psi_wst_true)
  A_sst <- simulate_osbm_toy(n, z_true, eta_true, kappa_true, regime = "SST", psi = psi_sst_true)

  z_draws <- matrix(rep(z_true, each = S), nrow = S, byrow = TRUE)
  eta_draws <- replicate(S, pmax(eta_true + rnorm(n, 0, 0.03), 1e-6), simplify = "matrix")
  eta_draws <- t(eta_draws)

  kappa_draws <- array(0, dim = c(S, K, K))
  psi_wst_draws <- array(0, dim = c(S, K, K))
  psi_sst_draws <- matrix(0, nrow = S, ncol = K - 1L)
  for (s in seq_len(S)) {
    kap_s <- pmax(kappa_true + matrix(rnorm(K * K, 0, 0.05), K, K), 1e-6)
    kap_s[lower.tri(kap_s)] <- t(kap_s)[lower.tri(kap_s)]
    kappa_draws[s, , ] <- kap_s

    pw_s <- pmax(psi_wst_true + matrix(rnorm(K * K, 0, 0.03), K, K), 0)
    pw_s[lower.tri(pw_s)] <- t(pw_s)[lower.tri(pw_s)]
    diag(pw_s) <- 0
    psi_wst_draws[s, , ] <- pw_s

    ps_s <- pmax(psi_sst_true + rnorm(K - 1L, 0, 0.03), 1e-6)
    psi_sst_draws[s, ] <- cummax(ps_s)
  }

  di_wst <- build_dyad_index_local(nrow(A_wst))
  out_wst <- list(z = z_draws, eta = eta_draws, kappa = kappa_draws, psi = psi_wst_draws)
  ll_fast_wst <- loglik_matrix_wst(A_wst, out_wst, di_wst)
  ll_slow_wst <- slow_loglik_matrix_osbm(A_wst, out_wst, regime = "WST", dyad_index = di_wst)

  out_wst_ragged <- list(
    z = lapply(seq_len(S), function(s) as.integer(z_draws[s, ])),
    eta = lapply(seq_len(S), function(s) as.numeric(eta_draws[s, ])),
    kappa = lapply(seq_len(S), function(s) kappa_draws[s, , ]),
    psi = lapply(seq_len(S), function(s) psi_wst_draws[s, , ])
  )
  ll_fast_wst_mod <- loglik_matrix_modular(A_wst, out_wst_ragged, regime = "WST", dyad_index = di_wst)
  ll_slow_wst_mod <- slow_loglik_matrix_osbm_ragged(A_wst, out_wst_ragged, regime = "WST", dyad_index = di_wst)

  di_sst <- build_dyad_index_local(nrow(A_sst))
  out_sst <- list(z = z_draws, eta = eta_draws, kappa = kappa_draws, psi = psi_sst_draws)
  ll_fast_sst <- loglik_matrix_sst(A_sst, out_sst, di_sst)
  ll_slow_sst <- slow_loglik_matrix_osbm(A_sst, out_sst, regime = "SST", dyad_index = di_sst)

  # Structural check: SST is nested in WST via psi_{k,l} = psi_{|k-l|}.
  psi_wst_equiv <- matrix(0, K, K)
  for (k in seq_len(K)) for (l in seq_len(K)) {
    if (k != l) psi_wst_equiv[k, l] <- psi_sst_true[abs(k - l)]
  }
  out_sst_1draw <- list(
    z = matrix(z_true, nrow = 1),
    eta = matrix(eta_true, nrow = 1),
    kappa = array(kappa_true, dim = c(1, K, K)),
    psi = matrix(psi_sst_true, nrow = 1)
  )
  out_wst_equiv_1draw <- list(
    z = matrix(z_true, nrow = 1),
    eta = matrix(eta_true, nrow = 1),
    kappa = array(kappa_true, dim = c(1, K, K)),
    psi = array(psi_wst_equiv, dim = c(1, K, K))
  )
  ll_sst_equiv <- loglik_matrix_sst(A_sst, out_sst_1draw, di_sst)
  ll_wst_equiv <- loglik_matrix_wst(A_sst, out_wst_equiv_1draw, di_sst)
  err_nested <- max(abs(ll_sst_equiv - ll_wst_equiv))

  out_sst_ragged <- list(
    z = lapply(seq_len(S), function(s) as.integer(z_draws[s, ])),
    eta = lapply(seq_len(S), function(s) as.numeric(eta_draws[s, ])),
    kappa = lapply(seq_len(S), function(s) kappa_draws[s, , ]),
    psi = lapply(seq_len(S), function(s) as.numeric(psi_sst_draws[s, ]))
  )
  ll_fast_sst_mod <- loglik_matrix_modular(A_sst, out_sst_ragged, regime = "SST", dyad_index = di_sst)
  ll_slow_sst_mod <- slow_loglik_matrix_osbm_ragged(A_sst, out_sst_ragged, regime = "SST", dyad_index = di_sst)

  theta_out_true <- runif(n, 0.7, 1.3)
  theta_in_true <- runif(n, 0.7, 1.3)
  lambda_true <- matrix(rgamma(K * K, shape = 1.8, rate = 1.2), K, K)
  A_dc <- simulate_dcsbm_toy(n, z_true, theta_out_true, theta_in_true, lambda_true)

  lambda_list <- vector("list", S)
  theta_out_draws <- matrix(0, nrow = S, ncol = n)
  theta_in_draws <- matrix(0, nrow = S, ncol = n)
  for (s in seq_len(S)) {
    lambda_list[[s]] <- pmax(lambda_true + matrix(rnorm(K * K, 0, 0.05), K, K), 1e-6)
    theta_out_draws[s, ] <- pmax(theta_out_true + rnorm(n, 0, 0.03), 1e-6)
    theta_in_draws[s, ] <- pmax(theta_in_true + rnorm(n, 0, 0.03), 1e-6)
  }
  out_dc <- list(
    z = z_draws,
    theta_out = theta_out_draws,
    theta_in = theta_in_draws,
    lambda = lambda_list
  )
  di_dc <- build_dyad_index_local(nrow(A_dc))
  ll_fast_dc <- loglik_matrix_dcsbm(A_dc, out_dc, di_dc)
  ll_slow_dc <- slow_loglik_matrix_dcsbm(A_dc, out_dc, di_dc)

  err_wst <- max(abs(ll_fast_wst - ll_slow_wst))
  err_sst <- max(abs(ll_fast_sst - ll_slow_sst))
  err_wst_mod <- max(abs(ll_fast_wst_mod - ll_slow_wst_mod))
  err_sst_mod <- max(abs(ll_fast_sst_mod - ll_slow_sst_mod))
  err_dc <- max(abs(ll_fast_dc - ll_slow_dc))

  res <- data.frame(
    model = c(
      "fast_helper_OSBM-WST",
      "fast_helper_OSBM-SST",
      "fast_application_OSBM-WST-modular",
      "fast_application_OSBM-SST-modular",
      "fast_application_DC-SBM",
      "nested_equivalence_SST_vs_WST"
    ),
    max_abs_diff = c(err_wst, err_sst, err_wst_mod, err_sst_mod, err_dc, err_nested),
    tolerance = tol,
    pass = c(
      err_wst <= tol,
      err_sst <= tol,
      err_wst_mod <= tol,
      err_sst_mod <= tol,
      err_dc <= tol,
      err_nested <= tol
    )
  )

  print(res, row.names = FALSE)
  if (!all(res$pass)) {
    stop("Fast-vs-slow log-likelihood mismatch exceeds tolerance.", call. = FALSE)
  }

  invisible(
    list(
      summary = res,
      toy = list(A_wst = A_wst, A_sst = A_sst, A_dc = A_dc),
      ll = list(
        fast_wst = ll_fast_wst, slow_wst = ll_slow_wst,
        fast_sst = ll_fast_sst, slow_sst = ll_slow_sst,
        fast_wst_mod = ll_fast_wst_mod, slow_wst_mod = ll_slow_wst_mod,
        fast_sst_mod = ll_fast_sst_mod, slow_sst_mod = ll_slow_sst_mod,
        fast_dc = ll_fast_dc, slow_dc = ll_slow_dc
      )
    )
  )
}

run_loglik_toy_checks <- function(seed = 20260304L, n = 10L, K = 3L, S = 8L, tol = 1e-10) {
  run_application_loglik_checks(seed = seed, n = n, K = K, S = S, tol = tol)
}

# Example run (adds this check to the existing battery style of tester scripts).
if (sys.nframe() == 0L) {
  run_application_loglik_checks()
  cat("Log-likelihood toy checks (application path + helper path) passed.\n")
}
