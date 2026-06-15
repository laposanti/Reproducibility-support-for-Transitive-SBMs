#!/usr/bin/env Rscript

## Robustness checks for the OSBM MCMC sampler.
## Exercises small a_eta values that have caused occasional failures.

suppressPackageStartupMessages({
  library(Matrix)
})

source("helper_folder/helper.R")
source("core/modular mcmc.R")

simulate_osbm_data <- function(n, K, psi_constraint = c("WST", "SST"), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  psi_constraint <- match.arg(psi_constraint)
  z <- sample.int(K, n, replace = TRUE)
  eta <- rgamma(n, shape = 2, rate = 2)
  eta <- eta * n / sum(eta)
  kappa <- matrix(1, K, K)
  kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]
  if (psi_constraint == "WST") {
    psi <- matrix(0, K, K)
    psi[upper.tri(psi)] <- seq_len(K * (K - 1) / 2) / (K * (K - 1) / 2)
  } else {
    psi <- seq(1, K, length.out = K - 1)
  }
  
  A <- Matrix(0, n, n, sparse = TRUE)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      lam <- eta[i] * eta[j] * kappa[z[i], z[j]]
      Nij <- rpois(1, lam)
      if (Nij > 0) {
        if (psi_constraint == "WST") {
          pr <- if (z[i] < z[j]) {
            logistic(psi[z[i], z[j]])
          } else if (z[i] > z[j]) {
            logistic(-psi[z[j], z[i]])
          } else {
            0.5
          }
        } else {
          d <- abs(z[i] - z[j])
          pr <- if (d == 0) {
            0.5
          } else if (z[i] < z[j]) {
            logistic(psi[d])
          } else {
            logistic(-psi[d])
          }
        }
        fwd <- rbinom(1, Nij, pr)
        A[i, j] <- fwd
        A[j, i] <- Nij - fwd
      }
    }
  }
  list(A = A, K = K)
}

check_mcmc_output <- function(out, K) {
  stopifnot(is.matrix(out$z))
  if (any(!is.finite(out$eta))) stop("Non-finite eta encountered.")
  if (any(!is.finite(out$kappa))) stop("Non-finite kappa encountered.")
  if (any(!is.finite(out$psi))) stop("Non-finite psi encountered.")
  if (any(out$z < 1 | out$z > K)) stop("Invalid z labels encountered.")
}

run_robustness_suite <- function() {
  a_eta_values <- c(0.1, 0.25, 0.4, 0.5, 1)
  settings <- expand.grid(
    psi_constraint = c("WST", "SST"),
    a_eta = a_eta_values,
    seed = c(11, 37, 91),
    stringsAsFactors = FALSE
  )
  
  for (row in seq_len(nrow(settings))) {
    cfg <- settings[row, ]
    sim <- simulate_osbm_data(n = 30, K = 3, psi_constraint = cfg$psi_constraint,
                              seed = cfg$seed)
    out <- modular_osbm_sampler(
      A = sim$A,
      K = sim$K,
      n_iter = 80,
      burn = 20,
      thin = 2,
      psi_constraint = cfg$psi_constraint,
      hyper = list(a_eta = cfg$a_eta, b_eta = 1),
      seed = cfg$seed
    )
    check_mcmc_output(out, sim$K)
  }
  message("Robustness checks completed without non-finite values.")
}

run_robustness_suite()
