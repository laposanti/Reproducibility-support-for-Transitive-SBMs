#!/usr/bin/env Rscript
# Quick smoke test for the ordered-SBM samplers on tiny synthetic datasets.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", script_arg[1L])))
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
  setwd(repo_root)
}

source("helper_folder/simulation/simulation_study_helpers.R")
source("helper_folder/models/ordered_sbm/sst_helpers.R")
source("helper_folder/models/ordered_sbm/wst_helpers.R")
source("helper_folder/config/hyperparameter_setup.R")
source("core/transitive_sbm_sampler.R")

set.seed(42)
n <- 30; K_true <- 3
z_true <- rep(1:K_true, length.out = n)
kappa <- matrix(runif(K_true^2, 0.5, 2), K_true, K_true)
kappa <- (kappa + t(kappa)) / 2
eta <- runif(n, 0.8, 1.2)

# WST test
psi_wst <- matrix(0, K_true, K_true)
psi_wst[1,2] <- 1.5; psi_wst[1,3] <- 2.0; psi_wst[2,3] <- 1.0
A <- as.matrix(simulate_osbm(n, K_true, z_true, eta, kappa, psi_wst, "WST"))

hypers <- get_corollary_calibrated_hypers(A, K_expected = 1,
  ordering_prior_mode = "equivalence_class", a_kappa = 2, a_eta = 2,
  mu0 = 1.0, gamma_bounds = c(0.3, 0.7))

cat("WST quick run (200 iter, use_mixing_moves=TRUE)...\n")
t0 <- proc.time()[3]
out <- modular_osbm_sampler(
  A = A, K = 5, n_iter = 200, burn = 50, thin = 1,
  psi_constraint = "WST" , partition_prior = "GN",
  gamma_gn = hypers$gamma, a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
  a_eta = hypers$a_eta, b_eta = hypers$b_eta,
  mu0 = hypers$mu0, sigma0 = hypers$sigma0, tau0 = hypers$tau0,
  use_mixing_moves = TRUE, seed = 42, verbose = FALSE)
t1 <- proc.time()[3]

z_chain <- do.call(rbind, lapply(out$z, as.integer))
K_trace <- apply(z_chain, 1, function(x) length(unique(x)))
z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L)
ari <- fossil::adj.rand.index(z_hat, z_true)
cat(sprintf("  Done in %.1fs. ARI=%.3f, K_mode=%d, K_range=[%d,%d]\n",
            t1 - t0, ari, as.integer(names(which.max(table(K_trace)))),
            min(K_trace), max(K_trace)))

# SST test
psi_sst <- cumsum(c(0.8, 0.6))
A_sst <- as.matrix(simulate_osbm(n, K_true, z_true, eta, kappa, psi_sst, "SST"))
cat("SST quick run (200 iter, use_mixing_moves=TRUE)...\n")
t0 <- proc.time()[3]
out2 <- modular_osbm_sampler(
  A = A_sst, K = 5, n_iter = 200, burn = 50, thin = 1,
  psi_constraint = "SST" , partition_prior = "GN",
  gamma_gn = hypers$gamma, a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
  a_eta = hypers$a_eta, b_eta = hypers$b_eta,
  mu0 = hypers$mu0, sigma0 = hypers$sigma0, tau0 = hypers$tau0,
  use_mixing_moves = TRUE, seed = 42, verbose = FALSE)
t1 <- proc.time()[3]

z_chain2 <- do.call(rbind, lapply(out2$z, as.integer))
K_trace2 <- apply(z_chain2, 1, function(x) length(unique(x)))
z_hat2 <- salso::salso(z_chain2, loss = salso::VI(), nRuns = 1L)
ari2 <- fossil::adj.rand.index(z_hat2, z_true)
cat(sprintf("  Done in %.1fs. ARI=%.3f, K_mode=%d, K_range=[%d,%d]\n",
            t1 - t0, ari2, as.integer(names(which.max(table(K_trace2)))),
            min(K_trace2), max(K_trace2)))
cat("All quick checks passed.\n")
