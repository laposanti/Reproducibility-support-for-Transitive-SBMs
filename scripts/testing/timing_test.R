#!/usr/bin/env Rscript
# Quick timing test
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")
source("./helper_folder/helper.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./helper_folder/Hyper_setup.R")
source("./core/my_best_try_so_far.R")

matrix_files <- list.files("./data/ShizukaMcDonald_Data",
                           full.names = TRUE, pattern = "[.]csv$")
n_each <- vapply(matrix_files, function(f) nrow(read.csv(f, row.names = 1)),
                 FUN.VALUE = integer(1))
A <- as.matrix(read.csv(matrix_files[which.max(n_each)],
                        row.names = 1, check.names = FALSE))
cat("n =", nrow(A), "\n")

t0 <- proc.time()
out <- modular_osbm_sampler(
  A = A, K = 5,
  n_iter = 100, burn = 1, thin = 1,
  verbose = 50,
  psi_constraint = "SST" ,
  partition_prior = "GN",
  gamma_gn = 0.96,
  a_kappa = 1, b_kappa = 22.5,
  a_eta = 7.7, b_eta = 7.7,
  mu0 = 0, sigma0 = 0.2, tau0 = 0.2,
  use_mixing_moves = TRUE,
  seed = 42
)
elapsed <- (proc.time() - t0)[["elapsed"]]
Kt <- out$K_trace
Kt <- Kt[!is.na(Kt)]
cat("\n100 iterations in", round(elapsed, 1), "seconds\n")
cat("K trace:", Kt, "\n")
cat("Per-iteration:", round(elapsed / 100, 3), "seconds\n")
