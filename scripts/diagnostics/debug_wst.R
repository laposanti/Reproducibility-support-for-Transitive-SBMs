## Quick debug script for WST crash
suppressPackageStartupMessages({
  library(BayesLogit); library(truncnorm); library(Matrix); library(Rcpp)
})
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")
source("./helper_folder/helper.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./core/my_best_try_so_far.R")

set.seed(42)
K_true <- 3; n <- 30
z_true <- rep(1:K_true, each = n / K_true)
kappa_true <- matrix(0.5, K_true, K_true); diag(kappa_true) <- 3
eta_true <- rep(1, n)
psi_true_sst <- c(1.0, 1.5)
A <- matrix(0L, n, n)
for (i in 1:(n-1)) for (j in (i+1):n) {
  ki <- z_true[i]; kj <- z_true[j]
  N_ij <- rpois(1, eta_true[i]*eta_true[j]*kappa_true[ki,kj])
  if (N_ij > 0) {
    d <- abs(ki-kj)
    rho <- if(d==0) 0.5 else 1/(1+exp(-sign(kj-ki)*psi_true_sst[d]))
    a_ij <- rbinom(1, N_ij, rho); A[i,j] <- a_ij; A[j,i] <- N_ij-a_ij
  }
}

options(error = function() { traceback(2); q(status=1) })

cat("Running WST sampler...\n")
res <- modular_osbm_sampler(
  A=A, K=5, free=c("psi","kappa","eta","z"),
  n_iter=50, burn=10, thin=1, verbose=TRUE,
  psi_constraint="WST", seed=456, gamma_gn=0.8,
  mu0=1, sigma0=2, tau0=0.15
)
cat("Done. K range:", range(res$K_trace), "\n")
