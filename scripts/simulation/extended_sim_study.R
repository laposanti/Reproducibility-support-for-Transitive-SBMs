# --- Required libraries -----------------------------------------------------
library(Matrix)
library(fossil)
library(mcclust)
library(coda)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggside)
library(mcclust.ext)
library(RColorBrewer)
library(scales)
library(loo)   # for WAIC / LOO
library(digest) # for deterministic per-cell seeds

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(filelock)
  library(data.table)
  library(progressr)
})

# --- Project setup & sources ------------------------------------------------
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")
source("./core/modular mcmc.R")
source("./helper_folder/helper.R")
source("./helper_folder/plugin_stream_save.R")

## =========================================================
## Utilities for simulation design
## =========================================================
logistic <- function(x) 1/(1+exp(-x))

## Generate psi by model, then enforce model’s support if needed
gen_psi <- function(K, model = c("SST","WST","NONE")) {
  model <- match.arg(model)
  stopifnot(K > 1)
  D <- K - 1
  base <- seq(1, 5, length.out = D)
  if (model == "SST") {
    psi_t <- base
  } else if (model == "WST") {
    psi_t <- rev(base)           # weak: shuffled positive increasing values
  } else { # "NONE"
    signs <- rep(c(1, -1), length.out = D)
    psi_t <- base * signs           # unconstrained signs
  }
  return(psi_t)
}

# Scale factor to hit a target mean rate
calibrate_tau <- function(z, eta, kappa, target_mean = 2) {
  n <- length(z)
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  lam_base <- eta[idx[,1]] * eta[idx[,2]] * kappa[cbind(z[idx[,1]], z[idx[,2]])]
  target_mean / mean(lam_base)
}

gen_kappa <- function(K, difficulty = c("dense","medium","sparse")) {
  difficulty <- match.arg(difficulty)
  disp <- switch(difficulty, dense = 1.0, medium = 0.6, sparse = 0.3)
  M <- matrix(rexp(K*K, rate = 1), K, K)
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  diag(M) <- pmax(diag(M), 0.5)
  M^disp
}

gen_eta <- function(n, z, K) {
  eta <- runif(n, 0.8, 1.2)
  n_k <- as.integer(tabulate(z, nbins = K))
  for (k in seq_len(K)) {
    idx <- which(z == k)
    if (length(idx) == 0L) next
    s_k <- sum(eta[idx])
    if (s_k > 0) {
      eta[idx] <- n_k[k] * eta[idx] / s_k
    } else {
      eta[idx] <- n_k[k] / length(idx)
    }
  }
  eta
}

## Data generator 
generate_osbm_data <- function(n, K,
                               gen_model = c("SST","WST","NONE"),
                               difficulty = c("dense","medium","sparse")) {
  gen_model <- match.arg(gen_model)
  difficulty <- match.arg(difficulty)
  
  z_t   <- rep(1:K, length.out = n)
  eta_t <- gen_eta(n, z_t, K)
  psi_t <- gen_psi(K, model = gen_model)
  kappa_t <- gen_kappa(K, difficulty)
  
  tau <- calibrate_tau(z_t, eta_t, kappa_t, 
                       target_mean = switch(difficulty, dense=3, medium=2, sparse=1))
  kappa_t <- tau * kappa_t
  
  A <- Matrix::Matrix(0, n, n, sparse = TRUE)
  for (i in 1:(n-1)) for (j in (i+1):n) {
    lam <- eta_t[i]*eta_t[j]*kappa_t[z_t[i], z_t[j]]
    N   <- rpois(1, lam)
    if (N>0) {
      d  <- abs(z_t[i]-z_t[j])
      pr <- if (d==0) 0.5 else logistic(if (z_t[i]<z_t[j]) psi_t[d] else -psi_t[d])
      fwd <- rbinom(1, N, pr)
      A[i,j] <- fwd; A[j,i] <- N - fwd
    }
  }
  truth <- list(z = z_t, eta = eta_t, kappa = kappa_t, psi = psi_t)
  list(A=A, truth=truth, meta=list(gen_model=gen_model, difficulty=difficulty))
}

# --- Example run -------------------------------------------------------------
res <- run_simulation_study_streaming(
  K_vals = c(3,5,7),
  n_vals = 80,
  gen_model = c('SST','WST',"NONE"),
  difficulties = c("dense","sparse"),
  n_reps = 1,
  n_iter = 5000, burn = 1000, thin = 2,
  seed = 2021,
  out_path = "./sim_results_stream_another_try5000.csv",
  workers = 4,
  plan_strategy = "multisession",
  verbose = T,
  use_progress = T
)

# Inspect results
df <- read_stream_results("sim_results_stream_another_try5000.csv")
conf <- confusion_from_stream(df)
print(head(df))
print(conf)

df %>%
  dplyr::filter(chosen, difficulty == 'sparse') |>
  dplyr::group_by(gen_model, difficulty) |>
  dplyr::summarise(mean_ari = mean(ari), .groups = "drop")

df %>%
  dplyr::filter(chosen, difficulty == 'sparse') |>
  dplyr::group_by(gen_model, fit_model) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
  pivot_wider(names_from = fit_model, values_from = n)

df %>%
  ggplot(aes(x = fit_model,y = ari))+
  geom_boxplot()




df%>%
  group_by(difficulty)%>%
  filter(fit_model == 'SST',gen_model == 'SST')%>%
  select(K,mae_kappa)



