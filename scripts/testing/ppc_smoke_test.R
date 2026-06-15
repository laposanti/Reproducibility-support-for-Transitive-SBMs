#!/usr/bin/env Rscript
# Minimal PPC smoke-test: fit 3 models on mountain_goats with tiny chains,
# then run PPC and save to posterior_predictive_checks/

cat("=== PPC smoke test ===\n")
gamma_gn <- 0.9
cat("Using GN prior gamma =", gamma_gn, "for all models\n")

suppressPackageStartupMessages({
  library(Matrix); library(mcclust); library(mcclust.ext); library(coda)
  library(dplyr); library(tidyr); library(loo); library(here); library(fs)
  library(readr); library(salso); library(glue); library(label.switching)
  library(purrr); library(bayesplot); library(ggplot2)
})

source("core/ppc_checks.R")
source("helper_folder/sim_study_helper.R")
source("helper_folder/transitivity_check_helper.R")
source("core/DCSBM_varK.R")
source("helper_folder/helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("core/my_best_try_so_far.R")

# ---- Load mountain_goats data ----
matrix_files <- list.files("data/ShizukaMcDonald_Data", full.names = TRUE,
                           pattern = "[.]csv$")
n_each <- vapply(matrix_files,
                 function(f) nrow(read.csv(f, row.names = 1)),
                 integer(1))
A <- as.matrix(read.csv(matrix_files[which.max(n_each)],
                        row.names = 1, check.names = FALSE))
cat("Mountain goats A:", nrow(A), "x", ncol(A), "\n")

# ---- Set up hyperparameters ----
hypers <- suppressWarnings(get_corollary_calibrated_hypers(
  A, K_expected = 1, ordering_prior_mode = "equivalence_class",
  a_kappa = 2, a_eta = 2, mu0 = 1.0, gamma_bounds = c(0.3, 0.7)
))
cat("Hypers OK (calibrated gamma ignored; using fixed gamma =", gamma_gn, ")\n")

# ---- Output directory ----
ppc_dir_base <- here::here("output", "application", "ppc", "mountain_goats")
fs::dir_create(ppc_dir_base)

n_iter <- 300L
burn <- 50L
thin <- 2L
seed <- 102L

# ================ SST ================
cat("\n--- Fitting SST ---\n")
out_sst <- tryCatch(
  modular_osbm_sampler(
    A = A, K = nrow(A), n_iter = n_iter, burn = burn, thin = thin,
    verbose = TRUE,
    a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
    a_eta = hypers$a_eta, b_eta = hypers$b_eta,
    psi_constraint = "SST" , partition_prior = "GN",
    gamma_gn = gamma_gn, mu0 = hypers$mu0,
    sigma0 = hypers$sigma0, tau0 = hypers$tau0, seed = seed
  ),
  error = function(e) { cat("SST ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(out_sst)) {
  cat("SST fit OK, draws:", length(out_sst$z), "\n")
  ppc_dir_sst <- file.path(ppc_dir_base, "SST")
  fs::dir_create(ppc_dir_sst)
  ppc_sst <- tryCatch(
    .run_ppc_model(A, out_sst, "SST",
                   n_draws = 50L, m_triples = 500L,
                   seed = 999L, out_dir = ppc_dir_sst,
                   tag = "_mountain_goats"),
    error = function(e) { cat("SST PPC ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(ppc_sst)) {
    cat("SST PPC completed:\n")
    print(ppc_sst$pvals)
  }
} else {
  cat("SST fit FAILED\n")
}

# ================ WST ================
cat("\n--- Fitting WST ---\n")
out_wst <- tryCatch(
  modular_osbm_sampler(
    A = A, K = nrow(A), n_iter = n_iter, burn = burn, thin = thin,
    verbose = TRUE,
    a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
    a_eta = hypers$a_eta, b_eta = hypers$b_eta,
    psi_constraint = "WST" , partition_prior = "GN",
    gamma_gn = gamma_gn, mu0 = hypers$mu0,
    sigma0 = hypers$sigma0, tau0 = hypers$tau0, seed = seed
  ),
  error = function(e) { cat("WST ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(out_wst)) {
  cat("WST fit OK, draws:", length(out_wst$z), "\n")
  ppc_dir_wst <- file.path(ppc_dir_base, "WST")
  fs::dir_create(ppc_dir_wst)
  ppc_wst <- tryCatch(
    .run_ppc_model(A, out_wst, "WST",
                   n_draws = 50L, m_triples = 500L,
                   seed = 999L, out_dir = ppc_dir_wst,
                   tag = "_mountain_goats"),
    error = function(e) { cat("WST PPC ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(ppc_wst)) {
    cat("WST PPC completed:\n")
    print(ppc_wst$pvals)
  }
} else {
  cat("WST fit FAILED\n")
}

# ================ DCSBM ================
cat("\n--- Fitting DCSBM ---\n")
out_dc <- tryCatch(
  fit_dcsbm_gibbs_gnedin(
    as.matrix(A), K_init = nrow(A),
    priors = list(
      a_eta = hypers$a_eta, b_eta = hypers$b_eta,
      a_lambda = hypers$a_kappa, b_lambda = hypers$b_kappa,
      gamma_gnedin = gamma_gn
    ),
    iters = n_iter, burn_in = burn, thin = thin,
    verbose = 200, seed = seed
  ),
  error = function(e) { cat("DCSBM ERROR:", conditionMessage(e), "\n"); NULL }
)

if (!is.null(out_dc)) {
  cat("DCSBM fit OK, draws:", nrow(out_dc$z), "\n")

  # Relabel
  out_dc_relab <- tryCatch(
    relabel_by_ECR(mcmc_out = out_dc, model = "DCSBM"),
    error = function(e) { cat("DCSBM relabel ERROR:", conditionMessage(e), "\n"); out_dc }
  )

  ppc_dir_dc <- file.path(ppc_dir_base, "DCSBM")
  fs::dir_create(ppc_dir_dc)
  ppc_dc <- tryCatch(
    .run_ppc_model(A, out_dc_relab, "DCSBM",
                   n_draws = 50L, m_triples = 500L,
                   seed = 999L, out_dir = ppc_dir_dc,
                   tag = "_mountain_goats"),
    error = function(e) { cat("DCSBM PPC ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(ppc_dc)) {
    cat("DCSBM PPC completed:\n")
    print(ppc_dc$pvals)
  }
} else {
  cat("DCSBM fit FAILED\n")
}

cat("\n=== Done. Check posterior_predictive_checks/mountain_goats/ ===\n")
