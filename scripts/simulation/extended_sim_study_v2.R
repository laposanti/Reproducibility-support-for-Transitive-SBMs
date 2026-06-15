#!/usr/bin/env Rscript
# =============================================================================
# Extended Simulation Study v2 — K recovery, mixing, model choice
#
# AXES OF COMPARISON:
#   1. Hyperparameter calibration: default | oracle Rule-A | spectral Rule-A+
#      | calibrate_osbm_hypers | hierarchical b_kappa
#   2. Mixing moves: with vs without
#   3. Model choice: WST vs SST vs DCSBM (cross-fit for LOOIC)
#
# SIMULATED DATA SCENARIOS:
#   - gen_model ∈ {WST, SST}
#   - K_true   ∈ {3, 5, 7}
#   - signal   ∈ {"strong", "weak"}
#   - n        = 60 (default, env-configurable)
#
# REAL DATA:
#   - mountain_goats, citations_data, macaques_data, high_school
#   - Fit: WST, SST, DCSBM (each with hyper variants)
#
# METRICS:
#   K_hat, K_post_mean, K_post_mode, ARI, VI, LOOIC, ESS_K, b_kappa_trace
#
# USAGE:
#   Rscript extended_sim_study_v2.R                     # full run
#   SMOKE=1 Rscript extended_sim_study_v2.R             # quick local check
#   sbatch extended_sim_study_launcher.sh                # SLURM cluster
# =============================================================================

## ---- Environment knobs (overridable via env vars) -----------------------
SMOKE       <- as.integer(Sys.getenv("SMOKE", unset = "0")) == 1L
N_REP       <- as.integer(Sys.getenv("N_REP",   unset = if (SMOKE) "1" else "5"))
N_ITER      <- as.integer(Sys.getenv("N_ITER",  unset = if (SMOKE) "600" else "8000"))
BURN        <- as.integer(Sys.getenv("BURN",    unset = if (SMOKE) "200" else "2000"))
THIN        <- as.integer(Sys.getenv("THIN",    unset = if (SMOKE) "1"  else "2"))
N_SIM       <- as.integer(Sys.getenv("N_SIM",   unset = if (SMOKE) "60" else "60"))
BASE_SEED   <- as.integer(Sys.getenv("BASE_SEED", unset = "2026"))

## ---- Libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(BayesLogit)
  library(truncnorm)
  library(salso)
  library(fossil)
  library(mcclust)
  library(coda)
  library(loo)
  library(dplyr)
  library(tidyr)
})

## ---- Source project code -------------------------------------------------
local_dir <- "/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/"
if (dir.exists(local_dir)) setwd(local_dir)
source("helper_folder/sim_study_helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("helper_folder/transitivity_check_helper.R")
source("helper_folder/helper.R")
source("core/my_best_try_so_far.R")
source("core/DCSBM_varK.R")

## ---- Output setup --------------------------------------------------------
out_dir <- file.path("output", "simulation", "raw", paste0("ext_sim_v2_",
                     format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
log_con  <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)

cat(sprintf("=== Extended simulation study v2 started: %s ===\n", Sys.time()))
cat(sprintf("SMOKE=%s  N_REP=%d  N_ITER=%d  BURN=%d  THIN=%d  N_SIM=%d  SEED=%d\n",
            SMOKE, N_REP, N_ITER, BURN, THIN, N_SIM, BASE_SEED))

sgn      <- function(x) ifelse(x > 0, 1L, ifelse(x < 0, -1L, 0L))
logistic <- function(x) 1 / (1 + exp(-x))

# =============================================================================
# Helper: set up hyperparameters for a given method
# =============================================================================
setup_hypers <- function(A, method, gen_model, K_true) {
  n <- nrow(A)
  # Defaults shared across most methods
  a_kappa <- 1
  b_kappa <- 1
  a_eta   <- 1
  b_eta   <- 1
  mu0     <- 0
  sigma0  <- 2
  tau0    <- 0.15
  gamma_gn <- 0.5
  sample_b   <- FALSE
  alpha0_bk  <- 1.0
  beta0_bk   <- 0.01
  use_moves  <- TRUE

  if (method == "default") {
    # Naive default: b_kappa=1, gamma from K_true
    gamma_gn <- choose_gamma_from_K_expected(n, K_expected = K_true,
                  ordering_prior_mode = "equivalence_class")

  } else if (method == "oracle_ruleA") {
    # Rule A with true K
    b_kappa  <- 5 * n / K_true
    gamma_gn <- choose_gamma_from_K_expected(n, K_expected = K_true,
                  ordering_prior_mode = "equivalence_class")

  } else if (method == "spectral_ruleA") {
    # Rule A+ with spectral K
    K_spec   <- estimate_K_spectral(A)
    b_kappa  <- 5 * n / max(K_spec, 2)
    gamma_gn <- choose_gamma_from_K_expected(n, K_expected = max(K_spec, 2),
                  ordering_prior_mode = "equivalence_class")

  } else if (method == "calibrate_full") {
    # calibrate_osbm_hypers (data-driven, no true K needed)
    cal <- calibrate_osbm_hypers(A, regime = gen_model)
    a_kappa  <- cal$a_kappa;  b_kappa <- cal$b_kappa
    a_eta    <- cal$a_eta;    b_eta   <- cal$b_eta
    mu0      <- cal$mu0;      sigma0  <- cal$sigma0
    tau0     <- cal$tau0
    gamma_gn <- cal$gamma

  } else if (method == "hierarchical") {
    # Hierarchical b_kappa (learns b during MCMC)
    sample_b <- TRUE
    b_kappa  <- 5 * n / K_true   # warm-start at oracle
    gamma_gn <- choose_gamma_from_K_expected(n, K_expected = K_true,
                  ordering_prior_mode = "equivalence_class")

  } else if (method == "no_moves") {
    # Same as oracle but without mixing moves
    b_kappa  <- 5 * n / K_true
    gamma_gn <- choose_gamma_from_K_expected(n, K_expected = K_true,
                  ordering_prior_mode = "equivalence_class")
    use_moves <- FALSE

  } else {
    stop("Unknown method: ", method)
  }

  list(
    a_kappa = a_kappa, b_kappa = b_kappa,
    a_eta = a_eta, b_eta = b_eta,
    mu0 = mu0, sigma0 = sigma0, tau0 = tau0,
    gamma_gn = gamma_gn,
    sample_b_kappa = sample_b,
    alpha0_bkappa  = alpha0_bk,
    beta0_bkappa   = beta0_bk,
    use_mixing_moves = use_moves
  )
}

# =============================================================================
# Helper: run OSBM sampler with given hypers, return summary row
# =============================================================================
run_osbm_fit <- function(A, K_true, z_true, fit_model, hypers, seed,
                         n_iter = N_ITER, burn = BURN, thin_ = THIN) {

  n <- nrow(A)
  K_init <- max(K_true + 3L, 6L)

  t0 <- proc.time()[3]
  out <- tryCatch(
    modular_osbm_sampler(
      A = A, K = K_init,
      n_iter = n_iter, burn = burn, thin = thin_,
      psi_constraint = fit_model ,
      partition_prior = "GN",
      gamma_gn       = hypers$gamma_gn,
      a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
      a_eta   = hypers$a_eta,   b_eta   = hypers$b_eta,
      mu0     = hypers$mu0,     sigma0  = hypers$sigma0,
      tau0    = hypers$tau0,
      use_mixing_moves = hypers$use_mixing_moves,
      sample_b_kappa   = hypers$sample_b_kappa,
      alpha0_bkappa    = hypers$alpha0_bkappa,
      beta0_bkappa     = hypers$beta0_bkappa,
      seed    = seed,
      verbose = FALSE
    ),
    error = function(e) { message("  OSBM ERROR: ", e$message); NULL }
  )
  elapsed <- proc.time()[3] - t0

  if (is.null(out)) {
    return(data.frame(
      K_hat = NA, K_post_mean = NA, K_post_mode = NA,
      ari = NA, vi = NA, elpd_loo = NA, looic = NA,
      ess_K = NA, bk_final = NA, time_sec = elapsed,
      stringsAsFactors = FALSE
    ))
  }

  # Extract K trace and z chain
  z_chain <- do.call(rbind, lapply(out$z, as.integer))
  K_trace <- out$K_trace
  S <- length(K_trace)

  # Partition estimate
  K_occ <- apply(z_chain, 1, function(z) length(unique(z)))
  salso_cap <- as.integer(max(3L, min(ncol(z_chain),
                  quantile(K_occ, 0.99, type = 1) + 1L)))
  z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L,
                         maxNClusters = salso_cap)
  K_hat <- length(unique(z_hat))
  K_post_mean <- mean(K_trace)
  K_post_mode <- as.integer(names(which.max(table(K_trace))))

  # ARI, VI (only when z_true available)
  ari <- if (!is.null(z_true)) fossil::adj.rand.index(z_hat, z_true) else NA_real_
  vi  <- if (!is.null(z_true)) mcclust::vi.dist(z_hat, z_true) else NA_real_

  # LOO-IC
  DI <- build_dyad_index(n)
  LL <- tryCatch(
    loglik_matrix_modular(A, out, regime = fit_model, dyad_index = DI),
    error = function(e) NULL
  )
  elpd_loo <- looic_val <- NA_real_
  if (!is.null(LL)) {
    loo_fit <- tryCatch(loo::loo(LL), error = function(e) NULL)
    if (!is.null(loo_fit)) {
      elpd_loo  <- loo_fit$estimates["elpd_loo", "Estimate"]
      looic_val <- loo_fit$estimates["looic", "Estimate"]
    }
  }

  # Mixing: ESS on K_trace
  ess_K <- tryCatch(
    as.numeric(coda::effectiveSize(coda::mcmc(K_trace))),
    error = function(e) NA_real_
  )

  # b_kappa final value
  bk_final <- if (hypers$sample_b_kappa) mean(tail(out$b_kappa_trace, S %/% 2))
              else hypers$b_kappa

  data.frame(
    K_hat = K_hat, K_post_mean = round(K_post_mean, 2),
    K_post_mode = K_post_mode,
    ari = round(ari, 4), vi = round(vi, 4),
    elpd_loo = round(elpd_loo, 2), looic = round(looic_val, 2),
    ess_K = round(ess_K, 1), bk_final = round(bk_final, 2),
    time_sec = round(elapsed, 1),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Helper: run DCSBM fit and return summary row
# =============================================================================
run_dcsbm_fit <- function(A, K_true, z_true, gamma_gn, seed,
                          n_iter = N_ITER, burn = BURN, thin_ = THIN) {
  n <- nrow(A)
  t0 <- proc.time()[3]
  out <- tryCatch(
    fit_dcsbm_gibbs_gnedin(
      A = as.matrix(A), K_init = 15L,
      priors = list(a_eta = 1, b_eta = 1, a_lambda = 1, b_lambda = 1,
                    gamma_gnedin = gamma_gn),
      iters = n_iter, burn_in = burn, thin = thin_,
      verbose = 0, seed = seed
    ),
    error = function(e) { message("  DCSBM ERROR: ", e$message); NULL }
  )
  elapsed <- proc.time()[3] - t0

  if (is.null(out)) {
    return(data.frame(
      K_hat = NA, K_post_mean = NA, K_post_mode = NA,
      ari = NA, vi = NA, elpd_loo = NA, looic = NA,
      ess_K = NA, bk_final = NA, time_sec = elapsed,
      stringsAsFactors = FALSE
    ))
  }

  # z chain
  if (is.matrix(out$z)) {
    z_chain <- out$z
  } else {
    z_chain <- do.call(rbind, lapply(out$z, as.integer))
  }
  K_occ <- apply(z_chain, 1, function(z) length(unique(z)))
  K_trace <- K_occ
  S <- nrow(z_chain)

  salso_cap <- as.integer(max(3L, min(ncol(z_chain),
                  quantile(K_occ, 0.99, type = 1) + 1L)))
  z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L,
                         maxNClusters = salso_cap)
  K_hat <- length(unique(z_hat))
  K_post_mean <- mean(K_trace)
  K_post_mode <- as.integer(names(which.max(table(K_trace))))

  ari <- if (!is.null(z_true)) fossil::adj.rand.index(z_hat, z_true) else NA_real_
  vi  <- if (!is.null(z_true)) mcclust::vi.dist(z_hat, z_true) else NA_real_

  DI <- build_dyad_index(n)
  LL <- tryCatch(
    loglik_matrix_dcsbm(A = as.matrix(A), dcsbm_out = out, dyad_index = DI),
    error = function(e) NULL
  )
  elpd_loo <- looic_val <- NA_real_
  if (!is.null(LL)) {
    loo_fit <- tryCatch(loo::loo(LL), error = function(e) NULL)
    if (!is.null(loo_fit)) {
      elpd_loo  <- loo_fit$estimates["elpd_loo", "Estimate"]
      looic_val <- loo_fit$estimates["looic", "Estimate"]
    }
  }

  ess_K <- tryCatch(
    as.numeric(coda::effectiveSize(coda::mcmc(K_trace))),
    error = function(e) NA_real_
  )

  data.frame(
    K_hat = K_hat, K_post_mean = round(K_post_mean, 2),
    K_post_mode = K_post_mode,
    ari = round(ari, 4), vi = round(vi, 4),
    elpd_loo = round(elpd_loo, 2), looic = round(looic_val, 2),
    ess_K = round(ess_K, 1), bk_final = NA_real_,
    time_sec = round(elapsed, 1),
    stringsAsFactors = FALSE
  )
}


# #########################################################################
# PART 1: SIMULATED DATA
# #########################################################################

cat("\n", strrep("=", 70), "\n")
cat("PART 1: SIMULATED DATA\n")
cat(strrep("=", 70), "\n")

## ---- Scenario grid -------------------------------------------------------
K_vals     <- if (SMOKE) c(3, 5) else c(3, 5, 7)
gen_models <- c("WST", "SST")
strengths  <- c("strong", "weak")

## Hyper / mixing methods to compare
osbm_methods <- c("default", "oracle_ruleA", "spectral_ruleA",
                   "calibrate_full", "hierarchical", "no_moves")

sim_grid <- expand.grid(
  gen_model = gen_models,
  K_true    = K_vals,
  strength  = strengths,
  rep       = seq_len(N_REP),
  stringsAsFactors = FALSE
)

cat(sprintf("Simulation grid: %d scenarios × %d reps = %d datasets\n",
            nrow(sim_grid) / N_REP, N_REP, nrow(sim_grid)))
cat(sprintf("Methods per dataset: %d OSBM variants + 1 DCSBM = %d fits\n",
            length(osbm_methods) * 2, length(osbm_methods) * 2 + 1))

set.seed(BASE_SEED)
sim_seeds <- sample.int(1e6, nrow(sim_grid))

all_sim_rows <- list()

for (i in seq_len(nrow(sim_grid))) {
  gm  <- sim_grid$gen_model[i]
  Kt  <- sim_grid$K_true[i]
  str <- sim_grid$strength[i]
  rep <- sim_grid$rep[i]
  seed_data <- sim_seeds[i]

  cat(sprintf("\n--- [%d/%d] gen=%s K=%d signal=%s rep=%d seed=%d ---\n",
              i, nrow(sim_grid), gm, Kt, str, rep, seed_data))

  ## ---- Generate data -----------------------------------------------------
  set.seed(seed_data)
  z_true  <- rep(1:Kt, length.out = N_SIM)
  kappa_t <- make_kappa(Kt, mean = 2.0, var = 1.0)
  eta_t   <- runif(N_SIM, 0.8, 1.2)

  if (gm == "WST") {
    psi_t <- make_psi_wst(Kt, strength = str)
    A <- as.matrix(simulate_osbm(N_SIM, Kt, z_true, eta_t, kappa_t, psi_t, "WST"))
  } else {
    psi_t <- make_psi_sst(Kt, strength = str)
    A <- as.matrix(simulate_osbm(N_SIM, Kt, z_true, eta_t, kappa_t, psi_t, "SST"))
  }

  cat(sprintf("  sum(A)=%d  density=%.3f\n", sum(A), sum(A > 0) / (N_SIM * (N_SIM - 1))))

  ## ---- OSBM methods: fit with BOTH gen_model regime AND cross-fit --------
  fit_models_osbm <- c("WST", "SST")

  for (method in osbm_methods) {
    for (fm in fit_models_osbm) {
      hypers <- setup_hypers(A, method, fm, Kt)
      seed_fit <- seed_data + match(method, osbm_methods) * 100L + match(fm, fit_models_osbm)

      cat(sprintf("  %s / fit=%s  b_kappa=%.1f gamma=%.3f moves=%s hier=%s ... ",
                  method, fm, hypers$b_kappa, hypers$gamma_gn,
                  hypers$use_mixing_moves, hypers$sample_b_kappa))
      flush.console()

      res <- run_osbm_fit(A, Kt, z_true, fm, hypers, seed_fit)

      cat(sprintf("K=%d ARI=%.3f LOOIC=%.0f (%.1fs)\n",
                  res$K_hat, res$ari, res$looic, res$time_sec))

      row <- cbind(
        data.frame(
          part = "sim", dataset = NA_character_,
          gen_model = gm, fit_model = fm, K_true = Kt,
          strength = str, rep = rep,
          method = method, seed = seed_data,
          stringsAsFactors = FALSE
        ),
        res
      )
      all_sim_rows[[length(all_sim_rows) + 1L]] <- row

      # Checkpoint after every single fit
      write.csv(dplyr::bind_rows(all_sim_rows),
                file.path(out_dir, "sim_results_streaming.csv"), row.names = FALSE)
    }
  }

  ## ---- DCSBM cross-fit ---------------------------------------------------
  gamma_dc <- choose_gamma_from_K_expected(N_SIM, K_expected = Kt,
                ordering_prior_mode = "equivalence_class")
  seed_dc <- seed_data + 999L

  cat(sprintf("  DCSBM  gamma=%.3f ... ", gamma_dc))
  flush.console()

  res_dc <- run_dcsbm_fit(A, Kt, z_true, gamma_dc, seed_dc)
  cat(sprintf("K=%d ARI=%.3f LOOIC=%.0f (%.1fs)\n",
              res_dc$K_hat, res_dc$ari, res_dc$looic, res_dc$time_sec))

  row_dc <- cbind(
    data.frame(
      part = "sim", dataset = NA_character_,
      gen_model = gm, fit_model = "DCSBM", K_true = Kt,
      strength = str, rep = rep,
      method = "dcsbm_varK", seed = seed_data,
      stringsAsFactors = FALSE
    ),
    res_dc
  )
  all_sim_rows[[length(all_sim_rows) + 1L]] <- row_dc

  # Checkpoint after DCSBM fit
  write.csv(dplyr::bind_rows(all_sim_rows),
            file.path(out_dir, "sim_results_streaming.csv"), row.names = FALSE)
}

sim_df <- dplyr::bind_rows(all_sim_rows)
write.csv(sim_df, file.path(out_dir, "sim_results.csv"), row.names = FALSE)
cat(sprintf("\nSimulated data results: %d rows saved\n", nrow(sim_df)))


# #########################################################################
# PART 2: REAL DATA
# #########################################################################

cat("\n", strrep("=", 70), "\n")
cat("PART 2: REAL DATA\n")
cat(strrep("=", 70), "\n")

datasets_real <- c("mountain_goats", "citations_data", "macaques_data", "high_school")
if (SMOKE) datasets_real <- c("mountain_goats", "macaques_data")

## Methods for real data (no z_true, no oracle, no "no_moves" needed)
real_methods <- c("default", "spectral_ruleA", "calibrate_full", "hierarchical")

all_real_rows <- list()

for (ds in datasets_real) {
  cat(sprintf("\n======== Dataset: %s ========\n", ds))

  A <- if (ds == "mountain_goats") {
    fls <- list.files("./data/ShizukaMcDonald_Data", full.names = TRUE, pattern = "\\.csv$")
    nn <- vapply(fls, function(f) nrow(read.csv(f, row.names = 1)), integer(1))
    as.matrix(read.csv(fls[which.max(nn)], row.names = 1, check.names = FALSE))
  } else if (ds == "citations_data") {
    tmp <- as.matrix(read.csv("./data/Citations_application/cross-citation-matrix.csv",
                              row.names = 1, header = TRUE, check.names = FALSE))
    diag(tmp) <- 0; tmp
  } else if (ds == "macaques_data") {
    el <- read.table("./data/macaques/out.moreno.txt")
    nodes <- sort(unique(c(el[[1]], el[[2]])))
    tmp <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (r in seq_len(nrow(el))) tmp[el[r, 1], el[r, 2]] <- el[r, "V3"]
    tmp
  } else if (ds == "high_school") {
    edges <- read.csv("./data/high-school/edges.csv", header = FALSE,
                      comment.char = "#", strip.white = TRUE)
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    nn <- max(c(edges$source, edges$target))
    tmp <- matrix(0L, nn, nn)
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]; j <- edges$target[r]; w <- edges$weight[r]
      if (w > 0L) tmp[i, j] <- tmp[i, j] + w
    }
    diag(tmp) <- 0L; tmp
  } else stop("Unknown dataset: ", ds)

  A <- as.matrix(A)
  n <- nrow(A)
  cat(sprintf("  n=%d  sum(A)=%d  density=%.3f\n",
              n, sum(A), sum(A > 0) / (n * (n - 1))))

  K_spec <- estimate_K_spectral(A)
  cat(sprintf("  K_spectral=%d\n", K_spec))

  ## ---- OSBM methods × fit_model ∈ {WST, SST} ----------------------------
  for (method in real_methods) {
    for (fm in c("WST", "SST")) {
      hypers <- setup_hypers(A, method, fm, K_true = max(K_spec, 2))
      seed_fit <- BASE_SEED + match(ds, datasets_real) * 1000L +
                  match(method, real_methods) * 10L + match(fm, c("WST", "SST"))

      cat(sprintf("  %s / fit=%s  b=%.1f gamma=%.3f hier=%s ... ",
                  method, fm, hypers$b_kappa, hypers$gamma_gn,
                  hypers$sample_b_kappa))
      flush.console()

      res <- run_osbm_fit(A, K_true = NA, z_true = NULL, fm, hypers, seed_fit)

      cat(sprintf("K=%d LOOIC=%.0f ESS_K=%.0f (%.1fs)\n",
                  res$K_hat, res$looic, res$ess_K, res$time_sec))

      row <- cbind(
        data.frame(
          part = "real", dataset = ds,
          gen_model = NA_character_, fit_model = fm,
          K_true = NA_integer_,
          strength = NA_character_, rep = NA_integer_,
          method = method, seed = seed_fit,
          stringsAsFactors = FALSE
        ),
        res
      )
      all_real_rows[[length(all_real_rows) + 1L]] <- row

      # Checkpoint after every single fit
      write.csv(dplyr::bind_rows(all_real_rows),
                file.path(out_dir, "real_results_streaming.csv"), row.names = FALSE)
    }
  }

  ## ---- DCSBM fit ---------------------------------------------------------
  gamma_dc <- choose_gamma_from_K_expected(n, K_expected = max(K_spec, 2),
                ordering_prior_mode = "equivalence_class")
  seed_dc <- BASE_SEED + match(ds, datasets_real) * 1000L + 500L

  cat(sprintf("  DCSBM  gamma=%.3f ... ", gamma_dc))
  flush.console()

  res_dc <- run_dcsbm_fit(A, K_true = NA, z_true = NULL, gamma_dc, seed_dc)

  cat(sprintf("K=%d LOOIC=%.0f ESS_K=%.0f (%.1fs)\n",
              res_dc$K_hat, res_dc$looic, res_dc$ess_K, res_dc$time_sec))

  row_dc <- cbind(
    data.frame(
      part = "real", dataset = ds,
      gen_model = NA_character_, fit_model = "DCSBM",
      K_true = NA_integer_,
      strength = NA_character_, rep = NA_integer_,
      method = "dcsbm_varK", seed = seed_dc,
      stringsAsFactors = FALSE
    ),
    res_dc
  )
  all_real_rows[[length(all_real_rows) + 1L]] <- row_dc

  # Checkpoint after DCSBM fit
  write.csv(dplyr::bind_rows(all_real_rows),
            file.path(out_dir, "real_results_streaming.csv"), row.names = FALSE)
}

real_df <- dplyr::bind_rows(all_real_rows)
write.csv(real_df, file.path(out_dir, "real_results.csv"), row.names = FALSE)
cat(sprintf("\nReal data results: %d rows saved\n", nrow(real_df)))


# #########################################################################
# PART 3: SUMMARY TABLES
# #########################################################################

cat("\n", strrep("=", 70), "\n")
cat("SUMMARY TABLES\n")
cat(strrep("=", 70), "\n")

## ---- A. K recovery (simulated data) --------------------------------------
if (nrow(sim_df) > 0) {
  cat("\n--- A. K RECOVERY (simulated, correct-model fit) ---\n")
  k_recovery <- sim_df %>%
    filter(gen_model == fit_model | fit_model == "DCSBM") %>%
    group_by(gen_model, K_true, strength, method) %>%
    summarise(
      K_mean = mean(K_hat, na.rm = TRUE),
      K_mode_mean = mean(K_post_mode, na.rm = TRUE),
      K_within1 = mean(abs(K_hat - K_true) <= 1, na.rm = TRUE),
      ARI_mean = mean(ari, na.rm = TRUE),
      VI_mean  = mean(vi, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(gen_model, K_true, strength, method)

  print(as.data.frame(k_recovery), row.names = FALSE)
  write.csv(k_recovery, file.path(out_dir, "summary_K_recovery.csv"), row.names = FALSE)

  ## ---- B. Mixing diagnostics -----------------------------------------------
  cat("\n--- B. MIXING (ESS on K trace) ---\n")
  mixing_table <- sim_df %>%
    filter(gen_model == fit_model | fit_model == "DCSBM") %>%
    group_by(gen_model, K_true, strength, method) %>%
    summarise(
      ESS_K_mean = mean(ess_K, na.rm = TRUE),
      ESS_K_min  = min(ess_K, na.rm = TRUE),
      time_mean  = mean(time_sec, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(gen_model, K_true, strength, method)

  print(as.data.frame(mixing_table), row.names = FALSE)
  write.csv(mixing_table, file.path(out_dir, "summary_mixing.csv"), row.names = FALSE)

  ## ---- C. Model choice (cross-fit LOOIC) ----------------------------------
  cat("\n--- C. MODEL CHOICE (cross-fit LOOIC, oracle method only) ---\n")
  model_choice <- sim_df %>%
    filter(method %in% c("oracle_ruleA", "dcsbm_varK")) %>%
    group_by(gen_model, K_true, strength, fit_model) %>%
    summarise(
      LOOIC_mean = mean(looic, na.rm = TRUE),
      LOOIC_sd   = sd(looic, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(gen_model, K_true, strength, LOOIC_mean)

  print(as.data.frame(model_choice), row.names = FALSE)
  write.csv(model_choice, file.path(out_dir, "summary_model_choice.csv"), row.names = FALSE)

  ## ---- D. Hierarchical b_kappa learning -----------------------------------
  cat("\n--- D. HIERARCHICAL b_kappa FINAL VALUES ---\n")
  bk_table <- sim_df %>%
    filter(method %in% c("hierarchical", "oracle_ruleA", "default")) %>%
    group_by(gen_model, K_true, strength, method) %>%
    summarise(
      bk_mean  = mean(bk_final, na.rm = TRUE),
      bk_sd    = sd(bk_final, na.rm = TRUE),
      K_mean   = mean(K_hat, na.rm = TRUE),
      ARI_mean = mean(ari, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(gen_model, K_true, strength, method)

  print(as.data.frame(bk_table), row.names = FALSE)
  write.csv(bk_table, file.path(out_dir, "summary_bkappa.csv"), row.names = FALSE)
}

## ---- E. Real data comparison ----------------------------------------------
if (nrow(real_df) > 0) {
  cat("\n--- E. REAL DATA ---\n")
  real_summary <- real_df %>%
    group_by(dataset, fit_model, method) %>%
    summarise(
      K_hat  = mean(K_hat, na.rm = TRUE),
      LOOIC  = mean(looic, na.rm = TRUE),
      ESS_K  = mean(ess_K, na.rm = TRUE),
      bk     = mean(bk_final, na.rm = TRUE),
      time_s = mean(time_sec, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(dataset, fit_model, method)

  print(as.data.frame(real_summary), row.names = FALSE)
  write.csv(real_summary, file.path(out_dir, "summary_real.csv"), row.names = FALSE)

  # Best model per dataset by LOOIC
  cat("\n--- Best model per dataset (lowest LOOIC) ---\n")
  best_model <- real_df %>%
    group_by(dataset) %>%
    filter(!is.na(looic)) %>%
    slice_min(looic, n = 1) %>%
    select(dataset, fit_model, method, K_hat, looic) %>%
    ungroup()
  print(as.data.frame(best_model), row.names = FALSE)
}

# ---- Combined output -------------------------------------------------------
full_df <- dplyr::bind_rows(sim_df, real_df)
write.csv(full_df, file.path(out_dir, "all_results.csv"), row.names = FALSE)
saveRDS(full_df, file.path(out_dir, "all_results.rds"))

cat(sprintf("\n=== Study complete: %s ===\n", Sys.time()))
cat(sprintf("Output directory: %s\n", out_dir))
cat(sprintf("Total results: %d rows (%d sim + %d real)\n",
            nrow(full_df), nrow(sim_df), nrow(real_df)))

sink(type = "message")
sink(type = "output")
close(log_con)
