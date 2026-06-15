#!/usr/bin/env Rscript
# ======================================================================
# Hyperparameter sensitivity study
#
# Goal: Understand how (a_kappa, b_kappa, a_eta, b_eta) affect
#       cluster creation across 4 real datasets. Compare with DC-SBM.
#
# Approach:
#   1. Load each dataset, compute network summaries
#   2. Vary c_kappa in {1, 3, 5, 10, 20} (holding a_kappa=1)
#   3. Test hierarchical b_kappa (sample_b_kappa=TRUE) as alternative
#   4. Test a_eta matched to CV(deg) vs fixed a_eta=2
#   5. Run DC-SBM for comparison
#   6. Save and display results
# ======================================================================
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

# --- Dedicated output folder with timestamp ---
run_tag  <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir  <- file.path("output", "hyperparameter_sensitivity", run_tag)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

# --- Tee all output to a log file ---
log_path <- file.path(run_dir, "run_log.txt")
log_con  <- file(log_path, open = "wt")
sink(log_con, type = "output",  split = TRUE)
sink(log_con, type = "message", append = TRUE)

cat("=== Hyperparameter Sensitivity Study ===\n")
cat("Start:", format(Sys.time()), "\n")
cat("Log:  ", log_path, "\n\n")

suppressPackageStartupMessages({
  library(coda)
  library(mcclust)
  library(mcclust.ext)
  library(salso)
  library(dplyr)
})

source("./helper_folder/helper.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./helper_folder/Hyper_setup.R")
source("./core/my_best_try_so_far.R")
source("./core/DCSBM_varK.R")

# ======================================================================
# Data loading (from application code)
# ======================================================================
choose_dataset <- function(dataset = c("mountain_goats","citations_data",
                                       "macaques_data","high_school")) {
  dataset <- match.arg(dataset)
  if (dataset == "mountain_goats") {
    matrix_files <- list.files("./data/ShizukaMcDonald_Data",
                               full.names = TRUE, pattern = "\\.csv$")
    n_each <- vapply(matrix_files, function(f) nrow(read.csv(f, row.names = 1)),
                     FUN.VALUE = integer(1))
    A <- as.matrix(read.csv(matrix_files[which.max(n_each)],
                            row.names = 1, check.names = FALSE))
  } else if (dataset == "citations_data") {
    A <- read.csv("./data/Citations_application/cross-citation-matrix.csv",
                  row.names = 1, header = TRUE, check.names = FALSE)
    diag(A) <- 0
  } else if (dataset == "macaques_data") {
    edge_list <- read.table("./data/macaques/out.moreno.txt")
    nodes <- sort(unique(c(edge_list[[1]], edge_list[[2]])))
    A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (i in seq_len(nrow(edge_list))) {
      A[edge_list[i,1], edge_list[i,2]] <- edge_list[i,"V3"]
    }
  } else if (dataset == "high_school") {
    edges <- read.csv("./data/high-school/edges.csv",
                      header = FALSE, comment.char = "#", strip.white = TRUE)
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    n_nodes <- max(c(edges$source, edges$target))
    node_ids <- as.character(seq_len(n_nodes))
    A <- matrix(0L, n_nodes, n_nodes, dimnames = list(node_ids, node_ids))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]; j <- edges$target[r]; w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w
    }
    diag(A) <- 0L
  }
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  colnames(A) <- rownames(A)
  A
}

# ======================================================================
# Network summaries
# ======================================================================
network_summary <- function(A) {
  n <- nrow(A)
  N <- A + t(A); diag(N) <- 0
  edges_total <- sum(N[upper.tri(N)])
  density <- edges_total / (n*(n-1))
  deg <- rowSums(A) + colSums(A)
  cv_deg <- sd(deg)/mean(deg)
  K_spectral <- estimate_K_spectral(A)
  cat(sprintf("  n=%d, edges=%d, density=%.4f, CV(deg)=%.3f, K_spectral=%d\n",
              n, edges_total, density, cv_deg, K_spectral))
  list(n = n, density = density, cv_deg = cv_deg, K_spectral = K_spectral,
       mean_deg = mean(deg), edges_total = edges_total)
}

# ======================================================================
# Run one OSBM fit and extract summary
# ======================================================================
run_osbm_quick <- function(A, psi_mode, a_kappa, b_kappa, a_eta, b_eta,
                           gamma_gn, mu0, sigma0, tau0,
                           sample_b_kappa = FALSE,
                           n_iter = 3000, burn = 1000, thin = 1,
                           seed = 42, label = "osbm",
                           K_init = NULL) {
  n <- nrow(A)
  if (is.null(K_init)) K_init <- min(n, max(5L, 2L * estimate_K_spectral(A)))
  t0 <- proc.time()
  out <- tryCatch({
    modular_osbm_sampler(
      A = A, K = K_init,
      n_iter = n_iter, burn = burn, thin = thin,
      verbose = TRUE,
      psi_constraint = psi_mode ,
      partition_prior = "GN",
      gamma_gn = gamma_gn,
      a_kappa = a_kappa, b_kappa = b_kappa,
      a_eta = a_eta, b_eta = b_eta,
      mu0 = mu0, sigma0 = sigma0, tau0 = tau0,
      use_mixing_moves = TRUE,
      sample_b_kappa = sample_b_kappa,
      seed = seed
    )
  }, error = function(e) { message("  ERROR: ", e$message); NULL })
  elapsed <- (proc.time() - t0)[["elapsed"]]

  if (is.null(out)) {
    return(data.frame(label = label, K_mean = NA, K_median = NA,
                      K_mode = NA, K_lo = NA, K_hi = NA,
                      b_kappa_final = b_kappa, elapsed = elapsed,
                      stringsAsFactors = FALSE))
  }

  K_post <- out$K_trace
  K_post <- K_post[!is.na(K_post)]

  # b_kappa trace if hierarchical
  b_final <- if (sample_b_kappa && !is.null(out$b_kappa_trace)) {
    tail_b <- out$b_kappa_trace[!is.na(out$b_kappa_trace)]
    if (length(tail_b) > 0) round(median(tail_b), 2) else b_kappa
  } else b_kappa

  data.frame(
    label    = label,
    K_mean   = round(mean(K_post), 2),
    K_median = round(median(K_post), 0),
    K_mode   = as.integer(names(sort(table(K_post), decreasing = TRUE))[1]),
    K_lo     = quantile(K_post, 0.025),
    K_hi     = quantile(K_post, 0.975),
    b_kappa_final = round(b_final, 2),
    elapsed  = round(elapsed, 1),
    stringsAsFactors = FALSE
  )
}

# ======================================================================
# Run DC-SBM and extract summary
# ======================================================================
run_dcsbm_quick <- function(A, n_iter = 2000, burn = 500, thin = 1,
                            seed = 42) {
  n <- nrow(A)
  K_init <- min(n, max(10L, 2L * estimate_K_spectral(A)))
  t0 <- proc.time()
  out <- tryCatch({
    fit_dcsbm_gibbs_gnedin(
      A = as.matrix(A), K_init = K_init,
      priors = list(a_eta = 1, b_eta = 1, a_lambda = 1, b_lambda = 1,
                    gamma_gnedin = 0.95),
      iters = n_iter, burn_in = burn, thin = thin,
      verbose = 50, seed = seed
    )
  }, error = function(e) { message("  DCSBM ERROR: ", e$message); NULL })
  elapsed <- (proc.time() - t0)[["elapsed"]]

  if (is.null(out)) {
    return(data.frame(label = "DCSBM", K_mean = NA, K_median = NA,
                      K_mode = NA, K_lo = NA, K_hi = NA,
                      b_kappa_final = NA, elapsed = elapsed,
                      stringsAsFactors = FALSE))
  }

  K_post <- out$K_trace
  K_post <- K_post[!is.na(K_post)]

  data.frame(
    label    = "DCSBM",
    K_mean   = round(mean(K_post), 2),
    K_median = round(median(K_post), 0),
    K_mode   = as.integer(names(sort(table(K_post), decreasing = TRUE))[1]),
    K_lo     = quantile(K_post, 0.025),
    K_hi     = quantile(K_post, 0.975),
    b_kappa_final = NA,
    elapsed  = round(elapsed, 1),
    stringsAsFactors = FALSE
  )
}

# ======================================================================
# Main study
# ======================================================================
datasets <- c("mountain_goats", "macaques_data")
c_kappas <- c(1, 5, 20)

n_iter <- 1000
burn   <- 200
thin   <- 1
seed   <- 42

all_results <- list()

for (ds in datasets) {
  cat("\n========================================\n")
  cat("Dataset:", ds, "\n")
  cat("========================================\n")

  tryCatch({

  A <- choose_dataset(ds)
  s <- network_summary(A)
  n <- s$n

  # Use SST as the main model (consistent with production)
  psi_mode <- "SST"

  # Calibrate directional + eta using the data-driven calibrator
  cal <- calibrate_osbm_hypers(A, regime = psi_mode)
  K_exp <- cal$meta$K_expected
  gamma_gn <- cal$gamma
  mu0     <- cal$mu0
  sigma0  <- cal$sigma0
  tau0    <- cal$tau0

  # Eta: data-matched
  a_eta_matched <- cal$a_eta
  b_eta_matched <- cal$b_eta

  cat(sprintf("  K_expected=%d, gamma=%.4f, a_eta_matched=%.3f\n",
              K_exp, gamma_gn, a_eta_matched))
  cat(sprintf("  mu0=%.3f, sigma0=%.3f, tau0=%.3f\n", mu0, sigma0, tau0))

  results_ds <- list()

  # --- Sweep over c_kappa ---
  for (cc in c_kappas) {
    b_k <- cc * n / K_exp
    lab <- sprintf("c_kappa=%d", cc)
    cat(sprintf("  Running %s (b_kappa=%.1f) ... ", lab, b_k))
    res <- run_osbm_quick(
      A, psi_mode,
      a_kappa = 1, b_kappa = b_k,
      a_eta = a_eta_matched, b_eta = b_eta_matched,
      gamma_gn = gamma_gn, mu0 = mu0, sigma0 = sigma0, tau0 = tau0,
      n_iter = n_iter, burn = burn, thin = thin,
      seed = seed, label = lab
    )
    res$dataset <- ds
    cat(sprintf("K_mode=%s, K_mean=%.1f [%.0f, %.0f] in %.0fs\n",
                res$K_mode, res$K_mean, res$K_lo, res$K_hi, res$elapsed))
    results_ds[[length(results_ds) + 1]] <- res
  }

  # --- Hierarchical b_kappa ---
  cat("  Running hierarchical b_kappa ... ")
  b_init <- 5 * n / K_exp  # start at c=5
  res <- run_osbm_quick(
    A, psi_mode,
    a_kappa = 1, b_kappa = b_init,
    a_eta = a_eta_matched, b_eta = b_eta_matched,
    gamma_gn = gamma_gn, mu0 = mu0, sigma0 = sigma0, tau0 = tau0,
    sample_b_kappa = TRUE,
    n_iter = n_iter, burn = burn, thin = thin,
    seed = seed, label = "hierarchical_b"
  )
  res$dataset <- ds
  cat(sprintf("K_mode=%s, K_mean=%.1f, b_final=%.1f in %.0fs\n",
              res$K_mode, res$K_mean, res$b_kappa_final, res$elapsed))
  results_ds[[length(results_ds) + 1]] <- res

  # --- Fixed a_eta = 2 (non-adaptive) vs matched ---
  cat("  Running a_eta=2 (fixed) ... ")
  b_k <- 5 * n / K_exp
  res <- run_osbm_quick(
    A, psi_mode,
    a_kappa = 1, b_kappa = b_k,
    a_eta = 2, b_eta = 2,
    gamma_gn = gamma_gn, mu0 = mu0, sigma0 = sigma0, tau0 = tau0,
    n_iter = n_iter, burn = burn, thin = thin,
    seed = seed, label = "a_eta=2_fixed"
  )
  res$dataset <- ds
  cat(sprintf("K_mode=%s, K_mean=%.1f in %.0fs\n",
              res$K_mode, res$K_mean, res$elapsed))
  results_ds[[length(results_ds) + 1]] <- res

  # --- DC-SBM baseline ---
  cat("  Running DC-SBM ... ")
  res <- run_dcsbm_quick(A, n_iter = n_iter, burn = burn, thin = thin, seed = seed)
  res$dataset <- ds
  cat(sprintf("K_mode=%s, K_mean=%.1f in %.0fs\n",
              res$K_mode, res$K_mean, res$elapsed))
  results_ds[[length(results_ds) + 1]] <- res

  all_results[[ds]] <- do.call(rbind, results_ds)
  }, error = function(e) {
    cat(sprintf("\n  *** ERROR in dataset %s: %s ***\n", ds, e$message))
  })
}

# ======================================================================
# Print summary table
# ======================================================================
cat("\n\n")
cat("================================================================\n")
cat("RESULTS SUMMARY\n")
cat("================================================================\n\n")

final <- do.call(rbind, all_results)
rownames(final) <- NULL

for (ds in datasets) {
  cat(ds, ":\n")
  sub <- final[final$dataset == ds, c("label", "K_mode", "K_mean", "K_lo",
                                       "K_hi", "b_kappa_final", "elapsed")]
  print(sub, row.names = FALSE)
  cat("\n")
}

# Save results
out_file <- file.path(run_dir, "hyperparameter_sensitivity.rds")
saveRDS(final, out_file)
cat("Results saved to:", out_file, "\n")

# ======================================================================
# Diagnostic: birth probability decomposition for one dataset
# ======================================================================
cat("\n\n")
cat("================================================================\n")
cat("BIRTH PROBABILITY DECOMPOSITION (mountain_goats, node 1)\n")
cat("================================================================\n\n")

A <- choose_dataset("mountain_goats")
n <- nrow(A)
cal <- calibrate_osbm_hypers(A, regime = "SST")
K_exp <- cal$meta$K_expected

cat("Decomposing birth kappa score for c_kappa = {1, 5, 20}:\n\n")

for (cc in c(1, 5, 20)) {
  b_k <- cc * n / K_exp
  # Compute gp_marginal for a typical node
  # Node 1's edge counts to "all other nodes" (treating as one block for simplicity)
  N_mat <- A + t(A); diag(N_mat) <- 0
  R_i <- sum(N_mat[1, ])  # total edge count for node 1
  T_i <- n - 1            # exposure (eta=1 assumption)

  score_new <- gp_marginal(R_i, T_i, 1, b_k)
  score_0   <- gp_marginal(0, 0, 1, b_k)  # no edges (reference)

  cat(sprintf("  c=%2d  b_kappa=%8.1f  R_i=%d  T_i=%d  score_new=%.3f  -R*log(b)=%.3f\n",
              cc, b_k, R_i, T_i, score_new, -R_i * log(b_k)))
}

cat("\nDone.\n")
cat("End:", format(Sys.time()), "\n")

# Close log
sink(type = "message")
sink(type = "output")
close(log_con)
