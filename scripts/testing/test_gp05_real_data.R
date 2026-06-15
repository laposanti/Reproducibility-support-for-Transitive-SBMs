#!/usr/bin/env Rscript
# ======================================================================
# Test GP05 prior on real datasets: mountain goats & citations
#
# Sweeps (alpha_gp05, gamma_gn) to verify they control K.
# Short runs (1000 iter, 300 burn) for quick turnaround.
# ======================================================================
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

cat("=== GP05 prior: real-data K-control test ===\n")
cat("Start:", format(Sys.time()), "\n\n")
t0_all <- proc.time()[3]

suppressPackageStartupMessages({
  library(Matrix)
  library(BayesLogit)
  library(salso)
  library(fossil)
})

source("helper_folder/helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("helper_folder/sim_study_helper.R")
source("core/my_best_try_so_far.R")

# ======================================================================
# Load datasets
# ======================================================================
cat("Loading datasets...\n")

# Mountain goats
matrix_files <- list.files("./data/ShizukaMcDonald_Data",
                           full.names = TRUE, pattern = "[.]csv$")
n_each <- vapply(matrix_files, function(f) nrow(read.csv(f, row.names = 1)),
                 FUN.VALUE = integer(1))
A_goats <- as.matrix(read.csv(matrix_files[which.max(n_each)],
                               row.names = 1, check.names = FALSE))
cat(sprintf("  Mountain goats: n=%d, sum(A)=%d, density=%.3f\n",
            nrow(A_goats), sum(A_goats),
            sum(A_goats > 0) / (nrow(A_goats) * (nrow(A_goats) - 1))))

# Citations
A_cit <- as.matrix(read.csv("./data/Citations_application/cross-citation-matrix.csv",
                             row.names = 1, header = TRUE, check.names = FALSE))
diag(A_cit) <- 0
cat(sprintf("  Citations:      n=%d, sum(A)=%d, density=%.3f\n",
            nrow(A_cit), sum(A_cit),
            sum(A_cit > 0) / (nrow(A_cit) * (nrow(A_cit) - 1))))

# ======================================================================
# Hyperparameter grid
# ======================================================================
# gamma controls new-block penalty; alpha controls composition structure
# Expect: small gamma => fewer blocks, large gamma => more blocks
#         alpha has a secondary effect
hyper_grid <- expand.grid(
  alpha_gp05 = c(0.3, 0.5, 0.7),
  gamma_gn   = c(0.3, 0.6, 0.9),
  stringsAsFactors = FALSE
)

n_iter <- 1000
burn   <- 300
thin   <- 1
seed   <- 42

# ======================================================================
# Run function
# ======================================================================
run_one <- function(A, dataset_name, alpha_gp05, gamma_gn, regime = "WST") {
  n <- nrow(A)
  cat(sprintf("  [%s] alpha=%.1f gamma=%.1f %s ... ",
              dataset_name, alpha_gp05, gamma_gn, regime))
  t0 <- proc.time()[3]

  out <- tryCatch(
    modular_osbm_sampler(
      A = A, K = n,
      n_iter = n_iter, burn = burn, thin = thin,
      psi_constraint = regime ,
      partition_prior = "GP05",
      gamma_gn   = gamma_gn,
      alpha_gp05 = alpha_gp05,
      a_kappa = 0.1, b_kappa = 2,
      a_eta   = 0.1, b_eta   = 1,
      mu0 = 0.1, sigma0 = 1, tau0 = 1,
      use_mixing_moves = TRUE,
      seed = seed, verbose = FALSE
    ),
    error = function(e) {
      cat(sprintf("ERROR: %s\n", e$message))
      return(NULL)
    }
  )

  elapsed <- proc.time()[3] - t0

  if (is.null(out)) {
    return(data.frame(
      dataset = dataset_name, regime = regime,
      alpha = alpha_gp05, gamma = gamma_gn,
      K_mode = NA, K_mean = NA, K_median = NA,
      K_min = NA, K_max = NA, time_s = elapsed,
      stringsAsFactors = FALSE
    ))
  }

  z_chain <- do.call(rbind, lapply(out$z, as.integer))
  K_trace <- apply(z_chain, 1, function(x) length(unique(x)))
  K_mode  <- as.integer(names(which.max(table(K_trace))))

  cat(sprintf("K_mode=%d K_mean=%.1f [%d,%d] (%.0fs)\n",
              K_mode, mean(K_trace), min(K_trace), max(K_trace), elapsed))

  data.frame(
    dataset = dataset_name, regime = regime,
    alpha = alpha_gp05, gamma = gamma_gn,
    K_mode = K_mode, K_mean = round(mean(K_trace), 2),
    K_median = median(K_trace),
    K_min = min(K_trace), K_max = max(K_trace),
    time_s = round(elapsed, 1),
    stringsAsFactors = FALSE
  )
}

# ======================================================================
# Run sweep
# ======================================================================
results <- list()
idx <- 1L

for (ds_name in c("mountain_goats", "citations_data")) {
  A <- if (ds_name == "mountain_goats") A_goats else A_cit
  cat(sprintf("\n--- %s (n=%d) ---\n", ds_name, nrow(A)))

  for (i in seq_len(nrow(hyper_grid))) {
    res <- run_one(
      A, ds_name,
      alpha_gp05 = hyper_grid$alpha_gp05[i],
      gamma_gn   = hyper_grid$gamma_gn[i],
      regime = "WST"
    )
    results[[idx]] <- res
    idx <- idx + 1L
  }
}

# ======================================================================
# Summary
# ======================================================================
tab <- do.call(rbind, results)

cat("\n\n========================================\n")
cat("     GP05 PRIOR: K-CONTROL SUMMARY\n")
cat("========================================\n\n")

for (ds in unique(tab$dataset)) {
  sub <- tab[tab$dataset == ds, ]
  cat(sprintf("Dataset: %s\n", ds))
  cat(sprintf("  %-8s %-8s  %6s  %6s  %8s  %6s\n",
              "alpha", "gamma", "K_mode", "K_mean", "K_range", "time"))
  cat(sprintf("  %s\n", paste(rep("-", 52), collapse = "")))
  for (r in seq_len(nrow(sub))) {
    cat(sprintf("  %-8.1f %-8.1f  %6d  %6.1f  [%2d, %2d]  %5.0fs\n",
                sub$alpha[r], sub$gamma[r],
                sub$K_mode[r], sub$K_mean[r],
                sub$K_min[r], sub$K_max[r],
                sub$time_s[r]))
  }
  cat("\n")
}

# Check monotonicity: does gamma increasing => K increasing?
cat("--- Monotonicity check (gamma -> K) ---\n")
for (ds in unique(tab$dataset)) {
  for (a in unique(tab$alpha)) {
    sub <- tab[tab$dataset == ds & tab$alpha == a, ]
    sub <- sub[order(sub$gamma), ]
    K_means <- sub$K_mean
    monotone <- all(diff(K_means) >= -0.5)  # allow slight noise
    cat(sprintf("  %s alpha=%.1f: gamma=%s => K_mean=%s  %s\n",
                ds, a,
                paste(sub$gamma, collapse=" -> "),
                paste(K_means, collapse=" -> "),
                ifelse(monotone, "OK", "NOT MONOTONE")))
  }
}

total_time <- proc.time()[3] - t0_all
cat(sprintf("\nTotal time: %.0fs (%.1f min)\n", total_time, total_time / 60))
cat("=== Done ===\n")
