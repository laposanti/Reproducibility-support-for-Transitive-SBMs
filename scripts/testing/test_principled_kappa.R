#!/usr/bin/env Rscript
# Quick end-to-end test of the principled kappa calibration:
#   a_kappa = 1,  b_kappa = c · n / K_expected
# Runs 100 MCMC iters per dataset and shows K trace.

suppressPackageStartupMessages({
  library(Rcpp)
  library(truncnorm)
})

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")
source("core/my_best_try_so_far.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/Hyper_setup.R")

# --- Dataset loader (same as application.R choose_dataset()) ---
load_ds <- function(name) {
  switch(name,
    mountain_goats = {
      files <- list.files("./data/ShizukaMcDonald_Data",
                          full.names = TRUE, pattern = "\\.csv$")
      n_each <- vapply(files, function(f) nrow(read.csv(f, row.names = 1)),
                        FUN.VALUE = integer(1))
      as.matrix(read.csv(files[which.max(n_each)],
                          row.names = 1, check.names = FALSE))
    },
    citations_data = {
      A <- as.matrix(read.csv("./data/Citations_application/cross-citation-matrix.csv",
                               row.names = 1, header = TRUE, check.names = FALSE))
      diag(A) <- 0; A
    },
    macaques_data = {
      el <- read.table("./data/macaques/out.moreno.txt")
      nodes <- sort(unique(c(el[[1]], el[[2]])))
      A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
      for (i in seq_len(nrow(el)))
        A[el[i,1], el[i,2]] <- el[i, "V3"]
      A
    },
    high_school = {
      edges <- read.csv("./data/high-school/edges.csv",
                         header = FALSE, comment.char = "#", strip.white = TRUE)
      names(edges)[1:3] <- c("source","target","weight")
      edges$source <- as.integer(edges$source)
      edges$target <- as.integer(edges$target)
      edges$weight <- as.integer(edges$weight)
      if (min(edges$source, edges$target) == 0L) {
        edges$source <- edges$source + 1L
        edges$target <- edges$target + 1L
      }
      nn <- max(c(edges$source, edges$target))
      A <- matrix(0L, nn, nn)
      for (i in seq_len(nrow(edges)))
        A[edges$source[i], edges$target[i]] <-
          A[edges$source[i], edges$target[i]] + edges$weight[i]
      A
    }
  )
}

# ====================================================================
cat("=== Principled kappa calibration: sampler test ===\n\n")
cat("Rule: a_kappa = 1, b_kappa = c · n / K_expected\n\n")

datasets  <- c("mountain_goats", "citations_data", "macaques_data", "high_school")
K_exp     <- 5L
c_values  <- c(1, 3, 5, 10)   # test several concentrations
n_iter    <- 100L
burn      <- 20L

results <- data.frame()

for (ds in datasets) {
  A <- load_ds(ds)
  n <- nrow(A)
  N_mat <- A + t(A); diag(N_mat) <- 0
  mean_N <- mean(N_mat[upper.tri(N_mat)])

  # Get principled hypers (direction, eta, gamma from existing calibration)
  hyp <- get_principled_hypers_v2(A, K_expected = K_exp, c_kappa = 5)

  cat(sprintf("--- %s (n=%d, mean_N=%.2f) ---\n", ds, n, mean_N))
  cat(sprintf("  gamma=%.3f, tau0=%.3f, a_eta=%.2f\n",
              hyp$gamma, hyp$tau0, hyp$a_eta))

  for (cc in c_values) {
    b_val <- cc * n / K_exp
    cat(sprintf("  c=%2d  b_kappa=%6.1f  ...", cc, b_val))

    out <- tryCatch(
      modular_osbm_sampler(
        A = A, K = K_exp, truth = NA,
        free = c("psi","kappa","eta","z"),
        n_iter = n_iter, burn = burn, thin = 1,
        verbose = FALSE,
        psi_constraint = "SST",
        seed = 42 ,
        a_kappa = 1, b_kappa = b_val,
        a_eta = hyp$a_eta, b_eta = hyp$b_eta,
        mu0 = hyp$mu0, sigma0 = hyp$sigma0,
        tau0 = hyp$tau0,
        gamma_gn = hyp$gamma
      ),
      error = function(e) list(error = conditionMessage(e))
    )

    if (!is.null(out$error)) {
      cat(sprintf(" ERROR: %s\n", out$error))
      next
    }

    K_trace <- out$K_trace
    K_trace <- K_trace[K_trace > 0]  # drop unused slots
    K_final <- tail(K_trace, 1)
    K_range <- range(K_trace)
    post_burn <- K_trace[seq(max(1, burn+1), length(K_trace))]
    K_med   <- median(post_burn)

    cat(sprintf(" K_final=%d, K_median=%.0f, K_range=[%d,%d]\n",
                K_final, K_med, K_range[1], K_range[2]))

    results <- rbind(results, data.frame(
      dataset = ds, n = n, mean_N = mean_N,
      c = cc, b_kappa = b_val, a_kappa = 1,
      K_final = K_final, K_median = K_med,
      K_min = K_range[1], K_max = K_range[2],
      stringsAsFactors = FALSE
    ))
  }
  cat("\n")
}

# ====================================================================
cat("\n=== SUMMARY TABLE ===\n\n")
cat(sprintf("%15s  %4s  %8s  %8s  %8s  %8s  %8s\n",
            "dataset", "c", "b_kappa", "K_final", "K_med", "K_min", "K_max"))
for (i in seq_len(nrow(results))) {
  r <- results[i,]
  cat(sprintf("%15s  %4d  %8.1f  %8d  %8.0f  %8d  %8d\n",
              r$dataset, r$c, r$b_kappa, r$K_final, r$K_median, r$K_min, r$K_max))
}

# Comparison with baseline (a=1, b=1)
cat("\n=== BASELINE COMPARISON (a=1, b=1) ===\n\n")
for (ds in datasets) {
  A <- load_ds(ds)
  n <- nrow(A)
  hyp <- get_principled_hypers_v2(A, K_expected = K_exp, c_kappa = 5)
  cat(sprintf("  %s (b=1): ", ds))
  out <- tryCatch(
    modular_osbm_sampler(
      A = A, K = K_exp, truth = NA,
      free = c("psi","kappa","eta","z"),
      n_iter = n_iter, burn = burn, thin = 1,
      verbose = FALSE, psi_constraint = "SST", seed = 42 ,
      a_kappa = 1, b_kappa = 1,
      a_eta = hyp$a_eta, b_eta = hyp$b_eta,
      mu0 = hyp$mu0, sigma0 = hyp$sigma0,
      tau0 = hyp$tau0, gamma_gn = hyp$gamma
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  if (!is.null(out$error)) { cat(sprintf("ERROR: %s\n", out$error)); next }
  K_trace <- out$K_trace
  K_trace <- K_trace[K_trace > 0]
  post_burn <- K_trace[seq(max(1,burn+1), length(K_trace))]
  cat(sprintf("K_final=%d, K_med=%.0f, K_range=[%d,%d]\n",
              tail(K_trace,1), median(post_burn),
              min(K_trace), max(K_trace)))
}

cat("\nDone.\n")
