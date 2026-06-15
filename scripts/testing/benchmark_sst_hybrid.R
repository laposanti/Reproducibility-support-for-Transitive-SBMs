#!/usr/bin/env Rscript
# =============================================================================
# Benchmark: SST z-update — collapsed vs full-augmented vs hybrid
# =============================================================================
# 3 variants:
#   (1) collapsed  – exact Binom-logistic kernel + numerical integration
#   (2) augmented  – full PG kernel (edge-level omega, quadratic for all)
#   (3) hybrid     – exact kernel for existing blocks & known distances;
#                    PG closed-form integral ONLY for the extreme distance
# =============================================================================

out_dir <- "output/simulation/benchmark_sst_hybrid"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(out_dir, paste0("benchmark_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat("=== Benchmark: SST collapsed vs augmented vs hybrid ===\n")
cat("Started:", format(Sys.time()), "\n")
t0_script <- proc.time()[3]

suppressPackageStartupMessages({
  library(Matrix)
  library(BayesLogit)
  library(mcclust)
})

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")

source("helper_folder/helper.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("helper_folder/mixing_moves.R")
Rcpp::sourceCpp("core/counts_by_block_exact_cpp.cpp")
Rcpp::sourceCpp("core/block_totals_for_poisson_cpp.cpp")
source("core/my_best_try_so_far.R")

cat("  Sources loaded.\n")

# =============================================================================
# SWEEP FUNCTIONS
# =============================================================================

sweep_z_sst_collapsed <- function(z, A, eta, kappa, psi,
                                  i_idx, j_idx, N_edge, edge_by_node,
                                  hyper, psi_hyper) {
  n <- length(z); K <- nrow(kappa)
  bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  slot_rad <- K + 1L
  for (i in sample.int(n)) {
    z_old_i <- z[i]; K_old <- K
    res <- sst_update_i_with_birth_LOO(
      i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
      Rkl = bt$Rkl, Tkl = bt$Tkl,
      i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
      gamma_gn = psi_hyper$gamma_gn,
      slot_radius = slot_rad,
      tau0 = psi_hyper$tau0, mu0 = psi_hyper$mu0,
      sig2_0 = psi_hyper$sigma0
    )
    z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K
    if (K != K_old || z[i] != z_old_i || isTRUE(res$changed))
      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  }
  list(z = z, kappa = kappa, psi = psi, K = nrow(kappa))
}

sweep_z_sst_augmented <- function(z, A, eta, kappa, psi,
                                  i_idx, j_idx, N_edge, edge_by_node,
                                  hyper, psi_hyper) {
  n <- length(z); K <- nrow(kappa)
  omega_edge <- draw_edge_omega(z, psi, i_idx, j_idx, N_edge, mode = "distance")
  bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  slot_rad <- K + 1L
  for (i in sample.int(n)) {
    z_old_i <- z[i]; K_old <- K
    res <- sst_update_i_with_birth_LOO_pg(
      i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
      Rkl = bt$Rkl, Tkl = bt$Tkl,
      i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
      gamma_gn = psi_hyper$gamma_gn,
      omega_edge = omega_edge,
      slot_radius = slot_rad,
      tau0 = psi_hyper$tau0, mu0 = psi_hyper$mu0,
      sig2_0 = psi_hyper$sigma0
    )
    z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K
    if (K != K_old || z[i] != z_old_i || isTRUE(res$changed))
      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  }
  list(z = z, kappa = kappa, psi = psi, K = nrow(kappa))
}

sweep_z_sst_hybrid <- function(z, A, eta, kappa, psi,
                               i_idx, j_idx, N_edge, edge_by_node,
                               hyper, psi_hyper) {
  n <- length(z); K <- nrow(kappa)
  bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  slot_rad <- K + 1L
  for (i in sample.int(n)) {
    z_old_i <- z[i]; K_old <- K
    res <- sst_update_i_with_birth_LOO_hybrid(
      i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
      Rkl = bt$Rkl, Tkl = bt$Tkl,
      i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
      gamma_gn = psi_hyper$gamma_gn,
      slot_radius = slot_rad,
      tau0 = psi_hyper$tau0, mu0 = psi_hyper$mu0,
      sig2_0 = psi_hyper$sigma0
    )
    z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K
    if (K != K_old || z[i] != z_old_i || isTRUE(res$changed))
      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
  }
  list(z = z, kappa = kappa, psi = psi, K = nrow(kappa))
}

# =============================================================================
# SHARED PARAMETER UPDATES
# =============================================================================

update_kappa_block <- function(z, eta, kappa, i_idx, j_idx, N_edge, hyper) {
  K <- nrow(kappa)
  zi <- z[i_idx]; zj <- z[j_idx]
  R <- matrix(0, K, K); T_mat <- matrix(0, K, K)
  for (e in seq_along(N_edge)) {
    p <- min(zi[e], zj[e]); q <- max(zi[e], zj[e])
    R[p, q] <- R[p, q] + N_edge[e]
    T_mat[p, q] <- T_mat[p, q] + eta[i_idx[e]] * eta[j_idx[e]]
  }
  for (k in seq_len(K)) for (l in k:K)
    kappa[k, l] <- kappa[l, k] <- rgamma(1,
      shape = hyper$a_kappa + R[k, l], rate = hyper$b_kappa + T_mat[k, l])
  kappa
}

update_eta_block <- function(z, eta, kappa, A, hyper) {
  n <- length(z); K <- nrow(kappa)
  G_i <- Matrix::rowSums(A) + Matrix::colSums(A)
  E_k <- tapply(eta, factor(z, levels = 1:K), sum); E_k[is.na(E_k)] <- 0
  for (i in seq_len(n)) {
    k <- z[i]; E_k[k] <- E_k[k] - eta[i]
    rate_i <- hyper$b_eta + sum(kappa[k, ] * E_k)
    eta[i] <- rgamma(1, shape = hyper$a_eta + G_i[i], rate = max(rate_i, 1e-10))
    E_k[k] <- E_k[k] + eta[i]
  }
  n_k <- as.integer(tabulate(z, nbins = K))
  for (k in seq_len(K)) {
    idx_k <- which(z == k); if (!length(idx_k)) next
    s_k <- sum(eta[idx_k])
    if (s_k > 0) eta[idx_k] <- n_k[k] * eta[idx_k] / s_k
    else         eta[idx_k] <- n_k[k] / length(idx_k)
  }
  eta
}

update_psi_sst_block <- function(z, psi_vec, A, i_idx, j_idx, N_edge, hyper) {
  K <- max(z); D <- K - 1L
  if (D < 1L) return(numeric(0))
  agg <- aggregate_by_distance(K, z_i = z[i_idx], z_j = z[j_idx],
                               A_ij = A[cbind(i_idx, j_idx)], N_edge = N_edge)
  B_d <- distance_totals(K, z_i = z[i_idx], z_j = z[j_idx], N_edge = N_edge)
  bar_omega <- draw_omega_bar(B_d = B_d, psi = psi_vec)
  update_psi_sst(K = K, bar_y = agg$bar_y, bar_omega = bar_omega,
                 psi_curr = psi_vec,
                 mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2,
                 tau2_0 = hyper$tau0^2, n_inner_sweeps = 4L)
}

# =============================================================================
# MINI-SAMPLER
# =============================================================================

run_chain <- function(A, z_init, kappa_init, psi_init, eta_init,
                      i_idx, j_idx, N_edge, edge_by_node,
                      hyper, psi_hyper,
                      method = c("collapsed", "augmented", "hybrid"),
                      n_iter = 500, burn = 150, seed = 42,
                      truth_z = NULL) {
  method <- match.arg(method)
  set.seed(seed)

  z <- z_init; kappa <- kappa_init; psi <- psi_init; eta <- eta_init

  z_sweep_fn <- switch(method,
    collapsed = sweep_z_sst_collapsed,
    augmented = sweep_z_sst_augmented,
    hybrid    = sweep_z_sst_hybrid
  )

  K_trace   <- integer(n_iter)
  ari_trace <- numeric(n_iter)
  z_time    <- numeric(n_iter)

  for (it in seq_len(n_iter)) {
    t1 <- proc.time()[3]
    res <- z_sweep_fn(z, A, eta, kappa, psi,
                      i_idx, j_idx, N_edge, edge_by_node,
                      hyper, psi_hyper)
    z_time[it] <- proc.time()[3] - t1
    z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K

    psi <- update_psi_sst_block(z, psi, A, i_idx, j_idx, N_edge, hyper)
    kappa <- update_kappa_block(z, eta, kappa, i_idx, j_idx, N_edge, hyper)
    eta   <- update_eta_block(z, eta, kappa, A, hyper)

    K_trace[it] <- K
    ari_trace[it] <- if (!is.null(truth_z)) mcclust::arandi(z, truth_z) else NA

    if (it %% 100 == 0)
      cat(sprintf("  [SST-%s] it=%d/%d K=%d z_sweep=%.3fs\n",
                  method, it, n_iter, K, z_time[it]))
  }

  post_idx <- (burn + 1):n_iter
  list(
    K_trace     = K_trace,
    ari_trace   = ari_trace,
    z_time      = z_time,
    K_mean      = mean(K_trace[post_idx]),
    ari_mean    = mean(ari_trace[post_idx], na.rm = TRUE),
    z_time_mean = mean(z_time[post_idx]),
    z_time_sd   = sd(z_time[post_idx]),
    z_final     = z
  )
}

# =============================================================================
# DATA GENERATION
# =============================================================================

logistic <- function(x) 1 / (1 + exp(-x))

generate_sst_network <- function(n, K, psi_true, kappa_true, z_true, eta_true) {
  A <- Matrix(0, n, n, sparse = TRUE)
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    lam <- eta_true[i] * eta_true[j] * kappa_true[z_true[i], z_true[j]]
    N <- rpois(1, lam)
    if (N > 0) {
      d <- abs(z_true[i] - z_true[j])
      th <- if (d == 0) 0 else psi_true[d]
      pr <- if (z_true[i] <= z_true[j]) logistic(th) else logistic(-th)
      fwd <- rbinom(1, N, pr)
      A[i, j] <- fwd; A[j, i] <- N - fwd
    }
  }
  A
}

build_edges <- function(A, n) {
  A_sym <- A + t(A)
  idx   <- which(A_sym > 0 & upper.tri(A_sym), arr.ind = TRUE)
  i_idx <- idx[, 1]; j_idx <- idx[, 2]
  N_edge <- A_sym[idx]
  edge_by_node <- vector("list", n)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], e)
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], e)
  }
  list(i_idx = i_idx, j_idx = j_idx, N_edge = N_edge, edge_by_node = edge_by_node)
}

# =============================================================================
# RUN BENCHMARK
# =============================================================================

set.seed(2026)

n <- 100; K_true <- 4
z_true <- rep(1:K_true, each = n / K_true)
eta_true <- rgamma(n, 10, 10)
n_k <- tabulate(z_true, nbins = K_true)
for (k in seq_len(K_true)) {
  idx_k <- which(z_true == k)
  eta_true[idx_k] <- n_k[k] * eta_true[idx_k] / sum(eta_true[idx_k])
}

psi_true_sst <- c(0.8, 1.5, 2.2)
kappa_true <- matrix(2, K_true, K_true); diag(kappa_true) <- 4

cat("\n--- Generating SST network (n=", n, ", K=", K_true, ") ---\n")
A_sst <- generate_sst_network(n, K_true, psi_true_sst, kappa_true, z_true, eta_true)
A_sst <- as.matrix(A_sst)
cat("  sum(A)=", sum(A_sst), " density=", round(sum(A_sst > 0) / (n * (n - 1)), 3), "\n")

edges <- build_edges(A_sst, n)

hyper <- list(a_kappa = 1, b_kappa = 1, a_eta = 2, b_eta = 2,
              mu0 = 0.5, sigma0 = 1.5, tau0 = 0.5)
psi_hyper <- list(gamma_gn = 0.5, mu0 = hyper$mu0,
                  sigma0 = hyper$sigma0, tau0 = hyper$tau0)

N_ITER <- 500
BURN   <- 150

cat("\n====== SST BENCHMARK (3 variants) ======\n")
cat("Running", N_ITER, "iterations (burn=", BURN, "), COLD start from K=1\n")

z_init <- rep(1L, n)
kappa_init_1 <- matrix(rgamma(1, 2, 1), 1, 1)
psi_init_1   <- numeric(0)
eta_init     <- rep(1, n)

methods <- c("collapsed", "augmented", "hybrid")
results <- list()

for (m in methods) {
  results[[m]] <- run_chain(
    A_sst, z_init, kappa_init_1, psi_init_1, eta_init,
    edges$i_idx, edges$j_idx, edges$N_edge, edges$edge_by_node,
    hyper, psi_hyper, method = m,
    n_iter = N_ITER, burn = BURN, seed = 42, truth_z = z_true
  )
}

# =============================================================================
# RESULTS
# =============================================================================

ess_bm <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 10) return(NA)
  b <- floor(sqrt(n)); a <- floor(n / b)
  bm <- sapply(seq_len(a), function(i) mean(x[((i-1)*b+1):(i*b)]))
  n / (1 + 2 * sum(abs(acf(bm, plot = FALSE)$acf[-1])))
}

cat("\n=============================================================================\n")
cat("                SST BENCHMARK RESULTS (3-way)\n")
cat("=============================================================================\n")

post_idx <- (BURN + 1):N_ITER

tab <- data.frame(
  method    = character(0),
  K_mean    = numeric(0),
  ARI       = numeric(0),
  z_time_ms = numeric(0),
  z_time_sd = numeric(0),
  ESS_K     = numeric(0),
  ESS_per_s = numeric(0),
  stringsAsFactors = FALSE
)

for (m in methods) {
  r <- results[[m]]
  ess <- ess_bm(r$K_trace[post_idx])
  tab <- rbind(tab, data.frame(
    method    = m,
    K_mean    = round(r$K_mean, 2),
    ARI       = round(r$ari_mean, 3),
    z_time_ms = round(r$z_time_mean * 1000, 1),
    z_time_sd = round(r$z_time_sd * 1000, 1),
    ESS_K     = round(ess, 1),
    ESS_per_s = round(ess / (sum(r$z_time[post_idx])), 3),
    stringsAsFactors = FALSE
  ))
}

print(tab, row.names = FALSE)

cat("\n--- z-sweep timing ---\n")
for (m in methods) {
  r <- results[[m]]
  cat(sprintf("  %-12s  %6.1f +/- %5.1f ms/sweep\n",
              m, r$z_time_mean * 1000, r$z_time_sd * 1000))
}

# Speedup relative to collapsed
t_col <- results[["collapsed"]]$z_time_mean
cat(sprintf("\n  Speedup augmented vs collapsed: %.2fx\n",
            t_col / results[["augmented"]]$z_time_mean))
cat(sprintf("  Speedup hybrid    vs collapsed: %.2fx\n",
            t_col / results[["hybrid"]]$z_time_mean))

cat("\n--- K trace (post-burn) ---\n")
for (m in methods) {
  r <- results[[m]]
  ess <- ess_bm(r$K_trace[post_idx])
  cat(sprintf("  SST %-12s: K_mean=%.2f  ESS(K)=%d  ARI=%.3f\n",
              m, r$K_mean, round(ess), r$ari_mean))
}

# --- Plots ---
pdf(file.path(out_dir, "benchmark_K_trace.pdf"), width = 10, height = 4)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
cols <- c("steelblue", "tomato", "seagreen")
for (idx in seq_along(methods)) {
  m <- methods[idx]; r <- results[[m]]
  plot(r$K_trace, type = "l", col = cols[idx],
       main = paste("SST", m), xlab = "iteration", ylab = "K",
       ylim = c(1, max(sapply(results, function(x) max(x$K_trace)))))
  abline(h = K_true, lty = 2, col = "grey40")
  abline(v = BURN, lty = 3, col = "grey60")
}
dev.off()

pdf(file.path(out_dir, "benchmark_timing.pdf"), width = 8, height = 4)
par(mfrow = c(1, 1), mar = c(5, 4, 3, 1))
bp <- barplot(
  sapply(methods, function(m) results[[m]]$z_time_mean * 1000),
  names.arg = methods, col = cols,
  main = "SST z-sweep time (ms/sweep)", ylab = "ms"
)
dev.off()

# Save RDS
saveRDS(results, file.path(out_dir, "benchmark_results.rds"))

cat("\nPlots saved to:", out_dir, "\n")
cat("Total script time:", round(proc.time()[3] - t0_script, 1), "s\n")
cat("Finished:", format(Sys.time()), "\n")

sink(type = "message"); sink()
close(log_con)
