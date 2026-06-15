#!/usr/bin/env Rscript
# =============================================================================
# Benchmark: Collapsed vs PG-Augmented z-update for SST (with mixing moves)
# =============================================================================
# Compares 4 SST variants:
#   (1) collapsed          – exact Binom-logistic kernel + numerical integration
#   (2) augmented          – quadratic PG kernel + closed-form new-block integral
#   (3) collapsed + mix    – (1) + adjacent_block_swap + split_merge
#   (4) augmented + mix    – (2) + adjacent_block_swap + split_merge
#
# Saves results to output/simulation/benchmark_sst_augmented/
# =============================================================================

# ---- Set up output directory and log file ----
out_dir <- "output/simulation/benchmark_sst_augmented"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(out_dir, paste0("benchmark_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat("=== Benchmark: SST collapsed vs PG-augmented z-update ===\n")
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

# --- SST z-sweep: collapsed ---
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
    if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
    }
  }
  K_out <- nrow(kappa)
  list(z = z, kappa = kappa, psi = psi, K = K_out)
}

# --- SST z-sweep: augmented ---
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
    if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
    }
  }
  K_out <- nrow(kappa)
  list(z = z, kappa = kappa, psi = psi, K = K_out)
}

# --- SST z-sweep + mixing moves: collapsed ---
sweep_z_sst_collapsed_mix <- function(z, A, eta, kappa, psi,
                                      i_idx, j_idx, N_edge, edge_by_node,
                                      hyper, psi_hyper) {
  res <- sweep_z_sst_collapsed(z, A, eta, kappa, psi,
                               i_idx, j_idx, N_edge, edge_by_node,
                               hyper, psi_hyper)
  z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K

  if (K >= 2L) {
    swap_res <- adjacent_block_swap_move(
      z = z, kappa = kappa, psi = psi, eta = eta,
      A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node, psi_mode = "distance"
    )
    z <- swap_res$z; kappa <- swap_res$kappa; psi <- swap_res$psi
  }

  sm_res <- split_merge_move(
    z = z, kappa = kappa, psi = psi, eta = eta,
    A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
    edge_by_node = edge_by_node,
    a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
    gamma_gn = psi_hyper$gamma_gn, psi_mode = "distance",
    hyper_psi = list(mu0 = hyper$mu0, sigma0 = hyper$sigma0, tau0 = hyper$tau0),
    n_restricted_scans = 3L
  )
  z <- sm_res$z; kappa <- sm_res$kappa; psi <- sm_res$psi
  K <- nrow(kappa)

  list(z = z, kappa = kappa, psi = psi, K = K)
}

# --- SST z-sweep + mixing moves: augmented ---
sweep_z_sst_augmented_mix <- function(z, A, eta, kappa, psi,
                                      i_idx, j_idx, N_edge, edge_by_node,
                                      hyper, psi_hyper) {
  res <- sweep_z_sst_augmented(z, A, eta, kappa, psi,
                               i_idx, j_idx, N_edge, edge_by_node,
                               hyper, psi_hyper)
  z <- res$z; kappa <- res$kappa; psi <- res$psi; K <- res$K

  if (K >= 2L) {
    swap_res <- adjacent_block_swap_move(
      z = z, kappa = kappa, psi = psi, eta = eta,
      A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node, psi_mode = "distance"
    )
    z <- swap_res$z; kappa <- swap_res$kappa; psi <- swap_res$psi
  }

  sm_res <- split_merge_move(
    z = z, kappa = kappa, psi = psi, eta = eta,
    A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
    edge_by_node = edge_by_node,
    a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
    gamma_gn = psi_hyper$gamma_gn, psi_mode = "distance",
    hyper_psi = list(mu0 = hyper$mu0, sigma0 = hyper$sigma0, tau0 = hyper$tau0),
    n_restricted_scans = 3L
  )
  z <- sm_res$z; kappa <- sm_res$kappa; psi <- sm_res$psi
  K <- nrow(kappa)

  list(z = z, kappa = kappa, psi = psi, K = K)
}

# =============================================================================
# HELPER BLOCKS (kappa, eta, psi updates; same as WST benchmark)
# =============================================================================

update_kappa_block <- function(z, eta, kappa, i_idx, j_idx, N_edge, hyper) {
  K <- nrow(kappa)
  E_k <- tapply(eta, factor(z, levels = 1:K), sum); E_k[is.na(E_k)] <- 0
  zi <- z[i_idx]; zj <- z[j_idx]
  R <- matrix(0, K, K); T_mat <- matrix(0, K, K)
  for (e in seq_along(N_edge)) {
    p <- min(zi[e], zj[e]); q <- max(zi[e], zj[e])
    R[p, q] <- R[p, q] + N_edge[e]
    T_mat[p, q] <- T_mat[p, q] + eta[i_idx[e]] * eta[j_idx[e]]
  }
  for (k in seq_len(K)) for (l in k:K) {
    kappa[k, l] <- kappa[l, k] <- rgamma(1,
      shape = hyper$a_kappa + R[k, l],
      rate  = hyper$b_kappa + T_mat[k, l])
  }
  kappa
}

update_eta_block <- function(z, eta, kappa, A, hyper) {
  n <- length(z); K <- nrow(kappa)
  G_i <- Matrix::rowSums(A) + Matrix::colSums(A)
  E_k <- tapply(eta, factor(z, levels = 1:K), sum); E_k[is.na(E_k)] <- 0

  for (i in seq_len(n)) {
    k <- z[i]; E_k[k] <- E_k[k] - eta[i]
    rate_i <- hyper$b_eta + sum(kappa[k, ] * E_k)
    eta[i] <- rgamma(1, shape = hyper$a_eta + G_i[i],
                     rate = max(rate_i, 1e-10))
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
  K <- max(z)
  D <- K - 1L
  if (D < 1L) return(numeric(0))
  agg <- aggregate_by_distance(K, z_i = z[i_idx], z_j = z[j_idx],
                               A_ij = A[cbind(i_idx, j_idx)], N_edge = N_edge)
  B_d <- distance_totals(K, z_i = z[i_idx], z_j = z[j_idx], N_edge = N_edge)
  bar_omega <- draw_omega_bar(B_d = B_d, psi = psi_vec)
  psi_new <- update_psi_sst(K = K, bar_y = agg$bar_y, bar_omega = bar_omega,
                            psi_curr = psi_vec,
                            mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2,
                            tau2_0 = hyper$tau0^2, n_inner_sweeps = 4L)
  psi_new
}

# =============================================================================
# MINI-SAMPLER
# =============================================================================

run_chain <- function(A, z_init, kappa_init, psi_init, eta_init,
                      i_idx, j_idx, N_edge, edge_by_node,
                      hyper, psi_hyper,
                      method = c("collapsed", "augmented",
                                 "collapsed_mix", "augmented_mix"),
                      n_iter = 500, burn = 100, seed = 42,
                      truth_z = NULL) {
  method <- match.arg(method)
  set.seed(seed)

  z     <- z_init
  kappa <- kappa_init
  psi   <- psi_init
  eta   <- eta_init
  n     <- length(z)

  z_sweep_fn <- switch(method,
    collapsed      = sweep_z_sst_collapsed,
    augmented      = sweep_z_sst_augmented,
    collapsed_mix  = sweep_z_sst_collapsed_mix,
    augmented_mix  = sweep_z_sst_augmented_mix
  )

  K_trace   <- integer(n_iter)
  ari_trace <- numeric(n_iter)
  z_time    <- numeric(n_iter)
  psi_time  <- numeric(n_iter)

  for (it in seq_len(n_iter)) {
    t1 <- proc.time()[3]
    res <- z_sweep_fn(z, A, eta, kappa, psi,
                      i_idx, j_idx, N_edge, edge_by_node,
                      hyper, psi_hyper)
    z_time[it] <- proc.time()[3] - t1

    z     <- res$z
    kappa <- res$kappa
    psi   <- res$psi
    K     <- res$K

    t2 <- proc.time()[3]
    psi <- update_psi_sst_block(z, psi, A, i_idx, j_idx, N_edge, hyper)
    psi_time[it] <- proc.time()[3] - t2

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
    K_trace   = K_trace,
    ari_trace = ari_trace,
    z_time    = z_time,
    psi_time  = psi_time,
    K_mean    = mean(K_trace[post_idx]),
    ari_mean  = mean(ari_trace[post_idx], na.rm = TRUE),
    z_time_mean = mean(z_time[post_idx]),
    z_time_sd   = sd(z_time[post_idx]),
    z_final   = z,
    kappa_final = kappa,
    psi_final  = psi,
    eta_final  = eta
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
      A[i, j] <- fwd
      A[j, i] <- N - fwd
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
  list(i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
       edge_by_node = edge_by_node)
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
kappa_true <- matrix(2, K_true, K_true)
diag(kappa_true) <- 4

cat("\n--- Generating SST network (n=", n, ", K=", K_true, ") ---\n")
A_sst <- generate_sst_network(n, K_true, psi_true_sst, kappa_true, z_true, eta_true)
A_sst <- as.matrix(A_sst)
cat("  sum(A)=", sum(A_sst), " density=", round(sum(A_sst > 0) / (n * (n - 1)), 3), "\n")

edges <- build_edges(A_sst, n)

hyper <- list(
  a_kappa = 1, b_kappa = 1,
  a_eta = 2, b_eta = 2,
  mu0 = 0.5, sigma0 = 1.5,
  tau0 = 0.5
)
psi_hyper <- list(
  gamma_gn = 0.5,
  mu0 = hyper$mu0,
  sigma0 = hyper$sigma0,
  tau0 = hyper$tau0
)

N_ITER <- 800
BURN   <- 200

cat("\n====== SST BENCHMARK (4 variants) ======\n")
cat("Running", N_ITER, "iterations (burn=", BURN, "), COLD start from K=1\n")

z_init <- rep(1L, n)
kappa_init_1 <- matrix(rgamma(1, 2, 1), 1, 1)
psi_init_sst_1 <- numeric(0)
eta_init <- rep(1, n)

res_collapsed <- run_chain(
  A_sst, z_init, kappa_init_1, psi_init_sst_1, eta_init,
  edges$i_idx, edges$j_idx, edges$N_edge, edges$edge_by_node,
  hyper, psi_hyper, method = "collapsed",
  n_iter = N_ITER, burn = BURN, seed = 42, truth_z = z_true
)

res_augmented <- run_chain(
  A_sst, z_init, kappa_init_1, psi_init_sst_1, eta_init,
  edges$i_idx, edges$j_idx, edges$N_edge, edges$edge_by_node,
  hyper, psi_hyper, method = "augmented",
  n_iter = N_ITER, burn = BURN, seed = 42, truth_z = z_true
)

res_collapsed_mix <- run_chain(
  A_sst, z_init, kappa_init_1, psi_init_sst_1, eta_init,
  edges$i_idx, edges$j_idx, edges$N_edge, edges$edge_by_node,
  hyper, psi_hyper, method = "collapsed_mix",
  n_iter = N_ITER, burn = BURN, seed = 42, truth_z = z_true
)

res_augmented_mix <- run_chain(
  A_sst, z_init, kappa_init_1, psi_init_sst_1, eta_init,
  edges$i_idx, edges$j_idx, edges$N_edge, edges$edge_by_node,
  hyper, psi_hyper, method = "augmented_mix",
  n_iter = N_ITER, burn = BURN, seed = 42, truth_z = z_true
)

# =============================================================================
# RESULTS
# =============================================================================

ess_bm <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 10) return(NA)
  b <- floor(sqrt(n))
  a <- floor(n / b)
  batch_means <- sapply(seq_len(a), function(i) mean(x[((i-1)*b+1):(i*b)]))
  var_grand <- var(x)
  var_batch <- var(batch_means)
  if (var_batch <= 0 || var_grand <= 0) return(n)
  pmin(n, n * var_grand / (b * var_batch))
}

cat("\n")
cat("=", rep("=", 76), "\n", sep = "")
cat("                SST BENCHMARK RESULTS\n")
cat("=", rep("=", 76), "\n", sep = "")

results <- data.frame(
  method  = c("collapsed", "augmented", "collapsed+mix", "augmented+mix"),
  K_mean  = c(res_collapsed$K_mean, res_augmented$K_mean,
              res_collapsed_mix$K_mean, res_augmented_mix$K_mean),
  ARI     = c(res_collapsed$ari_mean, res_augmented$ari_mean,
              res_collapsed_mix$ari_mean, res_augmented_mix$ari_mean),
  z_time_ms = c(res_collapsed$z_time_mean, res_augmented$z_time_mean,
                res_collapsed_mix$z_time_mean, res_augmented_mix$z_time_mean) * 1000,
  z_time_sd = c(res_collapsed$z_time_sd, res_augmented$z_time_sd,
                res_collapsed_mix$z_time_sd, res_augmented_mix$z_time_sd) * 1000,
  ESS_K   = c(ess_bm(res_collapsed$K_trace[(BURN+1):N_ITER]),
              ess_bm(res_augmented$K_trace[(BURN+1):N_ITER]),
              ess_bm(res_collapsed_mix$K_trace[(BURN+1):N_ITER]),
              ess_bm(res_augmented_mix$K_trace[(BURN+1):N_ITER])),
  stringsAsFactors = FALSE
)
results$speedup <- NA
results$speedup[2] <- results$z_time_ms[1] / results$z_time_ms[2]
results$speedup[4] <- results$z_time_ms[3] / results$z_time_ms[4]
results$ESS_per_s <- results$ESS_K / (results$z_time_ms / 1000 * (N_ITER - BURN))

print(results, digits = 3)

cat("\n--- SST z-sweep timing ---\n")
cat(sprintf("  Collapsed:      %.1f +/- %.1f ms/sweep\n",
            res_collapsed$z_time_mean * 1000,
            res_collapsed$z_time_sd * 1000))
cat(sprintf("  Augmented:      %.1f +/- %.1f ms/sweep\n",
            res_augmented$z_time_mean * 1000,
            res_augmented$z_time_sd * 1000))
cat(sprintf("  Collapsed+mix:  %.1f +/- %.1f ms/sweep\n",
            res_collapsed_mix$z_time_mean * 1000,
            res_collapsed_mix$z_time_sd * 1000))
cat(sprintf("  Augmented+mix:  %.1f +/- %.1f ms/sweep\n",
            res_augmented_mix$z_time_mean * 1000,
            res_augmented_mix$z_time_sd * 1000))
cat(sprintf("  Speedup (no mix):   %.2fx\n",
            res_collapsed$z_time_mean / res_augmented$z_time_mean))
cat(sprintf("  Speedup (with mix): %.2fx\n",
            res_collapsed_mix$z_time_mean / res_augmented_mix$z_time_mean))

cat("\n--- K trace (post-burn) ---\n")
all_res <- list(res_collapsed, res_augmented,
                res_collapsed_mix, res_augmented_mix)
all_labels <- c("collapsed", "augmented", "collapsed+mix", "augmented+mix")
for (j in seq_along(all_res)) {
  r <- all_res[[j]]
  cat(sprintf("  SST %-16s: K_mean=%.2f  ESS(K)=%.0f  ARI=%.3f\n",
              all_labels[j], r$K_mean,
              ess_bm(r$K_trace[(BURN+1):N_ITER]),
              r$ari_mean))
}

# ---- Save diagnostic plots ----
pdf(file.path(out_dir, "K_traces.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(res_collapsed$K_trace, type = "l", col = "steelblue",
     main = "SST collapsed: K trace", xlab = "iteration", ylab = "K")
abline(h = K_true, lty = 2, col = "red")

plot(res_augmented$K_trace, type = "l", col = "darkorange",
     main = "SST augmented: K trace", xlab = "iteration", ylab = "K")
abline(h = K_true, lty = 2, col = "red")

plot(res_collapsed_mix$K_trace, type = "l", col = "steelblue",
     main = "SST collapsed+mix: K trace", xlab = "iteration", ylab = "K")
abline(h = K_true, lty = 2, col = "red")

plot(res_augmented_mix$K_trace, type = "l", col = "darkorange",
     main = "SST augmented+mix: K trace", xlab = "iteration", ylab = "K")
abline(h = K_true, lty = 2, col = "red")

dev.off()

pdf(file.path(out_dir, "z_sweep_times.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

boxplot(list(collapsed = res_collapsed$z_time * 1000,
             augmented = res_augmented$z_time * 1000),
        main = "SST (no mix): z-sweep time (ms)", ylab = "ms",
        col = c("steelblue", "darkorange"))

boxplot(list("collapsed+mix" = res_collapsed_mix$z_time * 1000,
             "augmented+mix" = res_augmented_mix$z_time * 1000),
        main = "SST (with mix): z-sweep time (ms)", ylab = "ms",
        col = c("steelblue", "darkorange"))

dev.off()

# ---- Save results RDS ----
saveRDS(results, file.path(out_dir, "benchmark_results.rds"))

cat("\nPlots saved to:", out_dir, "\n")
cat("Total script time:", round(proc.time()[3] - t0_script, 1), "s\n")
cat("Finished:", format(Sys.time()), "\n")

# Close log
sink(type = "message")
sink(type = "output")
close(log_con)

cat("SST benchmark complete. Log saved to:", log_file, "\n")
