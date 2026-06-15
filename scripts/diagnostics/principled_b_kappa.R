#!/usr/bin/env Rscript
# ============================================================
# principled_b_kappa.R
#
# Derive and verify: a_kappa = 1,  b_kappa = c · n / K_expected
#
# Theoretical prediction:
#   birth_advantage_kappa ≈ −R_i · ln(c)
# where R_i = total dyad count for node i, c = b·K/n.
#
# Three candidate rules are compared:
#   Rule A: a = 1           (constant, our proposal)
#   Rule B: a = λ̄ · b       (mean-matched)
#   Rule C: a = λ̄           (density-calibrated constant)
# ============================================================

suppressPackageStartupMessages({
  library(Rcpp)
})

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")
sourceCpp("core/block_totals_for_poisson_cpp.cpp")

# ---- gp_marginal (vectorised) ----
gp_marginal <- function(R, T, a, b) {
  out <- rep(-Inf, length(R))
  ok0 <- (T == 0)
  out[ok0 & (R == 0)] <- 0
  ok <- !ok0
  if (any(ok))
    out[ok] <- lgamma(a + R[ok]) - lgamma(a) +
      a * log(b) - (a + R[ok]) * log(b + T[ok])
  out
}

# ---- Load datasets (mirrors choose_dataset() in application.R) ----
load_dataset <- function(name) {
  switch(name,
    mountain_goats = {
      files <- list.files("./data/ShizukaMcDonald_Data",
                          full.names = TRUE, pattern = "\\.csv$")
      n_each <- vapply(files, function(f) nrow(read.csv(f, row.names = 1)),
                        FUN.VALUE = integer(1))
      A <- as.matrix(read.csv(files[which.max(n_each)],
                               row.names = 1, check.names = FALSE))
      storage.mode(A) <- "integer"; A
    },
    citations_data = {
      A <- as.matrix(read.csv("./data/Citations_application/cross-citation-matrix.csv",
                               row.names = 1, header = TRUE, check.names = FALSE))
      diag(A) <- 0; storage.mode(A) <- "integer"; A
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
    },
    stop("Unknown dataset: ", name)
  )
}

# ---- Compute per-node sufficient statistics -----------------------
# Given adjacency A and partition z (labels 1..K), for each node i
# compute: r_add[ℓ], t_add[ℓ], R_minus[g,h], T_minus[g,h]
node_suffstats <- function(A, z, i) {
  n   <- nrow(A)
  K   <- max(z)
  N   <- A + t(A);  diag(N) <- 0       # dyad-total matrix

  old_k <- z[i]
  eta   <- rep(1, n)                     # simplify: all eta = 1

  # block sizes and sizes-minus-i
  n_b <- tabulate(z, nbins = K)
  n_b_minus <- n_b;  n_b_minus[old_k] <- n_b_minus[old_k] - 1L

  # --- r_add: total dyad count between i and each block (excluding self)
  r_add <- numeric(K)
  for (j in seq_len(n)) {
    if (j == i) next
    r_add[z[j]] <- r_add[z[j]] + N[i, j]
  }

  # --- t_add: exposure of i to each block (eta = 1 → t_add = n_b_minus)
  t_add <- as.numeric(n_b_minus)

  # --- R full and T full (upper-triangular block totals)
  bt <- block_totals_for_poisson_cpp(as.integer(z), eta,
    integer(0), integer(0), numeric(0), K)
  # build R and T from observed dyads
  R_full <- matrix(0, K, K);  T_full <- matrix(0, K, K)
  # Rkl from edge list
  upper_idx <- which(upper.tri(N, diag = FALSE), arr.ind = TRUE)
  for (e in seq_len(nrow(upper_idx))) {
    u <- upper_idx[e, 1]; v <- upper_idx[e, 2]
    g <- z[u]; h <- z[v]
    p <- min(g, h); q <- max(g, h)
    R_full[p, q] <- R_full[p, q] + N[u, v]
  }
  # Tkl (eta = 1)
  for (g in 1:K) for (h in g:K) {
    if (g == h) T_full[g, h] <- n_b[g] * (n_b[g] - 1) / 2
    else        T_full[g, h] <- n_b[g] * n_b[h]
  }
  # symmetrise
  R_full <- R_full + t(R_full) - diag(diag(R_full))
  T_full <- T_full + t(T_full) - diag(diag(T_full))

  # --- R_minus, T_minus: subtract i's contribution from its old block
  R_minus <- R_full;  T_minus <- T_full
  for (ell in 1:K) {
    p <- min(old_k, ell); q <- max(old_k, ell)
    R_minus[p, q] <- R_minus[p, q] - r_add[ell]
    R_minus[q, p] <- R_minus[p, q]
    T_minus[p, q] <- T_minus[p, q] - t_add[ell]
    T_minus[q, p] <- T_minus[p, q]
    R_minus[p, q] <- max(R_minus[p, q], 0)
    R_minus[q, p] <- R_minus[p, q]
    T_minus[p, q] <- max(T_minus[p, q], 0)
    T_minus[q, p] <- T_minus[p, q]
  }

  list(r_add = r_add, t_add = t_add,
       R_minus = R_minus, T_minus = T_minus,
       R_i = sum(r_add), old_k = old_k, K = K)
}

# ---- Compute birth advantage from kappa for a single node --------
kappa_birth_advantage <- function(ss, a, b) {
  K <- ss$K

  # Birth: full marginals for K new block pairs
  birth_kappa <- sum(gp_marginal(ss$r_add, ss$t_add, a, b))

  # Join: best existing block
  join_scores <- numeric(K)
  for (k in 1:K) {
    acc <- 0
    for (ell in 1:K) {
      if (ss$r_add[ell] == 0 && ss$t_add[ell] == 0) next
      p <- min(k, ell); q <- max(k, ell)
      acc <- acc + gp_marginal(ss$R_minus[p,q] + ss$r_add[ell],
                                ss$T_minus[p,q] + ss$t_add[ell], a, b) -
                   gp_marginal(ss$R_minus[p,q], ss$T_minus[p,q], a, b)
    }
    join_scores[k] <- acc
  }

  birth_kappa - max(join_scores)
}

# ============================================================
# MAIN ANALYSIS
# ============================================================
datasets <- c("mountain_goats", "citations_data", "macaques_data", "high_school")
K_exp    <- 5L
set.seed(42)

# b grid (in units of c = b·K/n, so b = c·n/K)
c_grid <- c(0.1, 0.2, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 50)

cat("\n========== PRINCIPLED b_kappa ANALYSIS ==========\n\n")
cat("Theoretical prediction: BA_kappa ≈ -R_i · ln(c),  c = b·K/n\n")
cat("  → crossover at c = 1 (b = n/K)\n")
cat("  → birth suppressed for c > 1\n\n")

all_results <- list()

for (ds in datasets) {
  A <- load_dataset(ds)
  n <- nrow(A)
  N_mat <- A + t(A); diag(N_mat) <- 0
  mean_N <- mean(N_mat[upper.tri(N_mat)])

  # Create a balanced random partition
  z <- rep(seq_len(K_exp), length.out = n)
  z <- sample(z)  # shuffle

  cat(sprintf("=== %s (n=%d, mean_N=%.2f) ===\n", ds, n, mean_N))

  # Precompute sufficient statistics for all nodes
  ss_all <- lapply(seq_len(n), function(i) node_suffstats(A, z, i))
  R_all  <- sapply(ss_all, function(s) s$R_i)

  cat(sprintf("  Node degrees (R_i): median=%.1f, range=[%.0f, %.0f]\n",
              median(R_all), min(R_all), max(R_all)))

  # For each c value, compute BA under three rules
  for (ci in seq_along(c_grid)) {
    cc <- c_grid[ci]
    b_val <- cc * n / K_exp

    # Rule A: a = 1
    ba_A <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = 1, b = b_val))
    # Rule B: a = mean_N * b_val  (mean-matched)
    a_B <- mean_N * b_val
    ba_B <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = a_B, b = b_val))
    # Rule C: a = mean_N  (density constant)
    ba_C <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = mean_N, b = b_val))

    # Theoretical prediction: BA ≈ -R_i * ln(c)
    ba_theory <- -R_all * log(cc)

    all_results[[length(all_results) + 1L]] <- data.frame(
      dataset  = ds,
      n        = n,
      mean_N   = mean_N,
      c        = cc,
      b        = b_val,
      rule     = rep(c("A_const", "B_mean_match", "C_dens_const", "theory"), each = n),
      BA       = c(ba_A, ba_B, ba_C, ba_theory),
      R_i      = rep(R_all, 4),
      stringsAsFactors = FALSE
    )
  }

  # ---- Summary table for Rule A ----
  cat("\n  Rule A (a=1): median BA_kappa across nodes\n")
  cat(sprintf("  %6s  %8s  %8s  %8s  %10s\n",
              "c", "b", "med_BA", "theory", "frac_birth"))
  for (cc in c_grid) {
    b_val <- cc * n / K_exp
    ba <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = 1, b = b_val))
    med_ba <- median(ba)
    med_theory <- -median(R_all) * log(cc)
    frac <- mean(ba > 0)
    cat(sprintf("  %6.1f  %8.1f  %8.2f  %8.2f  %10.1f%%\n",
                cc, b_val, med_ba, med_theory, 100 * frac))
  }

  cat("\n  Comparison of rules at c=5:\n")
  b5 <- 5 * n / K_exp
  ba_A5 <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = 1, b = b5))
  ba_B5 <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = mean_N * b5, b = b5))
  ba_C5 <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = mean_N, b = b5))
  cat(sprintf("    Rule A (a=1):         median BA = %.2f, frac_birth = %.1f%%\n",
              median(ba_A5), 100 * mean(ba_A5 > 0)))
  cat(sprintf("    Rule B (a=λ̄·b):      median BA = %.2f, frac_birth = %.1f%%\n",
              median(ba_B5), 100 * mean(ba_B5 > 0)))
  cat(sprintf("    Rule C (a=λ̄):        median BA = %.2f, frac_birth = %.1f%%\n",
              median(ba_C5), 100 * mean(ba_C5 > 0)))
  cat(sprintf("    Theory (-R̃·ln5):     median BA = %.2f\n\n",
              -median(R_all) * log(5)))
}

# ---- Cross-dataset summary ----
cat("\n========== CROSS-DATASET SUMMARY ==========\n\n")
cat("The crossover c* (where median BA_kappa ≈ 0) should be near c=1 for all datasets.\n")
cat("Proposed rule: a_kappa = 1, b_kappa = c · n / K_expected,  default c = 5\n\n")

for (ds in datasets) {
  A <- load_dataset(ds)
  n <- nrow(A)
  N_mat <- A + t(A); diag(N_mat) <- 0
  mean_N <- mean(N_mat[upper.tri(N_mat)])
  z <- rep(seq_len(K_exp), length.out = n)
  set.seed(42); z <- sample(z)
  ss_all <- lapply(seq_len(n), function(i) node_suffstats(A, z, i))

  # Find the crossover c where median BA crosses 0
  f_cross <- function(logc) {
    cc <- exp(logc)
    b_val <- cc * n / K_exp
    ba <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = 1, b = b_val))
    median(ba)
  }
  # Quick scan
  logc_grid <- log(c(0.1, 0.5, 1, 2, 5, 10, 50))
  vals <- sapply(logc_grid, f_cross)

  # Find approximate crossover
  sign_changes <- which(diff(sign(vals)) != 0)
  if (length(sign_changes) > 0) {
    lo <- logc_grid[sign_changes[1]]
    hi <- logc_grid[sign_changes[1] + 1]
    root <- tryCatch(
      uniroot(f_cross, interval = c(lo, hi))$root,
      error = function(e) NA
    )
    c_star <- if (!is.na(root)) exp(root) else NA
  } else {
    c_star <- NA
  }

  b_default <- 5 * n / K_exp
  ba_default <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = 1, b = b_default))

  cat(sprintf("  %15s: n=%d, λ̄=%.2f, c*=%.2f, b*=%.1f | default: b=%.0f, med_BA=%.1f, birth%%=%.0f%%\n",
              ds, n, mean_N,
              ifelse(is.na(c_star), NA, c_star),
              ifelse(is.na(c_star), NA, c_star * n / K_exp),
              b_default, median(ba_default), 100 * mean(ba_default > 0)))
}

cat("\n========== MONOTONICITY CHECK ==========\n")
cat("Verify that median BA is monotonically decreasing in c for Rule A:\n\n")

for (ds in datasets) {
  A <- load_dataset(ds)
  n <- nrow(A)
  z <- rep(seq_len(K_exp), length.out = n)
  set.seed(42); z <- sample(z)
  ss_all <- lapply(seq_len(n), function(i) node_suffstats(A, z, i))

  prev <- Inf; monotone <- TRUE
  cat(sprintf("  %s: ", ds))
  for (cc in c_grid) {
    b_val <- cc * n / K_exp
    ba <- sapply(ss_all, function(s) kappa_birth_advantage(s, a = 1, b = b_val))
    med <- median(ba)
    if (med > prev + 0.01) { monotone <- FALSE }
    prev <- med
  }
  cat(if (monotone) "MONOTONE ✓" else "NOT MONOTONE ✗", "\n")
}

cat("\n========== DONE ==========\n")
