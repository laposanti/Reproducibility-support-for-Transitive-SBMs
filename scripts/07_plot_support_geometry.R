#!/usr/bin/env Rscript

# =============================================================================
# Simulate support coverage in a 3D cube (K = 3)
# Supports compared:
#   - WST
#   - SST (generic strong-transitivity-compatible construction)
#   - Toeplitz-SST (distance-based)
#   - Bradley-Terry (LST)
#
# The script:
#   1) Samples psi = (psi12, psi13, psi23) under each support.
#   2) Varies hyperparameters across random settings.
#   3) Voxelises [0, cube_max]^3 and computes average support coverage.
#   4) Produces an interactive 2x2 3D plot of average support coverage.
#   5) Writes summary metrics of relative restrictiveness.
#
# Usage:
#   Rscript scripts/07_plot_support_geometry.R
#
# Optional env vars:
#   OUT_DIR           (default: output/diagnostics/support_geometry)
#   SEED              (default: 123)
#   K_SIM             (default: 3)
#   K_GRID            (default: 3,4,5,6)
#   CUBE_MAX          (default: 4)
#   N_BINS            (default: 26)
#   N_PARAM_SETTINGS  (default: 36)
#   N_SAMPLES         (default: 5000)
#   OCC_THRESHOLD     (default: 0.10)
# =============================================================================

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
if (length(file_arg)) {
  script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg[1L])),
                               winslash = "/", mustWork = TRUE)
  setwd(normalizePath(file.path(dirname(script_path), ".."),
                      winslash = "/", mustWork = TRUE))
}

suppressPackageStartupMessages({
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Install with install.packages('plotly').")
  }
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    stop("Package 'htmlwidgets' is required. Install with install.packages('htmlwidgets').")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with install.packages('ggplot2').")
  }
})

# -----------------------------
# Configuration
# -----------------------------
out_dir <- Sys.getenv("OUT_DIR", unset = "output/diagnostics/support_geometry")
seed <- as.integer(Sys.getenv("SEED", unset = "123"))
K_sim <- as.integer(Sys.getenv("K_SIM", unset = "3"))
k_grid_raw <- Sys.getenv("K_GRID", unset = "3,4,5,6")
cube_max <- as.numeric(Sys.getenv("CUBE_MAX", unset = "4"))
n_bins <- as.integer(Sys.getenv("N_BINS", unset = "26"))
n_param_settings <- as.integer(Sys.getenv("N_PARAM_SETTINGS", unset = "36"))
n_samples <- as.integer(Sys.getenv("N_SAMPLES", unset = "5000"))
occ_threshold <- as.numeric(Sys.getenv("OCC_THRESHOLD", unset = "0.10"))

if (is.na(seed)) seed <- 123L
if (is.na(K_sim) || K_sim < 3) stop("K_SIM must be an integer >= 3")
k_grid <- as.integer(trimws(strsplit(k_grid_raw, ",", fixed = TRUE)[[1]]))
if (length(k_grid) == 0 || any(is.na(k_grid)) || any(k_grid < 3)) {
  stop("K_GRID must be a comma-separated list of integers >= 3 (example: 3,4,5,6)")
}
k_grid <- sort(unique(k_grid))
if (is.na(cube_max) || cube_max <= 0) stop("CUBE_MAX must be > 0")
if (is.na(n_bins) || n_bins < 6) stop("N_BINS must be >= 6")
if (is.na(n_param_settings) || n_param_settings < 4) stop("N_PARAM_SETTINGS must be >= 4")
if (is.na(n_samples) || n_samples < 1000) stop("N_SAMPLES must be >= 1000")
if (is.na(occ_threshold) || occ_threshold <= 0 || occ_threshold >= 1) {
  stop("OCC_THRESHOLD must be in (0,1)")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

# -----------------------------
# Helpers
# -----------------------------

# Positive-truncated normal via rejection sampling.
rtnorm_pos <- function(n, mean, sd) {
  if (sd <= 0) return(rep(max(mean, 1e-8), n))
  x <- rnorm(n, mean = mean, sd = sd)
  bad <- which(x <= 0)
  while (length(bad) > 0) {
    x[bad] <- rnorm(length(bad), mean = mean, sd = sd)
    bad <- which(x <= 0)
  }
  x
}

clip_cube <- function(df, cube_max) {
  df[df$psi12 <= cube_max & df$psi13 <= cube_max & df$psi23 <= cube_max, , drop = FALSE]
}

lin_index_3d <- function(ix, iy, iz, n_bins) {
  ix + (iy - 1L) * n_bins + (iz - 1L) * n_bins * n_bins
}

voxel_occupancy <- function(points, breaks) {
  n_bins <- length(breaks) - 1L
  n_vox <- n_bins^3
  occ <- integer(n_vox)
  if (nrow(points) == 0) return(occ)

  ix <- findInterval(points$psi12, breaks, all.inside = TRUE)
  iy <- findInterval(points$psi13, breaks, all.inside = TRUE)
  iz <- findInterval(points$psi23, breaks, all.inside = TRUE)

  idx <- unique(lin_index_3d(ix, iy, iz, n_bins))
  occ[idx] <- 1L
  occ
}

random_hyper <- function(family) {
  switch(
    family,
    WST = list(
      mu = runif(1, 0.15, 1.50),
      sd = runif(1, 0.10, 1.00)
    ),
    SST = list(
      mu_base = runif(1, 0.15, 1.25),
      sd_base = runif(1, 0.10, 0.90),
      mu_gap = runif(1, 0.05, 1.00),
      sd_gap = runif(1, 0.05, 0.80)
    ),
    Toeplitz_SST = list(
      mu_phi1 = runif(1, 0.15, 1.25),
      sd_phi1 = runif(1, 0.08, 0.80),
      mu_delta = runif(1, 0.05, 1.00),
      sd_delta = runif(1, 0.05, 0.80)
    ),
    LST = list(
      mu_a = runif(1, 0.15, 1.25),
      sd_a = runif(1, 0.08, 0.80),
      mu_b = runif(1, 0.15, 1.25),
      sd_b = runif(1, 0.08, 0.80)
    ),
    stop("Unknown family: ", family)
  )
}

sample_support <- function(family, n, hyper, cube_max, K_sim) {
  # Coordinates are projected psi = (psi12, psi13, psi23), anchored to ranks (1,2,3).
  # For K_sim > 3, this remains a fixed low-dimensional projection.
  if (family == "WST") {
    psi12 <- rtnorm_pos(n, hyper$mu, hyper$sd)
    psi13 <- rtnorm_pos(n, hyper$mu, hyper$sd)
    psi23 <- rtnorm_pos(n, hyper$mu, hyper$sd)
    pts <- data.frame(psi12 = psi12, psi13 = psi13, psi23 = psi23)
    return(clip_cube(pts, cube_max))
  }

  if (family == "SST") {
    # Generic SST-compatible construction for K=3:
    # psi12 > 0, psi23 > 0, psi13 >= max(psi12, psi23).
    psi12 <- rtnorm_pos(n, hyper$mu_base, hyper$sd_base)
    psi23 <- rtnorm_pos(n, hyper$mu_base, hyper$sd_base)
    gap <- rtnorm_pos(n, hyper$mu_gap, hyper$sd_gap)
    psi13 <- pmax(psi12, psi23) + gap
    pts <- data.frame(psi12 = psi12, psi13 = psi13, psi23 = psi23)
    return(clip_cube(pts, cube_max))
  }

  if (family == "Toeplitz_SST") {
    # Toeplitz distance-based SST for K=3:
    # psi12 = phi1, psi23 = phi1, psi13 = phi2, with 0 < phi1 <= phi2.
    phi1 <- rtnorm_pos(n, hyper$mu_phi1, hyper$sd_phi1)
    delta <- rtnorm_pos(n, hyper$mu_delta, hyper$sd_delta)
    phi2 <- phi1 + delta
    pts <- data.frame(psi12 = phi1, psi13 = phi2, psi23 = phi1)
    return(clip_cube(pts, cube_max))
  }

  if (family == "LST") {
    # LST for K=3:
    # psi12 = a > 0, psi23 = b > 0, psi13 = a + b.
    a <- rtnorm_pos(n, hyper$mu_a, hyper$sd_a)
    b <- rtnorm_pos(n, hyper$mu_b, hyper$sd_b)
    pts <- data.frame(psi12 = a, psi13 = a + b, psi23 = b)
    return(clip_cube(pts, cube_max))
  }

  stop("Unknown family: ", family)
}

build_centers_df <- function(breaks, n_bins) {
  mids <- (breaks[-1] + breaks[-length(breaks)]) / 2
  grid <- expand.grid(ix = seq_len(n_bins), iy = seq_len(n_bins), iz = seq_len(n_bins))
  grid$psi12 <- mids[grid$ix]
  grid$psi13 <- mids[grid$iy]
  grid$psi23 <- mids[grid$iz]
  grid$lin <- lin_index_3d(grid$ix, grid$iy, grid$iz, n_bins)
  grid
}

build_pair_index <- function(K) {
  pairs <- utils::combn(K, 2)
  idx_mat <- matrix(NA_integer_, nrow = K, ncol = K)
  for (col in seq_len(ncol(pairs))) {
    i <- pairs[1, col]
    j <- pairs[2, col]
    idx_mat[i, j] <- col
    idx_mat[j, i] <- col
  }
  list(
    pair_i = pairs[1, ],
    pair_j = pairs[2, ],
    idx_mat = idx_mat
  )
}

build_triple_index <- function(K, idx_mat) {
  triples <- utils::combn(K, 3)
  n_t <- ncol(triples)
  triple_cols <- matrix(NA_integer_, nrow = n_t, ncol = 3)
  for (t in seq_len(n_t)) {
    i <- triples[1, t]
    l <- triples[2, t]
    m <- triples[3, t]
    triple_cols[t, 1] <- idx_mat[i, l]
    triple_cols[t, 2] <- idx_mat[i, m]
    triple_cols[t, 3] <- idx_mat[l, m]
  }
  triple_cols
}

sample_psi_matrix_allpairs <- function(family, n, K, hyper, pair_i, pair_j) {
  pK <- length(pair_i)

  if (family == "WST") {
    out <- matrix(rtnorm_pos(n * pK, hyper$mu, hyper$sd), nrow = n, ncol = pK)
    return(out)
  }

  if (family == "SST") {
    # Build from adjacent pairs upward using conditional gap sampling so that
    # SST holds by construction and values stay moderate inside the cube.
    # psi[k, k+1] (adjacent):  drawn from TN+(mu_base, sd_base).
    # psi[k, k+d] for d >= 2:  SST requires psi[k,k+d] >= max over m of
    #   max(psi[k,m], psi[m,k+d]).  We set psi[k,k+d] = lower_bound + gap,
    #   gap ~ TN+(mu_gap, sd_gap), so SST holds exactly and values grow slowly.
    psi_full <- array(0, dim = c(n, K, K))
    for (k in seq_len(K - 1L)) {
      psi_full[, k, k + 1L] <- rtnorm_pos(n, hyper$mu_base, hyper$sd_base)
    }
    if (K > 2L) {
      for (len in 2L:(K - 1L)) {
        for (i in seq_len(K - len)) {
          j <- i + len
          lb <- rep(0, n)
          for (m in (i + 1L):(j - 1L)) {
            lb <- pmax(lb, psi_full[, i, m], psi_full[, m, j])
          }
          psi_full[, i, j] <- lb + rtnorm_pos(n, hyper$mu_gap, hyper$sd_gap)
        }
      }
    }
    out <- matrix(0, nrow = n, ncol = pK)
    for (c in seq_len(pK)) {
      out[, c] <- psi_full[, pair_i[c], pair_j[c]]
    }
    return(out)
  }

  if (family == "Toeplitz_SST") {
    # Draw K-1 distance-specific log-odds phi1 <= phi2 <= ... <= phi_{K-1}.
    phi1 <- rtnorm_pos(n, hyper$mu_phi1, hyper$sd_phi1)
    increments <- matrix(0, nrow = n, ncol = K - 1L)
    increments[, 1] <- phi1
    if (K > 2) {
      for (d in 2:(K - 1L)) {
        increments[, d] <- increments[, d - 1L] + rtnorm_pos(n, hyper$mu_delta, hyper$sd_delta)
      }
    }
    # phi[d] = increments[, d].  psi_{ij} = phi[j - i].
    out <- matrix(0, nrow = n, ncol = pK)
    for (c in seq_len(pK)) {
      d <- pair_j[c] - pair_i[c]   # always > 0 since pair_i < pair_j
      out[, c] <- increments[, d]
    }
    return(out)
  }

  if (family == "LST") {
    # Positive increments induce ordered utilities u1 > u2 > ... > uK.
    increments <- matrix(rtnorm_pos(n * (K - 1L), hyper$mu_a, hyper$sd_a), nrow = n, ncol = K - 1L)
    u <- matrix(0, nrow = n, ncol = K)
    for (k in 2:K) {
      u[, k] <- u[, k - 1L] - increments[, k - 1L]
    }
    out <- matrix(0, nrow = n, ncol = pK)
    for (c in seq_len(pK)) {
      i <- pair_i[c]
      j <- pair_j[c]
      out[, c] <- u[, i] - u[, j]
    }
    return(out)
  }

  stop("Unknown family: ", family)
}

compute_mc_volume_fractions_k_grid <- function(k_grid, cube_max, n_mc = 400000) {
  # For each K, estimate P(u in support) for u ~ Uniform([0, cube_max]^{p_K}).
  # This is the correct dimension-aware measure: it lives in the full p_K-dimensional
  # parameter space, NOT a 3D projection of a fixed triplet.
  #
  # WST:          1 analytically (sampling space = positive orthant = WST support).
  # SST:          Monte Carlo estimate (fraction satisfying all C(K,3) max constraints).
  # LST:           0 exactly (K-1 dimensional manifold, zero Lebesgue measure in R^{p_K}).
  # Toeplitz_SST: 0 exactly (K-1 dimensional manifold, zero Lebesgue measure in R^{p_K}).
  rows <- list()
  for (K in k_grid) {
    pK      <- K * (K - 1L) / 2L
    pairs   <- utils::combn(K, 2L)
    triples <- utils::combn(K, 3L)
    n_triples <- ncol(triples)

    idx_mat <- matrix(NA_integer_, K, K)
    for (c in seq_len(pK)) idx_mat[pairs[1L, c], pairs[2L, c]] <- c

    cat(sprintf("  -> K=%d  p_K=%d  n_triples=%d\n", K, pK, n_triples))

    u <- matrix(stats::runif(as.numeric(n_mc) * pK, 0, cube_max),
                nrow = n_mc, ncol = pK)

    sst_ok <- rep(TRUE, n_mc)
    for (t in seq_len(n_triples)) {
      i <- triples[1L, t]; l <- triples[2L, t]; m <- triples[3L, t]
      sst_ok <- sst_ok &
        (u[, idx_mat[i, m]] >= u[, idx_mat[i, l]]) &
        (u[, idx_mat[i, m]] >= u[, idx_mat[l, m]])
    }
    sst_frac <- mean(sst_ok)
    sst_se   <- sqrt(sst_frac * (1 - sst_frac) / n_mc)

    cat(sprintf("     WST=100%% (analytic)  SST=%.2f%% +/- %.3f%%  [LST, Toeplitz-SST: measure-zero]\n",
                100 * sst_frac, 100 * sst_se))

    rows <- c(rows, list(
      data.frame(K = K, family = "WST",          volume_fraction = 1.0,      pK = pK, intrinsic_dim = pK),
      data.frame(K = K, family = "SST",          volume_fraction = sst_frac, pK = pK, intrinsic_dim = pK),
      data.frame(K = K, family = "Toeplitz_SST", volume_fraction = 0.0,      pK = pK, intrinsic_dim = K - 1L),
      data.frame(K = K, family = "LST",           volume_fraction = 0.0,      pK = pK, intrinsic_dim = K - 1L)
    ))
  }
  do.call(rbind, rows)
}

compute_bt_rmse <- function(u, K, pair_i, pair_j, cube_max) {
  # Least-squares projection onto LST manifold psi_{ij} = v_i - v_j.
  # With complete pair graph, normal equations simplify to v = (u B) / K.
  pK <- ncol(u)
  B <- matrix(0, nrow = pK, ncol = K)
  for (c in seq_len(pK)) {
    B[c, pair_i[c]] <- 1
    B[c, pair_j[c]] <- -1
  }
  v_hat <- (u %*% B) / K
  u_hat <- v_hat %*% t(B)
  sqrt(rowMeans((u - u_hat)^2)) / cube_max
}

weighted_pava <- function(y, w) {
  # Exact weighted isotonic regression (nondecreasing) via PAVA.
  means <- as.numeric(y)
  weights <- as.numeric(w)
  starts <- seq_along(y)
  ends <- seq_along(y)
  nb <- length(means)
  i <- 1L
  while (i < nb) {
    if (means[i] <= means[i + 1L]) {
      i <- i + 1L
    } else {
      new_w <- weights[i] + weights[i + 1L]
      new_m <- (weights[i] * means[i] + weights[i + 1L] * means[i + 1L]) / new_w
      weights[i] <- new_w
      means[i] <- new_m
      ends[i] <- ends[i + 1L]
      if (i + 1L < nb) {
        weights[(i + 1L):(nb - 1L)] <- weights[(i + 2L):nb]
        means[(i + 1L):(nb - 1L)] <- means[(i + 2L):nb]
        starts[(i + 1L):(nb - 1L)] <- starts[(i + 2L):nb]
        ends[(i + 1L):(nb - 1L)] <- ends[(i + 2L):nb]
      }
      nb <- nb - 1L
      if (i > 1L) i <- i - 1L
    }
  }

  out <- numeric(length(y))
  for (b in seq_len(nb)) {
    out[starts[b]:ends[b]] <- means[b]
  }
  out
}

compute_toep_rmse <- function(u, K, pair_i, pair_j, cube_max) {
  # Projection onto monotone distance-based Toeplitz support:
  # psi_{ij} = phi_{j-i}, with 0 <= phi_1 <= ... <= phi_{K-1}.
  pK <- ncol(u)
  d_vec <- pair_j - pair_i

  # Mean by distance class for each row.
  means_d <- matrix(0, nrow = nrow(u), ncol = K - 1L)
  for (d in seq_len(K - 1L)) {
    cols <- which(d_vec == d)
    if (length(cols) == 1L) {
      means_d[, d] <- u[, cols]
    } else {
      means_d[, d] <- rowMeans(u[, cols, drop = FALSE])
    }
  }

  # Exact weighted isotonic projection row-by-row (weights = distance-class counts).
  w_d <- as.numeric(tabulate(d_vec, nbins = K - 1L))
  phi_hat <- matrix(0, nrow = nrow(u), ncol = K - 1L)
  for (r in seq_len(nrow(u))) {
    phi_hat[r, ] <- weighted_pava(means_d[r, ], w_d)
  }
  phi_hat <- pmax(phi_hat, 0)

  u_hat <- matrix(0, nrow = nrow(u), ncol = pK)
  for (c in seq_len(pK)) {
    u_hat[, c] <- phi_hat[, d_vec[c]]
  }
  sqrt(rowMeans((u - u_hat)^2)) / cube_max
}

build_avg_from_neighborhood_indicator <- function(u, indicator, breaks) {
  n_bins <- length(breaks) - 1L
  n_vox <- n_bins^3
  ix <- findInterval(u[, 1], breaks, all.inside = TRUE)
  iy <- findInterval(u[, 2], breaks, all.inside = TRUE)
  iz <- findInterval(u[, 3], breaks, all.inside = TRUE)
  lin <- lin_index_3d(ix, iy, iz, n_bins)
  counts <- tabulate(lin, nbins = n_vox)
  sums <- numeric(n_vox)
  tmp <- tapply(as.numeric(indicator), lin, sum)
  sums[as.integer(names(tmp))] <- as.numeric(tmp)
  avg <- numeric(n_vox)
  keep <- counts > 0
  avg[keep] <- sums[keep] / counts[keep]
  avg
}

compute_unified_k3_visuals <- function(cube_max, breaks, eps_vis = 0.05, n_mc_vis = 400000) {
  K <- 3L
  pair_info <- build_pair_index(K)
  u <- matrix(stats::runif(as.numeric(n_mc_vis) * 3L, 0, cube_max), nrow = n_mc_vis, ncol = 3L)

  dist_wst <- rep(0, n_mc_vis)
  dist_sst <- compute_sst_violation(u, K, pair_info$pair_i, pair_info$pair_j) / cube_max
  dist_toep <- compute_toep_rmse(u, K, pair_info$pair_i, pair_info$pair_j, cube_max)
  dist_bt <- compute_bt_rmse(u, K, pair_info$pair_i, pair_info$pair_j, cube_max)

  avg_cov_list <- list(
    WST = build_avg_from_neighborhood_indicator(u, dist_wst <= eps_vis, breaks),
    SST = build_avg_from_neighborhood_indicator(u, dist_sst <= eps_vis, breaks),
    Toeplitz_SST = build_avg_from_neighborhood_indicator(u, dist_toep <= eps_vis, breaks),
    LST = build_avg_from_neighborhood_indicator(u, dist_bt <= eps_vis, breaks)
  )

  list(avg_cov_list = avg_cov_list, eps_vis = eps_vis, n_mc_vis = n_mc_vis)
}

compute_sst_violation <- function(u, K, pair_i, pair_j) {
  # For each row, return max normalized SST deficit:
  # deficit(i,l,m) = max(0, max(psi_il, psi_lm) - psi_im).
  pK <- ncol(u)
  idx_mat <- matrix(NA_integer_, nrow = K, ncol = K)
  for (c in seq_len(pK)) idx_mat[pair_i[c], pair_j[c]] <- c

  triples <- utils::combn(K, 3L)
  max_def <- rep(0, nrow(u))
  for (t in seq_len(ncol(triples))) {
    i <- triples[1L, t]
    l <- triples[2L, t]
    m <- triples[3L, t]
    def_t <- pmax(0, pmax(u[, idx_mat[i, l]], u[, idx_mat[l, m]]) - u[, idx_mat[i, m]])
    max_def <- pmax(max_def, def_t)
  }
  max_def
}

compute_epsilon_neighborhood_k_grid <- function(k_grid, cube_max, n_mc = 100000,
                                                eps_rel = c(0.02, 0.05, 0.10)) {
  # Compare restrictiveness via epsilon-neighborhood volume under a common law:
  # U ~ Uniform([0, cube_max]^{p_K}). For each model M, estimate
  #   P( dist(U, M) / cube_max <= eps_rel ).
  # This yields non-zero, comparable quantities for SST, Toeplitz-SST, and LST.
  rows <- list()
  idx <- 1L
  for (K in k_grid) {
    pair_info <- build_pair_index(K)
    pK <- K * (K - 1L) / 2L
    cat(sprintf("  -> epsilon-neighborhoods: K=%d  p_K=%d\n", K, pK))

    u <- matrix(stats::runif(as.numeric(n_mc) * pK, 0, cube_max), nrow = n_mc, ncol = pK)

    sst_dist <- compute_sst_violation(u, K, pair_info$pair_i, pair_info$pair_j) / cube_max
    toep_dist <- compute_toep_rmse(u, K, pair_info$pair_i, pair_info$pair_j, cube_max)
    bt_dist <- compute_bt_rmse(u, K, pair_info$pair_i, pair_info$pair_j, cube_max)

    dist_map <- list(SST = sst_dist, Toeplitz_SST = toep_dist, LST = bt_dist)
    for (fam in names(dist_map)) {
      d <- dist_map[[fam]]
      for (eps in eps_rel) {
        n_in <- sum(d <= eps)
        # Jeffreys-style smoothing avoids exact zeros in rare-event regimes.
        share_smooth <- (n_in + 0.5) / (length(d) + 1)
        rows[[idx]] <- data.frame(
          K = K,
          pK = pK,
          family = fam,
          epsilon_rel = eps,
          neighborhood_share = share_smooth,
          median_dist = as.numeric(stats::quantile(d, 0.50)),
          q90_dist = as.numeric(stats::quantile(d, 0.90)),
          mean_dist = mean(d)
        )
        idx <- idx + 1L
      }
    }
  }
  do.call(rbind, rows)
}

summarise_restrictiveness <- function(avg_cov, threshold) {
  c(
    mean_coverage = mean(avg_cov),
    q95_coverage = as.numeric(stats::quantile(avg_cov, probs = 0.95)),
    support_share = mean(avg_cov >= threshold)
  )
}

build_slice_df <- function(centers, avg_cov, n_bins, cube_max, slice_targets = NULL) {
  # Build 2D cross-sections on psi13 = constant planes.
  if (is.null(slice_targets)) {
    slice_targets <- c(0.20, 0.45, 0.70) * cube_max
  }

  y_vals <- sort(unique(centers$psi13))
  out <- list()
  for (s in slice_targets) {
    y_pick <- y_vals[which.min(abs(y_vals - s))]
    idx <- which(abs(centers$psi13 - y_pick) < 1e-12)
    df <- centers[idx, c("psi12", "psi23")]
    df$avg_cov <- avg_cov[idx]
    df$slice_value <- y_pick
    df$slice_label <- sprintf("psi13 = %.2f", y_pick)
    out[[length(out) + 1L]] <- df
  }
  do.call(rbind, out)
}

make_plotly_panel <- function(centers, avg_cov, title, threshold) {
  keep <- avg_cov >= threshold
  df <- centers[keep, c("psi12", "psi13", "psi23")]
  if (nrow(df) == 0) {
    # Fallback if threshold is too high.
    top_idx <- order(avg_cov, decreasing = TRUE)[seq_len(min(1200, length(avg_cov)))]
    df <- centers[top_idx, c("psi12", "psi13", "psi23")]
    vals <- avg_cov[top_idx]
  } else {
    vals <- avg_cov[keep]
  }

  plotly::plot_ly(
    data = df,
    x = ~psi12, y = ~psi13, z = ~psi23,
    type = "scatter3d", mode = "markers",
    marker = list(
      size = 2.3,
      opacity = 0.72,
      color = vals,
      colorscale = "Viridis",
      cmin = 0,
      cmax = 1,
      colorbar = list(title = "Coverage")
    ),
    hovertemplate = paste0(
      "psi12=%{x:.2f}<br>",
      "psi13=%{y:.2f}<br>",
      "psi23=%{z:.2f}<br>",
      "coverage=%{marker.color:.3f}<extra></extra>"
    )
  ) |>
    plotly::layout(
      title = list(text = title, y = 0.96),
      scene = list(
        xaxis = list(title = "psi12", range = c(0, cube_max)),
        yaxis = list(title = "psi13", range = c(0, cube_max)),
        zaxis = list(title = "psi23", range = c(0, cube_max)),
        aspectmode = "cube"
      ),
      margin = list(l = 0, r = 0, b = 0, t = 32)
    )
}

build_support_points_for_overlay <- function(centers, avg_cov, family, threshold) {
  keep <- avg_cov >= threshold
  df <- centers[keep, c("psi12", "psi13", "psi23")]
  if (nrow(df) == 0) {
    top_idx <- order(avg_cov, decreasing = TRUE)[seq_len(min(1200, length(avg_cov)))]
    df <- centers[top_idx, c("psi12", "psi13", "psi23")]
    cov_vals <- avg_cov[top_idx]
  } else {
    cov_vals <- avg_cov[keep]
  }
  df$coverage <- cov_vals
  df$family <- family
  df
}

save_static_3d_overlay <- function(overlay_df, file, cube_max) {
  family_levels <- c("WST", "SST", "Toeplitz_SST", "LST")
  family_colors <- c(
    WST = "#4C78A8",
    SST = "#F58518",
    Toeplitz_SST = "#54A24B",
    LST = "#B279A2"
  )

  # Subsample for readability and speed.
  set.seed(seed + 101)
  max_per_family <- 2200
  sub_df <- do.call(rbind, lapply(family_levels, function(fam) {
    d <- overlay_df[overlay_df$family == fam, , drop = FALSE]
    if (nrow(d) <= max_per_family) return(d)
    d[sample.int(nrow(d), max_per_family), , drop = FALSE]
  }))

  # Perspective transform.
  pmat <- graphics::persp(
    x = c(0, cube_max),
    y = c(0, cube_max),
    z = matrix(c(0, cube_max, cube_max, 0), nrow = 2, ncol = 2),
    zlim = c(0, cube_max),
    theta = 38, phi = 24, r = 2.2, expand = 0.95,
    col = "#f8f8f8", border = "#d0d0d0",
    ticktype = "detailed",
    xlab = "psi12", ylab = "psi23", zlab = "psi13"
  )

  # Draw points from far to near for cleaner occlusion.
  depth_key <- sub_df$psi12 + sub_df$psi13 + sub_df$psi23
  ord <- order(depth_key, decreasing = FALSE)
  sub_df <- sub_df[ord, , drop = FALSE]

  proj <- grDevices::trans3d(
    x = sub_df$psi12,
    y = sub_df$psi23,
    z = sub_df$psi13,
    pmat = pmat
  )

  base_col <- unname(family_colors[as.character(sub_df$family)])
  alpha_vec <- pmin(0.75, pmax(0.20, 0.15 + sub_df$coverage))
  point_col <- mapply(
    function(cl, a) grDevices::adjustcolor(cl, alpha.f = a),
    cl = base_col,
    a = alpha_vec,
    USE.NAMES = FALSE
  )

  points(proj$x, proj$y, pch = 16, cex = 0.40, col = point_col)

  legend(
    "topright",
    legend = family_levels,
    col = family_colors[family_levels],
    pch = 16,
    pt.cex = 1.1,
    bty = "n",
    title = "Support"
  )
  title(main = "Static 3D support coverage overlay")
  mtext("Each point = cube cell covered in at least threshold fraction of hyperparameter settings", side = 1, line = 3, cex = 0.85)
}

save_static_3d_small_multiples <- function(overlay_df, cube_max, support_share_map = NULL) {
  family_levels <- c("WST", "SST", "Toeplitz_SST", "LST")
  family_labels <- c(
    WST = "WST",
    SST = "SST",
    Toeplitz_SST = "Toepliz SST",
    LST = "LST"
  )
  family_colors <- c(
    WST = "#4C78A8",
    SST = "#F58518",
    Toeplitz_SST = "#54A24B",
    LST = "#D81B60"
  )

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 2), mar = c(2.1, 2.2, 2.4, 0.8), oma = c(0.5, 0.5, 1.8, 0.5))

  for (fam in family_levels) {
    d <- overlay_df[overlay_df$family == fam, , drop = FALSE]
    if (nrow(d) == 0) next

    set.seed(seed + 202)
    if (nrow(d) > 2500) d <- d[sample.int(nrow(d), 2500), , drop = FALSE]

    # Per-panel viewing angle: lift phi for LST to reveal the triangular face.
    theta_use <- 36
    phi_use   <- if (fam == "LST") 35 else 22

    pmat <- graphics::persp(
      x = c(0, cube_max),
      y = c(0, cube_max),
      z = matrix(c(0, cube_max, cube_max, 0), nrow = 2, ncol = 2),
      zlim = c(0, cube_max),
      theta = theta_use, phi = phi_use, r = 2.2, expand = 0.95,
      col = "#ffffff", border = NA,
      ticktype = "simple",
      xlab = "psi12", ylab = "psi23", zlab = "psi13"
    )

    proj <- grDevices::trans3d(d$psi12, d$psi23, d$psi13, pmat)
    alpha_use <- if (fam == "LST") 0.62 else 0.48
    cex_use <- if (fam == "LST") 0.48 else 0.42
    base_col <- grDevices::adjustcolor(family_colors[fam], alpha.f = alpha_use)
    points(proj$x, proj$y, pch = 16, cex = cex_use, col = base_col)

    # --- Wireframe edges and vertex annotations ---
    B <- cube_max
    vtx3d <- function(px, py, pz) grDevices::trans3d(px, py, pz, pmat)
    draw_edge3d <- function(x0, y0, z0, x1, y1, z1,
                            col = "#444444", lty = 1, lwd = 1.15) {
      p0 <- vtx3d(x0, y0, z0)
      p1 <- vtx3d(x1, y1, z1)
      lines(c(p0$x, p1$x), c(p0$y, p1$y), col = col, lty = lty, lwd = lwd)
    }
    mark_vtx <- function(px, py, pz, lbl, adj = c(0.5, 0.5), cex = 0.58) {
      p <- vtx3d(px, py, pz)
      points(p$x, p$y, pch = 21, bg = "white", col = "#333333", cex = 0.80)
      text(p$x, p$y, labels = lbl, adj = adj, cex = cex, font = 2)
    }
    ec <- grDevices::adjustcolor(family_colors[fam], alpha.f = 0.90)

    if (fam == "WST") {
      # Whole cube: just anchor the two extremes.
      mark_vtx(0, 0, 0, "(0,0,0)", adj = c(1.15, 1.30), cex = 0.55)
      mark_vtx(B, B, B, "(B,B,B)", adj = c(-0.15, -0.25), cex = 0.55)

    } else if (fam == "SST") {
      # SST polytope vertices: (0,0,0),(0,0,B),(B,0,B),(0,B,B),(B,B,B).
      # Boundary edges of the pyramid.
      draw_edge3d(0, 0, 0, 0, 0, B, ec, lwd = 1.2)          # vertical edge
      draw_edge3d(0, 0, 0, B, 0, B, ec, lwd = 1.2)          # ridge psi13=psi12
      draw_edge3d(0, 0, 0, 0, B, B, ec, lwd = 1.2)          # ridge psi13=psi23
      draw_edge3d(0, 0, B, B, 0, B, ec, lty = 2, lwd = 1.0) # top face
      draw_edge3d(0, 0, B, 0, B, B, ec, lty = 2, lwd = 1.0) # top face
      draw_edge3d(B, 0, B, B, B, B, ec, lwd = 1.2)
      draw_edge3d(0, B, B, B, B, B, ec, lwd = 1.2)
      mark_vtx(0, 0, 0, "(0,0,0)", adj = c(1.15, 1.30), cex = 0.55)
      mark_vtx(0, 0, B, "(0,0,B)", adj = c(1.20, 0.50), cex = 0.55)
      mark_vtx(B, 0, B, "(B,0,B)", adj = c(-0.15, 0.50), cex = 0.55)
      mark_vtx(0, B, B, "(0,B,B)", adj = c(1.15, 0.50), cex = 0.55)
      mark_vtx(B, B, B, "(B,B,B)", adj = c(-0.15, -0.25), cex = 0.55)

    } else if (fam == "Toeplitz_SST") {
      # Triangle: (0,0,0), (0,0,B), (B,B,B).
      draw_edge3d(0, 0, 0, 0, 0, B, ec, lwd = 1.5)
      draw_edge3d(0, 0, B, B, B, B, ec, lwd = 1.5)
      draw_edge3d(0, 0, 0, B, B, B, ec, lwd = 1.5)
      mark_vtx(0, 0, 0, "(0,0,0)", adj = c(1.15, 1.30), cex = 0.58)
      mark_vtx(0, 0, B, "(0,0,B)", adj = c(1.20, 0.50), cex = 0.58)
      mark_vtx(B, B, B, "(B,B,B)", adj = c(-0.15, -0.25), cex = 0.58)

    } else if (fam == "LST") {
      # Triangle: (0,0,0), (B,0,B), (0,B,B).
      draw_edge3d(0, 0, 0, B, 0, B, ec, lwd = 1.5)
      draw_edge3d(0, 0, 0, 0, B, B, ec, lwd = 1.5)
      draw_edge3d(B, 0, B, 0, B, B, ec, lwd = 1.5) # top edge = simplex at psi13=B
      mark_vtx(0, 0, 0, "(0,0,0)", adj = c(1.15, 1.30), cex = 0.58)
      mark_vtx(B, 0, B, "(B,0,B)", adj = c(-0.15, 0.50), cex = 0.58)
      mark_vtx(0, B, B, "(0,B,B)", adj = c(1.15, 0.50), cex = 0.58)
    }

    title(main = family_labels[fam])

    if (!is.null(support_share_map) && fam %in% names(support_share_map)) {
      pct <- 100 * support_share_map[[fam]]
      usr <- par("usr")
      graphics::text(
        x = usr[1] + 0.03 * (usr[2] - usr[1]),
        y = usr[4] - 0.05 * (usr[4] - usr[3]),
        labels = sprintf("Covered: %.1f%%", pct),
        adj = c(0, 1),
        cex = 0.92,
        col = "#202020",
        font = 2
      )
    }
  }

  mtext("Support Geometry In Psi-Space By Model Assumption", outer = TRUE, line = 0.3, cex = 1.15, font = 2)
}

# Clean geometric version: shaded filled polygons, no scatter points, phi=35 for all panels.
save_static_3d_shaded_geometry <- function(cube_max) {
  B <- cube_max
  family_levels <- c("WST", "SST", "Toeplitz_SST", "LST")
  family_labels  <- c(
    WST = "WST",
    SST = "SST",
    Toeplitz_SST = "Toeplitz SST",
    LST = "LST"
  )
  family_colors <- c(
    WST          = "#4C78A8",
    SST          = "#F58518",
    Toeplitz_SST = "#54A24B",
    LST           = "#D81B60"
  )

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 2), mar = c(2.1, 2.2, 2.8, 0.8), oma = c(0.5, 0.5, 1.8, 0.5))

  for (fam in family_levels) {
    col_edge  <- grDevices::adjustcolor(family_colors[fam], alpha.f = 0.95)
    col_fill  <- grDevices::adjustcolor(family_colors[fam], alpha.f = 0.28)
    col_fill2 <- grDevices::adjustcolor(family_colors[fam], alpha.f = 0.14)

    pmat <- graphics::persp(
      x = c(0, B), y = c(0, B),
      z = matrix(c(0, B, B, 0), nrow = 2),
      zlim = c(0, B),
      theta = 36, phi = 35, r = 2.2, expand = 0.95,
      col = NA, border = NA,
      ticktype = "simple",
      xlab = "psi12", ylab = "psi23", zlab = "psi13"
    )

    p3d <- function(x, y, z) grDevices::trans3d(x, y, z, pmat)

    fill_face <- function(xs, ys, zs, fill = col_fill, lwd = 1.4) {
      pts <- p3d(xs, ys, zs)
      graphics::polygon(pts$x, pts$y, col = fill, border = col_edge, lwd = lwd)
    }
    draw_edge3d <- function(x0, y0, z0, x1, y1, z1, lty = 1, lwd = 1.4) {
      p0 <- p3d(x0, y0, z0); p1 <- p3d(x1, y1, z1)
      graphics::lines(c(p0$x, p1$x), c(p0$y, p1$y), col = col_edge, lty = lty, lwd = lwd)
    }
    mark_vtx <- function(px, py, pz, lbl, adj = c(0.5, 0.5), cex = 0.60) {
      p <- p3d(px, py, pz)
      graphics::points(p$x, p$y, pch = 21, bg = "white", col = "#333333", cex = 0.85)
      graphics::text(p$x, p$y, labels = lbl, adj = adj, cex = cex, font = 2)
    }

    if (fam == "WST") {
      # Three pairs of visible/back cube faces, back-to-front.
      fill_face(c(0, 0, 0, 0), c(0, B, B, 0), c(0, 0, B, B), fill = col_fill2) # x=0
      fill_face(c(0, B, B, 0), c(0, 0, 0, 0), c(0, 0, B, B), fill = col_fill2) # y=0
      fill_face(c(0, B, B, 0), c(0, 0, B, B), c(0, 0, 0, 0), fill = col_fill2) # z=0
      fill_face(c(B, B, B, B), c(0, B, B, 0), c(0, 0, B, B))                   # x=B
      fill_face(c(0, B, B, 0), c(B, B, B, B), c(0, 0, B, B))                   # y=B
      fill_face(c(0, B, B, 0), c(0, 0, B, B), c(B, B, B, B))                   # z=B (top)
      mark_vtx(0, 0, 0, "(0,0,0)", adj = c(1.15, 1.30))
      mark_vtx(B, B, B, "(B,B,B)", adj = c(-0.15, -0.25))

    } else if (fam == "SST") {
      # 5 faces drawn back-to-front.
      fill_face(c(0, 0, B),    c(0, 0, 0),    c(0, B, B),    fill = col_fill2) # front wall y=0
      fill_face(c(0, 0, 0),    c(0, 0, B),    c(0, B, B),    fill = col_fill2) # left wall  x=0
      fill_face(c(0, B, B),    c(0, 0, B),    c(0, B, B))                      # right slope z=x
      fill_face(c(0, 0, B),    c(0, B, B),    c(0, B, B))                      # back slope  z=y
      fill_face(c(0, B, B, 0), c(0, 0, B, B), c(B, B, B, B))                  # top (z=B)
      mark_vtx(0, 0, 0, "(0,0,0)", adj = c(1.15, 1.30))
      mark_vtx(0, 0, B, "(0,0,B)", adj = c(1.20, 0.55))
      mark_vtx(B, 0, B, "(B,0,B)", adj = c(-0.15, 0.55))
      mark_vtx(0, B, B, "(0,B,B)", adj = c(1.15, 0.55))
      mark_vtx(B, B, B, "(B,B,B)", adj = c(-0.15, -0.25))

    } else if (fam == "Toeplitz_SST") {
      # Triangle: (0,0,0), (0,0,B), (B,B,B).
      fill_face(c(0, 0, B), c(0, 0, B), c(0, B, B))
      mark_vtx(0, 0, 0, "(0,0,0)", adj = c(1.15, 1.30))
      mark_vtx(0, 0, B, "(0,0,B)", adj = c(1.20, 0.55))
      mark_vtx(B, B, B, "(B,B,B)", adj = c(-0.15, -0.25))

    } else if (fam == "LST") {
      # Triangle: (0,0,0), (B,0,B), (0,B,B).
      fill_face(c(0, B, 0), c(0, 0, B), c(0, B, B))
      # Redraw edges so they sit on top of fill.
      draw_edge3d(0, 0, 0, B, 0, B)
      draw_edge3d(0, 0, 0, 0, B, B)
      draw_edge3d(B, 0, B, 0, B, B)
      mark_vtx(0, 0, 0, "(0,0,0)", adj = c(1.15, 1.30))
      mark_vtx(B, 0, B, "(B,0,B)", adj = c(-0.15, 0.55))
      mark_vtx(0, B, B, "(0,B,B)", adj = c(1.15, 0.55))
    }

    title(main = family_labels[fam], cex.main = 1.05)
  }

  mtext("Support Geometry In Psi-Space By Model Assumption", outer = TRUE, line = 0.3, cex = 1.15, font = 2)
}

build_projection_long <- function(overlay_df) {
  rbind(
    data.frame(
      projection = "psi12 vs psi13",
      x = overlay_df$psi12,
      y = overlay_df$psi13,
      family = overlay_df$family,
      coverage = overlay_df$coverage
    ),
    data.frame(
      projection = "psi12 vs psi23",
      x = overlay_df$psi12,
      y = overlay_df$psi23,
      family = overlay_df$family,
      coverage = overlay_df$coverage
    ),
    data.frame(
      projection = "psi13 vs psi23",
      x = overlay_df$psi13,
      y = overlay_df$psi23,
      family = overlay_df$family,
      coverage = overlay_df$coverage
    )
  )
}

analyse_numeric_tendencies <- function(centers, avg_cov_list, threshold) {
  families_local <- names(avg_cov_list)
  out <- list()
  for (fam in families_local) {
    w <- avg_cov_list[[fam]]
    keep <- w > 0
    df <- centers[keep, c("psi12", "psi13", "psi23")]
    ww <- w[keep]

    # Weighted correlations and linear tendencies.
    cor_xz <- stats::cov.wt(cbind(df$psi12, df$psi23), wt = ww, cor = TRUE)$cor[1, 2]
    cor_xy <- stats::cov.wt(cbind(df$psi12, df$psi13), wt = ww, cor = TRUE)$cor[1, 2]
    cor_zy <- stats::cov.wt(cbind(df$psi23, df$psi13), wt = ww, cor = TRUE)$cor[1, 2]

    lm_xz <- stats::lm(psi23 ~ psi12, data = df, weights = ww)
    lm_xy <- stats::lm(psi13 ~ psi12, data = df, weights = ww)

    # Structural diagnostics implied by each assumption.
    d_bt <- df$psi13 - (df$psi12 + df$psi23)
    d_toep <- df$psi12 - df$psi23
    d_sst <- df$psi13 - pmax(df$psi12, df$psi23)

    wmean_abs <- function(x, w) sum(abs(x) * w) / sum(w)
    wmean_signed <- function(x, w) sum(x * w) / sum(w)

    out[[fam]] <- data.frame(
      family = fam,
      support_share = mean(w >= threshold),
      corr_psi12_psi23 = cor_xz,
      corr_psi12_psi13 = cor_xy,
      corr_psi23_psi13 = cor_zy,
      slope_psi23_on_psi12 = unname(stats::coef(lm_xz)[2]),
      slope_psi13_on_psi12 = unname(stats::coef(lm_xy)[2]),
      mean_abs_bt_gap = wmean_abs(d_bt, ww),
      mean_signed_bt_gap = wmean_signed(d_bt, ww),
      mean_abs_toep_gap = wmean_abs(d_toep, ww),
      mean_abs_sst_gap = wmean_abs(d_sst, ww),
      share_bt_like = sum((abs(d_bt) <= 0.20) * ww) / sum(ww),
      share_toep_like = sum((abs(d_toep) <= 0.15 & d_sst >= -1e-10) * ww) / sum(ww)
    )
  }
  do.call(rbind, out)
}

# -----------------------------
# Main simulation
# -----------------------------
cat("Running support-geometry cube simulation...\n")
cat(sprintf("seed=%d | K_sim=%d | cube_max=%.2f | n_bins=%d | n_param_settings=%d | n_samples=%d\n",
            seed, K_sim, cube_max, n_bins, n_param_settings, n_samples))

families <- c("WST", "SST", "Toeplitz_SST", "LST")
breaks <- seq(0, cube_max, length.out = n_bins + 1L)
centers <- build_centers_df(breaks, n_bins)

avg_cov_list <- list()
summary_rows <- list()

for (fam in families) {
  cat(sprintf("  -> Family: %s\n", fam))
  occ_sum <- integer(n_bins^3)
  n_nonempty <- 0L

  for (s in seq_len(n_param_settings)) {
    hyper <- random_hyper(fam)
    pts <- sample_support(fam, n_samples, hyper, cube_max, K_sim)
    if (nrow(pts) > 0) {
      occ_sum <- occ_sum + voxel_occupancy(pts, breaks)
      n_nonempty <- n_nonempty + 1L
    }
  }

  if (n_nonempty == 0L) {
    avg_cov <- rep(0, n_bins^3)
  } else {
    avg_cov <- occ_sum / n_nonempty
  }

  avg_cov_list[[fam]] <- avg_cov
  stats_vec <- summarise_restrictiveness(avg_cov, occ_threshold)
  summary_rows[[fam]] <- data.frame(
    family = fam,
    n_param_settings = n_param_settings,
    n_nonempty_settings = n_nonempty,
    mean_coverage = unname(stats_vec["mean_coverage"]),
    q95_coverage = unname(stats_vec["q95_coverage"]),
    support_share = unname(stats_vec["support_share"])
  )
}

summary_df <- do.call(rbind, summary_rows)
summary_df$family <- factor(summary_df$family, levels = families)
summary_df <- summary_df[order(summary_df$support_share, decreasing = TRUE), ]

# -----------------------------
# Outputs: tables and plots
# -----------------------------
summary_csv <- file.path(out_dir, "support_restrictiveness_summary.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)

# Bar plot of support share.
bar_png <- file.path(out_dir, "support_share_barplot.png")
bar_png_legacy <- file.path(out_dir, "support_effective_volume_barplot.png")
g <- ggplot2::ggplot(
  summary_df,
  ggplot2::aes(x = reorder(family, -support_share), y = support_share, fill = family)
) +
  ggplot2::geom_col(width = 0.72, alpha = 0.9) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f%%", 100 * support_share)),
    vjust = -0.3,
    size = 3.8
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, min(1, max(summary_df$support_share) * 1.12)),
    labels = function(x) sprintf("%.0f%%", 100 * x)
  ) +
  ggplot2::labs(
    title = sprintf("Average support coverage in [0, %.1f]^3 (projection, K=%d)", cube_max, K_sim),
    subtitle = sprintf("Support share = fraction of cube cells covered in at least %.0f%% of parameter settings", 100 * occ_threshold),
    x = NULL,
    y = "Support share"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(bar_png, g, width = 8.0, height = 5.0, dpi = 160)
ggplot2::ggsave(bar_png_legacy, g, width = 8.0, height = 5.0, dpi = 160)

# Build 2x2 interactive 3D panel.
plt_wst <- make_plotly_panel(centers, avg_cov_list[["WST"]], "WST", occ_threshold)
plt_sst <- make_plotly_panel(centers, avg_cov_list[["SST"]], "SST", occ_threshold)
plt_toep <- make_plotly_panel(centers, avg_cov_list[["Toeplitz_SST"]], "Toeplitz-SST", occ_threshold)
plt_bt <- make_plotly_panel(centers, avg_cov_list[["LST"]], "LST", occ_threshold)

panel <- plotly::subplot(
  plt_wst, plt_sst, plt_toep, plt_bt,
  nrows = 2,
  margin = 0.03,
  titleX = TRUE,
  titleY = TRUE
) |>
  plotly::layout(
    title = list(
      text = sprintf("Average 3D support coverage (projection, K = %d)", K_sim),
      x = 0.5
    ),
    showlegend = FALSE
  )

html_file <- file.path(out_dir, "average_3d_supports_plotly.html")
has_pandoc <- FALSE
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  has_pandoc <- isTRUE(rmarkdown::pandoc_available())
}
htmlwidgets::saveWidget(panel, html_file, selfcontained = has_pandoc)

# Also save the raw averaged coverage arrays for reproducibility.
rds_file <- file.path(out_dir, "average_coverage_arrays.rds")
saveRDS(
  list(
    config = list(
      seed = seed,
      K_sim = K_sim,
      cube_max = cube_max,
      n_bins = n_bins,
      n_param_settings = n_param_settings,
      n_samples = n_samples,
      occ_threshold = occ_threshold
    ),
    voxel_centers = centers[, c("psi12", "psi13", "psi23")],
    avg_coverage = avg_cov_list,
    summary = summary_df
  ),
  rds_file
)

# Static cross-sectional plot designed to better distinguish supports.
slice_df_list <- list()
for (fam in families) {
  tmp <- build_slice_df(
    centers = centers,
    avg_cov = avg_cov_list[[fam]],
    n_bins = n_bins,
    cube_max = cube_max
  )
  tmp$family <- fam
  slice_df_list[[fam]] <- tmp
}
slice_df <- do.call(rbind, slice_df_list)
slice_df$family <- factor(slice_df$family, levels = families)
slice_df$family_label <- factor(
  as.character(slice_df$family),
  levels = c("WST", "SST", "Toeplitz_SST", "LST"),
  labels = c("WST", "SST", "Toeplitz SST", "LST")
)

panel_share <- data.frame(
  family = as.character(summary_df$family),
  support_share = summary_df$support_share,
  stringsAsFactors = FALSE
)
slice_df <- merge(slice_df, panel_share, by = "family", all.x = TRUE)
slice_df$panel_note <- sprintf("Covered: %.1f%%", 100 * slice_df$support_share)

slice_png <- file.path(out_dir, "support_cross_sections_static.png")
g_slice <- ggplot2::ggplot(slice_df, ggplot2::aes(x = psi12, y = psi23, fill = avg_cov)) +
  ggplot2::geom_raster(interpolate = FALSE) +
  ggplot2::geom_contour(
    ggplot2::aes(z = avg_cov),
    breaks = c(0.15, 0.35, 0.55, 0.75),
    color = "white",
    linewidth = 0.24,
    alpha = 0.78
  ) +
  ggplot2::geom_text(
    data = unique(slice_df[, c("family_label", "slice_label", "panel_note")]),
    mapping = ggplot2::aes(x = 0.10 * cube_max, y = 0.93 * cube_max, label = panel_note),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    color = "white",
    fontface = "bold"
  ) +
  ggplot2::facet_grid(slice_label ~ family_label) +
  ggplot2::scale_fill_gradientn(
    colours = c("#fffdf5", "#fee8a8", "#fdd072", "#fca85e", "#f67e4b", "#e34a33"),
    limits = c(0, 1),
    name = "Coverage"
  ) +
  ggplot2::scale_x_continuous(breaks = 0:4, minor_breaks = NULL, expand = c(0, 0)) +
  ggplot2::scale_y_continuous(breaks = 0:4, minor_breaks = NULL, expand = c(0, 0)) +
  ggplot2::coord_equal(xlim = c(0, cube_max), ylim = c(0, cube_max), expand = FALSE) +
  ggplot2::labs(
    title = "Support Geometry Cross-Sections In Psi-Space",
    subtitle = "Rows fix psi13; columns are model assumptions. Brighter cells indicate higher support coverage.",
    x = "psi12",
    y = "psi23"
  ) +
  ggplot2::theme_minimal(base_size = 11.5) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    strip.background = ggplot2::element_rect(fill = "#f2f2f2", color = "#d9d9d9", linewidth = 0.3),
    panel.border = ggplot2::element_rect(color = "#e2e2e2", fill = NA, linewidth = 0.35),
    axis.text = ggplot2::element_text(color = "#3a3a3a"),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(color = "#4a4a4a"),
    panel.spacing = grid::unit(0.4, "lines")
  )

ggplot2::ggsave(slice_png, g_slice, width = 12, height = 8.5, dpi = 180)

# Static non-interactive 3D overlay (all supports in one cube).
overlay_df <- do.call(rbind, lapply(families, function(fam) {
  build_support_points_for_overlay(centers, avg_cov_list[[fam]], fam, occ_threshold)
}))
overlay_png <- file.path(out_dir, "support_3d_overlay_static.png")
grDevices::png(overlay_png, width = 1900, height = 1300, res = 180)
par(mar = c(2.5, 2.8, 2.8, 1.2))
save_static_3d_overlay(overlay_df, overlay_png, cube_max)
dev.off()

# Readability approach 1: 3D small-multiples (one support per panel).
overlay_small_png <- file.path(out_dir, "support_3d_small_multiples_static.png")
grDevices::png(overlay_small_png, width = 1900, height = 1500, res = 180)
share_map <- stats::setNames(summary_df$support_share, as.character(summary_df$family))
save_static_3d_small_multiples(overlay_df, cube_max, support_share_map = share_map)
dev.off()

# Geometric shaded version: filled polygons only, no scatter points.
geom_clean_png <- file.path(out_dir, "support_3d_shaded_geometry.png")
grDevices::png(geom_clean_png, width = 1900, height = 1500, res = 180)
save_static_3d_shaded_geometry(cube_max)
dev.off()
message(" - ", geom_clean_png)

# Readability approach 2: 2D projection contour overlays.
proj_df <- build_projection_long(overlay_df)
proj_png <- file.path(out_dir, "support_projection_contours_static.png")
g_proj <- ggplot2::ggplot(
  proj_df,
  ggplot2::aes(x = x, y = y, color = family, group = family)
) +
  ggplot2::stat_density_2d(
    geom = "path",
    contour_var = "ndensity",
    breaks = c(0.20, 0.40, 0.60, 0.80),
    linewidth = 0.70,
    alpha = 0.95
  ) +
  ggplot2::coord_equal(xlim = c(0, cube_max), ylim = c(0, cube_max), expand = FALSE) +
  ggplot2::facet_wrap(~projection, nrow = 1) +
  ggplot2::scale_color_manual(values = c(
    WST = "#4C78A8",
    SST = "#F58518",
    Toeplitz_SST = "#54A24B",
    LST = "#B279A2"
  )) +
  ggplot2::labs(
    title = "Support coverage contours across 2D projections",
    subtitle = "Readability view: same color = same assumption across projections",
    x = NULL,
    y = NULL,
    color = "Assumption"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(proj_png, g_proj, width = 12.5, height = 4.9, dpi = 190)

# Numeric tendencies for interpretation and reporting.
tendency_df <- analyse_numeric_tendencies(centers, avg_cov_list, occ_threshold)
tendency_csv <- file.path(out_dir, "support_numeric_tendencies.csv")
write.csv(tendency_df, tendency_csv, row.names = FALSE)

cat("\nNumeric tendency summary:\n")
print(data.frame(
  family = tendency_df$family,
  support_share = sprintf("%.2f%%", 100 * tendency_df$support_share),
  corr_psi12_psi23 = sprintf("%.3f", tendency_df$corr_psi12_psi23),
  slope_psi23_on_psi12 = sprintf("%.3f", tendency_df$slope_psi23_on_psi12),
  mean_abs_bt_gap = sprintf("%.3f", tendency_df$mean_abs_bt_gap),
  mean_abs_toep_gap = sprintf("%.3f", tendency_df$mean_abs_toep_gap)
), row.names = FALSE)

# Barplot content table: print, save as CSV, and save as LaTeX.
bar_table <- summary_df[, c("family", "support_share")]
names(bar_table) <- c("Model", "SupportShare")
bar_table$SupportSharePct <- 100 * bar_table$SupportShare
bar_table <- bar_table[, c("Model", "SupportSharePct")]

cat("\nBarplot content table (support share):\n")
print(data.frame(
  Model = bar_table$Model,
  SupportSharePct = sprintf("%.2f", bar_table$SupportSharePct)
), row.names = FALSE)

bar_table_csv <- file.path(out_dir, "support_share_barplot_table.csv")
write.csv(
  data.frame(
    Model = bar_table$Model,
    SupportSharePct = round(bar_table$SupportSharePct, 3)
  ),
  bar_table_csv,
  row.names = FALSE
)

bar_table_tex <- file.path(out_dir, "support_share_barplot_table.tex")
tex_lines <- c(
  "\\begin{tabular}{lr}",
  "\\hline",
  "Model & Support share (\\%) \\\\",
  "\\hline"
)
for (i in seq_len(nrow(bar_table))) {
  model_label <- gsub("_", "\\\\_", as.character(bar_table$Model[i]), fixed = TRUE)
  tex_lines <- c(
    tex_lines,
    sprintf("%s & %.2f \\\\", model_label, bar_table$SupportSharePct[i])
  )
}
tex_lines <- c(
  tex_lines,
  "\\hline",
  "\\end{tabular}"
)
writeLines(tex_lines, con = bar_table_tex)

# K-evolution: Monte Carlo volume fraction in the full p_K-dimensional hypercube.
# This is the correct computation: for K=3 the space is R^3, for K=4 it is R^6, etc.
cat("\nRunning Monte Carlo volume fraction experiment across K...\n")
cat(sprintf("  (n_mc = 400000 uniform draws from [0, %.1f]^{p_K} per K)\n", cube_max))
k_trend_df <- compute_mc_volume_fractions_k_grid(k_grid = k_grid, cube_max = cube_max)

k_trend_df$family <- factor(k_trend_df$family,
                             levels = c("WST", "SST", "Toeplitz_SST", "LST"))
k_trend_csv <- file.path(out_dir, "support_volume_fraction_k_trend.csv")
write.csv(k_trend_df, k_trend_csv, row.names = FALSE)

family_colors_trend <- c(
  WST = "#4C78A8", SST = "#F58518", Toeplitz_SST = "#54A24B", LST = "#B279A2"
)

# Panel A: volume fractions for positive-volume families (WST and SST)
k_vol_df <- k_trend_df[k_trend_df$family %in% c("WST", "SST"), ]
k_vol_df$vol_plot <- pmax(k_vol_df$volume_fraction, 1e-5)  # avoid log(0)
g_vol <- ggplot2::ggplot(
  k_vol_df,
  ggplot2::aes(x = K, y = vol_plot, color = family, group = family)
) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::annotate(
    "text", x = max(k_vol_df$K) - 0.1, y = 0.04,
    label = "LST, Toeplitz-SST: zero volume (dim = K-1 << p_K)",
    hjust = 1, size = 3.3, color = "#555555", fontface = "italic"
  ) +
  ggplot2::annotate(
    "text", x = 3, y = 0.5,
    label = "K=3: SST = 1/3 (analytic)", hjust = 0, size = 3.2,
    color = "#F58518", fontface = "italic"
  ) +
  ggplot2::scale_x_continuous(breaks = sort(unique(k_vol_df$K))) +
  ggplot2::scale_y_log10(
    limits = c(1e-5, 1.5),
    breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1.0),
    labels = c("0%", "0.01%", "0.1%", "1%", "10%", "100%")
  ) +
  ggplot2::scale_color_manual(values = family_colors_trend) +
  ggplot2::labs(
    title = "(a) Lebesgue volume fraction in [0, cube_max]^{p_K}  (log scale)",
    subtitle = "WST = 100% (analytic). SST by Monte Carlo (n=400k). LST/Toeplitz-SST: measure zero.",
    x = NULL, y = "Volume fraction (log scale)", color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"))

# Panel B: intrinsic dimension / p_K for all four families
k_dim_df <- k_trend_df
k_dim_df$dim_fraction <- k_dim_df$intrinsic_dim / k_dim_df$pK
g_dim <- ggplot2::ggplot(
  k_dim_df,
  ggplot2::aes(x = K, y = dim_fraction, color = family, group = family,
               linetype = family %in% c("LST", "Toeplitz_SST"))
) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::scale_x_continuous(breaks = sort(unique(k_dim_df$K))) +
  ggplot2::scale_y_continuous(limits = c(0, 1.05),
                               labels = function(x) sprintf("%.0f%%", 100 * x)) +
  ggplot2::scale_color_manual(values = family_colors_trend) +
  ggplot2::scale_linetype_manual(values = c(`TRUE` = "dashed", `FALSE` = "solid"),
                                  guide = "none") +
  ggplot2::labs(
    title = "(b) Intrinsic dimension / p_K",
    subtitle = "WST, SST: full-dim (= 1). LST, Toeplitz-SST: dim = K-1, fraction = 2/K (dashed).",
    x = "Number of blocks K", y = "dim / p_K", color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "none",
                 panel.grid.minor = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = "bold"))

k_trend_png <- file.path(out_dir, "support_volume_fraction_k_trend.png")
if (requireNamespace("patchwork", quietly = TRUE)) {
  ggplot2::ggsave(k_trend_png,
                  patchwork::wrap_plots(g_vol, g_dim, ncol = 1),
                  width = 8.0, height = 10.0, dpi = 190)
} else {
  # fallback: save panel A only and panel B separately
  ggplot2::ggsave(k_trend_png, g_vol, width = 8.0, height = 5.0, dpi = 190)
  ggplot2::ggsave(sub(".png", "_dim.png", k_trend_png, fixed = TRUE),
                  g_dim, width = 8.0, height = 5.0, dpi = 190)
}

cat("\nMonte Carlo volume fraction summary:\n")
print(data.frame(
  K             = k_trend_df$K,
  family        = as.character(k_trend_df$family),
  pK            = k_trend_df$pK,
  intrinsic_dim = k_trend_df$intrinsic_dim,
  vol_fraction  = sprintf("%.4f", k_trend_df$volume_fraction)
), row.names = FALSE)

# Alternative large-K restrictiveness metric for measure-zero families:
# epsilon-neighborhood volume around each support in ambient space.
cat("\nRunning epsilon-neighborhood restrictiveness analysis...\n")
eps_rel <- c(0.02, 0.05, 0.10)
eps_df <- compute_epsilon_neighborhood_k_grid(
  k_grid = k_grid,
  cube_max = cube_max,
  n_mc = 100000,
  eps_rel = eps_rel
)
eps_df$family <- factor(eps_df$family, levels = c("SST", "Toeplitz_SST", "LST"))

eps_csv <- file.path(out_dir, "support_epsilon_neighborhood_k_trend.csv")
write.csv(eps_df, eps_csv, row.names = FALSE)

# Relative restrictiveness versus Toeplitz-SST (same epsilon, same K).
eps_ref <- eps_df[eps_df$family == "Toeplitz_SST", c("K", "epsilon_rel", "neighborhood_share")]
names(eps_ref)[3] <- "toep_share"
eps_rel_df <- merge(eps_df, eps_ref, by = c("K", "epsilon_rel"), all.x = TRUE)
eps_rel_df$rel_to_toep <- eps_rel_df$neighborhood_share / pmax(eps_rel_df$toep_share, 1e-12)
eps_rel_csv <- file.path(out_dir, "support_epsilon_relative_to_toep.csv")
write.csv(eps_rel_df, eps_rel_csv, row.names = FALSE)

eps_palette <- c(SST = "#F58518", Toeplitz_SST = "#54A24B", LST = "#B279A2")

g_eps <- ggplot2::ggplot(
  eps_df,
  ggplot2::aes(x = K, y = neighborhood_share, color = family, group = family)
) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::facet_wrap(~epsilon_rel, ncol = 1,
                      labeller = ggplot2::labeller(epsilon_rel = function(x) sprintf("epsilon = %.2f", as.numeric(x)))) +
  ggplot2::scale_x_continuous(breaks = sort(unique(eps_df$K))) +
  ggplot2::scale_y_log10(
    limits = c(1e-6, 1.05),
    breaks = c(1e-6, 1e-4, 1e-2, 1),
    labels = c("1e-6", "1e-4", "1e-2", "1")
  ) +
  ggplot2::scale_color_manual(values = eps_palette) +
  ggplot2::labs(
    title = "Epsilon-neighborhood volume across K",
    subtitle = "U ~ Uniform([0, cube_max]^{p_K}); y-axis is P(dist(U, support) / cube_max <= epsilon)",
    x = "Number of blocks K",
    y = "Neighborhood share (log scale)",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

g_eps_rel <- ggplot2::ggplot(
  eps_rel_df[eps_rel_df$family %in% c("SST", "LST"), ],
  ggplot2::aes(x = K, y = rel_to_toep, color = family, group = family)
) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#666666", linewidth = 0.5) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::facet_wrap(~epsilon_rel, ncol = 1,
                      labeller = ggplot2::labeller(epsilon_rel = function(x) sprintf("epsilon = %.2f", as.numeric(x)))) +
  ggplot2::scale_x_continuous(breaks = sort(unique(eps_rel_df$K))) +
  ggplot2::scale_y_log10() +
  ggplot2::scale_color_manual(values = eps_palette[c("SST", "LST")]) +
  ggplot2::labs(
    title = "Relative neighborhood volume vs Toeplitz-SST",
    subtitle = "Values > 1 indicate less restrictive support than Toeplitz-SST",
    x = "Number of blocks K",
    y = "share(model) / share(Toeplitz-SST)  (log scale)",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

eps_png <- file.path(out_dir, "support_epsilon_neighborhood_k_trend.png")
eps_rel_png <- file.path(out_dir, "support_epsilon_relative_to_toep.png")
ggplot2::ggsave(eps_png, g_eps, width = 8.0, height = 10.0, dpi = 190)
ggplot2::ggsave(eps_rel_png, g_eps_rel, width = 8.0, height = 10.0, dpi = 190)

# Unified K=3 visuals using the same epsilon-neighborhood method.
cat("\nBuilding K=3 visuals with unified epsilon-neighborhood method...\n")
unified_k3 <- compute_unified_k3_visuals(
  cube_max = cube_max,
  breaks = breaks,
  eps_vis = 0.05,
  n_mc_vis = 400000
)
avg_cov_list_unified <- unified_k3$avg_cov_list

summary_rows_unified <- list()
for (fam in families) {
  stats_vec <- summarise_restrictiveness(avg_cov_list_unified[[fam]], occ_threshold)
  summary_rows_unified[[fam]] <- data.frame(
    family = fam,
    mean_coverage = unname(stats_vec["mean_coverage"]),
    q95_coverage = unname(stats_vec["q95_coverage"]),
    support_share = unname(stats_vec["support_share"]),
    eps_vis = unified_k3$eps_vis,
    n_mc_vis = unified_k3$n_mc_vis
  )
}
summary_unified_df <- do.call(rbind, summary_rows_unified)
summary_unified_csv <- file.path(out_dir, "support_unified_k3_summary.csv")
write.csv(summary_unified_df, summary_unified_csv, row.names = FALSE)

slice_df_list_unified <- list()
for (fam in families) {
  tmp <- build_slice_df(
    centers = centers,
    avg_cov = avg_cov_list_unified[[fam]],
    n_bins = n_bins,
    cube_max = cube_max
  )
  tmp$family <- fam
  slice_df_list_unified[[fam]] <- tmp
}
slice_df_u <- do.call(rbind, slice_df_list_unified)
slice_df_u$family <- factor(slice_df_u$family, levels = families)
slice_df_u$family_label <- factor(
  as.character(slice_df_u$family),
  levels = c("WST", "SST", "Toeplitz_SST", "LST"),
  labels = c("WST", "SST", "Toeplitz SST", "LST")
)

panel_share_u <- data.frame(
  family = as.character(summary_unified_df$family),
  support_share = summary_unified_df$support_share,
  stringsAsFactors = FALSE
)
slice_df_u <- merge(slice_df_u, panel_share_u, by = "family", all.x = TRUE)
slice_df_u$panel_note <- sprintf("Covered: %.1f%%", 100 * slice_df_u$support_share)

# Overwrite manuscript figure paths with unified-method visuals.
slice_png <- file.path(out_dir, "support_cross_sections_static.png")
g_slice_u <- ggplot2::ggplot(slice_df_u, ggplot2::aes(x = psi12, y = psi23, fill = avg_cov)) +
  ggplot2::geom_raster(interpolate = FALSE) +
  ggplot2::geom_contour(
    ggplot2::aes(z = avg_cov),
    breaks = c(0.15, 0.35, 0.55, 0.75),
    color = "white",
    linewidth = 0.24,
    alpha = 0.78
  ) +
  ggplot2::geom_text(
    data = unique(slice_df_u[, c("family_label", "slice_label", "panel_note")]),
    mapping = ggplot2::aes(x = 0.10 * cube_max, y = 0.93 * cube_max, label = panel_note),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    color = "white",
    fontface = "bold"
  ) +
  ggplot2::facet_grid(slice_label ~ family_label) +
  ggplot2::scale_fill_gradientn(
    colours = c("#fffdf5", "#fee8a8", "#fdd072", "#fca85e", "#f67e4b", "#e34a33"),
    limits = c(0, 1),
    name = "Neighborhood\nshare"
  ) +
  ggplot2::scale_x_continuous(breaks = 0:4, minor_breaks = NULL, expand = c(0, 0)) +
  ggplot2::scale_y_continuous(breaks = 0:4, minor_breaks = NULL, expand = c(0, 0)) +
  ggplot2::coord_equal(xlim = c(0, cube_max), ylim = c(0, cube_max), expand = FALSE) +
  ggplot2::labs(
    title = "Unified Geometry Cross-Sections In Psi-Space",
    subtitle = sprintf("All panels use P(dist(psi, support) <= %.2f * cube_max) from uniform ambient draws", unified_k3$eps_vis),
    x = "psi12",
    y = "psi23"
  ) +
  ggplot2::theme_minimal(base_size = 11.5) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    strip.background = ggplot2::element_rect(fill = "#f2f2f2", color = "#d9d9d9", linewidth = 0.3),
    panel.border = ggplot2::element_rect(color = "#e2e2e2", fill = NA, linewidth = 0.35),
    axis.text = ggplot2::element_text(color = "#3a3a3a"),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(color = "#4a4a4a"),
    panel.spacing = grid::unit(0.4, "lines")
  )

ggplot2::ggsave(slice_png, g_slice_u, width = 12, height = 8.5, dpi = 180)

overlay_df_u <- do.call(rbind, lapply(families, function(fam) {
  build_support_points_for_overlay(centers, avg_cov_list_unified[[fam]], fam, occ_threshold)
}))
overlay_small_png <- file.path(out_dir, "support_3d_small_multiples_static.png")
grDevices::png(overlay_small_png, width = 1900, height = 1500, res = 180)
share_map_u <- stats::setNames(summary_unified_df$support_share, as.character(summary_unified_df$family))
save_static_3d_small_multiples(overlay_df_u, cube_max, support_share_map = share_map_u)
dev.off()

cat("\nEpsilon-neighborhood summary (selected rows):\n")
print(data.frame(
  K = eps_df$K,
  family = as.character(eps_df$family),
  epsilon_rel = sprintf("%.2f", eps_df$epsilon_rel),
  neighborhood_share = sprintf("%.6f", eps_df$neighborhood_share),
  mean_dist = sprintf("%.4f", eps_df$mean_dist)
), row.names = FALSE)

cat("\nDone. Files written:\n")
cat(sprintf(" - %s\n", summary_csv))
cat(sprintf(" - %s\n", bar_png))
cat(sprintf(" - %s\n", html_file))
cat(sprintf(" - %s\n", rds_file))
cat(sprintf(" - %s\n", slice_png))
cat(sprintf(" - %s\n", overlay_png))
cat(sprintf(" - %s\n", overlay_small_png))
cat(sprintf(" - %s\n", proj_png))
cat(sprintf(" - %s\n", bar_table_csv))
cat(sprintf(" - %s\n", bar_table_tex))
cat(sprintf(" - %s\n", tendency_csv))
cat(sprintf(" - %s\n", k_trend_csv))
cat(sprintf(" - %s\n", k_trend_png))
cat(sprintf(" - %s\n", eps_csv))
cat(sprintf(" - %s\n", eps_rel_csv))
cat(sprintf(" - %s\n", eps_png))
cat(sprintf(" - %s\n", eps_rel_png))
cat(sprintf(" - %s\n", summary_unified_csv))
