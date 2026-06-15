# ======================================================================
# Posterior Predictive Checks for OSBM (WST / SST) + DCSBM
# ======================================================================
#
# Key exported function:
#   .run_ppc_model(A_obs, out_relab, regime, ...)
#
# Generates:
#   - Network-level summary statistics (obs vs. replicated)
#   - Bayesian p-values
#   - bayesplot:: ppc_dens_overlay, ppc_stat, ppc_bars
#   - Multi-page PDF with base-R histograms
#   - CSV / RDS outputs
# ======================================================================

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

.gini <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  if (all(x == 0))     return(0)
  x <- sort(x)
  n <- length(x)
  (2 * sum(x * seq_len(n)) - (n + 1) * sum(x)) / (n * sum(x))
}

# ------------------------------------------------------------------
# Extract a single MCMC draw from sampler output
# Handles both list-of-draws and matrix/array storage.
# ------------------------------------------------------------------
.extract_draw <- function(out, s) {
  get_s <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) return(x[[s]])
    if (is.matrix(x)) return(x[s, , drop = TRUE])
    if (!is.null(dim(x)) && length(dim(x)) == 3L)
      return(x[s, , , drop = TRUE])
    stop("Cannot extract draw from object of class: ",
         paste(class(x), collapse = ","))
  }
  list(
    z      = get_s(out$z),
    eta    = get_s(out$eta),
    theta  = get_s(out$theta),
    kappa  = get_s(out$kappa),
    psi    = get_s(out$psi),
    lambda = get_s(out$lambda)
  )
}

# ------------------------------------------------------------------
# Convert a WST psi draw to an upper-triangular K x K matrix
# ------------------------------------------------------------------
.psi_to_upper_matrix_wst <- function(psi_draw, K) {
  M <- matrix(0, K, K)
  if (is.matrix(psi_draw)) {
    nr <- min(nrow(psi_draw), K)
    nc <- min(ncol(psi_draw), K)
    sub <- psi_draw[seq_len(nr), seq_len(nc), drop = FALSE]
    if (nr == K && nc == K) {
      M[upper.tri(M)] <- sub[upper.tri(sub)]
    } else {
      # Pad: use what we have
      M[seq_len(nr), seq_len(nc)] <- sub
      M[lower.tri(M)] <- 0
      diag(M) <- 0
    }
    return(M)
  }
  psi_draw <- as.numeric(psi_draw)
  if (length(psi_draw) >= choose(K, 2)) {
    M[upper.tri(M)] <- psi_draw[seq_len(choose(K, 2))]
    return(M)
  }
  stop("Unrecognised WST psi representation (need KxK matrix or ",
       "length >= choose(K,2) vector).")
}

# ------------------------------------------------------------------
# Simulate one replicated network from a single posterior draw
# ------------------------------------------------------------------
.simulate_A_from_draw <- function(draw, regime) {
  stopifnot(regime %in% c("WST", "SST", "DCSBM"))

  z <- as.integer(draw$z)
  n <- length(z)
  labs <- sort(unique(z))
  K <- length(labs)

  # --- degree-correction parameters ---
  if (regime %in% c("WST", "SST")) {
    if (is.null(draw$eta) || length(draw$eta) < n)
      stop("eta draw missing or wrong length.")
    eta <- as.numeric(draw$eta)
  } else {
    if (is.null(draw$theta) || length(draw$theta) < n)
      stop("theta draw missing or wrong length.")
    theta <- as.numeric(draw$theta)
  }

  # Remap labels to 1..K (handles gaps from variable-K)
  z_remap <- match(z, labs)

  # --- block-level parameters ---
  if (regime %in% c("WST", "SST")) {
    kappa <- draw$kappa
    if (!is.matrix(kappa))
      stop("kappa draw is not a matrix.")
    # Subset to active labels if raw matrix is larger
    if (nrow(kappa) >= max(labs) && ncol(kappa) >= max(labs)) {
      kappa <- kappa[labs, labs, drop = FALSE]
    } else if (nrow(kappa) >= K && ncol(kappa) >= K) {
      kappa <- kappa[seq_len(K), seq_len(K), drop = FALSE]
    } else {
      stop("kappa matrix too small for active labels (K=", K, ").")
    }
    kappa <- 0.5 * (kappa + t(kappa))

    if (regime == "WST") {
      psiU <- .psi_to_upper_matrix_wst(draw$psi, K)
    } else {
      psi_vec <- as.numeric(draw$psi)
      if (length(psi_vec) < (K - 1L))
        stop("SST psi vector too short for K=", K, ".")
      psi_vec <- psi_vec[seq_len(K - 1L)]
    }

  } else {
    lambda <- draw$lambda
    if (!is.matrix(lambda))
      stop("lambda draw is not a matrix for DCSBM.")
    if (nrow(lambda) >= max(labs) && ncol(lambda) >= max(labs)) {
      lambda <- lambda[labs, labs, drop = FALSE]
    } else if (nrow(lambda) >= K && ncol(lambda) >= K) {
      lambda <- lambda[seq_len(K), seq_len(K), drop = FALSE]
    } else {
      stop("lambda matrix too small for active labels (K=", K, ").")
    }
  }

  # --- simulate edges ---
  Arep <- matrix(0L, n, n)

  for (i in 1:(n - 1L)) {
    for (j in (i + 1L):n) {
      ki <- z_remap[i]; kj <- z_remap[j]

      if (regime %in% c("WST", "SST")) {
        lam_tot <- eta[i] * eta[j] * kappa[ki, kj]
        Nij <- rpois(1L, lam_tot)
        if (Nij == 0L) next

        if (ki == kj) {
          p <- 0.5
        } else {
          sgn <- sign(kj - ki)
          if (regime == "WST") {
            psi_val <- psiU[min(ki, kj), max(ki, kj)]
          } else {
            d <- abs(ki - kj)
            psi_val <- psi_vec[d]
          }
          p <- plogis(sgn * psi_val)
        }
        Aij <- rbinom(1L, size = Nij, prob = p)
        Aji <- Nij - Aij

      } else {
        # DCSBM: two independent Poisson
        lam_ij <- theta[i] * theta[j] * lambda[ki, kj]
        lam_ji <- theta[j] * theta[i] * lambda[kj, ki]
        Aij <- rpois(1L, lam_ij)
        Aji <- rpois(1L, lam_ji)
      }

      Arep[i, j] <- Aij
      Arep[j, i] <- Aji
    }
  }

  diag(Arep) <- 0L
  Arep
}

# ------------------------------------------------------------------
# Triad cycle rate (random sample of triples)
# ------------------------------------------------------------------
.triad_cycle_rate <- function(A, m_triples = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(A)
  if (n < 3L) return(NA_real_)

  cycles <- 0L; used <- 0L

  for (t in seq_len(m_triples)) {
    idx <- sample.int(n, 3L, replace = FALSE)
    i <- idx[1]; j <- idx[2]; k <- idx[3]

    s_ij <- sign(A[i, j] - A[j, i])
    s_jk <- sign(A[j, k] - A[k, j])
    s_ki <- sign(A[k, i] - A[i, k])

    if (s_ij == 0L || s_jk == 0L || s_ki == 0L) next
    used <- used + 1L

    if ((s_ij ==  1L && s_jk ==  1L && s_ki ==  1L) ||
        (s_ij == -1L && s_jk == -1L && s_ki == -1L)) {
      cycles <- cycles + 1L
    }
  }

  if (used == 0L) return(NA_real_)
  cycles / used
}

# ------------------------------------------------------------------
# Compute a vector of network summary statistics
# ------------------------------------------------------------------
.ppc_stats <- function(A, seed_triples = 1L, m_triples = 2000L) {
  n <- nrow(A)

  out_strength <- rowSums(A)
  in_strength  <- colSums(A)

  N <- A + t(A)
  Nij <- N[upper.tri(N)]

  D <- abs(A - t(A))
  Dij <- D[upper.tri(D)]
  asym <- ifelse(Nij > 0, Dij / Nij, NA_real_)

  recip <- if (sum(A) > 0) sum(pmin(A, t(A))) / sum(A) else NA_real_

  # Score-based order conformity
  score    <- out_strength - in_strength
  ord      <- order(score, decreasing = TRUE)
  rank_vec <- integer(n); rank_vec[ord] <- seq_len(n)

  i_idx <- row(N)[upper.tri(N)]
  j_idx <- col(N)[upper.tri(N)]
  valid <- Nij > 0 & (A[cbind(i_idx, j_idx)] != A[cbind(j_idx, i_idx)])
  if (any(valid)) {
    winner_is_i <- A[cbind(i_idx, j_idx)] > A[cbind(j_idx, i_idx)]
    higher_is_i <- rank_vec[i_idx] < rank_vec[j_idx]
    conform     <- mean(winner_is_i[valid] == higher_is_i[valid])
  } else {
    conform <- NA_real_
  }

  cycle_rate <- .triad_cycle_rate(A, m_triples = m_triples,
                                  seed = seed_triples)

  c(
    total_mass         = sum(A),
    edge_nz_rate       = mean(A[row(A) != col(A)] > 0),
    dyad_zero_rate     = mean(Nij == 0),
    dyad_mean          = mean(Nij),
    dyad_var           = stats::var(Nij),
    dyad_max           = max(Nij),
    out_gini           = .gini(out_strength),
    in_gini            = .gini(in_strength),
    reciprocity        = recip,
    asym_mean          = mean(asym, na.rm = TRUE),
    asym_q95           = stats::quantile(asym, 0.95, na.rm = TRUE,
                                         names = FALSE),
    cycle_rate         = cycle_rate,
    node_order_conform = conform
  )
}

# ======================================================================
# Main PPC driver
# ======================================================================

.run_ppc_model <- function(A_obs, out_relab, regime,
                           n_draws   = 200L,
                           m_triples = 2000L,
                           seed      = 1L,
                           out_dir   = ".",
                           tag       = "") {
  stopifnot(regime %in% c("WST", "SST", "DCSBM"))

  # --- Guard: ensure bayesplot + ggplot2 are available ----
  if (!requireNamespace("bayesplot", quietly = TRUE))
    stop("Install bayesplot:  install.packages('bayesplot')")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Install ggplot2:  install.packages('ggplot2')")

  set.seed(seed)
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Build z-chain matrix (S x n) to index draws ----
  z_chain <- if (is.matrix(out_relab$z)) {
    out_relab$z
  } else {
    do.call(rbind, out_relab$z)
  }
  S <- nrow(z_chain)
  draw_ids <- sort(sample.int(S, size = min(n_draws, S), replace = FALSE))

  cat(sprintf("[PPC] regime=%s  S=%d  n_draws=%d  tag=%s\n",
              regime, S, length(draw_ids), tag))

  # --- Observed statistics ----
  obs_stats <- .ppc_stats(A_obs, seed_triples = seed + 11L,
                          m_triples = m_triples)

  # --- Simulate replicated networks & collect statistics ----
  sim_stats <- matrix(NA_real_, nrow = length(draw_ids),
                      ncol = length(obs_stats))
  colnames(sim_stats) <- names(obs_stats)

  # Dyad-total samples for bayesplot overlay
  N_obs     <- A_obs + t(A_obs)
  Nij_obs   <- N_obs[upper.tri(N_obs)]
  m_dyads   <- min(1500L, length(Nij_obs))
  dyad_idx  <- sample.int(length(Nij_obs), size = m_dyads, replace = FALSE)
  y_dyad    <- Nij_obs[dyad_idx]
  yrep_dyad <- matrix(NA_real_, nrow = length(draw_ids), ncol = m_dyads)

  # Edge-presence indicator (0/1) for ppc_bars
  edge_obs  <- as.integer(Nij_obs[dyad_idx] > 0)
  yrep_edge <- matrix(NA_integer_, nrow = length(draw_ids), ncol = m_dyads)

  n_ok <- 0L; n_fail <- 0L

  for (r in seq_along(draw_ids)) {
    s  <- draw_ids[r]
    dr <- .extract_draw(out_relab, s)

    A_rep <- tryCatch(
      .simulate_A_from_draw(dr, regime = regime),
      error = function(e) {
        if (n_fail < 5L)
          message(sprintf("[PPC] draw %d failed: %s", s, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(A_rep)) { n_fail <- n_fail + 1L; next }
    n_ok <- n_ok + 1L

    sim_stats[r, ] <- .ppc_stats(A_rep,
                                 seed_triples = seed + 1000L + r,
                                 m_triples    = m_triples)
    N_rep          <- A_rep + t(A_rep)
    Nij_rep        <- N_rep[upper.tri(N_rep)]
    yrep_dyad[r, ] <- Nij_rep[dyad_idx]
    yrep_edge[r, ] <- as.integer(Nij_rep[dyad_idx] > 0)
  }

  cat(sprintf("[PPC] %d/%d draws simulated successfully (%d failed)\n",
              n_ok, length(draw_ids), n_fail))

  if (n_ok == 0L)
    stop("[PPC] All draws failed for regime=", regime)

  # Drop rows where simulation failed
  good      <- !is.na(sim_stats[, 1])
  sim_stats <- sim_stats[good, , drop = FALSE]
  yrep_dyad <- yrep_dyad[good, , drop = FALSE]
  yrep_edge <- yrep_edge[good, , drop = FALSE]
  sim_df    <- as.data.frame(sim_stats)
  sim_df$draw_id <- draw_ids[good]

  # --- Bayesian p-values (upper tail) ----
  pvals <- vapply(names(obs_stats), function(nm)
    mean(sim_df[[nm]] >= obs_stats[[nm]], na.rm = TRUE), numeric(1))
  pval_df <- tibble::tibble(
    statistic = names(obs_stats),
    observed  = as.numeric(obs_stats),
    p_upper   = as.numeric(pvals)
  )

  # ==================================================================
  # PLOTS
  # ==================================================================

  # ---- 1. Multi-page PDF: base-R histograms ----
  pdf_file <- file.path(out_dir,
                        paste0("ppc_histograms_", regime, tag, ".pdf"))
  grDevices::pdf(pdf_file, width = 8, height = 6)

  key_stats <- c("dyad_zero_rate", "dyad_mean", "dyad_max",
                 "reciprocity", "asym_mean", "cycle_rate",
                 "node_order_conform")
  key_stats <- intersect(key_stats, names(obs_stats))

  for (nm in key_stats) {
    vals <- sim_df[[nm]]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 2L) next
    hist(vals, main = paste0("PPC: ", nm, "  (", regime, ")"),
         xlab = nm, breaks = 30, col = "grey85", border = "grey60")
    abline(v = obs_stats[[nm]], col = "red", lwd = 2, lty = 2)
    legend("topright",
           legend = sprintf("obs = %.4f  (p = %.3f)",
                            obs_stats[[nm]], pvals[[nm]]),
           bty = "n", cex = 0.9)
  }
  grDevices::dev.off()
  cat(sprintf("[PPC] Saved histogram PDF: %s\n", pdf_file))

  # ---- 2. bayesplot: density overlay (dyad totals) ----
  max_overlay <- min(50L, nrow(yrep_dyad))
  yrep_plot   <- yrep_dyad[sample.int(nrow(yrep_dyad), max_overlay), ,
                           drop = FALSE]

  p1 <- bayesplot::ppc_dens_overlay(y = y_dyad, yrep = yrep_plot) +
    ggplot2::ggtitle(paste0("PPC density overlay (dyad totals) — ", regime))

  png1 <- file.path(out_dir,
                    paste0("ppc_dens_overlay_", regime, tag, ".png"))
  ggplot2::ggsave(filename = png1, plot = p1,
                  width = 8, height = 5, dpi = 150)
  cat(sprintf("[PPC] Saved bayesplot density overlay: %s\n", png1))

  # ---- 3. bayesplot: stat check — mean of dyad totals ----
  p2 <- bayesplot::ppc_stat(y = y_dyad, yrep = yrep_dyad, stat = "mean") +
    ggplot2::ggtitle(paste0("PPC T(y) = mean  (dyad totals) — ", regime))

  png2 <- file.path(out_dir,
                    paste0("ppc_stat_mean_", regime, tag, ".png"))
  ggplot2::ggsave(filename = png2, plot = p2,
                  width = 8, height = 5, dpi = 150)
  cat(sprintf("[PPC] Saved bayesplot stat-mean: %s\n", png2))

  # ---- 4. bayesplot: stat check — sd of dyad totals ----
  p3 <- bayesplot::ppc_stat(y = y_dyad, yrep = yrep_dyad, stat = "sd") +
    ggplot2::ggtitle(paste0("PPC T(y) = sd  (dyad totals) — ", regime))

  png3 <- file.path(out_dir,
                    paste0("ppc_stat_sd_", regime, tag, ".png"))
  ggplot2::ggsave(filename = png3, plot = p3,
                  width = 8, height = 5, dpi = 150)
  cat(sprintf("[PPC] Saved bayesplot stat-sd: %s\n", png3))

  # ---- 5. bayesplot: bars (edge presence 0/1) ----
  p4 <- bayesplot::ppc_bars(y = edge_obs, yrep = yrep_edge) +
    ggplot2::ggtitle(paste0("PPC edge presence — ", regime))

  png4 <- file.path(out_dir,
                    paste0("ppc_bars_edge_", regime, tag, ".png"))
  ggplot2::ggsave(filename = png4, plot = p4,
                  width = 7, height = 5, dpi = 150)
  cat(sprintf("[PPC] Saved bayesplot bars: %s\n", png4))

  # ==================================================================
  # Return
  # ==================================================================
  list(
    obs_stats  = obs_stats,
    draw_stats = sim_df,
    pvals      = pval_df,
    pdf        = pdf_file,
    bayesplot  = c(density_overlay = png1, stat_mean = png2,
                   stat_sd = png3, bars_edge = png4)
  )
}
