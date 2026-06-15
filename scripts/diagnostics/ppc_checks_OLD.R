# ======================================================================
# Posterior predictive checks for OSBM (WST/SST) + DCSBM
# ======================================================================

.gini <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (all(x == 0)) return(0)
  x <- sort(x)
  n <- length(x)
  (2 * sum(x * seq_len(n)) - (n + 1) * sum(x)) / (n * sum(x))
}

.extract_draw <- function(out, s) {
  get_s <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.list(x)) return(x[[s]])
    if (is.matrix(x)) return(x[s, , drop = TRUE])
    if (!is.null(dim(x)) && length(dim(x)) == 3) return(x[s, , , drop = TRUE])
    stop("Don't know how to extract draw s from object of class: ", paste(class(x), collapse = ","))
  }
  list(
    z     = get_s(out$z),
    eta   = get_s(out$eta),
    theta = get_s(out$theta),
    kappa = get_s(out$kappa),
    psi   = get_s(out$psi),
    lambda = get_s(out$lambda)
  )
}

.psi_to_upper_matrix_wst <- function(psi_draw, K) {
  # Goal: return KxK matrix with upper-triangular entries psi_{k<l}, diag 0
  # Accepts: KxK matrix, length choose(K,2) vector, or list element already formatted.
  M <- matrix(0, K, K)
  if (is.matrix(psi_draw)) {
    if (!all(dim(psi_draw) == c(K, K))) stop("WST psi matrix has wrong dims.")
    # If full matrix, assume upper triangle contains positive parameters
    M[upper.tri(M)] <- psi_draw[upper.tri(psi_draw)]
    return(M)
  }
  psi_draw <- as.numeric(psi_draw)
  if (length(psi_draw) == choose(K, 2)) {
    M[upper.tri(M)] <- psi_draw
    return(M)
  }
  stop("Unrecognised WST psi representation (need KxK or length choose(K,2)).")
}

.simulate_A_from_draw <- function(draw, regime) {
  stopifnot(regime %in% c("WST","SST", "DCSBM"))
  z <- as.integer(draw$z)
  n <- length(z)
  labs <- sort(unique(z))
  K <- length(labs)

  if (regime %in% c("WST", "SST")) {
    if (is.null(draw$eta) || length(draw$eta) != n) stop("eta draw missing or wrong length.")
    eta <- as.numeric(draw$eta)
  } else {
    if (is.null(draw$theta) || length(draw$theta) != n) stop("theta draw missing or wrong length.")
    theta <- as.numeric(draw$theta)
  }
  
  # Recode labels to 1..K (defensive, avoids gaps)
  z <- match(z, labs)
  
  if (regime %in% c("WST", "SST")) {
    # kappa: want KxK symmetric
    kappa <- draw$kappa
    if (is.matrix(kappa)) {
      if (nrow(kappa) < max(labs) || ncol(kappa) < max(labs)) {
        stop("kappa matrix dims are too small for active labels in this draw.")
      }
      if (!all(dim(kappa) == c(K, K))) kappa <- kappa[labs, labs, drop = FALSE]
    } else {
      stop("kappa draw is not a matrix; adapt extractor to your output structure.")
    }
    kappa <- 0.5 * (kappa + t(kappa))  # enforce symmetry

    # psi
    if (regime == "WST") {
      psiU <- .psi_to_upper_matrix_wst(draw$psi, K)
    } else {
      psi_vec <- as.numeric(draw$psi)
      if (length(psi_vec) < (K - 1)) stop("SST psi vector too short for this draw's K.")
      psi_vec <- psi_vec[seq_len(K - 1)]
    }
  } else {
    lambda <- draw$lambda
    if (!is.matrix(lambda)) stop("lambda draw is not a matrix for DCSBM draw.")
    if (nrow(lambda) < max(labs) || ncol(lambda) < max(labs)) {
      stop("lambda matrix dims are too small for active labels in this draw.")
    }
    if (!all(dim(lambda) == c(K, K))) lambda <- lambda[labs, labs, drop = FALSE]
  }
  
  Arep <- matrix(0L, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      ki <- z[i]; kj <- z[j]
      if (regime %in% c("WST", "SST")) {
        lam_tot <- eta[i] * eta[j] * kappa[ki, kj]
        Nij <- rpois(1L, lam_tot)

        if (Nij == 0L) next

        if (ki == kj) {
          p <- 0.5
        } else {
          sgn <- sign(kj - ki)  # +1 if i is "stronger block" (smaller index) -> forward i->j
          if (regime == "WST") {
            psi <- psiU[min(ki,kj), max(ki,kj)]
          } else {
            d <- abs(ki - kj)
            psi <- psi_vec[d]
          }
          p <- plogis(sgn * psi)
        }

        Aij <- rbinom(1L, size = Nij, prob = p)
        Aji <- Nij - Aij
      } else {
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

.triad_cycle_rate <- function(A, m_triples = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(A)
  if (n < 3) return(NA_real_)
  
  pick3 <- function() sample.int(n, 3L, replace = FALSE)
  
  cycles <- 0L
  used <- 0L
  
  for (t in seq_len(m_triples)) {
    idx <- pick3()
    i <- idx[1]; j <- idx[2]; k <- idx[3]
    
    # Pairwise direction by sign(Aij - Aji); drop ties and zero-total pairs
    s_ij <- sign(A[i,j] - A[j,i])
    s_jk <- sign(A[j,k] - A[k,j])
    s_ki <- sign(A[k,i] - A[i,k])
    
    if (s_ij == 0L || s_jk == 0L || s_ki == 0L) next  # ignore ties / no signal
    
    used <- used + 1L
    # A 3-cycle exists if directions are consistent around the loop
    # (i beats j, j beats k, k beats i) OR the reverse.
    if ((s_ij ==  1L && s_jk ==  1L && s_ki ==  1L) ||
        (s_ij == -1L && s_jk == -1L && s_ki == -1L)) {
      cycles <- cycles + 1L
    }
  }
  
  if (used == 0L) return(NA_real_)
  cycles / used
}

.ppc_stats <- function(A, seed_triples = 1L, m_triples = 2000L) {
  n <- nrow(A)
  off <- which(row(A) != col(A), arr.ind = TRUE)
  
  out_strength <- rowSums(A)
  in_strength  <- colSums(A)
  
  # dyad totals
  N <- A + t(A)
  Nij <- N[upper.tri(N)]
  
  # asymmetry (only where Nij>0)
  D <- abs(A - t(A))
  Dij <- D[upper.tri(D)]
  asym <- ifelse(Nij > 0, Dij / Nij, NA_real_)
  
  # weighted reciprocity: sum(min(Aij,Aji)) / sum(Aij)
  recip <- if (sum(A) > 0) sum(pmin(A, t(A))) / sum(A) else NA_real_
  
  # simple score-based order conformity at node level: score = out - in
  score <- out_strength - in_strength
  ord <- order(score, decreasing = TRUE)
  rank <- integer(n); rank[ord] <- seq_len(n)
  # for dyads with Nij>0, check whether higher-rank tends to "win"
  i_idx <- row(N)[upper.tri(N)]
  j_idx <- col(N)[upper.tri(N)]
  valid <- Nij > 0 & (A[cbind(i_idx, j_idx)] != A[cbind(j_idx, i_idx)])
  if (any(valid)) {
    winner_is_i <- A[cbind(i_idx, j_idx)] > A[cbind(j_idx, i_idx)]
    higher_is_i <- rank[i_idx] < rank[j_idx]
    conform <- mean(winner_is_i[valid] == higher_is_i[valid])
  } else {
    conform <- NA_real_
  }
  
  cycle_rate <- .triad_cycle_rate(A, m_triples = m_triples, seed = seed_triples)
  
  c(
    total_mass = sum(A),
    edge_nz_rate = mean(A[off] > 0),
    dyad_zero_rate = mean(Nij == 0),
    dyad_mean = mean(Nij),
    dyad_var = stats::var(Nij),
    dyad_max = max(Nij),
    out_gini = .gini(out_strength),
    in_gini  = .gini(in_strength),
    reciprocity = recip,
    asym_mean = mean(asym, na.rm = TRUE),
    asym_q95  = stats::quantile(asym, 0.95, na.rm = TRUE, names = FALSE),
    cycle_rate = cycle_rate,
    node_order_conform = conform
  )
}

.run_ppc_model <- function(A_obs, out_relab, regime,
                           n_draws = 200L, m_triples = 2000L,
                           seed = 1L, out_dir = ".", tag = "") {
  stopifnot(regime %in% c("WST", "SST", "DCSBM"))
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    stop("Package 'bayesplot' is required for PPC plotting. Install with install.packages('bayesplot').")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for PPC plotting. Install with install.packages('ggplot2').")
  }

  set.seed(seed)

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Determine number of saved draws
  z_chain <- if (is.matrix(out_relab$z)) out_relab$z else do.call(rbind, out_relab$z)
  S <- nrow(z_chain)

  draw_ids <- sample.int(S, size = min(n_draws, S), replace = FALSE)

  obs_stats <- .ppc_stats(A_obs, seed_triples = seed + 11L, m_triples = m_triples)

  sim_stats <- matrix(NA_real_, nrow = length(draw_ids), ncol = length(obs_stats))
  colnames(sim_stats) <- names(obs_stats)

  N_obs <- A_obs + t(A_obs)
  Nij_obs <- N_obs[upper.tri(N_obs)]
  m_dyads <- min(1500L, length(Nij_obs))
  dyad_idx <- sample.int(length(Nij_obs), size = m_dyads, replace = FALSE)
  y_dyad <- Nij_obs[dyad_idx]
  yrep_dyad <- matrix(NA_real_, nrow = length(draw_ids), ncol = m_dyads)

  for (r in seq_along(draw_ids)) {
    s <- draw_ids[r]
    dr <- .extract_draw(out_relab, s)
    A_rep <- .simulate_A_from_draw(dr, regime = regime)
    sim_stats[r, ] <- .ppc_stats(A_rep, seed_triples = seed + 1000L + r, m_triples = m_triples)

    N_rep <- A_rep + t(A_rep)
    Nij_rep <- N_rep[upper.tri(N_rep)]
    yrep_dyad[r, ] <- Nij_rep[dyad_idx]
  }

  sim_df <- as.data.frame(sim_stats)
  sim_df$draw_id <- draw_ids

  # Bayesian p-values (two-sided-ish is your choice; here: upper-tail)
  pvals <- vapply(names(obs_stats), function(nm) mean(sim_df[[nm]] >= obs_stats[[nm]], na.rm = TRUE), numeric(1))
  pval_df <- tibble::tibble(
    statistic = names(obs_stats),
    observed  = as.numeric(obs_stats),
    p_upper   = as.numeric(pvals)
  )
  
  # Plot a few key stats (pdf)
  pdf_file <- file.path(out_dir, paste0("ppc_", regime, tag, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"))
  grDevices::pdf(pdf_file, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)

  key_stats <- c("dyad_zero_rate","dyad_mean","dyad_max","reciprocity","asym_mean","cycle_rate","node_order_conform")
  key_stats <- intersect(key_stats, names(obs_stats))

  for (nm in key_stats) {
    hist(sim_df[[nm]], main = paste0("PPC: ", nm, " (", regime, ")"),
         xlab = nm, breaks = 30)
    abline(v = obs_stats[[nm]], lwd = 2)
  }

  yrep_plot <- yrep_dyad
  max_overlay <- min(50L, nrow(yrep_plot))
  if (max_overlay < nrow(yrep_plot)) {
    keep <- sample.int(nrow(yrep_plot), size = max_overlay, replace = FALSE)
    yrep_plot <- yrep_plot[keep, , drop = FALSE]
  }

  p1 <- bayesplot::ppc_dens_overlay(y = y_dyad, yrep = yrep_plot) +
    ggplot2::ggtitle(paste0("PPC density overlay (dyad totals) - ", regime))

  p2 <- bayesplot::ppc_stat(y = y_dyad, yrep = yrep_dyad, stat = "mean") +
    ggplot2::ggtitle(paste0("PPC mean statistic (dyad totals) - ", regime))

  png1 <- file.path(out_dir, paste0("ppc_", regime, tag, "_dens_overlay.png"))
  png2 <- file.path(out_dir, paste0("ppc_", regime, tag, "_stat_mean.png"))
  ggplot2::ggsave(filename = png1, plot = p1, width = 8, height = 5, dpi = 140)
  ggplot2::ggsave(filename = png2, plot = p2, width = 8, height = 5, dpi = 140)

  list(
    obs_stats = obs_stats,
    draw_stats = sim_df,
    pvals = pval_df,
    pdf = pdf_file,
    bayesplot = c(density_overlay = png1, stat_mean = png2)
  )
}