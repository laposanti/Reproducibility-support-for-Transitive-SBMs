# Variable-K ordered SBM sampler.
#
# This file contains the main Gibbs/MH loop for the ordered SBM models.
# Sampler-local utilities now live under helper_folder/models/ordered_sbm/.

Rcpp::sourceCpp("core/block_totals_for_poisson_cpp.cpp")
source("helper_folder/load_sampler_helpers.R", chdir = FALSE)

modular_osbm_sampler <- function(
    A, K, truth = NA,
    free   = c("psi","kappa","eta","z"),
    n_iter = 4000, burn = 500, thin = 1,
    verbose = FALSE,
    psi_constraint = c("WST","SST"),
    seed = NULL,
    eta_identifiability = c("none", "block_sum"),
    # NEW: controls
    shrink_when = c("after_z_sweep","end_of_iter","never"),
    refresh_pg_after_birth = TRUE,
    a_kappa = 1, b_kappa = 1,
    a_eta = 1, b_eta = 1,
    mu0 = 1, sigma0 = 2, tau0 = 0.15,
    alpha0 = 1.0, discount = 0.0,
    partition_prior = 'OCRP',
    gamma_gn = 0.8,
    alpha_gp05 = 0.5,
    theta_ocrp = 1.0,
    sst_birth_score_mode = c("exact_nonlocal", "local_approx"),
    slot_radius_after_burnin = NULL,
    sample_b_kappa = FALSE,
    alpha0_bkappa  = 1.0,    # hyperprior shape: b_kappa ~ Gamma(alpha0, beta0)
    beta0_bkappa   = 0.01,   # hyperprior rate:  prior mean = alpha0 / beta0 = 100
    # Mixing moves: set FALSE to reproduce the results in the original paper (no mixing moves)
    use_mixing_moves = FALSE,
    # DEBUG: dump sampler state to disk on NaN / invariant failure
    debug_dump_dir = NULL
){
  ## ---- hygiene -----------------------------------------------------------
  op <- options(
    warnPartialMatchArgs   = TRUE,
    warnPartialMatchDollar = TRUE,
    warnPartialMatchAttr   = TRUE
  )
  on.exit(options(op), add = TRUE)
  
  if (!is.null(seed)) set.seed(as.integer(seed))
  psi_constraint <- match.arg(psi_constraint)
  sst_birth_score_mode <- match.arg(sst_birth_score_mode)
  psi_mode <- if (psi_constraint == "WST") "pair" else "distance"
  eta_identifiability <- match.arg(eta_identifiability)
  shrink_when <- match.arg(shrink_when)
  partition_prior <- .normalize_partition_prior(partition_prior)
  
  if (!all(free %in% c("psi","kappa","eta","z")))
    stop("`free` must be subset of {'psi','kappa','eta','z'}.", call. = FALSE)
  
  needs <- setdiff(c("psi","kappa","eta","z"), free)
  if (length(needs) > 0) {
    if (exists(".assert_truth", mode = "function", inherits = TRUE)) {
      .assert_truth(truth, needs)
    } else {
      if (is.null(truth) || is.na(truth)) {
        stop("`truth` must be a named list when parameters are fixed.", call. = FALSE)
      }
      miss <- setdiff(needs, names(truth))
      if (length(miss)) {
        stop("`truth` is missing: ", paste(miss, collapse = ", "), call. = FALSE)
      }
    }
  }
  
  .assert_scalar_int(K, "K")
  .assert_scalar_int(n_iter, "n_iter")
  .assert_scalar_int(burn, "burn")
  .assert_scalar_int(thin, "thin")
  if (burn >= n_iter) stop("`burn` must be < `n_iter`.", call. = FALSE)
  if (thin < 1L)      stop("`thin` must be >= 1.", call. = FALSE)
  
  n <- nrow(A)
  .assert_scalar_int(n, "n (nrow(A))")
  .assert_matrix(A, n, "A")
  
  if (!requireNamespace("BayesLogit", quietly = TRUE))
    stop("Package 'BayesLogit' is required.", call. = FALSE)
  if (psi_constraint == "SST" &&
      !exists("update_psi_sst", mode = "function", inherits = TRUE))
    stop("`psi_constraint='SST'` requested but `update_psi_sst()` is not visible.", call. = FALSE)
  
  ## ---- hyper-parameters --------------------------------------------------
  hyper <- list(
    a_kappa = a_kappa, b_kappa = b_kappa,
    a_eta   = a_eta,   b_eta   = b_eta,
    mu0     = mu0,     sigma0  = sigma0,
    tau0    = tau0,
    partition_prior = partition_prior,
    gamma_gn = gamma_gn,
    alpha_gp05 = alpha_gp05,
    theta_ocrp = theta_ocrp
  )
  # 
  # hyper <- calibrate_osbm_hypers(A,
  #                              rho0 = 0.60,
  #                              r0 = 0.65,
  #                              a_eta = 2.0,
  #                              v_shape = 20,
  #                              gn_gamma = 0.95)
  # after hyper is created, once per run
  # hyper$alpha_b <- 10
  # hyper$beta_b  <- hyper$alpha_b / hyper$b_kappa  # prior mean = initial b_kappa
  # 
  psi_hyper <- list(
    mu0 = hyper$mu0, sigma0 = hyper$sigma0, tau0 = hyper$tau0,
    alpha0 = hyper$alpha0, discount = hyper$discount,
    partition_prior = hyper$partition_prior,
    gamma_gn = hyper$gamma_gn,
    alpha_gp05 = hyper$alpha_gp05,
    theta_ocrp = hyper$theta_ocrp
  )
  # right after you build hyper / psi_hyper
  
  
  
  prior_type <- partition_prior
  
  ## ---- flags -------------------------------------------------------------
  psi_est   <- "psi"   %in% free
  z_est     <- "z"     %in% free
  eta_est   <- "eta"   %in% free
  kappa_est <- "kappa" %in% free
  
  ## ---- edge list ---------------------------------------------------------
  A <- as.matrix(A)
  
  # total matches per dyad
  N_mat <- A + t(A)
  diag(N_mat) <- 0
  
  idx   <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  i_idx <- idx[,1]
  j_idx <- idx[,2]
  N_edge <- as.numeric(N_mat[idx])   # total matches
  
  # keep A as wins matrix for directional use
  
  edge_by_node <- replicate(n, integer(0), simplify = FALSE)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], as.integer(e))
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], as.integer(e))
  }
  
  
  ## ---- initial values ----------------------------------------------------
  # Save the requested initial K before any override
  K_init <- as.integer(K)
  # z: use the passed K_init to create a balanced starting partition
  z <- if (z_est) {
    as.integer(sample(seq_len(K_init), n, replace = TRUE))
  } else {
    .assert_len(truth$z, n, "truth$z"); as.integer(truth$z)
  }
  # K is now the number of *occupied* blocks (may be < K_init by chance)
  K <- max(tabulate(z, nbins = K_init) > 0L) * K_init  # at most K_init
  K <- length(unique(z))  # exact occupied count
  # Relabel z to 1..K (contiguous)
  uv <- sort(unique(z))
  if (!identical(uv, seq_len(K))) {
    remap <- integer(max(uv)); remap[uv] <- seq_len(K)
    z <- remap[z]
  }
  # ψ
  if (psi_est) {
    if (psi_mode == "pair") {
      psi <- matrix(0, K, K); psi[upper.tri(psi)] <- 0.25
    } else {
      psi <- cumsum(rep(0.25, max(K-1L,0)))
    }
  } else {
    if (psi_mode == "pair") {
      .assert_matrix(truth$psi, K, "truth$psi (WST pair)"); psi <- truth$psi
    } else {
      .assert_len(truth$psi, K-1, "truth$psi (SST distance)"); psi <- as.numeric(truth$psi)
    }
  }
  
  # η
  if (eta_est) {
    eta <- rgamma(n, hyper$a_eta, hyper$b_eta)
    if (identical(eta_identifiability, "block_sum")) {
      eta <- eta_rescale_by_block(eta, z, K)
    }
  } else {
    .assert_len(truth$eta, n, "truth$eta"); eta <- truth$eta
  }
  
  # κ
  if (kappa_est) {
    p <- pmin(z[i_idx], z[j_idx]); q <- pmax(z[i_idx], z[j_idx])
    p <- factor(p, levels = 1:K);  q <- factor(q, levels = 1:K)
    Rkl <- as.matrix(xtabs(N_edge ~ p + q, drop.unused.levels = FALSE))
    
    E_k   <- tapply(eta,   factor(z, levels = 1:K), sum);   E_k[is.na(E_k)] <- 0
    eta2k <- tapply(eta^2, factor(z, levels = 1:K), sum); eta2k[is.na(eta2k)] <- 0
    
    Tkl <- outer(E_k, E_k, `*`); diag(Tkl) <- pmax((E_k^2 - eta2k) / 2, 0)
    kappa <- (Rkl + hyper$a_kappa) / (Tkl + hyper$b_kappa)
    kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]
  } else {
    .assert_matrix(truth$kappa, K, "truth$kappa"); kappa <- truth$kappa
  }
  
  ## ---- storage (ragged -> lists) -----------------------------------------
  keep_seq <- seq(burn + 1L, n_iter, by = thin)
  n_keep   <- length(keep_seq)
  
  draws_z     <- vector("list", n_keep)
  draws_kappa <- vector("list", n_keep)
  draws_eta   <- vector("list", n_keep)
  draws_psi   <- vector("list", n_keep)
  K_trace     <- integer(n_keep)
  b_kappa_trace <- numeric(n_keep)   # b_kappa value at each saved iteration
  
  keep <- 0L
  t0 <- proc.time()[3]
  
  A_ij <- A[cbind(i_idx, j_idx)]
  
  ## ============================================================
  ## Main Gibbs loop (refactored SST variable-K using LOO updater)
  ## ============================================================
  
  # Optional: only keep omega_e if you truly need it elsewhere (e.g. diagnostics)
  need_omega_e <- FALSE  # set TRUE only if some later code uses omega_e
  
  .fast_sst_dimcheck <- function(z, kappa, psi) {
    K <- nrow(kappa)
    if (ncol(kappa) != K) stop("kappa not square.", call. = FALSE)
    if (length(psi) != max(K - 1L, 0L)) stop("psi length != K-1.", call. = FALSE)
    if (max(z) != K) stop("max(z) != K (labels not contiguous).", call. = FALSE)
    K
  }
  
for (it in seq_len(n_iter)) {

  ## --- (0) cheap coherence repair ONLY if needed --------------------------
  K <- nrow(kappa)

  if (psi_mode == "distance") {  # ===== SST =====
    # Goal: maintain the invariants required by the SST parameterisation:
    #   - labels contiguous 1..K
    #   - kappa is KxK
    #   - psi is length K-1 (distance vector)
    #
    # Why "only if needed"?
    # pack_state_sst() is not free (it relabels, shrinks, etc.).
    # Doing it every iter would be wasted work and can hide bugs.
    if (ncol(kappa) != K || length(psi) != max(K - 1L, 0L) || max(z) != K) {
      packed <- pack_state_sst(z, kappa, psi)
      z <- packed$z; kappa <- packed$kappa; psi <- packed$psi
      K <- packed$K
    }

    # One cheap check per iteration, not per node i.
    # If this fails, you want it to fail early rather than 200 lines later.
    K <- .fast_sst_dimcheck(z, kappa, psi)

  } else {
    # ===== WST =====
    # In WST, psi is a KxK matrix (pairwise).
    # There's no "psi length must be K-1" constraint.
    # You still want contiguous labels, but if you're using the state-returning
    # updater + pack_state_wst() after sweeps, you can often skip packing here.
    K <- nrow(kappa)
  }

  # omega_e is an OPTIONAL cache of per-edge Polya–Gamma variables.
  #
  # omega_dirty tracks whether the cached omega_e no longer matches (z, psi).
  # We reset to FALSE at the start of each iteration because:
  #   - omega_e, if it exists, is assumed correct at this moment
  #   - we will set omega_dirty = TRUE if we change z or psi later in the iter
  #
  # Important: omega_dirty is meaningless if we do not keep omega_e at all.
  omega_dirty <- FALSE


  ## ---- (A) psi block ------------------------------------------------------
  if (psi_est) {

    if (psi_mode == "pair") {  # ===== WST ψ update =====

      K <- nrow(kappa)

      # Packing here ensures the WST invariants before updating psi:
      #   - psi must be KxK
      #   - max(z) must equal K
      #
      # Why do it here rather than only at the top of the iter?
      # Because the previous iteration might have changed K via births/deaths,
      # and you want psi update to see a coherent (z, psi, kappa) state.
      if (nrow(psi) != K || ncol(psi) != K || max(z) != K) {
        packed <- pack_state_wst(z, kappa, psi)
        z <- packed$z; kappa <- packed$kappa; psi <- packed$psi; K <- packed$K
      }

      # This function changes psi (and only psi) by sampling ψ_{kℓ} | rest.
      psi <- update_psi_wst_pair(
        K = K, z = z,
        i_idx = i_idx, j_idx = j_idx,
        A_ij = A_ij, N_edge = N_edge,
        psi_curr = psi,
        mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2
      )

      # Why set omega_dirty = TRUE after psi changes?
      # Because ω_e ~ PG(N_e, θ_e) and θ_e depends on ψ (and z).
      # If you store ω_e and intend to reuse it later in this iteration,
      # it is now out-of-date.
      #
      # Why guard with if (need_omega_e)?
      # Because if you do NOT store ω_e, there is nothing to refresh.
      if (need_omega_e) omega_dirty <- TRUE

    } else {  # ===== SST ψ update =====

      # In SST you update the distance vector ψ_d using aggregated-by-distance
      # statistics (bar_y, bar_omega). This does NOT require ω_e per-edge;
      # it uses bar_omega ~ PG(B_d, ψ_d) (distance-pooled omega).
      #
      # Even though this update doesn't use ω_e, changing ψ still invalidates
      # any cached ω_e, if you happen to keep ω_e for some other reason.
      agg   <- aggregate_by_distance(
        K,
        z_i = z[i_idx], z_j = z[j_idx],
        A_ij = A_ij, N_edge = N_edge
      )
      bar_y <- agg$bar_y

      B_d <- distance_totals(K, z_i = z[i_idx], z_j = z[j_idx], N_edge = N_edge)
      bar_omega <- draw_omega_bar(B_d = B_d, psi = psi)

      psi <- update_psi_sst(
        K = K, bar_y = bar_y, bar_omega = bar_omega,
        psi_curr = psi,
        mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2,
        tau2_0 = hyper$tau0^2, n_inner_sweeps = 4L
      )

      # Same reasoning: psi changed => cached omega_e is now stale *if it exists*.
      if (need_omega_e) omega_dirty <- TRUE
    }
  }


  ## ---- (B) z block (variable-K) ------------------------------------------
  if (z_est) {

    if (psi_mode == "pair") {  # ===== WST z-update (variable K) =====

      # Ensure coherent state before doing many single-site updates.
      K <- nrow(kappa)
      if (nrow(psi) != K || ncol(psi) != K || max(z) != K) {
        packed <- pack_state_wst(z, kappa, psi)
        z <- packed$z; kappa <- packed$kappa; psi <- packed$psi; K <- packed$K
      }

      # Draw edge-level PG latents ONCE for the entire z-sweep.
      # This is the PG-augmented approach: omega ~ PG(N_edge, |phi|).
      omega_edge <- draw_edge_omega(z, psi, i_idx, j_idx, N_edge, mode = "pair")

      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)

      slot_rad <- K + 1L

      node_order <- sample.int(n)
      for (i in node_order) {
        z_old_i <- z[i]
        K_old   <- K

        res <- wst_update_i_with_birth_LOO_pg(
          i = i,
          A = A, z = z, eta = eta, kappa = kappa, psi = psi,
          Rkl = bt$Rkl, Tkl = bt$Tkl,
          i_idx = i_idx, j_idx = j_idx, N_edge = N_edge, edge_by_node = edge_by_node,
          a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
          gamma_gn = psi_hyper$gamma_gn,
          mu0 = hyper$mu0, sigma0 = hyper$sigma0, sig2_0 = hyper$sigma0^2,
          omega_edge = omega_edge,
          slot_radius = slot_rad,
          partition_prior = psi_hyper$partition_prior,
          alpha_gp05 = psi_hyper$alpha_gp05,
          theta_ocrp = psi_hyper$theta_ocrp
        )

        z     <- res$z
        kappa <- res$kappa
        psi   <- res$psi
        K     <- res$K

        if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
          bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
          if (need_omega_e) omega_dirty <- TRUE
        }
      }

      packed <- pack_state_wst(z, kappa, psi)
      z <- packed$z; kappa <- packed$kappa; psi <- packed$psi; K <- packed$K

    } else if (psi_mode == "distance") {  # ===== SST z-update (variable K) =====

      # Draw edge-level PG latents ONCE for the entire z-sweep.
      omega_edge <- draw_edge_omega(z, psi, i_idx, j_idx, N_edge, mode = "distance")

      bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
      slot_rad <- K + 1

      node_order <- sample.int(n)
      for (i in node_order) {
        z_old_i <- z[i]
        K_old   <- K

        res <- sst_update_i_with_birth_LOO_pg(
          i = i,
          A = A, z = z, eta = eta, kappa = kappa, psi = psi,
          Rkl = bt$Rkl, Tkl = bt$Tkl,
          i_idx = i_idx, j_idx = j_idx, N_edge = N_edge, edge_by_node = edge_by_node,
          a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
          gamma_gn = psi_hyper$gamma_gn,
          omega_edge = omega_edge,
          slot_radius = slot_rad,
          tau0 = psi_hyper$tau0,
          mu0 = psi_hyper$mu0,
          sig2_0 = psi_hyper$sigma0^2,
          partition_prior = psi_hyper$partition_prior,
          alpha_gp05 = psi_hyper$alpha_gp05,
          theta_ocrp = psi_hyper$theta_ocrp,
          birth_score_mode = sst_birth_score_mode
        )

        z     <- res$z
        kappa <- res$kappa
        psi   <- res$psi
        K     <- res$K

        if (K != K_old || z[i] != z_old_i || isTRUE(res$changed)) {
          bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
          if (need_omega_e) omega_dirty <- TRUE
        }
      }

      # SST post-sweep invariant check
      K <- .fast_sst_dimcheck(z, kappa, psi)
    }
  }


  ## ---- (B2) Adjacent-block swap MH move -----------------------------------
  # Proposes swapping two adjacent blocks k <-> k+1.
  # Symmetric proposal; only directional likelihood changes.
  if (use_mixing_moves && z_est && K >= 2L) {
    swap_res <- adjacent_block_swap_move(
      z = z, kappa = kappa, psi = psi, eta = eta,
      A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      psi_mode = psi_mode
    )
    z <- swap_res$z; kappa <- swap_res$kappa; psi <- swap_res$psi
    if (isTRUE(swap_res$accepted) && need_omega_e) omega_dirty <- TRUE
  }


  ## ---- (B3) Split-merge MH move ------------------------------------------
  # Proposes splitting a block in two or merging two adjacent blocks.
  if (use_mixing_moves && z_est && K >= 1L) {
    sm_res <- split_merge_move(
      z = z, kappa = kappa, psi = psi, eta = eta,
      A = A, i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
      edge_by_node = edge_by_node,
      a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa,
      gamma_gn = psi_hyper$gamma_gn,
      psi_mode = psi_mode,
      hyper_psi = list(mu0 = hyper$mu0, sigma0 = hyper$sigma0, tau0 = hyper$tau0),
      n_restricted_scans = 3L,
      partition_prior = psi_hyper$partition_prior,
      theta_ocrp = psi_hyper$theta_ocrp
    )
    z <- sm_res$z; kappa <- sm_res$kappa; psi <- sm_res$psi
    if (isTRUE(sm_res$accepted)) {
      K <- if (psi_mode == "pair") nrow(kappa) else length(psi) + 1L
      if (need_omega_e) omega_dirty <- TRUE
    }
  }

  if (eta_est && identical(eta_identifiability, "block_sum")) {
    eta <- eta_rescale_by_block(eta, z, K)
  }


  ## ---- (C) omega_e (ONLY if you keep it) ---------------------------------
  # This is the *only* place where omega_e is actually refreshed.
  #
  # Why here?
  # Because you want to avoid expensive draws unless you have to,
  # and you want ω_e to correspond to the *final* (z, psi) state
  # after both psi-update and z-sweep in the current iteration.
  #
  # Why the AND condition?
  # - need_omega_e: only refresh if somebody will actually use omega_e later
  # - omega_dirty: only refresh if z/psi changed this iter
  if (need_omega_e && omega_dirty) {
    omega_e <- draw_pg_latents(i_idx, j_idx, N_edge, z, psi, mode = psi_mode)
    omega_dirty <- FALSE
  }


  ## ---- (D) kappa block ----------------------------------------------------
  if (kappa_est) {
    # κ depends on (z, eta) through Rkl, Tkl.
    # After the z-sweep (and possibly eta update later), you rebuild totals
    # and sample κ fresh.
    bt  <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)
    Rkl <- bt$Rkl
    Tkl <- bt$Tkl

    kappa <- matrix(0, K, K)
    # Vectorized kappa sampling from Gamma posteriors
    shape_mat <- hyper$a_kappa + Rkl
    rate_mat  <- hyper$b_kappa + Tkl
    # Fill upper triangle (including diagonal)
    ut_idx <- which(upper.tri(kappa, diag = TRUE), arr.ind = TRUE)
    kappa[ut_idx] <- rgamma(nrow(ut_idx),
                            shape = shape_mat[ut_idx],
                            rate  = rate_mat[ut_idx])
    # Symmetrize
    kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]
  }

  ## ---- (D2) b_kappa update (hierarchical, optional) ----------------------
  # Partially collapsed Gibbs: collapse kappa for z-updates (step B), then
  # un-collapse here to draw b_kappa | kappa ~ Gamma(alpha0 + P*a_kappa,
  #                                                    beta0  + sum(kappa_{ut}))
  # where P = K*(K+1)/2 (number of unique block pairs).
  # Self-regulation: K↑ => more kappa params => sum(kappa)↑ => b_kappa↑ => births suppressed.
  if (sample_b_kappa && kappa_est) {
    hyper$b_kappa <- update_b_kappa(kappa,
                                    a_kappa = hyper$a_kappa,
                                    alpha_b = alpha0_bkappa,
                                    beta_b  = beta0_bkappa)
  }


  ## ---- (E) eta block ------------------------------------------------------
  if (eta_est) {
    # η depends on κ and z; it’s sampled after κ in this design.
    eta <- update_eta_all_dyads(
      eta = eta, z = z, kappa = kappa,
      a_eta = hyper$a_eta, b_eta = hyper$b_eta,
      i_idx = i_idx, j_idx = j_idx, N_edge = N_edge
    )
    if (identical(eta_identifiability, "block_sum")) {
      eta <- eta_rescale_by_block(eta, z, K)
    }
  }


  ## ---- save draws ---------------------------------------------------------
  if (it %in% keep_seq) {
    keep <- keep + 1L
    draws_z[[keep]]     <- as.integer(z)
    draws_kappa[[keep]] <- kappa
    draws_eta[[keep]]   <- as.numeric(eta)
    draws_psi[[keep]]   <- psi
    K_trace[keep]       <- nrow(kappa)
    b_kappa_trace[keep] <- hyper$b_kappa
  }

  # DEBUG: check for NaN/Inf in psi after each iteration
  if (!is.null(debug_dump_dir) || verbose) {
    bad_psi <- if (psi_mode == "distance") {
      anyNA(psi) || any(!is.finite(psi))
    } else {
      anyNA(psi) || any(!is.finite(psi[upper.tri(psi)]))
    }
    if (bad_psi) {
      msg <- sprintf("[DEBUG|it=%d] psi contains NA/Inf/NaN (K=%d, mode=%s)",
                     it, nrow(kappa), psi_mode)
      if (!is.null(debug_dump_dir)) {
        dir.create(debug_dump_dir, showWarnings = FALSE, recursive = TRUE)
        fn <- file.path(debug_dump_dir,
                        sprintf("bad_psi_it%05d.rds", it))
        saveRDS(list(it = it, K = nrow(kappa), z = z, psi = psi,
                     kappa = kappa, eta = eta, mode = psi_mode), fn)
        message(msg, "  -> state dumped to ", fn)
      } else {
        warning(msg)
      }
      # Hard repair: reset bad psi entries to small positive
      if (psi_mode == "distance") {
        psi[!is.finite(psi)] <- 0.01
        psi <- pmax(psi, 0)
      } else {
        psi[!is.finite(psi)] <- 0
        psi <- pmax(psi, 0)
      }
    }
  }

  if (verbose && it %% 20 == 0) {
    cat("iter", it, "/", n_iter, " saved:", keep, "  K=", nrow(kappa), "  ",
        round(proc.time()[3] - t0, 2), "s \r")
  }
}

  
  
  
  invisible(list(
    z        = draws_z,
    psi      = draws_psi,
    kappa    = draws_kappa,
    eta      = draws_eta,
    K_trace  = K_trace,
    b_kappa_trace = b_kappa_trace,
    keep     = keep_seq,
    meta     = list(psi_mode = psi_mode, psi_constraint = psi_constraint,
                    hyper = hyper, free = free,
                    eta_identifiability = eta_identifiability,
                    shrink_when = shrink_when,
                    refresh_pg_after_birth = refresh_pg_after_birth,
                    sample_b_kappa = sample_b_kappa,
                    alpha0_bkappa = alpha0_bkappa,
                    beta0_bkappa  = beta0_bkappa)
  ))
}




# ===== Minimal assertion helpers (no extra packages) =========================
