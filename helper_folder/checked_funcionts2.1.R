exercise_all_slots_wst <- function(z, kappa, psi_mat) {
  K <- nrow(kappa); out <- list()
  for (r in seq_len(K + 1L)) {
    ins  <- insert_block_labels(z, K, r)
    init_kap <- rep(0.1, K)                 # placeholder: will be replaced by collapsed/posteriors in Step 3
    init_psi <- rep(0.5, K)                 # placeholder (WST half-Normal prior mean-ish)
    grown <- grow_params_wst(kappa, psi_mat, r, init_kappa_rowcol = init_kap, init_psi_row = init_psi)
    # set label for a hypothetical new node i at slot r
    z_new <- ins$z_new
    # invariants
    stopifnot(check_invariants(z_new, grown$kappa, grown$psi, "WST"))
    out[[r]] <- list(r = r, z = z_new, kappa = grown$kappa, psi = grown$psi)
  }
  out
}

exercise_all_slots_sst <- function(z, kappa, psi) {
  K <- nrow(kappa); out <- list()
  for (r in seq_len(K + 1L)) {
    ins   <- insert_block_labels(z, K, r)
    grown <- grow_sst_params(kappa, psi, r)  # placeholder
    stopifnot(check_invariants(ins$z_new, grown$kappa, grown$psi[-length(grown$psi)], "SST"))
    out[[r]] <- list(r = r, z = ins$z_new, kappa = grown$kappa, psi = grown$psi)
  }
  out
}
## -------- slots policy (burn-in: all; post burn-in: adjacent) ----------
slot_set <- function(K, k_prime = NULL, burnin = TRUE) {
  if (burnin || is.null(k_prime)) return(seq_len(K + 1L))
  c(max(1L, k_prime), min(K + 1L, k_prime + 1L))
}

## -------- relabel: insert empty block at slot r (labels only) ----------
# Returns list(z_new, K_new, perm_old2new) where perm maps old labels -> new
insert_block_labels <- function(z, K, r) {
  stopifnot(r >= 1L, r <= K + 1L)
  # shift labels >= r upward by +1
  z_new <- ifelse(z >= r, z + 1L, z)
  # identity on [1..r-1], then a hole at r, then +1 for the tail
  perm <- seq_len(K)
  perm_new <- integer(K)
  # old labels < r keep order
  keep_lt <- which(perm < r)
  if (length(keep_lt)) perm_new[keep_lt] <- perm[keep_lt]
  # old labels >= r shift by +1
  keep_ge <- which(perm >= r)
  if (length(keep_ge)) perm_new[keep_ge] <- perm[keep_ge] + 1L
  list(z_new = z_new, K_new = K + 1L, perm_old2new = perm_new)
}

## -------- utilities to insert row/col into KxK matrices ----------------
insert_rowcol_sym <- function(M, r, fill = 0) {
  K  <- nrow(M); stopifnot(ncol(M) == K)
  M2 <- matrix(fill, nrow = K + 1L, ncol = K + 1L)
  # top-left (r-1)x(r-1)
  if (r > 1L) M2[1:(r-1), 1:(r-1)] <- M[1:(r-1), 1:(r-1)]
  # bottom-right (K-r+1)x(K-r+1)
  if (r <= K) M2[(r+1):(K+1), (r+1):(K+1)] <- M[r:K, r:K]
  # top-right / bottom-left blocks
  if (r > 1L && r <= K) {
    M2[1:(r-1), (r+1):(K+1)] <- M[1:(r-1), r:K]
    M2[(r+1):(K+1), 1:(r-1)] <- t(M[1:(r-1), r:K])
  }
  # keep symmetry (diagonal left as `fill` for now)
  M2
}

## -------- WST: grow κ and ψ (pairwise ψ_{kℓ}) --------------------------
# psi_mat is KxK with meaningful upper-tri; κ is symmetric KxK
# init_kappa_rowcol: length-K numeric for new row/col (excluding diagonal)
# init_psi_row:      length-K numeric for ψ_{rℓ} at ℓ = 1..K (we’ll only keep ℓ>r)
grow_params_wst <- function(kappa, psi_mat, r,
                            init_kappa_rowcol = NULL,
                            init_psi_row = NULL,
                            diag_kappa = 0) {
  K <- nrow(kappa); stopifnot(ncol(kappa) == K, nrow(psi_mat) == K, ncol(psi_mat) == K)
  # expand κ
  kap2 <- insert_rowcol_sym(kappa, r, fill = diag_kappa)
  if (is.null(init_kappa_rowcol)) init_kappa_rowcol <- rep(NA_real_, K)
  stopifnot(length(init_kappa_rowcol) == K)
  # fill new row/col (excluding diagonal)
  if (r > 1L) {
    kap2[r, 1:(r-1)] <- init_kappa_rowcol[1:(r-1)]
    kap2[1:(r-1), r] <- init_kappa_rowcol[1:(r-1)]
  }
  if (r <= K) {
    kap2[r, (r+1):(K+1)] <- init_kappa_rowcol[r:K]
    kap2[(r+1):(K+1), r] <- init_kappa_rowcol[r:K]
  }
  
  # expand ψ (pairwise, upper-tri only meaningful)
  psi2 <- insert_rowcol_sym(psi_mat, r, fill = 0)
  if (is.null(init_psi_row)) init_psi_row <- rep(NA_real_, K)
  stopifnot(length(init_psi_row) == K)
  # set ψ_{min,max} with positivity on upper-tri
  # upper-tri entries at (min(r,ℓ), max(r,ℓ))
  for (ell in seq_len(K)) {
    if (ell == r) next
    a <- min(r, ifelse(ell >= r, ell + 1L, ell))
    b <- max(r, ifelse(ell >= r, ell + 1L, ell))
    if (!is.na(init_psi_row[ell])) psi2[a, b] <- max(init_psi_row[ell], 0)  # enforce WST sign
  }
  list(kappa = kap2, psi = psi2)
}

## -------- SST: grow κ and ψ (distance ψ_d) ------------------------------

## -------- invariants / sanity checks ------------------------------------
check_invariants <- function(z, kappa, psi, regime = c("WST","SST")) {
  regime <- match.arg(regime)
  K <- max(z)
  stopifnot(nrow(kappa) == K, ncol(kappa) == K)
  stopifnot(all(abs(kappa - t(kappa)) < 1e-10))  # symmetric κ
  if (regime == "WST") {
    stopifnot(is.matrix(psi), nrow(psi) == K, ncol(psi) == K)
    ut <- which(upper.tri(psi), arr.ind = TRUE)
    if (nrow(ut)) stopifnot(all(psi[ut] >= 0))     # half-Normal truncation
  } else {
    stopifnot(is.numeric(psi), length(psi) == K - 1L)
    # optional monotonic check (allow tiny jitter)
    if (length(psi)) stopifnot(all(diff(psi) >= -1e-10))
  }
  TRUE
}
## =======================
## Step 2 — DRY RUN
## =======================

## --- If you haven't sourced these helpers yet, keep them here ----------
slot_set <- function(K, k_prime = NULL, burnin = TRUE) {
  if (burnin || is.null(k_prime)) return(seq_len(K + 1L))
  c(max(1L, k_prime), min(K + 1L, k_prime + 1L))
}

insert_block_labels <- function(z, K, r) {
  stopifnot(r >= 1L, r <= K + 1L)
  z_new <- ifelse(z >= r, z + 1L, z)
  perm <- seq_len(K)
  perm_new <- integer(K)
  keep_lt <- which(perm < r)
  if (length(keep_lt)) perm_new[keep_lt] <- perm[keep_lt]
  keep_ge <- which(perm >= r)
  if (length(keep_ge)) perm_new[keep_ge] <- perm[keep_ge] + 1L
  list(z_new = z_new, K_new = K + 1L, perm_old2new = perm_new)
}

insert_rowcol_sym <- function(M, r, fill = 0) {
  K  <- nrow(M); stopifnot(ncol(M) == K)
  M2 <- matrix(fill, nrow = K + 1L, ncol = K + 1L)
  if (r > 1L) M2[1:(r-1), 1:(r-1)] <- M[1:(r-1), 1:(r-1)]
  if (r <= K) M2[(r+1):(K+1), (r+1):(K+1)] <- M[r:K, r:K]
  if (r > 1L && r <= K) {
    M2[1:(r-1), (r+1):(K+1)] <- M[1:(r-1), r:K]
    M2[(r+1):(K+1), 1:(r-1)] <- t(M[1:(r-1), r:K])
  }
  M2
}

grow_params_wst <- function(kappa, psi_mat, r,
                            init_kappa_rowcol = NULL,
                            init_psi_row = NULL,
                            diag_kappa = 0) {
  K <- nrow(kappa); stopifnot(ncol(kappa) == K, nrow(psi_mat) == K, ncol(psi_mat) == K)
  kap2 <- insert_rowcol_sym(kappa, r, fill = diag_kappa)
  if (is.null(init_kappa_rowcol)) init_kappa_rowcol <- rep(NA_real_, K)
  stopifnot(length(init_kappa_rowcol) == K)
  if (r > 1L) {
    kap2[r, 1:(r-1)] <- init_kappa_rowcol[1:(r-1)]
    kap2[1:(r-1), r] <- init_kappa_rowcol[1:(r-1)]
  }
  if (r <= K) {
    kap2[r, (r+1):(K+1)] <- init_kappa_rowcol[r:K]
    kap2[(r+1):(K+1), r] <- init_kappa_rowcol[r:K]
  }
  
  psi2 <- insert_rowcol_sym(psi_mat, r, fill = 0)
  if (is.null(init_psi_row)) init_psi_row <- rep(NA_real_, K)
  stopifnot(length(init_psi_row) == K)
  for (ell in seq_len(K)) {
    if (ell == r) next
    a <- min(r, ifelse(ell >= r, ell + 1L, ell))
    b <- max(r, ifelse(ell >= r, ell + 1L, ell))
    if (!is.na(init_psi_row[ell])) psi2[a, b] <- max(init_psi_row[ell], 0)
  }
  list(kappa = kap2, psi = psi2)
}

# --- SST growth: add one slot at position r and extend ψ monotonically ---
# kappa: KxK symmetric
# psi   : length K-1, nondecreasing
# r     : slot in {1, ..., K+1}
# returns list(kappa=K+1 x K+1, psi=length K)
grow_sst_params <- function(kappa, psi, r,
                            a_kappa = 1, b_kappa = 1,
                            delta_extreme = 1e-6) {
  stopifnot(is.matrix(kappa), nrow(kappa) == ncol(kappa),
            is.numeric(psi), length(psi) == nrow(kappa) - 1L)
  K  <- nrow(kappa)
  Kn <- K + 1L
  
  # ---- grow κ symmetrically with Gamma(a,b) prior for new row/col ----
  kap_new <- matrix(0, Kn, Kn)
  kap_new[1:K, 1:K] <- kappa
  # sample new interactions κ_{r,ℓ} (sym)
  for (ell in 1:Kn) {
    if (ell == r) next
    # small prior draw; in your actual sampler you’ll use collapsed or data-driven updates
    val <- rgamma(1, shape = a_kappa, rate = b_kappa)
    kap_new[r, ell] <- val
    kap_new[ell, r] <- val
  }
  # keep symmetry
  diag(kap_new) <- diag(kap_new) # noop but explicit
  
  # ---- grow ψ (distance-based) while preserving monotonicity ----
  # Old ψ has length K-1; new ψ should have length K (distances 1..K).
  psi_new <- numeric(Kn - 1L)
  
  # Distances 1..(K-1) reuse old ψ exactly (middle insertions don’t need new values here)
  psi_new[1:(K-1)] <- psi
  
  # New maximal distance: K_new-1 = K
  # If r at extremes -> create a strictly larger (or slightly larger) value;
  # If r in middle   -> reuse last value (equality is allowed and keeps monotone).
  is_extreme <- (r == 1L || r == Kn)
  inc <- if (is_extreme) abs(delta_extreme) else 0.0
  psi_new[K] <- psi[K-1] + inc
  
  list(kappa = kap_new, psi = psi_new)
}

# Utility: insert an empty block at slot r in labels (no data movement)
insert_empty_block <- function(z, r) {
  K <- max(z)
  stopifnot(r >= 1L, r <= K + 1L)
  z + as.integer(z >= r)
}


check_invariants <- function(z, kappa, psi, regime = c("WST","SST")) {
  regime <- match.arg(regime)
  
  # Parameter dimension is authoritative (may exceed max(z) if new empty block exists)
  K_par <- nrow(kappa)
  stopifnot(ncol(kappa) == K_par)
  
  # Labels must be within 1..K_par (allow empty blocks)
  stopifnot(length(z) >= 1L)
  stopifnot(all(z >= 1L & z <= K_par))
  
  # κ symmetry
  stopifnot(all(abs(kappa - t(kappa)) < 1e-10))
  
  if (regime == "WST") {
    # ψ must be K_par x K_par; enforce nonnegativity on upper triangle
    stopifnot(is.matrix(psi), nrow(psi) == K_par, ncol(psi) == K_par)
    ut <- which(upper.tri(psi), arr.ind = TRUE)
    if (nrow(ut)) stopifnot(all(psi[ut] >= -1e-12))
  } else { # SST
    # Accept ψ length K_par−1 (standard) or K_par (when carrying a ψ_K placeholder)
    stopifnot(is.numeric(psi))
    stopifnot(length(psi) %in% c(K_par - 1L, K_par))
    v <- if (length(psi) == K_par) psi[-length(psi)] else psi
    if (length(v)) stopifnot(all(diff(v) >= -1e-10))  # monotone nondecreasing on used distances
  }
  TRUE
}


## --- Tiny toy states ---------------------------------------------------
make_toy_state <- function(K = 3L, n = 15L, regime = c("WST","SST"), seed = 1L) {
  regime <- match.arg(regime)
  set.seed(seed)
  z <- sample.int(K, n, replace = TRUE)
  # κ symmetric, positive
  B <- matrix(rexp(K*K, rate = 1), K, K); B <- (B + t(B))/2
  diag(B) <- diag(B) + 0.5
  if (regime == "WST") {
    psi <- matrix(0, K, K)
    ut  <- which(upper.tri(psi), arr.ind = TRUE)
    psi[ut] <- pmax(rnorm(nrow(ut), mean = 0.6, sd = 0.2), 0) # half-Normal-ish
  } else {
    psi <- cumsum(pmax(rnorm(K-1, mean = 0.5, sd = 0.15), 0.05)) # monotone
  }
  list(z = z, kappa = B, psi = psi)
}

## --- Dry run per regime: try ALL slots and check -----------------------
dry_run_step2_wst <- function(state, verbose = TRUE) {
  z <- state$z; kappa <- state$kappa; psi <- state$psi
  K <- nrow(kappa)
  if (verbose) cat("WST dry-run over slots 1..", K+1, "\n", sep = "")
  res <- vector("list", K + 1L)
  for (r in seq_len(K + 1L)) {
    ins  <- insert_block_labels(z, K, r)
    init_kap <- rep(0.1, K)
    init_psi <- rep(0.4, K)
    grown <- grow_params_wst(kappa, psi, r, init_kappa_rowcol = init_kap, init_psi_row = init_psi)
    ok <- check_invariants(ins$z_new, grown$kappa, grown$psi, "WST")
    if (verbose) {
      nz <- tabulate(ins$z_new, nbins = K + 1L)
      ut <- which(upper.tri(grown$psi), arr.ind = TRUE)
      minpsi <- if (nrow(ut)) min(grown$psi[ut]) else 0
      cat(sprintf("  slot r=%d: OK | K->%d | block sizes: %s | min ψ(upper)=%.3f\n",
                  r, K+1, paste(nz, collapse = " "), minpsi))
    }
    res[[r]] <- list(r = r, z = ins$z_new, kappa = grown$kappa, psi = grown$psi)
  }
  invisible(res)
}

dry_run_step2_sst <- function(state, verbose = TRUE) {
  z <- state$z; kappa <- state$kappa; psi <- state$psi
  K <- nrow(kappa)
  if (verbose) cat("SST dry-run over slots 1..", K+1, "\n", sep = "")
  res <- vector("list", K + 1L)
  for (r in seq_len(K + 1L)) {
    ins   <- insert_block_labels(z, K, r)
    grown <- grow_sst_params(kappa, psi, r)
    # For invariants we check monotonicity on the *existing* distances only
    ok <- check_invariants(ins$z_new, grown$kappa, grown$psi, "SST")
    if (verbose) {
      nz <- tabulate(ins$z_new, nbins = K + 1L)
      # report monotonicity on first K-1 entries
      v <- grown$psi
      msg_extra <- if (r %in% c(1, K+1)) " (ψ_K placeholder present)" else " (middle slot: ψ_K unused)"
      cat(sprintf("  slot r=%d: OK | K->%d | block sizes: %s | ψ[1..%d] monotone%s\n",
                  r, K+1, paste(nz, collapse = " "), K-1, msg_extra))
    }
    res[[r]] <- list(r = r, z = ins$z_new, kappa = grown$kappa, psi = grown$psi)
  }
  invisible(res)
}

## --- One-liner to run both --------------------------------------------
run_step2_dryrun <- function(K = 3L, n = 15L, seed = 1L) {
  st_w <- make_toy_state(K, n, "WST", seed)
  st_s <- make_toy_state(K, n, "SST", seed + 1L)
  cat("== DRY RUN: WST ==\n")
  out_w <- dry_run_step2_wst(st_w, verbose = TRUE)
  cat("== DRY RUN: SST ==\n")
  out_s <- dry_run_step2_sst(st_s, verbose = TRUE)
  invisible(list(WST = out_w, SST = out_s))
}
out <- run_step2_dryrun(K = 4, n = 20, seed = 42)
# == DRY RUN: WST ==
# WST dry-run over slots 1..5
#   slot r=1: OK | K->5 | block sizes: 6 6 3 5 0 | min ψ(upper)=0.000
#   slot r=2: OK | K->5 | block sizes: 7 0 6 3 4 | min ψ(upper)=0.000
#   ...
# == DRY RUN: SST ==
# SST dry-run over slots 1..5
#   slot r=1: OK | K->5 | block sizes: ... | ψ[1..3] monotone (ψ_K placeholder present)
#   slot r=2: OK | K->5 | block sizes: ... | ψ[1..3] monotone (middle slot: ψ_K unused)
#   ...
