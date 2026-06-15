#!/usr/bin/env Rscript
# ============================================================================
# diagnose_sst_births.R  —  Expose WHY K proliferates in SST
#
# Strategy:
#   1. GN prior test: verify weights & confirm gamma_gn effect
#   2. Collapsed-kappa (volume) asymmetry test
#   3. Full birth/join decomposition: log each component per node
#   4. Run on toy + real datasets, tabulate which component drives births
# ============================================================================

cat("========== SST BIRTH DIAGNOSTICS ==========\n\n")

suppressPackageStartupMessages({
  source("helper_folder/helper.R")
  source("helper_folder/SST_helpers.R")
  source("helper_folder/WST_helpers.R")
  source("helper_folder/SST_tester.R")
  source("core/my_best_try_so_far.R")
})

# ============================================================================
# TEST 1:  GN prior weights — does gamma actually penalize births?
# ============================================================================
cat("============================================================\n")
cat("TEST 1: GN prior weights for different gamma and H\n")
cat("============================================================\n\n")

for (gamma in c(0.3, 0.5, 0.8, 0.9, 0.95, 0.99)) {
  cat(sprintf("  gamma = %.2f\n", gamma))
  for (H in c(2, 3, 5, 10, 15, 20)) {
    # v_minus: H equally-sized blocks of 3 nodes each (minus node i)
    v_minus <- rep(3, H)
    pw <- gn_log_weights_packed(v_minus, gamma_gn = gamma)

    # Probability of birth = softmax(new) / (softmax(new) + softmax(exist))
    # Here we just compare raw log-weights
    w_ex_total <- lse(pw$exist)
    w_new <- pw$new

    # Birth probability (prior only, no likelihood)
    p_birth <- exp(w_new - lse(c(w_new, w_ex_total)))
    cat(sprintf("    H=%2d: log(w_new)=%7.2f  log(sum_w_exist)=%7.2f  P(birth|prior)=%.4f\n",
                H, w_new, w_ex_total, p_birth))
  }
  cat("\n")
}

# Key observation: w_new = H^2 - H*gamma grows as H^2.
# Even gamma=0.99 only subtracts ~0.99*H.
# So P(birth|prior) INCREASES with H — the prior ENCOURAGES births when K is big!
cat("CONCLUSION: GN new-block weight = H^2 - H*gamma grows like H^2.\n")
cat("  => P(birth|prior) INCREASES with H (# occupied blocks).\n")
cat("  => gamma only gives a LINEAR correction; quadratic growth dominates.\n")
cat("  => The prior CANNOT prevent K proliferation once K is moderately large.\n\n")


# ============================================================================
# TEST 2:  Collapsed kappa asymmetry — new vs existing
# ============================================================================
cat("============================================================\n")
cat("TEST 2: Collapsed kappa (volume) — new vs existing asymmetry\n")
cat("============================================================\n\n")

# The new-block kappa score is:
#   lp_kappa_new = sum_{ell} gp_marginal(r_add[ell], t_add[ell], a, b)
#
# This uses (R=0, T=0) as the base for the new block pair, i.e. it evaluates
# the GP marginal from scratch for each (node_i, block_ell) pair.
#
# The existing-block kappa score is:
#   lp_kappa_exist[k] = sum_{ell} [gp_marginal(R_minus[k,ell] + r_add[ell], ...)
#                                 - gp_marginal(R_minus[k,ell], ...)]
#
# KEY ASYMMETRY: For existing blocks, the score is a DELTA (marginal gain of
# adding node i). For the new block, the score is the FULL marginal (not a delta).
#
# If a_kappa is small (e.g. 1) and b_kappa is small (e.g. 0.1):
#   gp_marginal(R, T, 1, 0.1) = lgamma(1+R) - lgamma(1) + 1*log(0.1) - (1+R)*log(0.1+T)
#                               = log(R!) - (1+R)*log(0.1+T) + log(0.1)
#
# For a new block with ONE node: r_add might be small, t_add = eta_i * E_ell.
# gp_marginal(r_add, t_add, a, b) includes the FULL lgamma(a+R) - lgamma(a) + a*log(b) term.
# For existing blocks, the +a*log(b) cancels in the delta!
#
# This means lp_kappa_new gets a BONUS of K_minus * a_kappa * log(b_kappa) that
# existing blocks do NOT get (it cancels in the delta).
# When a_kappa * log(b_kappa) is negative (b_kappa < 1), this is a penalty.
# When b_kappa >= 1, it's a bonus.
# But importantly, the lgamma(a+R)-lgamma(a) terms also don't cancel for new blocks.

cat("Illustrating the asymmetry with typical parameters:\n\n")
for (a_k in c(1, 5, 50)) {
  for (b_k in c(0.1, 1, 10)) {
    # Simulate: 5 existing blocks, node i has moderate connections
    K_test <- 5
    r_add <- c(3, 5, 2, 8, 1)  # edges from i to each block
    t_add <- c(2, 3, 1.5, 4, 0.5) # exposure from i to each block

    # Base R/T for existing blocks (moderate counts)
    R_base <- matrix(20, K_test, K_test)
    T_base <- matrix(10, K_test, K_test)

    # Existing block scores (pick block 1)
    kappa_exist <- numeric(K_test)
    for (kc in 1:K_test) {
      acc <- 0
      for (ell in 1:K_test) {
        p <- min(kc, ell); q <- max(kc, ell)
        acc <- acc + (gp_marginal(R_base[p,q] + r_add[ell], T_base[p,q] + t_add[ell], a_k, b_k) -
                        gp_marginal(R_base[p,q], T_base[p,q], a_k, b_k))
      }
      kappa_exist[kc] <- acc
    }

    # New block score
    kappa_new <- sum(gp_marginal(r_add, t_add, a_k, b_k))

    cat(sprintf("  a_kappa=%.0f, b_kappa=%.1f: best_exist=%7.2f  new_block=%7.2f  diff(new-best)=%7.2f\n",
                a_k, b_k, max(kappa_exist), kappa_new, kappa_new - max(kappa_exist)))
  }
}

cat("\nCONCLUSION: When a_kappa is small and b_kappa is small,\n")
cat("  the new-block kappa score can systematically exceed the best existing score.\n")
cat("  This is because the 'full marginal' includes lgamma(a+R) bonuses\n")
cat("  that cancel in the 'delta' computation for existing blocks.\n\n")


# ============================================================================
# TEST 3: SST directional tester (existing fast-vs-brute checks)
# ============================================================================
cat("============================================================\n")
cat("TEST 3: SST directional score correctness (fast vs brute)\n")
cat("============================================================\n\n")

toy <- simulate_sst_toy(n = 14L, K = 4L, lambda_N = 3, seed = 42L)
cat("Toy psi:", paste(round(toy$psi, 3), collapse = ", "), "\n")
tst <- run_fast_bruteforce_tests(
  A = toy$A, N = toy$N, z = toy$z, psi_vec = toy$psi,
  tau0 = 0.15, m_inc = 0, tol = 1e-8, verbose = TRUE
)


# ============================================================================
# TEST 4: Full birth/join decomposition on a real dataset
# ============================================================================
cat("\n============================================================\n")
cat("TEST 4: Full birth/join decomposition during z-sweep\n")
cat("============================================================\n\n")

# Instrumented version of the SST z-update that logs all components
sst_update_i_instrumented <- function(
    i,
    A, z, eta, kappa, psi,
    Rkl, Tkl,
    i_idx, j_idx, N_edge, edge_by_node,
    a_kappa, b_kappa,
    gamma_gn,
    tau0, mu0 = 0, sig2_0 = 1
) {

  K_full <- nrow(kappa)
  oldk  <- z[i]
  eta_i <- eta[i]

  n_minus_full <- tabulate(z[-i], nbins = K_full)
  keep_full <- sort(which(n_minus_full > 0L))
  K_minus <- length(keep_full)
  if (K_minus < 1L) return(NULL)

  map <- integer(K_full); map[keep_full] <- seq_len(K_minus)
  z_packed <- map[z]
  oldk_pos <- match(oldk, keep_full)

  kappa_minus <- kappa[keep_full, keep_full, drop = FALSE]
  psi_minus   <- reindex_psi_sst_keep(psi_old = psi, K_old = K_full, keep = keep_full)

  C_i_full <- counts_by_block_exact_cpp(
    i = i, A = A, z = as.integer(z),
    i_idx = as.integer(i_idx), j_idx = as.integer(j_idx),
    N_edge = as.numeric(N_edge),
    edge_by_node = edge_by_node, K = K_full
  )
  c_plus <- as.numeric(C_i_full$c_plus)[keep_full]
  N_tot  <- as.numeric(C_i_full$N_tot )[keep_full]

  E_excl_full <- E_by_block_excluding_i(i, z, eta, K_full)
  E_excl <- as.numeric(E_excl_full[keep_full])

  r_add <- N_tot
  t_add <- eta_i * E_excl

  R_minus <- Rkl[keep_full, keep_full, drop = FALSE]
  T_minus <- Tkl[keep_full, keep_full, drop = FALSE]

  if (!is.na(oldk_pos)) {
    for (ell_pos in seq_len(K_minus)) {
      subR <- N_tot[ell_pos]; subT <- eta_i * E_excl[ell_pos]
      p <- min(oldk_pos, ell_pos); q <- max(oldk_pos, ell_pos)
      R_minus[p, q] <- R_minus[p, q] - subR
      T_minus[p, q] <- T_minus[p, q] - subT
      R_minus[q, p] <- R_minus[p, q]; T_minus[q, p] <- T_minus[p, q]
      if (R_minus[p, q] < 0) R_minus[p, q] <- R_minus[q, p] <- 0
      if (T_minus[p, q] < 0) T_minus[p, q] <- T_minus[q, p] <- 0
    }
  }

  r_set <- seq_len(K_minus + 1L)

  # (A) directional
  lp_dir_exist <- vapply(seq_len(K_minus), function(kc)
    dir_exact_SST_existing_counts(k = kc, c_plus = c_plus, N_tot = N_tot, psi_vec = psi_minus),
    numeric(1))
  lp_dir_new <- rep(-Inf, K_minus + 1L)
  lp_dir_new[r_set] <- dir_exact_SST_new_vec_counts(
    c_plus = c_plus, N_tot = N_tot, psi_vec = psi_minus,
    r_set = r_set, tau0 = tau0)

  # (B) collapsed kappa
  lp_kappa_exist <- numeric(K_minus)
  for (kc in seq_len(K_minus)) {
    acc <- 0
    for (ell in seq_len(K_minus)) {
      if (r_add[ell] == 0 && t_add[ell] == 0) next
      p <- min(kc, ell); q <- max(kc, ell)
      acc <- acc + (gp_marginal(R_minus[p,q] + r_add[ell], T_minus[p,q] + t_add[ell], a_kappa, b_kappa) -
                      gp_marginal(R_minus[p,q], T_minus[p,q], a_kappa, b_kappa))
    }
    lp_kappa_exist[kc] <- acc
  }
  lp_kappa_new <- sum(gp_marginal(r_add, t_add, a_kappa, b_kappa))

  # (C) GN prior
  v_minus <- n_minus_full[keep_full]
  pw <- gn_log_weights_packed(v_minus, gamma_gn = gamma_gn)
  lp_prior_exist <- pw$exist
  lp_prior_new   <- pw$new

  # totals
  total_exist <- lp_dir_exist + lp_kappa_exist + lp_prior_exist
  total_new   <- lp_dir_new   + lp_kappa_new   + lp_prior_new

  logW_join <- lse(total_exist)
  logW_new  <- lse(total_new[r_set] - log(length(r_set)))

  p_birth <- exp(logW_new - lse(c(logW_new, logW_join)))

  # best existing
  best_k <- which.max(total_exist)

  list(
    i = i, K_minus = K_minus, oldk = oldk,
    p_birth = p_birth,
    # best existing block scores
    dir_best_exist = lp_dir_exist[best_k],
    kappa_best_exist = lp_kappa_exist[best_k],
    prior_best_exist = lp_prior_exist[best_k],
    total_best_exist = total_exist[best_k],
    # new block scores (best slot)
    best_new_r = which.max(total_new),
    dir_best_new = max(lp_dir_new[r_set]),
    kappa_new = lp_kappa_new,
    prior_new = lp_prior_new,
    total_best_new = max(total_new[r_set]),
    # means for comparing
    dir_new_mean = mean(lp_dir_new[r_set]),
    # block sizes
    v_minus = v_minus,
    N_tot = N_tot, r_add = r_add, t_add = t_add
  )
}


# --- Run on a real dataset ---
run_decomposition_test <- function(dataset_name, A, K_init, n_warmup, gamma_gn,
                                    a_kappa, b_kappa, a_eta, b_eta,
                                    tau0 = 0.15, mu0 = 0, sigma0 = 0.5) {
  cat(sprintf("\n--- Dataset: %s  (n=%d, K_init=%d, gamma=%.2f, a_kappa=%.1f, b_kappa=%.1f) ---\n",
              dataset_name, nrow(A), K_init, gamma_gn, a_kappa, b_kappa))

  n <- nrow(A)

  # Short warmup fit to get a reasonable state
  out <- modular_osbm_sampler(
    A = A, K = K_init,
    n_iter = n_warmup, burn = 1, thin = 1,
    verbose = FALSE,
    psi_constraint = "SST" ,
    partition_prior = "GN",
    gamma_gn = gamma_gn,
    a_kappa = a_kappa, b_kappa = b_kappa,
    a_eta = a_eta, b_eta = b_eta,
    mu0 = mu0, sigma0 = sigma0, tau0 = tau0,
    seed = 1
  )

  cat(sprintf("  K_trace after %d iters: %s\n", n_warmup,
              paste(out$K_trace, collapse=",")))

  # Take the LAST draw as our state
  S <- length(out$z)
  z <- out$z[[S]]
  kappa <- out$kappa[[S]]
  psi <- out$psi[[S]]
  eta <- out$eta[[S]]
  K <- max(z)

  cat(sprintf("  Final state: K=%d, block_sizes=%s\n", K, paste(tabulate(z, K), collapse=",")))
  cat(sprintf("  psi (length %d): %s\n", length(psi),
              paste(round(psi, 3), collapse=", ")))

  # Build edge list
  N_mat <- A + t(A); diag(N_mat) <- 0
  idx   <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  i_idx <- idx[,1]; j_idx <- idx[,2]
  N_edge <- as.numeric(N_mat[idx])
  edge_by_node <- replicate(n, integer(0), simplify = FALSE)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], as.integer(e))
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], as.integer(e))
  }

  bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K)

  # Instrument a sample of nodes
  sample_nodes <- sample.int(n, min(n, 20))
  results <- lapply(sample_nodes, function(i) {
    tryCatch(
      sst_update_i_instrumented(
        i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
        Rkl = bt$Rkl, Tkl = bt$Tkl,
        i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
        edge_by_node = edge_by_node,
        a_kappa = a_kappa, b_kappa = b_kappa,
        gamma_gn = gamma_gn, tau0 = tau0, mu0 = mu0, sig2_0 = sigma0
      ),
      error = function(e) {
        cat(sprintf("  ERROR at node %d: %s\n", i, conditionMessage(e)))
        NULL
      }
    )
  })
  results <- Filter(Nonnull <- function(x) !is.null(x), results)

  if (length(results) == 0) {
    cat("  No nodes could be evaluated!\n")
    return(invisible(NULL))
  }

  # Tabulate
  df <- data.frame(
    node       = sapply(results, `[[`, "i"),
    K_minus    = sapply(results, `[[`, "K_minus"),
    p_birth    = sapply(results, `[[`, "p_birth"),
    dir_exist  = sapply(results, `[[`, "dir_best_exist"),
    dir_new    = sapply(results, `[[`, "dir_best_new"),
    dir_delta  = sapply(results, function(r) r$dir_best_new - r$dir_best_exist),
    kap_exist  = sapply(results, `[[`, "kappa_best_exist"),
    kap_new    = sapply(results, `[[`, "kappa_new"),
    kap_delta  = sapply(results, function(r) r$kappa_new - r$kappa_best_exist),
    prior_exist= sapply(results, `[[`, "prior_best_exist"),
    prior_new  = sapply(results, `[[`, "prior_new"),
    prior_delta= sapply(results, function(r) r$prior_new - r$prior_best_exist)
  )

  cat("\n  Per-node birth decomposition (new - best_exist for each component):\n")
  cat(sprintf("  %-5s  %6s  %8s |  %8s  %8s  %8s |  %8s  %8s  %8s |  %8s  %8s  %8s\n",
              "node", "K-", "P(birth)",
              "dir_ex", "dir_new", "dir_d",
              "kap_ex", "kap_new", "kap_d",
              "pri_ex", "pri_new", "pri_d"))
  for (r in seq_len(nrow(df))) {
    cat(sprintf("  %-5d  %6d  %8.4f |  %8.2f  %8.2f  %8.2f |  %8.2f  %8.2f  %8.2f |  %8.2f  %8.2f  %8.2f\n",
                df$node[r], df$K_minus[r], df$p_birth[r],
                df$dir_exist[r], df$dir_new[r], df$dir_delta[r],
                df$kap_exist[r], df$kap_new[r], df$kap_delta[r],
                df$prior_exist[r], df$prior_new[r], df$prior_delta[r]))
  }

  cat(sprintf("\n  SUMMARY over %d nodes:\n", nrow(df)))
  cat(sprintf("    Mean P(birth)       = %.4f\n", mean(df$p_birth)))
  cat(sprintf("    Mean dir  delta     = %.2f  (new-exist: + favors birth)\n", mean(df$dir_delta)))
  cat(sprintf("    Mean kappa delta    = %.2f  (new-exist: + favors birth)\n", mean(df$kap_delta)))
  cat(sprintf("    Mean prior delta    = %.2f  (new-exist: + favors birth)\n", mean(df$prior_delta)))
  cat(sprintf("    Births mostly driven by: %s\n",
              if (abs(mean(df$kap_delta)) > abs(mean(df$dir_delta)) &&
                  abs(mean(df$kap_delta)) > abs(mean(df$prior_delta))) "KAPPA (volume)"
              else if (abs(mean(df$dir_delta)) > abs(mean(df$prior_delta))) "DIRECTION (psi)"
              else "PRIOR (GN)"))

  invisible(df)
}


# ============================================================================
# TEST 5: Apply to real data
# ============================================================================
cat("\n============================================================\n")
cat("TEST 5: Real-data birth decomposition\n")
cat("============================================================\n")

# --- Mountain goats ---
source("scripts/application/application.R", local = TRUE, chdir = FALSE)
mg <- tryCatch(choose_dataset("mountain_goats"), error = function(e) NULL)
if (!is.null(mg)) {
  cat("\n========== MOUNTAIN GOATS ==========\n")

  # Baseline hyperparams
  df1 <- run_decomposition_test("mountain_goats (baseline)",
    A = mg$A, K_init = 3, n_warmup = 50,
    gamma_gn = 0.9,
    a_kappa = 50, b_kappa = 0.1,
    a_eta = 0.1, b_eta = 1)

  # Higher b_kappa (user says this helps)
  df2 <- run_decomposition_test("mountain_goats (b_kappa=10)",
    A = mg$A, K_init = 3, n_warmup = 50,
    gamma_gn = 0.9,
    a_kappa = 50, b_kappa = 10,
    a_eta = 0.1, b_eta = 1)

  # Very high gamma
  df3 <- run_decomposition_test("mountain_goats (gamma=0.99)",
    A = mg$A, K_init = 3, n_warmup = 50,
    gamma_gn = 0.99,
    a_kappa = 50, b_kappa = 0.1,
    a_eta = 0.1, b_eta = 1)
}

# --- High school ---
hs <- tryCatch(choose_dataset("high_school"), error = function(e) NULL)
if (!is.null(hs)) {
  cat("\n========== HIGH SCHOOL ==========\n")
  df4 <- run_decomposition_test("high_school (baseline)",
    A = hs$A, K_init = 3, n_warmup = 50,
    gamma_gn = 0.9,
    a_kappa = 50, b_kappa = 0.1,
    a_eta = 0.1, b_eta = 1)
}

# --- Macaques ---
mc <- tryCatch(choose_dataset("macaques_data"), error = function(e) NULL)
if (!is.null(mc)) {
  cat("\n========== MACAQUES ==========\n")
  df5 <- run_decomposition_test("macaques (baseline)",
    A = mc$A, K_init = 3, n_warmup = 50,
    gamma_gn = 0.9,
    a_kappa = 50, b_kappa = 0.1,
    a_eta = 0.1, b_eta = 1)
}


# ============================================================================
# TEST 6: How many K_minus+1 new slots vs K_minus existing — counting asymmetry
# ============================================================================
cat("\n============================================================\n")
cat("TEST 6: Slot-counting asymmetry in birth/join decision\n")
cat("============================================================\n")
cat("\nThe birth/join decision uses:\n")
cat("  logW_join = lse(total_exist)                        — K_minus terms\n")
cat("  logW_new  = lse(total_new[r_set] - log(|r_set|))   — K_minus+1 terms\n")
cat("\nThe new-block score is the LSE over K+1 slots, each getting the SAME kappa_new.\n")
cat("So the new-block log-weight gets a 'bonus' of log(K+1) from the slot averaging.\n")
cat("Meanwhile, each new slot also gets the SAME prior_new.\n")
cat("So the effective new-block weight is:\n")
cat("  lse(dir_new[r] + kappa_new + prior_new for r in 1..K+1) - log(K+1)\n")
cat("  = lse(dir_new[r]) + kappa_new + prior_new - log(K+1)\n")
cat("This means the directional scores over different slots get to 'pool'\n")
cat("while the kappa_new and prior_new are constants added to ALL slots.\n\n")

cat("KEY INSIGHT:\n")
cat("  lp_kappa_new is a CONSTANT added to all K+1 slots, so it enters\n")
cat("  the LSE as: lse(dir_new[r]) + kappa_new + prior_new - log(K+1)\n")
cat("  vs existing: max(dir_exist[k] + kappa_exist[k] + prior_exist[k])\n\n")
cat("  If kappa_new >> max(kappa_exist[k]), this systematically favors births.\n\n")


cat("\n========== ALL TESTS COMPLETE ==========\n")
