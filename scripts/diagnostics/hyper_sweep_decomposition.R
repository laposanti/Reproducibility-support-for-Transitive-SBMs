#!/usr/bin/env Rscript
# ============================================================================
# hyper_sweep_decomposition.R
#
# Systematic sweep: run SST sampler under different prior hyperparameters
# and different datasets, then decompose the birth/join decision into
#    prior (GN)  +  kappa (volume)  +  direction (psi)
# for every node at a snapshot state.
#
# Strategy:
#   1.  Run short warmup (n_warmup iterations) under each hyperparameter set
#   2.  At the final state, instrument a z-sweep over all nodes
#   3.  Record per-node: prior/kappa/direction contributions to birth vs best-exist
#   4.  Aggregate into a tidy CSV for analysis
# ============================================================================

cat("====================================================================\n")
cat("  HYPERPARAMETER SWEEP — birth decomposition\n")
cat("====================================================================\n\n")

suppressPackageStartupMessages({
  library(Rcpp)
  source("helper_folder/helper.R")
  source("helper_folder/SST_helpers.R")
  source("helper_folder/WST_helpers.R")
  source("helper_folder/Hyper_setup.R")
  source("core/my_best_try_so_far.R")
})

# ---- Dataset loader (extracted from application.R, no side effects) ----
load_dataset <- function(dataset) {
  if (dataset == "mountain_goats") {
    fs <- list.files("./data/ShizukaMcDonald_Data", full.names = TRUE, pattern = "\\.csv$")
    n_each <- vapply(fs, function(f) nrow(read.csv(f, row.names = 1)), integer(1))
    A <- as.matrix(read.csv(fs[which.max(n_each)], row.names = 1, check.names = FALSE))
  } else if (dataset == "citations_data") {
    A <- as.matrix(read.csv("./data/Citations_application/cross-citation-matrix.csv",
                            row.names = 1, header = TRUE, check.names = FALSE))
    diag(A) <- 0
  } else if (dataset == "macaques_data") {
    el <- read.table("./data/macaques/out.moreno.txt")
    nodes <- sort(unique(c(el[[1]], el[[2]])))
    A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (i in seq_len(nrow(el))) A[el[i,1], el[i,2]] <- el[i, "V3"]
  } else if (dataset == "high_school") {
    edges <- read.csv("./data/high-school/edges.csv", header = FALSE,
                      comment.char = "#", strip.white = TRUE)
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    nn <- max(c(edges$source, edges$target))
    A <- matrix(0L, nn, nn,
                dimnames = list(as.character(seq_len(nn)), as.character(seq_len(nn))))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]; j <- edges$target[r]; w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w
    }
    diag(A) <- 0L
  } else stop("Unknown dataset: ", dataset)
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  A
}

# ----- Instrumented per-node decomposition (from diagnose_sst_births.R) -----
sst_decompose_node <- function(
    i, A, z, eta, kappa, psi,
    Rkl, Tkl,
    i_idx, j_idx, N_edge, edge_by_node,
    a_kappa, b_kappa, gamma_gn,
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
  z_packed  <- map[z]
  oldk_pos  <- match(oldk, keep_full)
  kappa_minus <- kappa[keep_full, keep_full, drop = FALSE]
  psi_minus   <- reindex_psi_sst_keep(psi_old = psi, K_old = K_full, keep = keep_full)

  C_i_full <- counts_by_block_exact_cpp(
    i = i, A = A, z = as.integer(z),
    i_idx = as.integer(i_idx), j_idx = as.integer(j_idx),
    N_edge = as.numeric(N_edge),
    edge_by_node = edge_by_node, K = K_full)
  c_plus <- as.numeric(C_i_full$c_plus)[keep_full]
  N_tot  <- as.numeric(C_i_full$N_tot )[keep_full]

  E_excl_full <- E_by_block_excluding_i(i, z, eta, K_full)
  E_excl <- as.numeric(E_excl_full[keep_full])
  r_add  <- N_tot
  t_add  <- eta_i * E_excl

  R_minus <- Rkl[keep_full, keep_full, drop = FALSE]
  T_minus <- Tkl[keep_full, keep_full, drop = FALSE]
  if (!is.na(oldk_pos)) {
    for (ell_pos in seq_len(K_minus)) {
      subR <- N_tot[ell_pos]; subT <- eta_i * E_excl[ell_pos]
      p <- min(oldk_pos, ell_pos); q <- max(oldk_pos, ell_pos)
      R_minus[p,q] <- R_minus[p,q] - subR
      T_minus[p,q] <- T_minus[p,q] - subT
      R_minus[q,p] <- R_minus[p,q]; T_minus[q,p] <- T_minus[p,q]
      if (R_minus[p,q] < 0) R_minus[p,q] <- R_minus[q,p] <- 0
      if (T_minus[p,q] < 0) T_minus[p,q] <- T_minus[q,p] <- 0
    }
  }

  r_set <- seq_len(K_minus + 1L)

  # (A) direction -------------------------------------------------------
  lp_dir_exist <- vapply(seq_len(K_minus), function(kc)
    dir_exact_SST_existing_counts(k = kc, c_plus = c_plus, N_tot = N_tot,
                                  psi_vec = psi_minus), numeric(1))
  lp_dir_new <- rep(-Inf, K_minus + 1L)
  lp_dir_new[r_set] <- dir_exact_SST_new_vec_counts(
    c_plus = c_plus, N_tot = N_tot, psi_vec = psi_minus,
    r_set = r_set, tau0 = tau0)

  # (B) collapsed kappa -------------------------------------------------
  lp_kappa_exist <- numeric(K_minus)
  for (kc in seq_len(K_minus)) {
    acc <- 0
    for (ell in seq_len(K_minus)) {
      if (r_add[ell] == 0 && t_add[ell] == 0) next
      p <- min(kc, ell); q <- max(kc, ell)
      acc <- acc + (gp_marginal(R_minus[p,q] + r_add[ell],
                                T_minus[p,q] + t_add[ell], a_kappa, b_kappa) -
                      gp_marginal(R_minus[p,q], T_minus[p,q], a_kappa, b_kappa))
    }
    lp_kappa_exist[kc] <- acc
  }
  lp_kappa_new <- sum(gp_marginal(r_add, t_add, a_kappa, b_kappa))

  # (C) GN prior --------------------------------------------------------
  v_minus <- n_minus_full[keep_full]
  pw <- gn_log_weights_packed(v_minus, gamma_gn = gamma_gn)
  lp_prior_exist <- pw$exist
  lp_prior_new   <- pw$new

  # (D) aggregate -------------------------------------------------------
  total_exist <- lp_dir_exist + lp_kappa_exist + lp_prior_exist
  total_new   <- lp_dir_new   + lp_kappa_new   + lp_prior_new

  logW_join <- lse(total_exist)
  logW_new  <- lse(total_new[r_set] - log(length(r_set)))
  p_birth   <- exp(logW_new - lse(c(logW_new, logW_join)))

  best_k    <- which.max(total_exist)
  best_r    <- r_set[which.max(total_new[r_set])]

  list(
    i = i, K_minus = K_minus, oldk = oldk, p_birth = p_birth,
    dir_exist   = lp_dir_exist[best_k],
    dir_new     = lp_dir_new[best_r],
    kappa_exist = lp_kappa_exist[best_k],
    kappa_new   = lp_kappa_new,
    prior_exist = lp_prior_exist[best_k],
    prior_new   = lp_prior_new
  )
}


# ============================================================================
# Core: run one (dataset x hypers) configuration
# ============================================================================
run_one_config <- function(
    dataset_name, A, K_init, n_warmup, seed,
    a_kappa, b_kappa, a_eta, b_eta,
    gamma_gn, tau0, mu0 = 0, sigma0 = 0.5,
    tag = ""
) {
  n <- nrow(A)
  label <- sprintf("%s  K0=%d  a_k=%.2f b_k=%.2f  a_e=%.2f b_e=%.2f  gam=%.2f tau0=%.3f%s",
                   dataset_name, K_init, a_kappa, b_kappa, a_eta, b_eta,
                   gamma_gn, tau0, if (nzchar(tag)) paste0("  [", tag, "]") else "")
  cat(sprintf("\n--- %s ---\n", label))

  out <- tryCatch(
    modular_osbm_sampler(
      A = A, K = K_init,
      n_iter = n_warmup, burn = 1, thin = 1,
      verbose = FALSE,
      psi_constraint = "SST" ,
      partition_prior = "GN",
      gamma_gn = gamma_gn,
      a_kappa = a_kappa, b_kappa = b_kappa,
      a_eta = a_eta, b_eta = b_eta,
      mu0 = mu0, sigma0 = sigma0, tau0 = tau0,
      seed = seed
    ),
    error = function(e) {
      cat("  SAMPLER ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(out)) return(NULL)

  # K trajectory
  K_trace <- out$K_trace
  cat(sprintf("  K_trace: %s\n", paste(K_trace, collapse=",")))

  # Grab final state
  S <- length(out$z)
  z <- out$z[[S]]; kappa <- out$kappa[[S]]
  psi <- out$psi[[S]]; eta <- out$eta[[S]]
  K_final <- max(z)
  blk_sizes <- tabulate(z, K_final)
  cat(sprintf("  K_final=%d  sizes=%s\n", K_final, paste(blk_sizes, collapse=",")))

  # Edge list
  N_mat <- A + t(A); diag(N_mat) <- 0
  idx   <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  i_idx <- idx[,1]; j_idx <- idx[,2]; N_edge <- as.numeric(N_mat[idx])
  edge_by_node <- replicate(n, integer(0), simplify = FALSE)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], as.integer(e))
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], as.integer(e))
  }
  bt <- block_totals_for_poisson_cpp(z, eta, i_idx, j_idx, N_edge, K_final)

  # Instrument ALL nodes
  results <- lapply(seq_len(n), function(i) {
    tryCatch(
      sst_decompose_node(
        i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
        Rkl = bt$Rkl, Tkl = bt$Tkl,
        i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
        edge_by_node = edge_by_node,
        a_kappa = a_kappa, b_kappa = b_kappa,
        gamma_gn = gamma_gn, tau0 = tau0, mu0 = mu0, sig2_0 = sigma0),
      error = function(e) NULL)
  })
  results <- Filter(Negate(is.null), results)
  if (length(results) == 0L) { cat("  No nodes evaluated!\n"); return(NULL) }

  # Build data.frame
  df <- data.frame(
    dataset   = dataset_name,
    tag       = tag,
    n         = n,
    K_init    = K_init,
    K_final   = K_final,
    a_kappa   = a_kappa,
    b_kappa   = b_kappa,
    a_eta     = a_eta,
    b_eta     = b_eta,
    gamma_gn  = gamma_gn,
    tau0      = tau0,
    node      = sapply(results, `[[`, "i"),
    K_minus   = sapply(results, `[[`, "K_minus"),
    p_birth   = sapply(results, `[[`, "p_birth"),
    dir_exist = sapply(results, `[[`, "dir_exist"),
    dir_new   = sapply(results, `[[`, "dir_new"),
    kap_exist = sapply(results, `[[`, "kappa_exist"),
    kap_new   = sapply(results, `[[`, "kappa_new"),
    pri_exist = sapply(results, `[[`, "prior_exist"),
    pri_new   = sapply(results, `[[`, "prior_new"),
    stringsAsFactors = FALSE
  )
  df$dir_delta   <- df$dir_new   - df$dir_exist
  df$kap_delta   <- df$kap_new   - df$kap_exist
  df$pri_delta   <- df$pri_new   - df$pri_exist
  df$total_delta <- df$dir_delta + df$kap_delta + df$pri_delta

  # Summary line
  cat(sprintf("  mean P(birth)=%.4f  dir_d=%+.2f  kap_d=%+.2f  pri_d=%+.2f  total_d=%+.2f\n",
              mean(df$p_birth), mean(df$dir_delta), mean(df$kap_delta),
              mean(df$pri_delta), mean(df$total_delta)))
  df
}


# ============================================================================
# Load datasets
# ============================================================================
cat("Loading datasets...\n")
datasets <- list()
for (ds in c("mountain_goats", "macaques_data", "high_school", "citations_data")) {
  datasets[[ds]] <- tryCatch(load_dataset(ds), error = function(e) {
    cat(sprintf("  SKIP %s: %s\n", ds, conditionMessage(e))); NULL
  })
}
for (ds in names(datasets)) {
  if (!is.null(datasets[[ds]]))
    cat(sprintf("  %-18s n=%d  density=%.3f  mean_N=%.2f\n",
                ds, nrow(datasets[[ds]]),
                mean(datasets[[ds]] > 0),
                compute_mean_edges(datasets[[ds]])))
}

# ============================================================================
# Define sweep grid
# ============================================================================
# We vary one thing at a time from a common "baseline" to isolate effects.

baseline <- list(K_init = 3, n_warmup = 60, seed = 42,
                 a_kappa = 1, b_kappa = 1,
                 a_eta = 1, b_eta = 1,
                 gamma_gn = 0.9, tau0 = 0.15)

sweep_configs <- list(
  # --- Baseline ---
  list(tag = "baseline",
       a_kappa = 1, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 0.15),

  # --- Vary a_kappa ---
  list(tag = "a_kappa=0.1",
       a_kappa = 0.1, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 0.15),
  list(tag = "a_kappa=10",
       a_kappa = 10, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 0.15),
  list(tag = "a_kappa=50",
       a_kappa = 50, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 0.15),

  # --- Vary b_kappa ---
  list(tag = "b_kappa=0.1",
       a_kappa = 1, b_kappa = 0.1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 0.15),
  list(tag = "b_kappa=10",
       a_kappa = 1, b_kappa = 10, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 0.15),
  list(tag = "b_kappa=50",
       a_kappa = 1, b_kappa = 50, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 0.15),

  # --- Vary gamma_gn ---
  list(tag = "gamma=0.3",
       a_kappa = 1, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.3, tau0 = 0.15),
  list(tag = "gamma=0.5",
       a_kappa = 1, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.5, tau0 = 0.15),
  list(tag = "gamma=0.99",
       a_kappa = 1, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.99, tau0 = 0.15),

  # --- Vary tau0 ---
  list(tag = "tau0=0.01",
       a_kappa = 1, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 0.01),
  list(tag = "tau0=1.0",
       a_kappa = 1, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 1.0),
  list(tag = "tau0=5.0",
       a_kappa = 1, b_kappa = 1, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 5.0),

  # --- Vary a_eta, b_eta ---
  list(tag = "a_eta=0.1,b_eta=0.1",
       a_kappa = 1, b_kappa = 1, a_eta = 0.1, b_eta = 0.1,
       gamma_gn = 0.9, tau0 = 0.15),
  list(tag = "a_eta=10,b_eta=10",
       a_kappa = 1, b_kappa = 1, a_eta = 10, b_eta = 10,
       gamma_gn = 0.9, tau0 = 0.15),

  # --- Joint: principled-style a_kappa,b_kappa combos ---
  list(tag = "strong_kappa(a=50,b=10)",
       a_kappa = 50, b_kappa = 10, a_eta = 1, b_eta = 1,
       gamma_gn = 0.9, tau0 = 0.15),
  list(tag = "informative(a=5,b=5)",
       a_kappa = 5, b_kappa = 5, a_eta = 2, b_eta = 2,
       gamma_gn = 0.9, tau0 = 0.15)
)

# ============================================================================
# Run sweep
# ============================================================================
all_results <- list()
cfg_id <- 0L

for (ds_name in names(datasets)) {
  A <- datasets[[ds_name]]
  if (is.null(A)) next

  for (cfg in sweep_configs) {
    cfg_id <- cfg_id + 1L
    df <- tryCatch(
      run_one_config(
        dataset_name = ds_name, A = A,
        K_init  = baseline$K_init,
        n_warmup = baseline$n_warmup,
        seed    = baseline$seed,
        a_kappa = cfg$a_kappa,
        b_kappa = cfg$b_kappa,
        a_eta   = cfg$a_eta,
        b_eta   = cfg$b_eta,
        gamma_gn = cfg$gamma_gn,
        tau0    = cfg$tau0,
        tag     = cfg$tag),
      error = function(e) {
        cat(sprintf("  CONFIG ERROR: %s\n", conditionMessage(e)))
        NULL
      })
    if (!is.null(df)) all_results[[cfg_id]] <- df
  }
}

# ============================================================================
# Combine and save
# ============================================================================
if (length(all_results) > 0) {
  full_df <- do.call(rbind, all_results)
  rownames(full_df) <- NULL

  # Save detailed per-node results
  write.csv(full_df, "hyper_sweep_per_node.csv", row.names = FALSE)
  cat(sprintf("\nSaved %d rows to hyper_sweep_per_node.csv\n", nrow(full_df)))

  # ---- Aggregated summary ----
  agg <- aggregate(
    cbind(p_birth, dir_delta, kap_delta, pri_delta, total_delta, K_minus) ~ 
      dataset + tag + n + K_init + K_final + a_kappa + b_kappa + a_eta + b_eta + gamma_gn + tau0,
    data = full_df, FUN = mean)
  names(agg)[names(agg) == "K_minus"] <- "mean_K_minus"

  # dominant driver
  agg$driver <- ifelse(
    abs(agg$kap_delta) >= abs(agg$dir_delta) & abs(agg$kap_delta) >= abs(agg$pri_delta),
    "kappa",
    ifelse(abs(agg$dir_delta) >= abs(agg$pri_delta), "direction", "prior"))

  # sort by dataset, then by mean total_delta (most birth-prone first)
  agg <- agg[order(agg$dataset, -agg$total_delta), ]

  write.csv(agg, "hyper_sweep_summary.csv", row.names = FALSE)
  cat(sprintf("Saved %d configs to hyper_sweep_summary.csv\n\n", nrow(agg)))

  # ---- Pretty-print ----
  cat("====================================================================\n")
  cat("                   AGGREGATED SUMMARY TABLE\n")
  cat("====================================================================\n\n")
  cat(sprintf("%-18s %-28s  K0->Kf  %8s  %8s  %8s  %8s  %8s  %-9s\n",
              "dataset", "config", "P(brth)", "dir_d", "kap_d", "pri_d", "tot_d", "driver"))
  cat(strrep("-", 130), "\n")

  for (r in seq_len(nrow(agg))) {
    cat(sprintf("%-18s %-28s  %d->%-3d  %8.4f  %+8.2f  %+8.2f  %+8.2f  %+8.2f  %-9s\n",
                agg$dataset[r], agg$tag[r],
                agg$K_init[r], agg$K_final[r],
                agg$p_birth[r],
                agg$dir_delta[r], agg$kap_delta[r],
                agg$pri_delta[r], agg$total_delta[r],
                agg$driver[r]))
  }

  # ---- Cross-dataset insight: which hyperparameter matters most? ----
  cat("\n====================================================================\n")
  cat("  WHICH HYPERPARAMETER HAS THE LARGEST EFFECT ON P(birth)?\n")
  cat("====================================================================\n\n")

  # For each dataset, compare the range of P(birth) across configs
  for (ds in unique(agg$dataset)) {
    sub <- agg[agg$dataset == ds, ]
    cat(sprintf("  %s (n=%d):\n", ds, sub$n[1]))
    cat(sprintf("    P(birth) range: [%.4f, %.4f]\n", min(sub$p_birth), max(sub$p_birth)))
    cat(sprintf("    K_final range:  [%d, %d]\n", min(sub$K_final), max(sub$K_final)))

    # Which varied parameter yields the widest P(birth) swing?
    # Group by param varied
    bp <- sub$p_birth; names(bp) <- sub$tag
    baseline_pb <- bp["baseline"]

    # Sensitivity = |P(birth) - baseline| for each config
    sens <- abs(bp - baseline_pb)
    top3 <- head(sort(sens, decreasing = TRUE), 5)
    cat("    Top-5 configs by |delta P(birth)| from baseline:\n")
    for (i in seq_along(top3)) {
      idx <- which(sub$tag == names(top3)[i])
      cat(sprintf("      %s: P(birth)=%.4f  K=%d  (kap_d=%+.1f  dir_d=%+.1f  pri_d=%+.1f)\n",
                  names(top3)[i], sub$p_birth[idx], sub$K_final[idx],
                  sub$kap_delta[idx], sub$dir_delta[idx], sub$pri_delta[idx]))
    }
    cat("\n")
  }

} else {
  cat("No results collected — check for errors above.\n")
}

cat("\n========== SWEEP COMPLETE ==========\n")
