#!/usr/bin/env Rscript
# ======================================================================
# Diagnostic: hyperparameter sensitivity to network density
#
# This script investigates WHY the OSBM birth/death dynamics are
# sensitive to network density and evaluates principled fixes.
#
# Three analyses:
#   1) Prior predictive analysis (analytical — no MCMC)
#   2) Birth-move Bayes factor decomposition
#   3) Short MCMC comparison: current rule vs sample_b_kappa
# ======================================================================
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

suppressPackageStartupMessages({
  library(dplyr)
})

source("./helper_folder/helper.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./helper_folder/Hyper_setup.R")
source("./helper_folder/sim_study_helper.R")
source("./core/ppc_checks.R")
source("./core/my_best_try_so_far.R")
if (!exists("log1pexp")) log1pexp <- function(x) ifelse(x < 35, log1p(exp(x)), x)

# === Re-use the same data loader from the application script ===
DATASET_CHOICES <- c(
  "moreno_sheep", "strauss_2019b", "mountain_goats",
  "citations_data", "macaques_data", "high_school", "hiv1_data"
)

choose_dataset <- function(dataset) {
  dataset <- match.arg(dataset, choices = DATASET_CHOICES)
  if (dataset == "mountain_goats") {
    matrix_files <- list.files("./data/ShizukaMcDonald_Data",
                               full.names = TRUE, pattern = "[.]csv$")
    n_each <- vapply(matrix_files, function(f) nrow(read.csv(f, row.names = 1)),
                     FUN.VALUE = integer(1))
    A <- as.matrix(read.csv(matrix_files[which.max(n_each)],
                            row.names = 1, check.names = FALSE))
  } else if (dataset == "citations_data") {
    A <- read.csv("./data/Citations_application/cross-citation-matrix.csv",
                  row.names = 1, header = TRUE, check.names = FALSE)
    diag(A) <- 0
  } else if (dataset == "macaques_data") {
    edge_list <- read.table("./data/macaques/out.moreno.txt")
    nodes <- sort(unique(c(edge_list[[1]], edge_list[[2]])))
    A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (i in seq_len(nrow(edge_list)))
      A[edge_list[i, 1], edge_list[i, 2]] <- edge_list[i, "V3"]
  } else if (dataset == "high_school") {
    edges <- read.csv("./data/high-school/edges.csv",
                      header = FALSE, comment.char = "#", strip.white = TRUE)
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    n_nodes <- max(c(edges$source, edges$target))
    A <- matrix(0L, n_nodes, n_nodes,
                dimnames = list(seq_len(n_nodes), seq_len(n_nodes)))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]; j <- edges$target[r]; w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w
    }
    diag(A) <- 0L
  } else if (dataset == "hiv1_data") {
    A <- read.csv("./data/HIV1_collapsed_adjacency_matrix.csv",
                  row.names = 1, header = TRUE, check.names = FALSE)
    diag(A) <- 0
  } else if (dataset == "moreno_sheep") {
    edges <- read.csv("./data/moreno_sheep/edges.csv",
                      comment.char = "#", strip.white = TRUE)
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    n_nodes <- max(c(edges$source, edges$target))
    A <- matrix(0L, n_nodes, n_nodes,
                dimnames = list(seq_len(n_nodes), seq_len(n_nodes)))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]; j <- edges$target[r]; w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w
    }
    diag(A) <- 0L
  } else if (dataset == "strauss_2019b") {
    edges <- read.csv("./data/Strauss_2019b/edges.csv",
                      comment.char = "#", strip.white = TRUE)
    names(edges)[1:2] <- c("source", "target")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- 1L
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    n_nodes <- max(c(edges$source, edges$target))
    A <- matrix(0L, n_nodes, n_nodes,
                dimnames = list(seq_len(n_nodes), seq_len(n_nodes)))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]; j <- edges$target[r]; w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w
    }
    diag(A) <- 0L
  }
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  colnames(A) <- rownames(A)
  A
}

# ======================================================================
# PART 1: Prior predictive analysis (pure analytical)
#
# Under the prior kappa ~ Gamma(a, b) with K equal blocks:
#   - Typical exposure per block pair: T_typ = (n/K)^2
#   - Prior predictive for block-pair total: R ~ NegBin(a, b/(b+T))
#   - Prior predictive mean: a*T/b
#   - Prior predictive expected total edges: P * a * T / b
#
# We compare this to the observed total for each dataset.
# ======================================================================
cat("\n")
cat("================================================================\n")
cat("PART 1: PRIOR PREDICTIVE ANALYSIS\n")
cat("================================================================\n\n")

# Define hyperparameter settings to compare
hyper_settings <- list(
  "current_rule"  = function(n, K, mean_deg) list(a = 1, b = 1.0 * mean_deg / K),
  "a=1, b=1"      = function(n, K, mean_deg) list(a = 1, b = 1),
  "a=0.5, b=0.5"  = function(n, K, mean_deg) list(a = 0.5, b = 0.5),
  "a=1, b=0.01"   = function(n, K, mean_deg) list(a = 1, b = 0.01)
)

pp_rows <- list()
for (ds in DATASET_CHOICES) {
  A <- tryCatch(choose_dataset(ds), error = function(e) NULL)
  if (is.null(A)) next

  n <- nrow(A)
  total_edges <- sum(A)
  mean_deg <- total_edges / n
  density <- mean_deg / (n - 1)
  K_sp <- estimate_K_spectral(A)
  deg <- rowSums(A) + colSums(A)
  cv_deg <- sd(deg) / mean(deg)

  P <- K_sp * (K_sp + 1) / 2          # number of block pairs
  T_typ <- (n / K_sp)^2               # typical exposure per pair

  for (nm in names(hyper_settings)) {
    h <- hyper_settings[[nm]](n, K_sp, mean_deg)
    a <- h$a; b <- h$b

    # Prior predictive expected total edges
    pp_mean <- P * a * T_typ / b
    pp_var  <- P * a * T_typ * (b + T_typ) / b^2
    pp_sd   <- sqrt(pp_var)

    # z-score: how surprising is the data under the prior?
    z_score <- (total_edges - pp_mean) / pp_sd

    # Shrinkage weight: prior contribution to posterior mean of kappa
    w_kappa <- b / (b + T_typ)

    # Prior mean of kappa
    prior_mean_kappa <- a / b

    # Observed mean kappa (naive: total_edges / total_exposure)
    obs_mean_kappa <- total_edges / (P * T_typ)

    pp_rows[[length(pp_rows) + 1]] <- data.frame(
      dataset = ds,
      setting = nm,
      n = n, K_sp = K_sp, density = round(density, 3),
      total_edges = total_edges,
      a_kappa = a, b_kappa = round(b, 2),
      prior_mean_kappa = round(prior_mean_kappa, 4),
      obs_kappa = round(obs_mean_kappa, 4),
      pp_expected_edges = round(pp_mean, 0),
      pp_z_score = round(z_score, 2),
      shrinkage_w = round(w_kappa, 4),
      stringsAsFactors = FALSE
    )
  }
}

pp_df <- do.call(rbind, pp_rows)
cat("Prior predictive comparison across datasets and hyper settings:\n\n")

# Print per-dataset tables
for (ds in unique(pp_df$dataset)) {
  sub <- pp_df[pp_df$dataset == ds, ]
  cat(sprintf("--- %s (n=%d, K_sp=%d, density=%.3f, total_edges=%d) ---\n",
              ds, sub$n[1], sub$K_sp[1], sub$density[1], sub$total_edges[1]))
  cat(sprintf("    Observed mean kappa = %.4f\n\n", sub$obs_kappa[1]))
  cat(sprintf("    %-18s  a_kap  b_kap  prior_mean  pp_edges  z_score  shrinkage\n", "Setting"))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("    %-18s  %5.2f  %5.2f  %9.4f  %8d  %7.2f  %8.4f\n",
                sub$setting[i], sub$a_kappa[i], sub$b_kappa[i],
                sub$prior_mean_kappa[i], sub$pp_expected_edges[i],
                sub$pp_z_score[i], sub$shrinkage_w[i]))
  }
  cat("\n")
}

# ======================================================================
# PART 2: Birth-move Bayes factor decomposition
#
# Scenario: K_current = 2, equal-sized blocks.
# A single node i is considered for birth (creating block 3).
# The birth marginal score for the kappa component is:
#   sum_l gp_marginal(R_il, T_il, a, b)
# where R_il = edges from i to block l, T_il = eta_i * E_l.
#
# We compute this for the "median" node (by degree).
# ======================================================================
cat("\n")
cat("================================================================\n")
cat("PART 2: BIRTH-MOVE KAPPA SCORE (per-node, analytical)\n")
cat("================================================================\n\n")

cat("This measures how much the Gamma-Poisson marginal for a NEW block\n")
cat("favors or penalizes birth, as a function of (a_kappa, b_kappa).\n")
cat("Higher = easier to create new block. Lower = harder.\n\n")

# gp_marginal for scalar inputs
gp_marg <- function(R, T, a, b) {
  if (T == 0 && R == 0) return(0)
  lgamma(a + R) - lgamma(a) + a * log(b) - (a + R) * log(b + T)
}

birth_rows <- list()
for (ds in DATASET_CHOICES) {
  A <- tryCatch(choose_dataset(ds), error = function(e) NULL)
  if (is.null(A)) next

  n <- nrow(A)
  total_edges <- sum(A)
  mean_deg <- total_edges / n
  K_sp <- estimate_K_spectral(A)

  # Simulate K=2 equal blocks; pick median-degree node
  deg_out <- rowSums(A)
  med_node <- order(deg_out)[ceiling(n / 2)]  # median by out-degree

  # Assume 2 blocks, evenly split. Node med_node in block 1.
  block_size <- floor(n / 2)
  block1 <- seq_len(block_size)
  block2 <- setdiff(seq_len(n), block1)

  # Edges from med_node to block2 nodes
  R_to_b2 <- sum(A[med_node, block2]) + sum(A[block2, med_node])
  T_to_b2 <- length(block2)  # assuming eta = 1

  # Edges from med_node to block1 (excluding self)
  b1_other <- setdiff(block1, med_node)
  R_to_b1 <- sum(A[med_node, b1_other]) + sum(A[b1_other, med_node])
  T_to_b1 <- length(b1_other)

  for (nm in names(hyper_settings)) {
    h <- hyper_settings[[nm]](n, K_sp, mean_deg)
    a <- h$a; b <- h$b

    # Birth kappa score: marginals for new block vs (block1, block2)
    score_b1 <- gp_marg(R_to_b1, T_to_b1, a, b)
    score_b2 <- gp_marg(R_to_b2, T_to_b2, a, b)
    # Also the self-block (new block has just this node, no self-edges)
    score_self <- gp_marg(0, 0, a, b)

    birth_kappa_score <- score_b1 + score_b2 + score_self

    birth_rows[[length(birth_rows) + 1]] <- data.frame(
      dataset = ds, setting = nm,
      R_b1 = R_to_b1, R_b2 = R_to_b2, T_b1 = T_to_b1, T_b2 = T_to_b2,
      a_kappa = a, b_kappa = round(b, 2),
      score_b1 = round(score_b1, 2),
      score_b2 = round(score_b2, 2),
      birth_total = round(birth_kappa_score, 2),
      stringsAsFactors = FALSE
    )
  }
}

birth_df <- do.call(rbind, birth_rows)
for (ds in unique(birth_df$dataset)) {
  sub <- birth_df[birth_df$dataset == ds, ]
  cat(sprintf("--- %s (R_b1=%d, T_b1=%d, R_b2=%d, T_b2=%d) ---\n",
              ds, sub$R_b1[1], sub$T_b1[1], sub$R_b2[1], sub$T_b2[1]))
  cat(sprintf("    %-18s  a_kap  b_kap  score_b1  score_b2  total\n", "Setting"))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("    %-18s  %5.2f  %5.2f  %8.2f  %8.2f  %8.2f\n",
                sub$setting[i], sub$a_kappa[i], sub$b_kappa[i],
                sub$score_b1[i], sub$score_b2[i], sub$birth_total[i]))
  }
  cat("\n")
}

cat("\nINTERPRETATION:\n")
cat("The 'total' column is the log-marginal kappa score for a birth.\n")
cat("This enters the birth vs join comparison directly.\n")
cat("Differences of >5 log-units dominate the birth/death balance.\n\n")

# ======================================================================
# PART 3: Cross-dataset comparison of implicit birth penalty
#
# For each dataset and hyper setting, compute:
#   penalty = [birth_score under setting] - [birth_score under 'ideal']
# where 'ideal' sets b = a * T / R (i.e., prior mean = observed rate).
# This isolates the density-driven bias.
# ======================================================================
cat("================================================================\n")
cat("PART 3: DENSITY-DRIVEN BIRTH PENALTY (vs ideal prior)\n")
cat("================================================================\n\n")

penalty_rows <- list()
for (ds in unique(birth_df$dataset)) {
  sub <- birth_df[birth_df$dataset == ds, ]
  R1 <- sub$R_b1[1]; T1 <- sub$T_b1[1]
  R2 <- sub$R_b2[1]; T2 <- sub$T_b2[1]

  # "Ideal" b: set prior mean = observed rates
  b_ideal_1 <- if (R1 > 0) 1 * T1 / R1 else 1  # a=1, so b = T/R
  b_ideal_2 <- if (R2 > 0) 1 * T2 / R2 else 1
  ideal_score <- gp_marg(R1, T1, 1, b_ideal_1) +
                 gp_marg(R2, T2, 1, b_ideal_2)

  for (i in seq_len(nrow(sub))) {
    penalty_rows[[length(penalty_rows) + 1]] <- data.frame(
      dataset = ds,
      setting = sub$setting[i],
      birth_total = sub$birth_total[i],
      ideal_score = round(ideal_score, 2),
      penalty = round(sub$birth_total[i] - ideal_score, 2),
      stringsAsFactors = FALSE
    )
  }
}

pen_df <- do.call(rbind, penalty_rows)
cat(sprintf("%-18s  %-18s  birth    ideal    penalty\n", "Dataset", "Setting"))
for (i in seq_len(nrow(pen_df))) {
  cat(sprintf("%-18s  %-18s  %7.2f  %7.2f  %7.2f\n",
              pen_df$dataset[i], pen_df$setting[i],
              pen_df$birth_total[i], pen_df$ideal_score[i], pen_df$penalty[i]))
}

# ======================================================================
# PART 4: Short MCMC runs comparing strategies
#
# For 2 key datasets (one moderate-density, one high-density), run
# short MCMC (1000 iter) under:
#   A) Current rule: b_kappa = c * mean_deg / K_sp
#   B) sample_b_kappa = TRUE (hierarchical, vague hyperprior)
#   C) Fixed a=b=1
#
# Compare: K trace, final K, b_kappa trace (for B)
# ======================================================================
cat("\n\n================================================================\n")
cat("PART 4: SHORT MCMC COMPARISON\n")
cat("================================================================\n\n")

SHORT_ITER <- 1000
SHORT_BURN <- 200
SHORT_THIN <- 1
SEED <- 42

# Focus on 2 representative datasets
mcmc_datasets <- c("moreno_sheep", "strauss_2019b")

for (ds in mcmc_datasets) {
  A <- tryCatch(choose_dataset(ds), error = function(e) NULL)
  if (is.null(A)) next

  n <- nrow(A)
  mean_deg <- sum(A) / n
  K_sp <- estimate_K_spectral(A)
  density_val <- mean_deg / (n - 1)

  cat(sprintf("\n--- %s (n=%d, density=%.3f, mean_deg=%.1f, K_sp=%d) ---\n\n",
              ds, n, density_val, mean_deg, K_sp))

  deg <- rowSums(A) + colSums(A)
  cv_deg <- sd(deg) / mean(deg)
  if (!is.finite(cv_deg) || cv_deg < 0.1) cv_deg <- 0.1
  a_eta <- 1 / cv_deg^2
  sigma0 <- sqrt(4 / 13)
  tau0 <- qlogis(0.75) / (max(1, K_sp - 1) * sqrt(2 / pi))

  strategies <- list(
    # A: Current rule
    list(
      name = "A: current_rule",
      a_kappa = 1, b_kappa = 1.0 * mean_deg / K_sp,
      sample_b_kappa = FALSE
    ),
    # B: Hierarchical b_kappa
    list(
      name = "B: sample_b_kappa",
      a_kappa = 1, b_kappa = 1.0,  # initial value; will be learned
      sample_b_kappa = TRUE
    ),
    # C: Fixed vague
    list(
      name = "C: fixed a=b=1",
      a_kappa = 1, b_kappa = 1,
      sample_b_kappa = FALSE
    )
  )

  for (strat in strategies) {
    cat(sprintf("  Strategy: %s (a_kappa=%.1f, b_kappa=%.2f, sample_b_kappa=%s)\n",
                strat$name, strat$a_kappa, strat$b_kappa, strat$sample_b_kappa))

    t0 <- proc.time()
    out <- tryCatch(
      modular_osbm_sampler(
        A = A, K = n,
        n_iter = SHORT_ITER, burn = SHORT_BURN, thin = SHORT_THIN,
        verbose = FALSE,
        psi_constraint = "WST" ,
        partition_prior = "OCRP",
        theta_ocrp = 0.5,
        a_kappa = strat$a_kappa, b_kappa = strat$b_kappa,
        a_eta = a_eta, b_eta = a_eta,
        mu0 = 0, sigma0 = sigma0, tau0 = tau0,
        use_mixing_moves = TRUE,
        sample_b_kappa = strat$sample_b_kappa,
        seed = SEED
      ),
      error = function(e) { cat(sprintf("    ERROR: %s\n", e$message)); NULL }
    )
    elapsed <- (proc.time() - t0)[["elapsed"]]

    if (is.null(out)) next

    # Extract K trace
    z_mat <- if (is.matrix(out$z)) out$z else do.call(rbind, out$z)
    K_trace <- apply(z_mat, 1, function(z) length(unique(z[!is.na(z)])))
    S <- length(K_trace)

    cat(sprintf("    Time: %.1f s | Samples: %d\n", elapsed, S))
    cat(sprintf("    K: mean=%.2f, median=%d, range=[%d, %d]\n",
                mean(K_trace), median(K_trace), min(K_trace), max(K_trace)))
    cat(sprintf("    K (last 25%%): mean=%.2f, median=%d\n",
                mean(K_trace[ceiling(S*0.75):S]),
                median(K_trace[ceiling(S*0.75):S])))

    # b_kappa trace (if hierarchical)
    if (strat$sample_b_kappa && !is.null(out$b_kappa_trace)) {
      bk <- out$b_kappa_trace
      bk <- bk[!is.na(bk) & bk > 0]
      if (length(bk) > 0) {
        cat(sprintf("    b_kappa (learned): mean=%.3f, median=%.3f, range=[%.3f, %.3f]\n",
                    mean(bk), median(bk), min(bk), max(bk)))
        cat(sprintf("    b_kappa (last 25%%): mean=%.3f\n",
                    mean(bk[ceiling(length(bk)*0.75):length(bk)])))
      }
    }
    cat("\n")
  }
}

# ======================================================================
# PART 5: Summary and recommendations
# ======================================================================
cat("\n================================================================\n")
cat("SUMMARY AND ANALYSIS\n")
cat("================================================================\n\n")

cat("THE PROBLEM:\n")
cat("  The collapsed Gamma-Poisson marginal in the birth move is:\n")
cat("    log p(R | T, a, b) = lgamma(a+R) - lgamma(a) + a*log(b) - (a+R)*log(b+T)\n\n")
cat("  The term `a*log(b)` creates an asymmetry: it enters birth scores as a\n")
cat("  CONSTANT PENALTY that depends on b_kappa but NOT on the data.\n")
cat("  Meanwhile `(a+R)*log(b+T)` grows with R (data size).\n\n")
cat("  When b_kappa = c * mean_deg / K_sp (current rule):\n")
cat("    - Dense networks: b is large -> small prior mean a/b\n")
cat("    - The marginal penalizes any kappa >> a/b\n")
cat("    - Birth needs kappa values that match the observed rate R/T\n")
cat("    - If R/T >> a/b, birth is penalized -> model over-splits to\n")
cat("      create many small blocks where R is small enough\n\n")
cat("  This is backwards: dense networks get MORE births, not fewer.\n\n")

cat("THE FIX:\n")
cat("  Option A (recommended): sample_b_kappa = TRUE\n")
cat("    - Fully Bayesian hierarchical prior: b_kappa ~ Gamma(1, 0.01)\n")
cat("    - b_kappa learns the right scale from the current kappa values\n")
cat("    - NOT empirical Bayes: uncertainty in b_kappa propagates\n")
cat("    - Already implemented; just enable it\n\n")
cat("  Option B (simpler): Fixed a_kappa = 1, b_kappa = 1\n")
cat("    - Prior mean = 1, prior CV = 100%\n")
cat("    - Works well for moderate densities but can conflict for\n")
cat("      very sparse or very dense networks\n")
cat("    - No adaptation, but removes the density-proportional scaling\n\n")

cat("DONE.\n")
