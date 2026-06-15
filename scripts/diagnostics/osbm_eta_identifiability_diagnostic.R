#!/usr/bin/env Rscript

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

suppressPackageStartupMessages({
  source("./helper_folder/helper.R")
  source("./helper_folder/SST_helpers.R")
  source("./helper_folder/WST_helpers.R")
  source("./helper_folder/Hyper_setup.R")
  source("./core/my_best_try_so_far.R")
})

if (!exists("log1pexp")) log1pexp <- function(x) ifelse(x < 35, log1p(exp(x)), x)

N_ITER <- as.integer(Sys.getenv("N_ITER", unset = "900"))
BURN <- as.integer(Sys.getenv("BURN", unset = "1"))
THIN <- as.integer(Sys.getenv("THIN", unset = "1"))
SEED <- as.integer(Sys.getenv("SEED", unset = "42"))
THETA_OCRP <- as.numeric(Sys.getenv("THETA_OCRP", unset = "0.5"))
C_KAPPA <- as.numeric(Sys.getenv("C_KAPPA", unset = "1"))
A_KAPPA <- as.numeric(Sys.getenv("A_KAPPA", unset = "1"))

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path("output", "diagnostics", paste0("osbm_eta_identifiability_", run_stamp))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

cat("OSBM eta identifiability diagnostic\n")
cat("run_dir:", run_dir, "\n")
cat(sprintf("n_iter=%d burn=%d thin=%d seed=%d theta_ocrp=%.2f c_kappa=%.2f\n\n",
            N_ITER, BURN, THIN, SEED, THETA_OCRP, C_KAPPA))

if (BURN < 1L) {
  stop("BURN must be >= 1 for the current sampler implementation.", call. = FALSE)
}

DATASET_CHOICES <- c("moreno_sheep", "strauss_2019b", "citations_data")

choose_dataset <- function(dataset = DATASET_CHOICES[1]) {
  dataset <- match.arg(dataset, choices = DATASET_CHOICES)
  if (dataset == "citations_data") {
    A <- read.csv("./data/Citations_application/cross-citation-matrix.csv",
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
    A <- matrix(0L, n_nodes, n_nodes, dimnames = list(seq_len(n_nodes), seq_len(n_nodes)))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]
      j <- edges$target[r]
      w <- edges$weight[r]
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
    A <- matrix(0L, n_nodes, n_nodes, dimnames = list(seq_len(n_nodes), seq_len(n_nodes)))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]
      j <- edges$target[r]
      if (edges$weight[r] > 0L) A[i, j] <- A[i, j] + 1L
    }
    diag(A) <- 0L
  } else {
    stop("Unknown dataset.")
  }
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  colnames(A) <- rownames(A)
  A
}

get_ocrp_hypers <- function(A, c_kappa = 1, a_kappa = 1, theta_ocrp = 0.5) {
  n <- nrow(A)
  K_sp <- estimate_K_spectral(A)
  mean_deg <- sum(A) / n
  b_kappa <- c_kappa * mean_deg / K_sp

  deg <- rowSums(A) + colSums(A)
  cv_deg <- sd(deg) / mean(deg)
  if (!is.finite(cv_deg) || cv_deg < 0.1) cv_deg <- 0.1
  a_eta <- 1 / (cv_deg^2)

  sigma0 <- sqrt(4 / 13)
  tau0 <- qlogis(0.75) / (max(1, K_sp - 1) * sqrt(2 / pi))

  list(
    a_kappa = a_kappa,
    b_kappa = b_kappa,
    a_eta = a_eta,
    b_eta = a_eta,
    mu0 = 0,
    sigma0 = sigma0,
    tau0 = tau0,
    theta_ocrp = theta_ocrp,
    partition_prior = "OCRP",
    gamma_gn = NA_real_,
    alpha_gp05 = NA_real_,
    meta = list(K_sp = K_sp, mean_deg = mean_deg, cv_deg = cv_deg)
  )
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok]))
}

build_edge_cache <- function(A) {
  n <- nrow(A)
  N_mat <- A + t(A)
  diag(N_mat) <- 0
  idx <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  i_idx <- idx[, 1]
  j_idx <- idx[, 2]
  N_edge <- as.numeric(N_mat[idx])
  edge_by_node <- replicate(n, integer(0), simplify = FALSE)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], as.integer(e))
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], as.integer(e))
  }
  list(i_idx = i_idx, j_idx = j_idx, N_edge = N_edge, edge_by_node = edge_by_node)
}

partition_log_weights <- function(v_minus, partition_prior, gamma_gn, alpha_gp05, theta_ocrp) {
  if (identical(partition_prior, "OCRP")) {
    ocrp_log_weights_packed(v_minus, theta_ocrp = theta_ocrp)
  } else if (identical(partition_prior, "ROCRP")) {
    rocrp_log_weights_packed(v_minus, theta_ocrp = theta_ocrp)
  } else if (identical(partition_prior, "GP05")) {
    gp05_log_weights_packed(v_minus, alpha_gp05 = alpha_gp05, gamma_gp05 = gamma_gn)
  } else {
    gn_log_weights_packed(v_minus, gamma_gn = gamma_gn)
  }
}

node_birth_decomposition <- function(i, A, z, eta, kappa, psi, regime,
                                     Rkl, Tkl, cache, hyper, node_deg) {
  K_full <- nrow(kappa)
  oldk <- z[i]
  eta_i <- eta[i]

  n_minus_full <- tabulate(z[-i], nbins = K_full)
  keep_full <- sort(which(n_minus_full > 0L))
  K_minus <- length(keep_full)
  if (K_minus < 1L) return(NULL)

  map <- integer(K_full)
  map[keep_full] <- seq_len(K_minus)
  oldk_pos <- match(oldk, keep_full)

  if (regime == "WST") {
    psi_minus <- psi[keep_full, keep_full, drop = FALSE]
    diag(psi_minus) <- 0
    psi_minus[lower.tri(psi_minus)] <- t(psi_minus)[lower.tri(psi_minus)]
  } else {
    psi_minus <- reindex_psi_sst_keep(psi_old = psi, K_old = K_full, keep = keep_full)
  }

  C_i_full <- counts_by_block_exact_cpp(
    i = i,
    A = A,
    z = as.integer(z),
    i_idx = as.integer(cache$i_idx),
    j_idx = as.integer(cache$j_idx),
    N_edge = as.numeric(cache$N_edge),
    edge_by_node = cache$edge_by_node,
    K = K_full
  )
  c_plus <- as.numeric(C_i_full$c_plus)[keep_full]
  N_tot <- as.numeric(C_i_full$N_tot)[keep_full]

  E_excl_full <- E_by_block_excluding_i(i, z, eta, K_full)
  E_excl <- as.numeric(E_excl_full[keep_full])

  r_add <- N_tot
  t_add <- eta_i * E_excl

  R_minus <- Rkl[keep_full, keep_full, drop = FALSE]
  T_minus <- Tkl[keep_full, keep_full, drop = FALSE]

  if (!is.na(oldk_pos)) {
    for (ell_pos in seq_len(K_minus)) {
      subR <- N_tot[ell_pos]
      subT <- eta_i * E_excl[ell_pos]
      p <- min(oldk_pos, ell_pos)
      q <- max(oldk_pos, ell_pos)
      R_minus[p, q] <- max(R_minus[p, q] - subR, 0)
      T_minus[p, q] <- max(T_minus[p, q] - subT, 0)
      R_minus[q, p] <- R_minus[p, q]
      T_minus[q, p] <- T_minus[p, q]
    }
  }

  r_set <- seq_len(K_minus + 1L)

  if (regime == "WST") {
    lp_dir_exist <- vapply(seq_len(K_minus), function(kc) {
      dir_exact_WST_existing_counts(
        k = kc, c_plus = c_plus, N_tot = N_tot,
        psi_mat = psi_minus, include_within = TRUE
      )
    }, numeric(1))
    lp_dir_new <- rep(-Inf, K_minus + 1L)
    lp_dir_new[r_set] <- dir_exact_WST_new_vec_counts(
      c_plus = c_plus, N_tot = N_tot,
      mu0 = hyper$mu0, sigma0 = hyper$sigma0, r_set = r_set
    )
  } else {
    lp_dir_exist <- vapply(seq_len(K_minus), function(kc) {
      dir_exact_SST_existing_counts(k = kc, c_plus = c_plus, N_tot = N_tot, psi_vec = psi_minus)
    }, numeric(1))
    lp_dir_new <- rep(-Inf, K_minus + 1L)
    lp_dir_new[r_set] <- dir_exact_SST_new_vec_counts(
      c_plus = c_plus, N_tot = N_tot, psi_vec = psi_minus,
      r_set = r_set, tau0 = hyper$tau0
    )
  }

  lp_kappa_exist <- numeric(K_minus)
  for (kc in seq_len(K_minus)) {
    acc <- 0
    for (ell in seq_len(K_minus)) {
      if (r_add[ell] == 0 && t_add[ell] == 0) next
      p <- min(kc, ell)
      q <- max(kc, ell)
      acc <- acc + (
        gp_marginal(R_minus[p, q] + r_add[ell], T_minus[p, q] + t_add[ell],
                    hyper$a_kappa, hyper$b_kappa) -
          gp_marginal(R_minus[p, q], T_minus[p, q], hyper$a_kappa, hyper$b_kappa)
      )
    }
    lp_kappa_exist[kc] <- acc
  }
  lp_kappa_new <- sum(gp_marginal(r_add, t_add, hyper$a_kappa, hyper$b_kappa))

  v_minus <- n_minus_full[keep_full]
  pw <- partition_log_weights(
    v_minus = v_minus,
    partition_prior = hyper$partition_prior,
    gamma_gn = hyper$gamma_gn,
    alpha_gp05 = hyper$alpha_gp05,
    theta_ocrp = hyper$theta_ocrp
  )
  lp_prior_exist <- pw$exist
  lp_prior_new <- pw$new
  per_slot <- isTRUE(pw$per_slot)

  total_exist <- lp_dir_exist + lp_kappa_exist + lp_prior_exist
  total_new <- lp_dir_new + lp_kappa_new + lp_prior_new

  logW_join <- lse(total_exist)
  if (per_slot) {
    logW_new <- lse(total_new[r_set])
  } else {
    logW_new <- lse(total_new[r_set] - log(length(r_set)))
  }
  p_birth <- exp(logW_new - lse(c(logW_new, logW_join)))

  best_k <- which.max(total_exist)
  best_r <- r_set[which.max(total_new[r_set])]
  prior_new_best <- if (length(lp_prior_new) > 1L) lp_prior_new[best_r] else lp_prior_new
  new_block_total <- lp_dir_new[best_r] + lp_kappa_new + prior_new_best

  n_k <- as.integer(tabulate(z, nbins = K_full))
  E_k <- tapply(eta, factor(z, levels = seq_len(K_full)), sum)
  E_k[is.na(E_k)] <- 0

  data.frame(
    node = i,
    degree = node_deg[i],
    eta = eta_i,
    block = z[i],
    block_size = n_k[z[i]],
    block_mass = E_k[z[i]],
    p_birth = p_birth,
    dir_delta = lp_dir_new[best_r] - lp_dir_exist[best_k],
    kap_delta = lp_kappa_new - lp_kappa_exist[best_k],
    prior_delta = prior_new_best - lp_prior_exist[best_k],
    total_delta = new_block_total - total_exist[best_k],
    t_add_sum = sum(t_add),
    stringsAsFactors = FALSE
  )
}

analyse_fit <- function(out, A, regime, hyper, dataset, eta_identifiability, fit_label) {
  cache <- build_edge_cache(A)
  n <- nrow(A)
  node_deg <- rowSums(A) + colSums(A)
  n_save <- length(out$z)
  if (n_save < 1L) stop("No saved draws in fit.")

  iter_rows <- vector("list", n_save)
  node_rows <- vector("list", n_save)

  for (s in seq_len(n_save)) {
    z <- as.integer(out$z[[s]])
    eta <- as.numeric(out$eta[[s]])
    kappa <- out$kappa[[s]]
    psi <- out$psi[[s]]
    K <- length(unique(z))

    bt <- block_totals_for_poisson_cpp(z, eta, cache$i_idx, cache$j_idx, cache$N_edge, K)
    node_df <- do.call(
      rbind,
      lapply(seq_len(n), function(i) {
        node_birth_decomposition(
          i = i, A = A, z = z, eta = eta, kappa = kappa, psi = psi,
          regime = regime, Rkl = bt$Rkl, Tkl = bt$Tkl,
          cache = cache, hyper = hyper, node_deg = node_deg
        )
      })
    )
    node_df$iter <- s
    node_df$dataset <- dataset
    node_df$regime <- regime
    node_df$eta_identifiability <- eta_identifiability

    n_k <- as.integer(tabulate(z, nbins = K))
    E_k <- tapply(eta, factor(z, levels = seq_len(K)), sum)
    E_k[is.na(E_k)] <- 0
    mass_ratio <- E_k / pmax(n_k, 1)
    singleton_mask <- n_k[z] == 1L
    ut <- upper.tri(kappa, diag = TRUE)

    iter_rows[[s]] <- data.frame(
      iter = s,
      dataset = dataset,
      regime = regime,
      eta_identifiability = eta_identifiability,
      fit = fit_label,
      K = K,
      singleton_blocks = sum(n_k == 1L),
      mean_p_birth = mean(node_df$p_birth),
      median_p_birth = median(node_df$p_birth),
      max_p_birth = max(node_df$p_birth),
      mean_dir_delta = mean(node_df$dir_delta),
      mean_kap_delta = mean(node_df$kap_delta),
      mean_prior_delta = mean(node_df$prior_delta),
      mean_total_delta = mean(node_df$total_delta),
      eta_cv = stats::sd(eta) / mean(eta),
      eta_min = min(eta),
      eta_q10 = as.numeric(stats::quantile(eta, 0.10)),
      eta_q90 = as.numeric(stats::quantile(eta, 0.90)),
      block_mass_cv = if (length(mass_ratio) > 1L) stats::sd(mass_ratio) / mean(mass_ratio) else 0,
      block_mass_min = min(mass_ratio),
      block_mass_max = max(mass_ratio),
      mean_eta_singletons = if (any(singleton_mask)) mean(eta[singleton_mask]) else NA_real_,
      corr_pbirth_eta = safe_cor(node_df$p_birth, log(pmax(node_df$eta, 1e-12))),
      corr_pbirth_degree = safe_cor(node_df$p_birth, node_df$degree),
      corr_eta_degree = safe_cor(log(pmax(eta, 1e-12)), node_deg),
      mean_kappa = mean(kappa[ut]),
      stringsAsFactors = FALSE
    )
    node_rows[[s]] <- node_df
  }

  iter_df <- do.call(rbind, iter_rows)
  node_df <- do.call(rbind, node_rows)
  tail_start <- max(1L, floor(nrow(iter_df) / 2))
  iter_tail <- iter_df[tail_start:nrow(iter_df), , drop = FALSE]
  node_tail <- node_df[node_df$iter >= tail_start, , drop = FALSE]

  node_summary <- aggregate(
    cbind(p_birth, eta, block_size, block_mass) ~ dataset + regime + eta_identifiability + node + degree,
    data = node_tail,
    FUN = mean
  )
  singleton_prob <- aggregate((block_size == 1) ~ dataset + regime + eta_identifiability + node,
                              data = node_tail, FUN = mean)
  names(singleton_prob)[names(singleton_prob) == "(block_size == 1)"] <- "singleton_prob"
  node_summary <- merge(node_summary, singleton_prob,
                        by = c("dataset", "regime", "eta_identifiability", "node"),
                        all.x = TRUE, sort = FALSE)

  summary_df <- data.frame(
    dataset = dataset,
    regime = regime,
    eta_identifiability = eta_identifiability,
    fit = fit_label,
    n = nrow(A),
    K_mean_all = mean(iter_df$K),
    K_mean_tail = mean(iter_tail$K),
    K_median_tail = stats::median(iter_tail$K),
    K_max = max(iter_df$K),
    singleton_blocks_tail = mean(iter_tail$singleton_blocks),
    mean_p_birth_tail = mean(iter_tail$mean_p_birth),
    max_p_birth_tail = mean(iter_tail$max_p_birth),
    mean_total_delta_tail = mean(iter_tail$mean_total_delta),
    mean_kap_delta_tail = mean(iter_tail$mean_kap_delta),
    mean_dir_delta_tail = mean(iter_tail$mean_dir_delta),
    mean_prior_delta_tail = mean(iter_tail$mean_prior_delta),
    block_mass_cv_tail = mean(iter_tail$block_mass_cv),
    block_mass_min_tail = mean(iter_tail$block_mass_min),
    block_mass_max_tail = mean(iter_tail$block_mass_max),
    mean_eta_singletons_tail = mean(iter_tail$mean_eta_singletons, na.rm = TRUE),
    corr_K_birth = safe_cor(iter_df$K, iter_df$mean_p_birth),
    corr_K_block_mass_cv = safe_cor(iter_df$K, iter_df$block_mass_cv),
    corr_birth_block_mass_cv = safe_cor(iter_df$mean_p_birth, iter_df$block_mass_cv),
    corr_birth_eta_tail = mean(iter_tail$corr_pbirth_eta, na.rm = TRUE),
    corr_birth_degree_tail = mean(iter_tail$corr_pbirth_degree, na.rm = TRUE),
    corr_eta_degree_tail = mean(iter_tail$corr_eta_degree, na.rm = TRUE),
    top_birth_nodes = paste(
      node_summary$node[order(node_summary$p_birth, decreasing = TRUE)][seq_len(min(5L, nrow(node_summary)))],
      collapse = ","
    ),
    stringsAsFactors = FALSE
  )

  list(iter = iter_df, node = node_df, node_summary = node_summary, summary = summary_df)
}

run_fit <- function(A, regime, hyper, eta_identifiability, seed) {
  t0 <- proc.time()[3]
  out <- modular_osbm_sampler(
    A = A,
    K = nrow(A),
    n_iter = N_ITER,
    burn = BURN,
    thin = THIN,
    verbose = FALSE,
    psi_constraint = regime ,
    eta_identifiability = eta_identifiability,
    partition_prior = hyper$partition_prior,
    theta_ocrp = hyper$theta_ocrp,
    a_kappa = hyper$a_kappa,
    b_kappa = hyper$b_kappa,
    a_eta = hyper$a_eta,
    b_eta = hyper$b_eta,
    mu0 = hyper$mu0,
    sigma0 = hyper$sigma0,
    tau0 = hyper$tau0,
    use_mixing_moves = TRUE,
    sample_b_kappa = FALSE,
    seed = seed
  )
  elapsed <- proc.time()[3] - t0
  list(out = out, elapsed = elapsed)
}

screen_rows <- list()
all_iter <- list()
all_node_summary <- list()

datasets_main <- c("moreno_sheep", "strauss_2019b")
eta_modes <- c("none", "block_sum")
regimes <- c("WST", "SST")

for (dataset in datasets_main) {
  A <- choose_dataset(dataset)
  hyper <- get_ocrp_hypers(A, c_kappa = C_KAPPA, a_kappa = A_KAPPA, theta_ocrp = THETA_OCRP)
  cat(sprintf("Dataset=%s n=%d mean_deg=%.2f K_sp=%d b_kappa=%.3f a_eta=%.3f\n",
              dataset, nrow(A), hyper$meta$mean_deg, hyper$meta$K_sp,
              hyper$b_kappa, hyper$a_eta))

  for (regime in regimes) {
    for (eta_mode in eta_modes) {
      fit_label <- paste(dataset, regime, eta_mode, sep = "__")
      cat(sprintf("  fitting %s\n", fit_label))
      fit <- run_fit(A, regime, hyper, eta_mode, SEED)
      saveRDS(fit$out, file.path(run_dir, paste0(fit_label, "_fit.rds")))
      diag_out <- analyse_fit(
        out = fit$out, A = A, regime = regime, hyper = hyper,
        dataset = dataset, eta_identifiability = eta_mode, fit_label = fit_label
      )
      diag_out$summary$elapsed_sec <- fit$elapsed
      screen_rows[[length(screen_rows) + 1L]] <- diag_out$summary
      all_iter[[fit_label]] <- diag_out$iter
      all_node_summary[[fit_label]] <- diag_out$node_summary
      write.csv(diag_out$iter, file.path(run_dir, paste0(fit_label, "_iter.csv")), row.names = FALSE)
      write.csv(diag_out$node_summary, file.path(run_dir, paste0(fit_label, "_node_summary.csv")), row.names = FALSE)
    }
  }
  cat("\n")
}

screen_df <- do.call(rbind, screen_rows)
write.csv(screen_df, file.path(run_dir, "screen_summary.csv"), row.names = FALSE)

strauss_uncon <- subset(screen_df, dataset == "strauss_2019b" & eta_identifiability == "none")
focal_regime <- strauss_uncon$regime[which.max(strauss_uncon$K_mean_tail)]
block_gain <- merge(
  subset(screen_df, dataset == "strauss_2019b" & eta_identifiability == "none", c(regime, K_mean_tail, mean_p_birth_tail, block_mass_cv_tail)),
  subset(screen_df, dataset == "strauss_2019b" & eta_identifiability == "block_sum", c(regime, K_mean_tail, mean_p_birth_tail, block_mass_cv_tail)),
  by = "regime", suffixes = c("_none", "_block")
)
clear_pattern <- any((block_gain$K_mean_tail_none - block_gain$K_mean_tail_block) >= 2)

validation_df <- NULL
if (clear_pattern) {
  dataset <- "citations_data"
  A <- choose_dataset(dataset)
  hyper <- get_ocrp_hypers(A, c_kappa = C_KAPPA, a_kappa = A_KAPPA, theta_ocrp = THETA_OCRP)
  cat(sprintf("Validation on %s with focal regime %s\n", dataset, focal_regime))
  val_rows <- list()
  for (eta_mode in eta_modes) {
    fit_label <- paste(dataset, focal_regime, eta_mode, sep = "__")
    cat(sprintf("  fitting %s\n", fit_label))
    fit <- run_fit(A, focal_regime, hyper, eta_mode, SEED)
    saveRDS(fit$out, file.path(run_dir, paste0(fit_label, "_fit.rds")))
    diag_out <- analyse_fit(
      out = fit$out, A = A, regime = focal_regime, hyper = hyper,
      dataset = dataset, eta_identifiability = eta_mode, fit_label = fit_label
    )
    diag_out$summary$elapsed_sec <- fit$elapsed
    val_rows[[length(val_rows) + 1L]] <- diag_out$summary
    write.csv(diag_out$iter, file.path(run_dir, paste0(fit_label, "_iter.csv")), row.names = FALSE)
    write.csv(diag_out$node_summary, file.path(run_dir, paste0(fit_label, "_node_summary.csv")), row.names = FALSE)
  }
  validation_df <- do.call(rbind, val_rows)
  write.csv(validation_df, file.path(run_dir, "validation_summary.csv"), row.names = FALSE)
}

summary_md <- file.path(run_dir, "summary.md")
con <- file(summary_md, open = "wt")
writeLines("# OSBM Eta Identifiability Diagnostic", con)
writeLines("", con)
writeLines(sprintf("- Run directory: `%s`", run_dir), con)
writeLines(sprintf("- Main settings: `n_iter=%d`, `burn=%d`, `thin=%d`, `seed=%d`", N_ITER, BURN, THIN, SEED), con)
writeLines(sprintf("- OCRP hyper rule: `theta_ocrp=%.2f`, `c_kappa=%.2f`, `a_kappa=%.2f`", THETA_OCRP, C_KAPPA, A_KAPPA), con)
writeLines("", con)
writeLines("## Screen Summary", con)
write.table(screen_df, con, sep = ",", row.names = FALSE, col.names = TRUE)
if (!is.null(validation_df)) {
  writeLines("", con)
  writeLines("## Validation Summary", con)
  write.table(validation_df, con, sep = ",", row.names = FALSE, col.names = TRUE)
}
close(con)

cat("\nSaved screen summary to:", file.path(run_dir, "screen_summary.csv"), "\n")
if (!is.null(validation_df)) {
  cat("Saved validation summary to:", file.path(run_dir, "validation_summary.csv"), "\n")
}
cat("Saved markdown summary to:", summary_md, "\n")
