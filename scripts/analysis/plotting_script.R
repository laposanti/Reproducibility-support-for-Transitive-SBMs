# =========================
# All-fits reader + "best K per (dataset, model)" + plots
# =========================

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggside)
  library(glue)
  library(readr)
  library(mcclust)
  library(reshape2)
})

# --- model sources (paths as in your repo)
source("./core/modular mcmc.R")
source("./helper_folder/helper.R")
source("./helper_folder/sim_study_helper.R")
source("./helper_folder/transitivity_check_helper.R")
source("./core/DDCSBM.R")
source("./scripts/analysis/osbm_visualization.R")

# Override relabel_osbm with variable-K-compatible version
extract_z_matrix <- function(mcmc_out) {
  z_raw <- mcmc_out$z_chain %||% mcmc_out$z
  if (is.null(z_raw)) stop("relabel_osbm: cannot find z or z_chain in mcmc_out.")
  if (is.matrix(z_raw)) return(apply(z_raw, 2, as.integer))
  if (is.data.frame(z_raw)) return(as.matrix(z_raw))
  if (is.list(z_raw)) return(do.call(rbind, lapply(z_raw, as.integer)))
  stop("relabel_osbm: unsupported z format (", class(z_raw), ").")
}

relabel_osbm <- function(mcmc_out, A,
                         ordering = c("WST","SST","NONE"),
                         score = c("success","outdeg","indeg","netdeg")) {
  ordering <- match.arg(ordering)
  score    <- match.arg(score)
  draws_z  <- extract_z_matrix(mcmc_out)
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  S <- nrow(draws_z); n <- ncol(draws_z)
  if (ordering != "WST") return(list(z = draws_z, perm = NULL))
  outdeg  <- rowSums(A); matches <- A + t(A)
  success <- outdeg / pmax(rowSums(matches), 1)
  node_metric <- switch(score, outdeg = outdeg, indeg = colSums(A),
                        netdeg = outdeg - colSums(A), success = success)
  relab_z <- matrix(NA_integer_, S, n)
  for (s in seq_len(S)) {
    z_s <- as.integer(draws_z[s, ]); K_s <- max(z_s, na.rm = TRUE)
    means <- vapply(seq_len(K_s), function(k) {
      idx <- which(z_s == k); if (length(idx)) mean(node_metric[idx]) else -Inf
    }, numeric(1))
    ord_old <- order(means, decreasing = TRUE, na.last = TRUE)
    p <- integer(K_s); p[ord_old] <- seq_len(K_s)
    relab_z[s, ] <- p[z_s]
  }
  list(z = relab_z, perm = NULL)
}

dcsbm_relabel <- function(mcmc_out) {
  draws_z <- extract_z_matrix(mcmc_out)
  list(z = draws_z)
}

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
have_BT2 <- requireNamespace("BradleyTerry2", quietly = TRUE)

# same as in the main driver
choose_dataset <- function(dataset = c("mountain_goats","citations_data","macaques_data",
                                       "high_school","moreno_sheep","strauss_2019b")) {
  dataset <- match.arg(dataset)
  if (dataset == "mountain_goats") {
    matrix_files <- list.files("./data/ShizukaMcDonald_Data",
                               full.names = TRUE, pattern = "\\.csv$")
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
    for (i in seq_len(nrow(edge_list))) {
      A[edge_list[i,1], edge_list[i,2]] <- edge_list[i,"V3"]
    }
  } else if (dataset == "high_school") {
    edge_list <- read.csv("./data/high-school/edges.csv", comment.char = "#", header = FALSE)
    colnames(edge_list) <- c("source", "target", "weight")
    nodes <- sort(unique(c(edge_list$source, edge_list$target)))
    node_names <- as.character(nodes)
    A <- matrix(0, length(nodes), length(nodes), dimnames = list(node_names, node_names))
    for (i in seq_len(nrow(edge_list))) {
      A[as.character(edge_list$source[i]), as.character(edge_list$target[i])] <- edge_list$weight[i]
    }
  } else if (dataset == "moreno_sheep") {
    edge_list <- read.csv("./data/moreno_sheep/edges.csv", comment.char = "#", header = FALSE)
    colnames(edge_list) <- c("source", "target", "weight")
    nodes <- sort(unique(c(edge_list$source, edge_list$target)))
    node_names <- as.character(nodes)
    A <- matrix(0, length(nodes), length(nodes), dimnames = list(node_names, node_names))
    for (i in seq_len(nrow(edge_list))) {
      A[as.character(edge_list$source[i]), as.character(edge_list$target[i])] <- edge_list$weight[i]
    }
  } else if (dataset == "strauss_2019b") {
    edge_list <- read.csv("./data/Strauss_2019b/edges.csv", comment.char = "#", header = FALSE)
    colnames(edge_list) <- c("source", "target", "time", "date")
    edge_agg <- aggregate(time ~ source + target, data = edge_list, FUN = length)
    colnames(edge_agg) <- c("source", "target", "weight")
    nodes <- sort(unique(c(edge_agg$source, edge_agg$target)))
    node_names <- as.character(nodes)
    A <- matrix(0, length(nodes), length(nodes), dimnames = list(node_names, node_names))
    for (i in seq_len(nrow(edge_agg))) {
      A[as.character(edge_agg$source[i]), as.character(edge_agg$target[i])] <- edge_agg$weight[i]
    }
  } else stop("Unknown dataset.")
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  colnames(A) <- rownames(A)
  A
}

# ------------------------------------------------------------
# Bradley–Terry item “strength” (with no hard dependency on BT2)
# ------------------------------------------------------------
bt_item_strength <- function(A) {
  n <- nrow(A); items <- rownames(A)
  pairs <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  d <- tibble::tibble(i = pairs[,1], j = pairs[,2],
                      w_ij = A[cbind(i, j)], w_ji = A[cbind(j, i)]) |>
    dplyr::filter(w_ij + w_ji > 0)
  if (nrow(d) == 0) {
    return(tibble::tibble(item = items, bt = 0, bt_lo = 0, bt_hi = 0))
  }
  if (have_BT2) {
    contests <- dplyr::bind_rows(
      dplyr::transmute(d, winner = i, loser = j, wins = w_ij),
      dplyr::transmute(d, winner = j, loser = i, wins = w_ji)
    ) |>
      dplyr::filter(wins > 0) |>
      dplyr::mutate(
        winner = factor(winner, levels = seq_len(n)),
        loser  = factor(loser,  levels = seq_len(n))
      )
    fit <- BradleyTerry2::BTm(
      outcome = rep(1L, nrow(contests)),
      player1 = contests$winner,
      player2 = contests$loser,
      weights = contests$wins
    )
    ab <- as.data.frame(BradleyTerry2::BTabilities(fit))
    ab$item <- items[as.integer(rownames(ab))]
    ab$ability[!is.finite(ab$ability)] <- 0
    ab$s.e.[!is.finite(ab$s.e.) | ab$s.e. <= 0] <- 1e-6
    tibble::tibble(
      item  = ab$item,
      bt    = ab$ability,
      bt_lo = ab$ability - 1.96*ab$s.e.,
      bt_hi = ab$ability + 1.96*ab$s.e.
    )
  } else {
    # GLM fallback
    contests <- dplyr::bind_rows(
      dplyr::transmute(d, winner = i, loser = j, w = w_ij, l = w_ji),
      dplyr::transmute(d, winner = j, loser = i, w = w_ji, l = w_ij)
    ) |>
      dplyr::filter(w + l > 0)
    Xw <- model.matrix(~ factor(winner) - 1, contests)
    Xl <- model.matrix(~ factor(loser)  - 1, contests)
    M  <- Xw - Xl
    y  <- contests$w / (contests$w + contests$l)
    wt <- contests$w + contests$l
    suppressWarnings({
      fit <- glm(y ~ M - 1, family = binomial(), weights = wt)
    })
    b <- numeric(ncol(M)); b[seq_along(coef(fit))] <- coef(fit)
    b <- b - mean(b, na.rm = TRUE)
    V  <- try(vcov(fit), silent = TRUE)
    se <- rep(sd(b, na.rm = TRUE), length(b))
    if (!inherits(V, "try-error")) se[seq_len(ncol(V))] <- sqrt(pmax(0, diag(V)))
    se[!is.finite(se) | se <= 0] <- 1e-6
    tibble::tibble(
      item  = items[seq_along(b)],
      bt    = as.numeric(b),
      bt_lo = as.numeric(b - 1.96*se),
      bt_hi = as.numeric(b + 1.96*se)
    )
  }
}

# ------------------------------------------------------------
# Rank draws: BT vs “top-block probability”
# ------------------------------------------------------------
rank_draws_bt_top <- function(bt_df, s_top, Sdraw, n_draws = 2000L, seed = 99) {
  stopifnot(nrow(bt_df) == length(s_top))
  set.seed(seed)
  n  <- nrow(bt_df)
  mu <- bt_df$bt; mu[!is.finite(mu)] <- 0
  se <- (bt_df$bt_hi - bt_df$bt_lo) / (2*1.96)
  se[!is.finite(se) | se <= 0] <- 1e-6
  abil_mat <- matrix(
    rnorm(n_draws * n, mean = rep(mu, each = n_draws), sd = rep(se, each = n_draws)),
    nrow = n_draws, ncol = n, byrow = FALSE
  )
  abil_mat <- abil_mat + matrix(rnorm(n_draws * n, 0, 1e-9), nrow = n_draws)  # break ties
  a <- pmax(s_top + 0.5, 1e-6)
  b <- pmax(Sdraw - s_top + 0.5, 1e-6)
  ptop_mat <- matrix(
    rbeta(n_draws * n, rep(a, each = n_draws), rep(b, each = n_draws)),
    nrow = n_draws, ncol = n, byrow = FALSE
  )
  rank_bt  <- t(apply(abil_mat, 1, function(v) rank(-v, ties.method = "average")))
  rank_top <- t(apply(ptop_mat, 1, function(v) rank(-v, ties.method = "average")))
  sp <- kd <- numeric(n_draws)
  for (r in seq_len(n_draws)) {
    sp[r] <- suppressWarnings(cor(rank_bt[r,],  rank_top[r,], method = "spearman"))
    kd[r] <- suppressWarnings(cor(rank_bt[r,],  rank_top[r,], method = "kendall"))
  }
  list(rank_bt = rank_bt, rank_top = rank_top, sp = sp, kd = kd)
}

# ---------- Summarise per-item rank distributions ----------
summarise_rank_draws <- function(rank_mat) {
  tibble::tibble(
    med = apply(rank_mat, 2, median, na.rm = TRUE),
    lo  = apply(rank_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    hi  = apply(rank_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

# ---------- Rho[d] by distance (SST) ----------
posterior_rho_d <- function(mcmc_out, regime = c("SST","WST"), block_maps = NULL, K = NULL) {
  regime <- match.arg(regime)
  invlogit <- function(x) 1/(1+exp(-x))
  if (regime == "SST") {
    psi <- mcmc_out$psi %||% mcmc_out$Psi
    if (is.null(psi) && !is.null(mcmc_out$delta)) {
      del <- as.matrix(mcmc_out$delta)
      psi <- t(apply(del, 1, cumsum))
    }
    # Handle variable-K list-format psi: filter to draws matching target K
    if (is.list(psi) && !is.matrix(psi)) {
      target_len <- if (!is.null(K)) K - 1L else {
        lens <- vapply(psi, length, integer(1))
        as.integer(names(sort(table(lens), decreasing = TRUE))[1])
      }
      keep <- which(vapply(psi, length, integer(1)) == target_len)
      if (length(keep) == 0) {
        warning("No SST psi draws have length K-1 = ", target_len,
                "; using closest available length.")
        lens <- vapply(psi, length, integer(1))
        target_len <- as.integer(names(sort(table(lens), decreasing = TRUE))[1])
        keep <- which(lens == target_len)
      }
      psi <- do.call(rbind, psi[keep])
      if (is.null(K)) K <- target_len + 1L
    } else {
      psi <- as.matrix(psi)
      if (nrow(psi) < ncol(psi)) psi <- t(psi)
      if (is.null(K)) K <- ncol(psi) + 1L
    }
    if (ncol(psi) != (K - 1L)) {
      warning("psi ncol (", ncol(psi), ") != K-1 (", K-1L, "); truncating/padding.")
      d <- K - 1L
      if (ncol(psi) > d) psi <- psi[, 1:d, drop = FALSE]
      else {
        pad <- matrix(NA_real_, nrow(psi), d - ncol(psi))
        psi <- cbind(psi, pad)
      }
    }
    rho_draws <- invlogit(psi)  # S x (K-1)
    summ <- tibble::tibble(
      d  = 1:(K-1),
      mean = colMeans(rho_draws, na.rm = TRUE),
      lo   = apply(rho_draws, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
      hi   = apply(rho_draws, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
    )
    return(list(draws = rho_draws, summ = summ))
  }
  stop("WST branch of posterior_rho_d not needed here (we only call for SST).")
}

# shape summaries for rho[d] (posterior over shape)
summarize_rho_shape <- function(rho_draws) {
  # rho_draws: S x D (D = K-1)
  S <- nrow(rho_draws); D <- ncol(rho_draws)
  if (D < 2L) {
    return(tibble::tibble(
      slope_logit_mean = NA_real_, slope_logit_lo = NA_real_, slope_logit_hi = NA_real_,
      concavity_logit_mean = NA_real_, concavity_logit_lo = NA_real_, concavity_logit_hi = NA_real_,
      early_late_ratio_mean = NA_real_, early_late_ratio_lo = NA_real_, early_late_ratio_hi = NA_real_,
      auc_rho_mean = NA_real_, auc_rho_lo = NA_real_, auc_rho_hi = NA_real_
    ))
  }
  qlogit <- function(p) log(p) - log1p(-p)
  d <- 1:D
  slope   <- conc <- elr <- auc <- numeric(S)
  for (s in 1:S) {
    r  <- pmin(pmax(rho_draws[s,], 1e-6), 1-1e-6)
    lg <- qlogit(r)
    slope[s] <- coef(lm(lg ~ d))[2] %||% NA_real_
    if (D >= 3L) conc[s]  <- mean(diff(lg, differences = 2), na.rm = TRUE) else conc[s] <- NA_real_
    elr[s]   <- r[1] / r[D]
    auc[s]   <- mean(r)  # simple trapezoid-free average over d
  }
  qs <- function(x) c(mean = mean(x, na.rm = TRUE),
                      lo   = quantile(x, 0.025, na.rm = TRUE),
                      hi   = quantile(x, 0.975, na.rm = TRUE))
  s1 <- qs(slope); s2 <- qs(conc); s3 <- qs(elr); s4 <- qs(auc)
  tibble::tibble(
    slope_logit_mean = s1["mean"], slope_logit_lo = s1["lo"], slope_logit_hi = s1["hi"],
    concavity_logit_mean = s2["mean"], concavity_logit_lo = s2["lo"], concavity_logit_hi = s2["hi"],
    early_late_ratio_mean = s3["mean"], early_late_ratio_lo = s3["lo"], early_late_ratio_hi = s3["hi"],
    auc_rho_mean = s4["mean"], auc_rho_lo = s4["lo"], auc_rho_hi = s4["hi"]
  )
}

# ------------------------------------------------------------
# “Implicit” ranks from Z draws (block-first, strength tie-break)
# ------------------------------------------------------------
rank_draws_from_Z <- function(Z_relab, strength) {
  # Z_relab: S x n; strength: length n (fixed)
  S <- nrow(Z_relab); n <- ncol(Z_relab)
  eps <- (rank(strength, ties.method = "average") / (n + 1)) # stable tie-break
  R <- matrix(NA_real_, S, n)
  for (s in 1:S) {
    score <- -as.numeric(Z_relab[s, ]) + 1e-3 * eps
    R[s, ] <- rank(-score, ties.method = "average")
  }
  R
}

# posterior “expected top/second/third” rank for a deterministic figure
implicit_rank_lexi <- function(Z_relab, block_order = NULL, tie_digits = 10) {
  if (!is.matrix(Z_relab)) Z_relab <- as.matrix(Z_relab)
  S <- nrow(Z_relab); n <- ncol(Z_relab)
  K <- max(Z_relab, na.rm = TRUE)
  if (is.null(block_order)) block_order <- seq_len(K)  # 1..K, top-first
  
  # Pk: n x K with Pk[i, k] = Pr(z_i = k)
  Pk <- matrix(NA_real_, n, K)
  for (k in seq_len(K)) Pk[, k] <- colMeans(Z_relab == k, na.rm = TRUE)
  # reorder columns so block 1 priority first, then 2, ...
  Pk <- Pk[, block_order, drop = FALSE]
  
  # Lexicographic order: descending by col1, then col2, ...
  ord_args <- as.list(as.data.frame(-Pk))  # minus for descending
  idx <- do.call(order, c(ord_args, list(method = "radix")))
  
  # Position -> rank; average ranks for exact ties across ALL K probs
  pos <- integer(n); pos[idx] <- seq_len(n)
  keys <- apply(signif(Pk, tie_digits), 1L, function(r) paste(r, collapse = "|"))
  rk <- pos
  for (u in unique(keys)) {
    ix <- which(keys == u)
    rk[ix] <- mean(pos[ix])
  }
  list(rank = rk, Pk = Pk)
}


# ------------------------------------------------------------
# All-fits CSV reader: read from run dir all_results.csv (variable-K runs),
# falling back to old applications_results_summary.csv
# ------------------------------------------------------------
read_all_fits <- function() {
  run_dir <- Sys.getenv("APP_RUN_DIR",
    unset = "./output/application/raw/application_run_20260411_163055")
  # Prefer all_results.csv from the current run dir
  all_res <- file.path(run_dir, "all_results.csv")
  if (file.exists(all_res)) {
    csv_read <- readr::read_csv(all_res, show_col_types = FALSE)
    # Harmonise column names from variable-K format
    if ("K_hat" %in% names(csv_read) && !"K" %in% names(csv_read))
      csv_read <- dplyr::rename(csv_read, K = K_hat)
    if ("elpd_loo" %in% names(csv_read) && !"looic" %in% names(csv_read))
      csv_read <- dplyr::mutate(csv_read, looic = -2 * elpd_loo)
    if (!"K_hat_VI_a05" %in% names(csv_read))
      csv_read$K_hat_VI_a05 <- csv_read$K
    if ("pk_max" %in% names(csv_read) && !"pareto_k_max" %in% names(csv_read))
      csv_read <- dplyr::rename(csv_read, pareto_k_max = pk_max)
    if (!"elpd_loo" %in% names(csv_read) && "looic" %in% names(csv_read))
      csv_read$elpd_loo <- -csv_read$looic / 2
    csv_read <- csv_read |>
      dplyr::mutate(K = as.integer(K))
    message("read_all_fits: loaded ", nrow(csv_read), " rows from ", all_res)
    return(csv_read)
  }
  # Fallback: old summary CSV
  cand <- c("./output/application/raw/applications_results_summary.csv")
  path <- cand[file.exists(cand)][1]
  if (length(path) == 0) stop("Missing all-fits CSV: look for ", paste(cand, collapse = " or "))
  csv_read <- readr::read_csv(path, show_col_types = FALSE)
  csv_read |>
    dplyr::select(dataset, fit_model, K, K_hat_VI_a05, looic, pareto_k_max, dplyr::everything()) |>
    dplyr::mutate(K = as.integer(K))
}

derive_best_tables <- function(all_fits) {
  best_by_model <- all_fits |>
    dplyr::filter(is.finite(looic)) |>
    dplyr::group_by(dataset, fit_model) |>
    dplyr::slice_min(order_by = looic, with_ties = FALSE) |>
    dplyr::ungroup()
  best_overall_per_dataset <- all_fits |>
    dplyr::filter(is.finite(looic)) |>
    dplyr::group_by(dataset) |>
    dplyr::slice_min(order_by = looic, with_ties = FALSE) |>
    dplyr::ungroup()
  list(best_by_model = best_by_model, best_overall = best_overall_per_dataset)
}

save_best_k_csv <- function(best_by_model,
                            out_path = "./output/application/plots/best_k_by_model.csv",
                            mirror_path = "./output/application/plots/best_k_by_model.csv") {
  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(mirror_path), showWarnings = FALSE, recursive = TRUE)
  
  # Use elpd_loo if available, fall back to computing from looic
  bm <- best_by_model
  if (!"elpd_loo" %in% names(bm) && "looic" %in% names(bm))
    bm$elpd_loo <- -bm$looic / 2
  if (!"K_hat_VI_a05" %in% names(bm))
    bm$K_hat_VI_a05 <- bm$K
  if (!"pareto_k_max" %in% names(bm) && "pk_max" %in% names(bm))
    bm <- dplyr::rename(bm, pareto_k_max = pk_max)
  if (!"pareto_k_max" %in% names(bm))
    bm$pareto_k_max <- NA_real_
  
  out_tbl <- bm |>
    dplyr::select(dataset, fit_model, K,
                  dplyr::any_of(c("K_hat_VI_a05", "elpd_loo", "pareto_k_max"))) |>
    dplyr::group_by(dataset) |>
    dplyr::mutate(
      elpd_loo_max = max(elpd_loo, na.rm = TRUE),
      Delta_ELPD = elpd_loo - elpd_loo_max
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::any_of(c("K_hat_VI_a05", "elpd_loo", "pareto_k_max", "Delta_ELPD")),
        .fns  = ~ round(.x, 2)
      )
    ) |>
    dplyr::select(-elpd_loo_max)
  
  readr::write_csv(out_tbl, out_path)
  if (out_path != mirror_path) readr::write_csv(out_tbl, mirror_path)
  message("Saved best-K-by-(dataset,model) to: ", out_path)
}

# ------------------------------------------------------------
# Global “order violation” count (given final z_hat)
#
# Canonical definition (matches helper_folder/transitivity_check_helper.R
# `violation_rate_zhat` and scripts/analysis/audit_hierarchy_synopsis.R):
#   viol_rate = sum_{(i,j) : z_i > z_j} A_{ij}
#               / sum_{(i,j) : z_i != z_j} A_{ij}
# Both numerator and denominator are over *directed* edges.
#
# We do NOT relabel z_hat ex-post. For WST/SST the partition comes from
# `get_z_hat_from_draws(..., method = "draws")` (a single posterior draw
# with the prior's intrinsic 1 = strongest ordering); for DC-SBM every
# draw is already strength-relabelled inside `get_z_hat_from_draws`.
# Applying a second empirical reorder here would re-shuffle WST/SST labels
# against the model's own ordering.
# ------------------------------------------------------------
count_order_violations_global <- function(A, z_hat) {
  stopifnot(length(z_hat) == nrow(A), nrow(A) == ncol(A))
  Z <- outer(z_hat, z_hat, `-`)
  is_cross <- (Z != 0)
  is_viola <- (Z > 0)
  viol_count  <- sum(A[is_viola],  na.rm = TRUE)
  total_cross <- sum(A[is_cross],  na.rm = TRUE)
  viol_rate   <- if (total_cross > 0) viol_count / total_cross else NA_real_
  tibble::tibble(viol_count  = as.double(viol_count),
                 total_cross = as.double(total_cross),
                 viol_rate   = viol_rate)
}

write_violations_csv <- function(base, ds, fm, K, viol_row) {
  out <- viol_row |>
    dplyr::mutate(dataset = ds, fit_model = fm, K = K, .before = 1L)
  path <- paste0(base, "_violations.csv")
  readr::write_csv(out, path)
  path
}

# ------------------------------------------------------------
# Driver
# ------------------------------------------------------------
run_all <- function(mode = c("all", "best_by_model", "best_overall"),
                    OUT_DIR = "./output/application/plots",
                    verbose = FALSE) {
  mode <- match.arg(mode)
  OUT_DIR <- sub("/+$","", OUT_DIR)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Cached fits go to raw/, not plots/
  FIT_DIR <- "./output/application/raw/cached_fits"
  dir.create(FIT_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Run directory for loading existing fits
  RUN_DIR <- Sys.getenv("APP_RUN_DIR",
    unset = "./output/application/raw/application_run_20260411_163055")
  
  all_fits <- read_all_fits()
  bests    <- derive_best_tables(all_fits)
  save_best_k_csv(bests$best_by_model,
                  out_path   = file.path(OUT_DIR, "best_k_by_model.csv"),
                  mirror_path = file.path(OUT_DIR, "best_k_by_model.csv"))
  
  run_tbl <-
    if (mode == "all") {
      all_fits |> distinct(dataset, fit_model, K)
    } else if (mode == "best_by_model") {
      bests$best_by_model |> select(dataset, fit_model, K)
    } else {
      bests$best_overall |> select(dataset, fit_model, K)
    }
  
  rows <- lapply(seq_len(nrow(run_tbl)), function(i) {
    ds <- run_tbl$dataset[i]
    fm <- run_tbl$fit_model[i]
    K  <- as.integer(run_tbl$K[i])
    
    # Per-dataset output directory
    DS_OUT <- file.path(OUT_DIR, ds)
    dir.create(DS_OUT, showWarnings = FALSE, recursive = TRUE)
    
    tryCatch({
      message(glue("Processing {ds} | {fm} | K={K} \n"))
      
      A <- choose_dataset(ds)
      items <- rownames(A)
      # ----- FIT OR LOAD -----
      fit_path <- file.path(FIT_DIR, sprintf("%s__%s__K%d.rds", ds, fm, K))
      
      # Also check run dir for fits (preferred)
      fm_file <- if (fm %in% c("DC-SBM", "DCSBM")) "DCSBM" else fm
      run_fit_path <- file.path(RUN_DIR, paste0(ds, "_", fm_file, "_fit.rds"))
      
      if (file.exists(run_fit_path)) {
        message(sprintf("Loading fit from run dir: %s", basename(run_fit_path)))
        out <- readRDS(run_fit_path)
      } else if (file.exists(fit_path)) {
        message(sprintf("Loading cached fit: %s", basename(fit_path)))
        out <- readRDS(fit_path)
      } else {
        message(sprintf("No existing fit found at %s. Fitting model...", fit_path))
        
        if (fm %in% c("WST", "SST")) {
          out <- fit_osbm(
            A, K, truth = NULL, regime = fm,
            n_iter = 30000L, burn = 5000L, thin = 3L, verbose = verbose
          )
        } else if (fm %in% c("DC-SBM", "DCSBM")) {
          fm <- "DCSBM"
          out <- fit_dcsbm(
            as.matrix(A), K,
            iters = 30000L, burn = 5000L, thin = 3L,
            verbose = if (isTRUE(verbose)) 500 else FALSE,
            seed = 123L
          )
        } else {
          stop("Unknown fit_model: ", fm)
        }
        
        saveRDS(out, fit_path)
        message(sprintf("Model fit saved to %s.", fit_path))
      }
      # ----- RELABEL -----
      if (fm %in% c("WST","SST")) {
        out_relab <- relabel_osbm(A = A, mcmc_out = out, ordering = fm)
        Z_relab   <- out_relab$z
      } else {
        out_relab <- dcsbm_relabel(mcmc_out = out)
        Z_relab   <- out_relab$z
      }


      # Canonical point estimate: same `get_z_hat_from_draws` used by the
      # audit/synopsis pipeline. For WST/SST this returns one posterior draw
      # in the prior's intrinsic 1=strongest order (NO ex-post strength
      # reordering); for DC-SBM it relabels every draw by strength internally
      # before `minVI(method="draws")`. This guarantees plot violation rates
      # match `violation_rate_zhat` in `hierarchy_diagnostics_overview.csv`.
      model_tag <- if (fm %in% c("DC-SBM", "DCSBM")) "DCSBM" else fm
      z_hat <- get_z_hat_from_draws(out, A, model = model_tag)
      z_hat <- as.integer(z_hat)
      K_hat <- max(z_hat, na.rm = TRUE)
      K <- K_hat  # actual number of occupied blocks
      
      # ----- PLOTTING PREP -----
      den <- A + t(A)
      A_prob <- ifelse(den > 0, A/den, NA_real_)
      mwr <- rowSums(A) / pmax(1, rowSums(A) + rowSums(t(A)))
      meta <- tibble(item = items, cl = as.integer(z_hat), strength = mwr)
      
      Y_long <- reshape2::melt(A); colnames(Y_long) <- c("Winner","Loser","Win_Count")
      Y_long <- Y_long %>% mutate(Winner = as.character(Winner), Loser = as.character(Loser))
      Y_long$Matches_Count <- reshape2::melt(den)$value
      Y_long$Y_prob        <- reshape2::melt(A_prob)$value
      
      Y_long_plot <- Y_long %>%
        mutate(perc_success = Win_Count / pmax(1, Matches_Count)) %>%
        left_join(meta, by = c("Loser"  = "item")) %>% rename(row_cl = cl, marginal_row = strength) %>%
        left_join(meta, by = c("Winner" = "item")) %>% rename(col_cl = cl, marginal_col = strength) %>%
        mutate(
          Winner = factor(Winner, levels = unique(Winner[order(col_cl, -marginal_col, decreasing = TRUE)])),
          Loser  = factor(Loser,  levels = unique(Loser [order(row_cl,  -marginal_row )])),
          col_cl = factor(col_cl, ordered = TRUE),
          row_cl = factor(row_cl, ordered = TRUE)
        )
      
      v_lines_list <- Y_long_plot %>% group_by(row_cl) %>%
        summarise(x_break = max(as.numeric(Loser)), .groups = "drop") %>% pull(x_break)
      if (length(v_lines_list) > 0) v_lines_list <- v_lines_list[-length(v_lines_list)]
      h_lines_list <- Y_long_plot %>% group_by(col_cl) %>%
        summarise(y_break = min(as.numeric(Winner)), .groups = "drop") %>% pull(y_break)
      if (length(h_lines_list) > 0) h_lines_list <- h_lines_list[-length(h_lines_list)]
      
      geom_adjacency_row <- ggplot(Y_long_plot, aes(x = Loser, y = Winner)) +
        geom_tile(aes(fill = log1p(Win_Count)), color = 'grey40', show.legend = FALSE) +
        geom_ysidecol(aes(color = factor(col_cl), x = 1), show.legend = FALSE) +
        scale_color_viridis_d(option = "turbo") +
        scale_fill_gradient(low = "white", high = "black", na.value = "grey70") +
        geom_vline(xintercept = v_lines_list + 0.5, color = 'red', linewidth = 0.3) +
        geom_hline(yintercept = h_lines_list - 0.5, color = 'red', linewidth = 0.3) +
        labs(x = "Items (ordered by block)", y = "Items (ordered by block)",
             fill = "log1p(count)", color = "Block") +
        theme_minimal(base_size = 13) +
        theme(axis.text.y = element_text(size = 7), axis.text.x = element_blank(),
              panel.grid = element_blank(), legend.position = "left") +
        theme_ggside_void() + scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
        coord_fixed(ratio = 1)
      
      geom_adjacency_success <- ggplot(Y_long_plot, aes(x = Loser, y = Winner)) +
        geom_tile(aes(fill = Y_prob), color = 'grey40', show.legend = FALSE) +
        geom_ysidecol(aes(color = factor(col_cl), x = 1), show.legend = FALSE) +
        scale_color_viridis_d(option = "turbo") +
        scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue",
                             midpoint = 0.5, na.value = "grey40") +
        geom_vline(xintercept = v_lines_list + 0.5, color = 'red', linewidth = 0.3) +
        geom_hline(yintercept = h_lines_list - 0.5, color = 'red', linewidth = 0.3) +
        labs(x = "Items (ordered by block)", y = "Items (ordered by block)",
             fill = "π(i→j)", color = "Block") +
        theme_minimal(base_size = 13) +
        theme(axis.text.y = element_text(size = 7), axis.text.x = element_blank(),
              panel.grid = element_blank(), legend.position = "left") +
        theme_ggside_void() + scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
        coord_fixed(ratio = 1)
      
      psm_mat <- mcclust::comp.psm(as.matrix(Z_relab))
      dimnames(psm_mat) <- list(items, items)
      psm_long <- reshape2::melt(psm_mat, varnames = c("Winner","Loser"), value.name = "psm") %>%
        mutate(Winner = as.character(Winner), Loser = as.character(Loser)) %>%
        left_join(meta, by = c("Loser"  = "item")) %>% rename(row_cl = cl, marginal_row = strength) %>%
        left_join(meta, by = c("Winner" = "item")) %>% rename(col_cl = cl, marginal_col = strength) %>%
        mutate(
          Winner = factor(Winner, levels = unique(Winner[order(col_cl, -marginal_col, decreasing = TRUE)])),
          Loser  = factor(Loser,  levels = unique(Loser [order(row_cl,  -marginal_row )])),
          col_cl = factor(col_cl, ordered = TRUE)
        )
      geom_sim <- ggplot(psm_long, aes(x = Loser, y = Winner)) +
        geom_tile(aes(fill = psm), color = 'grey40', show.legend = FALSE) +
        geom_ysidecol(aes(color = factor(col_cl), x = 1), show.legend = FALSE) +
        scale_color_viridis_d(option = "turbo") +
        scale_fill_gradient2(low = "white", mid = "grey60", high = "black",
                             midpoint = 0.5, na.value = "grey90") +
        geom_vline(xintercept = v_lines_list + 0.5, color = 'red', linewidth = 0.3) +
        geom_hline(yintercept = h_lines_list - 0.5, color = 'red', linewidth = 0.3) +
        labs(x = "Items (ordered by block)", y = "Items (ordered by block)",
             fill = "Posterior similarity", color = "Block") +
        theme_minimal(base_size = 13) +
        theme(axis.text.y = element_text(size = 7), axis.text.x = element_blank(),
              panel.grid = element_blank(), legend.position = "left") +
        theme_ggside_void() + scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
        coord_fixed(ratio = 1)
      
      base <- file.path(DS_OUT, sprintf("%s__%s__K%d", ds, fm, K))
      ggsave(paste0(base, "_adjacency.png"),  geom_adjacency_row,    width = 7, height = 6, dpi = 300)
      ggsave(paste0(base, "_success.png"),    geom_adjacency_success, width = 7, height = 6, dpi = 300)
      ggsave(paste0(base, "_similarity.png"), geom_sim,               width = 7, height = 6, dpi = 300)
      
      # ----- SST-only analytics -----
      path_rho_d <- NA_character_
      path_rank_vs_rank <- NA_character_
      path_rank_vs_rank_simple <- NA_character_
      path_points <- NA_character_
      rho_shape_tbl <- tibble(
        slope_logit_mean = NA_real_, slope_logit_lo = NA_real_, slope_logit_hi = NA_real_,
        concavity_logit_mean = NA_real_, concavity_logit_lo = NA_real_, concavity_logit_hi = NA_real_,
        early_late_ratio_mean = NA_real_, early_late_ratio_lo = NA_real_, early_late_ratio_hi = NA_real_,
        auc_rho_mean = NA_real_, auc_rho_lo = NA_real_, auc_rho_hi = NA_real_
      )
      
      if (identical(fm, "SST")) {
        rho <- posterior_rho_d(out, regime = "SST", K = K)
        rho_shape_tbl <- summarize_rho_shape(rho$draws)
        rho_plot <- ggplot(rho$summ, aes(x = d, y = mean)) +
          geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
          geom_line() +
          scale_y_continuous(limits = c(0,1)) +
          scale_x_continuous(breaks = 1:(K-1)) +
          labs(x = "Block distance d", y = expression(rho[d]),
               title = glue("{ds} | {fm} | K={K}: distance curve ρ[d]"),
               subtitle = glue("logit-slope: {round(rho_shape_tbl$slope_logit_mean,3)}  ",
                               "[{round(rho_shape_tbl$slope_logit_lo,3)}, {round(rho_shape_tbl$slope_logit_hi,3)}]")) +
          theme_minimal(base_size = 12)
        path_rho_d <- file.path(DS_OUT, sprintf("%s__%s__K%d_rho_d.png", ds, fm, K))
        ggsave(path_rho_d, rho_plot, width = 6, height = 4.5, dpi = 300)
      }
      
      if (fm %in% c('SST','WST')){
        
        # 0) Sanity & convenience
        stopifnot(is.matrix(A) || is.data.frame(A))
        A <- as.matrix(A)
        items <- rownames(A)
        if (is.null(items)) items <- seq_len(nrow(A)) |> as.character()
        
        # 1) BT side
        bt <- bt_item_strength(A)
        s_top <- colSums(Z_relab == 1L)  # assumes block 1 = "top-block"
        # Ensure rows = draws for Z
        Z_use <- Z_relab
        if (ncol(Z_relab) == length(items) && nrow(Z_relab) < ncol(Z_relab)) {
          # looks like items in columns, draws in rows -> OK
        } else if (nrow(Z_relab) == length(items) && ncol(Z_relab) > nrow(Z_relab)) {
          # items in rows, draws in cols -> transpose
          Z_use <- t(Z_relab)
        } else if (! (ncol(Z_relab) == length(items) || nrow(Z_relab) == length(items)) ) {
          stop("Z_relab has incompatible shape: expected one dimension to match #items = ", length(items))
        }
        Sdraw <- nrow(Z_use)
        
        
        # 2) Rank draws from BT (top prob & bt rank)
        rr <- rank_draws_bt_top(bt_df = bt, s_top = s_top, Sdraw = Sdraw, n_draws = 2000L)
        
        # 3) OSBM ranks from Z draws (block-first, strength tie-break)
        #    If mwr is not defined, use BT strength as tie-breaker (higher -> better)
        if (!exists("mwr")) {
          # If bt$bt is on "bigger=stronger" scale, we negate when we later call rank()
          mwr <- as.numeric(bt$bt)
          names(mwr) <- items
        }
        # Ensure strength vector aligned to items
        if (is.null(names(mwr))) names(mwr) <- items
        mwr <- mwr[items]
        
        osbm_rank_draws <- rank_draws_from_Z(Z_use, strength = mwr)
        
        # 4) Summaries
        osbm_rank_summ <- summarise_rank_draws(osbm_rank_draws)
        bt_rank_summ   <- summarise_rank_draws(rr$rank_bt)
        
        # 5) Data for plot (robust join)
        meta <- tibble(item = items) |> dplyr::left_join(meta, by = "item")
        if (!"cl" %in% names(meta)) meta$cl <- 1L
        
        rank_ci_df <- tibble::tibble(
          item   = items,
          bt_med = bt_rank_summ$med, bt_lo = bt_rank_summ$lo, bt_hi = bt_rank_summ$hi,
          osbm_med = osbm_rank_summ$med, osbm_lo = osbm_rank_summ$lo, osbm_hi = osbm_rank_summ$hi
        ) |> dplyr::left_join(meta, by = "item")
        
        # 6) Plot: rank vs rank with uncertainty
        rank_ci_plot <- ggplot(rank_ci_df, aes(x = bt_med, y = osbm_med, color = factor(cl))) +
          geom_segment(aes(x = bt_lo, xend = bt_hi, y = osbm_med, yend = osbm_med), linewidth = 0.35, alpha = 0.9) +
          geom_segment(aes(x = bt_med, xend = bt_med, y = osbm_lo, yend = osbm_hi), linewidth = 0.35, alpha = 0.9) +
          geom_point(size = 1.8) +
          scale_color_viridis_d(option = "turbo", name = "Block") +
          geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, alpha = 0.7) +
          scale_x_reverse() + scale_y_reverse() + coord_equal() +
          labs(x = "BT rank (median, 95% CI)", y = "OSBM rank (median, 95% CI)",
               title = glue("{ds} | {fm} | K={K}: rank vs rank with uncertainty")) +
          theme_minimal(base_size = 12)
        
        path_rank_vs_rank <- paste0(base, "_rank_vs_rank.png")
        ggsave(path_rank_vs_rank, rank_ci_plot, width = 7.2, height = 5.6, dpi = 300)
        
        
        # 7) Deterministic rank vs rank (lexicographic by p(block1), p(block2), ...)
        rk_bt_det <- rank(-bt$bt, ties.method = "average")
        
        # Ensure Z_use has rows=draws, cols=items (you already set Z_use earlier)
        imp_lexi <- implicit_rank_lexi(Z_use, block_order = seq_len(K), tie_digits = 12)
        rk_top_det <- imp_lexi$rank
        
        n_items <- nrow(A)
        rank_simple <- tibble::tibble(
          item = items, bt_rank = rk_bt_det, top_rank = rk_top_det
        ) |> dplyr::left_join(meta, by = "item")
        
        sp_det <- suppressWarnings(cor(rank_simple$bt_rank, rank_simple$top_rank, method = "spearman"))
        
        kd_det <- suppressWarnings(cor(rank_simple$bt_rank, rank_simple$top_rank, method = "kendall"))
        
        rank_simple_plot <- ggplot(rank_simple, aes(x = bt_rank, y = top_rank, color = factor(cl))) +
          geom_point(size = 1.8, alpha = 0.95) +
          geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, alpha = 0.7) +
          scale_color_viridis_d(option = "turbo", name = "Block") +
          scale_x_reverse(limits = c(n_items, 1)) +
          scale_y_reverse(limits = c(n_items, 1)) +
          coord_equal() +
          labs(x = "BT rank",
               y = "Implicit rank (lexi: p(b1), p(b2), …)",
               title = glue("{ds} | {fm} | K={K}: Rank vs Rank (deterministic, lexicographic)"),
               subtitle = glue("Spearman={round(sp_det, 3)}; Kendall={round(kd_det, 3)}")) +
          theme_minimal(base_size = 12)
        
        
        
        
        path_rank_vs_rank_simple <- paste0(base, "_rank_vs_rank_simple.png")
        ggsave(path_rank_vs_rank_simple, rank_simple_plot, width = 7.2, height = 5.6, dpi = 300)
        
        # 8) Save points for post-processing
        top_rank_summ <- summarise_rank_draws(rr$rank_top)
        rank_df <- tibble::tibble(
          item   = items,
          bt_med = bt_rank_summ$med,  bt_lo = bt_rank_summ$lo,  bt_hi = bt_rank_summ$hi,
          top_med = top_rank_summ$med, top_lo = top_rank_summ$lo, top_hi = top_rank_summ$hi
        ) |> dplyr::left_join(meta, by = "item")
        
        path_points <- paste0(base, "_rank_vs_rank_points.csv")
        readr::write_csv(rank_df, path_points)
      }
      # ----- GLOBAL VIOLATIONS -----
      viol_row <- count_order_violations_global(A, z_hat)
      path_viol <- write_violations_csv(base, ds, fm, K, viol_row)
      
      # ----- METRICS ROW -----
      model_row <- tibble(
        dataset = ds, fit_model = fm, K = K,
        slope_logit_mean = rho_shape_tbl$slope_logit_mean,
        slope_logit_lo   = rho_shape_tbl$slope_logit_lo,
        slope_logit_hi   = rho_shape_tbl$slope_logit_hi,
        concavity_logit_mean = rho_shape_tbl$concavity_logit_mean,
        concavity_logit_lo   = rho_shape_tbl$concavity_logit_lo,
        concavity_logit_hi   = rho_shape_tbl$concavity_logit_hi,
        early_late_ratio_mean = rho_shape_tbl$early_late_ratio_mean,
        early_late_ratio_lo   = rho_shape_tbl$early_late_ratio_lo,
        early_late_ratio_hi   = rho_shape_tbl$early_late_ratio_hi,
        auc_rho_mean = rho_shape_tbl$auc_rho_mean,
        auc_rho_lo   = rho_shape_tbl$auc_rho_lo,
        auc_rho_hi   = rho_shape_tbl$auc_rho_hi,
        viol_count = viol_row$viol_count,
        viol_rate  = viol_row$viol_rate
      )
      
      tibble(
        dataset = ds, fit_model = fm, K = K,
        path_adjacency            = paste0(base, "_adjacency.png"),
        path_success              = paste0(base, "_success.png"),
        path_similarity           = paste0(base, "_similarity.png"),
        path_rho_d                = path_rho_d,
        path_rank_vs_rank         = path_rank_vs_rank,
        path_rank_vs_rank_simple  = path_rank_vs_rank_simple,
        path_rank_vs_rank_points  = path_points,
        path_violations           = path_viol,
        note = NA_character_,
        metrics = list(model_row)
      )
    }, error = function(e) {
      warning(glue("Failed {ds}|{fm}|K={K}: {conditionMessage(e)}"))
      tibble(dataset = ds, fit_model = fm, K = K,
             path_adjacency = NA_character_, path_success = NA_character_, path_similarity = NA_character_,
             path_rho_d = NA_character_, path_rank_vs_rank = NA_character_,
             path_rank_vs_rank_simple = NA_character_, path_rank_vs_rank_points = NA_character_,
             path_violations = NA_character_,
             note = conditionMessage(e), metrics = list(NULL))
    })
  })
  
  idx <- dplyr::bind_rows(rows)
  metric_tbl <- idx$metrics[!vapply(idx$metrics, is.null, logical(1))]
  metric_tbl <- if (length(metric_tbl)) dplyr::bind_rows(metric_tbl) else tibble()
  
  out_index <- file.path(OUT_DIR, "index_plots.csv")
  readr::write_csv(dplyr::select(idx, -metrics), out_index)
  if (nrow(metric_tbl)) {
    out_metrics <- file.path(OUT_DIR, "model_insights.csv")
    readr::write_csv(metric_tbl, out_metrics)
    message("Saved model insights: ", out_metrics)
  }
  message("Saved index: ", out_index)
}

# ---------------- Run ----------------
OUT_DIR <- "./output/application/plots/"
run_all("best_by_model", OUT_DIR)
# run_all("best_overall", OUT_DIR)
# run_all("all", OUT_DIR)
