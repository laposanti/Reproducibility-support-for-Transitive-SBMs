#!/usr/bin/env Rscript
# =============================================================================
# audit_ordered_model_diagnostics.R
#
# Single-file audit for hierarchy diagnostics in the application study.
#
# It recomputes, from the saved fits:
#   1. canonical partition-level backward mass and backward-mass rate;
#   2. same-K posterior empirical block WST/SST conformity;
#   3. same-K DC-SBM model-implied rho conformity;
#   4. corrected DC-SBM encompassing Bayes factors for WST/SST regions;
#   5. comparisons against the current hierarchy_diagnostics_overview.csv;
#   6. a compact research readout combining structure diagnostics with LOO.
#
# Run from the repository root:
#   APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#   EBF_PRIOR_MC=20000 \
#   Rscript scripts/analysis/audit_ordered_model_diagnostics.R
# =============================================================================

suppressPackageStartupMessages({
  source("scripts/analysis/osbm_visualization.R", chdir = FALSE)
  source("helper_folder/transitivity_check_helper.R", chdir = FALSE)
})

RUN_DIR <- Sys.getenv(
  "APP_RUN_DIR",
  unset = "output/application/raw/application_run_20260414_104327"
)
if (!dir.exists(RUN_DIR)) {
  stop("APP_RUN_DIR does not exist: ", RUN_DIR)
}

RUN_ID <- basename(normalizePath(RUN_DIR))
OUT_DIR <- Sys.getenv(
  "ORDER_DIAG_OUT_DIR",
  unset = file.path("output", "application", "tables")
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

PAPER_TABLE_DIR <- file.path("output", "paper", "tables", RUN_ID)
OVERVIEW_PATH <- file.path("output", "application", "tables",
                           "hierarchy_diagnostics_overview.csv")
ALL_RESULTS_PATH <- file.path(RUN_DIR, "all_results.csv")
MODEL_SELECTION_PATH <- file.path(PAPER_TABLE_DIR, "model_selection_paper.csv")

EBF_PRIOR_MC <- as.integer(Sys.getenv("EBF_PRIOR_MC", unset = "20000"))
EBF_SEED <- as.integer(Sys.getenv("EBF_SEED", unset = "20260516"))
ALPHA <- as.numeric(Sys.getenv("ORDER_DIAG_ALPHA", unset = "0.5"))
COMPARE_TOL <- as.numeric(Sys.getenv("ORDER_DIAG_COMPARE_TOL", unset = "1e-8"))

DATASETS <- c(
  "moreno_sheep", "strauss_2019b", "mountain_goats",
  "citations_data", "macaques_data", "high_school"
)
MODELS <- c("WST", "SST", "DCSBM")

DATASET_LABEL <- c(
  moreno_sheep = "Bighorn sheep",
  strauss_2019b = "Spotted hyenas",
  mountain_goats = "Mountain goats",
  citations_data = "Stat. journals",
  macaques_data = "Japanese macaques",
  high_school = "High school"
)

fmt <- function(value, digits = 3) {
  if (is.na(value) || !is.finite(value)) return("NA")
  formatC(value, format = "f", digits = digits)
}

fmt_pct <- function(value, digits = 1) {
  if (is.na(value) || !is.finite(value)) return("NA")
  paste0(formatC(100 * value, format = "f", digits = digits), "%")
}

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_csv <- function(data_frame, path) {
  utils::write.csv(data_frame, path, row.names = FALSE, na = "")
  cat("Wrote:", path, "\n")
}

summarise_numeric <- function(values) {
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(c(mean = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  c(
    mean = mean(values),
    lo = as.numeric(stats::quantile(values, 0.025, names = FALSE, type = 7)),
    hi = as.numeric(stats::quantile(values, 0.975, names = FALSE, type = 7))
  )
}

get_draw_count <- function(z_object) {
  if (is.matrix(z_object)) return(nrow(z_object))
  if (is.list(z_object)) return(length(z_object))
  stop("fit$z must be a matrix or list.")
}

get_z_draw <- function(z_object, draw_index) {
  if (is.matrix(z_object)) return(z_object[draw_index, ])
  z_object[[draw_index]]
}

occupied_block_count <- function(z_draw) {
  length(unique(z_draw[!is.na(z_draw)]))
}

load_adjacency <- function(dataset_id) {
  ppc_path <- file.path("output", "application", "ppc", dataset_id,
                        paste0(dataset_id, "_A_obs.rds"))
  adjacency <- if (file.exists(ppc_path)) {
    readRDS(ppc_path)
  } else {
    choose_dataset_local(dataset_id)
  }
  if (inherits(adjacency, "Matrix")) adjacency <- as.matrix(adjacency)
  as.matrix(adjacency)
}

backward_mass_stats <- function(adjacency, z_ranked) {
  edge_index <- which(adjacency > 0, arr.ind = TRUE)
  if (!nrow(edge_index)) {
    return(c(rate = NA_real_, backward_mass = 0, cross_mass = 0))
  }
  from_rank <- z_ranked[edge_index[, 1]]
  to_rank <- z_ranked[edge_index[, 2]]
  weights <- adjacency[edge_index]
  cross_edge <- from_rank != to_rank & !is.na(from_rank) & !is.na(to_rank)
  backward_edge <- cross_edge & from_rank > to_rank
  cross_mass <- sum(weights[cross_edge], na.rm = TRUE)
  backward_mass <- sum(weights[backward_edge], na.rm = TRUE)
  c(
    rate = if (cross_mass > 0) backward_mass / cross_mass else NA_real_,
    backward_mass = backward_mass,
    cross_mass = cross_mass
  )
}

ordered_rho <- function(rho_matrix,
                        method_order = c("mean", "identity"),
                        eps = 1e-8) {
  method_order <- match.arg(method_order)
  rho_matrix <- pmin(pmax(as.matrix(rho_matrix), eps), 1 - eps)
  diag(rho_matrix) <- 0.5
  if (method_order == "identity") return(rho_matrix)
  score <- rowMeans(rho_matrix)
  rho_matrix[order(-score, seq_along(score)), order(-score, seq_along(score)),
             drop = FALSE]
}

score_ordered_triples <- function(rho_matrix, eps = 1e-8) {
  block_count <- nrow(rho_matrix)
  if (block_count < 3L) {
    return(list(thetaW = NA_real_, thetaS = NA_real_, thetaW_all = NA_real_,
                thetaS_all = NA_real_, prem = 0L, coverage = 0))
  }

  rho_matrix <- pmin(pmax(rho_matrix, eps), 1 - eps)
  diag(rho_matrix) <- 0.5
  triples <- utils::combn(block_count, 3L)
  total <- ncol(triples)
  premise_count <- 0L
  ok_wst_given_premise <- 0L
  ok_sst_given_premise <- 0L

  for (triple_index in seq_len(total)) {
    block_a <- triples[1, triple_index]
    block_b <- triples[2, triple_index]
    block_c <- triples[3, triple_index]
    rho_ab <- rho_matrix[block_a, block_b]
    rho_bc <- rho_matrix[block_b, block_c]
    rho_ac <- rho_matrix[block_a, block_c]
    if (rho_ab >= 0.5 && rho_bc >= 0.5) {
      premise_count <- premise_count + 1L
      if (rho_ac >= 0.5) ok_wst_given_premise <- ok_wst_given_premise + 1L
      if (rho_ac >= max(rho_ab, rho_bc)) {
        ok_sst_given_premise <- ok_sst_given_premise + 1L
      }
    }
  }

  list(
    thetaW = if (premise_count > 0L) ok_wst_given_premise / premise_count else NA_real_,
    thetaS = if (premise_count > 0L) ok_sst_given_premise / premise_count else NA_real_,
    thetaW_all = ok_wst_given_premise / total,
    thetaS_all = ok_sst_given_premise / total,
    prem = premise_count,
    coverage = premise_count / total
  )
}

exact_region_check <- function(rho_matrix,
                               method_order = c("mean", "identity"),
                               tol = 1e-12) {
  rho_ordered <- ordered_rho(rho_matrix, method_order = method_order)
  block_count <- nrow(rho_ordered)
  if (block_count < 2L) return(list(wst = NA, sst = NA))
  wst_ok <- all(rho_ordered[upper.tri(rho_ordered)] + tol >= 0.5)
  if (block_count < 3L) return(list(wst = wst_ok, sst = NA))

  sst_ok <- wst_ok
  for (block_a in seq_len(block_count - 2L)) {
    for (block_b in (block_a + 1L):(block_count - 1L)) {
      for (block_c in (block_b + 1L):block_count) {
        if (!(rho_ordered[block_a, block_c] + tol >=
              max(rho_ordered[block_a, block_b], rho_ordered[block_b, block_c]))) {
          sst_ok <- FALSE
          break
        }
      }
      if (!sst_ok) break
    }
    if (!sst_ok) break
  }
  list(wst = wst_ok, sst = sst_ok)
}

conformity_from_rho <- function(rho_matrix,
                                method_order = c("mean", "identity")) {
  rho_ordered <- ordered_rho(rho_matrix, method_order = method_order)
  score_ordered_triples(rho_ordered)
}

block_flow_matrix <- function(adjacency, z_draw) {
  valid_node <- !is.na(z_draw)
  adjacency <- adjacency[valid_node, valid_node, drop = FALSE]
  compact_z <- match(z_draw[valid_node], sort(unique(z_draw[valid_node])))
  row_aggregate <- rowsum(adjacency, group = compact_z, reorder = TRUE)
  t(rowsum(t(row_aggregate), group = compact_z, reorder = TRUE))
}

empirical_conformity_one_draw <- function(adjacency, z_draw,
                                          method_order = c("mean", "identity"),
                                          alpha = 0.5) {
  method_order <- match.arg(method_order)
  if (occupied_block_count(z_draw) < 3L) {
    return(list(thetaW = NA_real_, thetaS = NA_real_, thetaW_all = NA_real_,
                thetaS_all = NA_real_, prem = 0L, coverage = 0))
  }
  flow_matrix <- block_flow_matrix(adjacency, z_draw)
  total_matrix <- flow_matrix + t(flow_matrix)
  rho_hat <- (flow_matrix + alpha) / pmax(total_matrix + 2 * alpha,
                                          .Machine$double.eps)
  diag(rho_hat) <- 0.5
  conformity_from_rho(rho_hat, method_order = method_order)
}

empirical_conformity_chain <- function(adjacency, z_object, target_block_count,
                                       model_id) {
  draw_count <- get_draw_count(z_object)
  method_order <- if (model_id == "DCSBM") "mean" else "identity"
  same_k <- rep(FALSE, draw_count)
  theta_w <- theta_s <- coverage <- prem <- rep(NA_real_, draw_count)

  for (draw_index in seq_len(draw_count)) {
    z_draw <- get_z_draw(z_object, draw_index)
    if (occupied_block_count(z_draw) != target_block_count) next
    same_k[draw_index] <- TRUE
    scores <- empirical_conformity_one_draw(
      adjacency = adjacency,
      z_draw = z_draw,
      method_order = method_order,
      alpha = ALPHA
    )
    theta_w[draw_index] <- scores$thetaW
    theta_s[draw_index] <- scores$thetaS
    coverage[draw_index] <- scores$coverage
    prem[draw_index] <- scores$prem
  }

  theta_w_summary <- summarise_numeric(theta_w)
  theta_s_summary <- summarise_numeric(theta_s)
  c(
    thetaW_block_emp_mean = theta_w_summary[["mean"]],
    thetaW_block_emp_lo = theta_w_summary[["lo"]],
    thetaW_block_emp_hi = theta_w_summary[["hi"]],
    thetaS_block_emp_mean = theta_s_summary[["mean"]],
    thetaS_block_emp_lo = theta_s_summary[["lo"]],
    thetaS_block_emp_hi = theta_s_summary[["hi"]],
    prem_block_emp_avg = mean(prem, na.rm = TRUE),
    coverage_block_emp_avg = mean(coverage, na.rm = TRUE),
    n_draws_K_match_emp = sum(same_k)
  )
}

get_theta_row <- function(theta_object, draw_index, node_count) {
  if (is.matrix(theta_object)) {
    if (ncol(theta_object) != node_count) {
      stop("theta matrix has ", ncol(theta_object), " columns, expected ", node_count)
    }
    return(theta_object[draw_index, ])
  }
  if (is.list(theta_object)) return(theta_object[[draw_index]])
  stop("theta must be a matrix or list.")
}

get_lambda_draw <- function(lambda_object, draw_index) {
  if (is.list(lambda_object)) return(lambda_object[[draw_index]])
  if (length(dim(lambda_object)) == 3L) return(lambda_object[draw_index, , ])
  stop("lambda must be a list or 3D array.")
}

align_lambda_to_occupied <- function(lambda_draw, occupied_labels, block_count) {
  lambda_draw <- as.matrix(lambda_draw)
  if (all(occupied_labels <= nrow(lambda_draw)) &&
      all(occupied_labels <= ncol(lambda_draw))) {
    return(lambda_draw[occupied_labels, occupied_labels, drop = FALSE])
  }
  if (all(dim(lambda_draw) == c(block_count, block_count))) return(lambda_draw)
  stop("Cannot align lambda draw to occupied labels.")
}

dcsbm_model_rho_chain <- function(fit, target_block_count, node_count) {
  z_object <- fit$z
  draw_count <- get_draw_count(z_object)
  theta_out_object <- if (!is.null(fit$theta_out)) fit$theta_out else fit$theta
  theta_in_object <- if (!is.null(fit$theta_in)) fit$theta_in else fit$theta
  lambda_object <- fit$lambda

  theta_w <- theta_s <- coverage <- prem <- rep(NA_real_, draw_count)
  region_w <- region_s <- bar_rho <- rep(NA_real_, draw_count)
  same_k <- rep(FALSE, draw_count)

  for (draw_index in seq_len(draw_count)) {
    z_draw <- get_z_draw(z_object, draw_index)
    occupied_labels <- sort(unique(z_draw[!is.na(z_draw)]))
    block_count <- length(occupied_labels)
    if (block_count < 3L || block_count != target_block_count) next
    same_k[draw_index] <- TRUE
    z_compact <- match(z_draw, occupied_labels)
    lambda_draw <- align_lambda_to_occupied(
      get_lambda_draw(lambda_object, draw_index),
      occupied_labels = occupied_labels,
      block_count = block_count
    )
    rho_matrix <- rho_from_dcsbm_draw(
      theta_out_draw = get_theta_row(theta_out_object, draw_index, node_count),
      theta_in_draw = get_theta_row(theta_in_object, draw_index, node_count),
      lambda_draw = lambda_draw,
      z_hat = z_compact,
      K = block_count,
      use_propensity = FALSE
    )

    conformity <- conformity_from_rho(rho_matrix, method_order = "mean")
    exact <- exact_region_check(rho_matrix, method_order = "mean")
    rho_ordered <- ordered_rho(rho_matrix, method_order = "mean")

    theta_w[draw_index] <- conformity$thetaW
    theta_s[draw_index] <- conformity$thetaS
    coverage[draw_index] <- conformity$coverage
    prem[draw_index] <- conformity$prem
    region_w[draw_index] <- as.numeric(isTRUE(exact$wst))
    region_s[draw_index] <- as.numeric(isTRUE(exact$sst))
    bar_rho[draw_index] <- mean(rho_ordered[upper.tri(rho_ordered)], na.rm = TRUE)
  }

  theta_w_summary <- summarise_numeric(theta_w)
  theta_s_summary <- summarise_numeric(theta_s)
  bar_rho_summary <- summarise_numeric(bar_rho)
  c(
    thetaW_block_model_mean = theta_w_summary[["mean"]],
    thetaW_block_model_lo = theta_w_summary[["lo"]],
    thetaW_block_model_hi = theta_w_summary[["hi"]],
    thetaS_block_model_mean = theta_s_summary[["mean"]],
    thetaS_block_model_lo = theta_s_summary[["lo"]],
    thetaS_block_model_hi = theta_s_summary[["hi"]],
    bar_rho_model_mean = bar_rho_summary[["mean"]],
    bar_rho_model_lo = bar_rho_summary[["lo"]],
    bar_rho_model_hi = bar_rho_summary[["hi"]],
    prem_block_model_avg = mean(prem, na.rm = TRUE),
    coverage_block_model_avg = mean(coverage, na.rm = TRUE),
    p_post_wst = mean(region_w, na.rm = TRUE),
    p_post_sst = mean(region_s, na.rm = TRUE),
    n_draws_K_match_model = sum(same_k)
  )
}

simulate_prior_region <- function(block_count, n_mc, seed,
                                  a_lambda = 1, b_lambda = 1) {
  if (!is.finite(block_count) || block_count < 3L || n_mc < 1L) {
    return(c(p_prior_wst = NA_real_, p_prior_sst = NA_real_,
             prior_wst_hits = NA_real_, prior_sst_hits = NA_real_,
             prior_mcse_wst = NA_real_, prior_mcse_sst = NA_real_))
  }
  set.seed(seed + block_count * 1009L)
  wst_hits <- logical(n_mc)
  sst_hits <- logical(n_mc)
  for (draw_index in seq_len(n_mc)) {
    lambda_matrix <- matrix(
      stats::rgamma(block_count * block_count, shape = a_lambda, rate = b_lambda),
      nrow = block_count,
      ncol = block_count
    )
    rho_matrix <- lambda_matrix / (lambda_matrix + t(lambda_matrix) + 1e-12)
    diag(rho_matrix) <- 0.5
    exact <- exact_region_check(rho_matrix, method_order = "mean")
    wst_hits[draw_index] <- isTRUE(exact$wst)
    sst_hits[draw_index] <- isTRUE(exact$sst)
  }
  p_wst <- mean(wst_hits)
  p_sst <- mean(sst_hits)
  c(
    p_prior_wst = p_wst,
    p_prior_sst = p_sst,
    prior_wst_hits = sum(wst_hits),
    prior_sst_hits = sum(sst_hits),
    prior_mcse_wst = sqrt(p_wst * (1 - p_wst) / n_mc),
    prior_mcse_sst = sqrt(p_sst * (1 - p_sst) / n_mc)
  )
}

bf_value <- function(posterior_prob, prior_prob) {
  if (!is.finite(posterior_prob) || !is.finite(prior_prob)) return(NA_real_)
  if (prior_prob == 0 && posterior_prob > 0) return(Inf)
  if (prior_prob == 0 && posterior_prob == 0) return(NA_real_)
  posterior_prob / prior_prob
}

bf_lower_bound <- function(posterior_prob, prior_hits, n_mc) {
  if (!is.finite(posterior_prob) || !is.finite(prior_hits)) return(NA_real_)
  if (prior_hits == 0 && posterior_prob > 0) return(posterior_prob * n_mc)
  NA_real_
}

compute_fit_row <- function(dataset_id, model_id, adjacency) {
  fit_path <- file.path(RUN_DIR, paste0(dataset_id, "_", model_id, "_fit.rds"))
  if (!file.exists(fit_path)) stop("Missing fit: ", fit_path)
  fit <- readRDS(fit_path)
  z_hat_object <- get_z_hat_from_draws(fit, adjacency, model = model_id)
  z_hat <- z_hat_object$z_hat
  target_block_count <- max(z_hat, na.rm = TRUE)
  backward <- backward_mass_stats(adjacency, z_hat)
  empirical <- empirical_conformity_chain(
    adjacency = adjacency,
    z_object = fit$z,
    target_block_count = target_block_count,
    model_id = model_id
  )

  row_values <- c(
    n = nrow(adjacency),
    K_hat = target_block_count,
    violation_rate_zhat = backward[["rate"]],
    violation_count_zhat = backward[["backward_mass"]],
    cross_count_zhat = backward[["cross_mass"]],
    empirical,
    thetaW_block_model_mean = NA_real_,
    thetaW_block_model_lo = NA_real_,
    thetaW_block_model_hi = NA_real_,
    thetaS_block_model_mean = NA_real_,
    thetaS_block_model_lo = NA_real_,
    thetaS_block_model_hi = NA_real_,
    bar_rho_model_mean = NA_real_,
    bar_rho_model_lo = NA_real_,
    bar_rho_model_hi = NA_real_,
    prem_block_model_avg = NA_real_,
    coverage_block_model_avg = NA_real_,
    n_draws_K_match_model = NA_real_,
    p_post_wst = NA_real_,
    p_prior_wst = NA_real_,
    bf_wst_0 = NA_real_,
    bf_wst_lower_bound = NA_real_,
    p_post_sst = NA_real_,
    p_prior_sst = NA_real_,
    bf_sst_0 = NA_real_,
    bf_sst_lower_bound = NA_real_,
    prior_wst_hits = NA_real_,
    prior_sst_hits = NA_real_,
    prior_mcse_wst = NA_real_,
    prior_mcse_sst = NA_real_,
    bf_prior_mc = NA_real_
  )

  if (model_id == "DCSBM") {
    model_diag <- dcsbm_model_rho_chain(
      fit = fit,
      target_block_count = target_block_count,
      node_count = nrow(adjacency)
    )
    prior_diag <- simulate_prior_region(
      block_count = target_block_count,
      n_mc = EBF_PRIOR_MC,
      seed = EBF_SEED
    )
    row_values[names(model_diag)] <- model_diag
    row_values[names(prior_diag)] <- prior_diag
    row_values[["bf_wst_0"]] <- bf_value(row_values[["p_post_wst"]],
                                         row_values[["p_prior_wst"]])
    row_values[["bf_sst_0"]] <- bf_value(row_values[["p_post_sst"]],
                                         row_values[["p_prior_sst"]])
    row_values[["bf_wst_lower_bound"]] <- bf_lower_bound(
      row_values[["p_post_wst"]], row_values[["prior_wst_hits"]], EBF_PRIOR_MC
    )
    row_values[["bf_sst_lower_bound"]] <- bf_lower_bound(
      row_values[["p_post_sst"]], row_values[["prior_sst_hits"]], EBF_PRIOR_MC
    )
    row_values[["bf_prior_mc"]] <- EBF_PRIOR_MC
  }

  data.frame(
    dataset = dataset_id,
    fit_model = model_id,
    as.list(row_values),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

comparison_rows_for_metric <- function(new_df, old_df, metric_name,
                                       dcsbm_only = FALSE) {
  if (is.null(old_df) || !(metric_name %in% names(old_df)) ||
      !(metric_name %in% names(new_df))) {
    return(NULL)
  }
  if (dcsbm_only) {
    new_df <- new_df[new_df$fit_model == "DCSBM", , drop = FALSE]
    old_df <- old_df[old_df$fit_model == "DCSBM", , drop = FALSE]
  }
  merged <- merge(
    new_df[, c("dataset", "fit_model", metric_name), drop = FALSE],
    old_df[, c("dataset", "fit_model", metric_name), drop = FALSE],
    by = c("dataset", "fit_model"), suffixes = c("_new", "_existing"),
    all.x = TRUE, sort = FALSE
  )
  new_name <- paste0(metric_name, "_new")
  old_name <- paste0(metric_name, "_existing")
  diff_value <- merged[[new_name]] - merged[[old_name]]
  both_na <- is.na(merged[[new_name]]) & is.na(merged[[old_name]])
  both_same_infinite <- is.infinite(merged[[new_name]]) &
    is.infinite(merged[[old_name]]) &
    sign(merged[[new_name]]) == sign(merged[[old_name]])
  finite_match <- is.finite(diff_value) & abs(diff_value) <= COMPARE_TOL
  data.frame(
    dataset = merged$dataset,
    fit_model = merged$fit_model,
    metric = metric_name,
    new_value = merged[[new_name]],
    existing_value = merged[[old_name]],
    diff = diff_value,
    matches_existing = both_na | both_same_infinite | finite_match,
    stringsAsFactors = FALSE
  )
}

make_comparison_table <- function(new_df, old_df) {
  metrics <- c(
    "K_hat", "violation_rate_zhat", "violation_count_zhat", "cross_count_zhat",
    "thetaW_block_emp_mean", "thetaS_block_emp_mean",
    "thetaW_block_model_mean", "thetaS_block_model_mean",
    "bar_rho_model_mean", "p_post_wst", "p_prior_wst", "bf_wst_0",
    "p_post_sst", "p_prior_sst", "bf_sst_0"
  )
  dcsbm_only_metrics <- c(
    "thetaW_block_model_mean", "thetaS_block_model_mean", "bar_rho_model_mean",
    "p_post_wst", "p_prior_wst", "bf_wst_0", "p_post_sst", "p_prior_sst",
    "bf_sst_0"
  )
  rows <- lapply(metrics, function(metric_name) {
    comparison_rows_for_metric(
      new_df, old_df, metric_name,
      dcsbm_only = metric_name %in% dcsbm_only_metrics
    )
  })
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

normalise_model_label <- function(label) {
  ifelse(grepl("Toeplitz|SST", label), "SST",
         ifelse(grepl("WST", label), "WST", "DCSBM"))
}

make_loo_readout <- function() {
  all_results <- read_csv_if_exists(ALL_RESULTS_PATH)
  if (is.null(all_results)) return(data.frame())
  model_selection <- read_csv_if_exists(MODEL_SELECTION_PATH)
  if (!is.null(model_selection)) {
    model_selection$fit_model <- normalise_model_label(model_selection$model)
  }

  rows <- list()
  for (dataset_id in DATASETS) {
    dataset_results <- all_results[all_results$dataset == dataset_id, , drop = FALSE]
    ordered_results <- dataset_results[dataset_results$fit_model %in% c("WST", "SST"), , drop = FALSE]
    dcsbm_result <- dataset_results[dataset_results$fit_model == "DCSBM", , drop = FALSE]
    if (!nrow(ordered_results) || !nrow(dcsbm_result)) next
    best_ordered <- ordered_results[which.max(ordered_results$elpd_loo), , drop = FALSE]
    delta_elpd <- best_ordered$elpd_loo - dcsbm_result$elpd_loo

    pair_se <- NA_real_
    if (!is.null(model_selection)) {
      selection_rows <- model_selection[model_selection$dataset == dataset_id, , drop = FALSE]
      best_ordered_row <- selection_rows[selection_rows$fit_model == best_ordered$fit_model, , drop = FALSE]
      dcsbm_row <- selection_rows[selection_rows$fit_model == "DCSBM", , drop = FALSE]
      if (nrow(best_ordered_row) && nrow(dcsbm_row)) {
        if (best_ordered_row$rank == 1L) pair_se <- dcsbm_row$se_delta_elpd
        if (dcsbm_row$rank == 1L) pair_se <- best_ordered_row$se_delta_elpd
      }
    }
    t_ratio <- if (is.finite(pair_se) && pair_se > 0) abs(delta_elpd / pair_se) else NA_real_
    conclusion <- if (!is.finite(t_ratio)) {
      "pairwise SE unavailable"
    } else if (delta_elpd > 0 && t_ratio >= 2) {
      "ordered predicts clearly better than DC-SBM"
    } else if (delta_elpd > 0) {
      "ordered predicts slightly better; difference is weak"
    } else if (delta_elpd < 0 && t_ratio >= 2) {
      "DC-SBM predicts clearly better than ordered"
    } else {
      "ordered and DC-SBM are predictively close"
    }
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = dataset_id,
      best_ordered_model = best_ordered$fit_model,
      elpd_best_ordered = best_ordered$elpd_loo,
      elpd_dcsbm = dcsbm_result$elpd_loo,
      delta_elpd_ordered_minus_dcsbm = delta_elpd,
      se_delta_elpd_ordered_vs_dcsbm = pair_se,
      abs_t_ordered_vs_dcsbm = t_ratio,
      predictive_conclusion = conclusion,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

bf_text <- function(bf_value, lower_bound) {
  if (is.finite(bf_value)) return(fmt(bf_value, 2))
  if (is.infinite(bf_value) && is.finite(lower_bound)) {
    return(paste0(">", fmt(lower_bound, 1)))
  }
  "NA"
}

structural_readout <- function(dcsbm_row) {
  wst_prob <- dcsbm_row$p_post_wst
  sst_prob <- dcsbm_row$p_post_sst
  if (!is.finite(wst_prob) || wst_prob == 0) {
    wst_text <- "no WST mass under the DC-SBM posterior"
  } else if (wst_prob >= 0.5) {
    wst_text <- "substantial WST mass under the DC-SBM posterior"
  } else if (wst_prob >= 0.1) {
    wst_text <- "partial WST mass under the DC-SBM posterior"
  } else {
    wst_text <- "weak WST mass under the DC-SBM posterior"
  }

  if (!is.finite(sst_prob) || sst_prob == 0) {
    sst_text <- "no exact SST mass"
  } else if (sst_prob >= 0.1) {
    sst_text <- "some exact SST mass"
  } else {
    sst_text <- "very little exact SST mass"
  }
  paste(wst_text, "and", sst_text)
}

write_markdown_report <- function(audit_df, comparison_df, loo_df, path) {
  dcsbm_df <- audit_df[audit_df$fit_model == "DCSBM", , drop = FALSE]
  mismatch_df <- comparison_df[!comparison_df$matches_existing %in% TRUE, , drop = FALSE]
  lines <- c(
    "# Ordered-model diagnostics audit",
    "",
    paste0("Run: `", RUN_ID, "`"),
    paste0("EBF prior simulations per dataset: `", EBF_PRIOR_MC, "`"),
    "",
    "## Definitions used here",
    "",
    "- Partition backward mass is `sum(A_ij)` over cross-block directed edges with `z_i > z_j`, divided by total cross-block mass for the percentage.",
    "- Empirical WST/SST conformity is a same-K posterior block-level conditional triple rate computed from observed block flow fractions.",
    "- DC-SBM model-implied conformity uses the same conditional triple score on each draw's model-implied directional probability matrix `rho`, after ordering blocks by mean success probability.",
    "- The EBF uses the exact ordered region: WST means all ordered upper-triangle entries are at least 0.5; SST means WST plus the monotone triple inequalities.",
    "- Predictive improvement is judged from LOO. The hierarchy diagnostics say whether the unordered posterior lies in an ordered region; they do not by themselves prove better prediction.",
    "",
    "## Comparison with current files",
    "",
    paste0("Metrics checked: `", length(unique(comparison_df$metric)), "`; mismatching rows: `", nrow(mismatch_df), "`."),
    ""
  )

  if (nrow(mismatch_df)) {
    mismatch_metrics <- paste(sort(unique(mismatch_df$metric)), collapse = ", ")
    lines <- c(lines, paste0("Mismatching metrics: ", mismatch_metrics, "."), "")
    ebf_metrics <- c("p_post_wst", "p_prior_wst", "bf_wst_0",
                     "p_post_sst", "p_prior_sst", "bf_sst_0")
    if (all(unique(mismatch_df$metric) %in% ebf_metrics)) {
      lines <- c(
        lines,
        "The partition-level backward mass, cross-block mass, empirical block conformity, and DC-SBM model-implied conformity all agree with the current overview. The remaining differences are the corrected EBF quantities.",
        ""
      )
    }
  } else {
    lines <- c(lines, "All recomputed metrics match the current overview within tolerance.", "")
  }

  lines <- c(
    lines,
    "## Dataset-level readout",
    "",
    "| Dataset | DCSBM WST post. | DCSBM SST post. | EBF WST | EBF SST | Ordered vs DC-SBM LOO | Reading |",
    "|---|---:|---:|---:|---:|---|---|"
  )

  for (dataset_id in DATASETS) {
    dcsbm_row <- dcsbm_df[dcsbm_df$dataset == dataset_id, , drop = FALSE]
    loo_row <- loo_df[loo_df$dataset == dataset_id, , drop = FALSE]
    if (!nrow(dcsbm_row) || !nrow(loo_row)) next
    loo_cell <- paste0(
      loo_row$best_ordered_model,
      " vs DC-SBM: Delta ELPD=", fmt(loo_row$delta_elpd_ordered_minus_dcsbm, 1),
      ", abs t=", fmt(loo_row$abs_t_ordered_vs_dcsbm, 2)
    )
    reading <- paste0(structural_readout(dcsbm_row), "; ", loo_row$predictive_conclusion, ".")
    table_row <- paste(
      DATASET_LABEL[[dataset_id]],
      fmt_pct(dcsbm_row$p_post_wst, 1),
      fmt_pct(dcsbm_row$p_post_sst, 1),
      bf_text(dcsbm_row$bf_wst_0, dcsbm_row$bf_wst_lower_bound),
      bf_text(dcsbm_row$bf_sst_0, dcsbm_row$bf_sst_lower_bound),
      loo_cell,
      reading,
      sep = " | "
    )
    lines <- c(lines, paste0("| ", table_row, " |"))
  }

  lines <- c(
    lines,
    "",
    "## Overall conclusion",
    "",
    "The corrected diagnostics separate two claims that should stay separate. First, a DCSBM posterior can contain ordered structure, measured by posterior mass in the WST/SST regions and the EBF. Second, an ordered model can predict better than the unordered DCSBM, measured by LOO. The first claim is structural; the second is predictive.",
    "",
    "The strongest predictive cases for ordered models are the sheep, goats, hyenas, and high school. The citations and macaques do not support a simple 'ordered beats unordered' conclusion: citations are best left to the DC-SBM predictively, and macaques retain a strong directional structure but the DC-SBM predicts better. Exact SST support inside the DCSBM posterior is generally scarce, so Toeplitz SST wins should be read mainly as useful regularisation or shape restriction unless the LOO comparison also clearly separates it from DC-SBM."
  )

  writeLines(lines, path)
  cat("Wrote:", path, "\n")
}

cat("Run id:", RUN_ID, "\n")
cat("Run dir:", RUN_DIR, "\n")
cat("Output dir:", OUT_DIR, "\n")
cat("EBF prior MC:", EBF_PRIOR_MC, "\n\n")

audit_rows <- list()
for (dataset_id in DATASETS) {
  adjacency <- load_adjacency(dataset_id)
  for (model_id in MODELS) {
    cat(sprintf("%-15s %-5s", dataset_id, model_id))
    fit_row <- compute_fit_row(dataset_id, model_id, adjacency)
    audit_rows[[length(audit_rows) + 1L]] <- fit_row
    cat(sprintf(" K=%2d  backward=%s\n",
                fit_row$K_hat, fmt_pct(fit_row$violation_rate_zhat)))
  }
}

audit_df <- do.call(rbind, audit_rows)
audit_path <- file.path(OUT_DIR, "ordered_model_diagnostics_audit.csv")
write_csv(audit_df, audit_path)

overview_df <- read_csv_if_exists(OVERVIEW_PATH)
comparison_df <- make_comparison_table(audit_df, overview_df)
comparison_path <- file.path(OUT_DIR, "ordered_model_diagnostics_comparison.csv")
write_csv(comparison_df, comparison_path)

mismatch_path <- file.path(OUT_DIR, "ordered_model_diagnostics_mismatches.csv")
mismatch_df <- comparison_df[!comparison_df$matches_existing %in% TRUE, , drop = FALSE]
write_csv(mismatch_df, mismatch_path)

loo_df <- make_loo_readout()
loo_path <- file.path(OUT_DIR, "ordered_model_vs_dcsbm_loo_readout.csv")
write_csv(loo_df, loo_path)

report_path <- file.path(OUT_DIR, "ordered_model_diagnostics_readout.md")
write_markdown_report(audit_df, comparison_df, loo_df, report_path)

cat("\nDone. Main report:", report_path, "\n")