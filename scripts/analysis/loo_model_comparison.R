#!/usr/bin/env Rscript
# ======================================================================
# LOO model comparison: pairwise elpd differences with SE
#
# For each dataset, reports best model, second-best, delta LOOIC,
# SE(delta), delta elpd per observation, and a verdict:
#   "clear winner"      |delta/SE| > 4
#   "weak edge"         2 < |delta/SE| <= 4
#   "indistinguishable" |delta/SE| <= 2
#
# Usage:
#   Rscript scripts/analysis/loo_model_comparison.R [RUN_DIR]
# If RUN_DIR is omitted, uses APP_RUN_DIR env var or finds the latest run.
# ======================================================================
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

suppressPackageStartupMessages({
  library(loo)
  library(dplyr)
  library(readr)
})

source("./helper_folder/helper.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./helper_folder/Hyper_setup.R")
source("./core/my_best_try_so_far.R")
source("./core/DCSBM_varK.R")
source("./helper_folder/sim_study_helper.R")
source("./helper_folder/transitivity_check_helper.R")

# ---- Resolve run directory ----
args <- commandArgs(trailingOnly = TRUE)
run_dir <- if (length(args) >= 1L) {
  args[1]
} else {
  Sys.getenv("APP_RUN_DIR", unset = "")
}

if (!nzchar(run_dir)) {
  raw_dir <- file.path("output", "application", "raw")
  dirs <- sort(list.dirs(raw_dir, full.names = TRUE, recursive = FALSE))
  dirs <- dirs[grepl("application_run_", basename(dirs))]
  if (length(dirs) == 0L) stop("No application runs found in ", raw_dir)
  run_dir <- dirs[length(dirs)]
}
cat("Run directory:", run_dir, "\n")

# ---- Data loader (from application_crp_ocrp.R) ----
choose_dataset <- function(dataset) {
  if (dataset == "mountain_goats") {
    fls <- list.files("./data/ShizukaMcDonald_Data", full.names = TRUE, pattern = "[.]csv$")
    ns <- vapply(fls, function(f) nrow(read.csv(f, row.names = 1)), integer(1))
    A <- as.matrix(read.csv(fls[which.max(ns)], row.names = 1, check.names = FALSE))
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
    edges$source <- as.integer(edges$source); edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L; edges$target <- edges$target + 1L }
    n_nodes <- max(c(edges$source, edges$target))
    A <- matrix(0L, n_nodes, n_nodes, dimnames = list(as.character(seq_len(n_nodes)), as.character(seq_len(n_nodes))))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]; j <- edges$target[r]; w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w }
    diag(A) <- 0L
  } else if (dataset == "moreno_sheep") {
    edges <- read.csv("./data/moreno_sheep/edges.csv", comment.char = "#", strip.white = TRUE)
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source); edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L; edges$target <- edges$target + 1L }
    n_nodes <- max(c(edges$source, edges$target))
    A <- matrix(0L, n_nodes, n_nodes, dimnames = list(as.character(seq_len(n_nodes)), as.character(seq_len(n_nodes))))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]; j <- edges$target[r]; w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w }
    diag(A) <- 0L
  } else if (dataset == "strauss_2019b") {
    edges <- read.csv("./data/Strauss_2019b/edges.csv", comment.char = "#", strip.white = TRUE)
    names(edges)[1:2] <- c("source", "target")
    edges$source <- as.integer(edges$source); edges$target <- as.integer(edges$target)
    edges$weight <- 1L
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L; edges$target <- edges$target + 1L }
    n_nodes <- max(c(edges$source, edges$target))
    A <- matrix(0L, n_nodes, n_nodes, dimnames = list(as.character(seq_len(n_nodes)), as.character(seq_len(n_nodes))))
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]; j <- edges$target[r]; w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w }
    diag(A) <- 0L
  } else stop("Unknown dataset: ", dataset)
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  colnames(A) <- rownames(A)
  A
}

# ---- Discover datasets from fit files ----
fit_files <- list.files(run_dir, pattern = "_fit\\.rds$", full.names = TRUE)
parsed <- regmatches(basename(fit_files),
                     regexec("^(.+)_(WST|SST|DCSBM)_fit\\.rds$", basename(fit_files)))
ds_model <- do.call(rbind, lapply(parsed, function(m) {
  if (length(m) == 3) data.frame(dataset = m[2], model = m[3], stringsAsFactors = FALSE)
  else NULL
}))
datasets <- unique(ds_model$dataset)
models   <- unique(ds_model$model)
cat(sprintf("Datasets: %s\n", paste(datasets, collapse = ", ")))
cat(sprintf("Models:   %s\n", paste(models, collapse = ", ")))

# ======================================================================
# Compute pointwise LOO for each (dataset, model)
# ======================================================================
compute_loo_object <- function(A, fit, model) {
  DI <- build_dyad_index(nrow(A))
  LL <- if (model %in% c("WST", "SST")) {
    loglik_matrix_modular(A, fit, regime = model, dyad_index = DI)
  } else {
    loglik_matrix_dcsbm(A = as.matrix(A), dcsbm_out = fit, dyad_index = DI)
  }
  loo::loo(LL)
}

loo_objects <- list()   # key: "dataset__model"
for (ds in datasets) {
  A <- tryCatch(choose_dataset(ds), error = function(e) {
    cat(sprintf("  WARNING: cannot load %s: %s\n", ds, e$message)); NULL
  })
  if (is.null(A)) next

  n_obs <- nrow(A) * (nrow(A) - 1)  # directed dyad count

  for (mod in models) {
    fit_path <- file.path(run_dir, sprintf("%s_%s_fit.rds", ds, mod))
    if (!file.exists(fit_path)) next
    cat(sprintf("  Computing LOO: %s / %s ... ", ds, mod))
    fit <- readRDS(fit_path)
    loo_obj <- tryCatch(compute_loo_object(A, fit, mod), error = function(e) {
      cat(sprintf("FAILED: %s\n", e$message)); NULL
    })
    if (!is.null(loo_obj)) {
      key <- paste0(ds, "__", mod)
      loo_objects[[key]] <- loo_obj
      cat("OK\n")
    }
    rm(fit); gc(verbose = FALSE)
  }
}

# ======================================================================
# Pairwise comparison table
# ======================================================================
comparison_rows <- list()

for (ds in datasets) {
  ds_keys <- paste0(ds, "__", models)
  ds_keys <- ds_keys[ds_keys %in% names(loo_objects)]
  if (length(ds_keys) < 2L) next

  ds_loos <- loo_objects[ds_keys]
  names(ds_loos) <- sub(paste0("^", ds, "__"), "", names(ds_loos))

  # loo_compare returns rows sorted best to worst
  comp <- loo::loo_compare(ds_loos)

  best_model    <- rownames(comp)[1]
  second_model  <- rownames(comp)[2]

  # Row 2 has delta elpd and SE relative to best
  delta_elpd    <- comp[2, "elpd_diff"]     # negative (worse than best)
  se_delta      <- comp[2, "se_diff"]

  delta_looic   <- -2 * delta_elpd          # positive = best wins by this much
  se_looic      <- 2 * se_delta

  # Number of observations (directed dyads)
  n_pw <- nrow(loo_objects[[ds_keys[1]]]$pointwise)
  delta_elpd_per_obs <- delta_elpd / n_pw

  # Verdict based on |delta/SE|
  z_ratio <- abs(delta_elpd / se_delta)
  verdict <- if (z_ratio > 4) {
    "clear winner"
  } else if (z_ratio > 2) {
    "weak edge"
  } else {
    "indistinguishable"
  }

  # Also compute third-best if 3+ models
  third_model <- if (nrow(comp) >= 3) rownames(comp)[3] else NA_character_
  delta_elpd_3rd <- if (nrow(comp) >= 3) comp[3, "elpd_diff"] else NA_real_
  se_delta_3rd   <- if (nrow(comp) >= 3) comp[3, "se_diff"]   else NA_real_

  comparison_rows[[ds]] <- data.frame(
    dataset         = ds,
    best_model      = best_model,
    second_model    = second_model,
    delta_elpd      = round(delta_elpd, 2),
    se_delta_elpd   = round(se_delta, 2),
    delta_looic     = round(delta_looic, 2),
    se_delta_looic  = round(se_looic, 2),
    z_ratio         = round(z_ratio, 2),
    delta_elpd_per_obs = round(delta_elpd_per_obs, 4),
    n_obs           = n_pw,
    verdict         = verdict,
    third_model     = third_model,
    delta_elpd_3rd  = round(delta_elpd_3rd, 2),
    se_delta_3rd    = round(se_delta_3rd, 2),
    stringsAsFactors = FALSE
  )
}

comp_table <- do.call(rbind, comparison_rows)
rownames(comp_table) <- NULL

# ---- Print ----
cat("\n================================================================\n")
cat("MODEL COMPARISON (best vs second-best, loo::loo_compare)\n")
cat("================================================================\n\n")

for (i in seq_len(nrow(comp_table))) {
  r <- comp_table[i, ]
  cat(sprintf("  %s:\n", r$dataset))
  cat(sprintf("    1st: %-6s   2nd: %-6s   delta_elpd = %.1f (SE = %.1f)  |z| = %.1f  =>  %s\n",
              r$best_model, r$second_model,
              r$delta_elpd, r$se_delta_elpd, r$z_ratio, r$verdict))
  cat(sprintf("    delta_LOOIC = %.1f (SE = %.1f)   delta_elpd/obs = %.4f   n_obs = %d\n",
              r$delta_looic, r$se_delta_looic, r$delta_elpd_per_obs, r$n_obs))
  if (!is.na(r$third_model)) {
    cat(sprintf("    3rd: %-6s   delta_elpd = %.1f (SE = %.1f)\n",
                r$third_model, r$delta_elpd_3rd, r$se_delta_3rd))
  }
  cat("\n")
}

# ---- Save ----
out_csv <- file.path(run_dir, "model_comparison_loo.csv")
readr::write_csv(comp_table, out_csv)
cat(sprintf("Saved: %s\n", out_csv))

out_rds <- file.path(run_dir, "model_comparison_loo.rds")
saveRDS(comp_table, out_rds)
cat(sprintf("Saved: %s\n", out_rds))
