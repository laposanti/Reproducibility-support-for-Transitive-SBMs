#!/usr/bin/env Rscript
# Run simulation plots + tables for a specific results CSV
# and save outputs into date+run specific subfolders.

args <- commandArgs(trailingOnly = TRUE)

results_csv <- if (length(args) >= 1L && nzchar(args[[1L]])) {
  args[[1L]]
} else {
  Sys.getenv("SIM_RESULTS_PATH", "")
}

if (!nzchar(results_csv)) {
  stop("Pass results CSV as first arg or set SIM_RESULTS_PATH.", call. = FALSE)
}

results_csv <- normalizePath(results_csv, winslash = "/", mustWork = TRUE)

base_name <- basename(results_csv)
run_id <- sub(
  "^full_simulation_crossfit_(final|progress)_(Demo(Kvar|OCRPvar)_run_[0-9]{8}_[0-9]{6})\\.csv$",
  "\\2",
  base_name
)
if (identical(run_id, base_name)) {
  run_id <- sub(
    "^snapshot_after_rep1_(Demo(Kvar|OCRPvar)_run_[0-9]{8}_[0-9]{6})\\.csv$",
    "\\1",
    base_name
  )
}
if (identical(run_id, base_name)) {
  run_id <- tools::file_path_sans_ext(base_name)
}

run_date <- sub("^.*_([0-9]{8})_[0-9]{6}$", "\\1", run_id)
if (!grepl("^[0-9]{8}$", run_date)) {
  run_date <- format(Sys.Date(), "%Y%m%d")
}

run_date_iso <- paste0(substr(run_date, 1, 4), "-", substr(run_date, 5, 6), "-", substr(run_date, 7, 8))
tag <- paste0(run_date_iso, "_", run_id)

plots_dir <- file.path("output", "simulation", "plots", tag)
tables_dir <- file.path("output", "simulation", "tables", tag)

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tables_dir, "appendix"), recursive = TRUE, showWarnings = FALSE)

Sys.setenv(SIM_RESULTS_PATH = results_csv)
Sys.setenv(RESULTS_PATH = results_csv)
Sys.setenv(COMPACT_RESULTS_PATH = results_csv)
Sys.setenv(SIM_PLOTS_OUTPUT_DIR = plots_dir)
Sys.setenv(SIM_TABLES_OUTPUT_DIR = tables_dir)

cat("============================================================\n")
cat(" Simulation analysis for run:\n")
cat("   Results: ", results_csv, "\n", sep = "")
cat("   Plots:   ", plots_dir, "\n", sep = "")
cat("   Tables:  ", tables_dir, "\n", sep = "")
cat("============================================================\n\n")

source("scripts/analysis/sim_visualization.R")
source("scripts/analysis/build_ocrp_snapshot_tables.R")
pointwise_dir <- file.path(dirname(results_csv), "loo_pointwise")
if (dir.exists(pointwise_dir) &&
    length(list.files(pointwise_dir, pattern = "\\.csv$", full.names = TRUE)) > 0L) {
  source("scripts/analysis/analyze_loo_pointwise_dyads.R")
}
source("scripts/analysis/read_process_save_sim_results.R")

cat("\nDone.\n")
