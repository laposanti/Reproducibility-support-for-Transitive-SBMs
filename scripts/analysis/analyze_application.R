#!/usr/bin/env Rscript
# =============================================================================
# analyze_application.R — One-stop application (empirical) analysis
#
# Generates all application-related outputs:
#   1. Model selection + hierarchy tables   → output/application/tables/
#   2. Adjacency / success / PSM plots      → output/application/plots/<dataset>/
#   3. OSBM network/heatmap visualizations  → output/application/plots/<dataset>/
#   4. Posterior K uncertainty plots + PSM   → output/application/plots/<dataset>/
#   5. Empirical LaTeX tables               → output/application/tables/
#
# Plots are grouped by dataset into subfolders under output/application/plots/.
#
# Usage:
#   # Analyse a specific run (the only safe default):
#   Rscript scripts/analysis/analyze_application.R --run-id=application_run_20260414_104327
#   APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#     Rscript scripts/analysis/analyze_application.R
#
#   # If neither is given, the script falls back to whichever run is currently
#   # blessed via the symlink output/paper/tables/current ; if no run is
#   # blessed, it errors out instead of guessing the latest folder.
#
#   # Bless a run as the paper-facing "current" one (idempotent):
#   APP_BLESS_RUN=1 APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#     Rscript scripts/analysis/analyze_application.R
#
# All scripts assume working directory = project root.
# =============================================================================

cat("============================================================\n")
cat("       APPLICATION ANALYSIS PIPELINE\n")
cat("============================================================\n\n")

# ---- Central run directory ------------------------------------------------
# The chosen run is explicit: either pass APP_RUN_DIR=<path>, or
# --run-id=<id> on the command line. We never silently pick "the latest";
# that is what produced the earlier mismatches between tables and plots.
parse_cli_run_id <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  rid <- NULL
  for (a in args) {
    m <- regmatches(a, regexec("^--run-id=(.+)$", a))[[1]]
    if (length(m) == 2L) rid <- m[[2]]
  }
  rid
}
cli_run_id <- parse_cli_run_id()

if (nzchar(Sys.getenv("APP_RUN_DIR", unset = ""))) {
  RUN_DIR <- Sys.getenv("APP_RUN_DIR")
} else if (!is.null(cli_run_id)) {
  RUN_DIR <- file.path("output", "application", "raw", cli_run_id)
} else {
  raw_base  <- "./output/application/raw"
  current_p <- file.path("output", "paper", "tables", "current")
  if (file.exists(current_p)) {
    # Resolve the run blessed as "current" in output/paper/tables/current
    target <- Sys.readlink(current_p)
    if (!nzchar(target)) target <- normalizePath(current_p, mustWork = FALSE)
    RUN_DIR <- file.path(raw_base, basename(target))
    cat(sprintf("Using blessed run from output/paper/tables/current -> %s\n",
                basename(RUN_DIR)))
  } else {
    stop("No run selected. Set APP_RUN_DIR=<path> or pass --run-id=<id>, ",
         "or bless a run by linking output/paper/tables/current.")
  }
}
stopifnot(dir.exists(RUN_DIR))

# Look for either the new all_results.csv or old summary CSV
all_res_file <- file.path(RUN_DIR, "all_results.csv")
sum_files <- list.files(RUN_DIR, pattern = "^applications_results_summary_.*\\.csv$",
                        full.names = TRUE)
if (file.exists(all_res_file)) {
  Sys.setenv(APP_SUMMARY_PATH = all_res_file)
} else if (length(sum_files) >= 1) {
  Sys.setenv(APP_SUMMARY_PATH = sum_files[1])
} else {
  stop("No results CSV found in run dir: ", RUN_DIR)
}
Sys.setenv(APP_RUN_DIR = RUN_DIR)
cat("Run directory:", RUN_DIR, "\n")
cat("Results CSV:  ", Sys.getenv("APP_SUMMARY_PATH"), "\n\n")

# ---- 1. Model selection + hierarchy diagnostics tables --------------------
cat(">>> Step 1/5: Model selection + hierarchy tables ...\n")
tryCatch({
  # app_var_analyze_new.R was moved to old_repo/deprecated_scripts/. Its
  # role is now covered by build_post_processing.R (the canonical cube) +
  # audit_hierarchy_synopsis.R (heavy block-level audits).
  legacy <- "old_repo/deprecated_scripts/app_var_analyze_new.R"
  if (file.exists(legacy)) {
    source(legacy)
  } else {
    cat("    SKIPPED: legacy script not found.\n")
  }
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

# ---- 2. Application adjacency / success / similarity plots ---------------
cat(">>> Step 2/5: Adjacency / success / similarity plots ...\n")
tryCatch({
  source("scripts/analysis/plotting_script.R")
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

# ---- 3. OSBM network / heatmap visualizations ----------------------------
cat(">>> Step 3/5: OSBM network / heatmap visualizations ...\n")
tryCatch({
  source("scripts/analysis/osbm_visualization.R")
  run_all_osbm_visualizations()
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

# ---- 3b. New experimental visualizations ----------------------------------
cat(">>> Step 3b/5: New experimental visualizations ...\n")
tryCatch({
  source("scripts/analysis/new_visualizations.R")
  run_all_new_visualizations()
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

# ---- 4. Posterior K uncertainty + PSM heatmaps ----------------------------
cat(">>> Step 4/5: Posterior K uncertainty plots ...\n")
tryCatch({
  source("scripts/analysis/posterior_K_uq.R")
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

# ---- 5. Empirical block/item LaTeX tables ---------------------------------
cat(">>> Step 5/5: Empirical LaTeX tables ...\n")
tryCatch({
  source("scripts/analysis/read_process_save_app_results.R")
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

# ---- 6. Paper-facing model-selection table (versioned) -------------------
# Writes output/paper/tables/<run_id>/model_selection_paper.tex (no verdict
# column; just K [95%], LOOIC, dELPD, SE(dELPD), |z|).
cat(">>> Step 6/6: Paper LOO table (output/paper/tables/<run_id>/) ...\n")
tryCatch({
  Sys.setenv(APP_PAPER_TABLES_DIR = file.path(
    "output", "paper", "tables", basename(RUN_DIR)
  ))
  source("scripts/analysis/build_paper_loo_table.R", local = TRUE)
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

# ---- 7. (Optional) bless this run as the paper "current" ------------------
# Set APP_BLESS_RUN=1 to update output/paper/{tables,figures}/current
# symlinks to point at the run we just analysed. Off by default so that
# experimental runs do not accidentally rewire the manuscript.
if (identical(tolower(Sys.getenv("APP_BLESS_RUN", unset = "")), "1") ||
    identical(tolower(Sys.getenv("APP_BLESS_RUN", unset = "")), "true")) {
  for (kind in c("tables", "figures")) {
    base <- file.path("output", "paper", kind)
    if (!dir.exists(file.path(base, basename(RUN_DIR)))) {
      cat(sprintf("    Skipping bless of %s: %s/%s does not exist yet.\n",
                  kind, base, basename(RUN_DIR)))
      next
    }
    cur <- file.path(base, "current")
    file.remove(cur) |> suppressWarnings()
    file.symlink(basename(RUN_DIR), cur)
    cat(sprintf("    Blessed: %s -> %s\n", cur, basename(RUN_DIR)))
  }
}

cat("============================================================\n")
cat("  Application analysis complete.\n")
cat("  Plots  → output/application/plots/<dataset>/\n")
cat("  Tables → output/application/tables/\n")
cat("  OSBM  → output/application/plots/osbm_visualizations/<dataset>/\n")
cat("  New   → output/application/plots/new_visualizations/<dataset>/\n")
cat("============================================================\n")
