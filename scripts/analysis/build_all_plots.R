#!/usr/bin/env Rscript
# =============================================================================
# build_all_plots.R -- DEFAULT one-shot regeneration of every plot in the repo.
#
# Calls, in order:
#   0. scripts/analysis/build_post_processing.R    (canonical cube, fast cache)
#   1. scripts/analysis/build_paper_plots.R        (paper figures, ~2 min)
#   2. scripts/analysis/build_exploratory_plots.R  (per-dataset library, slower)
#
# All children honour APP_RUN_DIR / output/paper/figures/current the same way.
# The cube under output/posterior_post_processing/<run_id>/ is the single
# source of truth for canonical z_hat / violation_rate / K_hat — every plot
# script reads from it via scripts/analysis/post_processing_helpers.R.
#
# Usage (typical):
#   APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#     Rscript scripts/analysis/build_all_plots.R
#
# Or, when output/paper/figures/current already points at the run you want:
#   Rscript scripts/analysis/build_all_plots.R
#
# Skip the exploratory half:   PAPER_ONLY=1   Rscript .../build_all_plots.R
# Skip the paper half:         EXPLORATORY_ONLY=1   Rscript .../build_all_plots.R
# Force cube rebuild:          FORCE_REBUILD=1
# Skip individual exploratory stages: BUILD_EXPL_SKIP="osbm,newviz"
# =============================================================================

t_total <- Sys.time()

stage <- function(label, script_path) {
  cat("############################################################\n")
  cat("# ", label, "\n", sep = "")
  cat("# script: ", script_path, "\n", sep = "")
  cat("############################################################\n")
  t0 <- Sys.time()
  ok <- tryCatch({
    source(script_path, local = TRUE); TRUE
  }, error = function(e) {
    cat("STAGE FAILED: ", conditionMessage(e), "\n", sep = ""); FALSE
  })
  cat(sprintf("# stage [%s, %.1fs]\n\n",
              if (ok) "ok" else "FAIL",
              as.numeric(Sys.time() - t0, units = "secs")))
  invisible(ok)
}

# Stage 0: cube — must come first; downstream plot scripts depend on it.
stage("Stage 0: post-processing cube (canonical z_hat / violations / K)",
      "scripts/analysis/build_post_processing.R")

if (!identical(tolower(Sys.getenv("EXPLORATORY_ONLY", "")), "1") &&
    !identical(tolower(Sys.getenv("EXPLORATORY_ONLY", "")), "true")) {
  stage("Stage A: paper figures", "scripts/analysis/build_paper_plots.R")
} else {
  cat("# Skipping paper stage (EXPLORATORY_ONLY set).\n\n")
}

if (!identical(tolower(Sys.getenv("PAPER_ONLY", "")), "1") &&
    !identical(tolower(Sys.getenv("PAPER_ONLY", "")), "true")) {
  stage("Stage B: exploratory plot library",
        "scripts/analysis/build_exploratory_plots.R")
} else {
  cat("# Skipping exploratory stage (PAPER_ONLY set).\n\n")
}

cat(sprintf("All plot pipelines complete in %.1f s.\n",
            as.numeric(Sys.time() - t_total, units = "secs")))
