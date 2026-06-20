#!/usr/bin/env Rscript
# =============================================================================
# build_simulation_tables_and_figures.R — One-stop simulation study rebuild
#
# Rebuilds the paper-facing simulation outputs from a bundled raw-results CSV:
#   1. ARI / VI boxplots                  → output/simulation/plots/
#   2. LaTeX tables (main + appendix)     → output/simulation/tables/
#   3. Cross-fit summary tables/CSV audits→ output/simulation/tables/
#
# Usage:
#   Rscript scripts/analysis/build_simulation_tables_and_figures.R
#   SIM_RESULTS_PATH=path/to/results.csv \
#     Rscript scripts/analysis/build_simulation_tables_and_figures.R
#
# All scripts assume working directory = project root.
# =============================================================================

cat("============================================================\n")
cat("       SIMULATION STUDY ANALYSIS PIPELINE\n")
cat("============================================================\n\n")

# ---- 1. ARI / VI Boxplots ------------------------------------------------
cat(">>> Step 1/3: ARI / VI boxplots ...\n")
tryCatch({
  source("scripts/analysis/sim_visualization.R")
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

# ---- 2. LaTeX tables (main + appendix) -----------------------------------
cat(">>> Step 2/3: LaTeX tables (main + appendix) ...\n")
tryCatch({
  source("scripts/analysis/read_process_save_sim_results.R")
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

# ---- 3. Cross-fit summary tables -------------------------------------------
cat(">>> Step 3/3: Cross-fit summary tables ...\n")
tryCatch({
  source("scripts/analysis/build_simulation_crossfit_tables.R")
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

cat("============================================================\n")
cat("  Simulation analysis complete.\n")
cat("  Plots  → output/simulation/plots/\n")
cat("  Tables → output/simulation/tables/\n")
cat("============================================================\n")
