#!/usr/bin/env Rscript
# =============================================================================
# analyze_simulation.R — One-stop simulation study analysis
#
# Generates all simulation-related outputs:
#   1. ARI / VI boxplots         → output/simulation/plots/
#   2. LaTeX tables (main+appx)  → output/simulation/tables/
#   3. Legacy partition/order tables + ARI plot (sim_var_analyze_new)
#
# Usage:
#   Rscript scripts/analysis/analyze_simulation.R                     # auto-detect latest
#   SIM_RESULTS_PATH=path/to/results.csv Rscript scripts/analysis/analyze_simulation.R
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

# ---- 3. Legacy partition/ELPD/order tables + ARI boxplot ------------------
cat(">>> Step 3/3: Partition / ELPD / order tables ...\n")
tryCatch({
  source("scripts/analysis/sim_var_analyze_new.R")
  cat("    Done.\n\n")
}, error = function(e) {
  cat("    FAILED:", conditionMessage(e), "\n\n")
})

cat("============================================================\n")
cat("  Simulation analysis complete.\n")
cat("  Plots  → output/simulation/plots/\n")
cat("  Tables → output/simulation/tables/\n")
cat("============================================================\n")
