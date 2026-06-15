#!/usr/bin/env Rscript
# =============================================================================
# build_exploratory_plots.R -- ONE entry point for the *exploratory* plot
# library that lives under output/application/plots/.
#
# This is the post-hoc visual diagnostics the user inspects manually, NOT the
# small set of figures embedded in the manuscript (those are produced by
# `Rscript scripts/analysis/build_paper_plots.R`).
#
# Outputs (per active run, resolved from APP_RUN_DIR or
# output/paper/figures/current):
#
#   plotting_script.R                  -> output/application/plots/<ds>/
#                                         {ds}_{model}_{adjacency,success,
#                                          similarity,rho_d,rank_vs_rank,
#                                          rank_vs_rank_simple}.png
#                                         (uses canonical get_z_hat_from_draws)
#
#   posterior_K_uq.R                   -> output/application/plots/<ds>/
#                                         K_posterior_pmf, K_trace, K_uq_panel
#
#   osbm_visualization.R::run_all_*    -> output/application/plots/osbm_visualizations/<ds>/
#
#   new_visualizations.R::run_all_*    -> output/application/plots/new_visualizations/<ds>/
#                                         arc_diagram, membership_heatmap,
#                                         rank_caterpillar, concentric_ring
#
# Usage:
#   APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#     Rscript scripts/analysis/build_exploratory_plots.R
#
# Or fall back to the run currently blessed via output/paper/figures/current:
#   Rscript scripts/analysis/build_exploratory_plots.R
#
# Skip a stage by setting BUILD_EXPL_SKIP="osbm,newviz" (comma list of any of
# perfit, posteriorK, osbm, newviz).
# =============================================================================

# ---- Run-directory resolution (shared with build_paper_plots.R) -------------
RUN_DIR <- Sys.getenv("APP_RUN_DIR", unset = "")
if (!nzchar(RUN_DIR)) {
  cur <- "output/paper/figures/current"
  if (file.exists(cur))
    RUN_DIR <- file.path("output/application/raw", basename(Sys.readlink(cur)))
}
if (!nzchar(RUN_DIR) || !dir.exists(RUN_DIR))
  stop("Set APP_RUN_DIR=output/application/raw/<run_id>; got '", RUN_DIR, "'.")

RUN_ID <- basename(RUN_DIR)
Sys.setenv(APP_RUN_DIR = RUN_DIR)

# Some downstream scripts also need APP_SUMMARY_PATH
all_res_file <- file.path(RUN_DIR, "all_results.csv")
sum_files <- list.files(RUN_DIR,
                        pattern = "^applications_results_summary_.*\\.csv$",
                        full.names = TRUE)
if (file.exists(all_res_file)) {
  Sys.setenv(APP_SUMMARY_PATH = all_res_file)
} else if (length(sum_files) >= 1L) {
  Sys.setenv(APP_SUMMARY_PATH = sum_files[1])
}

skip_raw  <- Sys.getenv("BUILD_EXPL_SKIP", unset = "")
SKIP <- if (nzchar(skip_raw))
  trimws(strsplit(skip_raw, ",", fixed = TRUE)[[1]]) else character(0)

cat("============================================================\n")
cat(" Exploratory plots pipeline\n")
cat(" Run id : ", RUN_ID, "\n")
cat(" Outputs: output/application/plots/\n")
if (length(SKIP)) cat(" Skipping: ", paste(SKIP, collapse = ", "), "\n")
cat("============================================================\n\n")

step <- function(tag, label, expr) {
  if (tag %in% SKIP) {
    cat(">>> ", label, " ... [skipped via BUILD_EXPL_SKIP]\n\n", sep = "")
    return(invisible(TRUE))
  }
  cat(">>> ", label, " ...\n", sep = "")
  t0 <- Sys.time()
  ok <- tryCatch({ force(expr); TRUE },
                 error = function(e) {
                   cat("    FAILED: ", conditionMessage(e), "\n", sep = "")
                   FALSE
                 })
  cat(sprintf("    [%s, %.1fs]\n\n",
              if (ok) "ok" else "FAIL",
              as.numeric(Sys.time() - t0, units = "secs")))
  invisible(ok)
}

# Stage 1: per-(dataset, model, K) plots (adjacency / success / similarity /
# rho_d / rank-vs-rank) using the canonical z_hat
step("perfit",
     "1/4 per-fit exploratory plots (plotting_script.R)",
     source("scripts/analysis/plotting_script.R", local = TRUE))

# Stage 2: posterior K uncertainty (PMF, trace, panel)
step("posteriorK",
     "2/4 posterior K uncertainty plots (posterior_K_uq.R)",
     source("scripts/analysis/posterior_K_uq.R", local = TRUE))

# Stage 3: OSBM network/heatmap visualisations
step("osbm",
     "3/4 OSBM visualisations (osbm_visualization.R)",
     {
       source("scripts/analysis/osbm_visualization.R", local = TRUE)
       run_all_osbm_visualizations()
     })

# Stage 4: experimental visualisations (arc / heatmap / caterpillar / ring)
step("newviz",
     "4/4 experimental visualisations (new_visualizations.R)",
     {
       source("scripts/analysis/new_visualizations.R", local = TRUE)
       run_all_new_visualizations()
     })

cat("Done. Exploratory plots in: output/application/plots/\n")
