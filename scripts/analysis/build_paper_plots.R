#!/usr/bin/env Rscript
# =============================================================================
# build_paper_plots.R -- ONE entry point that produces every figure file the
# manuscript reads via \includegraphics{...}.
#
# Manuscript figures (tex file/main.tex), grouped by source script:
#
#   regen_paper_network_figs.R
#     tex file/Figures/<dataset>_combined_block_networks_clean.pdf  (x6)
#     tex file/Figures/moreno_sheep_{DCSBM,SST}_network.png
#     (also output/paper/figures/<run>/ copies)
#
#   build_bt_delta_summary.R
#     output/paper/figures/<run>/bt_delta_wst_applications.{pdf,png}
#
#   scripts/diagnostics/simulate_support_geometries_cube.R
#     output/diagnostics/support_geometry/support_3d_shaded_geometry.png
#     (already produced; copied into the paper figures folder)
#
# This wrapper produces ONLY the small set of figures the paper actually
# \includegraphics{}-es, so that they always reflect the audit-canonical
# violation counts (strength-reordered z_hat;
# helper_folder/transitivity_check_helper.R::violation_rate_zhat).
#
# To also regenerate the global exploratory plot library
# (output/application/plots/<dataset>/...), use
#   Rscript scripts/analysis/build_exploratory_plots.R
# or, to refresh everything in one shot:
#   Rscript scripts/analysis/build_all_plots.R
#
# Usage:
#   APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#     Rscript scripts/analysis/build_paper_plots.R
# =============================================================================

RUN_DIR <- Sys.getenv("APP_RUN_DIR", unset = "")
if (!nzchar(RUN_DIR)) {
  cur <- "output/paper/figures/current"
  if (file.exists(cur))
    RUN_DIR <- file.path("output/application/raw", basename(Sys.readlink(cur)))
}
if (!nzchar(RUN_DIR) || !dir.exists(RUN_DIR))
  stop("Set APP_RUN_DIR=output/application/raw/<run_id>; got '", RUN_DIR, "'.")

RUN_ID  <- basename(RUN_DIR)
FIG_DIR <- Sys.getenv("APP_PAPER_FIGURES_DIR", unset = "")
if (!nzchar(FIG_DIR)) {
  FIG_DIR <- file.path("output/paper/figures", RUN_ID)
}
TAB_DIR <- Sys.getenv("APP_PAPER_TABLES_DIR", unset = "")
if (!nzchar(TAB_DIR)) {
  TAB_DIR <- file.path("output/paper/tables", RUN_ID)
}
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(APP_RUN_DIR = RUN_DIR,
           APP_PAPER_FIGURES_DIR = FIG_DIR,
           APP_PAPER_TABLES_DIR  = TAB_DIR)

cat("============================================================\n")
cat(" Paper plots pipeline\n")
cat(" Run id : ", RUN_ID, "\n")
cat(" Figures: ", FIG_DIR, "\n")
cat("============================================================\n\n")

step <- function(label, expr) {
  cat(">>> ", label, " ...\n", sep = "")
  t0 <- Sys.time()
  ok <- tryCatch({ force(expr); TRUE },
                 error = function(e) { cat("    FAILED: ", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("    [%s, %.1fs]\n\n", if (ok) "ok" else "FAIL",
              as.numeric(Sys.time() - t0, units = "secs")))
  invisible(ok)
}

step("1/3 combined block networks + sheep ordered networks (regen_paper_network_figs.R)",
     source("scripts/analysis/regen_paper_network_figs.R", local = TRUE))

step("2/3 bt_delta_wst_applications.{pdf,png} (build_bt_delta_summary.R)",
     source("scripts/analysis/build_bt_delta_summary.R", local = TRUE))

step("3/3 copy support_3d_shaded_geometry.png into paper figures",
     {
       src <- "output/diagnostics/support_geometry/support_3d_shaded_geometry.png"
       if (file.exists(src))
         file.copy(src, file.path(FIG_DIR, basename(src)), overwrite = TRUE)
       else
         message("    skipped: ", src, " not found; run scripts/diagnostics/simulate_support_geometries_cube.R first.")
     })

cat("Done. Paper figures in: ", FIG_DIR, "\n", sep = "")
