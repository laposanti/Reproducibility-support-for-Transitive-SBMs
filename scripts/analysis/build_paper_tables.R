#!/usr/bin/env Rscript
# =============================================================================
# build_paper_tables.R -- ONE entry point that produces every numeric artifact
# the manuscript reads via \input{...}.
#
# Manuscript inputs (tex file/main.tex):
#   ../output/paper/tables/current/model_selection_paper.tex
#   ../output/paper/tables/current/hierarchy_synopsis_paper.tex
#   ../output/paper/tables/current/violation_wst_region_paper.tex
#   ../output/paper/tables/current/violation_zeta_emp_full_paper.tex
#   ../output/paper/tables/current/bt_delta_summary_table.tex
#
# All three are written, together with the supporting CSV audits, under
#   output/paper/tables/<run_id>/
# and the `current/` symlink is repointed at <run_id>.
#
# Pipeline (canonical definitions in helper_folder/diagnostics/transitivity_diagnostics.R):
#   1. build_paper_loo_table.R     -> model_selection_paper.{tex,csv}
#   2. audit_hierarchy_synopsis.R  -> hierarchy_synopsis_audit.csv
#                                     violation_rates_by_model.csv
#   3. build_hierarchy_synopsis_table.R -> hierarchy_synopsis_paper.tex
#   4. build_violation_zeta_emp_table.R -> compact violation/WST-region table
#                                      and posterior empirical conformity table
#   5. build_bt_delta_summary.R    -> bt_delta_summary.{tex,csv}
#                                     bt_delta_wst_applications.{pdf,png}
#   6. build_application_supplement_tables.R
#                                  -> PSIS, cycle diagnostics, and high-school
#                                     cycle-inspection table for supplement
#
# Usage:
#   APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#     Rscript scripts/analysis/build_paper_tables.R
#
#   # Update output/paper/tables/current symlink to point at this run:
#   APP_BLESS_RUN=1 APP_RUN_DIR=... Rscript scripts/analysis/build_paper_tables.R
# =============================================================================

source("scripts/bundle_defaults.R", local = TRUE)
RUN_DIR <- bundle_resolve_application_run_dir(must_exist = TRUE)
RUN_ID  <- basename(RUN_DIR)
TAB_DIR <- file.path("output/paper/tables", RUN_ID)
FIG_DIR <- file.path("output/paper/figures", RUN_ID)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(APP_RUN_DIR = RUN_DIR,
           APP_PAPER_TABLES_DIR  = TAB_DIR,
           APP_PAPER_FIGURES_DIR = FIG_DIR)

cat("============================================================\n")
cat(" Paper tables pipeline\n")
cat(" Run id : ", RUN_ID, "\n")
cat(" Tables : ", TAB_DIR, "\n")
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

step("1/6 model_selection_paper.{tex,csv}  (build_paper_loo_table.R)",
     source("scripts/analysis/build_paper_loo_table.R", local = TRUE))

step("2/6 hierarchy_synopsis_audit.csv     (audit_hierarchy_synopsis.R)",
     source("scripts/analysis/audit_hierarchy_synopsis.R", local = TRUE))

step("3/6 hierarchy_synopsis_paper.tex     (build_hierarchy_synopsis_table.R)",
     source("scripts/analysis/build_hierarchy_synopsis_table.R", local = TRUE))

step("4/6 violation/WST-region tables    (build_violation_zeta_emp_table.R)",
     source("scripts/analysis/build_violation_zeta_emp_table.R", local = TRUE))

step("5/6 bt_delta_summary.{tex,csv,pdf}   (build_bt_delta_summary.R)",
     source("scripts/analysis/build_bt_delta_summary.R", local = TRUE))

step("6/6 application supplement tables    (build_application_supplement_tables.R)",
     source("scripts/analysis/build_application_supplement_tables.R", local = TRUE))

# Optional: repoint output/paper/{tables,figures}/current at this run
if (identical(tolower(Sys.getenv("APP_BLESS_RUN", unset = "")), "1") ||
    identical(tolower(Sys.getenv("APP_BLESS_RUN", unset = "")), "true")) {
  for (base in c("output/paper/tables", "output/paper/figures")) {
    cur <- file.path(base, "current")
    if (file.exists(cur) || file.exists(Sys.readlink(cur))) suppressWarnings(file.remove(cur))
    file.symlink(RUN_ID, cur)
    cat("Blessed: ", cur, " -> ", RUN_ID, "\n", sep = "")
  }
}

cat("Done. Paper tables in: ", TAB_DIR, "\n", sep = "")
