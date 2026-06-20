#!/usr/bin/env Rscript
# ----------------------------------------------------------------------
# Audit: recompute the hierarchy-synopsis numbers for the paper from the
# current MCMC run.
#
# Partition extraction (canonical, see get_z_hat_from_draws):
#   * WST / SST: kept as-is from the MCMC chain; minVI with
#                method = "draws" returns a single posterior sample, so
#                labels carry the prior's intrinsic 1 = strongest ordering.
#   * DC-SBM:    every draw is per-draw relabelled by empirical block
#                score (compute_block_scores), then minVI is taken on the
#                relabelled chain.
#
# No ex-post strength reordering is applied: the resulting partition is
# already canonical for each model, and a second empirical reorder would
# re-shuffle WST/SST labels against the model's own order.
#
# Per (dataset, model) we compute:
#   K            -- block count of z_hat.
#   viol_zhat    -- official violation rate: sum_{z_i>z_j} A_ij over
#                   sum_{z_i!=z_j} A_ij on directed edges.
# Posterior conformity summaries are produced by summarise_*_diagnostics:
#   WST/SST keep their intrinsic labels; DC-SBM is reordered by empirical
#   strength; all posterior means use draws with K equal to the displayed K.
#
# Outputs:
#   output/paper/tables/<run>/hierarchy_synopsis_audit.csv
#   output/paper/tables/<run>/violation_rates_by_model.csv
#   output/application/tables/hierarchy_diagnostics_<ds>.csv  (per dataset)
#   output/application/tables/hierarchy_diagnostics_overview.csv
# ----------------------------------------------------------------------

suppressPackageStartupMessages({
  source("scripts/analysis/osbm_visualization.R", chdir = FALSE)
  source("helper_folder/diagnostics/transitivity_diagnostics.R", chdir = FALSE)
})
source("scripts/bundle_defaults.R", local = TRUE)

RUN_DIR <- bundle_resolve_application_run_dir(must_exist = TRUE)
OUT_DIR <- file.path("output/paper/tables", basename(RUN_DIR))
APP_TBL_DIR <- "output/application/tables"
dir.create(OUT_DIR,     recursive = TRUE, showWarnings = FALSE)
dir.create(APP_TBL_DIR, recursive = TRUE, showWarnings = FALSE)

DATASETS <- c("moreno_sheep", "strauss_2019b", "mountain_goats",
              "citations_data", "macaques_data", "high_school")
MODELS   <- c("WST", "SST", "DCSBM")

block_count_matrix <- function(A, z) {
  K <- max(z, na.rm = TRUE)
  C <- matrix(0, K, K)
  for (k in seq_len(K)) {
    rows <- which(z == k)
    for (l in seq_len(K)) {
      cols <- which(z == l)
      if (length(rows) && length(cols)) C[k, l] <- sum(A[rows, cols])
    }
  }
  C
}

viol_rate_at <- function(A, z) {
  ij <- which(A > 0, arr.ind = TRUE)
  if (!nrow(ij)) return(c(rate = NA_real_, back = 0, cross = 0))
  zi <- z[ij[, 1]]; zj <- z[ij[, 2]]; w <- A[ij]
  cross <- zi != zj
  back  <- cross & (zi > zj)
  total_cross <- sum(w[cross])
  total_back  <- sum(w[back])
  c(rate = if (total_cross > 0) total_back / total_cross else NA_real_,
    back = total_back, cross = total_cross)
}

rows <- list()
rich_rows <- list()
for (ds in DATASETS) {
  ppc_path <- file.path("output/application/ppc", ds, paste0(ds, "_A_obs.rds"))
  A <- if (file.exists(ppc_path)) readRDS(ppc_path) else choose_dataset_local(ds)
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  n_nodes <- nrow(A)
  for (m in MODELS) {
    fp <- file.path(RUN_DIR, paste0(ds, "_", m, "_fit.rds"))
    if (!file.exists(fp)) next
    fit <- readRDS(fp)
    z_obj <- get_z_hat_from_draws(fit, A, model = m)
    z_hat <- z_obj$z_hat
    K_hat <- max(z_hat, na.rm = TRUE)

    v  <- viol_rate_at(A, z_hat)

    rows[[length(rows) + 1L]] <- data.frame(
      dataset         = ds,
      model           = m,
      n               = n_nodes,
      K               = K_hat,
      viol_zhat       = unname(v["rate"]),
      back_count      = unname(v["back"]),
      cross_count     = unname(v["cross"])
    )
    cat(sprintf(
      "%-15s %-5s K=%2d  viol=%.3f\n",
      ds, m, K_hat, v["rate"]))

    # --- Rich per-fit diagnostics (theta_W/S block, EBF, transitive triads,
    #     hierarchy/curl energy, etc.). The fits on disk are already
    #     relabelled by application.R (ECR for DCSBM; WST/SST untouched),
    #     so summarise_*_diagnostics is invoked on `fit` as-is, with the
    #     canonical z_hat (no ex-post strength reorder) for the *_zhat
    #     fields. See helper_folder/diagnostics/transitivity_diagnostics.R.
    rich_vec <- tryCatch({
      if (m %in% c("WST", "SST")) {
        summarise_osbm_diagnostics(
          out_relab = fit, regime = m, K_max_hint = K_hat,
          z_hat = z_hat, n = n_nodes, m_items = 2000L, alpha = 0.5,
          A = A, T_block = NULL, seed_block = 123L,
          method_order = "identity"
        )
      } else {
        summarise_dcsbm_diagnostics(
          fit = fit, z_hat = z_hat, K = K_hat, n = n_nodes,
          m_items = 2000L, alpha = 0.5, A = A,
          T_block = NULL, seed_block = 123L,
          method_order = "mean"
        )
      }
    }, error = function(e) {
      cat(sprintf("  [WARN] %s/%s rich-diagnostics failed: %s\n",
                  ds, m, conditionMessage(e)))
      .empty_diag_vec()
    })

    rich_rows[[length(rich_rows) + 1L]] <- cbind(
      data.frame(dataset = ds, fit_model = m,
                 n = n_nodes, K_hat = K_hat,
                 stringsAsFactors = FALSE),
      as.data.frame(as.list(rich_vec), stringsAsFactors = FALSE)
    )
  }
}
df <- do.call(rbind, rows)

readr::write_csv(df, file.path(OUT_DIR, "hierarchy_synopsis_audit.csv"))
cat("\nWrote:", file.path(OUT_DIR, "hierarchy_synopsis_audit.csv"), "\n")

slim <- df[, c("dataset", "model", "K", "viol_zhat",
               "back_count", "cross_count")]
names(slim) <- c("dataset", "model", "K", "viol_rate", "viol_count",
                 "cross_count")
readr::write_csv(slim, file.path(OUT_DIR, "violation_rates_by_model.csv"))
cat("Wrote:", file.path(OUT_DIR, "violation_rates_by_model.csv"), "\n")

diag_df <- df
names(diag_df)[names(diag_df) == "model"]        <- "fit_model"
names(diag_df)[names(diag_df) == "K"]            <- "K_hat"
names(diag_df)[names(diag_df) == "viol_zhat"]    <- "violation_rate_zhat"
names(diag_df)[names(diag_df) == "back_count"]   <- "violation_count_zhat"
names(diag_df)[names(diag_df) == "cross_count"]  <- "cross_count_zhat"
diag_df <- diag_df[, c("dataset", "fit_model", "n", "K_hat",
                       "violation_rate_zhat", "violation_count_zhat",
                       "cross_count_zhat")]

# Merge in rich per-fit diagnostics (theta_W/S, EBF, hierarchy energy, ...)
rich_df <- do.call(rbind, rich_rows)
# Drop any rich columns that duplicate the canonical empirical ones
# we just renamed (violation_rate_zhat / violation_count_zhat come from
# the helper in the rich frame as well; we keep the canonical-z_hat
# versions computed above).
drop_dupes <- intersect(
  names(rich_df),
  c("violation_rate_zhat", "violation_count_zhat")
)
if (length(drop_dupes)) rich_df <- rich_df[, !names(rich_df) %in% drop_dupes]

diag_df <- merge(diag_df, rich_df,
                 by = c("dataset", "fit_model", "n", "K_hat"),
                 all.x = TRUE, sort = FALSE)

readr::write_csv(diag_df,
  file.path(APP_TBL_DIR, "hierarchy_diagnostics_overview.csv"))
cat("Wrote:", file.path(APP_TBL_DIR, "hierarchy_diagnostics_overview.csv"), "\n")
for (ds in unique(diag_df$dataset)) {
  one <- diag_df[diag_df$dataset == ds, , drop = FALSE]
  readr::write_csv(one,
    file.path(APP_TBL_DIR, paste0("hierarchy_diagnostics_", ds, ".csv")))
}
cat("Wrote per-dataset hierarchy_diagnostics_*.csv in", APP_TBL_DIR, "\n")
