#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# build_violation_zeta_emp_table.R
#
# Builds a paper-facing LaTeX table reporting, for every dataset x model:
#   - K_hat                            (occupied blocks at the canonical z_hat)
#   - violation_rate_zhat              (single canonical partition)
#   - thetaW_block_emp  mean [lo, hi]  (posterior empirical block WST conformity)
#   - thetaS_block_emp  mean [lo, hi]  (posterior empirical block SST conformity)
#   - thetaW_item_emp   mean [lo, hi]  (posterior empirical item WST conformity)
#   - thetaS_item_emp   mean [lo, hi]  (posterior empirical item SST conformity)
# WST/SST use the intrinsic block order. DC-SBM is strength-reordered. All
# posterior summaries use only draws with K equal to the displayed K_hat.
#
# Source: `output/application/tables/hierarchy_diagnostics_overview.csv`,
# produced by `scripts/analysis/audit_hierarchy_synopsis.R`.
#
# Writes to:
#   output/paper/tables/<run-id>/violation_zeta_emp_full_paper.tex
#   output/paper/tables/<run-id>/violation_zeta_emp_full_paper.csv
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

OVERVIEW <- Sys.getenv(
  "HIERARCHY_OVERVIEW_CSV",
  unset = "output/application/tables/hierarchy_diagnostics_overview.csv"
)
source("scripts/bundle_defaults.R", local = TRUE)
APP_RUN_DIR <- bundle_resolve_application_run_dir(must_exist = TRUE)
run_id  <- basename(APP_RUN_DIR)
OUT_DIR <- Sys.getenv(
  "APP_PAPER_TABLES_DIR",
  unset = file.path("output/paper/tables", run_id)
)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(OVERVIEW))
df <- readr::read_csv(OVERVIEW, show_col_types = FALSE)

# Friendly dataset labels for the paper
ds_label <- c(
  moreno_sheep   = "Bighorn sheep",
  strauss_2019b  = "Spotted hyenas",
  mountain_goats = "Mountain goats",
  macaques_data  = "Macaques",
  citations_data = "Stat.\\ journals",
  high_school    = "High school"
)
ds_order <- c("moreno_sheep","strauss_2019b","mountain_goats",
              "citations_data","macaques_data","high_school")
m_order  <- c("WST","SST","DCSBM")

df <- df |>
  mutate(
    dataset   = factor(dataset,   levels = ds_order),
    fit_model = factor(fit_model, levels = m_order)
  ) |>
  arrange(dataset, fit_model)

fmt3 <- function(x) ifelse(is.na(x), "---", formatC(x, digits = 3, format = "f"))
fmt_int <- function(x) ifelse(is.na(x), "---", formatC(round(x), format = "d", big.mark = "{,}"))
fmt_pct <- function(x) ifelse(is.na(x), "---", paste0(formatC(100*x, digits = 1, format = "f"), "\\%"))
fmt_pct_count <- function(rate, back, cross) {
  ifelse(is.na(rate), "---", paste0(fmt_pct(rate), "\\,(", fmt_int(back), "/", fmt_int(cross), ")"))
}
fmt_prob_pct <- function(x) {
  ifelse(
    is.na(x), "---",
    ifelse(x > 0 & x < 0.001, "<0.1\\%", paste0(formatC(100*x, digits = 1, format = "f"), "\\%"))
  )
}
fmt_ci  <- function(m, lo, hi) {
  ifelse(is.na(m), "---",
         sprintf("%s\\,[%s,\\,%s]", fmt3(m), fmt3(lo), fmt3(hi)))
}

tbl <- df |>
  transmute(
    Dataset = ds_label[as.character(dataset)],
    Model   = ifelse(as.character(fit_model) == "DCSBM", "DC--SBM",
                     as.character(fit_model)),
    K       = K_hat,
    Viol    = fmt_pct_count(violation_rate_zhat, violation_count_zhat, cross_count_zhat),
    zWb     = fmt_ci(thetaW_block_emp_mean, thetaW_block_emp_lo, thetaW_block_emp_hi),
    zSb     = fmt_ci(thetaS_block_emp_mean, thetaS_block_emp_lo, thetaS_block_emp_hi),
    zWi     = fmt_ci(thetaW_item_emp_mean,  thetaW_item_emp_lo,  thetaW_item_emp_hi),
    zSi     = fmt_ci(thetaS_item_emp_mean,  thetaS_item_emp_lo,  thetaS_item_emp_hi)
  )

readr::write_csv(tbl, file.path(OUT_DIR, "violation_zeta_emp_full_paper.csv"))

dc_wst_prob <- df |>
  filter(as.character(fit_model) == "DCSBM") |>
  transmute(dataset_key = as.character(dataset), dcsbm_wst_post = p_post_wst)

compact_tbl <- df |>
  mutate(dataset_key = as.character(dataset)) |>
  left_join(dc_wst_prob, by = "dataset_key") |>
  transmute(
    Dataset = ds_label[as.character(dataset)],
    Model   = ifelse(as.character(fit_model) == "DCSBM", "DC--SBM",
                     as.character(fit_model)),
    K       = K_hat,
    Viol    = fmt_pct_count(violation_rate_zhat, violation_count_zhat, cross_count_zhat),
    DCSBMWSTPost = ifelse(as.character(fit_model) == "DCSBM", fmt_prob_pct(dcsbm_wst_post), "")
  )
readr::write_csv(compact_tbl, file.path(OUT_DIR, "violation_wst_region_paper.csv"))

# Build LaTeX body. We emit only the tabular core; the wrapping environment
# lives in main.tex (\begin{table}\input{...}\end{table}) so the caller can
# control captions/labels there.
body_rows <- character(0)
prev_ds <- NA_character_
for (i in seq_len(nrow(tbl))) {
  row <- tbl[i, ]
  ds_cell <- if (!identical(row$Dataset, prev_ds))
    sprintf("\\multirow{3}{*}{%s}", row$Dataset) else ""
  body_rows <- c(body_rows, sprintf(
    "%s & %s & %d & %s & $%s$ & $%s$ & $%s$ & $%s$ \\\\",
    ds_cell, row$Model, row$K, row$Viol,
    row$zWb, row$zSb, row$zWi, row$zSi
  ))
  if (i %% 3L == 0L && i < nrow(tbl)) {
    body_rows <- c(body_rows, "\\midrule")
  }
  prev_ds <- row$Dataset
}

# Output two flavours:
#  (1) a self-contained `table` environment (drop-in via \input{...}) and
#  (2) a tabular-only body, in case the user prefers to wrap manually.

tabular_only <- c(
  "% Auto-generated by scripts/analysis/build_violation_zeta_emp_table.R",
  "% Do NOT edit by hand.",
  "\\begin{tabular}{llrlcccc}",
  "\\toprule",
  paste(
    "Dataset & Model & $\\hat K$ & $\\zeta^{\\mathrm{viol}}_{\\hat z}$ (back/cross) &",
    "$\\bar{\\zeta}^{\\mathrm{blk}}_{W}$ &",
    "$\\bar{\\zeta}^{\\mathrm{blk}}_{S}$ &",
    "$\\bar{\\zeta}^{\\mathrm{item}}_{W}$ &",
    "$\\bar{\\zeta}^{\\mathrm{item}}_{S}$ \\\\"
  ),
  "\\midrule",
  body_rows,
  "\\bottomrule",
  "\\end{tabular}"
)
writeLines(tabular_only, file.path(OUT_DIR, "violation_zeta_emp_full_paper_tabular.tex"))

full_env <- c(
  "% Auto-generated by scripts/analysis/build_violation_zeta_emp_table.R",
  "% Do NOT edit by hand.",
  "\\begin{table}[t]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
    paste0("\\caption{Canonical violation rate and posterior empirical ",
      "conformity scores for every dataset $\\times$ model fit. ",
      "Violations are computed at the canonical partition ",
      "$\\hat z = \\mathrm{get\\_z\\_hat\\_from\\_draws}$ ",
      "(single posterior draw via $\\mathrm{minVI}$). WST and Toeplitz ",
      "SST keep the ordered-prior labels; DC--SBM draws are relabelled ",
      "by empirical block strength. $\\bar{\\zeta}_{W,S}^{\\mathrm{blk}}$ ",
      "and $\\bar{\\zeta}_{W,S}^{\\mathrm{item}}$ are empirical block-/item-level ",
      "conformity rates, reported as posterior means with 95\\% credible ",
      "intervals after filtering to draws with $K=\\hat K$. ",
      "The violation column also reports the corresponding backward ",
      "cross-block mass over total cross-block mass. ",
         "Values are produced by ",
         "\\texttt{summarise\\_osbm\\_diagnostics} (WST/SST) and ",
         "\\texttt{summarise\\_dcsbm\\_diagnostics} (DC-SBM) via ",
         "\\texttt{scripts/analysis/audit\\_hierarchy\\_synopsis.R}.}"),
  "\\label{tab:violation-zeta-emp-full}",
  "\\input{../output/paper/tables/current/violation_zeta_emp_full_paper_tabular.tex}",
  "\\end{table}"
)
writeLines(full_env, file.path(OUT_DIR, "violation_zeta_emp_full_paper.tex"))

compact_body <- character(0)
prev_ds <- NA_character_
for (i in seq_len(nrow(compact_tbl))) {
  row <- compact_tbl[i, ]
  is_new_ds <- !identical(row$Dataset, prev_ds)
  ds_cell <- if (is_new_ds) sprintf("\\multirow{3}{*}{%s}", row$Dataset) else ""
  post_cell <- if (nzchar(row$DCSBMWSTPost)) sprintf("$%s$", row$DCSBMWSTPost) else ""
  compact_body <- c(compact_body, sprintf(
    "%s & %s & %d & %s & %s \\\\",
    ds_cell, row$Model, row$K, row$Viol, post_cell
  ))
  if (i %% 3L == 0L && i < nrow(compact_tbl)) {
    compact_body <- c(compact_body, "\\addlinespace")
  }
  prev_ds <- row$Dataset
}

compact_env <- c(
  "% Auto-generated by scripts/analysis/build_violation_zeta_emp_table.R",
  "% Do NOT edit by hand.",
  "\\begin{table}[htpb]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{@{}llclc@{}}",
  "\\toprule",
  paste(
    "Dataset & Model & $K$ &",
    "$\\zeta^{\\mathrm{viol}}_{\\hat z}$ (back/cross) &",
    "\\makecell[c]{$\\Pr_{\\mathrm{DC}}\\{\\rho\\in\\mathcal C_{\\mathrm{WST}}\\mid A,N\\}$} \\\\"
  ),
  "\\midrule",
  compact_body,
  "\\bottomrule",
  "\\end{tabular}",
  paste0("\\caption{Canonical violation rate and exact WST-region posterior ",
         "probability. $K$ is the block count of the \\texttt{minVI} ",
         "partition. The violation column reports the weighted fraction of ",
         "cross-block edge mass running against the displayed block order, ",
         "with the absolute backward/cross-block counts in parentheses. ",
         "The last column is computed only from the unconstrained DC--SBM ",
         "posterior and is printed on the DC--SBM row: after relabelling each draw by block strength and ",
         "conditioning on the displayed $K$, it is the posterior probability ",
         "that all upper-triangular directional probabilities satisfy WST.}"),
  "\\label{tab:violation-wst-region}",
  "\\end{table}"
)
writeLines(compact_env, file.path(OUT_DIR, "violation_wst_region_paper.tex"))

cat("Wrote:\n",
    " ", file.path(OUT_DIR, "violation_zeta_emp_full_paper_tabular.tex"), "\n",
    " ", file.path(OUT_DIR, "violation_zeta_emp_full_paper.tex"), "\n",
    " ", file.path(OUT_DIR, "violation_zeta_emp_full_paper.csv"), "\n",
    " ", file.path(OUT_DIR, "violation_wst_region_paper.tex"), "\n",
    " ", file.path(OUT_DIR, "violation_wst_region_paper.csv"), "\n",
    sep = "")
