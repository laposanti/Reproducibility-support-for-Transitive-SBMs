#!/usr/bin/env Rscript
# =============================================================================
# build_paper_loo_table.R
#
# Builds the model-selection table that goes into the paper. The table reports,
# for every (dataset, model) pair: the posterior median number of blocks K with
# its 95% credible interval, the LOOIC, and -- only on the second-best row of
# each dataset -- the difference in expected log predictive density relative to
# the within-dataset winner together with its standard error and the |t_LOO|
# ratio delta_elpd / SE(delta_elpd). A dagger marks rows above the two-sided
# 95% t reference qt(0.975, n_obs - 1). No p-values or verdict column.
#
# Inputs (from a single run):
#   <run>/all_results.csv               (K posterior summaries, looic, elpd_loo)
#   <run>/model_comparison_loo.csv      (delta_elpd, SE, |t|, n_obs)
#
# Outputs:
#   <out_dir>/model_selection_paper.tex   LaTeX (booktabs, multirow)
#   <out_dir>/model_selection_paper.csv   tidy long form for record keeping
#
# Run:
#   APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#   APP_PAPER_TABLES_DIR=output/paper/tables/application_run_20260414_104327 \
#   Rscript scripts/analysis/build_paper_loo_table.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(glue)
})

suppressPackageStartupMessages({
  source("helper_folder/simulation/simulation_study_helpers.R", local = TRUE)
  source("scripts/analysis/osbm_visualization.R", local = TRUE)
})
source("scripts/bundle_defaults.R", local = TRUE)

# ---- 1. Resolve paths ------------------------------------------------------
run_arg <- ""
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) run_arg <- args[[1]]
RUN_DIR <- bundle_resolve_application_run_dir(run_dir = run_arg, must_exist = TRUE)

run_id <- basename(normalizePath(RUN_DIR))

OUT_DIR <- Sys.getenv("APP_PAPER_TABLES_DIR", unset = "")
if (!nzchar(OUT_DIR)) {
  OUT_DIR <- file.path("output", "paper", "tables", run_id)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Run id   :", run_id, "\n")
cat("Run dir  :", RUN_DIR, "\n")
cat("Out dir  :", OUT_DIR, "\n")

# ---- 2. Load inputs --------------------------------------------------------
all_res_path <- file.path(RUN_DIR, "all_results.csv")
loo_path     <- file.path(RUN_DIR, "model_comparison_loo.csv")
stopifnot(file.exists(all_res_path), file.exists(loo_path))

# all_results.csv has no header in some runs; sniff it
header_line <- readLines(all_res_path, n = 1L)
if (!grepl("dataset", header_line)) {
  res <- readr::read_csv(
    all_res_path, show_col_types = FALSE,
    col_names = c("dataset","fit_model","K_hat","K_mean","K_lo","K_hi",
                  "elpd_loo","p_loo","looic","pk_max","pk_bad",
                  "b_kappa","time_sec","note")
  )
} else {
  res <- readr::read_csv(all_res_path, show_col_types = FALSE)
}

cmp <- readr::read_csv(loo_path, show_col_types = FALSE)

count_unordered_dyads <- function(ds) {
  tryCatch({
    A <- choose_dataset_local(ds)
    if (inherits(A, "Matrix")) A <- as.matrix(A)
    DI <- build_dyad_index(nrow(A))
    if (is.list(DI) && !is.null(DI$i)) {
      length(DI$i)
    } else if (is.matrix(DI)) {
      nrow(DI)
    } else {
      length(DI)
    }
  }, error = function(e) NA_integer_)
}

n_obs_lookup <- vapply(as.character(unique(cmp$dataset)),
                       count_unordered_dyads,
                       integer(1))
cmp <- cmp %>%
  mutate(n_obs = dplyr::coalesce(unname(n_obs_lookup[dataset]),
                                 as.integer(n_obs)))

# ---- 3. Display labels & ordering -----------------------------------------
dataset_label <- c(
  moreno_sheep    = "Bighorn sheep",
  strauss_2019b   = "Spotted hyenas",
  mountain_goats  = "Mountain goats",
  citations_data  = "Stat.\\ journals",
  macaques_data   = "Macaques",
  high_school     = "High school"
)
dataset_order <- names(dataset_label)

model_label <- c(
  WST   = "WST--OSBM",
  SST   = "Toeplitz SST--OSBM",
  DCSBM = "DC--SBM"
)

# ---- 4. Within-dataset ranking by elpd_loo --------------------------------
res <- res %>%
  filter(dataset %in% dataset_order) %>%
  mutate(
    dataset    = factor(dataset, levels = dataset_order),
    model_disp = model_label[fit_model]
  ) %>%
  arrange(dataset, desc(elpd_loo)) %>%
  group_by(dataset) %>%
  mutate(
    rank          = row_number(),
    elpd_winner   = elpd_loo[rank == 1],
    looic_winner  = looic[rank == 1]
  ) %>%
  ungroup()

# attach delta_elpd / SE on the row that is "second_model" within each dataset
cmp_small <- cmp %>%
  transmute(
    dataset,
    fit_model      = second_model,
    delta_elpd_se  = se_delta_elpd,
    delta_elpd     = delta_elpd,
    z_ratio        = z_ratio,
    n_obs          = n_obs
  )

res <- res %>%
  left_join(cmp_small, by = c("dataset" = "dataset", "fit_model" = "fit_model"))

# For rows that are neither winner nor 2nd (the 3rd model), compute
# delta_elpd directly. SE for the 3rd model vs 1st is in cmp$se_delta_3rd.
cmp_third <- cmp %>%
  transmute(
    dataset,
    fit_model      = third_model,
    delta_elpd_3   = delta_elpd_3rd,
    se_delta_3     = se_delta_3rd,
    n_obs_3        = n_obs
  )
res <- res %>%
  left_join(cmp_third, by = c("dataset" = "dataset", "fit_model" = "fit_model")) %>%
  mutate(
    delta_elpd_eff = dplyr::case_when(
      rank == 1L ~ NA_real_,
      rank == 2L ~ delta_elpd,
      rank == 3L ~ delta_elpd_3,
      TRUE       ~ NA_real_
    ),
    se_eff = dplyr::case_when(
      rank == 1L ~ NA_real_,
      rank == 2L ~ delta_elpd_se,
      rank == 3L ~ se_delta_3,
      TRUE       ~ NA_real_
    ),
    n_obs_eff = dplyr::case_when(
      rank == 1L ~ NA_integer_,
      rank == 2L ~ as.integer(n_obs),
      rank == 3L ~ as.integer(n_obs_3),
      TRUE       ~ NA_integer_
    ),
    loo_df = ifelse(is.na(n_obs_eff), NA_integer_, n_obs_eff - 1L),
    z_eff = ifelse(is.na(se_eff) | se_eff == 0, NA_real_,
                   abs(delta_elpd_eff / se_eff)),
    t_crit_95 = ifelse(is.na(loo_df), NA_real_, stats::qt(0.975, df = loo_df)),
    clear_95 = !is.na(z_eff) & !is.na(t_crit_95) & z_eff >= t_crit_95
  )

# ---- 5. Format helpers -----------------------------------------------------
fmt_int  <- function(x) ifelse(is.na(x), "", formatC(x, format = "d",
                                                    big.mark = "{,}"))
fmt_num  <- function(x, d = 1) ifelse(is.na(x),
                                      "",
                                      formatC(round(x, d), format = "f",
                                              big.mark = "{,}", digits = d))
fmt_K    <- function(hat, lo, hi) {
  if (any(is.na(c(hat, lo, hi)))) {
    return("")
  }
  sprintf("$%s\\;[%s,%s]$",
          fmt_num(hat, 0),
          fmt_num(lo, 0),
          fmt_num(hi, 0))
}

write_model_selection_paper_outputs <- function(res, dataset_order, dataset_label,
                                                out_dir, run_id, run_dir) {
  tex_lines <- character()
  push <- function(x) tex_lines[[length(tex_lines) + 1L]] <<- x

  push("\\begin{tabular}{llcrrrr}")
  push("\\toprule")
  push(paste0(
    "Dataset & Model & $\\hat K\\;[95\\%\\text{ CrI}]$ & LOOIC & ",
    "$\\Delta\\text{ELPD}$ & SE($\\Delta\\text{ELPD}$) & $|t_{\\mathrm{LOO}}|$ \\\\"
  ))
  push("\\midrule")

  n_ds <- length(dataset_order)
  for (k in seq_along(dataset_order)) {
    ds <- dataset_order[[k]]
    rows <- res %>% filter(dataset == ds) %>% arrange(rank)
    if (!nrow(rows)) next
    ds_lab <- dataset_label[[ds]]

    for (i in seq_len(nrow(rows))) {
      r <- rows[i, ]
      ds_cell <- if (i == 1L) sprintf("\\multirow{%d}{*}{%s}", nrow(rows), ds_lab) else ""

      model_cell <- if (i == 1L) sprintf("\\textbf{%s}", r$model_disp) else r$model_disp
      K_cell     <- fmt_K(r$K_hat, r$K_lo, r$K_hi)
      looic_cell <- if (i == 1L) sprintf("$\\mathbf{%s}$", fmt_num(r$looic, 1))
                    else                  sprintf("$%s$",         fmt_num(r$looic, 1))

      if (i == 1L) {
        delta_cell <- "$0$"
        se_cell    <- "---"
        z_cell     <- "---"
      } else {
        delta_cell <- sprintf("$%s$", fmt_num(r$delta_elpd_eff, 1))
        se_cell    <- sprintf("$%s$", fmt_num(r$se_eff, 1))
        z_val      <- fmt_num(r$z_eff, 2)
        z_cell     <- if (isTRUE(r$clear_95)) {
          sprintf("$\\mathbf{%s}^{\\dagger}$", z_val)
        } else {
          sprintf("$%s$", z_val)
        }
      }

      push(sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
                   ds_cell, model_cell, K_cell, looic_cell,
                   delta_cell, se_cell, z_cell))
    }
    if (k < n_ds) push("\\addlinespace")
  }
  push("\\bottomrule")
  push("\\end{tabular}")

  tex_path <- file.path(out_dir, "model_selection_paper.tex")
  writeLines(tex_lines, tex_path)
  cat("Wrote:", tex_path, "\n")

  csv_out <- res %>%
    select(dataset, model = model_disp, rank,
           K_hat, K_lo, K_hi, elpd_loo, looic,
           delta_elpd = delta_elpd_eff,
           se_delta_elpd = se_eff,
           z_ratio = z_eff,
           t_ratio = z_eff,
           n_obs = n_obs_eff,
           loo_df,
           t_crit_95,
           clear_95)

  csv_path <- file.path(out_dir, "model_selection_paper.csv")
  readr::write_csv(csv_out, csv_path)
  cat("Wrote:", csv_path, "\n")

  manifest_path <- file.path(out_dir, "MANIFEST.txt")
  writeLines(c(
    paste0("run_id: ", run_id),
    paste0("run_dir: ", normalizePath(run_dir)),
    paste0("generated_by: scripts/analysis/build_paper_loo_table.R"),
    paste0("generated_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("R_version: ", R.version.string)
  ), manifest_path)
  cat("Wrote:", manifest_path, "\n")
}

# ---- 6. Write outputs ------------------------------------------------------
write_model_selection_paper_outputs(
  res = res,
  dataset_order = dataset_order,
  dataset_label = dataset_label,
  out_dir = OUT_DIR,
  run_id = run_id,
  run_dir = RUN_DIR
)
