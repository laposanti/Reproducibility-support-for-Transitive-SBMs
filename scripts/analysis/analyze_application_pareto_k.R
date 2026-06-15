#!/usr/bin/env Rscript
# Diagnose application PSIS-LOO Pareto-k tails at the dyad level.
#
# Outputs:
#   output/paper/tables/<run_id>/application_pareto_k_stress_summary.csv
#   output/paper/tables/<run_id>/application_pareto_k_stress_summary.tex
#   output/paper/tables/<run_id>/application_pareto_k_top_dyads.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(loo)
})

RUN_DIR <- Sys.getenv("APP_RUN_DIR", unset = "")
if (!nzchar(RUN_DIR)) {
  cur <- "output/paper/tables/current"
  if (file.exists(cur)) {
    RUN_DIR <- file.path("output/application/raw", basename(Sys.readlink(cur)))
  }
}
if (!nzchar(RUN_DIR) || !dir.exists(RUN_DIR)) {
  stop("Set APP_RUN_DIR=output/application/raw/<run_id>; got '", RUN_DIR, "'.",
       call. = FALSE)
}

RUN_ID <- basename(RUN_DIR)
OUT_DIR <- Sys.getenv(
  "APP_PAPER_TABLES_DIR",
  unset = file.path("output", "paper", "tables", RUN_ID)
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  source("helper_folder/helper.R", local = TRUE)
  source("helper_folder/SST_helpers.R", local = TRUE)
  source("helper_folder/WST_helpers.R", local = TRUE)
  source("helper_folder/Hyper_setup.R", local = TRUE)
  source("helper_folder/sim_study_helper.R", local = TRUE)
  source("helper_folder/transitivity_check_helper.R", local = TRUE)
  source("core/DCSBM_varK.R", local = TRUE)
  source("scripts/analysis/osbm_visualization.R", local = TRUE)
  source("scripts/analysis/post_processing_helpers.R", local = TRUE)
})

DATASET_LABEL <- c(
  moreno_sheep   = "Bighorn sheep",
  strauss_2019b  = "Spotted hyenas",
  mountain_goats = "Mountain goats",
  citations_data = "Stat.\\ journals",
  macaques_data  = "Macaques",
  high_school    = "High school"
)
MODEL_LABEL <- c(
  WST = "WST--OSBM",
  SST = "Toeplitz SST--OSBM",
  DCSBM = "DC--SBM"
)
DATASET_ORDER <- names(DATASET_LABEL)
MODEL_ORDER <- names(MODEL_LABEL)

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "---",
         formatC(x, format = "f", digits = digits, big.mark = "{,}"))
}
fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "---",
         paste0(formatC(100 * x, format = "f", digits = digits), "\\%"))
}

read_all_results <- function() {
  path <- file.path(RUN_DIR, "all_results.csv")
  if (!file.exists(path)) stop("Missing all_results.csv in ", RUN_DIR, call. = FALSE)
  first <- readLines(path, n = 1L)
  if (grepl("dataset", first)) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    readr::read_csv(
      path,
      col_names = c("dataset", "fit_model", "K_hat", "K_mean", "K_lo", "K_hi",
                    "elpd_loo", "p_loo", "looic", "pk_max", "pk_bad",
                    "b_kappa", "time_sec", "note"),
      show_col_types = FALSE
    )
  }
}

load_cube_or_null <- function(ds, model) {
  tryCatch(load_fit_cube(ds, model, read_cube_meta(RUN_DIR)), error = function(e) NULL)
}

dyad_frame <- function(A, dyad_index, z_hat = NULL) {
  i <- dyad_index[, 1]
  j <- dyad_index[, 2]
  out <- data.frame(
    i = i,
    j = j,
    Aij = A[cbind(i, j)],
    Aji = A[cbind(j, i)]
  )
  out$Nij <- out$Aij + out$Aji
  out$zero_dyad <- out$Nij == 0
  out$reciprocal_positive <- out$Aij > 0 & out$Aji > 0
  if (!is.null(z_hat)) {
    out$zi <- z_hat[i]
    out$zj <- z_hat[j]
    out$same_block <- out$zi == out$zj
    out$block_pair <- paste(pmin(out$zi, out$zj), pmax(out$zi, out$zj), sep = "-")
  }
  out
}

summarise_fit <- function(ds, model, all_res) {
  fit_path <- file.path(RUN_DIR, paste0(ds, "_", model, "_fit.rds"))
  if (!file.exists(fit_path)) return(NULL)

  cat(sprintf("[%s | %s]\n", ds, model))
  A <- choose_dataset_local(ds)
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  DI <- build_dyad_index(nrow(A))
  fit <- readRDS(fit_path)
  LL <- if (model %in% c("WST", "SST")) {
    loglik_matrix_modular(A, fit, regime = model, dyad_index = DI)
  } else {
    loglik_matrix_dcsbm(A = as.matrix(A), dcsbm_out = fit, dyad_index = DI)
  }
  loo_fit <- loo::loo(LL)
  k <- as.numeric(loo_fit$diagnostics$pareto_k)

  cube <- load_cube_or_null(ds, model)
  z_hat <- if (!is.null(cube)) cube$z_hat else NULL
  dyads <- dyad_frame(A, DI, z_hat)
  dyads$pareto_k <- k
  dyads$dataset <- ds
  dyads$model <- model

  bad <- k > 0.7
  zero <- dyads$zero_dyad
  positive <- !zero
  reciprocal <- dyads$reciprocal_positive
  res_row <- all_res |>
    filter(dataset == ds, fit_model == model) |>
    slice_head(n = 1)
  pk_bad_run <- if (nrow(res_row) && !is.na(res_row$pk_bad)) {
    res_row$pk_bad
  } else {
    mean(bad, na.rm = TRUE)
  }

  list(
    summary = data.frame(
      dataset = ds,
      model = model,
      n_dyads = length(k),
      n_positive_dyads = sum(positive),
      p_loo = if (nrow(res_row)) res_row$p_loo else loo_fit$estimates["p_loo", "Estimate"],
      pareto_k_max = max(k, na.rm = TRUE),
      pareto_k_gt_0_7 = pk_bad_run,
      pareto_k_gt_0_7_recomputed = mean(bad, na.rm = TRUE),
      pareto_k_gt_1 = mean(k > 1, na.rm = TRUE),
      bad_zero_share = if (sum(bad) > 0) mean(zero[bad]) else NA_real_,
      bad_positive_share = if (sum(bad) > 0) mean(positive[bad]) else NA_real_,
      zero_bad_rate = if (sum(zero) > 0) mean(bad[zero]) else NA_real_,
      positive_bad_rate = if (sum(positive) > 0) mean(bad[positive]) else NA_real_,
      reciprocal_bad_rate = if (sum(reciprocal) > 0) mean(bad[reciprocal]) else NA_real_,
      nonreciprocal_positive_bad_rate =
        if (sum(positive & !reciprocal) > 0) mean(bad[positive & !reciprocal]) else NA_real_,
      median_edges_bad = if (sum(bad) > 0) median(dyads$Nij[bad]) else NA_real_,
      p90_edges_bad = if (sum(bad) > 0) as.numeric(quantile(dyads$Nij[bad], 0.9)) else NA_real_,
      max_edges_bad = if (sum(bad) > 0) max(dyads$Nij[bad]) else NA_real_
    ),
    top = dyads |>
      arrange(desc(pareto_k)) |>
      slice_head(n = 10)
  )
}

all_res <- read_all_results()
pieces <- list()
for (ds in DATASET_ORDER) {
  for (model in MODEL_ORDER) {
    pieces[[length(pieces) + 1L]] <- summarise_fit(ds, model, all_res)
  }
}
pieces <- Filter(Negate(is.null), pieces)

summary_tbl <- bind_rows(lapply(pieces, `[[`, "summary")) |>
  mutate(
    dataset = factor(dataset, levels = DATASET_ORDER),
    model = factor(model, levels = MODEL_ORDER),
    p_loo_per_dyad = p_loo / n_dyads,
    p_loo_per_positive_dyad = p_loo / n_positive_dyads
  ) |>
  arrange(dataset, model)

top_tbl <- bind_rows(lapply(pieces, `[[`, "top")) |>
  mutate(
    dataset = factor(dataset, levels = DATASET_ORDER),
    model = factor(model, levels = MODEL_ORDER)
  ) |>
  arrange(dataset, model, desc(pareto_k))

readr::write_csv(summary_tbl, file.path(OUT_DIR, "application_pareto_k_stress_summary.csv"))
readr::write_csv(top_tbl, file.path(OUT_DIR, "application_pareto_k_top_dyads.csv"))

highlight <- summary_tbl |>
  filter(pareto_k_gt_0_7 >= 0.10 | pareto_k_max >= 5) |>
  arrange(desc(pareto_k_gt_0_7))

rows <- character(0)
for (i in seq_len(nrow(highlight))) {
  r <- highlight[i, ]
  rows <- c(rows, sprintf(
    "%s & %s & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\\\",
    DATASET_LABEL[[as.character(r$dataset)]],
    MODEL_LABEL[[as.character(r$model)]],
    fmt_pct(r$pareto_k_gt_0_7, 1),
    fmt_num(r$pareto_k_max, 2),
    fmt_num(r$p_loo_per_dyad, 2),
    fmt_pct(r$bad_zero_share, 1),
    fmt_pct(r$zero_bad_rate, 1),
    fmt_pct(r$positive_bad_rate, 1)
  ))
}

lines <- c(
  "% Auto-generated by scripts/analysis/analyze_application_pareto_k.R",
  "% Do NOT edit by hand.",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{@{}llrrrrrr@{}}",
  "\\toprule",
  "Dataset & Model & $\\Pr(k>0.7)$ & max $k$ & $p_{\\mathrm{loo}}/M$ & Bad zero share & Zero bad rate & Positive bad rate \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  paste0(
    "\\caption{Dyad-level decomposition of the largest PSIS--LOO Pareto-",
    "$k$ tails. Here $M$ is the number of unordered dyads. The displayed ",
    "$\\Pr(k>0.7)$ column is the value stored in the application run; the ",
    "zero/positive decomposition uses recomputed pointwise $k$ values in ",
    "order to classify dyads. `Bad zero share' is the fraction of dyads with ",
    "$k>0.7$ that have no observed interaction in either direction, while ",
    "the last two columns condition on zero and positive dyads respectively.}"
  ),
  "\\label{tab:application-pareto-k-stress}",
  "\\end{table}"
)
writeLines(lines, file.path(OUT_DIR, "application_pareto_k_stress_summary.tex"))

cat("Wrote:\n",
    "  ", file.path(OUT_DIR, "application_pareto_k_stress_summary.csv"), "\n",
    "  ", file.path(OUT_DIR, "application_pareto_k_stress_summary.tex"), "\n",
    "  ", file.path(OUT_DIR, "application_pareto_k_top_dyads.csv"), "\n",
    sep = "")
