#!/usr/bin/env Rscript
# Pointwise LOO diagnostics for saved dyad-level OSBM/DC-SBM comparisons.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(kableExtra)
})

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L && nzchar(args[[1L]])) {
  args[[1L]]
} else {
  Sys.getenv("SIM_RESULTS_PATH", "")
}

if (!nzchar(input_path)) {
  stop("Pass a simulation results CSV or run directory.", call. = FALSE)
}

input_path <- normalizePath(input_path, winslash = "/", mustWork = TRUE)
if (dir.exists(input_path)) {
  run_dir <- if (basename(input_path) == "loo_pointwise") dirname(input_path) else input_path
  run_id <- basename(run_dir)
} else {
  run_dir <- dirname(input_path)
  base_name <- basename(input_path)
  run_id <- sub(
    "^full_simulation_crossfit_(final|progress)_(Demo(Kvar|OCRPvar)_run_[0-9]{8}_[0-9]{6})\\.csv$",
    "\\2",
    base_name
  )
  if (identical(run_id, base_name)) {
    run_id <- sub(
      "^snapshot_after_rep1_(Demo(Kvar|OCRPvar)_run_[0-9]{8}_[0-9]{6})\\.csv$",
      "\\1",
      base_name
    )
  }
  if (identical(run_id, base_name)) run_id <- basename(run_dir)
}

loo_dir <- file.path(run_dir, "loo_pointwise")
if (!dir.exists(loo_dir)) {
  stop("No loo_pointwise directory found under: ", run_dir, call. = FALSE)
}

run_date <- sub("^.*_([0-9]{8})_[0-9]{6}$", "\\1", run_id)
if (!grepl("^[0-9]{8}$", run_date)) run_date <- format(Sys.Date(), "%Y%m%d")
run_date_iso <- paste0(substr(run_date, 1, 4), "-", substr(run_date, 5, 6), "-", substr(run_date, 7, 8))
tag <- paste0(run_date_iso, "_", run_id)

tables_root <- Sys.getenv("SIM_TABLES_OUTPUT_DIR", file.path("output", "simulation", "tables", tag))
plots_root <- Sys.getenv("SIM_PLOTS_OUTPUT_DIR", file.path("output", "simulation", "plots", tag))
tables_dir <- file.path(tables_root, "loo_pointwise")
plots_dir <- file.path(plots_root, "loo_pointwise")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(loo_dir, pattern = "\\.csv$", full.names = TRUE))
if (!length(files)) stop("No pointwise LOO CSV files found in: ", loo_dir, call. = FALSE)

message("Reading ", length(files), " pointwise LOO files from: ", loo_dir)

pointwise <- bind_rows(lapply(files, function(path) {
  readr::read_csv(path, show_col_types = FALSE) |>
    mutate(source_file = basename(path), .before = 1)
}))

delta_cols <- grep("^delta_elpd_", names(pointwise), value = TRUE)
if (!length(delta_cols)) stop("No delta_elpd_* columns found in pointwise files.", call. = FALSE)

fmt_num <- function(x, digits = 2) {
  ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), "--")
}

sig_2se <- function(delta, se) {
  ifelse(is.finite(delta) & is.finite(se) & se > 0 & abs(delta) > 2 * se, "\\checkmark", "")
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_=-]+", "_", x)
}

scenario_cols <- intersect(
  c("source_file", "gen_model", "K_true", "density", "hierch", "rep_id", "scenario_id",
    "theta_ocrp", "dgp_partition"),
  names(pointwise)
)

loo_se_sum <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sqrt(length(x) * stats::var(x))
}

long <- pointwise |>
  mutate(
    block_low = pmin(zi_true, zj_true),
    block_high = pmax(zi_true, zj_true),
    block_pair = paste(block_low, block_high, sep = "-"),
    block_distance = abs(zi_true - zj_true),
    total_edges = Aij + Aji,
    net_edges = Aij - Aji,
    observed_share_ij = ifelse(Nij > 0, Aij / Nij, NA_real_)
  ) |>
  pivot_longer(
    cols = all_of(delta_cols),
    names_to = "comparison",
    values_to = "pointwise_delta"
  ) |>
  mutate(
    comparison = sub("^delta_elpd_", "", comparison),
    first_model = sub("_minus_.*$", "", comparison),
    second_model = sub("^.*_minus_", "", comparison),
    favoured_model = ifelse(pointwise_delta >= 0, first_model, second_model)
  )

overall <- long |>
  group_by(across(all_of(c(scenario_cols, "comparison", "first_model", "second_model")))) |>
  summarise(
    n_dyads = n(),
    sum_delta = sum(pointwise_delta, na.rm = TRUE),
    se_delta = loo_se_sum(pointwise_delta),
    z_abs = ifelse(is.finite(se_delta) & se_delta > 0, abs(sum_delta / se_delta), NA_real_),
    significant_2se = is.finite(se_delta) & se_delta > 0 & abs(sum_delta) > 2 * se_delta,
    mean_delta = mean(pointwise_delta, na.rm = TRUE),
    sd_delta = sd(pointwise_delta, na.rm = TRUE),
    positive_fraction = mean(pointwise_delta > 0, na.rm = TRUE),
    max_abs_pointwise_delta = max(abs(pointwise_delta), na.rm = TRUE),
    mean_Nij = mean(Nij, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(source_file, comparison)

by_distance <- long |>
  group_by(across(all_of(c(scenario_cols, "comparison", "block_distance")))) |>
  summarise(
    n_dyads = n(),
    sum_delta = sum(pointwise_delta, na.rm = TRUE),
    se_delta = loo_se_sum(pointwise_delta),
    z_abs = ifelse(is.finite(se_delta) & se_delta > 0, abs(sum_delta / se_delta), NA_real_),
    significant_2se = is.finite(se_delta) & se_delta > 0 & abs(sum_delta) > 2 * se_delta,
    mean_delta = mean(pointwise_delta, na.rm = TRUE),
    positive_fraction = mean(pointwise_delta > 0, na.rm = TRUE),
    mean_Nij = mean(Nij, na.rm = TRUE),
    mean_total_edges = mean(total_edges, na.rm = TRUE),
    mean_net_edges = mean(net_edges, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(source_file, comparison, block_distance)

by_block_pair <- long |>
  group_by(across(all_of(c(scenario_cols, "comparison", "block_low", "block_high", "block_pair",
                           "block_distance")))) |>
  summarise(
    n_dyads = n(),
    sum_delta = sum(pointwise_delta, na.rm = TRUE),
    se_delta = loo_se_sum(pointwise_delta),
    z_abs = ifelse(is.finite(se_delta) & se_delta > 0, abs(sum_delta / se_delta), NA_real_),
    significant_2se = is.finite(se_delta) & se_delta > 0 & abs(sum_delta) > 2 * se_delta,
    mean_delta = mean(pointwise_delta, na.rm = TRUE),
    positive_fraction = mean(pointwise_delta > 0, na.rm = TRUE),
    mean_Nij = mean(Nij, na.rm = TRUE),
    mean_Aij = mean(Aij, na.rm = TRUE),
    mean_Aji = mean(Aji, na.rm = TRUE),
    mean_total_edges = mean(total_edges, na.rm = TRUE),
    mean_net_edges = mean(net_edges, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(source_file, comparison, desc(abs(sum_delta)))

top_dyads <- long |>
  group_by(across(all_of(c(scenario_cols, "comparison")))) |>
  arrange(desc(abs(pointwise_delta)), .by_group = TRUE) |>
  slice_head(n = 20) |>
  ungroup() |>
  select(any_of(c(
    scenario_cols, "comparison", "favoured_model", "dyad_id", "i", "j",
    "zi_true", "zj_true", "block_pair", "block_distance", "Aij", "Aji", "Nij",
    "total_edges", "net_edges", "observed_share_ij", "pointwise_delta"
  )))

readr::write_csv(overall, file.path(tables_dir, "loo_pointwise_overall.csv"))
readr::write_csv(by_distance, file.path(tables_dir, "loo_pointwise_by_distance.csv"))
readr::write_csv(by_block_pair, file.path(tables_dir, "loo_pointwise_by_block_pair.csv"))
readr::write_csv(top_dyads, file.path(tables_dir, "loo_pointwise_top_dyads.csv"))

overall_tex <- overall |>
  transmute(
    Gen = gen_model,
    Rep = rep_id,
    Scenario = scenario_id,
    Comparison = comparison,
    `sum Delta ELPD` = fmt_num(sum_delta, 1),
    `SE` = fmt_num(se_delta, 1),
    `Sig.` = sig_2se(sum_delta, se_delta),
    `mean Delta` = fmt_num(mean_delta, 3),
    `% positive` = fmt_num(100 * positive_fraction, 1)
  )

overall_kbl <- kableExtra::kbl(
  overall_tex,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Pointwise LOO decomposition by scenario. Positive deltas favor the first model in the comparison."
) |>
  kableExtra::kable_styling(full_width = FALSE, position = "center",
                            latex_options = c("hold_position"))
kableExtra::save_kable(overall_kbl, file.path(tables_dir, "loo_pointwise_overall.tex"))

top_tex <- top_dyads |>
  group_by(source_file, comparison) |>
  slice_head(n = 8) |>
  ungroup() |>
  transmute(
    Gen = gen_model,
    Rep = rep_id,
    Comparison = comparison,
    Dyad = paste(i, j, sep = "-"),
    Blocks = block_pair,
    `Aij/Aji/N` = paste(Aij, Aji, Nij, sep = "/"),
    Favours = favoured_model,
    `Delta ELPD` = fmt_num(pointwise_delta, 2)
  )

top_kbl <- kableExtra::kbl(
  top_tex,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Largest absolute dyad-level LOO deltas. Positive deltas favor the first model in the comparison."
) |>
  kableExtra::kable_styling(full_width = FALSE, position = "center",
                            latex_options = c("hold_position"))
kableExtra::save_kable(top_kbl, file.path(tables_dir, "loo_pointwise_top_dyads.tex"))

heatmap_data <- by_block_pair |>
  filter(comparison == "WST_minus_SST")

if (nrow(heatmap_data) > 0L) {
  max_abs <- max(abs(heatmap_data$sum_delta), na.rm = TRUE)
  for (sf in unique(heatmap_data$source_file)) {
    plot_data <- heatmap_data |> filter(source_file == sf)
    title <- paste0(
      unique(plot_data$gen_model), " rep ", unique(plot_data$rep_id),
      ": pointwise LOO WST-SST by true block pair"
    )
    p <- ggplot(plot_data, aes(x = factor(block_low), y = factor(block_high), fill = sum_delta)) +
      geom_tile(color = "white", linewidth = 0.4) +
      geom_text(aes(label = sprintf("%.1f", sum_delta)), size = 3) +
      scale_fill_gradient2(
        low = "#B2182B", mid = "white", high = "#2166AC",
        midpoint = 0, limits = c(-max_abs, max_abs),
        name = "sum Delta\nWST-SST"
      ) +
      coord_equal() +
      labs(title = title, x = "lower true block", y = "higher true block") +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank())

    out_base <- file.path(plots_dir, paste0("loo_heatmap_WST_minus_SST_", safe_name(tools::file_path_sans_ext(sf))))
    ggsave(paste0(out_base, ".pdf"), p, width = 6.5, height = 5.2)
    ggsave(paste0(out_base, ".png"), p, width = 6.5, height = 5.2, dpi = 200)
  }
}

message("Saved pointwise LOO tables to: ", tables_dir)
message("Saved pointwise LOO plots to:  ", plots_dir)
