#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: build_extended_sim_tables.R <input_csv> <output_dir>")
}

input_csv <- args[[1]]
output_dir <- args[[2]]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_csv(input_csv, show_col_types = FALSE) |>
  mutate(
    fit_model = case_when(
      fit_model == "DC-SBM" ~ "DC",
      fit_model == "SST" ~ "SST",
      fit_model == "WST" ~ "WST",
      TRUE ~ fit_model
    ),
    gen_model = as.character(gen_model),
    K_true = as.integer(K_true),
    kappa_mean = as.numeric(kappa_mean),
    psi_mean = as.numeric(psi_mean)
  )

fit_levels <- c("WST", "SST", "DC")

format_num <- function(x, digits = 3) sprintf(paste0("%.", digits, "f"), x)

write_overall_table <- function(df, out_file) {
  s <- df |>
    filter(fit_model %in% fit_levels) |>
    group_by(gen_model, K_true, fit_model) |>
    summarise(
      mean_looic = mean(looic, na.rm = TRUE),
      sd_looic = sd(looic, na.rm = TRUE),
      mean_vi = mean(vi, na.rm = TRUE),
      mean_ari = mean(ari, na.rm = TRUE),
      pr_keq = mean(K_hat == K_true, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    mutate(fit_model = factor(fit_model, levels = fit_levels)) |>
    arrange(gen_model, K_true, fit_model)

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Extended run summary by generating model and true block count. Values are means across all hierarchy and density settings in DemoKvar\\_run\\_20260604\\_165415 (4 replicates per cell). Lower LOOIC and VI are better; higher ARI and $\\Pr(\\hat K=K^\\star)$ are better.}",
    "\\label{tab:extended-overall-by-k}",
    "\\begin{tabular}{lllrrrrr}",
    "\\toprule",
    "Gen. & $K^\\star$ & Fit & Mean LOOIC & SD LOOIC & Mean VI & Mean ARI & $\\Pr(\\hat K=K^\\star)$ \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(s))) {
    row <- s[i, ]
    lines <- c(lines, sprintf(
      "%s & %d & %s & %s & %s & %s & %s & %s \\\\",
      row$gen_model,
      row$K_true,
      as.character(row$fit_model),
      format_num(row$mean_looic, 1),
      format_num(row$sd_looic, 1),
      format_num(row$mean_vi, 3),
      format_num(row$mean_ari, 3),
      format_num(row$pr_keq, 2)
    ))
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, out_file)
}

write_k8_focus_table <- function(df, out_file) {
  focus <- df |>
    filter(fit_model %in% fit_levels, K_true == 8, abs(kappa_mean - 0.75) < 1e-12, psi_mean %in% c(0.2, 1.3)) |>
    group_by(gen_model, psi_mean, fit_model) |>
    summarise(
      mean_looic = mean(looic, na.rm = TRUE),
      sd_looic = sd(looic, na.rm = TRUE),
      mean_vi = mean(vi, na.rm = TRUE),
      mean_ari = mean(ari, na.rm = TRUE),
      .groups = "drop"
    )

  winners <- focus |>
    group_by(gen_model, psi_mean) |>
    summarise(winner_looic = min(mean_looic, na.rm = TRUE), .groups = "drop")

  focus <- focus |>
    left_join(winners, by = c("gen_model", "psi_mean")) |>
    mutate(delta_looic = mean_looic - winner_looic) |>
    mutate(fit_model = factor(fit_model, levels = fit_levels)) |>
    arrange(gen_model, psi_mean, fit_model)

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Scenario check used in the manuscript update: $K^\\star=8$, $\\bar\\kappa=0.75$. $\\Delta$LOOIC is computed against the best mean LOOIC within each $(\\mathrm{Gen.},\\psi^\\star)$ block.}",
    "\\label{tab:extended-k8-focus}",
    "\\begin{tabular}{lllrrrr}",
    "\\toprule",
    "Gen. & $\\psi^\\star$ & Fit & Mean LOOIC & $\\Delta$LOOIC & Mean VI & Mean ARI \\\\",
    "\\midrule"
  )

  for (i in seq_len(nrow(focus))) {
    row <- focus[i, ]
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      row$gen_model,
      format_num(row$psi_mean, 1),
      as.character(row$fit_model),
      format_num(row$mean_looic, 1),
      format_num(row$delta_looic, 1),
      format_num(row$mean_vi, 3),
      format_num(row$mean_ari, 3)
    ))
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, out_file)
}

write_overall_table(d, file.path(output_dir, "extended_overall_by_k.tex"))
write_k8_focus_table(d, file.path(output_dir, "extended_k8_focus.tex"))

message("Wrote: ", file.path(output_dir, "extended_overall_by_k.tex"))
message("Wrote: ", file.path(output_dir, "extended_k8_focus.tex"))