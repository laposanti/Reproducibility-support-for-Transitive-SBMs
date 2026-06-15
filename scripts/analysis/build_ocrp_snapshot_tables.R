#!/usr/bin/env Rscript
# Snapshot-specific tables for DemoOCRPvar/DemoKvar-style one-repetition runs.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(kableExtra)
})

snapshot_results_path <- Sys.getenv("RESULTS_PATH", Sys.getenv("SIM_RESULTS_PATH", ""))
snapshot_tables_dir <- Sys.getenv("SIM_TABLES_OUTPUT_DIR", "output/simulation/tables")

if (!nzchar(snapshot_results_path) || !file.exists(snapshot_results_path)) {
  stop("Snapshot table builder needs RESULTS_PATH or SIM_RESULTS_PATH.", call. = FALSE)
}

dir.create(snapshot_tables_dir, recursive = TRUE, showWarnings = FALSE)

snap_raw <- readr::read_csv(snapshot_results_path, show_col_types = FALSE)
if (!"rep_id" %in% names(snap_raw)) snap_raw$rep_id <- 1L
if (!"density" %in% names(snap_raw)) {
  snap_raw$density <- ifelse(snap_raw$kappa_mean <= median(snap_raw$kappa_mean, na.rm = TRUE),
                             "sparse", "dense")
}
if (!"scenario_id" %in% names(snap_raw)) {
  snap_raw$scenario_id <- paste(snap_raw$gen_model, snap_raw$density, snap_raw$hierch, sep = "_")
}
if (!"K_target" %in% names(snap_raw)) {
  snap_raw$K_target <- snap_raw$K_true
}

fmt_num <- function(x, digits = 3) {
  out <- ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), "--")
  as.character(out)
}

fmt_int <- function(x) {
  out <- ifelse(is.finite(x), sprintf("%d", as.integer(round(x))), "--")
  as.character(out)
}

fmt_mean_sd <- function(mu, sigma, digits = 1) {
  paste0(fmt_num(mu, digits), " (", fmt_num(sigma, digits), ")")
}

sig_2se <- function(delta, se) {
  ifelse(is.finite(delta) & is.finite(se) & se > 0 & abs(delta) > 2 * se, "\\checkmark", "")
}

save_latex <- function(df, file, caption, label, align = NULL) {
  if (is.null(align)) align <- paste(rep("l", ncol(df)), collapse = "")
  tab <- kableExtra::kbl(
    df,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = paste0(caption, "\\\\label{", label, "}"),
    align = align,
    linesep = ""
  ) |>
    kableExtra::kable_styling(full_width = FALSE, position = "center",
                              latex_options = c("hold_position"))
  kableExtra::save_kable(tab, file = file)
  message("Saved: ", file)
}

scenario_cols <- c("scenario_id", "gen_model", "density", "hierch", "rep_id")

snap <- snap_raw |>
  mutate(
    fit_model_clean = ifelse(fit_model %in% c("DC-SBM", "DCSBM"), "DC-SBM", fit_model),
    nvi = ifelse(K_true > 1, vi / log2(K_true), NA_real_),
    delta_elpd_report = if ("delta_elpd_best" %in% names(snap_raw)) delta_elpd_best else NA_real_,
    delta_elpd_se_report = if ("delta_elpd_best_se" %in% names(snap_raw)) delta_elpd_best_se else NA_real_
  ) |>
  group_by(across(all_of(scenario_cols))) |>
  mutate(
    best_elpd = max(elpd, na.rm = TRUE),
    delta_elpd = ifelse(is.finite(delta_elpd_report), delta_elpd_report, elpd - best_elpd),
    delta_elpd_se = delta_elpd_se_report,
    elpd_winner = is.finite(elpd) & abs(elpd - best_elpd) < 1e-9,
    best_vi = min(vi, na.rm = TRUE),
    vi_winner = is.finite(vi) & abs(vi - best_vi) < 1e-9,
    best_ari = max(ari, na.rm = TRUE),
    ari_winner = is.finite(ari) & abs(ari - best_ari) < 1e-9
  ) |>
  ungroup()

partition_table <- snap |>
  arrange(scenario_id, gen_model, density, hierch, rep_id, fit_model_clean) |>
  transmute(
    Scenario = scenario_id,
    Gen = gen_model,
    `$K^\\star$` = K_true,
    `$K_{target}$` = ifelse(is.na(K_target), "--", fmt_int(K_target)),
    Density = density,
    Hierarchy = hierch,
    Rep = rep_id,
    Model = fit_model_clean,
    `$\\hat K$` = fmt_int(K_hat),
    `$E[K\\mid y]$` = fmt_num(`K_p.mean`, 2),
    `ARI` = fmt_num(ari, 3),
    `VI` = fmt_num(vi, 3),
    `NVI` = fmt_num(nvi, 3),
    `OCRP tries` = fmt_int(ocrp_attempts),
    `true block sizes` = block_sizes_true
  )

partition_csv <- file.path(snapshot_tables_dir, "ocrp_snapshot_partition.csv")
readr::write_csv(partition_table, partition_csv)
save_latex(
  partition_table,
  file.path(snapshot_tables_dir, "ocrp_snapshot_partition.tex"),
  "OCRP snapshot: partition recovery and true block sizes for the one-repetition local run.",
  "tab:ocrp-snapshot-partition"
)

predictive_table <- snap |>
  arrange(scenario_id, gen_model, density, hierch, rep_id, delta_elpd, fit_model_clean) |>
  transmute(
    Scenario = scenario_id,
    Gen = gen_model,
    `$K^\\star$` = K_true,
    Density = density,
    Hierarchy = hierch,
    Rep = rep_id,
    Model = fit_model_clean,
    `ELPD (SE)` = if ("elpd_se" %in% names(snap)) {
      paste0(fmt_num(elpd, 1), " (", fmt_num(elpd_se, 1), ")")
    } else {
      fmt_num(elpd, 1)
    },
    `$\\Delta$ELPD (SE)` = if (any(is.finite(delta_elpd_se))) {
      paste0(fmt_num(delta_elpd, 1), " (", fmt_num(delta_elpd_se, 1), ")")
    } else {
      fmt_num(delta_elpd, 1)
    },
    `$|\\Delta|>2SE$` = if (any(is.finite(delta_elpd_se))) {
      sig_2se(delta_elpd, delta_elpd_se)
    } else {
      ""
    },
    `LOOIC` = fmt_num(looic, 1),
    `PSIS $k>0.7$` = fmt_num(pk_bad, 2),
    `ELPD winner` = ifelse(elpd_winner, "\\checkmark", "")
  )

predictive_csv <- file.path(snapshot_tables_dir, "ocrp_snapshot_predictive.csv")
readr::write_csv(predictive_table, predictive_csv)
save_latex(
  predictive_table,
  file.path(snapshot_tables_dir, "ocrp_snapshot_predictive.tex"),
  "OCRP snapshot: predictive comparison by scenario. $\\Delta$ELPD is relative to the best fitted model in each scenario.",
  "tab:ocrp-snapshot-predictive"
)

predictive_summary <- snap |>
  group_by(scenario_id, gen_model, density, hierch, fit_model_clean) |>
  summarise(
    n_rep = n_distinct(rep_id),
    mean_K_true = mean(K_true, na.rm = TRUE),
    min_K_true = min(K_true, na.rm = TRUE),
    max_K_true = max(K_true, na.rm = TRUE),
    mean_elpd = mean(elpd, na.rm = TRUE),
    sd_elpd = sd(elpd, na.rm = TRUE),
    mean_delta = mean(delta_elpd, na.rm = TRUE),
    sd_delta = sd(delta_elpd, na.rm = TRUE),
    significant_losses = sum(is.finite(delta_elpd) & is.finite(delta_elpd_se) &
                               delta_elpd < 0 & abs(delta_elpd) > 2 * delta_elpd_se,
                             na.rm = TRUE),
    elpd_wins = sum(elpd_winner, na.rm = TRUE),
    mean_vi = mean(vi, na.rm = TRUE),
    mean_K_hat = mean(K_hat, na.rm = TRUE),
    mean_pk_bad = mean(pk_bad, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(scenario_id, gen_model, density, hierch, mean_delta, fit_model_clean) |>
  transmute(
    Scenario = scenario_id,
    Gen = gen_model,
    `$K^\\star$ mean [range]` = paste0(fmt_num(mean_K_true, 1), " [", fmt_int(min_K_true), ",", fmt_int(max_K_true), "]"),
    Density = density,
    Hierarchy = hierch,
    Model = fit_model_clean,
    Reps = n_rep,
    `ELPD mean (SD)` = fmt_mean_sd(mean_elpd, sd_elpd, 1),
    `$\\Delta$ELPD mean (SD)` = fmt_mean_sd(mean_delta, sd_delta, 1),
    `Sig. loss` = paste0(significant_losses, "/", n_rep),
    `ELPD wins` = paste0(elpd_wins, "/", n_rep),
    `mean VI` = fmt_num(mean_vi, 3),
    `mean $\\hat K$` = fmt_num(mean_K_hat, 2),
    `mean PSIS $k>0.7$` = fmt_num(mean_pk_bad, 3)
  )

predictive_summary_csv <- file.path(snapshot_tables_dir, "ocrp_snapshot_predictive_summary.csv")
readr::write_csv(predictive_summary, predictive_summary_csv)
save_latex(
  predictive_summary,
  file.path(snapshot_tables_dir, "ocrp_snapshot_predictive_summary.tex"),
  "OCRP snapshot: predictive comparison summarized across repetitions. Standard deviations are across repetitions.",
  "tab:ocrp-snapshot-predictive-summary"
)

pairwise_wide <- snap |>
  mutate(model_key = ifelse(fit_model_clean == "DC-SBM", "DC", fit_model_clean)) |>
  select(all_of(scenario_cols), K_true, model_key, elpd, vi, K_hat,
         any_of(c("elpd_wst_minus_sst", "elpd_wst_minus_sst_se", "elpd_wst_minus_sst_z",
                  "elpd_wst_minus_dc", "elpd_wst_minus_dc_se", "elpd_wst_minus_dc_z",
                  "elpd_sst_minus_dc", "elpd_sst_minus_dc_se", "elpd_sst_minus_dc_z"))) |>
  tidyr::pivot_wider(
    names_from = model_key,
    values_from = c(elpd, vi, K_hat),
    names_sep = "_"
  ) |>
  arrange(scenario_id, gen_model, density, hierch, rep_id)

has_wst_sst_se <- "elpd_wst_minus_sst_se" %in% names(pairwise_wide)
has_wst_dc_se <- "elpd_wst_minus_dc_se" %in% names(pairwise_wide)
has_sst_dc_se <- "elpd_sst_minus_dc_se" %in% names(pairwise_wide)

pairwise_table <- pairwise_wide |>
  transmute(
    Scenario = scenario_id,
    Gen = gen_model,
    `$K^\\star$` = K_true,
    Density = density,
    Hierarchy = hierch,
    Rep = rep_id,
    `ELPD WST$-$SST (SE)` = if (has_wst_sst_se) {
      paste0(fmt_num(elpd_wst_minus_sst, 1), " (", fmt_num(elpd_wst_minus_sst_se, 1), ")")
    } else {
      fmt_num(elpd_WST - elpd_SST, 1)
    },
    `Sig. WST$-$SST` = if (has_wst_sst_se) sig_2se(elpd_wst_minus_sst, elpd_wst_minus_sst_se) else "",
    `ELPD WST$-$DC (SE)` = if (has_wst_dc_se) {
      paste0(fmt_num(elpd_wst_minus_dc, 1), " (", fmt_num(elpd_wst_minus_dc_se, 1), ")")
    } else {
      fmt_num(elpd_WST - elpd_DC, 1)
    },
    `Sig. WST$-$DC` = if (has_wst_dc_se) sig_2se(elpd_wst_minus_dc, elpd_wst_minus_dc_se) else "",
    `ELPD SST$-$DC (SE)` = if (has_sst_dc_se) {
      paste0(fmt_num(elpd_sst_minus_dc, 1), " (", fmt_num(elpd_sst_minus_dc_se, 1), ")")
    } else {
      fmt_num(elpd_SST - elpd_DC, 1)
    },
    `Sig. SST$-$DC` = if (has_sst_dc_se) sig_2se(elpd_sst_minus_dc, elpd_sst_minus_dc_se) else "",
    `VI WST` = fmt_num(vi_WST, 3),
    `VI SST` = fmt_num(vi_SST, 3),
    `VI DC` = fmt_num(vi_DC, 3),
    `$\\hat K$ W/S/DC` = paste(fmt_int(K_hat_WST), fmt_int(K_hat_SST), fmt_int(K_hat_DC), sep = "/")
  )

pairwise_csv <- file.path(snapshot_tables_dir, "ocrp_snapshot_pairwise_deltas.csv")
readr::write_csv(pairwise_table, pairwise_csv)
save_latex(
  pairwise_table,
  file.path(snapshot_tables_dir, "ocrp_snapshot_pairwise_deltas.tex"),
  "OCRP snapshot: direct pairwise ELPD deltas by scenario. Positive ELPD differences favor the first model in the column name.",
  "tab:ocrp-snapshot-pairwise-deltas"
)

pairwise_for_summary <- pairwise_wide |>
  mutate(
    wst_sst_delta_report = if ("elpd_wst_minus_sst" %in% names(pairwise_wide)) {
      elpd_wst_minus_sst
    } else {
      elpd_WST - elpd_SST
    },
    wst_sst_abs_z_report = if ("elpd_wst_minus_sst_z" %in% names(pairwise_wide)) {
      elpd_wst_minus_sst_z
    } else {
      NA_real_
    },
    wst_dc_delta_report = if ("elpd_wst_minus_dc" %in% names(pairwise_wide)) {
      elpd_wst_minus_dc
    } else {
      elpd_WST - elpd_DC
    },
    wst_dc_abs_z_report = if ("elpd_wst_minus_dc_z" %in% names(pairwise_wide)) {
      elpd_wst_minus_dc_z
    } else {
      NA_real_
    },
    sst_dc_delta_report = if ("elpd_sst_minus_dc" %in% names(pairwise_wide)) {
      elpd_sst_minus_dc
    } else {
      elpd_SST - elpd_DC
    },
    sst_dc_abs_z_report = if ("elpd_sst_minus_dc_z" %in% names(pairwise_wide)) {
      elpd_sst_minus_dc_z
    } else {
      NA_real_
    },
    wst_sst_sig = if ("elpd_wst_minus_sst_se" %in% names(pairwise_wide)) {
      is.finite(wst_sst_delta_report) & is.finite(elpd_wst_minus_sst_se) &
        abs(wst_sst_delta_report) > 2 * elpd_wst_minus_sst_se
    } else {
      FALSE
    },
    wst_dc_sig = if ("elpd_wst_minus_dc_se" %in% names(pairwise_wide)) {
      is.finite(wst_dc_delta_report) & is.finite(elpd_wst_minus_dc_se) &
        abs(wst_dc_delta_report) > 2 * elpd_wst_minus_dc_se
    } else {
      FALSE
    },
    sst_dc_sig = if ("elpd_sst_minus_dc_se" %in% names(pairwise_wide)) {
      is.finite(sst_dc_delta_report) & is.finite(elpd_sst_minus_dc_se) &
        abs(sst_dc_delta_report) > 2 * elpd_sst_minus_dc_se
    } else {
      FALSE
    }
  )

pairwise_summary <- pairwise_for_summary |>
  group_by(scenario_id, gen_model, density, hierch) |>
  summarise(
    n_rep = n_distinct(rep_id),
    mean_K_true = mean(K_true, na.rm = TRUE),
    min_K_true = min(K_true, na.rm = TRUE),
    max_K_true = max(K_true, na.rm = TRUE),
    wst_sst_mean = mean(wst_sst_delta_report, na.rm = TRUE),
    wst_sst_sd = sd(wst_sst_delta_report, na.rm = TRUE),
    wst_sst_pos_sig = sum(wst_sst_sig & wst_sst_delta_report > 0, na.rm = TRUE),
    wst_sst_neg_sig = sum(wst_sst_sig & wst_sst_delta_report < 0, na.rm = TRUE),
    wst_dc_mean = mean(wst_dc_delta_report, na.rm = TRUE),
    wst_dc_sd = sd(wst_dc_delta_report, na.rm = TRUE),
    wst_dc_pos_sig = sum(wst_dc_sig & wst_dc_delta_report > 0, na.rm = TRUE),
    wst_dc_neg_sig = sum(wst_dc_sig & wst_dc_delta_report < 0, na.rm = TRUE),
    sst_dc_mean = mean(sst_dc_delta_report, na.rm = TRUE),
    sst_dc_sd = sd(sst_dc_delta_report, na.rm = TRUE),
    sst_dc_pos_sig = sum(sst_dc_sig & sst_dc_delta_report > 0, na.rm = TRUE),
    sst_dc_neg_sig = sum(sst_dc_sig & sst_dc_delta_report < 0, na.rm = TRUE),
    mean_vi_wst = mean(vi_WST, na.rm = TRUE),
    mean_vi_sst = mean(vi_SST, na.rm = TRUE),
    mean_vi_dc = mean(vi_DC, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(scenario_id, gen_model, density, hierch) |>
  transmute(
    Scenario = scenario_id,
    Gen = gen_model,
    `$K^\\star$ mean [range]` = paste0(fmt_num(mean_K_true, 1), " [", fmt_int(min_K_true), ",", fmt_int(max_K_true), "]"),
    Density = density,
    Hierarchy = hierch,
    Reps = n_rep,
    `ELPD WST$-$SST mean (SD)` = fmt_mean_sd(wst_sst_mean, wst_sst_sd, 1),
    `Sig. WST/SST` = paste0(wst_sst_pos_sig, "/", wst_sst_neg_sig),
    `ELPD WST$-$DC mean (SD)` = fmt_mean_sd(wst_dc_mean, wst_dc_sd, 1),
    `Sig. WST/DC` = paste0(wst_dc_pos_sig, "/", wst_dc_neg_sig),
    `ELPD SST$-$DC mean (SD)` = fmt_mean_sd(sst_dc_mean, sst_dc_sd, 1),
    `Sig. SST/DC` = paste0(sst_dc_pos_sig, "/", sst_dc_neg_sig),
    `mean VI W/S/DC` = paste(fmt_num(mean_vi_wst, 3), fmt_num(mean_vi_sst, 3),
                              fmt_num(mean_vi_dc, 3), sep = "/")
  )

pairwise_summary_csv <- file.path(snapshot_tables_dir, "ocrp_snapshot_pairwise_summary.csv")
readr::write_csv(pairwise_summary, pairwise_summary_csv)
save_latex(
  pairwise_summary,
  file.path(snapshot_tables_dir, "ocrp_snapshot_pairwise_summary.tex"),
  "OCRP snapshot: direct pairwise ELPD deltas summarized across repetitions. Positive ELPD differences favor the first model in the column name.",
  "tab:ocrp-snapshot-pairwise-summary"
)

order_table <- snap |>
  arrange(scenario_id, gen_model, density, hierch, rep_id, fit_model_clean) |>
  transmute(
    Scenario = scenario_id,
    Gen = gen_model,
    `$K^\\star$` = K_true,
    Density = density,
    Hierarchy = hierch,
    Rep = rep_id,
    Model = fit_model_clean,
    `$\\widehat\\zeta_W^{blk}$` = fmt_num(thetaW_block_emp_mean, 3),
    `$\\widehat\\zeta_S^{blk}$` = fmt_num(thetaS_block_emp_mean, 3),
    `Block coverage` = fmt_num(coverage_block_emp_avg, 3),
    `Violation mean` = fmt_num(violation_rate_mean, 3),
    `Violation $\\hat z$` = fmt_num(violation_rate_zhat, 3),
    `$p_{post}(WST)$` = fmt_num(p_post_wst, 3),
    `$p_{post}(SST)$` = fmt_num(p_post_sst, 3)
  )

order_csv <- file.path(snapshot_tables_dir, "ocrp_snapshot_order.csv")
readr::write_csv(order_table, order_csv)
save_latex(
  order_table,
  file.path(snapshot_tables_dir, "ocrp_snapshot_order.tex"),
  "OCRP snapshot: block-level hierarchy diagnostics and posterior order probabilities.",
  "tab:ocrp-snapshot-order"
)

pick_winner <- function(metric, better = c("max", "min")) {
  better <- match.arg(better)
  function(df) {
    vals <- df[[metric]]
    ok <- is.finite(vals)
    if (!any(ok)) return(NA_character_)
    target <- if (better == "max") max(vals[ok]) else min(vals[ok])
    paste(df$fit_model_clean[ok & abs(vals - target) < 1e-9], collapse = "/")
  }
}

winner_table <- snap |>
  group_by(across(all_of(scenario_cols))) |>
  summarise(
    K_true = first(K_true),
    `ELPD winner` = pick_winner("elpd", "max")(pick(everything())),
    `VI winner` = pick_winner("vi", "min")(pick(everything())),
    `ARI winner` = pick_winner("ari", "max")(pick(everything())),
    `Exact K models` = paste(fit_model_clean[K_hat == K_true], collapse = "/"),
    `Block sizes` = first(block_sizes_true),
    `theta OCRP` = first(theta_ocrp),
    .groups = "drop"
  ) |>
  arrange(scenario_id, gen_model, density, hierch, rep_id) |>
  transmute(
    Scenario = scenario_id,
    Gen = gen_model,
    `$K^\\star$` = K_true,
    Density = density,
    Hierarchy = hierch,
    Rep = rep_id,
    `ELPD winner`,
    `VI winner`,
    `ARI winner`,
    `Exact K models` = ifelse(nchar(`Exact K models`) > 0, `Exact K models`, "--"),
    `Block sizes`,
    `$\\vartheta$` = fmt_num(`theta OCRP`, 2)
  )

winner_csv <- file.path(snapshot_tables_dir, "ocrp_snapshot_winners.csv")
readr::write_csv(winner_table, winner_csv)
save_latex(
  winner_table,
  file.path(snapshot_tables_dir, "ocrp_snapshot_winners.tex"),
  "OCRP snapshot: scenario-level winners for predictive performance, VI, ARI, and exact $K$ recovery.",
  "tab:ocrp-snapshot-winners"
)

winner_counts <- snap |>
  group_by(across(all_of(scenario_cols))) |>
  mutate(
    elpd_win = elpd_winner,
    vi_win = vi_winner,
    ari_win = ari_winner,
    exact_K = K_hat == K_true
  ) |>
  ungroup() |>
  group_by(gen_model, fit_model_clean) |>
  summarise(
    `ELPD wins` = sum(elpd_win, na.rm = TRUE),
    `VI wins` = sum(vi_win, na.rm = TRUE),
    `ARI wins` = sum(ari_win, na.rm = TRUE),
    `Exact K` = sum(exact_K, na.rm = TRUE),
    `Scenarios` = n_distinct(paste(scenario_id, gen_model, density, hierch, rep_id)),
    .groups = "drop"
  ) |>
  arrange(gen_model, fit_model_clean) |>
  transmute(
    Gen = gen_model,
    Model = fit_model_clean,
    `ELPD wins`,
    `VI wins`,
    `ARI wins`,
    `Exact K`,
    Scenarios
  )

counts_csv <- file.path(snapshot_tables_dir, "ocrp_snapshot_winner_counts.csv")
readr::write_csv(winner_counts, counts_csv)
save_latex(
  winner_counts,
  file.path(snapshot_tables_dir, "ocrp_snapshot_winner_counts.tex"),
  "OCRP snapshot: winner counts across one-repetition scenarios.",
  "tab:ocrp-snapshot-winner-counts"
)

message("Snapshot table CSVs written to: ", snapshot_tables_dir)
