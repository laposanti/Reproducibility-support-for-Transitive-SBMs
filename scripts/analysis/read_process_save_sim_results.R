# ============================================================
# LaTeX tables (booktabs) for OSBM/DC-SBM simulations
# Main text + Appendix pack
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(glue)
  library(kableExtra)
  library(purrr)
  library(stringr)
})

`%+%` <- function(a, b) paste0(a, b)

# ---------- SETTINGS -----------------------------------------
compact_results_path <- Sys.getenv(
  "COMPACT_RESULTS_PATH",
  "./output/simulation/raw/full_simulation_crossfit_final_DemoKvar_run_20260302_153429.csv"
)
results_path <- Sys.getenv("RESULTS_PATH", compact_results_path)
out_dir_main <- Sys.getenv("SIM_TABLES_OUTPUT_DIR", "output/simulation/tables")
out_dir_appx <- file.path(out_dir_main, "appendix")
digits_small <- 3
digits_mid   <- 2

dir.create(out_dir_main, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_appx, showWarnings = FALSE, recursive = TRUE)

# ---------- HELPERS ------------------------------------------

kbl_ltx <- function(df, caption, label, align = NULL) {
  if (is.null(align)) align <- paste(rep("r", ncol(df)), collapse = "")
  kbl(df, format = "latex", booktabs = TRUE,
      caption = glue("{caption}\\\\\\label{{{label}}}"),
      align = align, escape = FALSE, linesep = "") |>
    kable_styling(full_width = FALSE, position = "center",
                  latex_options = c("hold_position"))
}

save_tab <- function(kbl_obj, file) {
  save_kable(kbl_obj, file = file)
  message("Saved: ", file)
}

pm <- function(x, d = 3) {
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  sprintf(paste0("%.", d, "f (%.", d, "f)"), m, s)
}

pct <- function(x, d = 1) {
  paste0(sprintf(paste0("%.", d, "f"), 100 * x), "\\%")
}

ci_str <- function(m, lo, hi, d = 3) {
  sprintf(paste0("%.", d, "f [%.", d, "f--%.", d, "f]"), m, lo, hi)
}

diag_flag <- function(ess, mcse, lo, hi, ess_thr = 400, mcse_rel_thr = 0.25) {
  width <- hi - lo
  too_small_ess <- !is.na(ess) & ess < ess_thr
  too_big_mcse  <- !is.na(mcse) & !is.na(width) & width > 0 & (mcse > mcse_rel_thr * width)
  flag <- ifelse(too_small_ess | too_big_mcse, "\\textsuperscript{\\dagger}", "")
  ifelse(is.na(flag), "", flag)
}

psisk_share_col <- function(df) {
  if ("pareto_k_share" %in% names(df)) {
    df$pareto_k_share
  } else if ("pareto_k_max" %in% names(df)) {
    as.numeric(df$pareto_k_max > 0.7)
  } else {
    NA_real_
  }
}

has_cols <- function(df, cols) all(cols %in% names(df))

fmt_wins <- function(x, n) {
  ifelse(abs(x - round(x)) < 1e-10,
         sprintf("%d/%d", as.integer(round(x)), as.integer(round(n))),
         sprintf("%.1f/%d", x, as.integer(round(n))))
}

# ---------- LOAD ---------------------------------------------
raw <- readr::read_csv(results_path, show_col_types = FALSE)

if (!"K" %in% names(raw)) {
  if ("K_true" %in% names(raw)) {
    raw <- raw |> mutate(K = K_true)
  } else {
    stop("Missing K column: expected `K` or `K_true` in results file: ", results_path)
  }
}

if (!"rep" %in% names(raw)) {
  rep_candidates <- c("rep_id", "rep_id.", "replicate")
  rep_col <- rep_candidates[rep_candidates %in% names(raw)][1]
  if (!is.na(rep_col)) {
    raw <- raw |> mutate(rep = .data[[rep_col]])
  } else {
    raw <- raw |> mutate(rep = row_number())
  }
}

raw <- raw |>
  mutate(
    gen_model = factor(gen_model, levels = c("WST","SST")),
    fit_model = factor(fit_model, levels = c("WST","SST","DC-SBM")),
    K = as.integer(K),
    rep = as.integer(rep)
  ) |>
  mutate(
    psis_tail_share = psisk_share_col(pick(everything())),
    vi_norm = vi / log2(K)
  )

# ---------- MAIN-TEXT TABLES ---------------------------------

make_table1 <- function(df_gen, gen) {
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  tab <- df_gen |>
    group_by(K, fit_model) |>
    summarise(
      `ARI mean (sd)`       = pm(ari, d = digits_small),
      `Pr(ARI ≥ 0.95)`      = pct(mean(ari >= 0.95, na.rm = TRUE), d = 1),
      `MAE $\\eta$`         = pm(mae_eta,   d = digits_small),
      `MAE $\\kappa$`       = pm(mae_kappa, d = digits_small),
      `MAE $\\psi$`         = pm(mae_psi,   d = digits_small),
      `Cov $\\eta$`         = sprintf("%.2f", mean(cov_eta,   na.rm = TRUE)),
      `Cov $\\kappa$`       = sprintf("%.2f", mean(cov_kappa, na.rm = TRUE)),
      `Cov $\\psi$`         = sprintf("%.2f", mean(cov_psi,   na.rm = TRUE)),
      .groups = "drop"
    ) |>
    arrange(K, fit_model)
  
  tab_kbl <- kbl_ltx(tab,
                     caption = glue("{gen_lab}: Partition and parameter recovery (ARI, MAE, Coverage)."),
                     label   = glue("tab:{tolower(gen)}-main-partition-params"),
                     align   = "l" %+% paste(rep("r", ncol(tab)-1), collapse = ""))
  
  save_tab(tab_kbl, file.path(out_dir_main, glue("{gen}_main_table1_partition_params.tex")))
}

make_table2 <- function(df_gen, gen) {
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  
  deltas <- df_gen |>
    group_by(K, rep) |>
    mutate(elpd_best_rep = max(elpd_loo, na.rm = TRUE),
           dELPD = elpd_loo - elpd_best_rep) |>
    ungroup()
  
  tab <- deltas |>
    group_by(K, fit_model) |>
    summarise(
      `ELPD\\_LOO mean (sd)` = pm(elpd_loo, d = digits_mid),
      `$\\Delta$ELPD vs best (sd)` = pm(dELPD, d = digits_mid),
      `$p_{\\text{loo}}$ mean` = sprintf("%.2f", mean(p_loo, na.rm = TRUE)),
      `PSIS tail share` = sprintf("%.2f", mean(psis_tail_share, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    arrange(K, fit_model)
  
  tab_kbl <- kbl_ltx(tab,
                     caption = glue("{gen_lab}: PSIS–LOO comparison. ΔELPD computed per replicate within K; PSIS tail share = fraction with Pareto-\\(k>0.7\\)."),
                     label   = glue("tab:{tolower(gen)}-main-psis-loo"),
                     align   = "l" %+% paste(rep("r", ncol(tab)-1), collapse = ""))
  
  save_tab(tab_kbl, file.path(out_dir_main, glue("{gen}_main_table2_psis_loo.tex")))
}

make_table3 <- function(df_gen, gen) {
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  
  cyc_has_diag <- has_cols(df_gen, c("cyc_block_ess","cyc_block_mcse"))
  
  tab <- df_gen |>
    mutate(
      thetaW_flag = diag_flag(thetaW_block_ess, thetaW_block_mcse,
                              thetaW_block_lo, thetaW_block_hi),
      thetaS_flag = diag_flag(thetaS_block_ess, thetaS_block_mcse,
                              thetaS_block_lo, thetaS_block_hi),
      cyc_flag = if (cyc_has_diag)
        diag_flag(cyc_block_ess, cyc_block_mcse, cyc_block_lo, cyc_block_hi) else "",
      thetaW_disp = paste0(ci_str(thetaW_block_mean, thetaW_block_lo, thetaW_block_hi, d = digits_small), thetaW_flag),
      thetaS_disp = paste0(ci_str(thetaS_block_mean, thetaS_block_lo, thetaS_block_hi, d = digits_small), thetaS_flag),
      cyc_disp    = paste0(ci_str(cyc_block_mean,    cyc_block_lo,    cyc_block_hi,    d = digits_small), cyc_flag),
      sw_gap = thetaS_block_mean - thetaW_block_mean
    ) |>
    group_by(K, fit_model) |>
    summarise(
      `WST (block)` = {
        vals <- unique(thetaW_disp); if (length(vals)==0) "" else vals[order(vals)][ceiling(length(vals)/2)]
      },
      `SST (block)` = {
        vals <- unique(thetaS_disp); if (length(vals)==0) "" else vals[order(vals)][ceiling(length(vals)/2)]
      },
      `Cycles (block)` = {
        vals <- unique(cyc_disp); if (length(vals)==0) "" else vals[order(vals)][ceiling(length(vals)/2)]
      },
      `Premise size (block)` = sprintf("%.1f", mean(prem_block_avg, na.rm = TRUE)),
      `S--W gap` = sprintf(paste0("%.", digits_small, "f"), mean(sw_gap, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    arrange(K, fit_model)
  
  cap <- glue("{gen_lab}: Order recovery at block level. Entries show posterior mean with 95\\% CI; ",
              "\\(\\dagger\\) flags small ESS or large MCSE (>25\\% of CI width). Higher WST/SST and lower cycles indicate stronger order; S--W gap = SST−WST.")
  
  tab_kbl <- kbl_ltx(tab,
                     caption = cap,
                     label   = glue("tab:{tolower(gen)}-main-order-block"),
                     align   = "l" %+% paste(rep("r", ncol(tab)-1), collapse = ""))
  
  save_tab(tab_kbl, file.path(out_dir_main, glue("{gen}_main_table3_order_block.tex")))
}

# ---------- APPENDIX TABLES ----------------------------------

# A) Full VI tables (mean (sd)); difficulty stratified if available
make_appx_vi <- function(df_gen, gen) {
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  by_cols <- c("K","fit_model", intersect("difficulty", names(df_gen)))
  tab <- df_gen |>
    group_by(across(all_of(by_cols))) |>
    summarise(
      `NVI mean (sd)` = pm(vi_norm, d = digits_small),
      .groups = "drop"
    ) |>
    arrange(across(all_of(by_cols)))
  
  cap <- glue("{gen_lab}: Normalised Variation of Information (NVI = VI$/\\log_2 K^\\star$) summary.",
              if ("difficulty" %in% names(df_gen)) " Stratified by difficulty." else "")
  lab <- glue("tab:{tolower(gen)}-appx-vi")
  file <- file.path(out_dir_appx, glue("{gen}_appx_A_vi.tex"))
  tab_kbl <- kbl_ltx(tab, caption = cap, label = lab,
                     align = "l" %+% paste(rep("r", ncol(tab)-1), collapse = ""))
  save_tab(tab_kbl, file)
}

# B) Detailed parameter tables: MAE/Bias/Coverage with uncertainty
make_appx_params <- function(df_gen, gen) {
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  by_cols <- c("K","fit_model", intersect("difficulty", names(df_gen)))
  
  summarise_param <- function(x) pm(x, d = digits_small)
  tab <- df_gen |>
    group_by(across(all_of(by_cols))) |>
    summarise(
      `MAE $\\eta$`   = summarise_param(mae_eta),
      `Bias $\\eta$`  = summarise_param(bias_eta),
      `Cov $\\eta$`   = summarise_param(cov_eta),
      `MAE $\\kappa$` = summarise_param(mae_kappa),
      `Bias $\\kappa$`= summarise_param(bias_kappa),
      `Cov $\\kappa$` = summarise_param(cov_kappa),
      `MAE $\\psi$`   = summarise_param(mae_psi),
      `Bias $\\psi$`  = summarise_param(bias_psi),
      `Cov $\\psi$`   = summarise_param(cov_psi),
      .groups = "drop"
    ) |>
    arrange(across(all_of(by_cols)))
  
  cap <- glue("{gen_lab}: Parameter recovery details (mean (sd) across replicates).",
              if ("difficulty" %in% names(df_gen)) " Stratified by difficulty." else "")
  lab <- glue("tab:{tolower(gen)}-appx-params")
  file <- file.path(out_dir_appx, glue("{gen}_appx_B_params.tex"))
  tab_kbl <- kbl_ltx(tab, caption = cap, label = lab,
                     align = "l" %+% paste(rep("r", ncol(tab)-1), collapse = ""))
  save_tab(tab_kbl, file)
}

# C) Winners (ARI and LOOIC), tie-aware (fractional credit)
make_appx_winners <- function(df_gen, gen) {
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  
  # ARI winners: handle ties by splitting credit equally
  ari_frac <- df_gen |>
    group_by(K, rep) |>
    mutate(ari_max = max(ari, na.rm = TRUE),
           is_win = abs(ari - ari_max) < 1e-12) |>
    mutate(split = ifelse(is_win, 1/sum(is_win), 0)) |>
    ungroup() |>
    group_by(K, fit_model) |>
    summarise(frac = sum(split, na.rm = TRUE),
              n_rep = n_distinct(rep),
              prop = frac / n_rep,
              .groups = "drop") |>
    arrange(K, fit_model)
  
  tab_ari <- ari_frac |>
    select(K, fit_model,
           `winners (frac)` = frac,
           `proportion` = prop)
  
  tab_ari_kbl <- kbl_ltx(tab_ari,
                         caption = glue("{gen_lab}: ARI winners by fitted model (tie-aware fractional credit)."),
                         label   = glue("tab:{tolower(gen)}-appx-ari-winners"),
                         align   = "lrrrr")
  
  save_tab(tab_ari_kbl, file.path(out_dir_appx, glue("{gen}_appx_C_ari_winners.tex")))
  
  # LOOIC winners (lower is better): tie-aware
  loo_frac <- df_gen |>
    group_by(K, rep) |>
    mutate(loo_min = min(looic, na.rm = TRUE),
           is_win = abs(looic - loo_min) < 1e-12) |>
    mutate(split = ifelse(is_win, 1/sum(is_win), 0)) |>
    ungroup() |>
    group_by(K, fit_model) |>
    summarise(frac = sum(split, na.rm = TRUE),
              n_rep = n_distinct(rep),
              prop = frac / n_rep,
              .groups = "drop") |>
    arrange(K, fit_model)
  
  tab_loo <- loo_frac |>
    select(K, fit_model,
           `winners (frac)` = frac,
           `proportion` = prop)
  
  tab_loo_kbl <- kbl_ltx(tab_loo,
                         caption = glue("{gen_lab}: LOOIC winners by fitted model (tie-aware fractional credit)."),
                         label   = glue("tab:{tolower(gen)}-appx-loo-winners"),
                         align   = "lrrrr")
  
  save_tab(tab_loo_kbl, file.path(out_dir_appx, glue("{gen}_appx_C_loo_winners.tex")))
}

# D) PSIS diagnostics detail
make_appx_psis <- function(df_gen, gen) {
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  
  # Build what we can from available columns
  cols_present <- names(df_gen)
  tab <- df_gen |>
    group_by(K, fit_model) |>
    summarise(
      `ELPD\\_LOO mean (sd)` = pm(elpd_loo, d = digits_mid),
      `LOOIC mean (sd)`      = pm(looic, d = digits_mid),
      `$p_{\\text{loo}}$ mean (sd)` = pm(p_loo, d = digits_mid),
      `max Pareto-k (mean)`  = if ("pareto_k_max" %in% cols_present) sprintf("%.2f", mean(pareto_k_max, na.rm = TRUE)) else NA_character_,
      `PSIS tail share (k>0.7)` = if ("pareto_k_share" %in% cols_present) sprintf("%.2f", mean(pareto_k_share, na.rm = TRUE)) else sprintf("%.2f", mean(as.numeric(pareto_k_max > 0.7), na.rm = TRUE)),
      .groups = "drop"
    ) |>
    arrange(K, fit_model)
  
  cap <- glue("{gen_lab}: Detailed PSIS–LOO diagnostics. Tail share = fraction of points with Pareto-\\(k>0.7\\).")
  lab <- glue("tab:{tolower(gen)}-appx-psis")
  file <- file.path(out_dir_appx, glue("{gen}_appx_D_psis_detail.tex"))
  tab_kbl <- kbl_ltx(tab, caption = cap, label = lab,
                     align = "l" %+% paste(rep("r", ncol(tab)-1), collapse = ""))
  save_tab(tab_kbl, file)
}

# E) Item-level order metrics
make_appx_item_order <- function(df_gen, gen) {
  need <- c("thetaW_item_mean","thetaW_item_lo","thetaW_item_hi","thetaW_item_ess","thetaW_item_mcse",
            "thetaS_item_mean","thetaS_item_lo","thetaS_item_hi","thetaS_item_ess","thetaS_item_mcse",
            "prem_item_avg")
  if (!has_cols(df_gen, need)) {
    message("Skipping item-level order (missing columns) for ", gen)
    return(invisible(NULL))
  }
  
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  
  tab <- df_gen |>
    mutate(
      w_flag = diag_flag(thetaW_item_ess, thetaW_item_mcse, thetaW_item_lo, thetaW_item_hi),
      s_flag = diag_flag(thetaS_item_ess, thetaS_item_mcse, thetaS_item_lo, thetaS_item_hi),
      w_disp = paste0(ci_str(thetaW_item_mean, thetaW_item_lo, thetaW_item_hi, d = digits_small), w_flag),
      s_disp = paste0(ci_str(thetaS_item_mean, thetaS_item_lo, thetaS_item_hi, d = digits_small), s_flag)
    ) |>
    group_by(K, fit_model) |>
    summarise(
      `WST (item)` = {v <- unique(w_disp); if (length(v)==0) "" else v[order(v)][ceiling(length(v)/2)]},
      `SST (item)` = {v <- unique(s_disp); if (length(v)==0) "" else v[order(v)][ceiling(length(v)/2)]},
      `Premise size (item)` = sprintf("%.1f", mean(prem_item_avg, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    arrange(K, fit_model)
  
  cap <- glue("{gen_lab}: Item-level order metrics. Entries are posterior mean with 95\\% CI; \\(\\dagger\\) flags small ESS or large MCSE.")
  lab <- glue("tab:{tolower(gen)}-appx-item-order")
  file <- file.path(out_dir_appx, glue("{gen}_appx_E_item_order.tex"))
  tab_kbl <- kbl_ltx(tab, caption = cap, label = lab,
                     align = "l" %+% paste(rep("r", ncol(tab)-1), collapse = ""))
  save_tab(tab_kbl, file)
}

# F) MCMC quality tables (ESS/MCSE summaries for order metrics)
make_appx_mcmc_quality <- function(df_gen, gen) {
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  
  # Collect available metric groups
  blocks <- list(
    W_block = c("thetaW_block_ess","thetaW_block_mcse","thetaW_block_lo","thetaW_block_hi"),
    S_block = c("thetaS_block_ess","thetaS_block_mcse","thetaS_block_lo","thetaS_block_hi")
  )
  if (has_cols(df_gen, c("cyc_block_lo","cyc_block_hi"))) {
    if (has_cols(df_gen, c("cyc_block_ess","cyc_block_mcse")))
      blocks$C_block <- c("cyc_block_ess","cyc_block_mcse","cyc_block_lo","cyc_block_hi")
  }
  if (has_cols(df_gen, c("thetaW_item_ess","thetaW_item_mcse","thetaW_item_lo","thetaW_item_hi")))
    blocks$W_item <- c("thetaW_item_ess","thetaW_item_mcse","thetaW_item_lo","thetaW_item_hi")
  if (has_cols(df_gen, c("thetaS_item_ess","thetaS_item_mcse","thetaS_item_lo","thetaS_item_hi")))
    blocks$S_item <- c("thetaS_item_ess","thetaS_item_mcse","thetaS_item_lo","thetaS_item_hi")
  
  if (length(blocks) == 0) {
    message("Skipping MCMC quality (no ESS/MCSE columns) for ", gen)
    return(invisible(NULL))
  }
  
  make_summary <- function(ess_col, mcse_col, lo_col, hi_col) {
    df_gen |>
      transmute(K, fit_model,
                ess = .data[[ess_col]],
                mcse = .data[[mcse_col]],
                width = .data[[hi_col]] - .data[[lo_col]],
                mcse_rel = ifelse(width > 0, mcse/width, NA_real_)) |>
      group_by(K, fit_model) |>
      summarise(
        `ESS (median)` = sprintf("%.0f", median(ess, na.rm = TRUE)),
        `ESS (min)`    = sprintf("%.0f", suppressWarnings(min(ess, na.rm = TRUE))),
        `MCSE (median)`= sprintf("%.4f", median(mcse, na.rm = TRUE)),
        `MCSE/width (median)` = sprintf("%.2f", median(mcse_rel, na.rm = TRUE)),
        .groups = "drop"
      ) |>
      arrange(K, fit_model)
  }
  
  tabs <- imap(blocks, function(cols, name) {
    make_summary(cols[1], cols[2], cols[3], cols[4]) |>
      mutate(Metric = name) |>
      relocate(Metric, .before = K)
  })
  
  tab <- bind_rows(tabs)
  
  cap <- glue("{gen_lab}: MCMC quality summaries for order metrics (ESS and MCSE). MCSE/width is MCSE scaled by CI width.")
  lab <- glue("tab:{tolower(gen)}-appx-mcmc-quality")
  file <- file.path(out_dir_appx, glue("{gen}_appx_F_mcmc_quality.tex"))
  tab_kbl <- kbl_ltx(tab, caption = cap, label = lab,
                     align = "l" %+% paste(rep("r", ncol(tab)-1), collapse = ""))
  save_tab(tab_kbl, file)
}

# G) Difficulty-stratified versions of Main Tables 1–3 (if difficulty present)
make_appx_difficulty_versions <- function(df_gen, gen) {
  if (!("difficulty" %in% names(df_gen))) {
    message("No difficulty column; skipping difficulty-stratified appendix for ", gen)
    return(invisible(NULL))
  }
  
  gen_lab <- if (gen=="WST") "WST-generated data" else "SST-generated data"
  
  # Table 1 stratified
  tab1 <- df_gen |>
    group_by(difficulty, K, fit_model) |>
    summarise(
      `ARI mean (sd)`  = pm(ari, d = digits_small),
      `Pr(ARI ≥ 0.95)` = pct(mean(ari >= 0.95, na.rm = TRUE), d = 1),
      `MAE $\\eta$`    = pm(mae_eta,   d = digits_small),
      `MAE $\\kappa$`  = pm(mae_kappa, d = digits_small),
      `MAE $\\psi$`    = pm(mae_psi,   d = digits_small),
      .groups = "drop"
    ) |>
    arrange(difficulty, K, fit_model)
  
  tab1_kbl <- kbl_ltx(tab1,
                      caption = glue("{gen_lab}: Partition and parameter summaries by difficulty."),
                      label   = glue("tab:{tolower(gen)}-appx-diff-partition"),
                      align   = "l" %+% paste(rep("r", ncol(tab1)-1), collapse = ""))
  save_tab(tab1_kbl, file.path(out_dir_appx, glue("{gen}_appx_G_diff_table1.tex")))
  
  # Table 2 stratified
  deltas <- df_gen |>
    group_by(difficulty, K, rep) |>
    mutate(elpd_best_rep = max(elpd_loo, na.rm = TRUE),
           dELPD = elpd_loo - elpd_best_rep) |>
    ungroup()
  
  tab2 <- deltas |>
    group_by(difficulty, K, fit_model) |>
    summarise(
      `ELPD\\_LOO mean (sd)` = pm(elpd_loo, d = digits_mid),
      `$\\Delta$ELPD vs best (sd)` = pm(dELPD, d = digits_mid),
      `$p_{\\text{loo}}$ mean` = sprintf("%.2f", mean(p_loo, na.rm = TRUE)),
      `PSIS tail share` = sprintf("%.2f", mean(psis_tail_share, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    arrange(difficulty, K, fit_model)
  
  tab2_kbl <- kbl_ltx(tab2,
                      caption = glue("{gen_lab}: PSIS–LOO by difficulty."),
                      label   = glue("tab:{tolower(gen)}-appx-diff-psis-loo"),
                      align   = "l" %+% paste(rep("r", ncol(tab2)-1), collapse = ""))
  save_tab(tab2_kbl, file.path(out_dir_appx, glue("{gen}_appx_G_diff_table2.tex")))
  
  # Table 3 stratified (block order)
  cyc_has_diag <- has_cols(df_gen, c("cyc_block_ess","cyc_block_mcse"))
  
  tab3 <- df_gen |>
    mutate(
      thetaW_flag = diag_flag(thetaW_block_ess, thetaW_block_mcse, thetaW_block_lo, thetaW_block_hi),
      thetaS_flag = diag_flag(thetaS_block_ess, thetaS_block_mcse, thetaS_block_lo, thetaS_block_hi),
      cyc_flag = if (cyc_has_diag) diag_flag(cyc_block_ess, cyc_block_mcse, cyc_block_lo, cyc_block_hi) else "",
      thetaW_disp = paste0(ci_str(thetaW_block_mean, thetaW_block_lo, thetaW_block_hi, d = digits_small), thetaW_flag),
      thetaS_disp = paste0(ci_str(thetaS_block_mean, thetaS_block_lo, thetaS_block_hi, d = digits_small), thetaS_flag),
      cyc_disp    = paste0(ci_str(cyc_block_mean,    cyc_block_lo,    cyc_block_hi,    d = digits_small), cyc_flag)
    ) |>
    group_by(difficulty, K, fit_model) |>
    summarise(
      `WST (block)` = {v <- unique(thetaW_disp); if (length(v)==0) "" else v[order(v)][ceiling(length(v)/2)]},
      `SST (block)` = {v <- unique(thetaS_disp); if (length(v)==0) "" else v[order(v)][ceiling(length(v)/2)]},
      `Cycles (block)` = {v <- unique(cyc_disp); if (length(v)==0) "" else v[order(v)][ceiling(length(v)/2)]},
      `Premise size (block)` = sprintf("%.1f", mean(prem_block_avg, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    arrange(difficulty, K, fit_model)
  
  tab3_kbl <- kbl_ltx(tab3,
                      caption = glue("{gen_lab}: Block-level order by difficulty (\\(\\dagger\\) flags small ESS/large MCSE)."),
                      label   = glue("tab:{tolower(gen)}-appx-diff-order"),
                      align   = "l" %+% paste(rep("r", ncol(tab3)-1), collapse = ""))
  save_tab(tab3_kbl, file.path(out_dir_appx, glue("{gen}_appx_G_diff_table3.tex")))
}

# H) Compact headline table (K recovery + ELPD choice + empirical violations)
make_compact_headline <- function(compact_df) {
  need <- c("gen_model", "fit_model", "K_true", "K_hat", "kappa_mean",
            "psi_mean", "vi", "violation_rate_mean", "violation_rate_lo.2.5.", "violation_rate_hi.97.5.")
  if (!has_cols(compact_df, need)) {
    missing <- setdiff(need, names(compact_df))
    stop("Compact headline table skipped: missing columns: ", paste(missing, collapse = ", "))
  }

  compact_df <- compact_df |>
    mutate(
      fit_model_clean = case_when(
        fit_model %in% c("DC-SBM", "DCSBM") ~ "DC",
        TRUE ~ as.character(fit_model)
      ),
      fit_model_clean = factor(fit_model_clean, levels = c("WST", "SST", "DC"))
    ) |>
    filter(kappa_mean %in% c(0.75, 1.5, 3),
           psi_mean %in% c(0.2, 1.3))

  # Requested model-wise summaries
  model_summary <- compact_df |>
    group_by(gen_model, kappa_mean, psi_mean, fit_model_clean) |>
    summarise(
      mean_vi = mean(vi / log2(K_true), na.rm = TRUE),
      pr_KeqKstar = mean(K_hat == K_true, na.rm = TRUE),
      viol_partition_mean = mean(violation_rate_mean, na.rm = TRUE),
      viol_partition_lo = mean(`violation_rate_lo.2.5.`, na.rm = TRUE),
      viol_partition_hi = mean(`violation_rate_hi.97.5.`, na.rm = TRUE),
      .groups = "drop"
    )

  # Save long-format CSV for full transparency
  long_file <- file.path(out_dir_main, "sim_headline_compact_long.csv")
  readr::write_csv(model_summary, long_file)
  message("Saved: ", long_file)

  # Wide multicolumn table for manuscript
  fmt_num_bold <- function(x, best, d = 2) {
    txt <- sprintf(paste0("%.", d, "f"), x)
    out <- ifelse(abs(x - best) < 1e-12, paste0("\\textbf{", txt, "}"), txt)
    ifelse(is.na(x) | is.na(best), "--", out)
  }

  fmt_int_bold <- function(m, lo, hi, best) {
    txt <- sprintf("%.3f$_{[%.3f,%.3f]}$", m, lo, hi)
    out <- ifelse(abs(m - best) < 1e-12, paste0("\\textbf{", txt, "}"), txt)
    ifelse(is.na(m) | is.na(lo) | is.na(hi) | is.na(best), "--", out)
  }

  compact_wide <- model_summary |>
    mutate(
      kappa = sprintf("%.2f", kappa_mean),
      psi = sprintf("%.1f", psi_mean)
    ) |>
    select(
      `Gen.` = gen_model,
      kappa,
      psi,
      fit_model_clean,
      mean_vi,
      pr_KeqKstar,
      viol_partition_mean,
      viol_partition_lo,
      viol_partition_hi
    ) |>
    pivot_wider(
      names_from = fit_model_clean,
      values_from = c(mean_vi, pr_KeqKstar, viol_partition_mean, viol_partition_lo, viol_partition_hi),
      names_glue = "{.value}_{fit_model_clean}"
    )

  needed_cols <- c(
    "mean_vi_WST", "mean_vi_SST", "mean_vi_DC",
    "pr_KeqKstar_WST", "pr_KeqKstar_SST", "pr_KeqKstar_DC",
    "viol_partition_mean_WST", "viol_partition_mean_SST", "viol_partition_mean_DC",
    "viol_partition_lo_WST", "viol_partition_lo_SST", "viol_partition_lo_DC",
    "viol_partition_hi_WST", "viol_partition_hi_SST", "viol_partition_hi_DC"
  )
  for (nm in needed_cols) {
    if (!nm %in% names(compact_wide)) compact_wide[[nm]] <- NA_real_
  }

  compact_tab <- compact_wide |>
    arrange(`Gen.`, as.numeric(kappa), as.numeric(psi)) |>
    mutate(
      best_vi = pmin(mean_vi_WST, mean_vi_SST, mean_vi_DC, na.rm = TRUE),
      best_pr = pmax(pr_KeqKstar_WST, pr_KeqKstar_SST, pr_KeqKstar_DC, na.rm = TRUE),
      best_part = pmin(viol_partition_mean_WST, viol_partition_mean_SST, viol_partition_mean_DC, na.rm = TRUE),

      vi_fmt_WST = fmt_num_bold(mean_vi_WST, best_vi, d = 3),
      vi_fmt_SST = fmt_num_bold(mean_vi_SST, best_vi, d = 3),
      vi_fmt_DC = fmt_num_bold(mean_vi_DC, best_vi, d = 3),

      pr_fmt_WST = fmt_num_bold(pr_KeqKstar_WST, best_pr, d = 2),
      pr_fmt_SST = fmt_num_bold(pr_KeqKstar_SST, best_pr, d = 2),
      pr_fmt_DC = fmt_num_bold(pr_KeqKstar_DC, best_pr, d = 2),

      part_fmt_WST = fmt_int_bold(viol_partition_mean_WST, viol_partition_lo_WST, viol_partition_hi_WST, best_part),
      part_fmt_SST = fmt_int_bold(viol_partition_mean_SST, viol_partition_lo_SST, viol_partition_hi_SST, best_part),
      part_fmt_DC = fmt_int_bold(viol_partition_mean_DC, viol_partition_lo_DC, viol_partition_hi_DC, best_part)
    ) |>
    select(
      `Gen.`, kappa, psi,
      vi_fmt_WST, vi_fmt_SST, vi_fmt_DC,
      pr_fmt_WST, pr_fmt_SST, pr_fmt_DC,
      part_fmt_WST, part_fmt_SST, part_fmt_DC
    )

  wide_file <- file.path(out_dir_main, "sim_headline_compact.csv")
  readr::write_csv(compact_tab, wide_file)
  message("Saved: ", wide_file)

  cap <- "Compact simulation headline summary stratified by generating model, density level (kappa), and hierarchy level (psi). Columns are grouped by metric with fitted models side-by-side (WST, SST, DC). NVI = VI$/\\log_2(K^\\star)$. Bold marks the best model per row for each metric (minimum NVI, maximum Pr(K=K*), minimum partition violation mean)."

  tab_kbl <- kbl(
    compact_tab,
    format = "latex",
    booktabs = TRUE,
    caption = glue("{cap}\\\\\\label{{tab:sim-headline-compact}}"),
    align = "l" %+% paste(rep("r", ncol(compact_tab) - 1), collapse = ""),
    col.names = c("Gen.", "kappa", "psi",
                  rep(c("WST", "SST", "DC"), 3)),
    escape = FALSE,
    linesep = ""
  ) |>
    add_header_above(c(" " = 3,
                       "Mean NVI" = 3,
                       "Pr(K=K*)" = 3,
                       "Partition violation rate" = 3), escape = FALSE) |>
    kable_styling(full_width = FALSE, position = "center",
                  latex_options = c("hold_position"))

  tex_file <- file.path(out_dir_main, "sim_headline_compact.tex")
  save_tab(tab_kbl, tex_file)

  message("\nCompact table preview:")
  print(compact_tab)

  invisible(list(long = model_summary, table = compact_tab))
}

# ---------- RUN (WST and SST) -------------------------------------------
legacy_required <- c(
  "mae_eta", "mae_kappa", "mae_psi",
  "cov_eta", "cov_kappa", "cov_psi",
  "elpd_loo", "p_loo", "looic",
  "thetaW_block_mean", "thetaS_block_mean", "cyc_block_mean",
  "thetaW_block_lo", "thetaW_block_hi", "thetaW_block_ess", "thetaW_block_mcse",
  "thetaS_block_lo", "thetaS_block_hi", "thetaS_block_ess", "thetaS_block_mcse",
  "cyc_block_lo", "cyc_block_hi", "prem_block_avg"
)

missing_legacy <- setdiff(legacy_required, names(raw))

if (length(missing_legacy) == 0) {
  for (gen in c("WST","SST")) {
    df_gen <- raw |> filter(gen_model == gen)

    # Main-text
    make_table1(df_gen, gen)
    make_table2(df_gen, gen)
    make_table3(df_gen, gen)

    # Appendix
    make_appx_vi(df_gen, gen)                 # A
    make_appx_params(df_gen, gen)             # B
    make_appx_winners(df_gen, gen)            # C
    make_appx_psis(df_gen, gen)               # D
    make_appx_item_order(df_gen, gen)         # E (skips if missing)
    make_appx_mcmc_quality(df_gen, gen)       # F (skips if missing)
    make_appx_difficulty_versions(df_gen, gen)# G (skips if missing)
  }
} else {
  message("Skipping legacy main/appendix tables for current results file: missing columns -> ",
          paste(missing_legacy, collapse = ", "))
}

# Compact headline summary from the current results file (or explicit override)
compact_source_path <- if (nzchar(Sys.getenv("COMPACT_RESULTS_PATH"))) compact_results_path else results_path
compact_raw <- readr::read_csv(compact_source_path, show_col_types = FALSE)
compact_headline <- make_compact_headline(compact_raw)

message("\nDone. Include main tables with, e.g.:")
message("\\input{", file.path(out_dir_main, "WST_main_table1_partition_params.tex"), "}")
message("\\input{", file.path(out_dir_main, "WST_main_table2_psis_loo.tex"), "}")
message("\\input{", file.path(out_dir_main, "WST_main_table3_order_block.tex"), "}")
message("\\input{", file.path(out_dir_main, "SST_main_table1_partition_params.tex"), "}")
message("\\input{", file.path(out_dir_main, "SST_main_table2_psis_loo.tex"), "}")
message("\\input{", file.path(out_dir_main, "SST_main_table3_order_block.tex"), "}")

message("\nAppendix suggestions:")
message("\\input{", file.path(out_dir_appx, "WST_appx_A_vi.tex"), "}")
message("\\input{", file.path(out_dir_appx, "WST_appx_B_params.tex"), "}")
message("\\input{", file.path(out_dir_appx, "WST_appx_C_ari_winners.tex"), "}")
message("\\input{", file.path(out_dir_appx, "WST_appx_C_loo_winners.tex"), "}")
message("\\input{", file.path(out_dir_appx, "WST_appx_D_psis_detail.tex"), "}")
message("\\input{", file.path(out_dir_appx, "WST_appx_E_item_order.tex"), "}")
message("\\input{", file.path(out_dir_appx, "WST_appx_F_mcmc_quality.tex"), "}")
message("\\input{", file.path(out_dir_appx, "WST_appx_G_diff_table1.tex"), "}")
message("\\input{", file.path(out_dir_appx, "WST_appx_G_diff_table2.tex"), "}")
message("\\input{", file.path(out_dir_appx, "WST_appx_G_diff_table3.tex"), "}")

message("\\input{", file.path(out_dir_appx, "SST_appx_A_vi.tex"), "}")
message("\\input{", file.path(out_dir_appx, "SST_appx_B_params.tex"), "}")
message("\\input{", file.path(out_dir_appx, "SST_appx_C_ari_winners.tex"), "}")
message("\\input{", file.path(out_dir_appx, "SST_appx_C_loo_winners.tex"), "}")
message("\\input{", file.path(out_dir_appx, "SST_appx_D_psis_detail.tex"), "}")
message("\\input{", file.path(out_dir_appx, "SST_appx_E_item_order.tex"), "}")
message("\\input{", file.path(out_dir_appx, "SST_appx_F_mcmc_quality.tex"), "}")
message("\\input{", file.path(out_dir_appx, "SST_appx_G_diff_table1.tex"), "}")
message("\\input{", file.path(out_dir_appx, "SST_appx_G_diff_table2.tex"), "}")
message("\\input{", file.path(out_dir_appx, "SST_appx_G_diff_table3.tex"), "}")

message("\nCompact headline table:")
message("\\input{", file.path(out_dir_main, "sim_headline_compact.tex"), "}")
