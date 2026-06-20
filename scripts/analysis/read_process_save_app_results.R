# applications_empirical_tables.R
# Produce LaTeX tables (empirical block + empirical item) with added columns:
#  - zeta_W (all), zeta_S (all), Prem. (avg), Cov. (avg)
# Bracketed 95% CI and per-cell dagger if ESS<200 or MCSE > 25% of CI width.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(glue)
  library(tidyr)
  library(purrr)
  library(knitr)
  library(kableExtra)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]),
                               winslash = "/", mustWork = TRUE)
  setwd(normalizePath(file.path(dirname(script_path), "../.."),
                      winslash = "/", mustWork = TRUE))
}
# ---------- I/O ----------
# Start with the well-established per-dataset summary files (richest columns)
in_paths <- c(
  citations   = "./output/application/raw/citations_data_applications_summary.csv",
  macaques    = "./output/application/raw/macaques_data_applications_summary.csv",
  mountain_goats = "./output/application/raw/mountain_goats_applications_summary.csv"
)
in_paths <- in_paths[file.exists(in_paths)]

# Auto-discover any additional per-dataset summary files in the raw dir
raw_dir <- "./output/application/raw"
extra_files <- list.files(raw_dir, pattern = "_applications_summary\\.csv$",
                           recursive = FALSE, full.names = TRUE)
extra_keys <- sub("_applications_summary\\.csv$", "",
                  basename(extra_files))
for (i in seq_along(extra_keys)) {
  if (!extra_keys[i] %in% names(in_paths)) {
    in_paths[extra_keys[i]] <- extra_files[i]
  }
}

# Also pull in datasets present in the hardcoded run dir that have no summary file yet
run_dir <- Sys.getenv("APP_RUN_DIR",
  unset = "./output/application/raw/application_run_20260411_163055")
if (dir.exists(run_dir)) {
  # Try old-format summary first
  run_sum_files <- list.files(run_dir,
                               pattern = "^applications_results_summary_.*\\.csv$",
                               full.names = TRUE)
  # Also try all_results.csv (new variable-K format)
  all_res_file <- file.path(run_dir, "all_results.csv")
  
  run_df <- NULL
  if (length(run_sum_files)) {
    run_df <- readr::read_csv(run_sum_files[1], show_col_types = FALSE)
  } else if (file.exists(all_res_file)) {
    run_df <- readr::read_csv(all_res_file, show_col_types = FALSE)
    # Harmonise column names for compatibility
    if ("K_hat" %in% names(run_df) && !"K" %in% names(run_df))
      run_df <- dplyr::rename(run_df, K = K_hat)
    if ("elpd_loo" %in% names(run_df) && !"looic" %in% names(run_df))
      run_df <- dplyr::mutate(run_df, looic = -2 * elpd_loo)
    if ("pk_max" %in% names(run_df) && !"pareto_k_max" %in% names(run_df))
      run_df <- dplyr::rename(run_df, pareto_k_max = pk_max)
  }
  
  if (!is.null(run_df)) {
    for (ds in unique(run_df$dataset)) {
      if (!ds %in% names(in_paths)) {
        # Write a per-dataset CSV so the downstream logic can read it uniformly
        ds_tmp <- file.path(raw_dir, paste0(ds, "_applications_summary.csv"))
        if (!file.exists(ds_tmp)) {
          readr::write_csv(run_df[run_df$dataset == ds, ], ds_tmp)
        }
        in_paths[ds] <- ds_tmp
      }
    }
  }
}

out_dir <- "./output/application/tables/"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- helpers ----------
fmt_cell <- function(m, lo, hi, ess, mcse, digits = 3) {
  # dagger rule: ESS < 200 OR MCSE > 25% of CI width
  width <- abs(hi - lo)
  bad <- (is.finite(ess) && ess < 200) ||
    (is.finite(width) && width > 0 && is.finite(mcse) && mcse > 0.25 * width)
  dag <- if (isTRUE(bad)) "\\textsuperscript{$\\dagger$}" else ""
  glue("{formatC(m, format='f', digits=digits)} [{formatC(lo, format='f', digits=digits)}--{formatC(hi, format='f', digits=digits)}]{dag}")
}

fmt_num <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

nice_fit <- function(x) recode(x,
                               "DCSBM" = "DC-SBM",
                               "DC-SBM" = "DC-SBM",
                               "SST" = "SST", "WST" = "WST",
                               .default = x
)

order_models <- function(x) factor(nice_fit(x), levels = c("SST","WST","DC-SBM"))

# safe column accessor (returns NA vector if missing)
col_or_na <- function(df, nm, n) if (nm %in% names(df)) df[[nm]] else rep(NA_real_, n)

# Build a single table (block or item) for a dataset df ----------------------
build_emp_table <- function(df, dataset_label, level = c("block","item")) {
  level <- match.arg(level)
  df <- df %>%
    group_by(fit_model) %>%
    slice_min(looic, with_ties = FALSE) %>%
    mutate(fit_model = order_models(fit_model)) %>%
    arrange(fit_model, K)
   
  
  n <- nrow(df)
  
  if (level == "block") {
    W_pref      <- "thetaW_block_emp_"
    S_pref      <- "thetaS_block_emp_"
    Wall_pref   <- "thetaW_block_emp_all_"
    Sall_pref   <- "thetaS_block_emp_all_"
    prem_col    <- "prem_block_emp_avg"
    cov_col     <- "coverage_block_emp_avg"
    cap <- glue("{dataset_label}: Empirical order at block level. Entries show posterior mean with 95\\% CI; $\\dagger$ flags small ESS ($<200$) or large MCSE ($>25\\%$ of CI width). “(all)” uses the unconditional success set (no premise); “Prem.” is the average number of premise-satisfying triples and “Cov.” the average premise coverage.")
    lab <- glue("tab:{str_replace_all(tolower(dataset_label),'[^a-z0-9]+','-')}-emp-block")
    zetaW_head <- "$\\zeta_{\\mathrm{W}}$ (block)"
    zetaS_head <- "$\\zeta_{\\mathrm{S}}$ (block)"
    zetaW_all_head <- "$\\zeta_{\\mathrm{W}}$ (block, all)"
    zetaS_all_head <- "$\\zeta_{\\mathrm{S}}$ (block, all)"
  } else {
    W_pref      <- "thetaW_item_emp_"
    S_pref      <- "thetaS_item_emp_"
    Wall_pref   <- "thetaW_item_emp_all_"
    Sall_pref   <- "thetaS_item_emp_all_"
    prem_col    <- "prem_item_emp_avg"
    cov_col     <- "coverage_item_emp_avg"
    cap <- glue("{dataset_label}: Empirical order at item level. Entries show posterior mean with 95\\% CI; $\\dagger$ flags small ESS ($<200$) or large MCSE ($>25\\%$ of CI width). “(all)” uses the unconditional success set (no premise); “Prem.” is the average number of premise-satisfying triples and “Cov.” the average premise coverage.")
    lab <- glue("tab:{str_replace_all(tolower(dataset_label),'[^a-z0-9]+','-')}-emp-item")
    zetaW_head <- "$\\zeta_{\\mathrm{W}}$ (item)"
    zetaS_head <- "$\\zeta_{\\mathrm{S}}$ (item)"
    zetaW_all_head <- "$\\zeta_{\\mathrm{W}}$ (item, all)"
    zetaS_all_head <- "$\\zeta_{\\mathrm{S}}$ (item, all)"
  }
  
  # pull vectors (or NAs) for safety
  Wm   <- col_or_na(df, paste0(W_pref,  "mean"), n)
  Wlo  <- col_or_na(df, paste0(W_pref,  "lo"),   n)
  Whi  <- col_or_na(df, paste0(W_pref,  "hi"),   n)
  Wess <- col_or_na(df, paste0(W_pref,  "ess"),  n)
  Wmc  <- col_or_na(df, paste0(W_pref,  "mcse"), n)
  
  Sm   <- col_or_na(df, paste0(S_pref,  "mean"), n)
  Slo  <- col_or_na(df, paste0(S_pref,  "lo"),   n)
  Shi  <- col_or_na(df, paste0(S_pref,  "hi"),   n)
  Sess <- col_or_na(df, paste0(S_pref,  "ess"),  n)
  Smc  <- col_or_na(df, paste0(S_pref,  "mcse"), n)
  
  WAm   <- col_or_na(df, paste0(Wall_pref,  "mean"), n)
  WAlo  <- col_or_na(df, paste0(Wall_pref,  "lo"),   n)
  WAhi  <- col_or_na(df, paste0(Wall_pref,  "hi"),   n)
  WAess <- col_or_na(df, paste0(Wall_pref,  "ess"),  n)
  WAmc  <- col_or_na(df, paste0(Wall_pref,  "mcse"), n)
  
  SAm   <- col_or_na(df, paste0(Sall_pref,  "mean"), n)
  SAlo  <- col_or_na(df, paste0(Sall_pref,  "lo"),   n)
  SAhi  <- col_or_na(df, paste0(Sall_pref,  "hi"),   n)
  SAess <- col_or_na(df, paste0(Sall_pref,  "ess"),  n)
  SAmc  <- col_or_na(df, paste0(Sall_pref,  "mcse"), n)
  
  Prem <- col_or_na(df, prem_col, n)
  Covg <- col_or_na(df, cov_col,  n)
  
  # Rowwise formatting so fmt_cell gets scalars
  tbl <- df %>%
    mutate(
      `K` = K,
      `fit model` = as.character(fit_model)
    ) %>%
    dplyr::select(`K`, `fit model`) %>%
    bind_cols(
      tibble(
        zetaW_fmt    = pmap_chr(list(Wm,  Wlo,  Whi,  Wess,  Wmc),  ~fmt_cell(..1,..2,..3,..4,..5)),
        zetaS_fmt    = pmap_chr(list(Sm,  Slo,  Shi,  Sess,  Smc),  ~fmt_cell(..1,..2,..3,..4,..5)),
        zetaW_allfmt = pmap_chr(list(WAm, WAlo, WAhi, WAess, WAmc), ~fmt_cell(..1,..2,..3,..4,..5)),
        zetaS_allfmt = pmap_chr(list(SAm, SAlo, SAhi, SAess, SAmc), ~fmt_cell(..1,..2,..3,..4,..5)),
        Prem_avg     = fmt_num(Prem, digits = 1),
        Cov_avg      = fmt_num(Covg, digits = 3)
      )
    ) %>%
    rename(
      !!zetaW_head     := zetaW_fmt,
      !!zetaS_head     := zetaS_fmt,
      !!zetaW_all_head := zetaW_allfmt,
      !!zetaS_all_head := zetaS_allfmt,
      `Prem. (avg)`    := Prem_avg,
      `Cov. (avg)`     := Cov_avg
    )
  
  # Build LaTeX table
  align <- c("l","l", rep("r", ncol(tbl) - 2))
  kt <- kable(
    tbl, format = "latex", booktabs = TRUE, align = align,
    caption = cap, escape = FALSE, label = lab, linesep = ""
  ) %>%
    kable_classic(full_width = FALSE, html_font = "Times") %>%
    add_header_above(c(" " = 2,
                       "Empirical (premise)" = 2,
                       "Empirical (all)" = 2,
                       "Sets" = 2)) %>%
    footnote(
      general = "Entries show posterior mean with 95% CI; \\dagger flags small ESS (<200) or large MCSE (>25% of CI width). “(all)” ignores the premise and scores unconditional success; “Prem.” is the average number of triples satisfying the premise; “Cov.” is the average premise coverage.",
      threeparttable = TRUE, escape = FALSE
    )
  
  list(table = kt, label = lab)
}

# ---------- read, build, write ----------
dfs <- lapply(names(in_paths), function(key) {
  df <- readr::read_csv(in_paths[[key]], show_col_types = FALSE)
  df$.__dataset__ <- key
  df
}) %>% bind_rows()

# Pretty labels: pre-defined for known datasets, auto-generated for new ones
key_to_label_base <- c(
  citations      = "Citations \u2014 Empirical",
  macaques       = "Macaques \u2014 Empirical",
  mountain_goats = "Mountain goats \u2014 Empirical",
  high_school    = "High school \u2014 Empirical",
  moreno_sheep   = "Sheep (Moreno) \u2014 Empirical",
  strauss_2019b  = "Hyenas (Strauss 2019b) \u2014 Empirical"
)
key_to_label <- key_to_label_base
for (k in names(in_paths)) {
  if (!k %in% names(key_to_label)) {
    key_to_label[k] <- paste0(tools::toTitleCase(gsub("_", " ", k)), " \u2014 Empirical")
  }
}

walk(names(in_paths), function(key) {
  dfk <- dfs %>% filter(.__dataset__ == key)
  if (nrow(dfk) == 0) {
    message(glue("Skipping {key}: no rows found"))
    return(invisible(NULL))
  }
  
  # BLOCK level
  res_block <- tryCatch(
    build_emp_table(dfk, dataset_label = key_to_label[[key]], level = "block"),
    error = function(e) { message("Block table error for ", key, ": ", conditionMessage(e)); NULL }
  )
  if (!is.null(res_block)) {
    block_path <- file.path(out_dir, glue("{key}_emp_block.tex"))
    cat(res_block$table, file = block_path)
    message(glue("Wrote: {block_path}  (\\ref{{{res_block$label}}})"))
  }
  
  # ITEM level
  res_item <- tryCatch(
    build_emp_table(dfk, dataset_label = key_to_label[[key]], level = "item"),
    error = function(e) { message("Item table error for ", key, ": ", conditionMessage(e)); NULL }
  )
  if (!is.null(res_item)) {
    item_path <- file.path(out_dir, glue("{key}_emp_item.tex"))
    cat(res_item$table, file = item_path)
    message(glue("Wrote: {item_path}   (\\ref{{{res_item$label}}})"))
  }
})
