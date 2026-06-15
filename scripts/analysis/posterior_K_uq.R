#!/usr/bin/env Rscript
############################################################
## Posterior K uncertainty quantification
##
## Given saved OSBM/DCSBM fit objects (from application.R or
## varK_app.R), produces:
##   1. Posterior PMF of K  (bar chart per dataset × model)
##   2. K trace plot        (line plot per dataset × model)
##   3. Summary table with posterior mean, median, mode,
##      95% HPD interval for K
##   4. Posterior similarity matrix (PSM) heatmaps
##
## Usage:
##   Rscript posterior_K_uq.R [fit_dir]
##
## fit_dir: directory containing *_fit.rds files (from
##          posterior_predictive_checks/<dataset>/<model>/
##          or figures/application_runs/<run_id>/)
##          Defaults to the latest run in figures/application_runs/
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(glue)
  library(fs)
})

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

# ============================================================
# 1. Locate fit files
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

# Scan run dir first (canonical, newer fits), ppc as fallback
find_fit_files <- function(base_dir = "output/application/ppc") {
  fits <- list()

  # Primary source: canonical run dir (flat layout)
  run_dir <- Sys.getenv("APP_RUN_DIR",
    unset = "./output/application/raw/application_run_20260411_163055")
  if (dir.exists(run_dir)) {
    rds_files <- list.files(run_dir, pattern = "_fit\\.rds$", full.names = TRUE)
    for (f in rds_files) {
      bn <- basename(f)
      m <- regmatches(bn, regexec("^(.+)_(WST|SST|DCSBM)_fit\\.rds$", bn, perl = TRUE))[[1]]
      if (length(m) < 3) next
      fits[[length(fits) + 1]] <- list(dataset = m[2], model = m[3], path = f)
    }
  }

  # Supplement with ppc dir fits for dataset/model combos not already found
  run_keys <- vapply(fits, function(x) paste0(x$dataset, "_", x$model), character(1))
  if (dir.exists(base_dir)) {
    for (ds_dir in list.dirs(base_dir, recursive = FALSE)) {
      dataset <- basename(ds_dir)
      for (model_dir in list.dirs(ds_dir, recursive = FALSE)) {
        model <- basename(model_dir)
        if (paste0(dataset, "_", model) %in% run_keys) next
        rds_files <- list.files(model_dir, pattern = "_fit\\.rds$", full.names = TRUE)
        for (f in rds_files) {
          fits[[length(fits) + 1]] <- list(
            dataset = dataset, model = model, path = f
          )
        }
      }
    }
  }

  fits
}

fit_meta <- find_fit_files()
if (length(fit_meta) == 0) {
  cat("No fit files found in output/application/ppc/ or latest run dir.\n")
  cat("Run application.R first to generate fits.\n")
  q(status = 0)
}

cat(sprintf("Found %d fit files:\n", length(fit_meta)))
for (fm in fit_meta) {
  cat(sprintf("  %s / %s : %s\n", fm$dataset, fm$model, basename(fm$path)))
}

# ============================================================
# 2. Extract K posterior from each fit
# ============================================================

extract_K_info <- function(fit_obj) {
  # Get z chain
  z_raw <- fit_obj$z
  if (is.null(z_raw)) z_raw <- fit_obj$out$z
  if (is.null(z_raw)) z_raw <- fit_obj$fit$z
  if (is.null(z_raw)) return(NULL)

  if (is.matrix(z_raw)) {
    z_chain <- z_raw
  } else if (is.list(z_raw)) {
    z_chain <- do.call(rbind, lapply(z_raw, as.integer))
  } else {
    return(NULL)
  }

  # K trace: number of distinct labels per draw
  K_trace <- apply(z_chain, 1, function(x) length(unique(x[!is.na(x)])))

  # Also check for pre-computed K_trace
  K_trace_stored <- fit_obj$K_trace
  if (is.null(K_trace_stored)) K_trace_stored <- fit_obj$out$K_trace
  if (!is.null(K_trace_stored) && length(K_trace_stored) > 0 &&
      all(K_trace_stored > 0, na.rm = TRUE)) {
    K_trace <- K_trace_stored[K_trace_stored > 0]
  }

  S <- length(K_trace)
  K_vals <- sort(unique(K_trace))
  K_pmf  <- table(K_trace) / S
  K_mean <- mean(K_trace)
  K_median <- median(K_trace)
  K_mode <- as.integer(names(which.max(table(K_trace))))
  K_sd   <- sd(K_trace)

  # 95% HPD-like interval (shortest interval containing 95% of mass)
  K_lo <- quantile(K_trace, probs = 0.025, type = 1)
  K_hi <- quantile(K_trace, probs = 0.975, type = 1)

  list(
    z_chain  = z_chain,
    K_trace  = K_trace,
    K_pmf    = K_pmf,
    K_mean   = K_mean,
    K_median = K_median,
    K_mode   = K_mode,
    K_sd     = K_sd,
    K_lo     = as.integer(K_lo),
    K_hi     = as.integer(K_hi),
    S        = S
  )
}

# ============================================================
# 3. Process all fits
# ============================================================

all_K_data <- list()
all_K_summaries <- list()
all_K_traces_df <- list()

for (fm in fit_meta) {
  cat(sprintf("Processing %s / %s ... ", fm$dataset, fm$model))
  fit_obj <- tryCatch(readRDS(fm$path), error = function(e) NULL)
  if (is.null(fit_obj)) { cat("FAILED to read\n"); next }

  ki <- extract_K_info(fit_obj)
  if (is.null(ki)) { cat("no z chain found\n"); next }

  cat(sprintf("S=%d, K_mode=%d, K_mean=%.1f [%d, %d]\n",
              ki$S, ki$K_mode, ki$K_mean, ki$K_lo, ki$K_hi))

  key <- paste0(fm$dataset, "_", fm$model)
  all_K_data[[key]] <- ki

  # Summary row
  all_K_summaries[[length(all_K_summaries) + 1]] <- data.frame(
    dataset  = fm$dataset,
    model    = fm$model,
    K_mode   = ki$K_mode,
    K_mean   = round(ki$K_mean, 2),
    K_median = ki$K_median,
    K_sd     = round(ki$K_sd, 2),
    K_lo_95  = ki$K_lo,
    K_hi_95  = ki$K_hi,
    S        = ki$S,
    stringsAsFactors = FALSE
  )

  # K trace data frame
  all_K_traces_df[[length(all_K_traces_df) + 1]] <- data.frame(
    dataset = fm$dataset,
    model   = fm$model,
    iter    = seq_along(ki$K_trace),
    K       = ki$K_trace,
    stringsAsFactors = FALSE
  )
}

if (length(all_K_summaries) == 0) {
  cat("No valid fits processed.\n")
  q(status = 0)
}

summary_df <- bind_rows(all_K_summaries)
traces_df  <- bind_rows(all_K_traces_df)

# ============================================================
# 4. Summary table
# ============================================================

cat("\n========== POSTERIOR K SUMMARY ==========\n")
print(as.data.frame(summary_df), row.names = FALSE)

out_dir <- "output/application/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(summary_df, file.path(out_dir, "K_posterior_summary.csv"), row.names = FALSE)

# ============================================================
# 5. Posterior PMF of K (bar charts)
# ============================================================

# Build PMF data
pmf_list <- list()
for (key in names(all_K_data)) {
  ki <- all_K_data[[key]]
  parts <- strsplit(key, "_")[[1]]
  model <- parts[length(parts)]
  dataset <- paste(parts[-length(parts)], collapse = "_")

  K_table <- as.data.frame(ki$K_pmf, stringsAsFactors = FALSE)
  names(K_table) <- c("K", "prob")
  K_table$K <- as.integer(as.character(K_table$K))
  K_table$dataset <- dataset
  K_table$model   <- model
  K_table$K_mode  <- ki$K_mode
  pmf_list[[key]] <- K_table
}
pmf_df <- bind_rows(pmf_list)

# Color palette for models
model_colors <- c("WST" = "#E69F00", "SST" = "#56B4E9",
                  "DCSBM" = "#009E73", "DC-SBM" = "#009E73")

p_pmf <- ggplot(pmf_df, aes(x = factor(K), y = prob, fill = model)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
  facet_wrap(~ dataset, scales = "free", ncol = 2) +
  scale_fill_manual(values = model_colors, name = "Model") +
  labs(x = "Number of blocks (K)",
       y = "Posterior probability",
       title = "Posterior distribution of K") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "K_posterior_pmf.pdf"), p_pmf,
       width = 10, height = 8)
ggsave(file.path(out_dir, "K_posterior_pmf.png"), p_pmf,
       width = 10, height = 8, dpi = 200)
cat(glue("Saved: {out_dir}/K_posterior_pmf.pdf\n"))

# ============================================================
# 6. K trace plots
# ============================================================

p_trace <- ggplot(traces_df, aes(x = iter, y = K, color = model)) +
  geom_line(alpha = 0.6, linewidth = 0.3) +
  facet_wrap(~ dataset, scales = "free_y", ncol = 2) +
  scale_color_manual(values = model_colors, name = "Model") +
  labs(x = "MCMC iteration (post burn-in, thinned)",
       y = "K (occupied blocks)",
       title = "Trace of K across MCMC iterations") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "K_trace.pdf"), p_trace,
       width = 10, height = 8)
ggsave(file.path(out_dir, "K_trace.png"), p_trace,
       width = 10, height = 8, dpi = 200)
cat(glue("Saved: {out_dir}/K_trace.pdf\n"))

# ============================================================
# 7. Per-dataset combined panels (PMF + trace side by side)
# ============================================================

datasets_present <- unique(summary_df$dataset)

for (ds in datasets_present) {
  ds_pmf <- pmf_df %>% filter(dataset == ds)
  ds_trace <- traces_df %>% filter(dataset == ds)
  ds_summ <- summary_df %>% filter(dataset == ds)

  if (nrow(ds_pmf) == 0 || nrow(ds_trace) == 0) next

  p1 <- ggplot(ds_pmf, aes(x = factor(K), y = prob, fill = model)) +
    geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
    scale_fill_manual(values = model_colors, name = "Model") +
    labs(x = "K", y = "P(K | data)", title = paste0(ds, ": Posterior of K")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  p2 <- ggplot(ds_trace, aes(x = iter, y = K, color = model)) +
    geom_line(alpha = 0.6, linewidth = 0.3) +
    scale_color_manual(values = model_colors, name = "Model") +
    labs(x = "Iteration", y = "K", title = paste0(ds, ": K trace")) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")

  p_combined <- p1 + p2 + plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")

  ds_out_dir <- file.path(out_dir, ds)
  dir.create(ds_out_dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(ds_out_dir, paste0(ds, "_K_uq_panel.pdf")), p_combined,
         width = 12, height = 5)
  ggsave(file.path(ds_out_dir, paste0(ds, "_K_uq_panel.png")), p_combined,
         width = 12, height = 5, dpi = 200)
  cat(glue("Saved: {ds_out_dir}/{ds}_K_uq_panel.pdf\n"))
}

# ============================================================
# 8. Posterior similarity matrix heatmaps (optional, if PSM
#    computation is feasible)
# ============================================================

# PSM plots are now generated by plotting_script.R (with block-coloured sidebars)
# so we skip the plain PSM generation here to avoid duplicates.

# ============================================================
# 9. LaTeX-ready summary table
# ============================================================

cat("\n=== LaTeX Summary Table ===\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%s & %s & %d & %.1f & [%d, %d] \\\\\n",
              r$dataset, r$model, r$K_mode, r$K_mean,
              r$K_lo_95, r$K_hi_95))
}

cat(sprintf("\nAll outputs saved to: %s/\n", out_dir))
cat("Done.\n")
