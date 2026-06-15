#!/usr/bin/env Rscript
# =============================================================================
# build_bt_delta_summary.R
#
# Computes Bradley--Terry additivity diagnostics from WST posterior psi draws
# across the six application datasets. For each dataset we compute
#   Delta_{k l m} = psi_{k m} - psi_{k l} - psi_{l m},   k < l < m
# and summarise the draw-level mean Delta with a 95% interval.
#
# Outputs:
#   <tables_dir>/bt_delta_summary.csv
#   <tables_dir>/bt_delta_summary_table.tex
#   <figures_dir>/bt_delta_wst_applications.pdf
#   <figures_dir>/bt_delta_wst_applications.png
#
# Usage:
#   APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#   APP_PAPER_TABLES_DIR=output/paper/tables/application_run_20260414_104327 \
#   APP_PAPER_FIGURES_DIR=output/paper/figures/application_run_20260414_104327 \
#   Rscript scripts/analysis/build_bt_delta_summary.R
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ---- Resolve paths ----------------------------------------------------------
RUN_DIR <- Sys.getenv("APP_RUN_DIR", unset = "")
if (!nzchar(RUN_DIR)) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1L) RUN_DIR <- args[[1]]
}
stopifnot(nzchar(RUN_DIR), dir.exists(RUN_DIR))

run_id <- basename(normalizePath(RUN_DIR))

TABLES_DIR <- Sys.getenv("APP_PAPER_TABLES_DIR", unset = "")
if (!nzchar(TABLES_DIR)) {
  TABLES_DIR <- file.path("output", "paper", "tables", run_id)
}
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

FIGURES_DIR <- Sys.getenv("APP_PAPER_FIGURES_DIR", unset = "")
if (!nzchar(FIGURES_DIR)) {
  FIGURES_DIR <- file.path("output", "paper", "figures", run_id)
}
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Run id      :", run_id, "\n")
cat("Run dir     :", RUN_DIR, "\n")
cat("Tables dir  :", TABLES_DIR, "\n")
cat("Figures dir :", FIGURES_DIR, "\n")

# ---- Dataset metadata -------------------------------------------------------
datasets <- tibble::tibble(
  dataset = c(
    "moreno_sheep", "strauss_2019b", "mountain_goats",
    "macaques_data", "citations_data", "high_school"
  ),
  dataset_label = c(
    "Bighorn sheep", "Spotted hyenas", "Mountain goats",
    "Macaques", "Stat. journals", "High school"
  ),
  x_order = 1:6
)

compute_delta_draw_mean <- function(psi_mat) {
  K <- nrow(psi_mat)
  vals <- numeric(0)
  for (k in 1:(K - 2)) {
    for (l in (k + 1):(K - 1)) {
      for (m in (l + 1):K) {
        vals <- c(vals, psi_mat[k, m] - psi_mat[k, l] - psi_mat[l, m])
      }
    }
  }
  mean(vals)
}

compute_delta_summary <- function(fit_path) {
  fit <- readRDS(fit_path)

  K_trace <- sapply(fit$z, function(z) length(unique(z)))
  K_modal <- as.integer(names(sort(table(K_trace), decreasing = TRUE))[1])
  idx_k <- which(K_trace == K_modal)

  # Keep draws whose psi matrix matches modal K.
  psi_ok <- sapply(idx_k, function(i) is.matrix(fit$psi[[i]]) && nrow(fit$psi[[i]]) == K_modal)
  idx <- idx_k[psi_ok]
  stopifnot(length(idx) > 10)

  draw_means <- sapply(idx, function(i) compute_delta_draw_mean(fit$psi[[i]]))

  # Triple-level range from posterior-mean psi, useful for interpretation.
  psi_mean <- Reduce("+", lapply(idx, function(i) fit$psi[[i]])) / length(idx)
  triple_vals <- numeric(0)
  for (k in 1:(K_modal - 2)) {
    for (l in (k + 1):(K_modal - 1)) {
      for (m in (l + 1):K_modal) {
        triple_vals <- c(triple_vals, psi_mean[k, m] - psi_mean[k, l] - psi_mean[l, m])
      }
    }
  }

  tibble::tibble(
    K_modal = K_modal,
    draws = length(idx),
    n_triples = length(triple_vals),
    delta_mean = mean(draw_means),
    delta_lo = as.numeric(stats::quantile(draw_means, 0.025)),
    delta_hi = as.numeric(stats::quantile(draw_means, 0.975)),
    triple_min = min(triple_vals),
    triple_max = max(triple_vals),
    pct_triples_neg = mean(triple_vals < 0)
  )
}

# ---- Compute summaries ------------------------------------------------------
res <- datasets %>%
  rowwise() %>%
  do({
    ds <- .$dataset
    fit_path <- file.path(RUN_DIR, paste0(ds, "_WST_fit.rds"))
    stopifnot(file.exists(fit_path))
    dplyr::bind_cols(tibble::tibble(dataset = ds), compute_delta_summary(fit_path))
  }) %>%
  ungroup() %>%
  left_join(datasets, by = "dataset") %>%
  arrange(x_order)

# ---- Write CSV --------------------------------------------------------------
csv_path <- file.path(TABLES_DIR, "bt_delta_summary.csv")
readr::write_csv(res, csv_path)
cat("Wrote:", csv_path, "\n")

# ---- Write LaTeX table ------------------------------------------------------
fmt1 <- function(x) formatC(round(x, 2), format = "f", digits = 2)
fmt0 <- function(x) formatC(round(100 * x, 0), format = "f", digits = 0)

tex <- c(
  "\\begin{tabular}{lrrrrc}",
  "\\toprule",
  "Dataset & $\\hat K$ & $\\bar\\Delta$ & $95\\%$ CI & Triple range & $\\%\\Delta<0$ \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  pct_str <- sprintf("%.0f", 100 * r$pct_triples_neg)
  tex <- c(tex, sprintf(
    "%s & %d & $%s$ & $[%s,\\,%s]$ & $[%s,\\,%s]$ & %s\\\\ ",
    r$dataset_label,
    r$K_modal,
    fmt1(r$delta_mean),
    fmt1(r$delta_lo), fmt1(r$delta_hi),
    fmt1(r$triple_min), fmt1(r$triple_max),
    pct_str
  ))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")

tex_path <- file.path(TABLES_DIR, "bt_delta_summary_table.tex")
writeLines(tex, tex_path)
cat("Wrote:", tex_path, "\n")

# ---- Plot: point + 95% interval bars ---------------------------------------
make_plot <- function(path, width = 10, height = 5, point_cex = 1.2) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  par(mar = c(8.5, 5, 2, 1), xpd = NA)

  x <- res$x_order
  y <- res$delta_mean
  lo <- res$delta_lo
  hi <- res$delta_hi

  ylim <- range(c(lo, hi, 0))

  plot(
    x, y,
    type = "p", pch = 16, cex = point_cex,
    xaxt = "n", xlab = "Dataset", ylab = expression(bar(Delta)~"(WST, with 95% interval)"),
    ylim = ylim
  )

  segments(x0 = x, y0 = lo, x1 = x, y1 = hi, lwd = 2)
  segments(x0 = x - 0.12, y0 = lo, x1 = x + 0.12, y1 = lo, lwd = 2)
  segments(x0 = x - 0.12, y0 = hi, x1 = x + 0.12, y1 = hi, lwd = 2)

  # Bradley--Terry benchmark.
  abline(h = 0, lty = 2, lwd = 1.4, col = "firebrick")

  axis(1, at = x, labels = res$dataset_label, las = 2, cex.axis = 0.9)
  mtext("Dashed line is the Bradley--Terry benchmark Delta = 0", side = 3, line = 0.2, cex = 0.85)
}

pdf_path <- file.path(FIGURES_DIR, "bt_delta_wst_applications.pdf")
pdf(pdf_path, width = 10, height = 5)
make_plot(pdf_path)
dev.off()
cat("Wrote:", pdf_path, "\n")

png_path <- file.path(FIGURES_DIR, "bt_delta_wst_applications.png")
png(png_path, width = 1800, height = 900, res = 180)
make_plot(png_path)
dev.off()
cat("Wrote:", png_path, "\n")

# ---- Manifest ---------------------------------------------------------------
manifest <- file.path(TABLES_DIR, "MANIFEST_bt_delta.txt")
writeLines(c(
  paste0("run_id: ", run_id),
  paste0("run_dir: ", normalizePath(RUN_DIR)),
  paste0("generated_by: scripts/analysis/build_bt_delta_summary.R"),
  paste0("generated_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R_version: ", R.version.string)
), manifest)
cat("Wrote:", manifest, "\n")
