#!/usr/bin/env Rscript
# Simulation Results Visualization
# Produces grids of ARI + VI boxplots: x=κ, rows=hierarchy, cols=K, colored by model

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(patchwork)
})
source("scripts/bundle_defaults.R", local = TRUE)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Path to simulation results
SIM_RESULTS_PATH <- Sys.getenv(
  "SIM_RESULTS_PATH",
  bundle_defaults$canonical_simulation_results_csv
)

# Output directory (override with SIM_PLOTS_OUTPUT_DIR)
OUTPUT_DIR <- Sys.getenv("SIM_PLOTS_OUTPUT_DIR", "output/simulation/plots")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading simulation results from:", SIM_RESULTS_PATH, "\n")
sim_data <- read.csv(SIM_RESULTS_PATH)

cat("Dimensions:", nrow(sim_data), "rows x", ncol(sim_data), "columns\n")
cat("Generating models:", paste(unique(sim_data$gen_model), collapse = ", "), "\n")
cat("Fitting models:", paste(unique(sim_data$fit_model), collapse = ", "), "\n")
cat("K values:", paste(sort(unique(sim_data$K_true)), collapse = ", "), "\n")
cat("Kappa values:", paste(sort(unique(sim_data$kappa_mean)), collapse = ", "), "\n")
cat("Psi values:", paste(sort(unique(sim_data$psi_mean)), collapse = ", "), "\n")

# =============================================================================
# DATA PREPARATION
# =============================================================================

# Create factor levels for proper ordering
sim_data <- sim_data %>%
  mutate(
    # Order kappa by increasing sparsity
    kappa_label = factor(format(kappa_mean, trim = TRUE),
                         levels = format(sort(unique(kappa_mean)), trim = TRUE)),
    # Order psi by increasing hierarchy (low to high psi)
    psi_label = factor(paste0("ψ=", psi_mean), 
                       levels = paste0("ψ=", sort(unique(psi_mean)))),
    # Hierarchy category
    hierch_label = factor(hierch, levels = c("weak", "medium", "strong")),
    # K label
    K_label = factor(paste0("K=", K_true), 
                     levels = paste0("K=", sort(unique(K_true)))),
    # Standardize fit_model names
    fit_model_clean = case_when(
      fit_model == "DC-SBM" ~ "DCSBM",
      TRUE ~ fit_model
    ),
    fit_model_clean = factor(fit_model_clean,
                             levels = c("WST", "SST", "DCSBM"),
                             labels = c("WST", "SST", "DCSBM"))
  )

# Facet labels with hierarchy + psi value
hierch_levels <- c("weak", "medium", "strong")
hierch_facet_labels <- sapply(hierch_levels, function(level_name) {
  psi_vals <- sort(unique(sim_data$psi_mean[sim_data$hierch == level_name]))
  psi_text <- if (length(psi_vals) > 0) format(psi_vals[1], trim = TRUE) else "NA"
  paste0("psi=", psi_text)
}, USE.NAMES = FALSE)

sim_data <- sim_data %>%
  mutate(
    hierch_facet = factor(
      hierch,
      levels = hierch_levels,
      labels = hierch_facet_labels
    )
  )

# =============================================================================
# MAIN VISUALIZATION FUNCTION
# =============================================================================

plot_ari_grid <- function(data, gen_model_filter = NULL, title_suffix = "") {
  # Filter by generating model if specified
  plot_data <- data
  if (!is.null(gen_model_filter)) {
    plot_data <- data %>% filter(gen_model == gen_model_filter)
    title_prefix <- paste0("Generated under ", gen_model_filter)
  } else {
    title_prefix <- "All Generating Models"
  }
  
  # Color palette for models
  model_colors <- c(
    "WST" = "#2166AC",    # Blue
    "SST" = "#B2182B",    # Red  
    "DCSBM" = "#4DAF4A"   # Green
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = kappa_label, y = ari, fill = fit_model_clean)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5,
                 position = position_dodge(width = 0.8)) +
    facet_grid(psi_label ~ K_label, scales = "fixed") +
    scale_fill_manual(values = model_colors, name = "Fitting Model") +
    scale_x_discrete() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = paste0("Adjusted Rand Index: ", title_prefix, title_suffix),
      x = expression(bar(kappa)),
      y = "Adjusted Rand Index (ARI)"
    ) +
    theme_bw(base_size = 18) +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 18),
      strip.text.y = element_text(face = "bold", size = 18),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 15),
      axis.title = element_text(size = 17),
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 15, color = "gray40")
    )
  
  return(p)
}

# =============================================================================
# ALTERNATIVE: HIERARCHY ON ROWS (as category labels)
# =============================================================================

plot_ari_grid_hierch <- function(data, gen_model_filter = NULL, title_suffix = "") {
  plot_data <- data
  if (!is.null(gen_model_filter)) {
    plot_data <- data %>% filter(gen_model == gen_model_filter)
    title_prefix <- paste0("Generated under ", gen_model_filter)
  } else {
    title_prefix <- "All Generating Models"
  }
  
  model_colors <- c(
    "WST" = "#2166AC",
    "SST" = "#B2182B",  
    "DCSBM" = "#4DAF4A"
  )
  
  p <- ggplot(plot_data, aes(x = kappa_label, y = ari, fill = fit_model_clean)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5,
                 position = position_dodge(width = 0.8)) +
    facet_grid(hierch_facet ~ K_label, scales = "fixed") +
    scale_fill_manual(values = model_colors, name = "Fitting Model") +
    scale_x_discrete() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = paste0("Adjusted Rand Index: ", title_prefix, title_suffix),
      x = expression(bar(kappa)),
      y = "Adjusted Rand Index (ARI)"
    ) +
    theme_bw(base_size = 18) +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 18),
      strip.text.y = element_text(face = "bold", size = 18),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 15),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 15),
      axis.title = element_text(size = 17),
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 15, color = "gray40")
    )
  
  return(p)
}

# =============================================================================
# GENERIC METRIC PLOTTING (e.g., VI)
# =============================================================================

plot_metric_grid_hierch <- function(
  data,
  metric_col,
  metric_title,
  metric_y_label,
  gen_model_filter = NULL,
  title_suffix = "",
  show_title = TRUE
) {
  plot_data <- data
  if (!is.null(gen_model_filter)) {
    plot_data <- data %>% filter(gen_model == gen_model_filter)
    title_prefix <- paste0("Generated under ", gen_model_filter)
  } else {
    title_prefix <- "All Generating Models"
  }

  plot_data <- plot_data %>%
    filter(is.finite(.data[[metric_col]]))

  model_colors <- c(
    "WST" = "#2166AC",
    "SST" = "#B2182B",
    "DCSBM" = "#4DAF4A"
  )

  y_max <- max(plot_data[[metric_col]], na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) {
    y_max <- 1
  }
  y_breaks <- pretty(c(0, y_max), n = 5)

  p <- ggplot(plot_data, aes(x = kappa_label, y = .data[[metric_col]], fill = fit_model_clean)) +
    geom_boxplot(outlier.size = 1.2, outlier.alpha = 0.5,
                 position = position_dodge(width = 0.8)) +
    facet_grid(hierch_facet ~ K_label, scales = "fixed") +
    scale_fill_manual(values = model_colors, name = "Fitting Model") +
    guides(fill = guide_legend(override.aes = list(alpha = 1, linewidth = 0.6))) +
    scale_x_discrete() +
    scale_y_continuous(limits = c(0, y_max), breaks = y_breaks) +
    labs(
      title = if (show_title) paste0(metric_title, ": ", title_prefix, title_suffix) else NULL,
      x = expression(bar(kappa)),
      y = metric_y_label
    ) +
    theme_bw(base_size = 24) +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 28),
      strip.text.y = element_text(face = "bold", size = 28),
      legend.position = "bottom",
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 23),
      legend.key.width = grid::unit(3.8, "lines"),
      legend.key.height = grid::unit(1.6, "lines"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
      axis.text.y = element_text(size = 23),
      axis.title = element_text(size = 25),
      plot.title = element_text(face = "bold", size = 24),
      plot.subtitle = element_text(size = 18, color = "gray40")
    )

  return(p)
}

# =============================================================================
# GENERATE PLOTS
# =============================================================================

cat("\nGenerating plots...\n")

# Readability subset for model-specific panels: drop middle hierarchy row and K=5 column.
readable_subset <- sim_data %>%
  filter(hierch %in% c("weak", "strong"), K_true %in% c(3, 8))

# 1. All generating models combined
p_all <- plot_ari_grid_hierch(sim_data)
ggsave(file.path(OUTPUT_DIR, "ari_boxplot_all.pdf"), p_all, 
       width = 12, height = 10)
ggsave(file.path(OUTPUT_DIR, "ari_boxplot_all.png"), p_all, 
       width = 12, height = 10, dpi = 300)
cat("  Saved: ari_boxplot_all.pdf/png\n")

# 2. SST-generated data only
if (any(sim_data$gen_model == "SST")) {
  p_sst <- plot_ari_grid_hierch(readable_subset, gen_model_filter = "SST")
  ggsave(file.path(OUTPUT_DIR, "ari_boxplot_SST_gen.pdf"), p_sst, 
         width = 12, height = 10)
  ggsave(file.path(OUTPUT_DIR, "ari_boxplot_SST_gen.png"), p_sst, 
         width = 12, height = 10, dpi = 300)
  cat("  Saved: ari_boxplot_SST_gen.pdf/png\n")
} else {
  cat("  Skipped: ari_boxplot_SST_gen (no SST-generated rows)\n")
}

# 3. WST-generated data only
if (any(sim_data$gen_model == "WST")) {
  p_wst <- plot_ari_grid_hierch(readable_subset, gen_model_filter = "WST")
  ggsave(file.path(OUTPUT_DIR, "ari_boxplot_WST_gen.pdf"), p_wst, 
         width = 12, height = 10)
  ggsave(file.path(OUTPUT_DIR, "ari_boxplot_WST_gen.png"), p_wst, 
         width = 12, height = 10, dpi = 300)
  cat("  Saved: ari_boxplot_WST_gen.pdf/png\n")
} else {
  cat("  Skipped: ari_boxplot_WST_gen (no WST-generated rows)\n")
}

# =============================================================================
# VI PLOTS (same layout as ARI)
# =============================================================================

cat("\nGenerating VI plots...\n")

p_vi_all <- plot_metric_grid_hierch(sim_data, metric_col = "vi",
          metric_title = "Variation of Information",
          metric_y_label = "VI distance")
ggsave(file.path(OUTPUT_DIR, "vi_boxplot_all.pdf"), p_vi_all,
  width = 16, height = 13)
ggsave(file.path(OUTPUT_DIR, "vi_boxplot_all.png"), p_vi_all,
  width = 16, height = 13, dpi = 300)
cat("  Saved: vi_boxplot_all.pdf/png\n")

if (any(sim_data$gen_model == "SST")) {
  p_vi_sst <- plot_metric_grid_hierch(readable_subset, metric_col = "vi",
            metric_title = "Variation of Information",
            metric_y_label = "VI distance",
            gen_model_filter = "SST",
            show_title = FALSE)
  ggsave(file.path(OUTPUT_DIR, "vi_boxplot_SST_gen.pdf"), p_vi_sst,
    width = 16, height = 13)
  ggsave(file.path(OUTPUT_DIR, "vi_boxplot_SST_gen.png"), p_vi_sst,
    width = 16, height = 13, dpi = 300)
  cat("  Saved: vi_boxplot_SST_gen.pdf/png\n")
} else {
  cat("  Skipped: vi_boxplot_SST_gen (no SST-generated rows)\n")
}

if (any(sim_data$gen_model == "WST")) {
  p_vi_wst <- plot_metric_grid_hierch(readable_subset, metric_col = "vi",
            metric_title = "Variation of Information",
            metric_y_label = "VI distance",
            gen_model_filter = "WST",
            show_title = FALSE)
  ggsave(file.path(OUTPUT_DIR, "vi_boxplot_WST_gen.pdf"), p_vi_wst,
    width = 16, height = 13)
  ggsave(file.path(OUTPUT_DIR, "vi_boxplot_WST_gen.png"), p_vi_wst,
    width = 16, height = 13, dpi = 300)
  cat("  Saved: vi_boxplot_WST_gen.pdf/png\n")
} else {
  cat("  Skipped: vi_boxplot_WST_gen (no WST-generated rows)\n")
}

# =============================================================================
# COMBINED SIDE-BY-SIDE PLOT (SST-gen | WST-gen)
# =============================================================================

if (all(c("SST", "WST") %in% unique(sim_data$gen_model))) {
  cat("\nGenerating combined side-by-side plot...\n")

  p_sst_compact <- plot_ari_grid_hierch(sim_data, gen_model_filter = "SST") +
    labs(title = "SST-Generated Data") +
    theme(legend.position = "none")

  p_wst_compact <- plot_ari_grid_hierch(sim_data, gen_model_filter = "WST") +
    labs(title = "WST-Generated Data")

  p_combined <- p_sst_compact + p_wst_compact + 
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Adjusted Rand Index by Generating Model",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, color = "gray40")
      )
    )

  ggsave(file.path(OUTPUT_DIR, "ari_boxplot_combined.pdf"), p_combined, 
         width = 20, height = 10)
  ggsave(file.path(OUTPUT_DIR, "ari_boxplot_combined.png"), p_combined, 
         width = 20, height = 10, dpi = 300)
  cat("  Saved: ari_boxplot_combined.pdf/png\n")
} else {
  cat("\nSkipped combined side-by-side ARI plot (requires both SST and WST generated rows).\n")
}

# =============================================================================
# COMBINED SIDE-BY-SIDE PLOT (VI: SST-gen | WST-gen)
# =============================================================================

if (all(c("SST", "WST") %in% unique(sim_data$gen_model))) {
  cat("\nGenerating combined side-by-side VI plot...\n")

  p_vi_sst_compact <- plot_metric_grid_hierch(readable_subset, metric_col = "vi",
                                              metric_title = "VI distance",
                                              metric_y_label = "VI distance",
                                              gen_model_filter = "SST") +
    labs(title = "SST-Generated Data") +
    theme(legend.position = "none")

  p_vi_wst_compact <- plot_metric_grid_hierch(readable_subset, metric_col = "vi",
                                              metric_title = "VI distance",
                                              metric_y_label = "VI distance",
                                              gen_model_filter = "WST") +
    labs(title = "WST-Generated Data")

  p_vi_combined <- p_vi_sst_compact + p_vi_wst_compact +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "VI distance by Generating Model",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 11, color = "gray40")
      )
    )

  ggsave(file.path(OUTPUT_DIR, "vi_boxplot_combined.pdf"), p_vi_combined,
      width = 22, height = 11)
  ggsave(file.path(OUTPUT_DIR, "vi_boxplot_combined.png"), p_vi_combined,
      width = 22, height = 11, dpi = 300)
  cat("  Saved: vi_boxplot_combined.pdf/png\n")
} else {
  cat("\nSkipped combined side-by-side VI plot (requires both SST and WST generated rows).\n")
}

# =============================================================================
# SUMMARY STATISTICS TABLE
# =============================================================================

cat("\nComputing summary statistics...\n")

summary_stats <- sim_data %>%
  group_by(gen_model, fit_model_clean, K_true, hierch_label, kappa_mean) %>%
  summarise(
    n = n(),
    ari_mean = mean(ari, na.rm = TRUE),
    ari_sd = sd(ari, na.rm = TRUE),
    ari_median = median(ari, na.rm = TRUE),
    ari_q25 = quantile(ari, 0.25, na.rm = TRUE),
    ari_q75 = quantile(ari, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_stats, file.path(OUTPUT_DIR, "ari_summary_stats.csv"), 
          row.names = FALSE)
cat("  Saved: ari_summary_stats.csv\n")

# VI summary stats
vi_summary_stats <- sim_data %>%
  group_by(gen_model, fit_model_clean, K_true, hierch_label, kappa_mean) %>%
  summarise(
    n = sum(is.finite(vi)),
    vi_mean = mean(vi, na.rm = TRUE),
    vi_sd = sd(vi, na.rm = TRUE),
    vi_median = median(vi, na.rm = TRUE),
    vi_q25 = quantile(vi, 0.25, na.rm = TRUE),
    vi_q75 = quantile(vi, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(vi_summary_stats, file.path(OUTPUT_DIR, "vi_summary_stats.csv"),
          row.names = FALSE)
cat("  Saved: vi_summary_stats.csv\n")

# Print brief summary
cat("\n=== ARI Summary by Fitting Model ===\n")
sim_data %>%
  group_by(fit_model_clean) %>%
  summarise(
    mean_ari = round(mean(ari, na.rm = TRUE), 3),
    sd_ari = round(sd(ari, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

cat("\n=== ARI Summary by Gen Model × Fit Model ===\n")
sim_data %>%
  group_by(gen_model, fit_model_clean) %>%
  summarise(
    mean_ari = round(mean(ari, na.rm = TRUE), 3),
    sd_ari = round(sd(ari, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

cat("\n=== VI Summary by Fitting Model ===\n")
sim_data %>%
  group_by(fit_model_clean) %>%
  summarise(
    mean_vi = round(mean(vi, na.rm = TRUE), 3),
    sd_vi = round(sd(vi, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

cat("\n=== VI Summary by Gen Model × Fit Model ===\n")
sim_data %>%
  group_by(gen_model, fit_model_clean) %>%
  summarise(
    mean_vi = round(mean(vi, na.rm = TRUE), 3),
    sd_vi = round(sd(vi, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

cat("\nAll simulation visualizations complete!\n")
cat("Output directory:", OUTPUT_DIR, "\n")
