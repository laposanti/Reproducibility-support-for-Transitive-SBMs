#!/usr/bin/env Rscript
# =============================================================================
# regen_paper_network_figs.R
#
# Regenerate only the network figures embedded in the paper after the
# SST relabel and the layout polish (no Block 1 ... K legend
# on the right, larger nodes, tighter margins). All MCMC fits are reused from
# the existing run dir; nothing is re-fit.
#
# Outputs (to output/paper/figures/<run_id>/):
#   <ds>_combined_block_networks_clean.{pdf,png}     for every dataset
#   moreno_sheep_SST_network.png                     for the side-by-side
#   moreno_sheep_DCSBM_network.png                   in the manuscript
# =============================================================================

suppressPackageStartupMessages({
  source("scripts/analysis/osbm_visualization.R", chdir = FALSE)
  source("scripts/analysis/post_processing_helpers.R", chdir = FALSE)
})
source("scripts/bundle_defaults.R", local = TRUE)

RUN_DIR <- bundle_resolve_application_run_dir(must_exist = TRUE)
stopifnot(dir.exists(RUN_DIR))
run_id <- basename(normalizePath(RUN_DIR))
PAPER_FIG_DIR <- Sys.getenv("APP_PAPER_FIGURES_DIR", unset = "")
if (!nzchar(PAPER_FIG_DIR)) {
  PAPER_FIG_DIR <- file.path("output", "paper", "figures", run_id)
}
fs::dir_create(PAPER_FIG_DIR)

# Single source of truth: read the canonical cube produced by
# scripts/analysis/build_post_processing.R. All z_hat / violation numbers
# rendered in the paper figures come from there.
cube_meta <- read_cube_meta(RUN_DIR)

DATASETS <- c("moreno_sheep", "strauss_2019b", "mountain_goats",
              "citations_data", "macaques_data", "high_school")

results_list <- list()
for (ds in DATASETS) {
  cat("==", ds, "==\n")
  cube_files <- file.path(cube_meta$per_fit_dir,
                          paste0(ds, "__", c("WST", "SST", "DCSBM"), ".rds"))
  if (!all(file.exists(cube_files))) {
    cat("  cube entries missing, skip\n"); next
  }
  fit_wst   <- load_fit_cube(ds, "WST",   cube_meta)
  fit_sst   <- load_fit_cube(ds, "SST",   cube_meta)
  fit_dcsbm <- load_fit_cube(ds, "DCSBM", cube_meta)
  A         <- fit_wst$A

  results_list[[ds]] <- list(
    A = A,
    z_hat_wst   = fit_wst$z_hat,
    z_hat_sst   = fit_sst$z_hat,
    z_hat_dcsbm = fit_dcsbm$z_hat,
    viol_wst    = get_violation_rate_cube(ds, "WST",   cube_meta),
    viol_sst    = get_violation_rate_cube(ds, "SST",   cube_meta),
    viol_dcsbm  = get_violation_rate_cube(ds, "DCSBM", cube_meta)
  )

  p_clean <- plot_combined_block_networks_clean(
    A, fit_wst$z_hat, fit_sst$z_hat, fit_dcsbm$z_hat, ds)
  ggplot2::ggsave(file.path(PAPER_FIG_DIR,
                            paste0(ds, "_combined_block_networks_clean.pdf")),
                  p_clean, width = 7.2, height = 2.95)
  ggplot2::ggsave(file.path(PAPER_FIG_DIR,
                            paste0(ds, "_combined_block_networks_clean.png")),
                  p_clean, width = 7.2, height = 2.95, dpi = 300)
  cat("  wrote combined block networks\n")
}

# Sheep individual ordered networks for the backward-edge-reduction figure
if (!is.null(results_list[["moreno_sheep"]])) {
  r <- results_list[["moreno_sheep"]]

  plot_simple_network_layout <- function(A, z_hat, layout_engine = c("fr", "stress"),
                                         subtitle_text = NULL, seed = 11,
                                         show_node_labels = TRUE) {
    layout_engine <- match.arg(layout_engine)
    if (inherits(A, "Matrix")) A <- as.matrix(A)
    diag(A) <- 0

    g <- igraph::graph_from_adjacency_matrix(
      A,
      mode = "directed",
      weighted = TRUE,
      diag = FALSE
    )
    V(g)$name <- if (!is.null(rownames(A))) rownames(A) else as.character(seq_len(nrow(A)))
    K <- max(z_hat, na.rm = TRUE)
    V(g)$tier <- factor(z_hat, levels = seq_len(K), labels = paste("Tier", seq_len(K)))

    g_und <- igraph::as_undirected(g, mode = "collapse", edge.attr.comb = list(weight = "sum"))
    if (layout_engine == "stress" && requireNamespace("graphlayouts", quietly = TRUE)) {
      lay <- graphlayouts::layout_with_stress(g_und, weights = E(g_und)$weight)
    } else {
      lay <- igraph::layout_with_fr(g_und, weights = E(g_und)$weight, niter = 2000)
    }
    layout_mat <- as.matrix(lay)

    ggraph::ggraph(g, layout = layout_mat) +
      ggraph::geom_edge_link(
        colour = "grey62",
        alpha = 0.38,
        linewidth = 0.34,
        arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed"),
        end_cap = ggraph::circle(2.7, "mm")
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(fill = tier),
        shape = 21,
        size = 3.6,
        stroke = 0.45,
        colour = "grey12"
      ) +
      {
        if (show_node_labels) {
          ggraph::geom_node_text(
            ggplot2::aes(label = name),
            size = 2.3,
            repel = TRUE,
            max.overlaps = 40,
            segment.color = NA,
            colour = "grey10"
          )
        }
      } +
      ggplot2::scale_fill_viridis_d(option = "C", name = "Tiered block", drop = FALSE) +
      ggplot2::labs(caption = subtitle_text) +
      ggplot2::coord_equal(clip = "off") +
      ggplot2::theme_void(base_size = 11) +
      ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(
          hjust = 0.5,
          size = 11,
          face = "bold",
          colour = "grey35",
          margin = ggplot2::margin(t = 2)
        ),
        plot.margin = ggplot2::margin(1, 1, 2, 1),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 10, face = "bold"),
        legend.text = ggplot2::element_text(size = 9),
        legend.margin = ggplot2::margin(t = -2),
        legend.box.margin = ggplot2::margin(t = -2)
      )
  }

  p_simple_fr <- plot_simple_network_layout(
    A = r$A,
    z_hat = r$z_hat_sst,
    layout_engine = "fr",
    subtitle_text = "Moreno sheep -- simple FR"
  )
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_simple_fr.png"),
                  p_simple_fr, width = 7.5, height = 8.5, dpi = 360, bg = "transparent")
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_simple_fr.pdf"),
                  p_simple_fr, width = 7.5, height = 8.5, bg = "transparent")

  p_simple_stress <- plot_simple_network_layout(
    A = r$A,
    z_hat = r$z_hat_sst,
    layout_engine = "stress",
    subtitle_text = "Moreno sheep -- simple stress"
  )
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_simple_stress.png"),
                  p_simple_stress, width = 7.5, height = 8.5, dpi = 360, bg = "transparent")
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_simple_stress.pdf"),
                  p_simple_stress, width = 7.5, height = 8.5, bg = "transparent")

  # Ordered stress/force layout variants requested for sheep.
  p_sst_stress <- plot_ordered_network_stress(
    A = r$A,
    z_hat = r$z_hat_sst,
    subtitle_text = "SST",
    layout_engine = "stress",
    tier_strength = 0.80,
    tier_spacing = 1.45,
    within_tier_jitter = 0.13,
    node_size_scale = 3.1,
    seed = 11,
    show_tier_guides = FALSE,
    show_node_labels = TRUE,
    precomputed_rate = r$viol_sst$rate,
    precomputed_backward = r$viol_sst$backward_mass,
    precomputed_cross = r$viol_sst$cross_mass
  )
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_SST_network_ordered_stress.png"),
                  p_sst_stress, width = 7.5, height = 8.5, dpi = 360, bg = "transparent")
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_SST_network_ordered_stress.pdf"),
                  p_sst_stress, width = 7.5, height = 8.5, bg = "transparent")

  p_sst_fr <- plot_ordered_network_stress(
    A = r$A,
    z_hat = r$z_hat_sst,
    subtitle_text = "SST",
    layout_engine = "fr",
    tier_strength = 0.80,
    tier_spacing = 1.45,
    within_tier_jitter = 0.13,
    node_size_scale = 3.1,
    seed = 11,
    show_tier_guides = FALSE,
    show_node_labels = TRUE,
    precomputed_rate = r$viol_sst$rate,
    precomputed_backward = r$viol_sst$backward_mass,
    precomputed_cross = r$viol_sst$cross_mass
  )
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_SST_network_ordered_fr.png"),
                  p_sst_fr, width = 7.5, height = 8.5, dpi = 360, bg = "transparent")
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_SST_network_ordered_fr.pdf"),
                  p_sst_fr, width = 7.5, height = 8.5, bg = "transparent")

  p_dc_stress <- plot_ordered_network_stress(
    A = r$A,
    z_hat = r$z_hat_dcsbm,
    subtitle_text = "DC-SBM",
    layout_engine = "stress",
    tier_strength = 0.82,
    tier_spacing = 1.55,
    within_tier_jitter = 0.16,
    node_size_scale = 3.1,
    seed = 11,
    show_tier_guides = FALSE,
    show_node_labels = TRUE,
    precomputed_rate = r$viol_dcsbm$rate,
    precomputed_backward = r$viol_dcsbm$backward_mass,
    precomputed_cross = r$viol_dcsbm$cross_mass
  )
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_DCSBM_network_ordered_stress.png"),
                  p_dc_stress, width = 7.5, height = 8.5, dpi = 360, bg = "transparent")
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_DCSBM_network_ordered_stress.pdf"),
                  p_dc_stress, width = 7.5, height = 8.5, bg = "transparent")

  p_dc_fr <- plot_ordered_network_stress(
    A = r$A,
    z_hat = r$z_hat_dcsbm,
    subtitle_text = "DC-SBM",
    layout_engine = "fr",
    tier_strength = 0.82,
    tier_spacing = 1.55,
    within_tier_jitter = 0.16,
    node_size_scale = 3.1,
    seed = 11,
    show_tier_guides = FALSE,
    show_node_labels = TRUE,
    precomputed_rate = r$viol_dcsbm$rate,
    precomputed_backward = r$viol_dcsbm$backward_mass,
    precomputed_cross = r$viol_dcsbm$cross_mass
  )
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_DCSBM_network_ordered_fr.png"),
                  p_dc_fr, width = 7.5, height = 8.5, dpi = 360, bg = "transparent")
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_DCSBM_network_ordered_fr.pdf"),
                  p_dc_fr, width = 7.5, height = 8.5, bg = "transparent")

  p_sst <- plot_ordered_network(
    r$A, r$z_hat_sst,
    layout_style = "fr_ordered",
    title = "Bighorn sheep -- SST",
    subtitle_text = "SST",
    precomputed_rate     = r$viol_sst$rate,
    precomputed_backward = r$viol_sst$backward_mass,
    precomputed_cross    = r$viol_sst$cross_mass
  )
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_SST_network.png"),
                  p_sst, width = 7, height = 8.5, dpi = 300)
  p_dc <- plot_ordered_network(
    r$A, r$z_hat_dcsbm,
    layout_style = "fr_ordered",
    title = "Bighorn sheep -- DC-SBM",
    subtitle_text = "DC-SBM",
    precomputed_rate     = r$viol_dcsbm$rate,
    precomputed_backward = r$viol_dcsbm$backward_mass,
    precomputed_cross    = r$viol_dcsbm$cross_mass
  )
  ggplot2::ggsave(file.path(PAPER_FIG_DIR, "moreno_sheep_DCSBM_network.png"),
                  p_dc, width = 7, height = 8.5, dpi = 300)

  # Match effective vertical span for tier_line so SST is not visually compressed.
  # Keep save dimensions identical between SST and DC-SBM tier_line exports.
  K_sst <- max(r$z_hat_sst, na.rm = TRUE)
  K_dc <- max(r$z_hat_dcsbm, na.rm = TRUE)
  base_tier_spacing <- 1.9
  target_tier_span <- base_tier_spacing * max(c(K_sst - 1, K_dc - 1, 1))
  full_tier_spacing_sst <- if (K_sst > 1) target_tier_span / (K_sst - 1) else base_tier_spacing
  full_tier_spacing_dc <- if (K_dc > 1) target_tier_span / (K_dc - 1) else base_tier_spacing
  tier_balance_blend <- 0.45
  tier_spacing_sst <- (base_tier_spacing + tier_balance_blend * (full_tier_spacing_sst - base_tier_spacing)) * 0.88
  tier_spacing_dc <- (base_tier_spacing + tier_balance_blend * (full_tier_spacing_dc - base_tier_spacing)) * 0.88
  tier_line_spacing_scale_sst <- 1.22
  tier_line_spacing_scale_dc <- 1.00
  alt_cmp_width <- 7.5
  alt_cmp_height <- 8.5
  tier_line_cmp_width <- 9.6
  tier_line_cmp_height <- 6.8

  # Standard-layout comparisons for the sheep network (same z_hat and violation stats).
  for (style in c("fr_ordered", "stress_ordered", "kk_ordered", "tier_line")) {
    p_sst_alt <- plot_ordered_network(
      r$A, r$z_hat_sst,
      layout_style = style,
      tier_spacing = if (identical(style, "tier_line")) tier_spacing_sst else base_tier_spacing,
      tier_line_spacing_scale = if (identical(style, "tier_line")) tier_line_spacing_scale_sst else 1.0,
      title = paste("Bighorn sheep -- SST (", style, ")", sep = ""),
      subtitle_text = paste("SST [", style, "]", sep = ""),
      precomputed_rate     = r$viol_sst$rate,
      precomputed_backward = r$viol_sst$backward_mass,
      precomputed_cross    = r$viol_sst$cross_mass
    )
    this_width <- if (identical(style, "tier_line")) tier_line_cmp_width else alt_cmp_width
    this_height <- if (identical(style, "tier_line")) tier_line_cmp_height else alt_cmp_height
    ggplot2::ggsave(file.path(PAPER_FIG_DIR, paste0("moreno_sheep_SST_network_", style, ".png")),
            p_sst_alt, width = this_width, height = this_height, dpi = 300)

    p_dc_alt <- plot_ordered_network(
      r$A, r$z_hat_dcsbm,
      layout_style = style,
      tier_spacing = if (identical(style, "tier_line")) tier_spacing_dc else base_tier_spacing,
      tier_line_spacing_scale = if (identical(style, "tier_line")) tier_line_spacing_scale_dc else 1.0,
      title = paste("Bighorn sheep -- DC-SBM (", style, ")", sep = ""),
      subtitle_text = paste("DC-SBM [", style, "]", sep = ""),
      precomputed_rate     = r$viol_dcsbm$rate,
      precomputed_backward = r$viol_dcsbm$backward_mass,
      precomputed_cross    = r$viol_dcsbm$cross_mass
    )
    ggplot2::ggsave(file.path(PAPER_FIG_DIR, paste0("moreno_sheep_DCSBM_network_", style, ".png")),
            p_dc_alt, width = this_width, height = this_height, dpi = 300)
  }
  cat("Wrote sheep individual ordered networks\n")
}

cat("Done.\n")
