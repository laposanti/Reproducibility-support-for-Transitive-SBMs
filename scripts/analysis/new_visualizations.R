# =============================================================================
# NEW EXPERIMENTAL VISUALIZATIONS
# =============================================================================
# scripts/analysis/new_visualizations.R
#
# 4 new plots per (dataset, model) combination:
#   1. Arc diagram — 1D hierarchical layout; forward above, backward below
#   2. Block membership probability heatmap — P(z_i = k) uncertainty
#   3. Posterior rank uncertainty caterpillar — median + 50%/95% CIs
#   4. Concentric ring network — radial layout, inner ring = strongest block
#
# Dependencies:
#   Requires osbm_visualization.R to have been sourced into the environment
#   (for find_fit_file_local, choose_dataset_local, get_z_hat_from_draws).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(scales)
  library(viridis)
})

# Optional canonical cube readers (read_cube_meta / load_fit_cube). Sourced
# from scripts/analysis/post_processing_helpers.R so that this script reads
# z_hat from the cube when available rather than recomputing it.
if (!exists("read_cube_meta", mode = "function")) {
  .pph <- "scripts/analysis/post_processing_helpers.R"
  if (file.exists(.pph)) source(.pph, chdir = FALSE)
}

# =============================================================================
# 1. ARC DIAGRAM
# =============================================================================
#' Nodes on a horizontal line ordered by block (strong→weak) then within-block
#' success. Forward edges (strong→weak) arch above in blue; backward edges
#' (violations) arch below in red; within-block edges are subtle gray above.
plot_arc_diagram <- function(A, z_hat, z_chain = NULL,
                             title = "Arc Diagram",
                             max_edges = 600,
                             show_labels = NULL) {
  n <- nrow(A)
  node_names <- rownames(A)
  if (is.null(node_names)) node_names <- as.character(seq_len(n))
  if (is.null(show_labels)) show_labels <- (n <= 40)

  node_success <- rowSums(A) / pmax(1, rowSums(A + t(A)))
  K <- max(z_hat, na.rm = TRUE)

  # Node ordering: block 1 (strongest) at left, block K at right,
  # within-block by descending success

  ordering <- integer(0)
  for (k in seq_len(K)) {
    idx <- which(z_hat == k)
    ordering <- c(ordering, idx[order(-node_success[idx])])
  }

  pos <- integer(n)
  pos[ordering] <- seq_len(n)

  # --- Edge data ---
  ij <- which(A > 0, arr.ind = TRUE)
  edge_df <- data.frame(
    x    = pos[ij[, 1]],
    xend = pos[ij[, 2]],
    weight    = A[ij],
    from_block = z_hat[ij[, 1]],
    to_block   = z_hat[ij[, 2]],
    stringsAsFactors = FALSE
  )
  edge_df$direction <- with(edge_df, ifelse(
    from_block < to_block, "forward",
    ifelse(from_block > to_block, "backward", "within")
  ))

  # Edge filtering for large networks: keep all backward, sample rest
  if (nrow(edge_df) > max_edges) {
    back <- edge_df[edge_df$direction == "backward", ]
    rest <- edge_df[edge_df$direction != "backward", ]
    n_rest <- max(0L, max_edges - nrow(back))
    if (nrow(rest) > n_rest) {
      set.seed(42)
      rest <- rest[sample.int(nrow(rest), n_rest), ]
    }
    edge_df <- rbind(back, rest)
  }

  fwd <- edge_df[edge_df$direction == "forward", , drop = FALSE]
  bwd <- edge_df[edge_df$direction == "backward", , drop = FALSE]
  wth <- edge_df[edge_df$direction == "within", , drop = FALSE]
  # Normalise within-block edges to left→right for consistent curvature
  if (nrow(wth) > 0) {
    swap <- wth$x > wth$xend
    tmp <- wth$x[swap]; wth$x[swap] <- wth$xend[swap]; wth$xend[swap] <- tmp
  }

  # --- Node data ---
  node_df <- data.frame(
    x = seq_len(n), block = factor(z_hat[ordering]),
    name = node_names[ordering], stringsAsFactors = FALSE
  )

  # --- Block boundaries / labels ---
  block_tbl   <- table(factor(z_hat[ordering], levels = seq_len(K)))
  block_bounds <- cumsum(block_tbl)
  block_starts <- c(0, block_bounds[-K])
  block_mids   <- (block_starts + block_bounds) / 2

  # Backward-edge summary
  total_cross <- sum(edge_df$weight[edge_df$direction != "within"])
  total_back  <- sum(edge_df$weight[edge_df$direction == "backward"])
  pct_back    <- if (total_cross > 0) 100 * total_back / total_cross else 0

  # Y-scaling reference for annotations
  y_ref <- n * 0.22

  # --- Block background rectangles ---
  rects <- data.frame(
    xmin = c(0.5, block_bounds[-K] + 0.5),
    xmax = block_bounds + 0.5,
    fill = factor(seq_len(K))
  )

  p <- ggplot() +
    geom_rect(data = rects,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
              alpha = 0.06, show.legend = FALSE) +
    scale_fill_viridis_d(option = "turbo", guide = "none")

  # Forward arcs — above line (blue)
  # For forward edges from_block < to_block, so their x < xend due to ordering.
  # Direction L→R, curvature < 0 ⇒ bend left of direction ⇒ above.
  if (nrow(fwd) > 0)
    p <- p + geom_curve(
      data = fwd, aes(x = x, xend = xend, y = 0, yend = 0, alpha = weight),
      curvature = -0.3, color = "#2166AC", linewidth = 0.22,
      show.legend = FALSE)

  # Within-block arcs — subtle above (gray), already normalised L→R
  if (nrow(wth) > 0)
    p <- p + geom_curve(
      data = wth, aes(x = x, xend = xend, y = 0, yend = 0),
      curvature = -0.15, color = "grey60", alpha = 0.06,
      linewidth = 0.15, show.legend = FALSE)

  # Backward arcs — below line (red)
  # For backward edges from_block > to_block, so x > xend.
  # Direction R→L, curvature < 0 ⇒ bend left of R→L direction ⇒ below.
  if (nrow(bwd) > 0)
    p <- p + geom_curve(
      data = bwd, aes(x = x, xend = xend, y = 0, yend = 0, alpha = weight),
      curvature = -0.4, color = "#B2182B", linewidth = 0.3,
      show.legend = FALSE)

  # Block boundary lines
  if (K > 1)
    p <- p + geom_vline(xintercept = block_bounds[-K] + 0.5,
                        linetype = "dashed", color = "grey40", linewidth = 0.4)

  # Nodes
  p <- p +
    geom_point(data = node_df, aes(x = x, y = 0, color = block),
               size = ifelse(n > 100, 1.4, 2.5), show.legend = TRUE) +
    scale_color_viridis_d(option = "turbo", name = "Block")

  # Node labels (small networks only)
  if (show_labels)
    p <- p + geom_text(data = node_df,
                       aes(x = x, y = -0.08 * y_ref, label = name),
                       angle = 90, hjust = 1,
                       size = max(1.5, 3 - 0.05 * n))

  # Block labels at top
  p <- p + annotate("text", x = block_mids, y = 0.42 * y_ref,
                    label = paste("Block", seq_len(K)),
                    size = 3, fontface = "bold", color = "grey30")

  # Backward-edge annotation
  p <- p + annotate("text", x = n / 2, y = -0.42 * y_ref,
                    label = sprintf("Backward edges: %.1f%%", pct_back),
                    size = 3.5, color = "#B2182B", fontface = "bold")

  p + scale_alpha_continuous(range = c(0.04, 0.55), guide = "none") +
    labs(title = title, x = "Node (strongest \u2192 weakest)", y = NULL) +
    theme_void(base_size = 11) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title.x  = element_text(size = 9, color = "grey40"),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA)
    )
}


# =============================================================================
# 2. BLOCK MEMBERSHIP PROBABILITY HEATMAP
# =============================================================================
#' Rows = nodes ordered by (z_hat block, certainty).
#' Columns = blocks 1..K. Cell colour = P(z_i = k).
plot_membership_heatmap <- function(z_chain, z_hat, A,
                                    title = "Block Membership Probabilities") {
  n <- ncol(z_chain)
  K <- max(z_hat, na.rm = TRUE)
  node_names <- rownames(A)
  if (is.null(node_names)) node_names <- as.character(seq_len(n))

  # P(z_i = k) matrix  (n × K)
  prob_mat <- matrix(0, n, K)
  for (k in seq_len(K)) prob_mat[, k] <- colMeans(z_chain == k)

  # Entropy per node (bits)
  node_entropy <- apply(prob_mat, 1, function(p) {
    p <- p[p > 0]; -sum(p * log2(p))
  })

  # Order: by z_hat block, within block by descending certainty
  node_certainty <- apply(prob_mat, 1, max)
  ordering <- order(z_hat, -node_certainty)

  # Long-form
  prob_df <- expand.grid(node_idx = seq_len(n), block = seq_len(K))
  prob_df$prob      <- as.vector(prob_mat)
  prob_df$node_name <- factor(node_names[prob_df$node_idx],
                              levels = node_names[ordering])
  prob_df$block     <- factor(prob_df$block)

  # Block boundaries for horizontal separator lines
  z_ordered    <- z_hat[ordering]
  block_bounds <- cumsum(table(factor(z_ordered, levels = seq_len(K))))

  # Entropy sidebar
  ent_df <- data.frame(
    node_name = factor(node_names[ordering], levels = node_names[ordering]),
    entropy   = node_entropy[ordering]
  )

  # --- Main heatmap ---
  p_main <- ggplot(prob_df, aes(x = block, y = node_name, fill = prob)) +
    geom_tile(color = "white", linewidth = 0.05) +
    scale_fill_viridis_c(option = "magma", direction = -1,
                         limits = c(0, 1), name = expression(P(z[i] == k))) +
    { if (K > 1) geom_hline(
        yintercept = n - block_bounds[-K] + 0.5,
        color = "white", linewidth = 0.7) } +
    labs(x = "Block", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y  = if (n > 50) element_blank() else element_text(size = 5),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(expand = FALSE)

  # --- Entropy sidebar ---
  p_ent <- ggplot(ent_df, aes(x = 1, y = node_name, fill = entropy)) +
    geom_tile(show.legend = TRUE) +
    scale_fill_gradient(low = "grey95", high = "#D7191C",
                        name = "Entropy\n(bits)") +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 11) +
    theme(
      axis.text.x  = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )

  # Combine with patchwork
  (p_main + p_ent +
     patchwork::plot_layout(widths = c(6, 1), guides = "collect")) +
    patchwork::plot_annotation(
      title = title,
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
}


# =============================================================================
# 3. POSTERIOR RANK CATERPILLAR
# =============================================================================
#' Per-node rank uncertainty: median point + 50% thick bar + 95% thin bar.
#' Rank is block-first (block 1 = rank 1-n1), with within-block success
#' tie-breaking.
plot_rank_caterpillar <- function(z_chain, z_hat, A,
                                 title = "Posterior Rank Uncertainty") {
  n <- ncol(z_chain)
  S <- nrow(z_chain)
  K <- max(z_hat, na.rm = TRUE)
  node_names <- rownames(A)
  if (is.null(node_names)) node_names <- as.character(seq_len(n))

  # Stable within-block tie-breaker from empirical success
  node_success <- rowSums(A) / pmax(1, rowSums(A + t(A)))
  eps <- rank(node_success, ties.method = "average") / (n + 1)

  # Rank for every MCMC draw
  rank_mat <- matrix(NA_integer_, S, n)
  for (s in seq_len(S)) {
    z_s   <- z_chain[s, ]
    score <- -as.numeric(z_s) + 1e-3 * eps
    rank_mat[s, ] <- rank(-score, ties.method = "average")
  }

  med  <- apply(rank_mat, 2, median)
  q025 <- apply(rank_mat, 2, quantile, 0.025)
  q975 <- apply(rank_mat, 2, quantile, 0.975)
  q25  <- apply(rank_mat, 2, quantile, 0.25)
  q75  <- apply(rank_mat, 2, quantile, 0.75)

  ord <- order(med)
  cat_df <- data.frame(
    x = seq_len(n), node = node_names[ord], block = factor(z_hat[ord]),
    med = med[ord], q025 = q025[ord], q975 = q975[ord],
    q25 = q25[ord], q75 = q75[ord], stringsAsFactors = FALSE
  )
  cat_df$node <- factor(cat_df$node, levels = cat_df$node)

  ggplot(cat_df, aes(x = x, y = med, color = block)) +
    # 95% CI (thin)
    geom_segment(aes(x = x, xend = x, y = q025, yend = q975),
                 linewidth = 0.3, alpha = 0.55) +
    # 50% CI (thick)
    geom_segment(aes(x = x, xend = x, y = q25, yend = q75),
                 linewidth = 1.0, alpha = 0.75) +
    # Median
    geom_point(size = 1.1) +
    # Perfect-agreement line
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                color = "grey55", linewidth = 0.3) +
    scale_color_viridis_d(option = "turbo", name = "Block") +
    scale_y_reverse() +
    labs(title = title,
         x = "Node (ordered by median rank)",
         y = "Rank  (1 = strongest)",
         subtitle = "Thin bar = 95% CI  \u2022  Thick bar = 50% CI") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
      axis.text.x   = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}


# =============================================================================
# 4. CONCENTRIC RING NETWORK
# =============================================================================
#' Block 1 (strongest) on the innermost ring, Block K on the outer ring.
#' Nodes evenly spaced within their ring, ordered by success.
#' Edges: blue = forward (inner→outer), red = backward (violation).
plot_concentric_network <- function(A, z_hat,
                                    title = "Concentric Ring Network",
                                    max_edges = 500) {
  n <- nrow(A)
  node_names <- rownames(A)
  if (is.null(node_names)) node_names <- as.character(seq_len(n))
  K <- max(z_hat, na.rm = TRUE)

  node_success <- rowSums(A) / pmax(1, rowSums(A + t(A)))

  # Compute positions: each block is a concentric ring
  node_x <- numeric(n)
  node_y <- numeric(n)
  for (k in seq_len(K)) {
    idx  <- which(z_hat == k)
    n_k  <- length(idx)
    if (n_k == 0) next
    radius <- k
    idx_ord <- idx[order(-node_success[idx])]
    angles  <- seq(0, 2 * pi, length.out = n_k + 1)[seq_len(n_k)]
    angles  <- angles + (k - 1) * pi / (2 * K)
    node_x[idx_ord] <- radius * cos(angles)
    node_y[idx_ord] <- radius * sin(angles)
  }

  # Build igraph
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE)
  V(g)$block <- z_hat
  V(g)$name  <- node_names

  # Classify edge direction
  el <- igraph::as_data_frame(g, what = "edges")
  from_idx <- match(el$from, node_names)
  to_idx   <- match(el$to,   node_names)
  direction <- ifelse(
    z_hat[from_idx] < z_hat[to_idx], "forward",
    ifelse(z_hat[from_idx] > z_hat[to_idx], "backward", "within")
  )
  E(g)$direction <- direction
  E(g)$edge_color <- dplyr::case_when(
    direction == "forward"  ~ "#2166AC",
    direction == "backward" ~ "#B2182B",
    TRUE                    ~ "#808080"
  )

  # Edge filtering
  if (igraph::ecount(g) > max_edges) {
    back_idx <- which(E(g)$direction == "backward")
    rest_idx <- setdiff(seq_len(igraph::ecount(g)), back_idx)
    n_keep   <- max(0L, max_edges - length(back_idx))
    if (length(rest_idx) > n_keep) {
      set.seed(42)
      rm_idx <- rest_idx[sample.int(length(rest_idx),
                                    length(rest_idx) - n_keep)]
      g <- igraph::delete_edges(g, rm_idx)
    }
  }

  layout_mat <- cbind(node_x, node_y)

  # Concentric guide circles
  circle_df <- do.call(rbind, lapply(seq_len(K), function(k) {
    theta <- seq(0, 2 * pi, length.out = 200)
    data.frame(x = k * cos(theta), y = k * sin(theta), ring = k)
  }))

  ggraph(g, layout = layout_mat) +
    # Ring guidelines
    geom_path(data = circle_df,
              aes(x = x, y = y, group = ring),
              color = "grey82", linewidth = 0.3, linetype = "dotted",
              inherit.aes = FALSE) +
    # Edges
    geom_edge_arc(aes(colour = I(edge_color), alpha = weight),
                  arrow   = arrow(length = unit(1.2, "mm"), type = "closed"),
                  end_cap = circle(2, "mm"),
                  strength = 0.15,
                  show.legend = FALSE) +
    # Nodes
    geom_node_point(aes(fill = factor(block)),
                    shape = 21, size = ifelse(n > 100, 1.5, 2.8),
                    stroke = 0.3) +
    # Labels for small networks
    { if (n <= 40)
      geom_node_text(aes(label = name), size = 2, repel = TRUE,
                     max.overlaps = 15) } +
    # Ring labels
    annotate("text",
             x = seq_len(K) + 0.35, y = 0.35,
             label = paste("B", seq_len(K)),
             size = 2.8, color = "grey40", fontface = "italic") +
    scale_fill_viridis_d(option = "turbo", name = "Block") +
    scale_edge_alpha(range = c(0.03, 0.45), guide = "none") +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed()
}


# =============================================================================
# ORCHESTRATOR: generate all 4 new plots for one dataset
# =============================================================================
generate_new_visualizations <- function(dataset_name,
                                        ppc_dir    = "output/application/ppc",
                                        output_dir = NULL,
                                        dataset_plots_dir = "output/application/plots") {
  cat("  [new_viz] Generating for:", dataset_name, "\n")

  # New layout: output/application/plots/<dataset>/new_visualizations/
  if (is.null(output_dir) || !nzchar(output_dir)) {
    out_path <- file.path(dataset_plots_dir, dataset_name, "new_visualizations")
  } else {
    out_path <- file.path(output_dir, dataset_name)
  }
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  # Load adjacency matrix (choose_dataset_local from osbm_visualization.R)
  A <- tryCatch(choose_dataset_local(dataset_name), error = function(e) NULL)
  if (is.null(A)) {
    cat("    ERROR loading adjacency for", dataset_name, "\n")
    return(invisible(NULL))
  }
  A <- as.matrix(A)
  n <- nrow(A)

  models <- c("WST", "SST", "DCSBM")

  # Prefer the canonical cube if available — single source of truth.
  cube_meta <- tryCatch(read_cube_meta(), error = function(e) NULL)

  for (model in models) {
    fit_path <- find_fit_file_local(dataset_name, model, ppc_dir)
    if (is.null(fit_path) || !file.exists(fit_path)) {
      cat("    ", model, "fit not found, skipping\n")
      next
    }

    cat("    Processing", model, "(", basename(fit_path), ") ...\n")
    fit <- readRDS(fit_path)

    cube_entry <- if (!is.null(cube_meta))
      tryCatch(load_fit_cube(dataset_name, model, cube_meta),
               error = function(e) NULL) else NULL
    if (!is.null(cube_entry)) {
      draws <- list(z_hat = cube_entry$z_hat,
                    z_chain = cube_entry$z_chain,
                    K = cube_entry$K_hat)
    } else {
      draws <- tryCatch(
        get_z_hat_from_draws(fit, A, model = model),
        error = function(e) {
          cat("      ERROR:", conditionMessage(e), "\n")
          NULL
        }
      )
    }
    if (is.null(draws)) next

    z_hat   <- draws$z_hat
    z_chain <- draws$z_chain
    K       <- draws$K

    cat("      K =", K, ", S =", nrow(z_chain), ", n =", n, "\n")

    prefix <- file.path(out_path, paste0(dataset_name, "_", model))

    # ---- 1. Arc Diagram ----
    cat("      [1/4] Arc diagram ...\n")
    tryCatch({
      p <- plot_arc_diagram(A, z_hat, z_chain,
                            title = paste(model, "Arc Diagram \u2014", dataset_name),
                            max_edges = if (n > 100) 800 else 600)
      ggsave(paste0(prefix, "_arc_diagram.png"), p,
             width = max(8, n * 0.04), height = 5, dpi = 300)
    }, error = function(e) cat("        FAILED:", conditionMessage(e), "\n"))

    # ---- 2. Block Membership Heatmap ----
    cat("      [2/4] Membership heatmap ...\n")
    tryCatch({
      p <- plot_membership_heatmap(z_chain, z_hat, A,
                                   title = paste(model, "Block Membership \u2014",
                                                 dataset_name))
      ggsave(paste0(prefix, "_membership_heatmap.png"), p,
             width = max(4.5, K * 0.8 + 2.5),
             height = max(5, n * 0.055), dpi = 300)
    }, error = function(e) cat("        FAILED:", conditionMessage(e), "\n"))

    # ---- 3. Rank Caterpillar ----
    cat("      [3/4] Rank caterpillar ...\n")
    tryCatch({
      p <- plot_rank_caterpillar(z_chain, z_hat, A,
                                 title = paste(model, "Rank Uncertainty \u2014",
                                               dataset_name))
      ggsave(paste0(prefix, "_rank_caterpillar.png"), p,
             width = max(7, n * 0.04), height = 5.5, dpi = 300)
    }, error = function(e) cat("        FAILED:", conditionMessage(e), "\n"))

    # ---- 4. Concentric Ring Network ----
    cat("      [4/4] Concentric ring ...\n")
    tryCatch({
      p <- plot_concentric_network(A, z_hat,
                                   title = paste(model, "Ring Layout \u2014",
                                                 dataset_name),
                                   max_edges = if (n > 100) 500 else 400)
      ggsave(paste0(prefix, "_concentric_ring.png"), p,
             width = 9, height = 9, dpi = 300)
    }, error = function(e) cat("        FAILED:", conditionMessage(e), "\n"))

    cat("      Done:", model, "\n")
  }

  cat("  [new_viz] Finished:", dataset_name, "\n")
  invisible(NULL)
}


# =============================================================================
# ENTRY POINT
# =============================================================================
run_all_new_visualizations <- function() {
  run_dir <- Sys.getenv("APP_RUN_DIR",
    unset = "./output/application/raw/application_run_20260411_163055")
  ppc_dir <- "output/application/ppc"

  # Collect dataset names from both sources
  datasets <- character(0)
  if (dir.exists(ppc_dir)) {
    datasets <- list.dirs(ppc_dir, recursive = FALSE, full.names = FALSE)
    datasets <- datasets[!datasets %in% c(".DS_Store")]
  }
  if (dir.exists(run_dir)) {
    run_fits <- list.files(run_dir, pattern = "_(WST|SST|DCSBM)_fit\\.rds$")
    run_ds   <- unique(sub("_(WST|SST|DCSBM)_fit\\.rds$", "", run_fits))
    datasets <- unique(c(datasets, run_ds))
  }

  if (length(datasets) == 0) {
    cat("[new_viz] No datasets found. Run application first.\n")
    return(invisible(NULL))
  }

  cat("[new_viz] Datasets:", paste(datasets, collapse = ", "), "\n\n")

  for (ds in datasets) {
    tryCatch(
      generate_new_visualizations(ds, ppc_dir = ppc_dir),
      error = function(e) cat("ERROR for", ds, ":", conditionMessage(e), "\n")
    )
  }

  cat("\n[new_viz] All done.\n")
  cat("Output \u2192 output/application/plots/<dataset>/new_visualizations/\n")
}
