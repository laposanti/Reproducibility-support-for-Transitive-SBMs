#!/usr/bin/env Rscript
# =============================================================================
# 4-panel comparison: mirrored OCRP vs standard CRP
#
# Panel (a): E[K] vs theta (both agree — same EPPF)
# Panel (b): Mean block size by relative position, theta = 0.5
# Panel (c): Mean block size by relative position, theta = 1
# Panel (d): Mean block size by relative position, theta = 2
#
# The CRP and mirrored OCRP share the same EPPF, so marginal
# partition statistics (K, block sizes) are identical.  The
# difference lies in the positional structure: the mirrored OCRP
# assigns larger blocks to lower positions (pyramid shape); the
# CRP, being exchangeable, produces a flat positional profile.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]),
                               winslash = "/", mustWork = TRUE)
  setwd(normalizePath(file.path(dirname(script_path), ".."),
                      winslash = "/", mustWork = TRUE))
}
source("helper_folder/helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("core/my_best_try_so_far.R")

out_dir <- "output/simulation/ocrp_tests_mirrored"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2024)

# --- Mirrored OCRP sampler ---
ocrp_mirrored_log_weights_packed <- function(v_minus, theta_ocrp) {
  v_minus <- as.numeric(v_minus)
  K <- length(v_minus)
  if (K < 1L) stop("K<1", call. = FALSE)
  if (theta_ocrp <= 0) stop("theta_ocrp must be > 0.", call. = FALSE)

  n <- sum(v_minus)
  S <- cumsum(v_minus)
  S_prev <- c(0, S)

  log_B <- numeric(K + 1L)
  if (K >= 1L) {
    for (j in K:1) {
      log_B[j] <- log_B[j + 1L] + log(S[j]) - log(S[j] + 1)
    }
  }

  lp_exist <- log(v_minus + 1) - log(theta_ocrp + n) + log_B[1:K]
  lp_new <- numeric(K + 1L)
  for (r in seq_len(K + 1L)) {
    lp_new[r] <- log(theta_ocrp) - log(theta_ocrp + n) -
      log(S_prev[r] + 1) + log_B[r]
  }

  list(exist = lp_exist, new = lp_new, per_slot = TRUE)
}

simulate_mirrored_ocrp <- function(n, theta) {
  blocks <- 1L
  for (i in 2:n) {
    K <- length(blocks)
    pw <- ocrp_mirrored_log_weights_packed(blocks, theta_ocrp = theta)
    probs <- c(exp(pw$exist), exp(pw$new))
    choice <- sample.int(length(probs), 1L, prob = probs)
    if (choice <= K) {
      blocks[choice] <- blocks[choice] + 1L
    } else {
      r <- choice - K
      if (r == 1L) {
        blocks <- c(1L, blocks)
      } else if (r == K + 1L) {
        blocks <- c(blocks, 1L)
      } else {
        blocks <- c(blocks[1:(r - 1L)], 1L, blocks[r:K])
      }
    }
  }
  blocks  # natural order — position matters
}

# --- Standard CRP sampler (exchangeable: no meaningful ordering) ---
simulate_crp <- function(n, alpha) {
  blocks <- 1L
  for (i in 2:n) {
    probs <- c(blocks, alpha)
    choice <- sample.int(length(probs), 1L, prob = probs)
    if (choice <= length(blocks)) {
      blocks[choice] <- blocks[choice] + 1L
    } else {
      blocks <- c(blocks, 1L)
    }
  }
  # Randomly permute: CRP partitions are exchangeable, so block
  # labels carry no positional information.
  blocks[sample.int(length(blocks))]
}

# =============================================================================
# Simulation
# =============================================================================

n_customers <- 100
n_sim_K     <- 3000   # for E[K] panel
n_sim_pos   <- 1500   # for positional profiles (costlier to store)

theta_grid_K   <- c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, 10)
selected_theta <- c(0.5, 1, 2)  # for positional barplots

# --- Part 1: E[K] across the full theta grid ---
cat("Part 1: simulating E[K] across theta grid...\n")

K_results <- list()
idx <- 0L
for (theta in theta_grid_K) {
  cat(sprintf("  theta = %.2f\n", theta))
  for (r in seq_len(n_sim_K)) {
    bm <- simulate_mirrored_ocrp(n_customers, theta)
    idx <- idx + 1L
    K_results[[idx]] <- data.frame(prior = "Mirrored OCRP",
                                   theta = theta, K = length(bm))
    bc <- simulate_crp(n_customers, theta)
    idx <- idx + 1L
    K_results[[idx]] <- data.frame(prior = "CRP",
                                   theta = theta, K = length(bc))
  }
}
K_data <- bind_rows(K_results)

K_summary <- K_data %>%
  group_by(prior, theta) %>%
  summarise(mean_K = mean(K), .groups = "drop")

analytic_EK <- data.frame(
  theta = theta_grid_K,
  EK = sapply(theta_grid_K, function(th) sum(th / (th + 0:(n_customers - 1))))
)

# --- Part 2: positional profiles for selected theta values ---
cat("Part 2: simulating positional profiles...\n")

breaks <- seq(0, 1, by = 0.2)
bin_labels <- c("[0, 0.2]", "(0.2, 0.4]", "(0.4, 0.6]",
                "(0.6, 0.8]", "(0.8, 1]")
bin_mids <- c(0.1, 0.3, 0.5, 0.7, 0.9)

pos_results <- list()
pidx <- 0L

for (theta in selected_theta) {
  cat(sprintf("  theta = %.2f\n", theta))
  for (prior_name in c("CRP", "Mirrored OCRP")) {
    for (r in seq_len(n_sim_pos)) {
      blocks <- if (prior_name == "CRP") {
        simulate_crp(n_customers, theta)
      } else {
        simulate_mirrored_ocrp(n_customers, theta)
      }
      K <- length(blocks)
      rel_pos <- seq_len(K) / K
      bin_idx <- findInterval(rel_pos, breaks, rightmost.closed = TRUE)
      bin_idx <- pmax(1L, pmin(bin_idx, length(bin_labels)))
      for (j in seq_len(K)) {
        pidx <- pidx + 1L
        pos_results[[pidx]] <- data.frame(
          prior = prior_name,
          theta = theta,
          bin = bin_labels[bin_idx[j]],
          size = blocks[j],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
pos_data <- bind_rows(pos_results)
pos_data$bin <- factor(pos_data$bin, levels = bin_labels)

pos_summary <- pos_data %>%
  group_by(prior, theta, bin) %>%
  summarise(mean_size = mean(size), .groups = "drop")

# =============================================================================
# Style
# =============================================================================

prior_colors <- c(
  "Mirrored OCRP" = "#B2182B",
  "CRP"           = "#2166AC"
)
prior_fills <- prior_colors

prior_shapes <- c("Mirrored OCRP" = 17, "CRP" = 16)

base_theme <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# Panel (a): E[K] vs theta
# =============================================================================

pa <- ggplot(K_summary, aes(x = theta, y = mean_K,
                             colour = prior, shape = prior)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_line(data = analytic_EK, aes(x = theta, y = EK),
            inherit.aes = FALSE, linetype = "dashed",
            colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = prior_colors) +
  scale_shape_manual(values = prior_shapes) +
  scale_x_continuous(breaks = c(0.5, 1, 2, 5, 10)) +
  labs(x = expression(theta), y = expression(E*"["*K[n]*"]"),
       title = expression(bold("(a)")~"Expected number of blocks")) +
  base_theme

# =============================================================================
# Panels (b)-(d): barplots of mean block size by relative position
# =============================================================================

make_pos_panel <- function(theta_val, label) {
  df <- pos_summary %>% filter(theta == theta_val)
  ggplot(df, aes(x = bin, y = mean_size, fill = prior)) +
    geom_col(position = position_dodge(width = 0.75),
             width = 0.65, alpha = 0.85) +
    scale_fill_manual(values = prior_fills) +
    labs(x = "Relative position  (j / K)",
         y = "Mean block size",
         title = bquote(bold((.(label)))~theta == .(theta_val))) +
    base_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1, size = 9))
}

pb <- make_pos_panel(0.5, "b")
pc <- make_pos_panel(1,   "c")
pd <- make_pos_panel(2,   "d")

# =============================================================================
# Combine
# =============================================================================

combined <- (pa | pb) / (pc | pd) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "crp_vs_mirrored_ocrp_4panel.pdf"),
       combined, width = 10, height = 7.5)

cat("Plot saved to:", file.path(out_dir, "crp_vs_mirrored_ocrp_4panel.pdf"), "\n")
cat("Done.\n")
