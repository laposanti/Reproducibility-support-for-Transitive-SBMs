#!/usr/bin/env Rscript
# =============================================================================
# Mirrored OCRP prior: comparison with the standard OCRP
# =============================================================================

cat("================================================================\n")
cat("  Mirrored OCRP Prior — Comparison with Standard OCRP\n")
cat("================================================================\n\n")

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")
source("helper_folder/helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("core/my_best_try_so_far.R")

out_dir <- "output/simulation/ocrp_tests_mirrored"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ocrp_mirrored_log_weights_packed <- function(v_minus, theta_ocrp) {
  v_minus <- as.numeric(v_minus)
  K <- length(v_minus)
  if (K < 1L) stop("ocrp_mirrored_log_weights_packed: K<1.", call. = FALSE)
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

log_ocrp_mirrored_composition <- function(n_k, theta) {
  K <- length(n_k)
  n <- sum(n_k)
  if (K == 0L || n == 0L) return(-Inf)
  S_j <- cumsum(n_k)
  K * log(theta) - sum(log(theta + 0:(n - 1L))) + lfactorial(n) - sum(log(S_j))
}

log_ocrp_standard_composition <- function(n_k, theta) {
  K <- length(n_k)
  n <- sum(n_k)
  if (K == 0L || n == 0L) return(-Inf)
  N_j <- rev(cumsum(rev(n_k)))
  K * log(theta) - sum(log(theta + 0:(n - 1L))) + lfactorial(n) - sum(log(N_j))
}

simulate_ocrp_variant <- function(n, theta, prior = c("standard", "mirrored"), seed = NULL) {
  prior <- match.arg(prior)
  if (!is.null(seed)) set.seed(seed)

  blocks <- 1L

  for (i in 2:n) {
    K <- length(blocks)
    pw <- if (prior == "standard") {
      ocrp_log_weights_packed(blocks, theta_ocrp = theta)
    } else {
      ocrp_mirrored_log_weights_packed(blocks, theta_ocrp = theta)
    }
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

  blocks
}

cat("Checking reversal identity of the predictive weights...\n")

theta_check <- c(0.5, 1, 2, 5)
configs_check <- list(
  c(5, 3, 2),
  c(10, 10, 10, 10),
  c(20, 15, 10, 5),
  c(4, 1, 7, 2, 3)
)

formula_checks <- data.frame(
  theta = numeric(),
  config = character(),
  max_exist_err = numeric(),
  max_new_err = numeric(),
  comp_err = numeric(),
  stringsAsFactors = FALSE
)

for (theta in theta_check) {
  for (v in configs_check) {
    K <- length(v)
    std_rev <- ocrp_log_weights_packed(rev(v), theta_ocrp = theta)
    mir <- ocrp_mirrored_log_weights_packed(v, theta_ocrp = theta)

    exist_target <- rev(exp(std_rev$exist))
    new_target <- rev(exp(std_rev$new))

    max_exist_err <- max(abs(exp(mir$exist) - exist_target))
    max_new_err <- max(abs(exp(mir$new) - new_target))
    comp_err <- abs(
      log_ocrp_mirrored_composition(v, theta) -
        log_ocrp_standard_composition(rev(v), theta)
    )

    formula_checks <- rbind(
      formula_checks,
      data.frame(
        theta = theta,
        config = paste(v, collapse = ","),
        max_exist_err = max_exist_err,
        max_new_err = max_new_err,
        comp_err = comp_err,
        stringsAsFactors = FALSE
      )
    )
  }
}

write.csv(formula_checks, file.path(out_dir, "formula_checks.csv"), row.names = FALSE)
print(formula_checks)

theta_sweep <- c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20)
n_customers <- 100
n_sim <- 2000

cat("\nSimulating standard and mirrored priors...\n")

sim_stats <- data.frame(
  prior = character(),
  theta = numeric(),
  rep = integer(),
  K = integer(),
  max_block = integer(),
  min_block = integer(),
  mean_block = numeric(),
  entropy = numeric(),
  first_block = integer(),
  last_block = integer(),
  stringsAsFactors = FALSE
)

for (prior in c("standard", "mirrored")) {
  cat(sprintf("  Prior: %s\n", prior))
  for (theta in theta_sweep) {
    for (r in seq_len(n_sim)) {
      blocks <- simulate_ocrp_variant(n_customers, theta, prior = prior)
      K <- length(blocks)
      p_bs <- blocks / sum(blocks)
      entropy <- -sum(p_bs * log(p_bs))

      sim_stats <- rbind(
        sim_stats,
        data.frame(
          prior = prior,
          theta = theta,
          rep = r,
          K = K,
          max_block = max(blocks),
          min_block = min(blocks),
          mean_block = mean(blocks),
          entropy = entropy,
          first_block = blocks[1L],
          last_block = blocks[K],
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

write.csv(sim_stats, file.path(out_dir, "mirror_vs_standard_raw_stats.csv"), row.names = FALSE)

summary_stats <- aggregate(
  cbind(K, max_block, min_block, mean_block, entropy, first_block, last_block) ~ prior + theta,
  data = sim_stats,
  FUN = mean
)
summary_stats$K_sd <- aggregate(K ~ prior + theta, data = sim_stats, FUN = sd)$K
write.csv(summary_stats, file.path(out_dir, "mirror_vs_standard_summary.csv"), row.names = FALSE)
print(summary_stats)

summary_std <- summary_stats[summary_stats$prior == "standard", ]
summary_mir <- summary_stats[summary_stats$prior == "mirrored", ]
summary_compare <- merge(
  summary_std,
  summary_mir,
  by = "theta",
  suffixes = c("_standard", "_mirrored")
)
summary_compare$delta_K <- summary_compare$K_standard - summary_compare$K_mirrored
summary_compare$delta_max_block <- summary_compare$max_block_standard - summary_compare$max_block_mirrored
summary_compare$delta_entropy <- summary_compare$entropy_standard - summary_compare$entropy_mirrored
summary_compare$delta_first_block <- summary_compare$first_block_standard - summary_compare$first_block_mirrored
summary_compare$delta_last_block <- summary_compare$last_block_standard - summary_compare$last_block_mirrored
write.csv(summary_compare, file.path(out_dir, "mirror_vs_standard_comparison.csv"), row.names = FALSE)

selected_theta <- c(0.5, 1, 2, 5)

# -----------------------------------------------------------------------------
# Plot 1: empirical pmf of K
# -----------------------------------------------------------------------------
pdf(file.path(out_dir, "mirror_vs_standard_K_pmf.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (theta in selected_theta) {
  sub_std <- sim_stats[sim_stats$theta == theta & sim_stats$prior == "standard", ]
  sub_mir <- sim_stats[sim_stats$theta == theta & sim_stats$prior == "mirrored", ]
  K_vals <- sort(unique(c(sub_std$K, sub_mir$K)))
  pmf_std <- table(factor(sub_std$K, levels = K_vals)) / nrow(sub_std)
  pmf_mir <- table(factor(sub_mir$K, levels = K_vals)) / nrow(sub_mir)

  plot(K_vals, as.numeric(pmf_std), type = "b", pch = 16, lwd = 2,
       col = "steelblue", ylim = range(c(pmf_std, pmf_mir)),
       xlab = "K", ylab = "Empirical probability",
       main = bquote(theta == .(theta)))
  lines(K_vals, as.numeric(pmf_mir), type = "b", pch = 17, lwd = 2,
        col = "firebrick")
  legend("topright", legend = c("Standard OCRP", "Mirrored prior"),
         col = c("steelblue", "firebrick"), pch = c(16, 17), lwd = 2, cex = 0.8)
}
dev.off()

# -----------------------------------------------------------------------------
# Plot 2: block size by relative position
# -----------------------------------------------------------------------------
position_profiles <- data.frame(
  prior = character(),
  theta = numeric(),
  bin = character(),
  bin_mid = numeric(),
  mean_size = numeric(),
  stringsAsFactors = FALSE
)

breaks <- seq(0, 1, 0.2)
bin_levels <- levels(cut(c(0.1, 0.3, 0.5, 0.7, 0.9), breaks = breaks, include.lowest = TRUE))
bin_midpoints <- c(0.1, 0.3, 0.5, 0.7, 0.9)

for (prior in c("standard", "mirrored")) {
  for (theta in selected_theta) {
    pos_size <- data.frame(rel_pos = numeric(), size = numeric())
    for (r in seq_len(500L)) {
      blocks <- simulate_ocrp_variant(n_customers, theta, prior = prior)
      K <- length(blocks)
      for (j in seq_len(K)) {
        pos_size <- rbind(
          pos_size,
          data.frame(rel_pos = j / K, size = blocks[j])
        )
      }
    }
    pos_size$bin <- cut(pos_size$rel_pos, breaks = breaks, include.lowest = TRUE)
    means <- aggregate(size ~ bin, data = pos_size, FUN = mean)
    mean_vals <- setNames(rep(NA_real_, length(bin_levels)), bin_levels)
    mean_vals[as.character(means$bin)] <- means$size

    position_profiles <- rbind(
      position_profiles,
      data.frame(
        prior = prior,
        theta = theta,
        bin = bin_levels,
        bin_mid = bin_midpoints,
        mean_size = as.numeric(mean_vals),
        stringsAsFactors = FALSE
      )
    )
  }
}

write.csv(position_profiles, file.path(out_dir, "mirror_vs_standard_position_profiles.csv"), row.names = FALSE)

pdf(file.path(out_dir, "mirror_vs_standard_position_profile.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (theta in selected_theta) {
  sub_std <- position_profiles[position_profiles$theta == theta & position_profiles$prior == "standard", ]
  sub_mir <- position_profiles[position_profiles$theta == theta & position_profiles$prior == "mirrored", ]
  yr <- range(c(sub_std$mean_size, sub_mir$mean_size), na.rm = TRUE)

  plot(sub_std$bin_mid, sub_std$mean_size, type = "b", pch = 16, lwd = 2,
       col = "steelblue", ylim = yr, xaxt = "n",
       xlab = "Relative position", ylab = "Mean block size",
       main = bquote(theta == .(theta)))
  axis(1, at = bin_midpoints, labels = bin_levels, las = 2, cex.axis = 0.8)
  lines(sub_mir$bin_mid, sub_mir$mean_size, type = "b", pch = 17, lwd = 2,
        col = "firebrick")
  legend("topright", legend = c("Standard OCRP", "Mirrored prior"),
         col = c("steelblue", "firebrick"), pch = c(16, 17), lwd = 2, cex = 0.8)
}
dev.off()

# Tutorial-style barplot version of the same position profile
pdf(file.path(out_dir, "mirror_vs_standard_position_barplots.pdf"), width = 10, height = 8)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
yr <- range(position_profiles$mean_size, na.rm = TRUE)

for (prior in c("standard", "mirrored")) {
  for (theta in selected_theta) {
    sub <- position_profiles[position_profiles$theta == theta & position_profiles$prior == prior, ]
    barplot(
      sub$mean_size,
      names.arg = sub$bin,
      ylim = yr,
      main = paste(if (prior == "standard") "Standard" else "Mirrored", ", theta =", theta),
      xlab = "Relative position",
      ylab = "Mean block size",
      col = "lightgreen",
      border = "darkgreen",
      las = 2,
      cex.names = 0.8
    )
  }
}
dev.off()

# -----------------------------------------------------------------------------
# Plot 3: exact existing-block probabilities under equal sizes
# -----------------------------------------------------------------------------
v_equal <- rep(10, 5)

pdf(file.path(out_dir, "mirror_vs_standard_equalsize_existing.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (theta in selected_theta) {
  std <- exp(ocrp_log_weights_packed(v_equal, theta_ocrp = theta)$exist)
  mir <- exp(ocrp_mirrored_log_weights_packed(v_equal, theta_ocrp = theta)$exist)
  yr <- range(c(std, mir))

  plot(seq_along(v_equal), std, type = "b", pch = 16, lwd = 2,
       col = "steelblue", ylim = yr,
       xlab = "Block position", ylab = expression(p[j]^old),
       main = bquote(theta == .(theta)))
  lines(seq_along(v_equal), mir, type = "b", pch = 17, lwd = 2,
        col = "firebrick")
  legend("topright", legend = c("Standard OCRP", "Mirrored prior"),
         col = c("steelblue", "firebrick"), pch = c(16, 17), lwd = 2, cex = 0.8)
}
dev.off()

# -----------------------------------------------------------------------------
# Plot 4: order-sensitive end-block summaries across theta
# -----------------------------------------------------------------------------
pdf(file.path(out_dir, "mirror_vs_standard_end_blocks_vs_theta.pdf"), width = 8, height = 5)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))

yr <- range(
  c(
    summary_stats$first_block,
    summary_stats$last_block
  ),
  na.rm = TRUE
)

sub_std <- summary_stats[summary_stats$prior == "standard", ]
sub_mir <- summary_stats[summary_stats$prior == "mirrored", ]

plot(sub_std$theta, sub_std$first_block, type = "b", pch = 16, lwd = 2,
     col = "steelblue", log = "x", ylim = yr,
     xlab = expression(theta), ylab = "Mean block size",
     main = "Mean first and last block size vs theta")
lines(sub_std$theta, sub_std$last_block, type = "b", pch = 1, lwd = 2,
      col = "steelblue", lty = 2)
lines(sub_mir$theta, sub_mir$first_block, type = "b", pch = 17, lwd = 2,
      col = "firebrick")
lines(sub_mir$theta, sub_mir$last_block, type = "b", pch = 2, lwd = 2,
      col = "firebrick", lty = 2)
legend(
  "topright",
  legend = c(
    "Standard: first block",
    "Standard: last block",
    "Mirrored: first block",
    "Mirrored: last block"
  ),
  col = c("steelblue", "steelblue", "firebrick", "firebrick"),
  pch = c(16, 1, 17, 2),
  lty = c(1, 2, 1, 2),
  lwd = 2,
  cex = 0.8
)
dev.off()

# -----------------------------------------------------------------------------
# Plot 5-7: theta-sweeps for invariant summaries
# -----------------------------------------------------------------------------
plot_summary_vs_theta <- function(stat_name, ylab, file_name, main_text) {
  pdf(file.path(out_dir, file_name), width = 8, height = 5)
  par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
  y_std <- sub_std[[stat_name]]
  y_mir <- sub_mir[[stat_name]]
  yr <- range(c(y_std, y_mir), na.rm = TRUE)

  plot(sub_std$theta, y_std, type = "b", pch = 16, lwd = 2,
       col = "steelblue", log = "x", ylim = yr,
       xlab = expression(theta), ylab = ylab, main = main_text)
  lines(sub_mir$theta, y_mir, type = "b", pch = 17, lwd = 2, col = "firebrick")
  legend("topright", legend = c("Standard OCRP", "Mirrored prior"),
         col = c("steelblue", "firebrick"), pch = c(16, 17), lwd = 2, cex = 0.8)
  dev.off()
}

plot_summary_vs_theta("K", "Mean K", "mirror_vs_standard_K_mean_vs_theta.pdf",
                      "Mean number of blocks vs theta")
plot_summary_vs_theta("max_block", "Mean largest block size",
                      "mirror_vs_standard_max_block_vs_theta.pdf",
                      "Mean largest block size vs theta")
plot_summary_vs_theta("entropy", "Mean entropy",
                      "mirror_vs_standard_entropy_vs_theta.pdf",
                      "Mean entropy vs theta")

cat(sprintf("\nOutputs saved to: %s/\n", out_dir))
