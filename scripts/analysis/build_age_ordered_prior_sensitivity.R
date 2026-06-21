#!/usr/bin/env Rscript
# Rebuild Supplement Figure 3: age-ordered prior sensitivity at n = 100.
#
# The figure reports:
# 1. the mean number of occupied blocks under the age-ordered prior across a
#    small theta grid, alongside the ordinary CRP benchmark with alpha = theta
# 2. the mean first-block and last-block sizes under the same age-ordered prior
#
# Outputs are written under output/diagnostics/age_ordered_prior/ and, when the
# manuscript asset folder exists, copied into output/all_figures/ so the local
# TeX bundle can pick up the regenerated figure directly.

ocrp_log_weights_packed_local <- function(v_minus, theta_ocrp) {
  v_minus <- as.numeric(v_minus)
  K <- length(v_minus)
  if (K < 1L) stop("ocrp_log_weights_packed_local: K<1.", call. = FALSE)
  if (!is.finite(theta_ocrp) || theta_ocrp <= 0) {
    stop("theta_ocrp must be a positive finite value.", call. = FALSE)
  }

  n <- sum(v_minus)
  S <- cumsum(v_minus)
  S_prev <- c(0, S)

  # B_j = prod_{i=j}^K S_i / (S_i + 1), with B_{K+1} = 1.
  log_B <- numeric(K + 1L)
  for (j in K:1L) {
    log_B[j] <- log_B[j + 1L] + log(S[j]) - log(S[j] + 1)
  }

  lp_exist <- log(v_minus + 1) - log(theta_ocrp + n) + log_B[1:K]

  lp_new <- numeric(K + 1L)
  for (r in seq_len(K + 1L)) {
    lp_new[r] <- log(theta_ocrp) - log(theta_ocrp + n) -
      log(S_prev[r] + 1) + log_B[r]
  }

  list(exist = lp_exist, new = lp_new)
}

simulate_age_ordered_partition <- function(n, theta_ocrp) {
  if (n < 1L) stop("n must be >= 1.", call. = FALSE)

  block_sizes <- 1L
  if (n == 1L) return(block_sizes)

  for (customer in 2:n) {
    weights <- ocrp_log_weights_packed_local(block_sizes, theta_ocrp = theta_ocrp)
    probs <- c(exp(weights$exist), exp(weights$new))
    choice <- sample.int(length(probs), size = 1L, prob = probs)
    K <- length(block_sizes)

    if (choice <= K) {
      block_sizes[choice] <- block_sizes[choice] + 1L
      next
    }

    slot <- choice - K
    if (slot == 1L) {
      block_sizes <- c(1L, block_sizes)
    } else if (slot == K + 1L) {
      block_sizes <- c(block_sizes, 1L)
    } else {
      block_sizes <- c(block_sizes[1:(slot - 1L)], 1L, block_sizes[slot:K])
    }
  }

  block_sizes
}

crp_expected_k <- function(n, alpha) {
  sum(alpha / (alpha + 0:(n - 1L)))
}

estimate_prior_summaries <- function(theta_grid, n, n_draws, seed) {
  set.seed(seed)

  k_summary <- data.frame(
    theta = theta_grid,
    age_ordered_mean_k = NA_real_,
    crp_mean_k = vapply(theta_grid, function(theta) crp_expected_k(n, theta), numeric(1)),
    stringsAsFactors = FALSE
  )
  block_summary <- data.frame(
    theta = theta_grid,
    top_block_mean = NA_real_,
    bottom_block_mean = NA_real_,
    stringsAsFactors = FALSE
  )

  for (idx in seq_along(theta_grid)) {
    theta <- theta_grid[idx]
    draws <- replicate(n_draws, simulate_age_ordered_partition(n, theta), simplify = FALSE)

    k_summary$age_ordered_mean_k[idx] <- mean(vapply(draws, length, integer(1)))
    block_summary$top_block_mean[idx] <- mean(vapply(draws, function(sizes) sizes[1L], numeric(1)))
    block_summary$bottom_block_mean[idx] <- mean(vapply(draws, function(sizes) sizes[length(sizes)], numeric(1)))
  }

  list(k_summary = k_summary, block_summary = block_summary)
}

draw_age_ordered_prior_sensitivity <- function(k_summary, block_summary, n, n_draws) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.5, 0.8), oma = c(0, 0, 0, 0))
  x_pos <- seq_len(nrow(k_summary))

  plot(
    x_pos,
    k_summary$age_ordered_mean_k,
    type = "o",
    pch = 16,
    lwd = 2,
    col = "#2C7FB8",
    xaxt = "n",
    xlab = expression(vartheta),
    ylab = expression("Mean " * K[n]),
    main = bquote("Mean " * K[n] * " vs " * vartheta),
    ylim = range(c(k_summary$age_ordered_mean_k, k_summary$crp_mean_k))
  )
  axis(1, at = x_pos, labels = format(k_summary$theta, trim = TRUE))
  lines(x_pos, k_summary$crp_mean_k, type = "o", pch = 1, lwd = 1.6, lty = 2, col = "grey45")
  legend(
    "topleft",
    legend = c("Age-ordered prior", expression(alpha == vartheta)),
    col = c("#2C7FB8", "grey45"),
    lty = c(1, 2),
    pch = c(16, 1),
    bty = "n",
    cex = 0.9
  )

  plot(
    x_pos,
    block_summary$top_block_mean,
    type = "o",
    pch = 16,
    lwd = 2,
    col = "#2C7FB8",
    xaxt = "n",
    xlab = expression(vartheta),
    ylab = "Mean block size",
    main = "End-block sizes vs vartheta",
    ylim = range(c(block_summary$top_block_mean, block_summary$bottom_block_mean))
  )
  axis(1, at = x_pos, labels = format(block_summary$theta, trim = TRUE))
  lines(x_pos, block_summary$bottom_block_mean, type = "o", pch = 17, lwd = 2, col = "#D7301F")
  legend(
    "topright",
    legend = c("Top block", "Bottom block"),
    col = c("#2C7FB8", "#D7301F"),
    lty = 1,
    pch = c(16, 17),
    bty = "n",
    cex = 0.9
  )

  mtext(
    sprintf("n = %d, %d prior draws per theta", n, n_draws),
    side = 1,
    outer = FALSE,
    line = 3.1,
    adj = 0.98,
    cex = 0.8,
    col = "grey35"
  )
}

save_age_ordered_prior_sensitivity <- function(
    output_dir = file.path("output", "diagnostics", "age_ordered_prior"),
    theta_grid = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0),
    n = 100L,
    n_draws = 500L,
    seed = 20240621L
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  summaries <- estimate_prior_summaries(
    theta_grid = theta_grid,
    n = as.integer(n),
    n_draws = as.integer(n_draws),
    seed = as.integer(seed)
  )

  utils::write.csv(
    summaries$k_summary,
    file.path(output_dir, "age_ordered_prior_theta_sensitivity_k_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    summaries$block_summary,
    file.path(output_dir, "age_ordered_prior_theta_sensitivity_block_summary.csv"),
    row.names = FALSE
  )

  pdf_path <- file.path(output_dir, "age_ordered_prior_theta_sensitivity.pdf")
  png_path <- file.path(output_dir, "age_ordered_prior_theta_sensitivity.png")

  grDevices::pdf(pdf_path, width = 9.2, height = 4.2)
  draw_age_ordered_prior_sensitivity(
    summaries$k_summary,
    summaries$block_summary,
    n = as.integer(n),
    n_draws = as.integer(n_draws)
  )
  grDevices::dev.off()

  grDevices::png(png_path, width = 2000, height = 900, res = 220, bg = "white")
  draw_age_ordered_prior_sensitivity(
    summaries$k_summary,
    summaries$block_summary,
    n = as.integer(n),
    n_draws = as.integer(n_draws)
  )
  grDevices::dev.off()

  manuscript_dir <- file.path("output", "all_figures")
  if (dir.exists(manuscript_dir)) {
    file.copy(pdf_path, file.path(manuscript_dir, basename(pdf_path)), overwrite = TRUE)
  }

  invisible(list(
    pdf = pdf_path,
    png = png_path,
    k_summary = summaries$k_summary,
    block_summary = summaries$block_summary
  ))
}

theta_grid <- as.numeric(strsplit(
  Sys.getenv("AGE_ORDERED_THETA_GRID", unset = "0.1,0.2,0.5,1,2,5"),
  ",",
  fixed = TRUE
)[[1L]])
theta_grid <- theta_grid[is.finite(theta_grid) & theta_grid > 0]
if (!length(theta_grid)) {
  stop("AGE_ORDERED_THETA_GRID must contain at least one positive numeric value.", call. = FALSE)
}

n <- as.integer(Sys.getenv("AGE_ORDERED_N", unset = "100"))
n_draws <- as.integer(Sys.getenv("AGE_ORDERED_N_DRAWS", unset = "500"))
seed <- as.integer(Sys.getenv("AGE_ORDERED_SEED", unset = "20240621"))

cat("============================================================\n")
cat(" Age-ordered prior sensitivity figure\n")
cat(" Output : output/diagnostics/age_ordered_prior/\n")
cat(" n      : ", n, "\n", sep = "")
cat(" draws  : ", n_draws, "\n", sep = "")
cat(" theta  : ", paste(theta_grid, collapse = ", "), "\n", sep = "")
cat("============================================================\n\n")

save_age_ordered_prior_sensitivity(
  output_dir = file.path("output", "diagnostics", "age_ordered_prior"),
  theta_grid = theta_grid,
  n = n,
  n_draws = n_draws,
  seed = seed
)

cat("Wrote:\n")
cat("  output/diagnostics/age_ordered_prior/age_ordered_prior_theta_sensitivity.pdf\n")
cat("  output/diagnostics/age_ordered_prior/age_ordered_prior_theta_sensitivity.png\n")
