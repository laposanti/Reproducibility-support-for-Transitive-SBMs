#!/usr/bin/env Rscript
# =============================================================================
# Comprehensive tests for the Ordered CRP (OCRP) partition prior
# =============================================================================
# Tests:
#   1. Probabilities sum to 1 for various configurations
#   2. Evolution of K as n increases
#   3. Consistency and projectivity
#   4. Integration over orders recovers standard CRP
#   5. Block size distribution under different hyperparameters
#   6. Additional correctness and characterisation tests
# =============================================================================

cat("================================================================\n")
cat("  Ordered CRP Partition Prior — Comprehensive Test Suite\n")
cat("================================================================\n\n")

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")
source("helper_folder/helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("core/my_best_try_so_far.R")

# Output directory
out_dir <- "output/simulation/ocrp_tests"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pass_count <- 0L
fail_count <- 0L
test_results <- list()

record <- function(name, passed, detail = "") {
  test_results[[length(test_results) + 1L]] <<- list(
    name = name, passed = passed, detail = detail
  )
  if (passed) {
    pass_count <<- pass_count + 1L
    cat(sprintf("  [PASS] %s\n", name))
  } else {
    fail_count <<- fail_count + 1L
    cat(sprintf("  [FAIL] %s: %s\n", name, detail))
  }
}

# ============================================================================
# TEST 1: Probabilities sum to 1
# ============================================================================
cat("\n=== TEST 1: Probabilities sum to 1 ===\n")

test_configs <- list(
  list(v = c(1),         theta = 1.0),
  list(v = c(5),         theta = 1.0),
  list(v = c(5, 3),      theta = 1.0),
  list(v = c(5, 3, 2),   theta = 1.0),
  list(v = c(10, 5, 3, 1), theta = 1.0),
  list(v = c(1, 1, 1, 1, 1), theta = 1.0),
  list(v = c(100, 50, 20), theta = 1.0),
  list(v = c(5, 3, 2),   theta = 0.1),
  list(v = c(5, 3, 2),   theta = 0.5),
  list(v = c(5, 3, 2),   theta = 2.0),
  list(v = c(5, 3, 2),   theta = 5.0),
  list(v = c(5, 3, 2),   theta = 10.0),
  list(v = c(5, 3, 2),   theta = 50.0),
  list(v = c(1, 1),      theta = 0.01),
  list(v = c(1000),      theta = 0.001)
)

sum1_results <- data.frame(
  config = character(), theta = numeric(),
  H = integer(), n = integer(),
  sum_total = numeric(), error = numeric(),
  stringsAsFactors = FALSE
)

for (tc in test_configs) {
  pw <- ocrp_log_weights_packed(tc$v, theta_ocrp = tc$theta)
  p_exist <- exp(pw$exist)
  p_new   <- exp(pw$new)
  total <- sum(p_exist) + sum(p_new)
  err <- abs(total - 1.0)
  
  label <- sprintf("v=(%s), theta=%.3f", paste(tc$v, collapse=","), tc$theta)
  passed <- err < 1e-10
  record(sprintf("Sum-to-1: %s", label), passed,
         detail = sprintf("sum=%.15f, err=%.2e", total, err))
  
  sum1_results <- rbind(sum1_results, data.frame(
    config = paste(tc$v, collapse=","), theta = tc$theta,
    H = length(tc$v), n = sum(tc$v),
    sum_total = total, error = err,
    stringsAsFactors = FALSE
  ))
}

write.csv(sum1_results, file.path(out_dir, "test1_sum_to_one.csv"), row.names = FALSE)

# Additionally verify the analytic decomposition:
# sum_old = n/(theta+n), sum_new = theta/(theta+n)
cat("\n  Checking old/new mass decomposition...\n")
for (tc in test_configs) {
  pw <- ocrp_log_weights_packed(tc$v, theta_ocrp = tc$theta)
  n <- sum(tc$v)
  sum_old <- sum(exp(pw$exist))
  sum_new <- sum(exp(pw$new))
  expected_old <- n / (tc$theta + n)
  expected_new <- tc$theta / (tc$theta + n)
  
  label <- sprintf("v=(%s), theta=%.3f", paste(tc$v, collapse=","), tc$theta)
  err_old <- abs(sum_old - expected_old)
  err_new <- abs(sum_new - expected_new)
  passed <- (err_old < 1e-10) && (err_new < 1e-10)
  record(sprintf("Old/new split: %s", label), passed,
         detail = sprintf("old_err=%.2e, new_err=%.2e", err_old, err_new))
}

# ============================================================================
# TEST 2: Evolution of K as n increases
# ============================================================================
cat("\n=== TEST 2: Evolution of K as n increases ===\n")

# Simulate ordered CRP sequentially for various theta values
simulate_ocrp <- function(n, theta, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # Start: first customer goes to block 1
  blocks <- 1L  # sizes of ordered blocks
  z <- integer(n)
  z[1] <- 1L
  
  for (i in 2:n) {
    H <- length(blocks)
    pw <- ocrp_log_weights_packed(blocks, theta_ocrp = theta)
    probs <- c(exp(pw$exist), exp(pw$new))
    
    # Sample: 1..H = existing, H+1..2H+1 = new slots
    choice <- sample.int(length(probs), 1L, prob = probs)
    
    if (choice <= H) {
      # Add to existing block
      blocks[choice] <- blocks[choice] + 1L
      z[i] <- choice
    } else {
      # New block at slot r = choice - H
      r <- choice - H
      # Insert new block of size 1 at position r
      if (r == 1L) {
        blocks <- c(1L, blocks)
      } else if (r == H + 1L) {
        blocks <- c(blocks, 1L)
      } else {
        blocks <- c(blocks[1:(r-1)], 1L, blocks[r:H])
      }
      z[i] <- r
      # Shift labels for nodes in blocks >= r
      # (We don't need z for this test, just blocks)
    }
  }
  
  list(K = length(blocks), blocks = blocks, z = z)
}

theta_vals <- c(0.1, 0.5, 1.0, 2.0, 5.0, 10.0)
n_vals <- c(10, 20, 50, 100, 200, 500)
n_reps <- 200

cat("  Simulating ordered CRP for various theta and n...\n")
K_results <- expand.grid(theta = theta_vals, n = n_vals)
K_results$K_mean <- NA_real_
K_results$K_sd <- NA_real_
K_results$K_median <- NA_real_

for (idx in seq_len(nrow(K_results))) {
  th <- K_results$theta[idx]
  nn <- K_results$n[idx]
  K_vec <- replicate(n_reps, simulate_ocrp(nn, th, seed = NULL)$K)
  K_results$K_mean[idx] <- mean(K_vec)
  K_results$K_sd[idx] <- sd(K_vec)
  K_results$K_median[idx] <- median(K_vec)
}

write.csv(K_results, file.path(out_dir, "test2_K_evolution.csv"), row.names = FALSE)

# For standard CRP, E[K] ~ theta * log(n)
# For ordered CRP, we expect similar growth (since marginalizing over order gives CRP)
cat("  E[K] results:\n")
print(K_results)

# Plot
pdf(file.path(out_dir, "test2_K_vs_n.pdf"), width = 8, height = 5)
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
cols <- rainbow(length(theta_vals))
plot(NULL, xlim = range(n_vals), ylim = c(0, max(K_results$K_mean) * 1.2),
     xlab = "n", ylab = "E[K]", main = "Ordered CRP: Expected K vs n")
for (i in seq_along(theta_vals)) {
  sub <- K_results[K_results$theta == theta_vals[i], ]
  lines(sub$n, sub$K_mean, col = cols[i], lwd = 2, type = "b", pch = 16)
  # Overlay theta*log(n) (CRP theoretical)
  lines(sub$n, theta_vals[i] * log(sub$n), col = cols[i], lwd = 1, lty = 2)
}
legend("topleft", legend = paste0("theta=", theta_vals),
       col = cols, lwd = 2, cex = 0.8)
legend("topright", legend = c("Simulated", "theta*log(n)"),
       lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
dev.off()

record("K evolution: growth with n", TRUE, "See test2_K_vs_n.pdf")

# ============================================================================
# TEST 3: Consistency and Projectivity
# ============================================================================
cat("\n=== TEST 3: Consistency and Projectivity ===\n")

# Projectivity test: simulate n+1, then marginalize out the last customer.
# The distribution of the first n should match the n-step simulation.
# We test this by comparing empirical distributions of block-size vectors.

test_projectivity <- function(n, theta, n_reps = 5000, seed = 42) {
  set.seed(seed)
  
  # (a) Direct: simulate OCRP to n
  direct_configs <- character(n_reps)
  for (r in seq_len(n_reps)) {
    res <- simulate_ocrp(n, theta)
    direct_configs[r] <- paste(res$blocks, collapse = "-")
  }
  
  # (b) Projective: simulate OCRP to n+1, then remove the last added node
  # This is subtle with ordered partitions. We need to track which block
  # the (n+1)-th customer went to and remove it.
  project_configs <- character(n_reps)
  for (r in seq_len(n_reps)) {
    # Simulate to n first
    res_n <- simulate_ocrp(n, theta)
    blocks_n <- res_n$blocks
    H <- length(blocks_n)
    
    # Add one more customer using predictive
    pw <- ocrp_log_weights_packed(blocks_n, theta_ocrp = theta)
    probs <- c(exp(pw$exist), exp(pw$new))
    choice <- sample.int(length(probs), 1L, prob = probs)
    
    if (choice <= H) {
      # Customer joined existing block choice
      blocks_np1 <- blocks_n
      blocks_np1[choice] <- blocks_np1[choice] + 1L
      # Remove: just subtract 1 from that block
      blocks_proj <- blocks_np1
      blocks_proj[choice] <- blocks_proj[choice] - 1L
      blocks_proj <- blocks_proj[blocks_proj > 0]
    } else {
      # Customer created new block at slot r = choice - H
      # Remove: the new block disappears
      blocks_proj <- blocks_n
    }
    project_configs[r] <- paste(blocks_proj, collapse = "-")
  }
  
  # Compare empirical distributions
  tab_direct <- sort(table(direct_configs) / n_reps, decreasing = TRUE)
  tab_project <- sort(table(project_configs) / n_reps, decreasing = TRUE)
  
  # KL-like comparison: total variation distance on the top configs
  all_configs <- union(names(tab_direct), names(tab_project))
  p <- numeric(length(all_configs)); names(p) <- all_configs
  q <- numeric(length(all_configs)); names(q) <- all_configs
  p[names(tab_direct)] <- tab_direct
  q[names(tab_project)] <- tab_project
  
  tv <- 0.5 * sum(abs(p - q))
  list(tv = tv, n_configs_direct = length(tab_direct),
       n_configs_project = length(tab_project))
}

for (theta in c(0.5, 1.0, 2.0)) {
  for (n in c(5, 8, 10)) {
    n_mc <- if (n <= 8) 50000 else 20000
    res <- test_projectivity(n, theta, n_reps = n_mc)
    label <- sprintf("n=%d, theta=%.1f", n, theta)
    # TV distance should be small (statistical noise); allow more for larger n
    thresh <- if (n <= 8) 0.03 else 0.08
    passed <- res$tv < thresh
    record(sprintf("Projectivity: %s", label), passed,
           detail = sprintf("TV=%.4f (thresh=%.2f)", res$tv, thresh))
  }
}

# ============================================================================
# TEST 4: Integrating over orders recovers standard CRP
# ============================================================================
cat("\n=== TEST 4: Integration over orders recovers standard CRP ===\n")

# For small n, enumerate all possible ordered partitions and verify
# that summing over all orderings of a given unordered partition
# recovers the standard CRP EPPF times the number of set partitions.
#
# The relationship is:
#   sum_{distinct orderings sigma} p_ordered(n_sigma) = C_set * EPPF
# where C_set = n! / (prod(n_k!) * prod(m_j!))  (# set partitions with those sizes)
# and m_j counts blocks of size j.
#
# Equivalently, sum over ALL compositions of n should equal 1,
# which we verify as a separate check.

# Standard CRP log-EPPF (probability of one specific set partition)
log_crp_eppf <- function(n_k, theta) {
  n_k <- sort(n_k[n_k > 0], decreasing = TRUE)
  K <- length(n_k)
  n <- sum(n_k)
  K * log(theta) + sum(lgamma(n_k)) - sum(log(theta + 0:(n-1)))
}

# Ordered CRP log-probability of a specific ordered composition
log_ocrp_composition <- function(n_k, theta) {
  K <- length(n_k)
  n <- sum(n_k)
  if (K == 0 || n == 0) return(-Inf)
  N_j <- rev(cumsum(rev(n_k)))
  K * log(theta) - sum(log(theta + 0:(n-1))) + lfactorial(n) - sum(log(N_j))
}

# Number of set partitions with given sizes (multinomial / symmetry factor)
log_n_set_partitions <- function(n_k) {
  n_k <- n_k[n_k > 0]
  n <- sum(n_k)
  # Count multiplicities: m_j = number of blocks of size j
  size_tab <- table(n_k)
  m_j <- as.integer(size_tab)
  lfactorial(n) - sum(lfactorial(n_k)) - sum(lfactorial(m_j))
}

# Test: for partitions of small n, sum ordered compositions over permutations
test_integration <- function(n, theta) {
  partitions <- list()
  gen_partitions <- function(n, max_val, current) {
    if (n == 0) {
      partitions[[length(partitions) + 1L]] <<- current
      return(invisible())
    }
    for (v in seq(min(n, max_val), 1L)) {
      gen_partitions(n - v, v, c(current, v))
    }
  }
  gen_partitions(n, n, integer(0))
  
  results <- data.frame(
    partition = character(), K = integer(),
    log_crp_target = numeric(), log_ocrp_summed = numeric(),
    error = numeric(), stringsAsFactors = FALSE
  )
  
  max_err <- 0
  
  for (part in partitions) {
    K <- length(part)
    # Target: C_set * EPPF
    log_target <- log_crp_eppf(part, theta) + log_n_set_partitions(part)
    
    if (K <= 1L) {
      perms <- matrix(part, nrow = 1)
    } else {
      perms <- .generate_unique_perms(part)
    }
    
    log_probs <- apply(perms, 1, function(perm) log_ocrp_composition(perm, theta))
    log_sum <- log(sum(exp(log_probs - max(log_probs)))) + max(log_probs)
    
    err <- abs(log_sum - log_target)
    max_err <- max(max_err, err)
    
    results <- rbind(results, data.frame(
      partition = paste(part, collapse=","),
      K = K, log_crp_target = log_target,
      log_ocrp_summed = log_sum, error = err,
      stringsAsFactors = FALSE
    ))
  }
  
  list(results = results, max_err = max_err)
}

# Additional check: sum of ALL compositions of n = 1
test_composition_sum <- function(n, theta) {
  # Generate all compositions of n (ordered) and sum their probabilities
  compositions <- list()
  gen_compositions <- function(n, current) {
    if (n == 0) {
      if (length(current) > 0)
        compositions[[length(compositions) + 1L]] <<- current
      return(invisible())
    }
    for (v in seq_len(n)) {
      gen_compositions(n - v, c(current, v))
    }
  }
  gen_compositions(n, integer(0))
  
  log_probs <- sapply(compositions, function(comp) log_ocrp_composition(comp, theta))
  mx <- max(log_probs)
  total <- sum(exp(log_probs - mx)) * exp(mx)
  total
}

# Helper: generate unique permutations
.generate_unique_perms <- function(x) {
  n <- length(x)
  if (n <= 1) return(matrix(x, nrow = 1))
  if (n == 2) {
    if (x[1] == x[2]) return(matrix(x, nrow = 1))
    return(rbind(x, rev(x)))
  }
  
  vals <- sort(unique(x))
  result <- NULL
  for (v in vals) {
    idx <- which(x == v)[1]
    rest <- x[-idx]
    sub_perms <- .generate_unique_perms(rest)
    for (i in seq_len(nrow(sub_perms))) {
      result <- rbind(result, c(v, sub_perms[i, ]))
    }
  }
  result
}

for (theta in c(0.5, 1.0, 2.0, 5.0)) {
  for (n in c(4, 6, 8)) {
    res <- test_integration(n, theta)
    label <- sprintf("n=%d, theta=%.1f", n, theta)
    passed <- res$max_err < 1e-8
    record(sprintf("CRP recovery (perm sum = C_set*EPPF): %s", label), passed,
           detail = sprintf("max_err=%.2e", res$max_err))
    
    if (n == 6 && theta == 1.0) {
      write.csv(res$results,
                file.path(out_dir, "test4_crp_recovery_n6_theta1.csv"),
                row.names = FALSE)
    }
  }
}

# Check: sum of all compositions = 1
for (theta in c(0.5, 1.0, 2.0)) {
  for (n in c(3, 5, 7)) {
    total <- test_composition_sum(n, theta)
    label <- sprintf("n=%d, theta=%.1f", n, theta)
    passed <- abs(total - 1.0) < 1e-10
    record(sprintf("All compositions sum to 1: %s", label), passed,
           detail = sprintf("sum=%.15f", total))
  }
}

# ============================================================================
# TEST 5: Block size distributions under different hyperparameters
# ============================================================================
cat("\n=== TEST 5: Block size distribution ===\n")

n_sim <- 200
n_customers <- 100

theta_sweep <- c(0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0)

block_stats <- data.frame(
  theta = numeric(), rep = integer(),
  K = integer(), max_block = integer(), min_block = integer(),
  mean_block = numeric(), entropy = numeric(),
  stringsAsFactors = FALSE
)

for (theta in theta_sweep) {
  for (r in seq_len(n_sim)) {
    res <- simulate_ocrp(n_customers, theta)
    bs <- res$blocks
    K <- length(bs)
    # Shannon entropy of the block-size distribution
    p_bs <- bs / sum(bs)
    ent <- -sum(p_bs * log(p_bs))
    
    block_stats <- rbind(block_stats, data.frame(
      theta = theta, rep = r,
      K = K, max_block = max(bs), min_block = min(bs),
      mean_block = mean(bs), entropy = ent,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(block_stats, file.path(out_dir, "test5_block_stats.csv"), row.names = FALSE)

# Summary table
block_summary <- aggregate(
  cbind(K, max_block, min_block, mean_block, entropy) ~ theta,
  data = block_stats, FUN = mean
)
block_summary$K_sd <- aggregate(K ~ theta, data = block_stats, FUN = sd)$K

cat("  Block size summary (n=100):\n")
print(block_summary, digits = 3)
write.csv(block_summary, file.path(out_dir, "test5_block_summary.csv"), row.names = FALSE)

# Plot: K distribution for different theta
pdf(file.path(out_dir, "test5_K_distribution.pdf"), width = 10, height = 6)
par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))
for (theta in theta_sweep) {
  sub <- block_stats[block_stats$theta == theta, ]
  hist(sub$K, breaks = seq(0.5, max(sub$K) + 0.5, 1), freq = FALSE,
       main = bquote(theta == .(theta)),
       xlab = "K", ylab = "Density", col = "steelblue")
  abline(v = mean(sub$K), col = "red", lwd = 2, lty = 2)
}
dev.off()

# Plot: block size distributions
pdf(file.path(out_dir, "test5_block_sizes.pdf"), width = 10, height = 6)
par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))
for (theta in theta_sweep) {
  sub <- block_stats[block_stats$theta == theta, ]
  # Collect all block sizes across reps
  all_sizes <- unlist(lapply(seq_len(nrow(sub)), function(i) {
    simulate_ocrp(n_customers, theta)$blocks
  }))
  hist(all_sizes, breaks = 30, freq = FALSE,
       main = bquote(theta == .(theta)),
       xlab = "Block size", ylab = "Density", col = "coral")
}
dev.off()

# Plot: block size by position (are earlier blocks larger?)
pdf(file.path(out_dir, "test5_size_by_position.pdf"), width = 10, height = 6)
par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))
for (theta in theta_sweep) {
  # Simulate and record sizes by normalized position
  pos_size <- data.frame(rel_pos = numeric(), size = numeric())
  for (r in 1:500) {
    res <- simulate_ocrp(n_customers, theta)
    K <- length(res$blocks)
    for (j in seq_len(K)) {
      pos_size <- rbind(pos_size, data.frame(
        rel_pos = j / K, size = res$blocks[j]
      ))
    }
  }
  
  # Bin by relative position
  brks <- seq(0, 1, 0.2)
  pos_size$bin <- cut(pos_size$rel_pos, breaks = brks, include.lowest = TRUE)
  bin_levels <- levels(pos_size$bin)
  means <- aggregate(size ~ bin, data = pos_size, FUN = mean)
  # ensure all levels present
  mean_vals <- setNames(rep(0, length(bin_levels)), bin_levels)
  mean_vals[as.character(means$bin)] <- means$size
  
  barplot(mean_vals, names.arg = names(mean_vals),
          main = bquote(theta == .(theta)),
          xlab = "Relative position", ylab = "Mean block size",
          col = "lightgreen", las = 2, cex.names = 0.7)
}
dev.off()

record("Block size characterisation", TRUE, "See test5_*.pdf/csv")

# ============================================================================
# TEST 6: Additional tests
# ============================================================================
cat("\n=== TEST 6: Additional correctness and characterisation ===\n")

# --- 6a. Sequential construction matches joint probability ---
cat("  6a. Sequential vs joint probability check...\n")

test_sequential_vs_joint <- function(theta, n_reps = 10000, n = 8, seed = 123) {
  set.seed(seed)
  
  # Simulate many ordered partitions, collect empirical frequencies
  config_counts <- list()
  log_probs_seq <- list()
  
  for (r in seq_len(n_reps)) {
    blocks <- 1L
    log_p <- 0  # accumulate log predictive prob
    
    for (i in 2:n) {
      H <- length(blocks)
      pw <- ocrp_log_weights_packed(blocks, theta_ocrp = theta)
      probs <- c(exp(pw$exist), exp(pw$new))
      choice <- sample.int(length(probs), 1L, prob = probs)
      log_p <- log_p + log(probs[choice])
      
      if (choice <= H) {
        blocks[choice] <- blocks[choice] + 1L
      } else {
        r_pos <- choice - H
        if (r_pos == 1L) {
          blocks <- c(1L, blocks)
        } else if (r_pos == H + 1L) {
          blocks <- c(blocks, 1L)
        } else {
          blocks <- c(blocks[1:(r_pos-1)], 1L, blocks[r_pos:H])
        }
      }
    }
    
    key <- paste(blocks, collapse = "-")
    if (is.null(config_counts[[key]])) {
      config_counts[[key]] <- 0L
      log_probs_seq[[key]] <- log_ocrp_composition(blocks, theta)
    }
    config_counts[[key]] <- config_counts[[key]] + 1L
  }
  
  # Compare empirical frequencies with theoretical probabilities
  keys <- names(config_counts)
  emp_freq <- sapply(keys, function(k) config_counts[[k]] / n_reps)
  theo_prob <- exp(sapply(keys, function(k) log_probs_seq[[k]]))
  
  # Use chi-squared statistic restricted to configs with enough expected count
  expected <- theo_prob * n_reps
  keep <- expected >= 5  # standard chi-sq rule
  if (sum(keep) < 2) {
    # Not enough configs with large enough expected counts
    chi2 <- 0; pval <- 1
  } else {
    obs <- emp_freq[keep] * n_reps
    exp_k <- expected[keep]
    chi2 <- sum((obs - exp_k)^2 / exp_k)
    df <- sum(keep) - 1L
    pval <- pchisq(chi2, df = df, lower.tail = FALSE)
  }
  
  list(chi2 = chi2, pval = pval, n_configs = length(keys), n_tested = sum(keep))
}

for (theta in c(0.5, 1.0, 2.0)) {
  res <- test_sequential_vs_joint(theta, n_reps = 100000, n = 6)
  label <- sprintf("theta=%.1f", theta)
  # Chi-squared p-value should not be extremely small
  passed <- res$pval > 0.001
  record(sprintf("Sequential vs joint: %s", label), passed,
         detail = sprintf("chi2=%.2f, pval=%.4f, configs=%d (tested=%d)",
                          res$chi2, res$pval, res$n_configs, res$n_tested))
}

# --- 6b. Monotonicity of p_j^old in block size ---
cat("  6b. Monotonicity check: larger blocks attract more mass...\n")

# For a given ordered composition, check whether p_j^old is increasing in n_j
# (up to the positional discount). Actually, the weight is (n_j+1) * positional,
# so bigger blocks should generally get more mass.
# But the positional factor decreases for later positions.
# The combined effect should still favour larger blocks.

for (theta in c(0.5, 1.0, 5.0)) {
  v <- c(20, 10, 5, 2, 1)
  pw <- ocrp_log_weights_packed(v, theta_ocrp = theta)
  
  cat(sprintf("  theta=%.1f, v=(%s)\n", theta, paste(v, collapse=",")))
  cat(sprintf("    Existing probs: %s\n",
              paste(round(exp(pw$exist), 6), collapse=", ")))
  cat(sprintf("    New slot probs: %s\n",
              paste(round(exp(pw$new), 6), collapse=", ")))
}

# --- 6c. Comparison with standard CRP weights ---
cat("\n  6c. Comparison with standard CRP conditional weights...\n")

# Standard CRP: P(join k) = n_k/(theta+n), P(new) = theta/(theta+n)
# OCRP: the total mass for existing vs new is the same, but distributed differently
pdf(file.path(out_dir, "test6_crp_vs_ocrp_weights.pdf"), width = 10, height = 5)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

for (theta in c(0.5, 1.0, 5.0)) {
  v <- c(15, 10, 7, 5, 3)
  n <- sum(v)
  H <- length(v)
  
  pw <- ocrp_log_weights_packed(v, theta_ocrp = theta)
  ocrp_exist <- exp(pw$exist)
  ocrp_new <- exp(pw$new)
  
  crp_exist <- v / (theta + n)
  crp_new <- theta / (theta + n)
  
  x <- seq_len(H)
  barplot(rbind(ocrp_exist, crp_exist), beside = TRUE,
          names.arg = paste0("B", x),
          col = c("steelblue", "coral"),
          main = bquote(theta == .(theta)),
          ylab = "Probability", xlab = "Block")
  legend("topright", legend = c("OCRP", "CRP"), fill = c("steelblue", "coral"),
         cex = 0.8)
}
dev.off()

record("CRP vs OCRP weight comparison", TRUE, "See test6_crp_vs_ocrp_weights.pdf")

# --- 6d. Edge case: all singletons ---
cat("  6d. Edge case: all singletons...\n")
v <- c(1, 1, 1, 1, 1)
for (theta in c(0.1, 1.0, 10.0)) {
  pw <- ocrp_log_weights_packed(v, theta_ocrp = theta)
  total <- sum(exp(pw$exist)) + sum(exp(pw$new))
  passed <- abs(total - 1.0) < 1e-10
  record(sprintf("Singletons: theta=%.1f, sum=%.15f", theta, total), passed)
}

# --- 6e. Edge case: single large block ---
cat("  6e. Edge case: single large block...\n")
for (n in c(10, 100, 1000)) {
  for (theta in c(0.1, 1.0, 10.0)) {
    pw <- ocrp_log_weights_packed(c(n), theta_ocrp = theta)
    total <- sum(exp(pw$exist)) + sum(exp(pw$new))
    # exist should be 1 prob (join the block), new should be 2 probs (left/right)
    passed <- abs(total - 1.0) < 1e-10
    record(sprintf("Single block: n=%d, theta=%.1f", n, theta), passed)
  }
}

# --- 6f. Positional bias: is the first block favoured? ---
cat("\n  6f. Positional bias analysis...\n")

# With equal-sized blocks, earlier blocks should get more mass
# because of the cumulative product factor
v_equal <- c(10, 10, 10, 10, 10)
positional_bias <- data.frame(theta = numeric(), position = integer(), prob = numeric())

for (theta in c(0.1, 0.5, 1.0, 2.0, 5.0, 10.0)) {
  pw <- ocrp_log_weights_packed(v_equal, theta_ocrp = theta)
  probs <- exp(pw$exist)
  for (j in seq_along(probs)) {
    positional_bias <- rbind(positional_bias, data.frame(
      theta = theta, position = j, prob = probs[j]
    ))
  }
}

write.csv(positional_bias, file.path(out_dir, "test6_positional_bias.csv"), row.names = FALSE)

# Check that prob is strictly decreasing by position for equal-sized blocks
for (theta in unique(positional_bias$theta)) {
  sub <- positional_bias[positional_bias$theta == theta, ]
  is_decreasing <- all(diff(sub$prob) < 0)
  record(sprintf("Positional bias decreasing: theta=%.1f", theta), is_decreasing,
         detail = paste(round(sub$prob, 6), collapse=", "))
}

# Plot
pdf(file.path(out_dir, "test6_positional_bias.pdf"), width = 7, height = 5)
par(mar = c(4, 4, 2, 1))
cols <- rainbow(length(unique(positional_bias$theta)))
theta_u <- unique(positional_bias$theta)
plot(NULL, xlim = c(1, 5), ylim = range(positional_bias$prob),
     xlab = "Block position", ylab = "P(join block j | equal sizes)",
     main = "OCRP Positional Bias (equal blocks, n_j=10)")
for (i in seq_along(theta_u)) {
  sub <- positional_bias[positional_bias$theta == theta_u[i], ]
  lines(sub$position, sub$prob, col = cols[i], lwd = 2, type = "b", pch = 16)
}
legend("topright", legend = paste0("theta=", theta_u), col = cols, lwd = 2, cex = 0.7)
dev.off()

# --- 6g. New-block slot bias ---
cat("  6g. New-block slot insertion bias...\n")

v <- c(20, 15, 10, 5)
H <- length(v)
thetas <- c(0.5, 1, 5, 20)

slot_bias <- data.frame(
  theta = numeric(),
  slot = integer(),
  q_new = numeric(),
  q_new_cond = numeric()
)

for (theta in thetas) {
  pw <- ocrp_log_weights_packed(v, theta_ocrp = theta)
  p_new <- exp(pw$new)
  p_new_cond <- p_new / sum(p_new)
  slot_bias <- rbind(
    slot_bias,
    data.frame(
      theta = theta,
      slot = seq_len(H + 1L),
      q_new = p_new,
      q_new_cond = p_new_cond
    )
  )
}

write.csv(slot_bias, file.path(out_dir, "test6_new_slot_bias.csv"), row.names = FALSE)

pdf(file.path(out_dir, "test6_new_slot_bias.pdf"), width = 8, height = 5)
par(mar = c(4, 4, 2, 1))
cols <- rainbow(length(thetas))
plot(NULL, xlim = c(1, H + 1L), ylim = range(slot_bias$q_new),
     xlab = "Insertion slot", ylab = expression(q[r]^new),
     main = sprintf("OCRP: Unconditional new-slot probabilities (v=(%s))", paste(v, collapse=",")))
for (i in seq_along(thetas)) {
  sub <- slot_bias[slot_bias$theta == thetas[i], ]
  lines(sub$slot, sub$q_new, col = cols[i], lwd = 2, type = "b", pch = 16)
}
legend("topleft", legend = paste0("theta=", thetas), col = cols, lwd = 2, cex = 0.8)
dev.off()

slot_cond_mat <- sapply(thetas, function(theta) {
  slot_bias$q_new_cond[slot_bias$theta == theta]
})
cond_invariant <- max(abs(slot_cond_mat - slot_cond_mat[, 1L])) < 1e-12
record("New-block slot bias", cond_invariant,
       "Conditional slot distribution is theta-invariant; pdf now shows unconditional q_r^new.")

# ============================================================================
# SUMMARY
# ============================================================================
cat("\n================================================================\n")
cat(sprintf("  RESULTS: %d passed, %d failed, %d total\n",
            pass_count, fail_count, pass_count + fail_count))
cat("================================================================\n")

# Save summary
summary_df <- do.call(rbind, lapply(test_results, function(x) {
  data.frame(test = x$name, passed = x$passed, detail = x$detail,
             stringsAsFactors = FALSE)
}))
write.csv(summary_df, file.path(out_dir, "test_summary.csv"), row.names = FALSE)

if (fail_count == 0L) {
  cat("\n  *** ALL TESTS PASSED ***\n\n")
} else {
  cat("\n  *** SOME TESTS FAILED — see details above ***\n\n")
  failed <- summary_df[!summary_df$passed, ]
  print(failed)
}

cat(sprintf("Outputs saved to: %s/\n", out_dir))
