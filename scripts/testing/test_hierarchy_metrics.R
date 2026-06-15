#!/usr/bin/env Rscript

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")
source("helper_folder/transitivity_check_helper.R")

expect_close <- function(value, target, tolerance = 1e-12, label = "value") {
  if (is.na(value) || abs(value - target) > tolerance) {
    stop(sprintf("%s: got %.12f, expected %.12f", label, value, target))
  }
}

make_three_block_rho <- function(rho_12, rho_23, rho_13) {
  rho_matrix <- matrix(0.5, 3L, 3L)
  rho_matrix[1L, 2L] <- rho_12
  rho_matrix[2L, 1L] <- 1 - rho_12
  rho_matrix[2L, 3L] <- rho_23
  rho_matrix[3L, 2L] <- 1 - rho_23
  rho_matrix[1L, 3L] <- rho_13
  rho_matrix[3L, 1L] <- 1 - rho_13
  rho_matrix
}

score_three_block_manually <- function(rho_matrix) {
  premise <- rho_matrix[1L, 2L] >= 0.5 && rho_matrix[2L, 3L] >= 0.5
  list(
    thetaW = if (premise) as.numeric(rho_matrix[1L, 3L] >= 0.5) else NA_real_,
    thetaS = if (premise) as.numeric(rho_matrix[1L, 3L] >= max(rho_matrix[1L, 2L], rho_matrix[2L, 3L])) else NA_real_
  )
}

score_block_adjacency <- function(adjacency_matrix, labels) {
  violation <- .violation_stats_from_ranked(adjacency_matrix, labels)
  block <- .empirical_block_rates_one(
    A = adjacency_matrix,
    z_vec = labels,
    alpha = 0.5,
    T_block = NULL,
    method_order = "identity"
  )
  list(violation = violation, block = block)
}

cat("=== hierarchy metric validation ===\n")

if (!"method_order" %in% names(formals(block_diag_rates_ext))) {
  stop("block_diag_rates_ext does not expose method_order")
}

perfect_rho <- make_three_block_rho(0.7, 0.6, 0.8)
perfect_helper <- block_diag_rates_ext(perfect_rho, T_block = NULL, method_order = "identity")
perfect_manual <- score_three_block_manually(perfect_rho)
expect_close(perfect_helper$thetaW, perfect_manual$thetaW, label = "perfect thetaW")
expect_close(perfect_helper$thetaS, perfect_manual$thetaS, label = "perfect thetaS")

wst_not_sst_rho <- make_three_block_rho(0.7, 0.7, 0.6)
wst_not_sst_helper <- block_diag_rates_ext(wst_not_sst_rho, T_block = NULL, method_order = "identity")
wst_not_sst_manual <- score_three_block_manually(wst_not_sst_rho)
expect_close(wst_not_sst_helper$thetaW, wst_not_sst_manual$thetaW, label = "WST thetaW")
expect_close(wst_not_sst_helper$thetaS, wst_not_sst_manual$thetaS, label = "SST thetaS")

permuted_rho <- make_three_block_rho(0.2, 0.9, 0.7)
identity_score <- block_diag_rates_ext(permuted_rho, T_block = NULL, method_order = "identity")
mean_score <- block_diag_rates_ext(permuted_rho, T_block = NULL, method_order = "mean")
if (!is.na(identity_score$thetaW)) {
  stop("identity ordering should have no active premise in the permuted example")
}
expect_close(mean_score$thetaW, 1, label = "mean-order thetaW")
expect_close(mean_score$thetaS, 1, label = "mean-order thetaS")

low_violation_bad_zeta <- matrix(0, 3L, 3L)
low_violation_bad_zeta[1L, 2L] <- 100
low_violation_bad_zeta[2L, 3L] <- 100
low_violation_bad_zeta[3L, 1L] <- 1
low_case <- score_block_adjacency(low_violation_bad_zeta, c(1L, 2L, 3L))
expect_close(low_case$violation$rate, 1 / 201, label = "low violation rate")
expect_close(low_case$block$thetaW, 0, label = "low violation thetaW")
expect_close(low_case$block$thetaS, 0, label = "low violation thetaS")

high_violation_good_zeta <- matrix(0, 4L, 4L)
for (block_from in seq_len(3L)) {
  for (block_to in seq.int(block_from + 1L, 4L)) {
    high_violation_good_zeta[block_from, block_to] <- 10
  }
}
high_violation_good_zeta[3L, 4L] <- 0
high_violation_good_zeta[4L, 3L] <- 1000
high_case <- score_block_adjacency(high_violation_good_zeta, c(1L, 2L, 3L, 4L))
expect_close(high_case$violation$rate, 1000 / 1050, label = "high violation rate")
expect_close(high_case$block$thetaW, 1, label = "high violation thetaW")
expect_close(high_case$block$thetaS, 1, label = "high violation thetaS")
expect_close(high_case$block$coverage, 0.5, label = "high violation coverage")

cat(sprintf("perfect rho: thetaW=%.1f thetaS=%.1f\n", perfect_helper$thetaW, perfect_helper$thetaS))
cat(sprintf("WST not SST: thetaW=%.1f thetaS=%.1f\n", wst_not_sst_helper$thetaW, wst_not_sst_helper$thetaS))
cat(sprintf("permuted labels: identity thetaW=%s, mean thetaW=%.1f thetaS=%.1f\n",
            ifelse(is.na(identity_score$thetaW), "NA", sprintf("%.1f", identity_score$thetaW)),
            mean_score$thetaW, mean_score$thetaS))
cat(sprintf("low violation / bad zeta: violation=%.4f thetaW=%.1f thetaS=%.1f\n",
            low_case$violation$rate, low_case$block$thetaW, low_case$block$thetaS))
cat(sprintf("high violation / good zeta: violation=%.4f thetaW=%.1f thetaS=%.1f coverage=%.1f\n",
            high_case$violation$rate, high_case$block$thetaW, high_case$block$thetaS,
            high_case$block$coverage))
cat("PASS\n")