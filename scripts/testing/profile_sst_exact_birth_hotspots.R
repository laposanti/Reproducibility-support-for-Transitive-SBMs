#!/usr/bin/env Rscript

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm")

out_dir <- "output/diagnostics/sst_profile"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(BayesLogit)
})

source("helper_folder/helper.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("core/my_best_try_so_far.R")
Rcpp::sourceCpp("core/counts_by_block_exact_cpp.cpp")
Rcpp::sourceCpp("core/block_totals_for_poisson_cpp.cpp")

read_mountain_goats <- function() {
  matrix_files <- list.files("data/ShizukaMcDonald_Data", full.names = TRUE, pattern = "[.]csv$")
  n_each <- vapply(matrix_files, function(f) nrow(read.csv(f, row.names = 1)), integer(1))
  A <- as.matrix(read.csv(matrix_files[which.max(n_each)], row.names = 1, check.names = FALSE))
  diag(A) <- 0
  storage.mode(A) <- "numeric"
  A
}

build_edge_data <- function(A) {
  n <- nrow(A)
  N_mat <- A + t(A)
  diag(N_mat) <- 0
  idx <- which(N_mat > 0 & upper.tri(N_mat), arr.ind = TRUE)
  i_idx <- idx[, 1]
  j_idx <- idx[, 2]
  N_edge <- as.numeric(N_mat[idx])
  edge_by_node <- replicate(n, integer(0), simplify = FALSE)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], as.integer(e))
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], as.integer(e))
  }
  list(i_idx = i_idx, j_idx = j_idx, N_edge = N_edge, edge_by_node = edge_by_node)
}

make_profile_state <- function(A, K_state = 24L, seed = 321L) {
  set.seed(seed)
  n <- nrow(A)
  K_state <- min(as.integer(K_state), n - 1L)
  z <- sample(rep(seq_len(K_state), length.out = n))
  eta <- rep(1, n)
  kappa <- matrix(1, K_state, K_state)
  psi <- cumsum(rep(0.25, K_state - 1L))
  edge_data <- build_edge_data(A)
  omega_edge <- draw_edge_omega(
    z, psi, edge_data$i_idx, edge_data$j_idx, edge_data$N_edge, mode = "distance"
  )
  c(list(A = A, z = as.integer(z), eta = eta, kappa = kappa, psi = psi,
         omega_edge = omega_edge), edge_data)
}

node_inputs <- function(state, i) {
  A <- state$A
  z <- state$z
  K_full <- nrow(state$kappa)
  n_minus_full <- tabulate(z[-i], nbins = K_full)
  keep_full <- sort(which(n_minus_full > 0L))
  K_minus <- length(keep_full)
  if (K_minus <= 1L) return(NULL)

  map <- integer(K_full)
  map[keep_full] <- seq_len(K_minus)
  z_packed <- map[z]
  psi_minus <- reindex_psi_sst_keep(state$psi, K_full, keep_full)

  C_i_full <- counts_by_block_exact_cpp(
    i = i, A = A, z = as.integer(z),
    i_idx = as.integer(state$i_idx), j_idx = as.integer(state$j_idx),
    N_edge = as.numeric(state$N_edge), edge_by_node = state$edge_by_node,
    K = K_full
  )
  c_plus <- as.numeric(C_i_full$c_plus)[keep_full]
  N_tot <- as.numeric(C_i_full$N_tot)[keep_full]
  Omega_blk <- aggregate_omega_by_block(
    i, z, state$omega_edge, state$i_idx, state$j_idx, state$edge_by_node, K_full
  )[keep_full]
  old_stats <- .sst_oldold_pair_stats_excluding_i(
    i = i, z_packed = z_packed, A = A,
    i_idx = state$i_idx, j_idx = state$j_idx, N_edge = state$N_edge,
    omega_edge = state$omega_edge, K_minus = K_minus
  )

  list(
    i = i, A = A, z_packed = z_packed,
    i_idx = state$i_idx, j_idx = state$j_idx, N_edge = state$N_edge,
    omega_edge = state$omega_edge, c_plus = c_plus, N_tot = N_tot,
    Omega_blk = Omega_blk, psi_vec = psi_minus,
    S_old = old_stats$S_old, O_old = old_stats$O_old,
    r_set = seq_len(K_minus + 1L), tau0 = 0.15
  )
}

# Independent reference copy of the corrected aggregated formula, kept only as
# a profiling sanity check against the production helper below.
dir_pg_SST_new_exact_nonlocal_corrected_proto <- function(c_plus, N_tot, Omega_blk, psi_vec,
                                                          S_old, O_old,
                                                          r_set = NULL, tau0,
                                                          m_inc = 0) {
  K_old <- length(N_tot)
  R <- K_old + 1L
  if (is.null(r_set)) r_set <- seq_len(R)
  lp <- rep(-Inf, length(r_set))
  if (K_old <= 0L) return(numeric(0))
  if (K_old == 1L) { lp[] <- 0; return(lp) }
  if (length(psi_vec) != (K_old - 1L)) {
    stop("prototype: psi_vec must have length K_old-1", call. = FALSE)
  }

  oldold_terms <- .sst_oldold_nonlocal_birth_terms(S_old, O_old, psi_vec)
  baseline_extreme <- rep(0, R)
  psiKm1 <- psi_vec[K_old - 1L]

  for (a in seq_len(K_old - 1L)) {
    b <- K_old
    if (FALSE) invisible(b)
  }
  for (a in seq_len(K_old - 1L)) {
    for (b in seq.int(a + 1L, K_old)) {
      if ((b - a) != (K_old - 1L)) next
      S_ab <- S_old[a, b]
      O_ab <- O_old[a, b]
      if (S_ab == 0 && O_ab == 0) next
      old_kernel <- S_ab * psiKm1 - 0.5 * O_ab * psiKm1^2
      lo <- a + 1L
      hi <- b
      baseline_extreme[lo:hi] <- baseline_extreme[lo:hi] + old_kernel
    }
  }

  for (idx in seq_along(r_set)) {
    r <- r_set[idx]
    val <- oldold_terms$corr_fixed[r] - baseline_extreme[r]
    YK_node <- 0
    OK_node <- 0

    for (ell in seq_len(K_old)) {
      n <- N_tot[ell]
      if (n == 0) next
      pos_ell <- if (ell < r) ell else ell + 1L
      d <- abs(pos_ell - r)
      if (d == 0L) next
      A_fwd <- if (pos_ell > r) c_plus[ell] else (n - c_plus[ell])
      Y_d <- A_fwd - 0.5 * n

      if (d <= (K_old - 1L)) {
        th <- psi_vec[d]
        val <- val + Y_d * th - 0.5 * Omega_blk[ell] * th^2
      } else if (d == K_old) {
        YK_node <- YK_node + Y_d
        OK_node <- OK_node + Omega_blk[ell]
      }
    }

    YK_tot <- YK_node + oldold_terms$YK_old[r]
    OK_tot <- OK_node + oldold_terms$OK_old[r]
    if (YK_tot != 0 || OK_tot != 0) {
      val <- val + log_int_pg_sst_extreme(YK_tot, OK_tot, psiKm1, tau0, m_inc)
    }
    lp[idx] <- val
  }
  lp
}

score_brute <- function(inp) {
  dir_pg_SST_new_exact_nonlocal_bruteforce(
    i = inp$i, z_packed = inp$z_packed, A = inp$A,
    i_idx = inp$i_idx, j_idx = inp$j_idx, N_edge = inp$N_edge,
    omega_edge = inp$omega_edge, psi_vec = inp$psi_vec,
    r_set = inp$r_set, tau0 = inp$tau0
  )
}

score_aggregated_current <- function(inp) {
  dir_pg_SST_new_exact_nonlocal(
    c_plus = inp$c_plus, N_tot = inp$N_tot, Omega_blk = inp$Omega_blk,
    psi_vec = inp$psi_vec, S_old = inp$S_old, O_old = inp$O_old,
    r_set = inp$r_set, tau0 = inp$tau0
  )
}

score_aggregated_corrected_proto <- function(inp) {
  dir_pg_SST_new_exact_nonlocal_corrected_proto(
    c_plus = inp$c_plus, N_tot = inp$N_tot, Omega_blk = inp$Omega_blk,
    psi_vec = inp$psi_vec, S_old = inp$S_old, O_old = inp$O_old,
    r_set = inp$r_set, tau0 = inp$tau0
  )
}

A <- read_mountain_goats()
state <- make_profile_state(A, K_state = 24L)
set.seed(2026)
nodes <- sample(seq_len(nrow(A)), 20L)
inputs <- Filter(Negate(is.null), lapply(nodes, function(i) node_inputs(state, i)))

cat(sprintf("dataset=mountain_goats n=%d edges=%d K_state=%d profiled_nodes=%d\n",
            nrow(A), length(state$N_edge), max(state$z), length(inputs)))

brute_vals <- lapply(inputs, score_brute)
current_vals <- lapply(inputs, score_aggregated_current)
proto_vals <- lapply(inputs, score_aggregated_corrected_proto)

max_abs_current <- max(unlist(Map(function(a, b) abs(a - b), brute_vals, current_vals)))
max_centered_current <- max(unlist(Map(function(a, b) abs((a - a[1]) - (b - b[1])), brute_vals, current_vals)))
max_abs_proto <- max(unlist(Map(function(a, b) abs(a - b), brute_vals, proto_vals)))
max_centered_proto <- max(unlist(Map(function(a, b) abs((a - a[1]) - (b - b[1])), brute_vals, proto_vals)))

reps <- 4L
time_call <- function(label, fun) {
  gc(FALSE)
  elapsed <- system.time({
    for (rr in seq_len(reps)) {
      invisible(lapply(inputs, fun))
    }
  })[["elapsed"]]
  data.frame(label = label, reps = reps, nodes = length(inputs), elapsed_sec = elapsed,
             calls_per_sec = reps * length(inputs) / elapsed)
}

timing <- rbind(
  time_call("bruteforce_current", score_brute),
  time_call("aggregated_production_exact", score_aggregated_current),
  time_call("aggregated_corrected_proto", score_aggregated_corrected_proto)
)
timing$speedup_vs_bruteforce <- timing$calls_per_sec / timing$calls_per_sec[timing$label == "bruteforce_current"]

accuracy <- data.frame(
  comparison = c("production_aggregated_vs_bruteforce", "corrected_proto_vs_bruteforce"),
  max_abs_error = c(max_abs_current, max_abs_proto),
  max_abs_error_centered = c(max_centered_current, max_centered_proto)
)

make_dense_stress_state <- function(n = 36L, K_state = 12L, seed = 99L) {
  set.seed(seed)
  z <- rep(seq_len(K_state), each = n / K_state)
  A <- matrix(0, n, n)
  for (u in seq_len(n - 1L)) {
    for (v in seq.int(u + 1L, n)) {
      N_uv <- 3L + ((u + v) %% 3L)
      fwd <- 1L + ((2L * u + v) %% max(1L, N_uv - 1L))
      A[u, v] <- fwd
      A[v, u] <- N_uv - fwd
    }
  }
  eta <- rep(1, n)
  kappa <- matrix(1, K_state, K_state)
  psi <- cumsum(rep(0.2, K_state - 1L))
  edge_data <- build_edge_data(A)
  omega_edge <- draw_edge_omega(
    z, psi, edge_data$i_idx, edge_data$j_idx, edge_data$N_edge, mode = "distance"
  )
  c(list(A = A, z = as.integer(z), eta = eta, kappa = kappa, psi = psi,
         omega_edge = omega_edge), edge_data)
}

stress_state <- make_dense_stress_state()
stress_inputs <- Filter(Negate(is.null), lapply(seq_len(nrow(stress_state$A)), function(i) {
  node_inputs(stress_state, i)
}))
stress_brute <- lapply(stress_inputs, score_brute)
stress_current <- lapply(stress_inputs, score_aggregated_current)
stress_proto <- lapply(stress_inputs, score_aggregated_corrected_proto)
stress_accuracy <- data.frame(
  comparison = c("stress_production_aggregated_vs_bruteforce", "stress_corrected_proto_vs_bruteforce"),
  max_abs_error = c(
    max(unlist(Map(function(a, b) abs(a - b), stress_brute, stress_current))),
    max(unlist(Map(function(a, b) abs(a - b), stress_brute, stress_proto)))
  ),
  max_abs_error_centered = c(
    max(unlist(Map(function(a, b) abs((a - a[1]) - (b - b[1])), stress_brute, stress_current))),
    max(unlist(Map(function(a, b) abs((a - a[1]) - (b - b[1])), stress_brute, stress_proto)))
  )
)
accuracy <- rbind(accuracy, stress_accuracy)

write.csv(timing, file.path(out_dir, "sst_birth_scorer_timing.csv"), row.names = FALSE)
write.csv(accuracy, file.path(out_dir, "sst_birth_scorer_accuracy.csv"), row.names = FALSE)

rprof_path <- file.path(out_dir, "sst_birth_bruteforce_Rprof.out")
Rprof(rprof_path, interval = 0.005)
for (rr in seq_len(reps)) invisible(lapply(inputs, score_brute))
Rprof(NULL)
prof <- summaryRprof(rprof_path)
by_total <- as.data.frame(prof$by.total)
by_total$function_name <- rownames(by_total)
by_total <- by_total[, c("function_name", setdiff(names(by_total), "function_name"))]
write.csv(by_total, file.path(out_dir, "sst_birth_bruteforce_Rprof_by_total.csv"), row.names = FALSE)

cat("\nTiming:\n")
print(timing, row.names = FALSE)
cat("\nAccuracy:\n")
print(accuracy, row.names = FALSE)
cat("\nTop Rprof by.total:\n")
print(head(by_total, 12), row.names = FALSE)
cat("\nWrote profiling outputs to ", out_dir, "\n", sep = "")
