#!/usr/bin/env Rscript

# Focused SST recovery diagnostic:
#   Does VI decrease as kappa increases under SST-generated data?
#
# The script compares two SST fit configurations:
#   - demo_current: mirrors DemoKvar's current simulation choices
#                   (K_init=n, K_expected=K_true, b_kappa=b_eta=1).
#   - truth_start: a truth-initialized validation path
#                  (K_init=K_true, K_expected=K_true, b_kappa=b_eta=1).
#
# It writes a CSV and a short Markdown report under
# output/diagnostics/sst_kappa_vi/.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(salso)
  library(mcclust)
  library(fossil)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
repo_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), "../.."),
                winslash = "/", mustWork = FALSE)
} else {
  getwd()
}
if (!file.exists(file.path(repo_root, "helper_folder", "SST_helpers.R"))) repo_root <- getwd()
setwd(repo_root)

source("helper_folder/helper.R")
source("helper_folder/sim_study_helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/mixing_moves.R")
source("helper_folder/Hyper_setup.R")
source("core/my_best_try_so_far.R")

parse_num_grid <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  vals <- suppressWarnings(as.numeric(strsplit(raw, ",")[[1]]))
  vals[is.finite(vals)]
}

parse_int <- function(name, default) {
  val <- suppressWarnings(as.integer(Sys.getenv(name, unset = "")))
  if (is.na(val)) default else val
}

out_dir <- "output/diagnostics/sst_kappa_vi"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_items <- parse_int("SST_KAPPA_DIAG_N", 36L)
K_true <- parse_int("SST_KAPPA_DIAG_K", 4L)
n_rep <- parse_int("SST_KAPPA_DIAG_REPS", 1L)
n_iter <- parse_int("SST_KAPPA_DIAG_N_ITER", 800L)
burn <- parse_int("SST_KAPPA_DIAG_BURN", 250L)
thin <- parse_int("SST_KAPPA_DIAG_THIN", 2L)
kappa_grid <- parse_num_grid("SST_KAPPA_DIAG_KAPPA", c(0.75, 1.5, 3, 6))

psi_true <- parse_num_grid("SST_KAPPA_DIAG_PSI", c(0.7, 1.15, 1.55))
if (length(psi_true) != K_true - 1L) {
  psi_true <- cumsum(rep(1.2 / max(K_true - 1L, 1L), K_true - 1L))
}

make_data <- function(kappa_mean, rep_id) {
  set.seed(7000L + rep_id + round(100 * kappa_mean))
  z_true <- rep(seq_len(K_true), length.out = n_items)
  z_true <- z_true[sample.int(n_items)]
  eta_true <- rep(1, n_items)
  kappa_true <- matrix(kappa_mean, K_true, K_true)
  A <- as.matrix(simulate_osbm(n_items, K_true, z_true, eta_true,
                               kappa_true, psi_true, "SST"))
  list(A = A, z_true = z_true)
}

fit_one <- function(A, z_true, config, seed) {
  if (identical(config, "demo_current")) {
    K_init <- nrow(A)
    K_expected <- K_true
    b_kappa <- 1
    b_eta <- 1
  } else if (identical(config, "truth_start")) {
    K_init <- K_true
    K_expected <- K_true
    b_kappa <- 1
    b_eta <- 1
  } else {
    stop("Unknown config: ", config)
  }

  hypers <- get_principled_hypers_v2(
    A = A,
    K_expected = K_expected,
    c_kappa = 3,
    ordering_prior_mode = "equivalence_class"
  )
  hypers$tau0 <- max(hypers$tau0, 0.2)

  fit <- modular_osbm_sampler(
    A = A,
    K = K_init,
    free = c("psi", "kappa", "eta", "z"),
    n_iter = n_iter,
    burn = burn,
    thin = thin,
    verbose = FALSE,
    psi_constraint = "SST",
    seed = seed,
    a_kappa = 1,
    b_kappa = b_kappa,
    a_eta = 1,
    b_eta = b_eta,
    mu0 = hypers$mu0,
    sigma0 = hypers$sigma0,
    tau0 = hypers$tau0,
    gamma_gn = hypers$gamma,
    partition_prior = "OCRP",
    theta_ocrp = 0.5,
    eta_identifiability = "block_sum",
    use_mixing_moves = TRUE,
    sst_birth_score_mode = "exact_nonlocal"
  )

  z_chain <- do.call(rbind, lapply(fit$z, as.integer))
  z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L)
  K_trace <- vapply(fit$z, function(z) length(unique(z)), integer(1))
  draw_vi <- vapply(fit$z, function(z) mcclust::vi.dist(z, z_true), numeric(1))

  data.frame(
    config = config,
    K_init = K_init,
    K_expected = K_expected,
    b_kappa = b_kappa,
    b_eta = b_eta,
    tau0 = hypers$tau0,
    gamma = hypers$gamma,
    K_hat = length(unique(z_hat)),
    K_mean = mean(K_trace),
    K_mode = as.integer(names(which.max(table(K_trace)))),
    ari = fossil::adj.rand.index(z_hat, z_true),
    vi = mcclust::vi.dist(z_hat, z_true),
    draw_vi_mean = mean(draw_vi),
    stringsAsFactors = FALSE
  )
}

rows <- list()
for (rep_id in seq_len(n_rep)) {
  for (kappa_mean in kappa_grid) {
    dat <- make_data(kappa_mean, rep_id)
    for (config in c("demo_current", "truth_start")) {
      cat(sprintf("rep=%d kappa=%.3f config=%s ... ", rep_id, kappa_mean, config))
      t0 <- Sys.time()
      row <- tryCatch(
        fit_one(dat$A, dat$z_true, config, seed = 9000L + rep_id * 100L + round(kappa_mean * 10)),
        error = function(e) data.frame(
          config = config, K_init = NA_integer_, K_expected = NA_integer_,
          b_kappa = NA_real_, b_eta = NA_real_, tau0 = NA_real_, gamma = NA_real_,
          K_hat = NA_integer_, K_mean = NA_real_, K_mode = NA_integer_,
          ari = NA_real_, vi = NA_real_, draw_vi_mean = NA_real_,
          error = conditionMessage(e), stringsAsFactors = FALSE
        )
      )
      row$rep_id <- rep_id
      row$kappa_mean <- kappa_mean
      row$n <- n_items
      row$K_true <- K_true
      rows[[length(rows) + 1L]] <- row
      cat(sprintf("done %.1fs\n", as.numeric(Sys.time() - t0, units = "secs")))
    }
  }
}

res <- bind_rows(rows) |>
  relocate(rep_id, kappa_mean, config, n, K_true)

csv_path <- file.path(out_dir, "sst_kappa_vi_diagnostic.csv")
readr::write_csv(res, csv_path)

summary <- res |>
  group_by(config, kappa_mean) |>
  summarise(
    n_runs = n(),
    vi_mean = mean(vi, na.rm = TRUE),
    ari_mean = mean(ari, na.rm = TRUE),
    exact_K_rate = mean(K_hat == K_true, na.rm = TRUE),
    K_mean = mean(K_mean, na.rm = TRUE),
    tau0 = mean(tau0, na.rm = TRUE),
    .groups = "drop"
  )

summary_path <- file.path(out_dir, "sst_kappa_vi_summary.csv")
readr::write_csv(summary, summary_path)

report_path <- file.path(out_dir, "sst_kappa_vi_report.md")
lines <- c(
  "# SST kappa-VI diagnostic",
  "",
  sprintf("- n=%d, K*=%d, reps=%d, n_iter=%d, burn=%d, thin=%d",
          n_items, K_true, n_rep, n_iter, burn, thin),
  sprintf("- psi_true=(%s)", paste(sprintf("%.3f", psi_true), collapse = ", ")),
  sprintf("- CSV: `%s`", csv_path),
  sprintf("- Summary CSV: `%s`", summary_path),
  "",
  "| config | kappa | VI | ARI | Pr(Khat=K*) | K_mean | tau0 |",
  "|---|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(summary))) {
  r <- summary[i, ]
  lines <- c(lines, sprintf(
    "| %s | %.3f | %.3f | %.3f | %.2f | %.2f | %.4f |",
    r$config, r$kappa_mean, r$vi_mean, r$ari_mean,
    r$exact_K_rate, r$K_mean, r$tau0
  ))
}

lines <- c(lines, "",
           "Interpretation:",
           "- `demo_current` checks the production-style overfitted start after the current calibration changes.",
           "- `truth_start` is the basic posterior-concentration check: VI should decline as kappa grows, modulo short-chain Monte Carlo error.",
           "- If `truth_start` still fails at large kappa, inspect the SST z-update and psi update before running larger studies.")
writeLines(lines, report_path)

cat("Wrote:\n")
cat(" - ", csv_path, "\n", sep = "")
cat(" - ", summary_path, "\n", sep = "")
cat(" - ", report_path, "\n", sep = "")
