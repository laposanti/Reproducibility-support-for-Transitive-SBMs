#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

get_repo_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grepl("^--file=", cmd)]
  if (!length(file_arg)) return(getwd())
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]),
                               winslash = "/", mustWork = TRUE)
  normalizePath(file.path(dirname(script_path), "../.."),
                winslash = "/", mustWork = TRUE)
}

setwd(get_repo_root())

args <- commandArgs(trailingOnly = TRUE)
sim_file <- if (length(args) >= 1L) args[[1L]] else NA_character_
out_dir <- if (length(args) >= 2L) args[[2L]] else {
  file.path("output", "diagnostics", paste0("ordered_edge_features_", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

choose_dataset <- function(dataset = c("moreno_sheep", "strauss_2019b", "mountain_goats",
                                       "citations_data", "macaques_data", "high_school")) {
  dataset <- match.arg(dataset)
  if (dataset == "moreno_sheep") {
    edges <- read.csv("data/moreno_sheep/edges.csv", header = FALSE,
                      comment.char = "#", strip.white = TRUE)
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    nodes <- sort(unique(c(edges$source, edges$target)))
    A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (r in seq_len(nrow(edges))) A[as.character(edges$source[r]), as.character(edges$target[r])] <- edges$weight[r]
  } else if (dataset == "strauss_2019b") {
    edges <- read.csv("data/Strauss_2019b/edges.csv", header = FALSE,
                      comment.char = "#", strip.white = TRUE)
    names(edges)[1:2] <- c("source", "target")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- 1L
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    nodes <- sort(unique(c(edges$source, edges$target)))
    A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (r in seq_len(nrow(edges))) A[as.character(edges$source[r]), as.character(edges$target[r])] <- A[as.character(edges$source[r]), as.character(edges$target[r])] + 1L
  } else if (dataset == "mountain_goats") {
    matrix_files <- list.files("data/ShizukaMcDonald_Data", full.names = TRUE, pattern = "\\.csv$")
    n_each <- vapply(matrix_files, function(f) nrow(read.csv(f, row.names = 1)), integer(1))
    A <- as.matrix(read.csv(matrix_files[which.max(n_each)], row.names = 1, check.names = FALSE))
  } else if (dataset == "citations_data") {
    A <- as.matrix(read.csv("data/Citations_application/cross-citation-matrix.csv",
                            row.names = 1, header = TRUE, check.names = FALSE))
  } else if (dataset == "macaques_data") {
    edge_list <- read.table("data/macaques/out.moreno.txt")
    nodes <- sort(unique(c(edge_list[[1]], edge_list[[2]])))
    A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (i in seq_len(nrow(edge_list))) A[edge_list[i, 1], edge_list[i, 2]] <- edge_list[i, "V3"]
  } else if (dataset == "high_school") {
    edges <- read.csv("data/high-school/edges.csv", header = FALSE, comment.char = "#", strip.white = TRUE)
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    A <- matrix(0L, max(c(edges$source, edges$target)), max(c(edges$source, edges$target)))
    for (r in seq_len(nrow(edges))) if (edges$weight[r] > 0L) A[edges$source[r], edges$target[r]] <- A[edges$source[r], edges$target[r]] + edges$weight[r]
  }
  A <- as.matrix(A)
  diag(A) <- 0L
  A
}

matrix_features <- function(A, label, source, K_ref = NA_integer_, best_model = NA_character_) {
  n <- nrow(A)
  dyads <- n * (n - 1) / 2
  i <- row(A)[upper.tri(A)]
  j <- col(A)[upper.tri(A)]
  N <- A[cbind(i, j)] + A[cbind(j, i)]
  nz <- N > 0
  abs_net <- abs(A[cbind(i, j)] - A[cbind(j, i)])
  data.frame(
    source = source,
    dataset = label,
    n = n,
    K_ref = K_ref,
    best_model = best_model,
    directed_edges = sum(A),
    unordered_dyads = dyads,
    mean_total_per_dyad = mean(N),
    nonzero_dyad_fraction = mean(nz),
    mean_total_nonzero_dyad = if (any(nz)) mean(N[nz]) else 0,
    reciprocity_share = sum(pmin(A[cbind(i, j)], A[cbind(j, i)])) / max(sum(A), 1),
    mean_abs_net_per_dyad = mean(abs_net),
    mean_abs_net_nonzero_dyad = if (any(nz)) mean(abs_net[nz]) else 0,
    degree_cv = sd(rowSums(A) + colSums(A)) / mean(rowSums(A) + colSums(A)),
    stringsAsFactors = FALSE
  )
}

app_selection_path <- file.path("output", "paper", "tables", "application_run_20260414_104327", "model_selection_paper.csv")
app_sel <- if (file.exists(app_selection_path)) {
  read.csv(app_selection_path, stringsAsFactors = FALSE) %>%
    group_by(dataset) %>%
    slice_min(rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(dataset, K_ref = K_hat, best_model = model)
} else {
  data.frame(dataset = character(), K_ref = integer(), best_model = character())
}

app_diag_path <- file.path("output", "posterior_post_processing",
                           "application_run_20260414_104327",
                           "per_fit_diagnostics.csv")
raw_app_features <- bind_rows(lapply(c("moreno_sheep", "strauss_2019b", "mountain_goats",
                                        "citations_data", "macaques_data", "high_school"), function(ds) {
  row <- app_sel[app_sel$dataset == ds, , drop = FALSE]
  tryCatch(
    matrix_features(
      choose_dataset(ds), ds, "application_raw",
      K_ref = if (nrow(row)) row$K_ref else NA_integer_,
      best_model = if (nrow(row)) row$best_model else NA_character_
    ),
    error = function(e) NULL
  )
}))
if (file.exists(app_diag_path)) {
  app_diag <- read.csv(app_diag_path, stringsAsFactors = FALSE)
  app_features <- app_diag %>%
    group_by(dataset) %>%
    slice(1) %>%
    ungroup() %>%
    left_join(app_sel, by = "dataset") %>%
    transmute(
      source = "application",
      dataset,
      n = n,
      K_ref = K_ref,
      best_model = best_model,
      directed_edges = n_edges_directed,
      unordered_dyads = n * (n - 1) / 2,
      mean_total_per_dyad = n_edges_directed / unordered_dyads,
      nonzero_dyad_fraction = NA_real_,
      mean_total_nonzero_dyad = NA_real_,
      reciprocity_share = NA_real_,
      mean_abs_net_per_dyad = NA_real_,
      mean_abs_net_nonzero_dyad = NA_real_,
      degree_cv = NA_real_
    )
} else {
  app_features <- raw_app_features
  app_features$source <- "application"
}

all_features <- app_features

if (!is.na(sim_file) && file.exists(sim_file)) {
  sim <- read.csv(sim_file, stringsAsFactors = FALSE)
  pointwise_dir <- file.path(dirname(sim_file), "loo_pointwise")
  pointwise_files <- if (dir.exists(pointwise_dir)) list.files(pointwise_dir, pattern = "\\.csv$", full.names = TRUE) else character()
  if (length(pointwise_files)) {
    sim_best <- sim %>%
      group_by(gen_model, rep_id) %>%
      summarise(
        best_model = first(elpd_best_model),
        K_ref = first(K_true),
        .groups = "drop"
      )
    sim_features <- bind_rows(lapply(pointwise_files, function(path) {
      pw <- read.csv(path, stringsAsFactors = FALSE)
      n <- max(c(pw$i, pw$j), na.rm = TRUE)
      deg <- numeric(n)
      for (r in seq_len(nrow(pw))) {
        deg[pw$i[r]] <- deg[pw$i[r]] + pw$Aij[r] + pw$Aji[r]
        deg[pw$j[r]] <- deg[pw$j[r]] + pw$Aij[r] + pw$Aji[r]
      }
      row <- sim_best[sim_best$gen_model == pw$gen_model[1] & sim_best$rep_id == pw$rep_id[1], , drop = FALSE]
      nz <- pw$Nij > 0
      data.frame(
        source = "simulation",
        dataset = paste0(pw$scenario_id[1], "_rep", pw$rep_id[1]),
        n = n,
        K_ref = if (nrow(row)) row$K_ref else pw$K_true[1],
        best_model = if (nrow(row)) row$best_model else NA_character_,
        directed_edges = sum(pw$Nij),
        unordered_dyads = nrow(pw),
        mean_total_per_dyad = mean(pw$Nij),
        nonzero_dyad_fraction = mean(nz),
        mean_total_nonzero_dyad = if (any(nz)) mean(pw$Nij[nz]) else 0,
        reciprocity_share = sum(pmin(pw$Aij, pw$Aji)) / max(sum(pw$Nij), 1),
        mean_abs_net_per_dyad = mean(abs(pw$Aij - pw$Aji)),
        mean_abs_net_nonzero_dyad = if (any(nz)) mean(abs(pw$Aij[nz] - pw$Aji[nz])) else 0,
        degree_cv = sd(deg) / mean(deg),
        stringsAsFactors = FALSE
      )
    }))
    all_features <- bind_rows(all_features, sim_features)
  }
}

write.csv(all_features, file.path(out_dir, "application_sim_feature_comparison.csv"), row.names = FALSE)
write.csv(app_features, file.path(out_dir, "application_effective_features.csv"), row.names = FALSE)
write.csv(raw_app_features, file.path(out_dir, "application_raw_features.csv"), row.names = FALSE)

app_by_best <- app_features %>%
  mutate(best_family = case_when(
    grepl("DC", best_model) ~ "DC-SBM",
    grepl("WST", best_model) ~ "WST",
    grepl("SST", best_model) ~ "SST",
    TRUE ~ best_model
  )) %>%
  group_by(best_family) %>%
  summarise(
    datasets = paste(dataset, collapse = ";"),
    n_mean = mean(n),
    K_mean = mean(K_ref, na.rm = TRUE),
    mean_total_per_dyad = mean(mean_total_per_dyad),
    nonzero_dyad_fraction = mean(nonzero_dyad_fraction),
    degree_cv = mean(degree_cv),
    .groups = "drop"
  )
write.csv(app_by_best, file.path(out_dir, "application_features_by_best_family.csv"), row.names = FALSE)

raw_app_by_best <- raw_app_features %>%
  mutate(best_family = case_when(
    grepl("DC", best_model) ~ "DC-SBM",
    grepl("WST", best_model) ~ "WST",
    grepl("SST", best_model) ~ "SST",
    TRUE ~ best_model
  )) %>%
  group_by(best_family) %>%
  summarise(
    datasets = paste(dataset, collapse = ";"),
    n_mean = mean(n),
    K_mean = mean(K_ref, na.rm = TRUE),
    raw_mean_total_per_dyad = mean(mean_total_per_dyad),
    raw_nonzero_dyad_fraction = mean(nonzero_dyad_fraction),
    raw_reciprocity_share = mean(reciprocity_share),
    raw_degree_cv = mean(degree_cv),
    .groups = "drop"
  )
write.csv(raw_app_by_best, file.path(out_dir, "application_raw_features_by_best_family.csv"), row.names = FALSE)

message("Saved feature comparison to: ", out_dir)
