#!/usr/bin/env Rscript
# =============================================================================
# build_post_processing.R -- canonical post-processing cube builder.
#
# Reads the raw MCMC fits from output/application/raw/<run>/ and writes
# the SINGLE source of truth for every downstream plot/table:
#
#   output/posterior_post_processing/<run>/
#     per_fit/<dataset>__<model>.rds        canonical z_hat, K_hat, A, K_trace,
#                                            psm, fit metadata
#     per_fit_diagnostics.csv               one row per (ds, model) with K_hat,
#                                            violation_rate_zhat (canonical),
#                                            violation_count_zhat, n_forward,
#                                            n_backward, n_within, total mass,
#                                            and (if available) heavy audit
#                                            columns (p_post_wst, BFs, theta*)
#                                            merged from
#                                            output/application/tables/
#                                            hierarchy_diagnostics_overview.csv
#     vi_pairs.csv                          pairwise VI between models within ds
#     K_posterior.csv                       K mode/mean/95% CI per (ds, model)
#
# Anything that needs z_hat or a violation rate downstream MUST read from
# this cube. Plot scripts no longer call get_z_hat_from_draws() or recompute
# backward_mass on their own.
#
# Usage:
#   APP_RUN_DIR=output/application/raw/application_run_20260414_104327 \
#     Rscript scripts/analysis/build_post_processing.R
#   # Or, with the blessed run already symlinked:
#   Rscript scripts/analysis/build_post_processing.R
#
#   # Force full rebuild even if per_fit/*.rds exist and are newer than the fit:
#   FORCE_REBUILD=1 Rscript scripts/analysis/build_post_processing.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

# --- Resolve run dir --------------------------------------------------------
RUN_DIR <- Sys.getenv("APP_RUN_DIR", unset = "")
if (!nzchar(RUN_DIR)) {
  cur <- "output/paper/figures/current"
  if (file.exists(cur))
    RUN_DIR <- file.path("output/application/raw", basename(Sys.readlink(cur)))
}
if (!nzchar(RUN_DIR) || !dir.exists(RUN_DIR))
  stop("Set APP_RUN_DIR=output/application/raw/<run_id>; got '", RUN_DIR, "'.")

RUN_ID  <- basename(RUN_DIR)
CUBE    <- file.path("output", "posterior_post_processing", RUN_ID)
PER_FIT <- file.path(CUBE, "per_fit")
dir.create(PER_FIT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(APP_RUN_DIR = RUN_DIR,
           APP_CUBE_DIR = CUBE)

FORCE <- identical(tolower(Sys.getenv("FORCE_REBUILD", "")), "1") ||
         identical(tolower(Sys.getenv("FORCE_REBUILD", "")), "true")

cat("============================================================\n")
cat(" Post-processing cube\n")
cat(" Run id:  ", RUN_ID, "\n")
cat(" Cube:    ", CUBE, "\n")
cat(" Force:   ", FORCE, "\n")
cat("============================================================\n\n")

# --- Source canonical machinery --------------------------------------------
source("scripts/analysis/osbm_visualization.R", local = TRUE) # for get_z_hat_from_draws, .compact_labels_row, compute_block_scores
source("helper_folder/transitivity_check_helper.R", local = TRUE) # violation_rate_zhat etc.
source("scripts/analysis/post_processing_helpers.R", local = TRUE)

# --- Locate fits + raw adjacency -------------------------------------------
fit_files <- list.files(RUN_DIR, pattern = "_(WST|SST|DCSBM)_fit\\.rds$",
                        full.names = TRUE)
if (length(fit_files) == 0L)
  stop("No *_(WST|SST|DCSBM)_fit.rds files in ", RUN_DIR)

parse_fit_filename <- function(path) {
  bn <- basename(path)
  m  <- regmatches(bn, regexec("^(.+)_(WST|SST|DCSBM)_fit\\.rds$", bn))[[1]]
  if (length(m) != 3L) return(NULL)
  list(dataset = m[[2]], model = m[[3]], fit_path = path)
}
manifest <- Filter(Negate(is.null), lapply(fit_files, parse_fit_filename))
cat("Found", length(manifest), "fits.\n\n")

# --- Per-fit processing -----------------------------------------------------
per_fit_rows <- list()
for (entry in manifest) {
  ds <- entry$dataset; mdl <- entry$model; fp <- entry$fit_path
  key <- paste0(ds, "__", mdl)
  out_rds <- file.path(PER_FIT, paste0(key, ".rds"))
  cat(sprintf("[%s | %s]\n", ds, mdl))

  needs_rebuild <- FORCE || !file.exists(out_rds) ||
                   file.mtime(out_rds) < file.mtime(fp)
  if (needs_rebuild) {
    cat("  computing z_hat (canonical) ...\n")
    fit <- readRDS(fp)
    A   <- load_adjacency_for_dataset(ds)
    if (is.null(A)) {
      cat("  WARNING: could not load adjacency for ", ds, "; skipping.\n")
      next
    }
    zres <- get_z_hat_from_draws(fit, A, model = mdl)
    K_trace <- extract_K_trace(fit)

    saveRDS(list(
      dataset   = ds,
      model     = mdl,
      fit_path  = fp,
      A         = A,
      z_hat     = zres$z_hat,
      z_chain   = zres$z_chain,
      K_hat     = zres$K,
      K_trace   = K_trace,
      n         = nrow(A),
      created   = Sys.time(),
      script    = "scripts/analysis/build_post_processing.R"
    ), out_rds)
  } else {
    cat("  cached (", basename(out_rds), ", up to date)\n", sep = "")
  }

  # Read back (cheap) for diagnostics row
  obj <- readRDS(out_rds)
  vstats <- compute_zhat_violation_stats(obj$A, obj$z_hat)
  Ksum   <- summarise_K_trace(obj$K_trace)
  n_edges <- sum(obj$A > 0)
  per_fit_rows[[length(per_fit_rows) + 1L]] <- tibble(
    dataset = ds, model = mdl, K_hat = obj$K_hat, n = obj$n,
    n_edges_directed = n_edges,
    violation_count_zhat = vstats$backward_mass,
    cross_count_zhat     = vstats$cross_mass,
    violation_rate_zhat  = vstats$rate,
    n_forward  = vstats$n_forward,
    n_backward = vstats$n_backward,
    n_within   = vstats$n_within,
    K_mode     = Ksum$mode,  K_mean = Ksum$mean,
    K_lo_95    = Ksum$lo,    K_hi_95 = Ksum$hi
  )
}

per_fit_df <- dplyr::bind_rows(per_fit_rows)

# --- Merge heavy audit columns if available --------------------------------
overview_csv <- "output/application/tables/hierarchy_diagnostics_overview.csv"
if (file.exists(overview_csv)) {
  audit <- readr::read_csv(overview_csv, show_col_types = FALSE)
  heavy_cols <- intersect(
    c("p_post_wst", "p_prior_wst", "bf_wst_0", "log_bf_wst_0",
      "p_post_sst", "p_prior_sst", "bf_sst_0", "log_bf_sst_0",
      "bf_sst_wst", "log_bf_sst_wst",
      "thetaW_block_emp_mean", "thetaW_block_emp_lo", "thetaW_block_emp_hi",
      "thetaS_block_emp_mean", "thetaS_block_emp_lo", "thetaS_block_emp_hi",
      "thetaW_item_emp_mean", "thetaW_item_emp_lo", "thetaW_item_emp_hi",
      "thetaS_item_emp_mean", "thetaS_item_emp_lo", "thetaS_item_emp_hi",
      "bar_rho_model_mean", "bar_rho_model_lo", "bar_rho_model_hi",
      "n_draws_K_match",
      "transitive_triads_mean", "transitive_triads_lo.2.5.",
      "transitive_triads_hi.97.5.",
      "hierarchy_energy_mean", "violation_rate_mean"),
    names(audit))
  if (length(heavy_cols)) {
    audit2 <- audit |>
      dplyr::rename(model = fit_model) |>
      dplyr::select(dataset, model, dplyr::all_of(heavy_cols))
    per_fit_df <- per_fit_df |>
      dplyr::left_join(audit2, by = c("dataset", "model"))
    cat("\nMerged ", length(heavy_cols),
        " heavy audit columns from hierarchy_diagnostics_overview.csv.\n",
        sep = "")
  }
} else {
  cat("\n(No hierarchy_diagnostics_overview.csv; per_fit_diagnostics will lack heavy columns.)\n")
}

readr::write_csv(per_fit_df, file.path(CUBE, "per_fit_diagnostics.csv"))
cat("Wrote: ", file.path(CUBE, "per_fit_diagnostics.csv"), "\n", sep = "")

# --- Pairwise VI between models for each dataset ---------------------------
vi_rows <- list()
for (ds in unique(per_fit_df$dataset)) {
  fits <- per_fit_df |> dplyr::filter(dataset == ds)
  if (nrow(fits) < 2) next
  for (i in seq_len(nrow(fits) - 1L)) {
    for (j in seq.int(i + 1L, nrow(fits))) {
      key_i <- paste0(ds, "__", fits$model[i])
      key_j <- paste0(ds, "__", fits$model[j])
      zi <- readRDS(file.path(PER_FIT, paste0(key_i, ".rds")))$z_hat
      zj <- readRDS(file.path(PER_FIT, paste0(key_j, ".rds")))$z_hat
      vi_rows[[length(vi_rows) + 1L]] <- tibble(
        dataset = ds,
        model1 = fits$model[i], K1 = fits$K_hat[i],
        model2 = fits$model[j], K2 = fits$K_hat[j],
        VI = vi_between(zi, zj),
        VI_norm = vi_between(zi, zj) / log(length(zi))
      )
    }
  }
}
if (length(vi_rows)) {
  vi_df <- dplyr::bind_rows(vi_rows)
  readr::write_csv(vi_df, file.path(CUBE, "vi_pairs.csv"))
  cat("Wrote: ", file.path(CUBE, "vi_pairs.csv"), "\n", sep = "")
}

# --- K posterior summary ----------------------------------------------------
K_rows <- list()
for (i in seq_len(nrow(per_fit_df))) {
  r <- per_fit_df[i, ]
  K_rows[[i]] <- tibble(
    dataset = r$dataset, model = r$model,
    K_hat = r$K_hat, K_mode = r$K_mode, K_mean = r$K_mean,
    K_lo_95 = r$K_lo_95, K_hi_95 = r$K_hi_95
  )
}
readr::write_csv(dplyr::bind_rows(K_rows),
                 file.path(CUBE, "K_posterior.csv"))
cat("Wrote: ", file.path(CUBE, "K_posterior.csv"), "\n", sep = "")

# --- README -----------------------------------------------------------------
writeLines(c(
  "# Posterior post-processing cube",
  "",
  "Single source of truth for every downstream plot and table.",
  "Generated by `scripts/analysis/build_post_processing.R`.",
  "",
  paste0("Run id: ", RUN_ID),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Files",
  "",
  "- `per_fit/<dataset>__<model>.rds`: list(z_hat, z_chain, K_hat, K_trace,",
  "  A, n, fit_path, dataset, model). z_hat is canonical:",
  "  WST/SST kept verbatim; DCSBM relabelled by block strength per draw.",
  "- `per_fit_diagnostics.csv`: one row per (dataset, model). Columns:",
  "  K_hat, n, violation_count_zhat, cross_count_zhat,",
  "  violation_rate_zhat (canonical), n_{forward,backward,within},",
  "  K_{mode,mean,lo_95,hi_95}, plus heavy audit columns merged from",
  "  output/application/tables/hierarchy_diagnostics_overview.csv when",
  "  that file exists (bf_*, p_post_*, theta*).",
  "- `vi_pairs.csv`: pairwise VI between models within each dataset.",
  "- `K_posterior.csv`: K mode/mean/95% CI summary."
), file.path(CUBE, "README.md"))

cat("\nDone.\n")
