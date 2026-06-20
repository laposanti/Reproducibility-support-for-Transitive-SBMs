#!/usr/bin/env Rscript
# =============================================================================
# post_processing_helpers.R -- small, canonical helpers used by
# build_post_processing.R AND by the consumer plot/table scripts. Sourcing
# this file is cheap: it loads no datasets and runs no models.
#
# Public API:
#   load_adjacency_for_dataset(ds)           -> A (numeric matrix)
#   extract_K_trace(fit)                     -> integer vector
#   summarise_K_trace(K_trace)               -> list(mode, mean, lo, hi)
#   compute_zhat_violation_stats(A, z_hat)   -> list(rate, backward_mass,
#                                                cross_mass, n_forward,
#                                                n_backward, n_within)
#   vi_between(z1, z2)                       -> numeric (variation of info)
#
#   read_cube_meta(run_dir = NULL)           -> list(cube_dir, per_fit_dir,
#                                                diag_df, runid)
#   load_fit_cube(ds, model, cube_meta)      -> the cached per-fit list
# =============================================================================

# load_adjacency_for_dataset: delegate to choose_dataset_local() defined in
# scripts/analysis/osbm_visualization.R. Returns NULL if it cannot find a
# loader for the dataset name.
load_adjacency_for_dataset <- function(ds) {
  if (!exists("choose_dataset_local", mode = "function")) {
    stop("choose_dataset_local() not found; source osbm_visualization.R first.")
  }
  tryCatch(choose_dataset_local(ds),
           error = function(e) {
             message("load_adjacency_for_dataset: ", conditionMessage(e))
             NULL
           })
}

# extract_K_trace mirrors the logic in scripts/analysis/posterior_K_uq.R:
# prefer fit$K_trace (or fit$out$K_trace) if present and valid; otherwise
# compute from z draws as the number of distinct labels per row.
extract_K_trace <- function(fit) {
  if (!is.null(fit$K_trace) && length(fit$K_trace) > 0 &&
      all(fit$K_trace > 0, na.rm = TRUE))
    return(as.integer(fit$K_trace[fit$K_trace > 0]))
  if (!is.null(fit$out$K_trace) && length(fit$out$K_trace) > 0 &&
      all(fit$out$K_trace > 0, na.rm = TRUE))
    return(as.integer(fit$out$K_trace[fit$out$K_trace > 0]))
  if (!exists("z_to_matrix", mode = "function"))
    stop("z_to_matrix() not found; source osbm_visualization.R first.")
  z <- z_to_matrix(fit$z)
  apply(z, 1, function(r) length(unique(r[!is.na(r)])))
}

summarise_K_trace <- function(K_trace) {
  K_trace <- K_trace[is.finite(K_trace)]
  if (!length(K_trace))
    return(list(mode = NA_integer_, mean = NA_real_,
                lo = NA_integer_, hi = NA_integer_))
  tbl <- table(K_trace)
  list(
    mode = as.integer(names(which.max(tbl))),
    mean = mean(K_trace),
    lo   = as.integer(quantile(K_trace, 0.025, type = 1)),
    hi   = as.integer(quantile(K_trace, 0.975, type = 1))
  )
}

# Canonical violation stats for a single z_hat against adjacency A.
# Mirrors helper_folder/diagnostics/transitivity_diagnostics.R: an edge is a violation
# iff z[i] > z[j] for a directed edge i -> j with z[i] != z[j].
compute_zhat_violation_stats <- function(A, z_hat) {
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  rows <- row(A); cols <- col(A); vals <- A
  z_i  <- z_hat[rows]; z_j <- z_hat[cols]
  pos  <- vals > 0
  fwd  <- pos & (z_i < z_j)
  bwd  <- pos & (z_i > z_j)
  win  <- pos & (z_i == z_j)
  cross_mass    <- sum(vals[fwd | bwd])
  backward_mass <- sum(vals[bwd])
  list(
    rate          = if (cross_mass > 0) backward_mass / cross_mass else 0,
    backward_mass = backward_mass,
    cross_mass    = cross_mass,
    n_forward     = sum(fwd),
    n_backward    = sum(bwd),
    n_within      = sum(win)
  )
}

# Variation of information between two partitions (label-invariant).
# Self-contained implementation to avoid hard mcclust dependency at read time.
vi_between <- function(z1, z2) {
  stopifnot(length(z1) == length(z2))
  n <- length(z1)
  tab <- table(z1, z2)
  p   <- tab / n
  p1  <- rowSums(p); p2 <- colSums(p)
  h1  <- -sum(p1[p1 > 0] * log(p1[p1 > 0]))
  h2  <- -sum(p2[p2 > 0] * log(p2[p2 > 0]))
  mi  <- sum(p[p > 0] * (log(p[p > 0]) -
           log(outer(p1, p2)[p > 0])))
  h1 + h2 - 2 * mi
}

# Resolve the cube directory for the active run and read the diagnostics CSV.
read_cube_meta <- function(run_dir = NULL) {
  if (is.null(run_dir) || !nzchar(run_dir))
    run_dir <- Sys.getenv("APP_RUN_DIR", unset = "")
  if (!nzchar(run_dir)) {
    source("scripts/bundle_defaults.R", local = TRUE)
    run_dir <- bundle_resolve_application_run_dir(must_exist = TRUE)
  }
  if (!nzchar(run_dir) || !dir.exists(run_dir))
    stop("read_cube_meta: cannot resolve run dir.")
  runid    <- basename(run_dir)
  cube_dir <- file.path("output", "posterior_post_processing", runid)
  diag_csv <- file.path(cube_dir, "per_fit_diagnostics.csv")
  if (!file.exists(diag_csv))
    stop("Cube not built yet: ", diag_csv,
         " — run Rscript scripts/analysis/build_post_processing.R first.")
  diag_df <- utils::read.csv(diag_csv, stringsAsFactors = FALSE,
                             check.names = FALSE)
  list(run_dir = run_dir, runid = runid,
       cube_dir = cube_dir,
       per_fit_dir = file.path(cube_dir, "per_fit"),
       diag_df = diag_df)
}

load_fit_cube <- function(ds, model, cube_meta = NULL) {
  if (is.null(cube_meta)) cube_meta <- read_cube_meta()
  fp <- file.path(cube_meta$per_fit_dir,
                  paste0(ds, "__", model, ".rds"))
  if (!file.exists(fp))
    stop("load_fit_cube: missing cube entry ", fp,
         " — rerun build_post_processing.R.")
  readRDS(fp)
}

# Convenience: get just the canonical z_hat / violation rate without loading
# the heavy A matrix.
get_violation_rate_cube <- function(ds, model, cube_meta = NULL) {
  if (is.null(cube_meta)) cube_meta <- read_cube_meta()
  row <- cube_meta$diag_df[cube_meta$diag_df$dataset == ds &
                           cube_meta$diag_df$model == model, ]
  if (nrow(row) == 0L)
    stop(sprintf("get_violation_rate_cube: %s/%s not in cube.", ds, model))
  list(rate          = as.numeric(row$violation_rate_zhat),
       backward_mass = as.numeric(row$violation_count_zhat),
       cross_mass    = as.numeric(row$cross_count_zhat),
       K_hat         = as.integer(row$K_hat))
}
