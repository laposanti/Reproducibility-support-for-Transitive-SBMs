#!/usr/bin/env Rscript
# Compute cross-block violation rates and re-emit agony summaries for each
# (dataset, model) pair. Used to verify the WST claim on citations.

suppressPackageStartupMessages({
  source("scripts/analysis/osbm_visualization.R", chdir = FALSE)
})

RUN_DIR <- "output/application/raw/application_run_20260414_104327"
DATASETS <- c("moreno_sheep", "strauss_2019b", "mountain_goats",
              "citations_data", "macaques_data", "high_school")

violation_summary <- function(A, z_hat) {
  z <- .rank_partition_by_strength(A = A, z_vec = z_hat, alpha = 0.5,
                                   order_direction = "strong_to_weak")
  K <- max(z, na.rm = TRUE)
  ij <- which(A > 0, arr.ind = TRUE)
  if (!nrow(ij)) return(c(viol_rate = NA_real_, viol_count = NA_real_,
                          cross = NA_real_))
  zi <- z[ij[, 1]]; zj <- z[ij[, 2]]
  w  <- A[ij]
  cross <- zi != zj
  back  <- cross & (zi > zj)
  total_cross <- sum(w[cross])
  total_back  <- sum(w[back])
  c(viol_rate = if (total_cross > 0) total_back / total_cross else NA_real_,
    viol_count = total_back,
    cross_count = total_cross,
    K = K)
}

cat(sprintf("%-18s %-7s %5s %8s %8s %8s\n",
            "dataset", "model", "K", "viol", "back", "cross"))
out <- list()
for (ds in DATASETS) {
  ppc_path <- file.path("output/application/ppc", ds, paste0(ds, "_A_obs.rds"))
  A <- if (file.exists(ppc_path)) readRDS(ppc_path) else choose_dataset_local(ds)
  for (m in c("WST", "SST", "DCSBM")) {
    fp <- file.path(RUN_DIR, paste0(ds, "_", m, "_fit.rds"))
    if (!file.exists(fp)) next
    fit <- readRDS(fp)
    z   <- get_z_hat_from_draws(fit, A, model = m)$z_hat
    s   <- violation_summary(A, z)
    out[[length(out) + 1L]] <- data.frame(
      dataset = ds, model = m, K = unname(s["K"]),
      viol_rate = unname(s["viol_rate"]),
      viol_count = unname(s["viol_count"]),
      cross_count = unname(s["cross_count"])
    )
    cat(sprintf("%-18s %-7s %5d %8.3f %8.0f %8.0f\n",
                ds, m, s["K"], s["viol_rate"], s["viol_count"],
                s["cross_count"]))
  }
}
df <- do.call(rbind, out)
out_path <- "output/paper/tables/application_run_20260414_104327/violation_rates_by_model.csv"
readr::write_csv(df, out_path)
cat("\nWrote:", out_path, "\n")
