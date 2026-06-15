#!/usr/bin/env Rscript
# Quick narrative stats for the citations WST fit.
suppressPackageStartupMessages({
  source("scripts/analysis/osbm_visualization.R", chdir = FALSE)
})
RUN <- "output/application/raw/application_run_20260414_104327"
A <- readRDS("output/application/ppc/citations_data/citations_data_A_obs.rds")
fit <- readRDS(file.path(RUN, "citations_data_WST_fit.rds"))
z   <- get_z_hat_from_draws(fit, A, model = "WST")$z_hat
K   <- max(z)
cat("WST minVI K=", K, "\n")

# posterior-mean rho_{kl} reordered by descending strength
psi_chain <- fit$psi_chain
if (is.null(psi_chain)) {
  cat("no psi_chain (using empirical forward share only)\n")
}
# psi is a list of K_t x K_t matrices over draws
S <- length(psi_chain)
# iterate, evaluate at z_hat partition is not direct since K varies. Compute
# block-pair forward share C_kl/(C_kl + C_lk) weighted on the minVI partition
ord <- .rank_partition_by_strength(A=A, z_vec=z, alpha=0.5,
                                   order_direction="strong_to_weak")
K <- max(ord)
C <- matrix(0, K, K)
for (i in seq_len(nrow(A))) for (j in seq_len(ncol(A))) {
  if (i==j) next
  C[ord[i], ord[j]] <- C[ord[i], ord[j]] + A[i, j]
}
rho <- matrix(NA_real_, K, K)
for (k in seq_len(K-1)) for (l in (k+1):K) {
  s <- C[k,l] + C[l,k]
  if (s > 0) rho[k,l] <- C[k,l]/s
}
rkl <- rho[upper.tri(rho)]
rkl <- rkl[!is.na(rkl)]
cat("Empirical forward share over", length(rkl), "ordered pairs:\n")
cat("  min=", round(min(rkl),3), " max=", round(max(rkl),3),
    " median=", round(median(rkl),3), " mean=", round(mean(rkl),3), "\n")
cat("  share with rho>=0.9:", round(mean(rkl>=0.9),2), "\n")
cat("  share with rho<=0.6:", round(mean(rkl<=0.6),2), "\n")
cat("  number of pairs with rho<=0.55:", sum(rkl<=0.55), "\n")
cat("Total directed cross-block edges:", sum(C[upper.tri(C)]+t(C)[upper.tri(C)]), "\n")
cat("Backward edges:", sum(C[lower.tri(C)]), " forward:", sum(C[upper.tri(C)]),"\n")
