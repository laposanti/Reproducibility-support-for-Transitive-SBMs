############################################################
## Minimal “cycle breaks DC-SBM” experiment (WST-gen)
## - Generate A from a WST model with a *cyclic* psi matrix
## - Fit DC-SBM, estimate Lambda (upper-tri, posterior mean)
## - Relabel DC-SBM z-chain (i) win-score ordering, (ii) ECR
## - Fit WST-OSBM, estimate Psi (posterior mean)
############################################################

## ---- 0) Packages ----
suppressPackageStartupMessages({
  library(Matrix)
  library(salso)
  library(mcclust)
  library(fossil)
  library(label.switching)  # implements ECR / ECR-iterative etc. :contentReference[oaicite:0]{index=0}
})

## ---- 1) Source *your* code (adjust paths) ----
## You said you already have these functions in your project:
## - simulate_osbm(...)
## - modular_osbm_sampler(...)
## - fit_dcsbm_gibbs_gnedin(...)
## If names differ slightly, just map them below.
setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")

source("./helper_folder/sim_study_helper.R")
source("./helper_folder/WST_helpers.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/Hyper_setup.R")
source("./core/DCSBM_varK.R")

source("./core/my_best_try_so_far.R")

logistic <- function(x) 1 / (1 + exp(-x))

## ---- 2) A tiny helper toolkit ----

# 2.1 Posterior mean of a list of KxK matrices
mat_list_mean <- function(M_list) {
  stopifnot(length(M_list) > 0)
  K <- nrow(M_list[[1]])
  acc <- matrix(0, K, K)
  for (s in seq_along(M_list)) acc <- acc + M_list[[s]]
  acc / length(M_list)
}

# 2.2 Upper-tri extraction (including diagonal) as a matrix
upper_tri_only <- function(M, diag = TRUE) {
  out <- matrix(0, nrow(M), ncol(M))
  out[upper.tri(out, diag = diag)] <- M[upper.tri(M, diag = diag)]
  out
}

# 2.3 Compute “average win score” per block for one allocation z
#     score_k = total wins by block k / total games involving block k
block_win_score <- function(A, z) {
  n <- nrow(A)
  K <- max(z)
  N <- A + t(A)  # total games (i,j)
  scores <- rep(NA_real_, K)
  
  for (k in 1:K) {
    Ik <- which(z == k)
    if (length(Ik) == 0) { scores[k] <- -Inf; next }
    
    wins  <- sum(A[Ik, , drop = FALSE])
    games <- sum(N[Ik, , drop = FALSE])
    scores[k] <- if (games > 0) wins / games else -Inf
  }
  scores
}

# 2.4 Permute labels in z (1..K) by a permutation vector "perm"
#     where perm[new_label] = old_label
apply_label_perm_z <- function(z, perm) {
  # invert mapping: old -> new
  inv <- integer(length(perm))
  inv[perm] <- seq_along(perm)
  inv[z]
}

# 2.5 Apply same permutation to a KxK matrix M (rows/cols are labels)
#     perm[new]=old
apply_label_perm_mat <- function(M, perm) {
  M[perm, perm, drop = FALSE]
}

# 2.6 Relabel one DC-SBM draw by win-score ordering
#     (ties broken by label index for stability)
relabel_one_by_winscore <- function(z, A, Lambda) {
  K <- max(z)
  sc <- block_win_score(A, z)
  # order blocks from “strongest” (highest win share) to weakest
  ord <- order(sc, decreasing = TRUE, na.last = TRUE)
  # ord gives new labels in terms of old labels? We want perm[new]=old:
  perm <- ord
  
  z_new <- apply_label_perm_z(z, perm)
  L_new <- apply_label_perm_mat(Lambda, perm)
  list(z = z_new, Lambda = L_new, perm = perm, score = sc[perm])
}

# 2.7 Relabel an entire chain by win-score (DC-SBM) and keep Lambda aligned
relabel_chain_by_winscore <- function(z_chain, A, Lambda_list) {
  S <- nrow(z_chain)
  K <- max(z_chain)
  
  z_out <- matrix(NA_integer_, S, ncol(z_chain))
  L_out <- vector("list", S)
  perms <- matrix(NA_integer_, S, K)
  
  for (s in 1:S) {
    tmp <- relabel_one_by_winscore(z_chain[s, ], A, Lambda_list[[s]])
    z_out[s, ] <- tmp$z
    L_out[[s]] <- tmp$Lambda
    perms[s, ] <- tmp$perm
  }
  list(z = z_out, Lambda = L_out, perms = perms)
}

# 2.8 Relabel by ECR-iterative-1 (label.switching package)
#     returns permutations m x K; apply them to z and Lambda
relabel_chain_by_ECR <- function(z_chain, Lambda_list) {
  S <- nrow(z_chain)
  K <- max(z_chain)
  
  # label.switching::ecr.iterative.1 expects an m x n array (m=S)
  perms <- label.switching::ecr.iterative.1(z = z_chain, K = K)$permutations
  stopifnot(all(dim(perms) == c(S, K)))
  
  z_out <- matrix(NA_integer_, S, ncol(z_chain))
  L_out <- vector("list", S)
  
  for (s in 1:S) {
    perm <- perms[s, ]  # perm[new]=old
    z_out[s, ] <- apply_label_perm_z(z_chain[s, ], perm)
    L_out[[s]] <- apply_label_perm_mat(Lambda_list[[s]], perm)
  }
  list(z = z_out, Lambda = L_out, perms = perms)
}

# 2.9 Point estimate for z from an S x n relabelled chain
z_point_estimate <- function(z_chain_relab) {
  salso::salso(z_chain_relab, loss = salso::VI(), nRuns = 1L)
}

# 2.10 “DC-SBM implied psi” from Lambda: psi = logit( lambda / (lambda + lambda^T) )
lambda_to_psi <- function(L) {
  vol <- L + t(L)
  rho <- matrix(0.5, nrow(L), ncol(L))
  mask <- vol > 0
  rho[mask] <- L[mask] / vol[mask]
  rho <- pmin(pmax(rho, 1e-6), 1 - 1e-6)
  qlogis(rho)
}

## ---- 3) Data generation: WST with a *cycle* in psi ----
set.seed(1)

K_true <- 3
n      <- 60

# balanced true clustering
z_true <- rep(seq_len(K_true), length.out = n)

# mild degree-correction
eta_true <- runif(n, 0.9, 1.1)

# --------- volumes: highly unbalanced ----------
# Make block 2 vs 3 very "high volume", 1 vs 3 low, 1 vs 2 medium.
kappa_true <- matrix(1, K_true, K_true)
diag(kappa_true) <- 1

kappa_true[1,2] <- kappa_true[2,1] <- 1.0
kappa_true[2,3] <- kappa_true[3,2] <- 6.0
kappa_true[1,3] <- kappa_true[3,1] <- 0.2

# --------- directionality: transitive (WST-consistent) ----------
p12 <- 0.60
p23 <- 0.65
p13 <- 0.55

psi_true <- matrix(0, K_true, K_true)
psi_true[1,2] <- qlogis(p12); psi_true[2,1] <- -qlogis(p12)
psi_true[2,3] <- qlogis(p23); psi_true[3,2] <- -qlogis(p23)
psi_true[1,3] <- qlogis(p13); psi_true[3,1] <- -qlogis(p13)

cat("True psi (transitive):\n"); print(round(psi_true, 2))

# Generate adjacency / win-counts
A <- simulate_osbm(n, K_true, z_true, eta_true, kappa_true, psi_true, regime = "WST")
A <- as.matrix(A)

cat("\nSanity checks:\n")
cat("Total wins:", sum(A), "\n")
cat("Diagonal mass:", sum(diag(A)), "\n")

## ---- 4) Fit DC-SBM (your sampler) ----
## Keep priors simple, and fix K=K_true for the didactic demo.
## If your function name/signature differs, adjust here.

n_iter <- 4000
burn   <- 1000
thin   <- 2
seed_dc <- 10

fit_dc <- fit_dcsbm_gibbs_gnedin(
  toy$A,
  K_init = nrow(A),
  priors = list(
    a_eta       = 1, b_eta = 1,
    a_lambda    = 1, b_lambda = 1,
    gamma_gnedin = 0.9
  ),
  iters   = n_iter,
  burn_in = burn,
  thin    = thin,
  verbose = 200,
  seed    = seed_dc
)

# Extract z-chain and Lambda-chain
# I’m assuming:
#   fit_dc$z is a matrix length Sxn
#   fit_dc$lambda is a list length S, each KxK matrix
z_dc <- fit_dc$z
Lambda_dc <- fit_dc$lambda
S <- nrow(z_dc)


unique_K = function(x) length(unique(x))

K_chain = apply(z_dc, unique_K, MARGIN = 1)
keep = which(K_chain ==3)
cat("\nDC-SBM saved draws:", S, "\n")

## ---- 5) Relabel DC-SBM: (i) win-score ordering ----
rel_win <- relabel_chain_by_winscore(z_dc[keep,], A, Lambda_dc[keep])
z_dc_win <- rel_win$z
Lambda_dc_win <- rel_win$Lambda

z_hat_dc_win <- z_point_estimate(z_dc_win)
ari_dc_win <- fossil::adj.rand.index(z_hat_dc_win, z_true)
vi_dc_win  <- mcclust::vi.dist(z_hat_dc_win, z_true)

cat("\nDC-SBM (win-score relabel) ARI:", round(ari_dc_win, 3),
    " VI:", round(vi_dc_win, 3), "\n")

Lambda_hat_dc_win <- mat_list_mean(Lambda_dc_win)
Lambda_hat_dc_win_ut <- upper_tri_only(Lambda_hat_dc_win, diag = TRUE)

cat("\nPosterior mean Lambda (DC-SBM, win-score relabel), upper-tri:\n")
print(round(Lambda_hat_dc_win_ut, 3))

# Optional: DC-implied psi estimate (mostly to show it struggles with cycles)
psi_hat_dc_win <- lambda_to_psi(Lambda_hat_dc_win)
cat("\nDC-implied psi from Lambda-hat (win-score relabel):\n")
print(round(psi_hat_dc_win, 2))

## ---- 6) Relabel DC-SBM: (ii) ECR (label.switching) ----
## This is the “more renowned algorithm” you mentioned – ECR and variants live here. :contentReference[oaicite:1]{index=1}
rel_ecr <- relabel_chain_by_ECR(z_chain = z_dc[keep,], Lambda_list = Lambda_dc[keep])
z_dc_ecr <- rel_ecr$z
Lambda_dc_ecr <- rel_ecr$Lambda

z_hat_dc_ecr <- z_point_estimate(z_dc_ecr)
ari_dc_ecr <- fossil::adj.rand.index(z_hat_dc_ecr, z_true)
vi_dc_ecr  <- mcclust::vi.dist(z_hat_dc_ecr, z_true)

cat("\nDC-SBM (ECR relabel) ARI:", round(ari_dc_ecr, 3),
    " VI:", round(vi_dc_ecr, 3), "\n")

Lambda_hat_dc_ecr <- mat_list_mean(Lambda_dc_ecr)
Lambda_hat_dc_ecr_ut <- upper_tri_only(Lambda_hat_dc_ecr, diag = TRUE)

cat("\nPosterior mean Lambda (DC-SBM, ECR relabel), upper-tri:\n")
print(round(Lambda_hat_dc_ecr_ut, 3))

psi_hat_dc_ecr <- lambda_to_psi(Lambda_hat_dc_ecr)
cat("\nDC-implied psi from Lambda-hat (ECR relabel):\n")
print(round(psi_hat_dc_ecr, 2))

## ---- 7) Fit WST model (your OSBM sampler) ----
## Here I force fixed K=K_true to keep it comparable and keep the demo short.
seed_wst <- 20

out_wst <- modular_osbm_sampler(
  A = toy$A,
  K = nrow(A),
  n_iter = n_iter,
  burn   = 20,
  thin   = thin,
  verbose = TRUE,
  a_kappa = 1, b_kappa = 1,
  a_eta   = 0.5, b_eta   = 1,
  psi_constraint = "WST" ,
  partition_prior = "GN",
  gamma_gn        = 0.8,
  mu0   = 0.1,
  sigma0 = 0.5,
  tau0   = 1,
  seed = 1
)

cat("\nWST psi posterior mean (reconstructed skew):\n")
print(round(psi_hat_wst_full, 3))

# out_wst$psi : list length S_wst, each KxK psi matrix (as in your pipeline)
# out_wst$z   : list length S_wst, each allocation vector length n
z_wst <- do.call(rbind, lapply(out_wst$z, as.integer))
psi_wst_list <- out_wst$psi
S_wst <- nrow(z_wst)
cat("\nWST saved draws:", S_wst, "\n")

K_unique = apply(z_wst, unique_K, MARGIN = 1)
keep_wst = which(K_unique == 3)

z_hat_wst <- z_point_estimate(z_wst)
ari_wst <- fossil::adj.rand.index(z_hat_wst, toy$z_true)
vi_wst  <- mcclust::vi.dist(z_hat_wst, toy$z_true)

cat("\nWST (win-score relabel) ARI:", round(ari_wst, 3),
    " VI:", round(vi_wst, 3), "\n")

psi_hat_wst <- mat_list_mean(psi_wst_list[keep_wst])
make_skew_from_upper <- function(PSI_upper) {
  PSI <- PSI_upper
  PSI[lower.tri(PSI)] <- -t(PSI)[lower.tri(PSI)]
  diag(PSI) <- 0
  PSI
}

psi_hat_wst_full <- make_skew_from_upper(psi_hat_wst)

cat("\nPosterior mean psi (WST fit, win-score relabel):\n")
print(round(psi_hat_wst_full, 2))

cat("\nTrue psi:\n")
print(round(psi_true, 2))

## ---- 8) A compact “did it fail?” summary ----
cat("\n================= SUMMARY =================\n")
cat("DC-SBM ARI (win-score):", round(ari_dc_win, 3), " VI:", round(vi_dc_win, 3), "\n")
cat("DC-SBM ARI (ECR):      ", round(ari_dc_ecr, 3), " VI:", round(vi_dc_ecr, 3), "\n")
cat("WST    ARI (win-score):", round(ari_wst, 3), " VI:", round(vi_wst, 3), "\n")

cat("\nHint for your write-up:\n")
cat("A cycle at the block level cannot be made consistent with any single latent ordering.\n")
cat("DC-SBM tries to absorb directionality into Lambda, but without a transitivity-aware psi structure,\n")
cat("it tends to smear the signal (or create awkward partitions) when the block tournament is non-transitive.\n")
cat("===========================================\n")



DI            <- build_dyad_index(nrow(A))


LL_dc  <- loglik_matrix_dcsbm(A, fit_dc, dyad_index = DI)
LL_wst <- loglik_matrix_modular(A, out_wst, regime = fit_model, dyad_index = DI)
loo_dc <- tryCatch(loo::loo(LL_dc),
                    error = function(e) list(estimates = matrix(NA, 3, 1),
                                             diagnostics = list(pareto_k = 0)))
loo_wst <- tryCatch(loo::loo(LL_wst),
                   error = function(e) list(estimates = matrix(NA, 3, 1),
                                            diagnostics = list(pareto_k = 0)))

elpd_dc   <- loo_dc$estimates["elpd_loo", "Estimate"]
looic_dc  <- loo_dc$estimates["looic",   "Estimate"]
pk_bad_dc <- mean(loo_dc$diagnostics$pareto_k > 0.7)


elpd_wst   <- loo_wst$estimates["elpd_loo", "Estimate"]
looic_wst  <- loo_wst$estimates["looic",   "Estimate"]
pk_bad_wst <- mean(loo_wst$diagnostics$pareto_k > 0.7)


model_comparison = data.frame(
  Model = c("DC-SBM", "WST-OSBM"),
  ELPD = c(elpd_dc, elpd_wst),
  LOOIC = c(looic_dc, looic_wst),
  Pk_bad = c(pk_bad_dc, pk_bad_wst),
  ari = c(ari_dc_win, ari_wst),
  vi = c(vi_dc_win, vi_wst)
) 


model_comparison%>%
  kable(format='latex', booktabs=T, digits=3, caption='Model comparison: DC-SBM vs WST-OSBM on cyclic data')%>%
  kable_styling(latex_options = c("hold_position"))


