# Shared ordered SBM utilities for math helpers, relabelling, and log-likelihood summaries.
#
###########################################################################
#  ---------------------    Helpers     -----------------------------------
#                       Lapo - July 2025
###########################################################################

##---------- helpers ----------
logistic      <- function(x) 1/(1+exp(-x)) #sigmoid/logistic function
#This is a wrapper for the Polya-Gamma sampler from the BayesLogit package.
rpg_vec       <- function(N, psi) BayesLogit::rpg(num = length(N), 
                                                  h = N,  #counts per dyad (N_ij)
                                                  z = psi) #the psi values at the appropriate distances

#To speed up the computation of the normalizing constant 
#of the Bernoulli/Binomial likelihood when computing 
#the directional contribution in the z update step:
#--> ...-N_ij log(1 + exp(s_ij psi_d))
#pos --> log(1+ exp(psi_d))
#neg --> log(1+ exp(psi_d))
make_log1p_tables <- function(psi){ 
  list(pos = log1pexp( psi),
       neg = log1pexp(-psi))}

#numerically stable version of log(1+exp(x)) that avoids overflow.
log1pexp <- function(x) ifelse(x > 50, x, log1p(exp(x))) #log1p(exp(x)) for better precision

## ---------- omega update ----------
draw_omega_bar <- function(B_d, psi) {
  if (length(B_d) != length(psi)) {
    stop(sprintf(paste0(
      "[draw_omega_bar] length(B_d)=%d implies K=%d, but length(psi)=%d. ",
      "Labels must be contiguous 1:K."
    ), length(B_d), length(B_d) + 1L, length(psi)), call. = FALSE)
  }
  bar_omega <- numeric(length(B_d))
  keep      <- B_d > 0   # only distances with edges (rpg fails with B_d=0)
  if (any(keep)) {
    if (anyNA(psi[keep]) || any(!is.finite(psi[keep]))) {
      stop("[draw_omega_bar] non-finite tilt in psi[keep]", call. = FALSE)
    }
    bar_omega[keep] <- BayesLogit::rpg(sum(keep),
                                       h = B_d[keep],
                                       z = abs(psi[keep]))
  }
  bar_omega                                     # zeros where B_d == 0
}

## ---------- psi update (Normal or half-Normal prior) -------------------------
## returns two length-(K-1) vectors
##  bar_y[d]     = sum S_{ij},      
##  bar_omega[d] = sum omega_{ij}   both at distance d.


#---------
#WST model 
#---------

# log(1+exp) lookup for pairwise ψ (upper-tri only is meaningful)
make_log1p_tables_pair <- function(psi_mat){
  stopifnot(is.matrix(psi_mat), nrow(psi_mat) == ncol(psi_mat))
  K <- nrow(psi_mat)
  pos <- matrix(0, K, K)
  neg <- matrix(0, K, K)
  ut  <- which(upper.tri(psi_mat), arr.ind = TRUE)
  if (nrow(ut)) {
    vals      <- psi_mat[ut]
    pos[ut]   <- log1pexp(vals)
    neg[ut]   <- log1pexp(-vals)
  }
  list(pos = pos, neg = neg)
}

aggregate_by_pair <- function(K, z_i, z_j, A_ij, N_edge) {
  # returns bar_y as KxK with upper-tri entries filled, else 0
  bar_y <- matrix(0, K, K)
  k_low  <- pmin(z_i, z_j)
  k_high <- pmax(z_i, z_j)
  forward <- z_i < z_j
  A_fwd   <- ifelse(forward, A_ij, N_edge - A_ij)
  S_vec   <- A_fwd - N_edge/2
  # sum S by (k_low, k_high)
  key <- paste(k_low, k_high, sep = "_")
  sums <- tapply(S_vec, key, sum)
  if (length(sums)) {
    idx <- do.call(rbind, strsplit(names(sums), "_", fixed = TRUE))
    idx <- cbind(as.integer(idx[,1]), as.integer(idx[,2]))
    for (t in seq_along(sums)) {
      k <- idx[t,1]; l <- idx[t,2]
      if (k < l) bar_y[k,l] <- sums[t]
    }
  }
  bar_y
}

pair_totals <- function(K, z_i, z_j, N_edge) {
  B <- matrix(0, K, K)
  k_low  <- pmin(z_i, z_j)
  k_high <- pmax(z_i, z_j)
  key <- paste(k_low, k_high, sep = "_")
  sums <- tapply(N_edge, key, sum)
  if (length(sums)) {
    idx <- do.call(rbind, strsplit(names(sums), "_", fixed = TRUE))
    idx <- cbind(as.integer(idx[,1]), as.integer(idx[,2]))
    for (t in seq_along(sums)) {
      k <- idx[t,1]; l <- idx[t,2]
      if (k < l) B[k,l] <- sums[t]
    }
  }
  B
}


draw_omega_pair <- function(B_mat, psi_mat) {
  bar_omega <- matrix(0, nrow(B_mat), ncol(B_mat))
  ut <- which(upper.tri(B_mat) & B_mat > 0, arr.ind = TRUE)
  if (nrow(ut)) {
    h <- B_mat[ut]
    z <- psi_mat[ut]
    bar_omega[ut] <- BayesLogit::rpg(nrow(ut), h = h, z = z)
  }
  bar_omega
}

update_psi_pair <- function(bar_y, bar_omega, mu0, sig2_0, trunc = TRUE) {
  K <- nrow(bar_y)
  psi_mat <- matrix(0, K, K)
  ut <- which(upper.tri(bar_y), arr.ind = TRUE)
  for (t in seq_len(nrow(ut))) {
    k <- ut[t,1]; l <- ut[t,2]
    prec   <- bar_omega[k,l] + 1/sig2_0
    mean_  <- (bar_y[k,l] + mu0/sig2_0)/prec
    sd_    <- 1/sqrt(prec)
    psi_mat[k,l] <- if (trunc)
      truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = mean_, sd = sd_)
    else
      rnorm(1, mean_, sd_)
  }
  psi_mat
}







## distance-wise aggregates: A signed count bar y_d
# which summarizes directional information about edges 
#between blocks at that distance.
# PATCH: aggregate_by_distance()
aggregate_by_distance <- function(K, z_i, z_j, A_ij, N_edge){ 
  d_vec <- abs(z_i - z_j)
  keep  <- d_vec > 0L
  bar_y <- numeric(K - 1)
  
  if (any(keep)) {
    d_use   <- d_vec[keep]
    # forward = lower block -> higher block
    forward <- z_i[keep] < z_j[keep]
    
    A_fwd   <- ifelse(forward, A_ij[keep], N_edge[keep] - A_ij[keep])
    S_vec   <- A_fwd - N_edge[keep]/2
    
    by_y <- tapply(S_vec, d_use, sum)
    bar_y[as.integer(names(by_y))] <- by_y
  }
  list(bar_y = bar_y)
}



## B_d = sum N_edge  for all dyads currently at distance d  (d = 1,…,K-1)
distance_totals <- function(K, z_i, z_j, N_edge)
{
  #z_i, z_j: block labels of the endpoints of each dyad.
  #d_vec: computes distances between each pair.
  #tapply(N_edge, factor(...)): sums the N_ij counts grouped by block distance d
  #levels = 1:(K-1): ensures missing distances get zero and code doen't break
  
  d_vec <- abs(z_i - z_j) 
  by_B <- tapply(N_edge, factor(d_vec, levels = 1:(K-1)), sum)
  by_B[is.na(by_B)] <- 0
  as.numeric(by_B)
  # Return a vector of length K-1, zeros where no dyad
  # entry d contains the total number of dyads at distance d
}

eta_rescale_by_block <- function(eta, z, K) {
  for (k in seq_len(K)) {
    idx_k <- which(z == k)
    nk    <- length(idx_k)
    if (nk) {
      s <- sum(eta[idx_k])
      eta[idx_k] <- nk * eta[idx_k] / s          #  Σ_i η_i = n_k
    }
  }
  eta
}

update_psi <- function(K, bar_y, bar_omega, mu0, sig2_0, trunc = TRUE)
{
  #This function samples the directional bias psi_d | rest
  # bar_y:     signed count of directional imbalances (from aggregate_by_distance)
  # bar_omega: total Polya-Gamma weights (from draw_omega_bar)
  # mu0, sig2_0 - prior mean and variance
  psi <- numeric(K-1)
  for (d in seq_len(K-1)) {
    prec   <- bar_omega[d] + 1 / sig2_0 
    mean_d <- (bar_y[d] + mu0 / sig2_0) / prec
    sd_d   <- 1 / sqrt(prec)
    
    psi[d] <- if (trunc)
      truncnorm::rtruncnorm(1, a = 0, b = Inf,
                            mean = mean_d, sd = sd_d)
    else
      rnorm(1, mean_d, sd_d)
  }
  #Return A length-K−1 vector with updated psi
  psi
}

# PATCH: log_post_fast_distance()
log_post_fast_distance <- function(i, k,
                                   N_blk, A_out_blk,
                                   n_minus,
                                   eta_sum_blk, eta_i,
                                   kappa, psi, alpha_vec,
                                   log_kappa, log1p_tab)
{
  log2 <- log(2)
  lp   <- log(alpha_vec[k] + n_minus[k])
  lambda_vec_k <- kappa %*% eta_sum_blk
  lp <- lp - eta_i * lambda_vec_k[k] + (log_kappa %*% N_blk)[k]
  
  idx        <- seq_along(alpha_vec)
  mask_cross <- idx != k
  d_vec      <- abs(idx - k)[mask_cross]
  
  sign_vec <- ifelse(idx[mask_cross] > k,  1, -1)
  
  A_dir      <- A_out_blk[mask_cross]
  lp <- lp +
    sum(sign_vec * A_dir * psi[d_vec]) -
    sum(N_blk[mask_cross] *
          ifelse(sign_vec == 1,
                 log1p_tab$pos[d_vec],
                 log1p_tab$neg[d_vec])) -
    N_blk[k] * log2
  lp
}

# Collapsed (per-move) Gamma–Poisson marginal for the κ part, SST (distance ψ).
# This *replaces* the κ-dependent terms in log_post_fast_distance():
#   -eta_i * (kappa %*% eta_sum_blk)[k]  +  (log_kappa %*% N_blk)[k]
#
# It keeps your directional/binomial part intact.
#
# Conventions:
# - forward = z_i < z_j  (stronger -> weaker), consistent with your DEMO.
# - For each candidate k, we consider pairs (min(k,l), max(k,l)) once.
# - ΔR_{kl}: counts contributed by node i to block-pair (k,l)
#     ΔR_{k,l != k} = N_blk[l]
#     ΔR_{k,k}      = N_blk[k]
# - ΔT_{kl}: exposure contributed by node i, using current eta sums excluding i
#     ΔT_{k,l != k} = eta_i * E_l
#     ΔT_{k,k}      = eta_i * E_k
#
# The collapsed log-factor for a single pair is:
#   log Γ(a + ΔR) - log Γ(a) - (a + ΔR) * log(b + ΔT)
# Up to additive constants common to all candidates (e.g., a*log b), which cancel.
compute_loo_waic <- function(log_lik_mat) {
  # r_eff helps PSIS-LOO diagnostics; with 1 chain this is approximate
  r_eff <- tryCatch(
    loo::relative_eff(exp(log_lik_mat), chain_id = rep(1, nrow(log_lik_mat))),
    error = function(e) 1
  )
  loo_obj  <- loo::loo(log_lik_mat, r_eff = r_eff)
  waic_obj <- loo::waic(log_lik_mat)
  
  list(
    loo_obj  = loo_obj,
    waic_obj = waic_obj,
    looic    = loo_obj$estimates["looic", "Estimate"],
    waic     = waic_obj$estimates["waic", "Estimate"]
  )
}


log_post_fast_distance_collapsed <- function(
    i, k,
    N_blk, A_out_blk,
    n_minus,
    eta_sum_blk, eta_i,
    psi, alpha_vec,
    log1p_tab,                     # from make_log1p_tables(psi)
    a_kappa, b_kappa               # Gamma prior hyper-params for κ
){
  K    <- length(alpha_vec)
  idx  <- seq_len(K)
  log2 <- log(2)
  
  # ----- Dirichlet–Multinomial prior on z -----
  lp <- log(alpha_vec[k] + n_minus[k])
  
  # ----- Directional/binomial piece (same as your distance scorer, but consistent sign) -----
  # forward = z_i < l (i.e., lower index -> higher index)
  mask_cross <- idx != k
  d_vec      <- abs(idx - k)[mask_cross]
  sign_vec   <- ifelse(idx[mask_cross] > k,  1, -1)  # +ψ when l > k, -ψ when l < k
  A_dir      <- A_out_blk[mask_cross]
  
  # Bernoulli/Binomial contribution with log(1+exp) caches
  lp <- lp +
    sum(sign_vec * A_dir * psi[d_vec]) -
    sum(N_blk[mask_cross] * ifelse(sign_vec == 1,
                                   log1p_tab$pos[d_vec],
                                   log1p_tab$neg[d_vec])) -
    N_blk[k] * log2
  
  # ----- Collapsed κ piece (Gamma–Poisson marginal per block-pair) -----
  # Build ΔR and ΔT for all l, then sum over unique (min,max) pairs.
  E_k_all <- eta_sum_blk           # sums over blocks excluding i (length K), NA -> 0
  E_k_all[is.na(E_k_all)] <- 0
  
  # ΔR / ΔT vectors indexed by l (block of the other endpoint)
  deltaR <- N_blk
  deltaT <- numeric(K)
  for (l in idx) {
    if (l == k) {
      deltaT[l] <- eta_i * E_k_all[l]
    } else {
      deltaT[l] <- eta_i * E_k_all[l]
    }
  }
  
  # Accumulate over pairs (min(k,l), max(k,l)) once
  # For l != k: one off-diagonal pair; for l == k: within-block pair.
  # Constants a*log(b) cancel across k, so we omit them.
  # We also skip l where ΔR=0 and ΔT=0 (no contribution).
  for (l in idx) {
    dR <- deltaR[l]; dT <- deltaT[l]
    if (dR > 0 || dT > 0) {
      lp <- lp + lgamma(a_kappa + dR) - lgamma(a_kappa) - (a_kappa + dR) * log(b_kappa + dT)
    }
  }
  
  lp
}
# Node-wise Gibbs update for z with collapsed κ (SST distance ψ).
# Drop-in replacement for z_update_osbm_distance().
#
# Usage inside your main loop:
#   z_update_fun <- z_update_osbm_distance_collapsed
#   # and pass a_kappa = hyper$a_kappa, b_kappa = hyper$b_kappa

z_update_osbm_distance_collapsed <- function(
    i, A, z, eta, psi, alpha_vec,
    log1p_tab,
    a_kappa, b_kappa
){
  K    <- length(alpha_vec)
  idxK <- seq_len(K)
  
  # ---- dyad summaries wrt current z (excluding i) ----
  e_out <- A[i, ];  e_in <- A[, i];  Nij <- e_out + e_in
  
  N_blk <- A_out_blk <- numeric(K)     # totals involving i toward each block l
  if (any(Nij > 0)) {
    ell <- z[which(Nij > 0)]
    N_blk_tmp     <- tapply(Nij[Nij > 0], ell, sum)
    A_out_blk_tmp <- tapply(e_out[Nij > 0], ell, sum)
    N_blk[as.integer(names(N_blk_tmp))]         <- N_blk_tmp
    A_out_blk[as.integer(names(A_out_blk_tmp))] <- A_out_blk_tmp
  }
  
  # Block sizes / eta sums excluding i
  n_minus     <- tabulate(z[-i], nbins = K)
  eta_sum_blk <- tapply(eta[-i], factor(z[-i], levels = idxK), sum)
  eta_sum_blk[is.na(eta_sum_blk)] <- 0
  
  # ---- score each candidate k with collapsed κ ----
  lp_vec <- sapply(
    idxK, log_post_fast_distance_collapsed, i = i,
    N_blk       = N_blk,
    A_out_blk   = A_out_blk,
    n_minus     = n_minus,
    eta_sum_blk = eta_sum_blk,
    eta_i       = eta[i],
    psi         = psi,
    alpha_vec   = alpha_vec,
    log1p_tab   = log1p_tab,
    a_kappa     = a_kappa,
    b_kappa     = b_kappa
  )
  
  p <- exp(lp_vec - max(lp_vec))
  sample.int(K, 1L, prob = p)
}




# Node-wise Gibbs update for z (distance ψ)
z_update_osbm_distance <- function(i, A, z, eta, kappa, psi, alpha_vec,
                                   log1p_tab, log_kappa,
                                   i_idx = NULL, j_idx = NULL, N_edge = NULL,
                                   edge_by_node = NULL)
{
  K <- length(alpha_vec)
  idxK <- seq_len(K)
  
  # dyad summaries wrt current z (same as your original code)
  e_out <- A[i, ];  e_in <- A[, i];  Nij <- e_out + e_in
  N_blk <- A_out_blk <- numeric(K)
  if (any(Nij > 0)) {
    ell <- z[which(Nij > 0)]
    N_blk_tmp     <- tapply(Nij[Nij > 0], ell, sum)
    A_out_blk_tmp <- tapply(e_out[Nij > 0], ell, sum)
    N_blk[as.integer(names(N_blk_tmp))]         <- N_blk_tmp
    A_out_blk[as.integer(names(A_out_blk_tmp))] <- A_out_blk_tmp
  }
  
  n_minus     <- tabulate(z[-i], nbins = K)
  eta_sum_blk <- tapply(eta[-i], factor(z[-i], levels = idxK), sum)
  eta_sum_blk[is.na(eta_sum_blk)] <- 0
  
  lp_vec <- sapply(
    idxK, log_post_fast_distance, i           = i,
    N_blk       = N_blk,
    A_out_blk   = A_out_blk,
    n_minus     = n_minus,
    eta_sum_blk = eta_sum_blk,
    eta_i       = eta[i],
    kappa       = kappa,
    psi         = psi,
    alpha_vec   = alpha_vec,
    log_kappa   = log_kappa,
    log1p_tab   = log1p_tab
  )
  
  p <- exp(lp_vec - max(lp_vec))
  sample.int(K, 1L, prob = p)
}


# fast log-posterior for candidate k (pairwise ψ)
# PATCH: log_post_fast_pair()
log_post_fast_pair <- function(i, k,
                               N_blk, A_out_blk,
                               n_minus,
                               eta_sum_blk, eta_i,
                               kappa, psi_mat, alpha_vec,
                               log_kappa, log1p_tab_pair)
{
  # Dirichlet prior mass (prefers non-empty blocks when alpha0>1)
  lp <- log(alpha_vec[k] + n_minus[k])
  
  # κ–η intensity part (Gamma–Poisson, non-directional)
  lambda_vec_k <- kappa %*% eta_sum_blk  # length K
  lp <- lp - eta_i * lambda_vec_k[k] + (log_kappa %*% N_blk)[k]
  
  # Pairwise directional part (Binomial with Pólya-Gamma augmentation cached)
  K   <- length(alpha_vec); idx <- seq_len(K); mask <- idx != k
  l   <- idx[mask]
  a   <- pmin(k, l); b <- pmax(k, l)
  sgn <- ifelse(l > k, 1, -1)   # +ψ for forward (l>k), -ψ otherwise
  
  
  psi_use   <- psi_mat[cbind(a, b)]
  log1p_pos <- log1p_tab_pair$pos[cbind(a, b)]
  log1p_neg <- log1p_tab_pair$neg[cbind(a, b)]
  
  A_dir <- A_out_blk[mask]
  lp <- lp +
    sum(sgn * A_dir * psi_use) -
    sum(N_blk[mask] * ifelse(sgn == 1, log1p_pos, log1p_neg)) -
    N_blk[k] * log(2)
  
  lp
}



# Node-wise Gibbs update for z (pairwise ψ)
z_update_osbm_pair <- function(i, A, z, eta, kappa, psi, alpha_vec,
                               log1p_tab, log_kappa,
                               i_idx = NULL, j_idx = NULL, N_edge = NULL,
                               edge_by_node = NULL)
{
  # here, `psi` is a KxK matrix; `log1p_tab` is from make_log1p_tables_pair()
  K <- length(alpha_vec)
  idxK <- seq_len(K)
  
  e_out <- A[i, ];  e_in <- A[, i];  Nij <- e_out + e_in
  N_blk <- A_out_blk <- numeric(K)
  if (any(Nij > 0)) {
    ell <- z[which(Nij > 0)]
    N_blk_tmp     <- tapply(Nij[Nij > 0], ell, sum)
    A_out_blk_tmp <- tapply(e_out[Nij > 0], ell, sum)
    N_blk[as.integer(names(N_blk_tmp))]         <- N_blk_tmp
    A_out_blk[as.integer(names(A_out_blk_tmp))] <- A_out_blk_tmp
  }
  
  n_minus     <- tabulate(z[-i], nbins = K)
  eta_sum_blk <- tapply(eta[-i], factor(z[-i], levels = idxK), sum)
  eta_sum_blk[is.na(eta_sum_blk)] <- 0
  
  lp_vec <- sapply(
    idxK, log_post_fast_pair, i           = i,
    N_blk       = N_blk,
    A_out_blk   = A_out_blk,
    n_minus     = n_minus,
    eta_sum_blk = eta_sum_blk,
    eta_i       = eta[i],
    kappa       = kappa,
    psi_mat     = psi,
    alpha_vec   = alpha_vec,
    log_kappa   = log_kappa,
    log1p_tab_pair = log1p_tab
  )
  
  p <- exp(lp_vec - max(lp_vec))
  sample.int(K, 1L, prob = p)
}
# Cumulative-sum operator mapping δ -> ψ
# ψ_1 = δ_1, ψ_r = ψ_{r-1} + δ_r  for r=2..(K-1)
make_T <- function(K){
  D <- K - 1
  if (D <= 0) return(matrix(0,0,0))
  T <- matrix(0, nrow = D, ncol = D)
  for (r in 1:D) T[r, 1:r] <- 1
  T
}


osbm_relabel <- function(mcmc_out,
                         regime = c("WST","SST"),
                         tie_break = c("size","label")) {
  regime    <- match.arg(regime)
  tie_break <- match.arg(tie_break)
  
  # expected fields
  Z  <- mcmc_out$z          # S x n
  KAP<- mcmc_out$kappa      # S x K x K
  PSI<- mcmc_out$psi        # WST: S x K x K (upper-tri meaningful); SST: S x (K-1)
  stopifnot(is.matrix(Z), length(dim(KAP))==3L, nrow(Z)==dim(KAP)[1])
  
  S <- nrow(Z)
  n <- ncol(Z)
  K <- dim(KAP)[2]
  
  relab_z   <- matrix(NA_integer_, S, n)
  relab_kap <- array(NA_real_, dim = c(S, K, K))
  relab_psi <- PSI # keep same shape
  perms     <- matrix(NA_integer_, S, K)
  
  if (regime == "SST") {
    # Intrinsic order = index order; only possible symmetry is global reversal.
    for (s in seq_len(S)) {
      z_s   <- as.integer(Z[s,])
      kap_s <- KAP[s,,]
      psi_s <- as.numeric(PSI[s,])
      
      if (sum(psi_s) < 0) {
        # enforce increasing distance-bias orientation
        p_old2new <- K:1
        relab_z[s,]      <- K + 1L - z_s
        relab_kap[s,,]   <- kap_s[K:1, K:1]
        relab_psi[s,]    <- -psi_s
      } else {
        p_old2new <- seq_len(K)
        relab_z[s,]      <- z_s
        relab_kap[s,,]   <- kap_s
        relab_psi[s,]    <- psi_s
      }
      perms[s,] <- p_old2new
    }
    
  } else { # WST (pairwise ψ)
    stopifnot(length(dim(PSI))==3L, dim(PSI)[2]==K, dim(PSI)[3]==K)
    
    for (s in seq_len(S)) {
      z_s   <- as.integer(Z[s,])
      kap_s <- KAP[s,,]
      psi_s <- PSI[s,,]
      
      # model-implied success probs ρ_{kℓ} from pairwise ψ:
      # for k<ℓ: ρ_{kℓ}=logit^{-1}(+ψ_{kℓ}); for k>ℓ: ρ_{kℓ}=logit^{-1}(-ψ_{ℓk})
      rho <- matrix(0.5, K, K)
      ut  <- which(upper.tri(psi_s), arr.ind = TRUE)
      if (nrow(ut)) {
        v <- psi_s[ut]
        p_up  <- 1/(1+exp(-v)) # k<l
        rho[ut] <- p_up
        # mirror for k>l with opposite sign
        rho[cbind(ut[,2], ut[,1])] <- 1 - p_up
      }
      
      # mean success per block (intrinsic strength)
      bar_rho <- rowMeans(rho) * K/(K-1)  # (since diag=0.5; factor rescales to mean over ℓ≠k)
      
      # deterministic tie-breaks
      n_k <- tabulate(z_s, nbins = K)
      if (tie_break == "size") {
        ord_old <- order(-bar_rho, -n_k, seq_len(K)) # higher ρ, larger size, then label
      } else {
        ord_old <- order(-bar_rho, seq_len(K))
      }
      
      p_old2new <- integer(K); p_old2new[ord_old] <- seq_len(K)
      
      # apply permutation
      relab_z[s,]      <- p_old2new[z_s]
      relab_kap[s,,]   <- kap_s[ord_old, ord_old, drop=FALSE]
      relab_psi[s,,]   <- psi_s[ord_old, ord_old, drop=FALSE]
      perms[s,]        <- p_old2new
    }
  }
  
  mcmc_out$z     <- relab_z
  mcmc_out$kappa <- relab_kap
  mcmc_out$psi   <- relab_psi
  mcmc_out$perm  <- perms
  attr(mcmc_out, "relabel_model") <- paste0(regime,"-OSBM")
  attr(mcmc_out, "relabel_order") <- "intrinsic"
  mcmc_out
}



make_loglik_matrix_osbm_edgewise <- function(A, fit_out) {
  A <- if (!is.matrix(A)) as.matrix(A) else A
  n <- nrow(A)
  stopifnot(ncol(A) == n)
  obs_idx <- which(row(A) != col(A), arr.ind = TRUE)
  N <- nrow(obs_idx)
  
  S <- dim(fit_out$eta)[1]
  log_lik <- matrix(NA_real_, nrow = S, ncol = N)
  
  psi_is_pair <- length(dim(fit_out$psi)) == 3
  
  for (s in 1:S) {
    z     <- fit_out$z[s, ]
    eta   <- fit_out$eta[s, ]
    kappa <- fit_out$kappa[s, , ]
    
    psi_use_vec <- if (!psi_is_pair) fit_out$psi[s, ] else NULL
    psi_use_mat <- if (psi_is_pair)  fit_out$psi[s, , ] else NULL
    
    # precompute p_ij on the fly
    for (m in 1:N) {
      i <- obs_idx[m,1]; j <- obs_idx[m,2]
      di <- z[i]; dj <- z[j]
      if (di == dj) {
        p_ij <- 0.5
      } else if (psi_is_pair) {
        a <- min(di, dj); b <- max(di, dj)
        sgn <- ifelse(di < dj, 1, -1)
        p_ij <- plogis(sgn * psi_use_mat[a, b])
      } else {
        d <- abs(di - dj)
        sgn <- ifelse(di > dj, 1, -1)
        p_ij <- plogis(sgn * psi_use_vec[d])
      }
      p_ij <- min(max(p_ij, 1e-12), 1 - 1e-12)
      
      lam_tot <- eta[i] * eta[j] * kappa[di, dj]
      lam_ij  <- lam_tot * p_ij
      
      # directed Poisson log-lik (keep zeros!)
      log_lik[s, m] <- dpois(A[i, j], lam_ij, log = TRUE)
    }
  }
  
  attr(log_lik, "obs_idx") <- obs_idx
  log_lik                   # S x N -> feed directly to loo()
}
# Apply a K-permutation p (old->new) to vector of labels
permute_labels <- function(z_vec, p) p[ z_vec ]

# Apply a K-permutation to a KxK matrix by simultaneous row/col permutation
permute_KxK <- function(M, p) M[p, p, drop = FALSE]

# Apply a K-permutation to an S x K x K array
permute_array_SK2 <- function(ARR, p) {
  ARRp <- ARR
  for (s in seq_len(dim(ARR)[1])) ARRp[s,,] <- permute_KxK(ARR[s,,], p)
  ARRp
}

# Apply a K-permutation to an S x n matrix of labels
permute_Z_SxN <- function(Z, p) {
  Zp <- Z
  for (s in seq_len(nrow(Z))) Zp[s,] <- permute_labels(Z[s,], p)
  Zp
}
# Reference partition from the chain (label-invariant)
make_reference_partition <- function(Z_SxN) {
  # SALSO VI is robust and fast; you already use it.
  as.integer(salso::salso(Z_SxN, loss = salso::VI(a = 0.5), nRuns = 1L, nCores = 1L))
}
# Returns a list with perm[s, ] for s = 1..S, and relabeled (Z, Lambda, etc.)
relabel_by_ECR <- function(mcmc_out, model = c("DCSBM","WST","SST")) {
  model <- match.arg(model)
  
  Z   <- mcmc_out$z     # S x n
  S   <- nrow(Z)
  n   <- ncol(Z)
  # Determine K from a block param present in the object
  if (model == "DCSBM") {
    LAM <- mcmc_out$lambda  # list of K_s x K_s or S x K x K
    if (is.list(LAM)) {
      K <- max(vapply(LAM, nrow, integer(1)), max(Z, na.rm = TRUE))
    } else {
      K <- dim(LAM)[2]
      if (is.null(K)) K <- max(Z, na.rm = TRUE)
    }
  } else {
    KAP <- mcmc_out$kappa   # S x K x K
    K   <- dim(KAP)[2]
  }
  
  # 1) Reference partition (label-invariant)
  z_ref <- make_reference_partition(Z)
  
  # 2) ECR permutations to match each draw z^(s) to z_ref
  #    label.switching::ecr expects a list of allocation vectors
  alloc_list <- lapply(seq_len(S), function(s) as.integer(Z[s,]))
  ecr_res <- label.switching::ecr(
    z = Z,
    K     = K,
    zpivot = as.integer(z_ref)
  )
  Perm <- ecr_res$permutations  # S x K, mapping old->new
  
  # 3) Apply permutations to Z and block-structured parameters
  Z_new <- Z
  for (s in 1:S) Z_new[s,] <- permute_labels(Z[s,], Perm[s,])
  
  if (model == "DCSBM") {
    pad_lambda <- function(lam, K) {
      if (!is.matrix(lam)) stop("lambda entries must be matrices.")
      if (nrow(lam) == K && ncol(lam) == K) return(lam)
      lam_pad <- matrix(0, K, K)
      lam_pad[seq_len(nrow(lam)), seq_len(ncol(lam))] <- lam
      lam_pad
    }
    
    if (is.list(LAM)) {
      LAM_new <- vector("list", S)
      for (s in seq_len(S)) {
        lam_s <- pad_lambda(LAM[[s]], K)
        LAM_new[[s]] <- permute_KxK(lam_s, Perm[s, ])
      }
    } else {
      LAM_new <- array(NA_real_, dim = c(S, K, K))
      for (s in seq_len(S)) {
        lam_s <- pad_lambda(LAM[s, , ], K)
        LAM_new[s, , ] <- permute_KxK(lam_s, Perm[s, ])
      }
    }
    
    # thetas are length-n and unaffected by K-permutation
    out <- mcmc_out
    out$z      <- Z_new
    out$lambda <- LAM_new
    out$perm   <- Perm
    attr(out, "relabel_model") <- "DCSBM"
    attr(out, "relabel_method") <- "ECR(z)"
    return(out)
  }
  
  # OSBM (WST/SST): relabel kappa; psi needs special care for SST orientation
  KAP_new <- array(NA_real_, dim = dim(mcmc_out$kappa))
  for (s in 1:S) KAP_new[s,,] <- permute_KxK(mcmc_out$kappa[s,,], Perm[s,])
  
  if (model == "WST") {
    PSI_new <- array(NA_real_, dim = dim(mcmc_out$psi))
    for (s in 1:S) PSI_new[s,,] <- permute_KxK(mcmc_out$psi[s,,], Perm[s,])
    out <- mcmc_out
    out$z     <- Z_new
    out$kappa <- KAP_new
    out$psi   <- PSI_new
    out$perm  <- Perm
    attr(out, "relabel_model") <- "WST-OSBM"
    attr(out, "relabel_method") <- "ECR(z)"
    return(out)
  }
  
  # SST: psi is length (K-1) monotone sequence up to global reversal
  # After ECR(z), fix the global reversal by a simple orientation check:
  PSI_new <- mcmc_out$psi  # S x (K-1)
  # Build a reference psi direction from posterior mean:
  psi_ref <- colMeans(PSI_new)
  sgn_ref <- sign(sum(psi_ref))  # >0 by convention
  
  for (s in 1:S) {
    # no KxK perm needed for psi vector; labels have fixed positions after ECR
    # Enforce same global sign as reference (reversal symmetry)
    if (sign(sum(PSI_new[s,])) != sgn_ref) {
      PSI_new[s,] <- -PSI_new[s,]
      # also reverse labels 1..K simultaneously to keep consistency
      revp <- K:1
      Perm[s,]    <- Perm[s, revp]         # compose with reversal
      Z_new[s,]   <- permute_labels(Z_new[s,], revp)
      KAP_new[s,,] <- permute_KxK(KAP_new[s,,], revp)
    }
  }
  out <- mcmc_out
  out$z     <- Z_new
  out$kappa <- KAP_new
  out$psi   <- PSI_new
  out$perm  <- Perm
  attr(out, "relabel_model") <- "SST-OSBM"
  attr(out, "relabel_method") <- "ECR(z)+orientation"
  out
}

dcsbm_relabel <- function(mcmc_out,
                          eps = 1e-12,
                          tie_break = c("size","label")) {
  tie_break <- match.arg(tie_break)
  
  Z   <- mcmc_out$z        # S x n
  LAM <- mcmc_out$lambda   # S x K x K
  stopifnot(is.matrix(Z), length(dim(LAM))==3L, nrow(Z)==dim(LAM)[1])
  
  S <- nrow(Z); n <- ncol(Z); K <- dim(LAM)[2]
  
  relab_z  <- matrix(NA_integer_, S, n)
  relab_la <- array(NA_real_, dim = c(S, K, K))
  perms    <- matrix(NA_integer_, S, K)
  
  for (s in seq_len(S)) {
    z_s <- as.integer(Z[s,])
    La  <- LAM[s,,]
    
    # ρ_{kℓ} = λ_{kℓ}/(λ_{kℓ}+λ_{ℓk})   (intrinsic, block-level)
    denom <- La + t(La) + eps
    rho   <- La / denom
    diag(rho) <- 0.5
    
    bar_rho <- rowMeans(rho) * K/(K-1)
    
    n_k <- tabulate(z_s, nbins = K)
    if (tie_break == "size") {
      ord_old <- order(-bar_rho, -n_k, seq_len(K))
    } else {
      ord_old <- order(-bar_rho, seq_len(K))
    }
    p_old2new <- integer(K); p_old2new[ord_old] <- seq_len(K)
    
    relab_z[s,]    <- p_old2new[z_s]
    relab_la[s,,]  <- La[ord_old, ord_old, drop=FALSE]
    perms[s,]      <- p_old2new
  }
  
  mcmc_out$z      <- relab_z
  mcmc_out$lambda <- relab_la
  mcmc_out$perm   <- perms
  attr(mcmc_out, "relabel_model") <- "DC-SBM"
  attr(mcmc_out, "relabel_order") <- "intrinsic"
  mcmc_out
}

upper_pairs <- function(n) which(matrix(TRUE,n,n)&upper.tri(diag(n)), arr.ind = TRUE)


#' Relabel SST draws to a common orientation and flag reversals
#'
#' SST has a global reversal symmetry: psi -> -psi, labels k -> K+1-k.
#' This function orients each draw to a common reference direction and
#' applies the corresponding label reversal when needed.
#'
#' @param mcmc_out list with components:
#'   - z    : S x n integer matrix of block labels per draw
#'   - psi  : S x (K-1) numeric matrix of distance parameters
#'   - kappa: (optional) S x K x K array of symmetric Poisson rates
#'   - (optionally) any other KxK arrays to be permuted can be passed via `extras`
#' @param enforce_monotone logical; if TRUE, projects each psi^(s) onto the
#'   nondecreasing cone (isotonic) after orientation (defaults to FALSE).
#' @param extras optional named list of S x K x K arrays in `mcmc_out` to be
#'   permuted together with kappa  (e.g. lambda, etc.). Each will be looked up
#'   by name in mcmc_out and permuted in place if present.
#' @return The input list with relabelled components, plus:
#'   - $reversed : logical length-S vector indicating which draws were flipped
#'   - $perm     : S x K integer matrix of permutations (old -> new labels)
#'   Attributes set:
#'   - attr(., "relabel_model")  = "SST-OSBM"
#'   - attr(., "relabel_method") = "orientation(+isotonic?)"
#'
relabel_sst_with_reversal <- function(mcmc_out)
{
  stopifnot(is.matrix(mcmc_out$z), is.matrix(mcmc_out$psi))
  Z   <- mcmc_out$z           # S x n
  PSI <- mcmc_out$psi         # S x (K-1)
  S   <- nrow(Z)
  n   <- ncol(Z)
  D   <- ncol(PSI)
  K   <- D + 1L
  if (!is.null(mcmc_out$kappa)) {
    stopifnot(length(dim(mcmc_out$kappa)) == 3L,
              dim(mcmc_out$kappa)[1] == S,
              dim(mcmc_out$kappa)[2] == K,
              dim(mcmc_out$kappa)[3] == K)
  }
  
  # Determine a reference direction from the posterior mean psi
  psi_ref <- colMeans(PSI, na.rm = TRUE)
  sgn_ref <- if (sum(psi_ref, na.rm = TRUE) >= 0) +1L else -1L
  
  # Helpers
  permute_labels <- function(z_vec, p) p[z_vec]
  permute_KxK    <- function(M, p)     M[p, p, drop = FALSE]
  
  # Storage
  reversed <- logical(S)
  Perm     <- matrix(NA_integer_, nrow = S, ncol = K)
  
  # Work copies
  Z_new   <- Z
  PSI_new <- PSI
  if (!is.null(mcmc_out$kappa)) {
    KAP_new <- mcmc_out$kappa
  } else {
    KAP_new <- NULL
  }
  

  # Main loop
  for (s in seq_len(S)) {
    psi_s <- as.numeric(PSI[s, ])
    z_s   <- as.integer(Z[s, ])
    
    # Decide orientation for this draw
    sgn_s <- if (sum(psi_s, na.rm = TRUE) >= 0) +1L else -1L
    do_flip <- (sgn_s != sgn_ref)
    
    if (do_flip) {
      reversed[s] <- TRUE
      # Reverse labels and parameters
      revp <- K:1
      Z_new[s, ]  <- permute_labels(z_s, revp)
      PSI_new[s,] <- -psi_s
      if (!is.null(KAP_new)) KAP_new[s,,] <- permute_KxK(KAP_new[s,,], revp)
      
      # Record permutation mapping: old->new
      P <- integer(K); P[revp] <- seq_len(K)
      Perm[s, ] <- P
    } else {
      # Keep as is; permutation is identity
      Z_new[s, ]  <- z_s
      PSI_new[s,] <- psi_s
      if (!is.null(KAP_new)) KAP_new[s,,] <- KAP_new[s,,]  # no-op
      
      Perm[s, ] <- seq_len(K)
    }
  }
  
  out <- mcmc_out
  out$z    <- Z_new
  out$psi  <- PSI_new
  if (!is.null(KAP_new)) out$kappa <- KAP_new
  out$reversed <- reversed
  out$perm     <- Perm
  attr(out, "relabel_model")  <- "SST-OSBM"
  attr(out, "relabel_method") <- "orientation"
  out
}


## ---------------------------
## DC-SBM: S x M joint log-lik
## ---------------------------
loglik_matrix_dcsbm <- function(A, dcsbm_out, dyad_index) {
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  stopifnot(is.matrix(dyad_index), ncol(dyad_index) == 2)
  stopifnot(all(dyad_index[,1] < dyad_index[,2]))     # upper-tri pairs

  Z   <- dcsbm_out$z           # either S x n matrix or list of length S
  THO <- dcsbm_out$theta_out   # typically S x n matrix, but we allow list
  THI <- dcsbm_out$theta_in    # same
  if (is.null(THO) || is.null(THI)) {
    TH <- dcsbm_out$theta
    if (is.null(TH)) {
      stop("dcsbm_out must provide theta_out/theta_in or a single theta matrix.")
    }
    THO <- TH
    THI <- TH
  }
  LAM <- dcsbm_out$lambda      # list of K_s x K_s, or 3D array S x K x K

  ## --- determine S and n from Z ---
  if (is.matrix(Z)) {
    S <- nrow(Z)
    n <- ncol(Z)
  } else if (is.list(Z)) {
    S <- length(Z)
    n <- length(Z[[1]])
  } else {
    stop("dcsbm_out$z must be either a matrix (S x n) or a list of length S.")
  }

  ## --- normalise lambda representation to: list of matrices ---
  if (is.list(LAM)) {
    if (length(LAM) != S) {
      stop("When lambda is a list, length(lambda) must equal the number of iterations inferred from z.")
    }
  } else if (length(dim(LAM)) == 3L) {
    # backward-compatible case: S x K x K array
    if (dim(LAM)[1] != S) {
      stop("First dimension of lambda array must match the number of iterations inferred from z.")
    }
    LAM <- lapply(seq_len(S), function(s) LAM[s, , ])
  } else {
    stop("dcsbm_out$lambda must be either a list of matrices or a 3D array (S x K x K).")
  }

  ## --- basic dimension checks ---
  stopifnot(nrow(A) == n)
  M <- nrow(dyad_index)
  i <- dyad_index[,1]; j <- dyad_index[,2]

  ## helper to pull theta for iteration s
  get_theta_row <- function(TH, s, n) {
    if (is.matrix(TH)) {
      if (ncol(TH) != n) {
        stop("theta_out/theta_in matrix must have n columns equal to the number of nodes.")
      }
      TH[s, ]
    } else if (is.list(TH)) {
      TH[[s]]
    } else {
      stop("theta_out/theta_in must be either a matrix (S x n) or a list of length S.")
    }
  }

  ## --- log-likelihood matrix ---
  LL <- matrix(NA_real_, nrow = S, ncol = M)

  for (s in seq_len(S)) {
    z_s <- if (is.matrix(Z)) Z[s, ] else Z[[s]]
    if (length(z_s) != n) {
      stop("Length of z at iteration ", s, " does not match n.")
    }

    lam_s <- LAM[[s]]
    if (!is.matrix(lam_s)) {
      stop("lambda[[", s, "]] is not a matrix.")
    }

    K_s <- max(z_s)
    if (!all(z_s %in% seq_len(K_s))) {
      stop("z labels at iteration ", s, " must be in {1, ..., K_s}.")
    }
    if (any(dim(lam_s) < K_s)) {
      stop("Dimensions of lambda[[", s, "]] are smaller than K_s x K_s implied by z.")
    }

    tho <- get_theta_row(THO, s, n)
    thi <- get_theta_row(THI, s, n)

    mu_ij <- tho[i] * thi[j] * lam_s[cbind(z_s[i], z_s[j])]
    mu_ji <- tho[j] * thi[i] * lam_s[cbind(z_s[j], z_s[i])]

    LL[s, ] <- dpois(A[cbind(i, j)], mu_ij, log = TRUE) +
      dpois(A[cbind(j, i)], mu_ji, log = TRUE)
  }

  LL
}

## ---------------------------------------
## WST: S x M joint log-lik on upper dyads
## ---------------------------------------
loglik_matrix_wst <- function(A, mcmc_out, dyad_index) {
  A  <- as.matrix(A)
  AT <- t(A)
  stopifnot(nrow(A) == ncol(A))
  stopifnot(is.matrix(dyad_index), ncol(dyad_index) == 2)
  stopifnot(all(dyad_index[,1] < dyad_index[,2]))
  
  Z    <- mcmc_out$z       # S x n
  ETA  <- mcmc_out$eta     # S x n
  KAP  <- mcmc_out$kappa   # S x K x K   (rate for total N_ij)
  PSI  <- mcmc_out$psi     # S x K x K   (pairwise "skill" > 0 off-diagonal)
  
  S <- nrow(Z)
  M <- nrow(dyad_index)
  i <- dyad_index[,1]; j <- dyad_index[,2]
  
  LL <- matrix(NA_real_, nrow = S, ncol = M)
  for (s in seq_len(S)) {
    z    <- Z[s,]; eta <- ETA[s,]
    kapp <- KAP[s,,]; psi <- PSI[s,,]
    
    Nij <- A[cbind(i, j)] + AT[cbind(i, j)]
    lam <- eta[i] * eta[j] * kapp[cbind(z[i], z[j])]
    
    same <- (z[i] == z[j])
    p    <- numeric(M)
    p[same]    <- 0.5
    if (any(!same)) {
      sgn        <- sign(z[j] - z[i])
      a      <- pmin(z[i], z[j])
      b      <- pmax(z[i], z[j])
      psi_ab <- psi[cbind(a, b)]
      p[!same] <- plogis(sgn[!same] * psi_ab[!same])
      
    }
    # Numerical guard (paranoia); dbinom handles 0/1 properly anyway
    eps <- 1e-12; p <- pmin(pmax(p, eps), 1 - eps)
    
    LL[s, ] <- dpois(Nij, lam, log = TRUE) + dbinom(A[cbind(i, j)], Nij, p, log = TRUE)
  }
  LL
}


relabel_osbm <- function(mcmc_out,
                         A,
                         ordering = c("WST","SST","NONE"),
                         score = c("success","outdeg","indeg","netdeg")) {
  ordering <- match.arg(ordering)
  score    <- match.arg(score)
  
  if (ordering != "WST") {
    # SST (distance psi_d) or NONE: do nothing
    mcmc_out$perm <- matrix(rep(seq_len(dim(mcmc_out$kappa)[2]), each = nrow(mcmc_out$z)),
                            nrow = nrow(mcmc_out$z))
    return(mcmc_out)
  }
  
  # --- WST relabelling by block score (A-driven) -----------------------
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  n <- nrow(A)
  
  # per-node metrics (fixed, from data)
  outdeg   <- rowSums(A)
  indeg    <- colSums(A)
  netdeg   <- outdeg - indeg
  matches  <- A + t(A)
  success <- outdeg / pmax(rowSums(matches), 1)  # row-wise success rate
  
  node_metric <- switch(score,
                        outdeg  = outdeg,
                        indeg   = indeg,
                        netdeg  = netdeg,
                        success = success)
  
  draws_z   <- mcmc_out$z                  # S x n
  draws_kap <- mcmc_out$kappa              # S x K x K
  draws_psi <- mcmc_out$psi                # WST: S x K x K (magnitudes)
  S         <- nrow(draws_z)
  K         <- dim(draws_kap)[2]
  
  relab_z   <- draws_z
  relab_kap <- array(NA_real_, dim = dim(draws_kap))
  relab_psi <- array(NA_real_, dim = dim(draws_psi))
  perms     <- matrix(NA_integer_, nrow = S, ncol = K)
  
  for (s in seq_len(S)) {
    z_s <- as.integer(draws_z[s,])
    
    # block means; empty blocks -> -Inf so they go to the tail
    means <- rep(-Inf, K)
    for (k in seq_len(K)) {
      idx <- which(z_s == k)
      if (length(idx) > 0) means[k] <- mean(node_metric[idx])
    }
    
    # ord_old: old labels sorted by decreasing block score (strongest first)
    ord_old <- order(means, decreasing = TRUE, na.last = TRUE)
    # mapping old -> new (p[old] = new)
    p_old2new <- integer(K); p_old2new[ord_old] <- seq_len(K)
    
    # apply to z
    relab_z[s, ] <- p_old2new[z_s]
    
    # apply to kappa and psi by reindexing rows/cols using ord_old
    relab_kap[s, , ] <- draws_kap[s, ord_old, ord_old, drop = FALSE]
    relab_psi[s, , ] <- draws_psi[s, ord_old, ord_old, drop = FALSE]
    
    perms[s, ] <- p_old2new
  }
  
  mcmc_out$z     <- relab_z
  mcmc_out$kappa <- relab_kap
  mcmc_out$psi   <- relab_psi
  mcmc_out$perm  <- perms
  attr(mcmc_out, "relabel_ordering") <- "WST"
  attr(mcmc_out, "relabel_score")    <- score
  mcmc_out
}


## ---------------------------------------
## SST: S x M joint log-lik on upper dyads
## (psi is length K-1 by distance |k-l|)
## ---------------------------------------
loglik_matrix_sst <- function(A, mcmc_out, dyad_index) {
  A  <- as.matrix(A)
  AT <- t(A)
  stopifnot(nrow(A) == ncol(A))
  stopifnot(is.matrix(dyad_index), ncol(dyad_index) == 2)
  stopifnot(all(dyad_index[,1] < dyad_index[,2]))
  
  Z    <- mcmc_out$z      # S x n
  ETA  <- mcmc_out$eta    # S x n
  KAP  <- mcmc_out$kappa  # S x K x K   (rate for total N_ij)
  PSI  <- mcmc_out$psi    # S x (K-1)   (by block-distance)
  
  S <- nrow(Z)
  M <- nrow(dyad_index)
  i <- dyad_index[,1]; j <- dyad_index[,2]
  
  LL <- matrix(NA_real_, nrow = S, ncol = M)
  for (s in seq_len(S)) {
    z    <- Z[s,]; eta <- ETA[s,]
    kapp <- KAP[s,,]; psi <- as.numeric(PSI[s,])  # length K-1
    
    Nij <- A[cbind(i, j)] + AT[cbind(i, j)]
    lam <- eta[i] * eta[j] * kapp[cbind(z[i], z[j])]
    
    same <- (z[i] == z[j])
    p    <- numeric(M)
    p[same] <- 0.5
    if (any(!same)) {
      sgn <- sign(z[j] - z[i])
      d   <- abs(z[i] - z[j])
      psi_d <- psi[d[!same]]  # distance-indexed
      p[!same] <- plogis(sgn[!same] * psi_d)
    }
    eps <- 1e-12; p <- pmin(pmax(p, eps), 1 - eps)
    
    LL[s, ] <- dpois(Nij, lam, log = TRUE) + dbinom(A[cbind(i, j)], Nij, p, log = TRUE)
  }
  LL
}



estimate_partition <- function(z_trace) {
  est_vi_K <- salso::salso(z_trace, loss= salso::VI(a = 0.5), nRuns=1, nCores=1)
  as.numeric(est_vi_K)
}



update_psi_sst <- function(K, bar_y, bar_omega, psi_curr,
                           mu0, sig2_0, tau2_0,
                           n_inner_sweeps = 2)
{
  D <- K - 1
  if (D == 0) return(numeric(0))
  
  # Build T, Ω, Q, g, and prior precision/shift
  Tmat  <- make_T(K)
  Omega <- diag(bar_omega, nrow = D, ncol = D)
  Q     <- crossprod(Tmat, Omega %*% Tmat)   # T^T Ω T
  g     <- crossprod(Tmat, bar_y)            # T^T \bar y   (note: no Ω)
  
  Vinv  <- diag(c(1/sig2_0, rep(1/tau2_0, max(0, D-1))), D, D)  # diag precision
  muvec <- c(mu0, rep(0, max(0, D-1)))
  
  Qstar <- Q + Vinv
  gstar <- g + Vinv %*% muvec
  
  # Start at current δ (invert the cumulative sum): δ0 = ψ1, δr = ψ_{r+1} - ψ_r
  delta <- numeric(D)
  if (length(psi_curr) == D) {
    delta[1] <- max(psi_curr[1], 1e-10)
    if (D >= 2) {
      incs <- pmax(psi_curr[-1] - psi_curr[-D], 1e-10)
      delta[2:D] <- incs
    }
  } else {
    # fall back if no current ψ provided
    delta[] <- 0.1
  }
  
  # Coordinate-wise truncated-Normal updates on δ_r > 0
  for (s in seq_len(n_inner_sweeps)) {
    for (r in seq_len(D)) {
      a_rr <- Qstar[r, r]
      b_r  <- as.numeric(gstar[r] - sum(Qstar[r, -r, drop = FALSE] * delta[-r]))
      mean_r <- b_r / a_rr
      sd_r   <- 1 / sqrt(a_rr)
      delta[r] <- truncnorm::rtruncnorm(1, a = 0, b = Inf, mean = mean_r, sd = sd_r)
      if (!is.finite(delta[r]) || delta[r] <= 0) delta[r] <- 1e-10
    }
  }
  
  # Map back to ψ = T δ (monotone & positive by construction)
  as.numeric(Tmat %*% delta)
}



# ===== Minimal assertion helpers (no extra packages) =========================
.assert_scalar_int <- function(x, name) {
  if (length(x) != 1L || is.na(x) || x != as.integer(x) || x < 1L)
    stop(sprintf("`%s` must be a positive integer scalar.", name), call. = FALSE)
  invisible(TRUE)
}
.assert_logical_scalar <- function(x, name){
  if (length(x) != 1L || is.na(x) || !is.logical(x))
    stop(sprintf("`%s` must be TRUE/FALSE.", name), call. = FALSE)
  invisible(TRUE)
}
.assert_len <- function(x, len, name){
  if (length(x) != len) stop(sprintf("`%s` must have length %d.", name, len), call. = FALSE)
  invisible(TRUE)
}
.assert_matrix <- function(A, n, name){
  if (!is.matrix(A) && !inherits(A,"Matrix"))
    stop(sprintf("`%s` must be a matrix or Matrix.", name), call. = FALSE)
  if (nrow(A) != n || ncol(A) != n)
    stop(sprintf("`%s` must be square %dx%d.", name, n, n), call. = FALSE)
  if (anyNA(A)) stop(sprintf("`%s` contains NA.", name), call. = FALSE)
  invisible(TRUE)
}
.assert_truth <- function(truth, needs){
  if (is.null(truth) || is.na(truth)) stop("`truth` must be a named list when parameters are fixed.", call. = FALSE)
  miss <- setdiff(needs, names(truth))
  if (length(miss)) stop("`truth` is missing: ", paste(miss, collapse=", "), call. = FALSE)
  invisible(TRUE)
}
.safe_log <- function(x) log(pmax(x, 1e-15))


# Funzione che decide se aggiungere ± sd oppure no
pretty_metric_summary <- function(x, name) {
  if (length(unique(x)) == 1) return(sprintf("%.3f", x[1]))
  if (name %in% c("MAE_eta", "MAE_kappa", "MAE_psi", "VI")) {
    return(sprintf("%.3f ± %.3f", mean(x), sd(x)))
  } else {
    return(sprintf("%.3f", mean(x)))  # ARI, coverage etc.
  }
}
