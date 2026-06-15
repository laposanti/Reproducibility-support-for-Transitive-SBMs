## ---- metrics + CI table for observed vs WST/SST ----------------------------
.validate_params <- function(n, K, z, eta, kappa, psi_scale,
                             align_z = c("resample","truncate","error")) {
  align_z <- match.arg(align_z)
  stopifnot(length(eta) == n)
  if (!is.matrix(kappa) || any(dim(kappa) != c(K,K))) {
    stop("kappa_t must be a K x K matrix.")
  }
  if (any(!is.finite(eta)))   stop("eta_t has non-finite values.")
  if (any(!is.finite(kappa))) stop("kappa_t has non-finite values.")
  if (!is.numeric(psi_scale) || length(psi_scale)!=1L || !is.finite(psi_scale)) {
    stop("psi_scale must be a finite scalar.")
  }
  
  # z labels must be in 1:K and length n
  if (length(z) != n) {
    if (align_z == "resample") {
      z <- sample(z, n, replace = TRUE)
    } else if (align_z == "truncate") {
      z <- z[seq_len(n)]
    } else {
      stop(sprintf("length(z_t)=%d differs from n=%d. Fix z_t or set align_z.", length(z), n))
    }
  }
  if (any(!is.finite(z))) stop("z_t has non-finite entries.")
  if (any(z < 1 | z > K)) stop("z_t contains labels outside 1:K.")
  
  z
}
# Compute simple directed graph metrics from an adjacency matrix.
# By default we binarize counts at > 0.
compute_metrics <- function(A, binarize = TRUE) {
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  n <- nrow(A)
  
  B <- if (binarize) (A > 0L) * 1L else A
  diag(B) <- 0L
  
  # For directed graphs w/out self loops:
  m <- sum(B)                          # number of arcs
  possible <- n * (n - 1)              # arcs possible
  density <- if (possible > 0) m / possible else NA_real_
  sparsity <- 1 - density
  
  outdeg <- rowSums(B)
  indeg  <- colSums(B)
  mean_outdeg <- mean(outdeg)
  mean_indeg  <- mean(indeg)
  
  # Reciprocity: fraction of arcs that are mutual
  # (each mutual dyad contributes 2 to the numerator)
  mutual_arcs <- sum(B * t(B))
  reciprocity <- if (m > 0) mutual_arcs / m else NA_real_
  
  data.frame(
    n = n,
    edges = m,
    density = density,
    sparsity = sparsity,
    mean_outdeg = mean_outdeg,
    mean_indeg  = mean_indeg,
    reciprocity = reciprocity,
    check.names = FALSE
  )
}

# Summarize a numeric vector into mean and central credible interval
cred_summ <- function(x, probs = c(0.025, 0.975)) {
  q <- stats::quantile(x, probs = probs, names = FALSE, type = 8)
  c(mean = mean(x), lwr = q[1], upr = q[2])
}

# Simulate nsample adjacency matrices under a regime and compute metrics each time.
# Assumes your functions and parameters exist: make_psi_wst/make_psi_sst, simulate_osbm.
simulate_metric_summaries <- function(regime = c("WST","SST"),
                                      nsample,
                                      n, K, z_t, eta_t, kappa_t,
                                      psi_scale, align_z = "resample") {
  regime <- match.arg(regime)
  
  dens <- spars <- mout <- min <- recip <- numeric(nsample)
  edges <- integer(nsample)
  # align/validate z and parameters
  z_use <- .validate_params(n, K, z_t, eta_t, kappa_t, psi_scale, align_z = align_z)
  
  for (s in seq_len(nsample)) {
    if (regime == "WST") {
      psi_t <- make_psi_wst(K, scale = psi_scale)
      A_s   <- simulate_osbm(n, K, z_use, eta_t, kappa_t, psi = psi_t, regime = "WST")
    } else {
      psi_t <- make_psi_sst(K, scale = psi_scale)
      A_s   <- simulate_osbm(n, K, z_use, eta_t, kappa_t, psi = psi_t, regime = "SST")
    }
    
    M <- compute_metrics(A_s)
    edges[s] <- M$edges
    dens[s]  <- M$density
    spars[s] <- M$sparsity
    mout[s]  <- M$mean_outdeg
    min[s]   <- M$mean_indeg
    recip[s] <- M$reciprocity
  }
  
  list(
    n             = n,
    edges         = cred_summ(edges),
    density       = cred_summ(dens),
    sparsity      = cred_summ(spars),
    mean_outdeg   = cred_summ(mout),
    mean_indeg    = cred_summ(min),
    reciprocity   = cred_summ(recip)
  )
}

# Format a mean and CI as a compact string: "m (l,u)"
fmt_ci <- function(mean, lwr, upr, digits = 3) {
  sprintf("%.*f [%.*f, %.*f]", digits, mean, digits, lwr, digits, upr)
}

## --- (Keep/Source the earlier helpers) --------------------------------------
## compute_metrics(), cred_summ(), fmt_ci(), 
## .validate_params(), .safe_simulate(), simulate_metric_summaries()

# If you don't have them in scope, paste them above or source your previous file.

## --- Build a single table with 3 observed datasets + 2 model rows -----------

build_all_tables <- function(datasets = c("mountain_goats","citations_data","macaques_data"),
                             # fixed n for simulations (independent of application n)
                             n_sim,
                             # model parameters
                             K, z_t, eta_t, kappa_t, psi_scale,
                             nsample = 500, seed = 123,
                             align_z = c("resample","truncate","error"),
                             model_dataset_label = "(fixed-n sim)") {
  align_z <- match.arg(align_z)
  set.seed(seed)
  
  # 1) Observed rows for requested datasets
  obs_rows <- lapply(datasets, function(ds) {
    A_obs <- choose_dataset(ds)
    M     <- compute_metrics(A_obs)
    data.frame(
      model = "Observed",
      dataset = ds,
      n = M$n,
      edges = M$edges,
      density = M$density,
      sparsity = M$sparsity,
      mean_outdeg = M$mean_outdeg,
      mean_indeg  = M$mean_indeg,
      reciprocity = M$reciprocity,
      check.names = FALSE
    )
  })
  obs_tbl <- do.call(rbind, obs_rows)
  rownames(obs_tbl) <- NULL
  
  # 2) Simulation rows at fixed n_sim (same for both regimes)
  #    Validate/align z once at the simulation n
  z_use <- .validate_params(n_sim, K, z_t, eta_t, kappa_t, psi_scale, align_z = align_z)
  
  wst <- simulate_metric_summaries("WST", nsample, n_sim, K, z_use, eta_t, kappa_t, psi_scale)
  sst <- simulate_metric_summaries("SST", nsample, n_sim, K, z_use, eta_t, kappa_t, psi_scale)
  
  wst_row <- data.frame(
    model = "WST",
    dataset = model_dataset_label,
    n = wst$n,
    edges = fmt_ci(wst$edges["mean"],   wst$edges["lwr"],   wst$edges["upr"]),
    density = fmt_ci(wst$density["mean"], wst$density["lwr"], wst$density["upr"]),
    sparsity = fmt_ci(wst$sparsity["mean"], wst$sparsity["lwr"], wst$sparsity["upr"]),
    mean_outdeg = fmt_ci(wst$mean_outdeg["mean"], wst$mean_outdeg["lwr"], wst$mean_outdeg["upr"]),
    mean_indeg  = fmt_ci(wst$mean_indeg["mean"],  wst$mean_indeg["lwr"],  wst$mean_indeg["upr"]),
    reciprocity = fmt_ci(wst$reciprocity["mean"], wst$reciprocity["lwr"], wst$reciprocity["upr"]),
    check.names = FALSE
  )
  
  sst_row <- data.frame(
    model = "SST",
    dataset = model_dataset_label,
    n = sst$n,
    edges = fmt_ci(sst$edges["mean"],   sst$edges["lwr"],   sst$edges["upr"]),
    density = fmt_ci(sst$density["mean"], sst$density["lwr"], sst$density["upr"]),
    sparsity = fmt_ci(sst$sparsity["mean"], sst$sparsity["lwr"], sst$sparsity["upr"]),
    mean_outdeg = fmt_ci(sst$mean_outdeg["mean"], sst$mean_outdeg["lwr"], sst$mean_outdeg["upr"]),
    mean_indeg  = fmt_ci(sst$mean_indeg["mean"],  sst$mean_indeg["lwr"],  sst$mean_indeg["upr"]),
    reciprocity = fmt_ci(sst$reciprocity["mean"], sst$reciprocity["lwr"], sst$reciprocity["upr"]),
    check.names = FALSE
  )
  
  # 3) Bind: 3 observed + 2 simulated rows
  out <- rbind(obs_tbl, wst_row, sst_row)
  rownames(out) <- NULL
  out
}

## --- Example usage -----------------------------------------------------------
K <- 4
n_sim <- 100     

# <-- your fixed simulation n
# ----- Generate truth & data -----
z_t   <- rep(1:K, length.out = n_sim)
eta_t <- runif(n_sim, 0.8, 1.2); eta_t <- n_sim * eta_t / sum(eta_t)
kappa_t <- make_kappa(K, rate_scale = 1)


tab_all <- build_all_tables(
  datasets = c("mountain_goats","citations_data","macaques_data"),
  n_sim    = n_sim,
  K        = K,
  z_t      = z_t,
  eta_t    = eta_t,
  kappa_t  = kappa_t,
  psi_scale = pars$psi_scale,
  nsample      = 50,
  seed      = 42,
  align_z   = "resample",               # if z_t length != n_sim, will be handled
  model_dataset_label = "(fixed-n sim)" # label for the two model rows
)
print(tab_all)


tab_all%>%
  kable(format = 'latex', booktabs = T,digits = 3)
