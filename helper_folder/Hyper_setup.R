### Hypers setup


#------------------------------------------------------------
# Helper: mean dyadic counts (N_ij = A_ij + A_ji)
#------------------------------------------------------------
compute_mean_edges <- function(A) {
  N_mat <- A + t(A)
  diag(N_mat) <- 0
  m <- mean(N_mat[upper.tri(N_mat)])
  if (!is.finite(m) || m <= 0) m <- 1e-6
  m
}

#------------------------------------------------------------
# Choose a_kappa from prior pseudo-counts
#
# K_target          : "typical" number of blocks you expect
# pseudo_factor     : how many prior pseudo-counts relative to a
#                     typical node–block count R_typ.
#                     - 0.5  => prior ~ half as informative as R_typ
#                     - 1    => prior ~ as informative as R_typ
#                     - 2    => prior ~ twice as informative as R_typ
# E_eta             : prior mean of eta (usually 1)
#------------------------------------------------------------
choose_a_kappa_from_pseudocount <- function(
    A,
    K_target      = 5,
    pseudo_factor = 0.5,
    E_eta         = 1
) {
  n <- nrow(A)
  if (K_target <= 0 || K_target > n) {
    stop("K_target must be in (0, n].")
  }
  
  mean_edges <- compute_mean_edges(A)
  
  # Typical block size and typical node–block total count
  m_typ  <- n / K_target               # typical block size
  R_typ  <- m_typ * mean_edges * E_eta^2
  
  # Prior pseudo-counts: a_kappa
  a_kappa <- pseudo_factor * R_typ
  
  list(
    a_kappa   = a_kappa,
    b_kappa   = a_kappa / mean_edges,
    R_typ     = R_typ,
    mean_edges = mean_edges
  )
}

#------------------------------------------------------------
# Gnedin prior utilities to set gamma from K_expected
#------------------------------------------------------------
gn_prior_K <- function(n, gamma,
                       ordering_prior_mode = c("equivalence_class", "labelled")) {
  ordering_prior_mode <- match.arg(ordering_prior_mode)
  k <- 1:n
  log_p <- lchoose(n, k) +
    (k - 1) * log1p(-gamma) +
    (n - k) * log(gamma) -
    (n - 1) * log1p(gamma)

  # Order-aware adjustment:
  # - equivalence_class: apply GN EPPF to partition classes (default)
  # - labelled: treat ordered labellings as distinct => multiplicity K!
  if (ordering_prior_mode == "labelled") {
    log_p <- log_p + lgamma(k + 1)
  }

  p <- exp(log_p - max(log_p))
  p <- p / sum(p)
  list(k = k, p = p)
}

expected_K_gn <- function(n, gamma,
                          ordering_prior_mode = c("equivalence_class", "labelled")) {
  ordering_prior_mode <- match.arg(ordering_prior_mode)
  kp <- gn_prior_K(n, gamma, ordering_prior_mode = ordering_prior_mode)
  sum(kp$k * kp$p)
}

gnedin_ordering_discrepancy <- function(n, gamma) {
  ek_equiv <- expected_K_gn(n, gamma, ordering_prior_mode = "equivalence_class")
  ek_label <- expected_K_gn(n, gamma, ordering_prior_mode = "labelled")
  rel_gap <- if (ek_equiv > 0) (ek_label - ek_equiv) / ek_equiv else NA_real_
  list(
    E_K_equivalence_class = ek_equiv,
    E_K_labelled = ek_label,
    relative_gap = rel_gap
  )
}

choose_gamma_from_K_expected <- function(
    n,
    K_expected,
    ordering_prior_mode = c("equivalence_class", "labelled"),
    lower = 1e-4,
    upper = 1 - 1e-4
) {
  ordering_prior_mode <- match.arg(ordering_prior_mode)
  f <- function(gamma) {
    expected_K_gn(n, gamma, ordering_prior_mode = ordering_prior_mode) - K_expected
  }
  f_low  <- f(lower)
  f_high <- f(upper)
  if (!is.finite(f_low) || !is.finite(f_high) || f_low * f_high > 0) {
    warning("Could not bracket a root for gamma; using gamma = 0.5 as fallback.")
    return(0.5)
  }
  uniroot(f, interval = c(lower, upper))$root
}

#------------------------------------------------------------
# Main: new get_principled_hypers()
#   - uses pseudo-count calibration for a_kappa
#   - uses Gnedin gamma to encode K_expected
#------------------------------------------------------------
get_principled_hypers <- function(
    A,
    K_expected      = 5,     # prior target for K
    pseudo_factor   = 0.5,   # prior pseudo-counts vs typical node-block count
    E_eta           = 1,     # prior mean of eta
    ordering_prior_mode = c("equivalence_class", "labelled"),
    dominance_target = 0.95  # reserved for psi calibration
) {
  ordering_prior_mode <- match.arg(ordering_prior_mode)
  n <- nrow(A)
  
  # 1. a_kappa & b_kappa from pseudo-count argument
  vol <- choose_a_kappa_from_pseudocount(
    A            = A,
    K_target     = K_expected,
    pseudo_factor = pseudo_factor,
    E_eta        = E_eta
  )
  a_kappa <- vol$a_kappa
  b_kappa <- vol$b_kappa
  
  # 2. Degree heterogeneity (as before)
  a_eta <- 2
  b_eta <- 2
  
  # 3. Directional priors (your defaults)
  tau0   <- 1.4   # SST increments
  sigma0 <- 0.84  # WST
  mu0    <- 0
  
  # 4. Gnedin gamma from K_expected
  gamma <- choose_gamma_from_K_expected(
    n,
    K_expected,
    ordering_prior_mode = ordering_prior_mode
  )
  gamma_diag <- gnedin_ordering_discrepancy(n, gamma)
  
  list(
    a_kappa = a_kappa,
    b_kappa = b_kappa,
    a_eta   = a_eta,
    b_eta   = b_eta,
    tau0    = tau0,
    sigma0  = sigma0,
    mu0     = mu0,
    gamma   = gamma,
    meta    = list(
      K_expected    = K_expected,
      ordering_prior_mode = ordering_prior_mode,
      pseudo_factor = pseudo_factor,
      R_typ         = vol$R_typ,
      mean_edges    = vol$mean_edges,
      gamma_discrepancy = gamma_diag
    )
  )
}

#------------------------------------------------------------
# Principled Gibbs/Bayesian calibration
#   - kappa prior weight set via conjugate shrinkage on a
#     "typical" block-pair exposure T_typ
#   - eta prior matched to empirical CV of degree
#   - psi priors from effective-sample-size targets (WST) or
#     dominance target (SST)
#------------------------------------------------------------
get_principled_hypers_gibbs <- function(
    A,
    K_expected        = 5,
    w_kappa_star      = 0.2,   # prior weight in posterior mean for kappa
    sparsity_factor   = 5.0,   # prior mean kappa = mean_edges / sparsity_factor
    N0_psi            = 13,    # WST: effective comparisons per block pair
    dominance_target  = 0.75,  # SST: P(extreme order wins) target
    ordering_prior_mode = c("equivalence_class", "labelled"),
    E_eta             = 1
) {
  ordering_prior_mode <- match.arg(ordering_prior_mode)
  if (w_kappa_star <= 0 || w_kappa_star >= 1) {
    stop("w_kappa_star must be in (0,1).")
  }
  n <- nrow(A)
  if (is.null(n) || n != ncol(A)) {
    stop("A must be a square adjacency / count matrix.")
  }

  mean_edges <- compute_mean_edges(A)
  mu_kappa <- mean_edges / sparsity_factor

  # Typical block size and block-pair exposure
  m_typ <- n / K_expected
  T_typ <- (m_typ^2) * (E_eta^2)

  # w_kappa_star = b / (b + T_typ)  => b = (w/(1-w)) * T_typ
  b_kappa <- (w_kappa_star / (1 - w_kappa_star)) * T_typ
  a_kappa <- mu_kappa * b_kappa

  degrees <- rowSums(A) + colSums(A)
  cv_deg <- sd(degrees) / mean(degrees)
  if (!is.finite(cv_deg) || cv_deg < 0.1) {
    cv_deg <- 0.1
  }
  a_eta <- 1 / (cv_deg^2)
  b_eta <- a_eta

  # WST: N0_psi = 4 / sigma0^2
  sigma0 <- sqrt(4 / N0_psi)
  mu0 <- 0

  # SST: dominance probability at max distance (K-1)
  logit_target <- qlogis(dominance_target)
  D_dist <- max(1, K_expected - 1)
  tau0 <- logit_target / (D_dist * sqrt(2 / pi))

  gamma <- choose_gamma_from_K_expected(
    n,
    K_expected,
    ordering_prior_mode = ordering_prior_mode
  )
  gamma_diag <- gnedin_ordering_discrepancy(n, gamma)

  list(
    a_kappa = a_kappa,
    b_kappa = b_kappa,
    a_eta   = a_eta,
    b_eta   = b_eta,
    tau0    = tau0,
    sigma0  = sigma0,
    mu0     = mu0,
    gamma   = gamma,
    meta    = list(
      K_expected      = K_expected,
      ordering_prior_mode = ordering_prior_mode,
      mean_edges      = mean_edges,
      sparsity_factor = sparsity_factor,
      T_typ           = T_typ,
      w_kappa_star    = w_kappa_star,
      cv_deg          = cv_deg,
      N0_psi          = N0_psi,
      dominance_target = dominance_target,
      gamma_discrepancy = gamma_diag
    )
  )
}

compute_dc_sbm_volume_priors <- function(A, a_lambda = 1) {
  
  # empirical mean intensity
  mean_lambda <- mean(A[A > 0])  # or: mean(A), or mean(A[A != 0])
  # If network is very sparse, you may want mean(A) instead.
  
  if (!is.finite(mean_lambda) || mean_lambda <= 0) {
    mean_lambda <- 1e-6
  }
  
  # scale-matching rate
  b_lambda <- a_lambda / mean_lambda
  
  list(
    a_lambda = a_lambda,
    b_lambda = b_lambda,
    mean_lambda = mean_lambda
  )
}
get_K_dcsbm <- function(dataset) {
  # dataset names must match those used in choose_dataset()
  switch(dataset,
         "citations_data"   = 7L,
         "macaques_data"    = 9L,
         "mountain_goats"   = 4L,
         8L)  # fallback / default if something new appears
}

#------------------------------------------------------------
# Principled kappa calibration  (February 2026)
#
# Derivation summary
# ------------------
# The birth kappa score for node i creating a new block is
#   B_kappa = Σ_ℓ gp_marginal(r_ℓ, t_ℓ, a, b)
#
# In the regime b >> t_ℓ this simplifies to
#   B_kappa ≈ Σ lg(a + r_ℓ) − K lg(a) − R_i ln(b) + O(t/b)
#
# where R_i = Σ r_ℓ is node i's total dyad count.
#
# The join score for an established block is approximately
# independent of (a, b), so the birth advantage is
#   BA_kappa ≈ −R_i · ln(b) + C(a, data)
#
# The intercept C scales as R_i · ln(n/K), giving crossover
#   b* ≈ n / K_expected       (density-independent!)
#
# Define the dimensionless concentration  c = b K / n .
# Then  BA_kappa ≈ −R_i · ln(c), i.e.
#   • c < 1 : births easy         (b < n/K)
#   • c = 1 : crossover           (b = n/K)
#   • c > 1 : births suppressed   (b > n/K)
# The slope is −R_i (node degree), so the penalty naturally
# scales with density — no explicit density correction needed.
#
# Numerical verification on 4 datasets (mountain_goats, citations,
# macaques, high_school) confirms:
#   1.  Monotone for c ≥ 1 across all datasets
#   2.  c = 5 suppresses births in ≥ 94 % of nodes universally
#   3.  Rule "a = λ̄ · b" (mean-matched) is ≈ 15× weaker
#
# Final rule:  a_kappa = 1,  b_kappa = c · n / K_expected
# Default c = 5  (moderate suppression; −1.6 R_i log-units).
#------------------------------------------------------------
get_principled_kappa <- function(
    n,
    K_expected = 5,
    c_kappa    = 5
) {
  if (K_expected < 1) K_expected <- 1
  if (K_expected > n) K_expected <- n
  b_kappa <- c_kappa * n / K_expected
  a_kappa <- 1
  list(a_kappa = a_kappa, b_kappa = b_kappa,
       meta = list(c_kappa = c_kappa, n = n, K_expected = K_expected,
                   crossover_b = n / K_expected))
}

#------------------------------------------------------------
# Spectral K estimator — Rule A+ (March 2026)
#
# Estimates K from the symmetrised adjacency A + t(A) via the
# eigenvalue-ratio criterion:  K_hat = argmax_{k} lambda_k / lambda_{k+1}
#
# Derivation idea
# ---------------
# Under an SBM with K blocks the symmetrised matrix A + t(A) has at most K
# non-trivial eigenvalues.  In practice the ratio lambda_k/lambda_{k+1}
# spikes at k = K_true.  We search over k = 1 .. K_max and return the
# argmax, clipped to [2, floor(n/3)] so that the rule stays conservative.
#
# No external package (irlba, etc.) is required — base::eigen() is used on
# the symmetrised matrix.
#
# Arguments:
#   A       square adjacency / count matrix (directed or undirected).
#   K_max   upper cap on the search range (default = min(30, floor(n/3))).
#   K_min   lower cap on the returned K (default = 2).
#
# Returns: scalar integer K_hat.
#------------------------------------------------------------
estimate_K_spectral <- function(A, K_max = NULL, K_min = 2L) {
  n <- nrow(A)
  if (is.null(n) || n != ncol(A)) stop("A must be square.", call. = FALSE)

  # Symmetrise and remove self-loops
  S <- A + t(A)
  diag(S) <- 0

  K_max <- if (is.null(K_max)) min(30L, floor(n / 3L)) else as.integer(K_max)
  K_max <- max(K_max, K_min + 1L)
  K_min <- max(1L, as.integer(K_min))

  # Compute top K_max+1 eigenvalues of S
  n_eig <- min(K_max + 1L, n)
  ev <- sort(abs(eigen(S, symmetric = TRUE, only.values = TRUE)$values),
             decreasing = TRUE)
  ev <- ev[seq_len(min(n_eig, length(ev)))]

  # Add a small floor so zero eigenvalues don't dominate the ratio
  ev_floor <- max(ev) * 1e-6 + 1e-12
  ev_safe  <- pmax(ev, ev_floor)

  # Eigenvalue-ratio criterion: ratio at position k uses ev[k]/ev[k+1]
  n_ratios <- length(ev_safe) - 1L
  if (n_ratios < 1L) return(as.integer(max(K_min, 1L)))

  ratios <- ev_safe[seq_len(n_ratios)] / ev_safe[seq_len(n_ratios) + 1L]
  k_hat  <- which.max(ratios)  # 1-based index = number of blocks

  # Clip to [K_min, K_max]
  as.integer(min(max(k_hat, K_min), K_max))
}

#------------------------------------------------------------
# Full principled calibration using get_principled_kappa
# When K_expected = NULL the spectral estimator is used (Rule A+).
#------------------------------------------------------------
get_principled_hypers_v2 <- function(
    A,
    K_expected        = NULL,    # NULL => auto-detect via spectral estimator
    c_kappa           = 5,       # kappa concentration (c = bK/n)
    dominance_target  = 0.75,    # SST: P(extreme order wins)
    N0_psi            = 13,      # WST: effective comparisons
    ordering_prior_mode = c("equivalence_class", "labelled")
) {
  ordering_prior_mode <- match.arg(ordering_prior_mode)
  n <- nrow(A)
  if (is.null(n) || n != ncol(A))
    stop("A must be a square adjacency / count matrix.")

  # ---- Auto-detect K_expected via spectral estimator when not provided ----
  if (is.null(K_expected)) {
    K_expected <- estimate_K_spectral(A)
  } else {
    K_expected <- max(2L, min(as.integer(K_expected), n))
  }

  # ---- Kappa: principled rule  a=1, b=c·n/K ----
  kp <- get_principled_kappa(n, K_expected, c_kappa)

  # ---- Eta: match degree CV ----
  degrees <- rowSums(A) + colSums(A)
  cv_deg <- sd(degrees) / mean(degrees)
  if (!is.finite(cv_deg) || cv_deg < 0.1) cv_deg <- 0.1
  a_eta <- 1 / (cv_deg^2)
  b_eta <- a_eta

  # ---- Directional (WST): sigma0 via effective sample size ----
  sigma0 <- sqrt(4 / N0_psi)
  mu0    <- 0

  # ---- Directional (SST): tau0 via dominance target ----
  logit_target <- qlogis(dominance_target)
  D_dist <- max(1, K_expected - 1)
  tau0 <- logit_target / (D_dist * sqrt(2 / pi))

  # ---- Gnedin gamma ----
  gamma <- choose_gamma_from_K_expected(
    n, K_expected, ordering_prior_mode = ordering_prior_mode
  )
  gamma_diag <- gnedin_ordering_discrepancy(n, gamma)

  list(
    a_kappa = kp$a_kappa,
    b_kappa = kp$b_kappa,
    a_eta   = a_eta,
    b_eta   = b_eta,
    tau0    = tau0,
    sigma0  = sigma0,
    mu0     = mu0,
    gamma   = gamma,
    meta    = list(
      K_expected    = K_expected,
      c_kappa       = c_kappa,
      crossover_b   = kp$meta$crossover_b,
      ordering_prior_mode = ordering_prior_mode,
      cv_deg        = cv_deg,
      gamma_discrepancy = gamma_diag,
      K_spectral    = K_expected   # set to spectral estimate if auto-detected
    )
  )
}

#------------------------------------------------------------
# Corollary-based practical calibration
#------------------------------------------------------------
get_corollary_calibrated_hypers <- function(
    A,
    K_expected = 5,
    ordering_prior_mode = c("equivalence_class", "labelled"),
    a_kappa = 2,
    a_eta = 2,
    mu0 = 1.0,
    gamma_bounds = c(0.3, 0.7)
) {
  ordering_prior_mode <- match.arg(ordering_prior_mode)
  n <- nrow(A)
  if (is.null(n) || n != ncol(A)) {
    stop("A must be a square adjacency / count matrix.")
  }

  mean_edges <- compute_mean_edges(A)

  a_kappa <- min(max(a_kappa, 1), 3)
  b_kappa <- a_kappa / mean_edges

  a_eta <- min(max(a_eta, 1), 3)
  b_eta <- a_eta

  mu0 <- max(mu0, 0.5)
  sigma0 <- max(2 * mu0, 1.0)

  D_dist <- max(1, K_expected - 1)
  tau0 <- mu0 / D_dist

  gamma_raw <- choose_gamma_from_K_expected(
    n,
    K_expected,
    ordering_prior_mode = ordering_prior_mode
  )
  gamma <- min(max(gamma_raw, gamma_bounds[1]), gamma_bounds[2])
  gamma_diag <- gnedin_ordering_discrepancy(n, gamma)

  list(
    a_kappa = a_kappa,
    b_kappa = b_kappa,
    a_eta   = a_eta,
    b_eta   = b_eta,
    tau0    = tau0,
    sigma0  = sigma0,
    mu0     = mu0,
    gamma   = gamma,
    meta    = list(
      K_expected = K_expected,
      ordering_prior_mode = ordering_prior_mode,
      mean_edges = mean_edges,
      gamma_raw = gamma_raw,
      gamma_bounds = gamma_bounds,
      gamma_discrepancy = gamma_diag
    )
  )
}


# ------------------------------------------------------------
# Data summaries that actually matter for OSBM calibration
# ------------------------------------------------------------
osbm_summaries <- function(A) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  n <- nrow(A)

  # dyad totals
  N <- A + t(A)
  diag(N) <- 0
  upper_idx <- which(upper.tri(N), arr.ind = TRUE)
  Nij <- N[upper.tri(N)]
  mean_N <- mean(Nij)

  Nij_pos <- Nij[Nij > 0]
  N_typ <- if (length(Nij_pos)) stats::median(Nij_pos) else 1

  # reversals / reciprocity among informative dyads
  i <- upper_idx[,1]; j <- upper_idx[,2]
  Aij <- A[cbind(i,j)]
  Aji <- A[cbind(j,i)]
  informative <- (Aij + Aji) > 0

  reversal_rate <- if (any(informative)) mean((Aij > 0) & (Aji > 0) & informative) else 0

  # dominance = max(Aij,Aji)/Nij on informative dyads
  dom <- if (any(informative)) (pmax(Aij, Aji)[informative] / (Aij + Aji)[informative]) else NA_real_
  dom_med <- if (length(dom)) stats::median(dom) else NA_real_

  # node “degree” totals for eta heterogeneity heuristic
  deg <- rowSums(A) + colSums(A)
  cv_deg <- stats::sd(deg) / mean(deg)
  if (!is.finite(cv_deg) || cv_deg < 0.1) cv_deg <- 0.1

  list(
    n = n,
    mean_N = if (is.finite(mean_N) && mean_N > 0) mean_N else 1e-6,
    N_typ = if (is.finite(N_typ) && N_typ > 0) N_typ else 1,
    reversal_rate = reversal_rate,
    dominance_median = dom_med,
    E_informative = sum(informative),
    cv_deg = cv_deg
  )
}

# ------------------------------------------------------------
# Solve for rho_* such that P(reversal | N_typ, rho_*) ~= target_rev
# reversal prob under Bin(N, rho): 1 - rho^N - (1-rho)^N
# ------------------------------------------------------------
rho_from_reversal_target <- function(N_typ, target_rev) {
  N_typ <- max(1, as.integer(round(N_typ)))
  target_rev <- min(max(target_rev, 0), 1 - 1e-12)

  f <- function(rho) 1 - rho^N_typ - (1 - rho)^N_typ - target_rev

  # If target_rev is too large, rho=0.5 is the max-reversal point
  max_rev_at_half <- 1 - 2 * (0.5^N_typ)
  if (target_rev >= max_rev_at_half) return(0.5)

  # otherwise root on [0.5, 1)
  stats::uniroot(f, interval = c(0.5, 1 - 1e-10))$root
}

# ------------------------------------------------------------
# Main calibrator: returns hypers + meta diagnostics
# ------------------------------------------------------------
calibrate_osbm_hypers <- function(
  A,
  regime = c("SST", "WST"),
  K_expected = NULL,
  ordering_prior_mode = c("equivalence_class", "labelled"),
  # volume controls
  a_kappa = 2,                  # keep modest: CV(kappa)=1/sqrt(a_kappa)
  birth_margin = 3,             # extra log-penalty beyond ~2 log n
  # direction controls
  target_rev_floor = 0.02,      # don't overreact to tiny samples
  target_rev_ceiling = 0.25,    # don't force psi to be small if reversals common
  psi_cap = 8,                  # logit^-1(8)=0.9997; beyond this it's usually model mismatch
  # eta controls
  a_eta_cap = c(0.5, 10),
  # info budget for default K
  c0 = 10
) {
  regime <- match.arg(regime)
  ordering_prior_mode <- match.arg(ordering_prior_mode)

  s <- osbm_summaries(A)
  n <- s$n

  # --- Choose K_expected if not provided.
  # Priority: (1) user-supplied, (2) spectral estimator, (3) info-budget fallback.
  K_spectral_used <- FALSE
  if (is.null(K_expected)) {
    K_spectral <- estimate_K_spectral(A)
    K_expected <- K_spectral
    K_spectral_used <- TRUE
  } else {
    K_spectral <- NA_integer_
    K_expected <- max(2, min(n, as.integer(K_expected)))
  }

  # --- Gnedin gamma to encode K_expected (consistent with your utilities)
  gamma <- choose_gamma_from_K_expected(
    n = n,
    K_expected = K_expected,
    ordering_prior_mode = ordering_prior_mode
  )

  # --- Volume mean m_kappa = a/b
  # Base match to empirical mean dyad total
  m_data <- s$mean_N

  # Also ensure singleton-birth penalty scale beats ~2 log n (neutral birth)
  # Using the crude approximation: penalty ~ m_kappa * (n-1)
  m_min <- (2 * log(n) + birth_margin) / max(1, (n - 1))
  m_kappa <- max(m_data, m_min)

  b_kappa <- a_kappa / m_kappa  # rate parameterisation

  # --- Eta: match degree CV roughly, but keep bounded
  a_eta <- 1 / (s$cv_deg^2)
  a_eta <- min(max(a_eta, a_eta_cap[1]), a_eta_cap[2])
  b_eta <- a_eta  # mean 1 under Gamma(shape=a_eta, rate=b_eta)

  # --- Directional calibration from reversal rate at typical N
  target_rev <- min(max(s$reversal_rate, target_rev_floor), target_rev_ceiling)
  rho_star <- rho_from_reversal_target(s$N_typ, target_rev)
  psi_star <- qlogis(rho_star)
  psi_star_capped <- min(psi_star, psi_cap)

  if (regime == "WST") {
    # Centre at psi_star, but avoid point-mass priors (sigma too small causes pathologies)
    mu0 <- max(0, psi_star_capped)
    sigma0 <- max(1.5, mu0 / 2)   # weakly informative; adjust if you like
    tau0 <- NA_real_
  } else {
    # SST: make increments comparable; choose tau0 so E[psi_{D}] ~ psi_star
    D <- max(1, K_expected - 1)
    # If increments ~ HalfNormal(scale=tau0): E[delta] = tau0*sqrt(2/pi)
    tau0 <- psi_star_capped / (D * sqrt(2 / pi))
    tau0 <- max(tau0, 0.2)        # prevent near-constant psi profile
    mu0 <- 0
    sigma0 <- tau0                # keep delta_0 on same scale as other increments
  }

  list(
    a_kappa = a_kappa,
    b_kappa = b_kappa,
    a_eta   = a_eta,
    b_eta   = b_eta,
    mu0     = mu0,
    sigma0  = sigma0,
    tau0    = tau0,
    gamma   = gamma,
    meta    = list(
      regime = regime,
      ordering_prior_mode = ordering_prior_mode,
      K_expected = K_expected,
      K_spectral = K_spectral,
      K_spectral_used = K_spectral_used,
      summaries = s,
      m_data = m_data,
      m_min_birth = m_min,
      m_kappa_used = m_kappa,
      target_rev = target_rev,
      rho_star = rho_star,
      psi_star = psi_star,
      psi_star_capped = psi_star_capped,
      warning = if (psi_star > psi_cap)
        "psi_star exceeded psi_cap; likely Binomial-direction mismatch (winner-takes-all dyads?)"
      else NULL
    )
  )
}