############################################################
## 0. Libraries & helpers
############################################################

# Reduced OCRP-DGP simulation:
# - default data partitions are sampled from OCRP(vartheta = 3), with a
#   minimum block-size rejection step to avoid unidentified singleton ranks
# - the realized K* is recorded for each replicate
# - fitted ordered OSBM uses the matching OCRP prior family
# - hypers can be fixed or calibrated; fixed hypers are the default
# - first-repetition snapshotting is optional and off by default

get_repo_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grepl("^--file=", cmd)]
  if (!length(file_arg)) {
    return(getwd())
  }
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]),
                               winslash = "/", mustWork = TRUE)
  normalizePath(file.path(dirname(script_path), ".."),
                winslash = "/", mustWork = TRUE)
}

setwd(get_repo_root())

library(stats)
library(salso)
library(BayesLogit)
library(truncnorm)
library(fossil)
library(loo)
library(mcclust)
library(dplyr)
library(tidyr)
library(ggplot2)
library(kableExtra)
library(Matrix)
library(parallel)
library(Rcpp)

# --- Load your existing helpers (unchanged) ---

source("./helper_folder/sim_study_helper.R")
source("./helper_folder/helper.R")
source("./helper_folder/transitivity_check_helper.R")
source("./core/DCSBM_varK.R")        # fit_dcsbm, dcsbm_relabel, etc.
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./helper_folder/Hyper_setup.R")
source("./core/my_best_try_so_far.R")

# sign and logistic helpers (just in case)
sgn      <- function(x) ifelse(x > 0, 1L, ifelse(x < 0, -1L, 0L))
logistic <- function(x) 1 / (1 + exp(-x))

############################################################
## 1. Generative priors: kappa, eta, psi_WST, psi_SST (mean/sd)
############################################################

# Symmetric kappa: Gamma with given mean and variance
sample_kappa_prior <- function(K, mean, var) {
  if (mean <= 0) stop("'mean' must be > 0 for Gamma.")
  if (var  <= 0) stop("'var'  must be > 0 for Gamma.")
  
  shape <- mean^2 / var    # alpha
  rate  <- mean / var      # beta
  
  KAP <- matrix(0, K, K)
  for (k in 1:K) {
    for (l in k:K) {
      KAP[k, l] <- rgamma(1L, shape = shape, rate = rate)
    }
  }
  KAP[lower.tri(KAP)] <- t(KAP)[lower.tri(KAP)]
  KAP
}

# Symmetric kappa with a more intuitive mean/sd parameterization.
sample_kappa_prior_sd <- function(K, mean, sd) {
  if (mean <= 0) stop("'mean' must be > 0 for Gamma.")
  if (sd   <= 0) stop("'sd'   must be > 0 for Gamma.")

  shape <- mean^2 / sd^2
  rate  <- mean / sd^2

  KAP <- matrix(0, K, K)
  for (k in 1:K) {
    for (l in k:K) {
      KAP[k, l] <- rgamma(1L, shape = shape, rate = rate)
    }
  }
  KAP[lower.tri(KAP)] <- t(KAP)[lower.tri(KAP)]
  KAP
}

sample_eta_prior <- function(n, mean = 1, sd = 0.1) {
  if (mean <= 0) stop("'mean' must be > 0 for Gamma.")
  if (sd   <= 0) stop("'sd'   must be > 0 for Gamma.")

  shape <- mean^2 / sd^2
  rate  <- mean / sd^2
  rgamma(n, shape = shape, rate = rate)
}

# DC-SBM truth: block-directed intensities with random pairwise asymmetry.
# The sign flips across block pairs so the induced directionality is generally
# non-transitive and therefore not ordered in the WST/SST sense.
sample_lambda_dcsbm_truth <- function(K, mean_kappa, sd_kappa,
                                      asym_mean, asym_sd) {
  if (K < 1L) stop("K must be positive.", call. = FALSE)
  if (mean_kappa <= 0 || sd_kappa <= 0) {
    stop("mean_kappa and sd_kappa must be positive.", call. = FALSE)
  }
  if (!is.finite(asym_mean) || asym_mean < 0 || !is.finite(asym_sd) || asym_sd <= 0) {
    stop("asym_mean must be >= 0 and asym_sd must be > 0.", call. = FALSE)
  }

  diag_mean <- mean_kappa / 4
  diag_sd   <- max(sd_kappa / 5, 0.04)
  off_mean  <- mean_kappa / 10
  off_sd    <- max(sd_kappa / 8, 0.02)

  lambda <- matrix(0, K, K)
  for (k in seq_len(K)) {
    lambda[k, k] <- rgamma(
      1L,
      shape = diag_mean^2 / diag_sd^2,
      rate  = diag_mean / diag_sd^2
    )
  }
  if (K >= 2L) {
    for (k in 1:(K - 1L)) {
      for (l in (k + 1L):K) {
        base <- rgamma(
          1L,
          shape = off_mean^2 / off_sd^2,
          rate  = off_mean / off_sd^2
        )
        log_ratio <- rnorm(1L, mean = asym_mean, sd = asym_sd)
        sign_flip <- sample(c(-1, 1), size = 1L)
        mult <- exp(sign_flip * log_ratio / 2)
        lambda[k, l] <- base * mult
        lambda[l, k] <- base / mult
      }
    }
  }
  lambda
}

# WST psi: block-pair specific, skew-symmetric
# psi_{kℓ} (k<ℓ) ~ N^+(mu_psi, sigma_psi^2)
sample_psi_wst_prior <- function(K, mu_psi, var_psi) {
  if (var_psi <= 0) stop("'var_psi' must be > 0")
  sigma_psi <- sqrt(var_psi)
  
  PSI <- matrix(0, K, K)
  PSI[upper.tri(PSI)] <- truncnorm::rtruncnorm(
    n    = K * (K - 1) / 2,
    a    = 0,
    mean = mu_psi,
    sd   = sigma_psi
  )
  PSI[lower.tri(PSI)] <- -t(PSI)[lower.tri(PSI)]
  PSI
}

# SST psi: distance-based, psi_d = sum_{i<=d} delta_i.
# In TOTAL mode, mu_psi is the approximate total span E[psi_{K-1}].
# In INCREMENT mode, mu_psi is the adjacent-rank increment mean.
sample_psi_sst_prior <- function(K, mu_psi, var_psi,
                                 mode = c("TOTAL", "INCREMENT")) {
  if (K <= 1) return(numeric(0))
  if (var_psi <= 0) stop("'var_psi' must be > 0")

  mode <- match.arg(toupper(mode), c("TOTAL", "INCREMENT"))
  if (mode == "TOTAL") {
    mu_delta  <- mu_psi / (K - 1)
    tau_delta <- sqrt(var_psi) / (K - 1)
  } else {
    mu_delta  <- mu_psi
    tau_delta <- sqrt(var_psi)
  }
  
  deltas <- truncnorm::rtruncnorm(
    n    = K - 1,
    a    = 0,
    mean = mu_delta,
    sd   = tau_delta
  )
  cumsum(deltas)
}

############################################################
## 2. Hyperparameters for OSBM (principled Gibbs calibration)
############################################################

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

demo_osbm_partition_prior <- toupper(trimws(Sys.getenv("DEMOOCRPVAR_OSBM_PARTITION_PRIOR", unset = "OCRP")))
if (!demo_osbm_partition_prior %in% c("OCRP", "ROCRP")) {
  stop("DEMOOCRPVAR_OSBM_PARTITION_PRIOR must be OCRP or ROCRP.", call. = FALSE)
}
demo_alpha_crp <- as.numeric(Sys.getenv("DEMOOCRPVAR_ALPHA_CRP", unset = ""))
demo_theta_override <- as.numeric(Sys.getenv("DEMOOCRPVAR_THETA_OCRP", unset = "3"))
demo_dgp_partition <- toupper(trimws(Sys.getenv("DEMOOCRPVAR_DGP_PARTITION", unset = "OCRP_MIN")))
if (demo_dgp_partition == "OCRP_CONDITIONED") demo_dgp_partition <- "OCRP_FIXED_K"
if (!demo_dgp_partition %in% c("OCRP", "OCRP_MIN", "OCRP_FIXED_K", "OCRP_FIXED_K_MIN", "BALANCED")) {
  stop("DEMOOCRPVAR_DGP_PARTITION must be OCRP, OCRP_MIN, OCRP_FIXED_K, OCRP_FIXED_K_MIN, or BALANCED.", call. = FALSE)
}
if (!is.finite(demo_theta_override) && demo_dgp_partition %in% c("OCRP", "OCRP_MIN")) {
  demo_theta_override <- 3
}
demo_min_block_size <- suppressWarnings(as.integer(Sys.getenv("DEMOOCRPVAR_MIN_BLOCK_SIZE", unset = "3")))
if (is.na(demo_min_block_size) || demo_min_block_size < 1L) demo_min_block_size <- 1L
demo_sst_psi_mode <- toupper(trimws(Sys.getenv("DEMOOCRPVAR_SST_PSI_MODE", unset = "INCREMENT")))
if (!demo_sst_psi_mode %in% c("TOTAL", "INCREMENT")) {
  stop("DEMOOCRPVAR_SST_PSI_MODE must be TOTAL or INCREMENT.", call. = FALSE)
}
demo_eta_mode <- toupper(trimws(Sys.getenv("DEMOOCRPVAR_ETA_MODE", unset = "GAMMA")))
if (!demo_eta_mode %in% c("GAMMA", "CONSTANT")) {
  stop("DEMOOCRPVAR_ETA_MODE must be GAMMA or CONSTANT.", call. = FALSE)
}
demo_eta_mean <- suppressWarnings(as.numeric(Sys.getenv("DEMOOCRPVAR_ETA_MEAN", unset = "1")))
demo_eta_sd <- suppressWarnings(as.numeric(Sys.getenv("DEMOOCRPVAR_ETA_SD", unset = "0.10")))
if (!is.finite(demo_eta_mean) || demo_eta_mean <= 0) demo_eta_mean <- 1
if (!is.finite(demo_eta_sd) || demo_eta_sd <= 0) demo_eta_sd <- 0.10
demo_save_loo_pointwise <- tolower(Sys.getenv("DEMOOCRPVAR_SAVE_LOO_POINTWISE", unset = "false")) %in% c("1", "true", "yes")
demo_order_dominance_target <- suppressWarnings(as.numeric(Sys.getenv("DEMOOCRPVAR_DOMINANCE_TARGET", unset = "")))
demo_order_N0_psi <- suppressWarnings(as.numeric(Sys.getenv("DEMOOCRPVAR_N0_PSI", unset = "")))
demo_tau0_override <- suppressWarnings(as.numeric(Sys.getenv("DEMOOCRPVAR_TAU0", unset = "")))
demo_sigma0_override <- suppressWarnings(as.numeric(Sys.getenv("DEMOOCRPVAR_SIGMA0", unset = "")))
demo_tau0_floor <- suppressWarnings(as.numeric(Sys.getenv("DEMOOCRPVAR_TAU0_FLOOR", unset = "0.2")))
if (!is.finite(demo_order_dominance_target) ||
    demo_order_dominance_target <= 0 ||
    demo_order_dominance_target >= 1) {
  demo_order_dominance_target <- 0.75
}
if (!is.finite(demo_order_N0_psi) || demo_order_N0_psi <= 0) {
  demo_order_N0_psi <- 13
}
if (!is.finite(demo_tau0_floor) || demo_tau0_floor < 0) {
  demo_tau0_floor <- 0.2
}
parse_numeric_vec_env <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(trimws(raw))) return(default)
  vals <- suppressWarnings(as.numeric(strsplit(raw, ",")[[1L]]))
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return(default)
  vals
}

parse_numeric_env <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(trimws(raw))) return(default)
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val)) return(default)
  val
}

parse_bool_env <- function(name, default = FALSE) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(trimws(raw))) return(default)
  tolower(trimws(raw)) %in% c("1", "true", "yes")
}

parse_char_vec_env <- function(name, default, allowed = NULL) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(trimws(raw))) return(default)
  vals <- toupper(trimws(strsplit(raw, ",")[[1L]]))
  vals <- vals[nzchar(vals)]
  if (!is.null(allowed)) vals <- vals[vals %in% allowed]
  if (!length(vals)) return(default)
  vals
}

parse_hierarchy_specs_env <- function(default) {
  raw <- Sys.getenv("DEMOOCRPVAR_HIERARCHY_SPECS", unset = "")
  if (!nzchar(trimws(raw))) return(default)
  pieces <- trimws(strsplit(raw, ",")[[1L]])
  rows <- lapply(pieces, function(piece) {
    parts <- trimws(strsplit(piece, ":")[[1L]])
    if (length(parts) < 2L) return(NULL)
    psi_mean <- suppressWarnings(as.numeric(parts[2L]))
    psi_sd <- if (length(parts) >= 3L) suppressWarnings(as.numeric(parts[3L])) else 0.30
    if (!is.finite(psi_mean) || !is.finite(psi_sd) || psi_sd <= 0) return(NULL)
    data.frame(hierch = parts[1L], psi_mean = psi_mean, psi_sd = psi_sd,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || !nrow(out)) default else out
}

demo_use_fixed_hypers <- parse_bool_env("DEMOOCRPVAR_FIXED_HYPERS", default = TRUE)
demo_first_rep_snapshot <- parse_bool_env("DEMOOCRPVAR_FIRST_REP_SNAPSHOT", default = FALSE)
demo_a_kappa_fixed <- parse_numeric_env("DEMOOCRPVAR_A_KAPPA", 1)
demo_b_kappa_fixed <- parse_numeric_env("DEMOOCRPVAR_B_KAPPA", 1)
demo_a_eta_fixed <- parse_numeric_env("DEMOOCRPVAR_A_ETA", 1)
demo_b_eta_fixed <- parse_numeric_env("DEMOOCRPVAR_B_ETA", 1)
demo_mu0_wst <- parse_numeric_env("DEMOOCRPVAR_MU0_WST", 0)
demo_sigma0_wst <- parse_numeric_env(
  "DEMOOCRPVAR_SIGMA0_WST",
  if (is.finite(demo_sigma0_override) && demo_sigma0_override > 0) demo_sigma0_override else 1.5
)
demo_mu0_sst <- parse_numeric_env("DEMOOCRPVAR_MU0_SST", 0)
demo_sigma0_sst <- parse_numeric_env(
  "DEMOOCRPVAR_SIGMA0_SST",
  if (is.finite(demo_sigma0_override) && demo_sigma0_override > 0) demo_sigma0_override else 1.5
)
demo_tau0_sst <- parse_numeric_env(
  "DEMOOCRPVAR_TAU0_SST",
  if (is.finite(demo_tau0_override) && demo_tau0_override > 0) demo_tau0_override else 0.5
)
demo_gamma_gn_fixed <- parse_numeric_env("DEMOOCRPVAR_GAMMA_GN", 0.8)

parse_scenario_specs_env <- function() {
  raw <- Sys.getenv("DEMOOCRPVAR_SCENARIO_SPECS", unset = "")
  if (!nzchar(trimws(raw))) return(NULL)
  pieces <- trimws(strsplit(raw, ",")[[1L]])
  rows <- lapply(seq_along(pieces), function(i) {
    parts <- trimws(strsplit(pieces[[i]], ":")[[1L]])
    if (!(length(parts) %in% c(6L, 7L, 8L))) {
      stop(
        paste(
          "DEMOOCRPVAR_SCENARIO_SPECS entries must be one of:",
          "label:gen_model:kappa_mean:hierch:psi_mean:psi_sd,",
          "label:gen_model:K:kappa_mean:hierch:psi_mean:psi_sd,",
          "label:gen_model:K:kappa_mean:kappa_sd:hierch:psi_mean:psi_sd."
        ),
        call. = FALSE
      )
    }
    gen_model <- toupper(parts[2L])
    if (!gen_model %in% c("SST", "WST", "DC-SBM")) {
      stop("Scenario gen_model must be SST, WST, or DC-SBM.", call. = FALSE)
    }
    if (length(parts) == 6L) {
      K_true <- NA_integer_
      kappa_mean <- suppressWarnings(as.numeric(parts[3L]))
      kappa_sd <- NA_real_
      hierch <- parts[4L]
      psi_mean <- suppressWarnings(as.numeric(parts[5L]))
      psi_sd <- suppressWarnings(as.numeric(parts[6L]))
    } else if (length(parts) == 7L && !is.na(suppressWarnings(as.integer(parts[3L])))) {
      K_true <- suppressWarnings(as.integer(parts[3L]))
      kappa_mean <- suppressWarnings(as.numeric(parts[4L]))
      kappa_sd <- NA_real_
      hierch <- parts[5L]
      psi_mean <- suppressWarnings(as.numeric(parts[6L]))
      psi_sd <- suppressWarnings(as.numeric(parts[7L]))
    } else {
      K_true <- suppressWarnings(as.integer(parts[3L]))
      kappa_mean <- suppressWarnings(as.numeric(parts[4L]))
      kappa_sd <- suppressWarnings(as.numeric(parts[5L]))
      hierch <- parts[6L]
      psi_mean <- suppressWarnings(as.numeric(parts[7L]))
      psi_sd <- suppressWarnings(as.numeric(parts[8L]))
    }
    if ((!is.na(K_true) && (!is.finite(K_true) || K_true < 1L)) ||
        !is.finite(kappa_mean) || kappa_mean <= 0 ||
        (!is.na(kappa_sd) && (!is.finite(kappa_sd) || kappa_sd <= 0)) ||
        !is.finite(psi_mean) || !is.finite(psi_sd) || psi_sd <= 0) {
      stop("Invalid numeric value in DEMOOCRPVAR_SCENARIO_SPECS.", call. = FALSE)
    }
    data.frame(
      scenario_label = parts[1L],
      gen_model = gen_model,
      K_true = K_true,
      kappa_mean = kappa_mean,
      kappa_sd = kappa_sd,
      hierch = hierch,
      psi_mean = psi_mean,
      psi_sd = psi_sd,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

sample_balanced_partition <- function(n, K_target, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sizes <- rep(floor(n / K_target), K_target)
  sizes[seq_len(n %% K_target)] <- sizes[seq_len(n %% K_target)] + 1L
  z <- rep(seq_len(K_target), times = sizes)
  z <- z[sample.int(n)]
  list(K = K_target, blocks = as.integer(tabulate(z, nbins = K_target)),
       z = as.integer(z), attempts = NA_integer_)
}

crp_expected_K <- function(n, theta) {
  sum(theta / (theta + seq.int(0L, n - 1L)))
}

choose_theta_for_expected_K <- function(n, K_expected) {
  if (!is.finite(K_expected) || K_expected < 1 || K_expected > n) {
    stop("K_expected must be in [1, n].", call. = FALSE)
  }
  if (K_expected <= 1 + 1e-8) return(1e-8)
  if (K_expected >= n - 1e-8) return(1e8)
  f <- function(theta) crp_expected_K(n, theta) - K_expected
  uniroot(f, lower = 1e-8, upper = 1e8, tol = 1e-10)$root
}

sample_ocrp_partition <- function(n, theta_ocrp, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (n < 1L) stop("n must be positive.", call. = FALSE)
  if (theta_ocrp <= 0) stop("theta_ocrp must be positive.", call. = FALSE)
  
  blocks <- 1L
  z <- integer(n)
  z[1L] <- 1L
  
  for (i in 2:n) {
    H <- length(blocks)
    pw <- ocrp_log_weights_packed(blocks, theta_ocrp = theta_ocrp)
    probs <- c(exp(pw$exist), exp(pw$new))
    choice <- sample.int(length(probs), 1L, prob = probs)
    
    if (choice <= H) {
      blocks[choice] <- blocks[choice] + 1L
      z[i] <- choice
    } else {
      r <- choice - H
      prev_idx <- seq_len(i - 1L)
      shift_idx <- prev_idx[z[prev_idx] >= r]
      z[shift_idx] <- z[shift_idx] + 1L
      if (r == 1L) {
        blocks <- c(1L, blocks)
      } else if (r == H + 1L) {
        blocks <- c(blocks, 1L)
      } else {
        blocks <- c(blocks[1:(r - 1L)], 1L, blocks[r:H])
      }
      z[i] <- r
    }
  }
  
  perm <- sample.int(n)
  z <- z[perm]
  list(
    K = length(blocks),
    blocks = as.integer(tabulate(z, nbins = length(blocks))),
    z = as.integer(z),
    attempts = 1L
  )
}

sample_ocrp_partition_min_blocks <- function(n, theta_ocrp, seed,
                                             min_block_size = 1L,
                                             max_attempts = 20000L) {
  for (attempt in seq_len(max_attempts)) {
    attempt_seed <- as.integer(((as.double(seed) + (attempt - 1) * 1000003) %% 2147483646) + 1)
    out <- sample_ocrp_partition(n, theta_ocrp, seed = attempt_seed)
    if (min(out$blocks) >= min_block_size) {
      out$attempts <- attempt
      return(out)
    }
  }
  stop(sprintf(
    "Failed to sample OCRP partition with min block size %d after %d attempts (n=%d, theta=%.4g).",
    min_block_size, max_attempts, n, theta_ocrp
  ), call. = FALSE)
}

sample_ocrp_partition_conditioned_K <- function(n, K_target, theta_ocrp, seed,
                                                min_block_size = 1L,
                                                max_attempts = 20000L) {
  if (K_target * min_block_size > n) {
    stop(sprintf(
      "Cannot sample K=%d blocks with min_block_size=%d when n=%d.",
      K_target, min_block_size, n
    ), call. = FALSE)
  }
  for (attempt in seq_len(max_attempts)) {
    attempt_seed <- as.integer(((as.double(seed) + (attempt - 1) * 1000003) %% 2147483646) + 1)
    out <- sample_ocrp_partition(n, theta_ocrp, seed = attempt_seed)
    if (out$K == K_target && min(out$blocks) >= min_block_size) {
      out$attempts <- attempt
      return(out)
    }
  }
  stop(sprintf(
    "Failed to sample OCRP partition with K=%d and min block size %d after %d attempts (n=%d, theta=%.4g).",
    K_target, min_block_size, max_attempts, n, theta_ocrp
  ), call. = FALSE)
}

############################################################
## 3. One replicate: generate + fit all three models
############################################################

run_replicate_crossfit <- function(
    gen_model,
    K_true,
    n,
    kappa_mean,
    kappa_var,
    psi_mean,
    psi_var,
    hierch_label,
    theta_ocrp,
    seed,
    rep_id = NA_integer_,
    scenario_id = NA_character_,
    density_label = NA_character_,
    n_iter = 5000,
    burn   = 1000,
    thin   = 2
) {
  set.seed(seed)
  K_target <- K_true
  
  # --- 1. Generate data from OSBM prior ---
  if (demo_dgp_partition == "OCRP") {
    z_draw <- sample_ocrp_partition(
      n = n,
      theta_ocrp = theta_ocrp,
      seed = seed + 100000L
    )
  } else if (demo_dgp_partition == "OCRP_MIN") {
    z_draw <- sample_ocrp_partition_min_blocks(
      n = n,
      theta_ocrp = theta_ocrp,
      seed = seed + 100000L,
      min_block_size = demo_min_block_size
    )
  } else if (demo_dgp_partition == "BALANCED") {
    if (is.na(K_target) || !is.finite(K_target) || K_target < 1L) {
      stop("BALANCED DGP needs a positive K target.", call. = FALSE)
    }
    z_draw <- sample_balanced_partition(
      n = n,
      K_target = K_target,
      seed = seed + 100000L
    )
  } else {
    if (is.na(K_target) || !is.finite(K_target) || K_target < 1L) {
      stop("Fixed-K OCRP DGP needs a positive K target.", call. = FALSE)
    }
    z_draw <- sample_ocrp_partition_conditioned_K(
      n = n,
      K_target = K_target,
      theta_ocrp = theta_ocrp,
      seed = seed + 100000L,
      min_block_size = if (demo_dgp_partition == "OCRP_FIXED_K_MIN") demo_min_block_size else 1L
    )
  }
  K_true <- z_draw$K
  z_true <- z_draw$z
  block_sizes_true <- paste(tabulate(z_true, nbins = K_true), collapse = ";")
  kappa_true <- sample_kappa_prior_sd(K_true, mean = kappa_mean, sd = sqrt(kappa_var))
  eta_true   <- if (demo_eta_mode == "CONSTANT") rep(demo_eta_mean, n) else sample_eta_prior(n, mean = demo_eta_mean, sd = demo_eta_sd)

  if (gen_model == "DC-SBM") {
    theta_true <- normalize_block_theta(z_true, eta_true)
    lambda_true <- sample_lambda_dcsbm_truth(
      K = K_true,
      mean_kappa = kappa_mean,
      sd_kappa = sqrt(kappa_var),
      asym_mean = psi_mean,
      asym_sd = sqrt(psi_var)
    )
    sim_dc <- simulate_dcsbm(
      n = n,
      K = K_true,
      z = z_true,
      theta = theta_true,
      lambda = lambda_true,
      seed = seed + 200000L
    )
    A <- sim_dc$A
  } else if (gen_model == "WST") {
    psi_true <- sample_psi_wst_prior(K_true, mu_psi = psi_mean, var_psi = psi_var)
    A        <- simulate_osbm(n, K_true, z_true, eta_true, kappa_true, psi_true, "WST")
  } else {
    psi_true <- sample_psi_sst_prior(
      K_true,
      mu_psi = psi_mean,
      var_psi = psi_var,
      mode = demo_sst_psi_mode
    )
    A        <- simulate_osbm(n, K_true, z_true, eta_true, kappa_true, psi_true, "SST")
  }
  A <- as.matrix(A)
  
  # --- 2. Model hypers ---
  K_init_fit <- nrow(A)
  hyper_source <- if (demo_use_fixed_hypers) "fixed" else "calibrated"

  if (demo_use_fixed_hypers) {
    hypers_wst <- list(
      a_kappa = demo_a_kappa_fixed,
      b_kappa = demo_b_kappa_fixed,
      a_eta = demo_a_eta_fixed,
      b_eta = demo_b_eta_fixed,
      mu0 = demo_mu0_wst,
      sigma0 = demo_sigma0_wst,
      tau0 = demo_tau0_sst,
      gamma = demo_gamma_gn_fixed
    )
    hypers_sst <- list(
      a_kappa = demo_a_kappa_fixed,
      b_kappa = demo_b_kappa_fixed,
      a_eta = demo_a_eta_fixed,
      b_eta = demo_b_eta_fixed,
      mu0 = demo_mu0_sst,
      sigma0 = demo_sigma0_sst,
      tau0 = demo_tau0_sst,
      gamma = demo_gamma_gn_fixed
    )
    gamma_gn_tuned <- demo_gamma_gn_fixed
  } else {
    K_expected_fit <- K_true
    hypers <- get_principled_hypers_v2(
      A = A,
      K_expected = K_expected_fit,
      c_kappa = 3,
      dominance_target = demo_order_dominance_target,
      N0_psi = demo_order_N0_psi,
      ordering_prior_mode = "equivalence_class"
    )
    hypers$tau0 <- max(hypers$tau0, demo_tau0_floor)
    if (is.finite(demo_tau0_override) && demo_tau0_override > 0) {
      hypers$tau0 <- demo_tau0_override
    }
    if (is.finite(demo_sigma0_override) && demo_sigma0_override > 0) {
      hypers$sigma0 <- demo_sigma0_override
    }
    gamma_gn_tuned <- hypers$gamma

    hypers_wst <- hypers
    hypers_wst$b_kappa <- demo_b_kappa_fixed
    hypers_wst$b_eta <- demo_b_eta_fixed
    hypers_wst$a_kappa <- demo_a_kappa_fixed
    hypers_wst$a_eta <- demo_a_eta_fixed

    hypers_sst <- hypers
    hypers_sst$b_kappa <- demo_b_kappa_fixed
    hypers_sst$b_eta <- demo_b_eta_fixed
    hypers_sst$a_kappa <- demo_a_kappa_fixed
    hypers_sst$a_eta <- demo_a_eta_fixed
  }
  
  # --- 3. Fit models: WST, SST, DC-SBM ---
  models_to_fit <- c("WST", "SST", "DC-SBM")
  results_list  <- list()
  loo_cache     <- list()
  DI            <- build_dyad_index(nrow(A))
  
  for (fit_model in models_to_fit) {
    loop_seed <- seed + which(models_to_fit == fit_model)
    set.seed(loop_seed)
    start_time <- proc.time()
    
    mcmc_out_for_diag <- NULL
    regime_for_diag   <- if (fit_model == "SST") "SST" else "WST"
    fit_sigma0_used   <- NA_real_
    fit_tau0_used     <- NA_real_
    fit_mu0_used      <- NA_real_
    fit_a_kappa_used  <- NA_real_
    fit_b_kappa_used  <- NA_real_
    fit_a_eta_used    <- NA_real_
    fit_b_eta_used    <- NA_real_
    
    # --- 3a. DC-SBM: initialize from K = n and learn the occupied K ---
    if (fit_model == "DC-SBM") {
      
      fit <- fit_dcsbm_gibbs_gnedin(
        as.matrix(A),
        K_init     = K_init_fit,
        priors = list(
          a_eta         = demo_a_eta_fixed,
          b_eta         = demo_b_eta_fixed,
          a_lambda      = demo_a_kappa_fixed,
          b_lambda      = demo_b_kappa_fixed,
          alpha_crp     = if (is.finite(demo_alpha_crp)) demo_alpha_crp else theta_ocrp,
          partition_prior = "CRP"
        ),
        iters   = n_iter,
        burn_in = burn,
        thin    = thin,
        verbose = 0,
        seed    = loop_seed
      )
      fit_a_kappa_used <- demo_a_kappa_fixed
      fit_b_kappa_used <- demo_b_kappa_fixed
      fit_a_eta_used <- demo_a_eta_fixed
      fit_b_eta_used <- demo_b_eta_fixed
      
      relab_out <- relabel_chain_by_block_score(fit$z, fit$lambda, A, A + t(A))
      Z_chain_relab    <- relab_out$Z      # matrix S x n or list of length S
      Lambda_chain_rel <- relab_out$Lambda # list length S, each K_s x K_s
      
      
      # Coerce Z_chain_relab to an S x n matrix for SALSO and diagnostics
      if (is.matrix(Z_chain_relab)) {
        z_chain <- Z_chain_relab              # S x n
      } else {
        z_chain <- do.call(rbind, Z_chain_relab)  # list of length S -> S x n
      }
      S <- nrow(z_chain)
      
      
      LL <- loglik_matrix_dcsbm(A, fit, dyad_index = DI)
      
      ## --- Pseudo-psi + ordered diagnostics (treat DC-SBM as WST on blocks) ---
      psi_list <- vector("list", S)
      z_list   <- vector("list", S)
      
      for (s in seq_len(S)) {
        kap   <- Lambda_chain_rel[[s]]          # K_s x K_s
        K_curr <- nrow(kap)
        vol   <- kap + t(kap)
        
        rho <- matrix(0.5, nrow = K_curr, ncol = K_curr)
        mask <- vol > 0
        rho[mask] <- kap[mask] / vol[mask]
        rho_clamped <- pmin(pmax(rho, 1e-6), 1 - 1e-6)
        
        psi_list[[s]] <- qlogis(rho_clamped)   # WST-style: full K_s x K_s logits
        z_list[[s]]   <- z_chain[s, ]
      }
      
      mcmc_out_for_diag <- list(psi = psi_list, z = z_list)
      
      
    } else {
      # --- 3b. WST / SST OSBM (unchanged sampler) ---

      # #check SST
      # res <- run_fast_bruteforce_tests(A = A, N = A+t(A), z = z_true, psi_vec = psi_true, 
      #                                  i_set = sample(seq_along(z_true), 5),)
      # 
      # #check WST
      # run_fast_bruteforce_tests_wst(
      #   A = A, N = A+t(A), z = z_true, psi_mat = psi_true,
      #   mu0 = 0, sigma0 = 1,
      #   tol = 1e-8, verbose = TRUE
      # )
      
      fail_on_osbm_error <- tolower(Sys.getenv("DEMOOCRPVAR_FAIL_ON_OSBM_ERROR", unset = "false")) %in% c("1", "true", "yes")

      fit_hypers <- if (fit_model == "SST") hypers_sst else hypers_wst
      fit_sigma0_used <- fit_hypers$sigma0
      fit_tau0_used <- if (fit_model == "SST") fit_hypers$tau0 else NA_real_
      fit_mu0_used <- fit_hypers$mu0
      fit_a_kappa_used <- fit_hypers$a_kappa
      fit_b_kappa_used <- fit_hypers$b_kappa
      fit_a_eta_used <- fit_hypers$a_eta
      fit_b_eta_used <- fit_hypers$b_eta
      out <- tryCatch({
        modular_osbm_sampler(
          A = A,
          K = K_init_fit,
          n_iter = n_iter,
          burn   = burn,
          thin   = thin,
          verbose = FALSE,
          a_kappa = fit_hypers$a_kappa,
          b_kappa = fit_hypers$b_kappa,
          a_eta   = fit_hypers$a_eta,
          b_eta   = fit_hypers$b_eta,
          psi_constraint = fit_model,     # "WST" or "SST"
          partition_prior = demo_osbm_partition_prior,
          theta_ocrp      = theta_ocrp,
          eta_identifiability = "block_sum",
          use_mixing_moves = TRUE,
          gamma_gn        = gamma_gn_tuned,
          mu0   = fit_hypers$mu0,
          sigma0 = fit_hypers$sigma0,
          tau0   = fit_hypers$tau0,
          seed = loop_seed
        )
      }, error = function(e) {
        message(sprintf(
          "OSBM fit failed | gen=%s fit=%s K_true=%d hierch=%s kappa_mean=%.3f psi_mean=%.3f seed=%d | %s",
          gen_model, fit_model, K_true, hierch_label, kappa_mean, psi_mean, loop_seed, conditionMessage(e)
        ))
        if (isTRUE(fail_on_osbm_error)) {
          stop(e)
        }
        NULL
      })
      
      if (is.null(out)) next
      
      z_chain <- do.call(rbind, lapply(out$z, as.integer))
      LL      <- loglik_matrix_modular(A, out, regime = fit_model, dyad_index = DI)
      mcmc_out_for_diag <- out
    }
    
    time_taken <- (proc.time() - start_time)["elapsed"]
    
    # --- 4. Clustering metrics ---
    z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L)
    ari   <- fossil::adj.rand.index(z_hat,z_true )
    vi    <- mcclust::vi.dist(z_hat, z_true)
    K_hat <- length(unique(z_hat))
    K_count = function(x)length(unique(x))
    K_draws <- apply(z_chain, MARGIN = 1, K_count)
    K_p.mean <- mean(K_draws)
    K_ci <- stats::quantile(K_draws, probs = c(0.025, 0.975), names = FALSE)
    
    
    loo_res <- tryCatch(loo::loo(LL),
                        error = function(e) NULL)
    
    if (is.null(loo_res)) {
      elpd   <- NA_real_
      elpd_se <- NA_real_
      p_loo <- NA_real_
      p_loo_se <- NA_real_
      looic  <- NA_real_
      looic_se <- NA_real_
      pk_bad <- NA_real_
      pk_bad_count <- NA_integer_
      pareto_k_max <- NA_real_
      pareto_k_mean <- NA_real_
      n_loo_obs <- ncol(LL)
    } else {
      elpd   <- loo_res$estimates["elpd_loo", "Estimate"]
      elpd_se <- loo_res$estimates["elpd_loo", "SE"]
      p_loo <- loo_res$estimates["p_loo", "Estimate"]
      p_loo_se <- loo_res$estimates["p_loo", "SE"]
      looic  <- loo_res$estimates["looic",   "Estimate"]
      looic_se <- loo_res$estimates["looic", "SE"]
      pk_bad <- mean(loo_res$diagnostics$pareto_k > 0.7)
      pk_bad_count <- sum(loo_res$diagnostics$pareto_k > 0.7)
      pareto_k_max <- max(loo_res$diagnostics$pareto_k, na.rm = TRUE)
      pareto_k_mean <- mean(loo_res$diagnostics$pareto_k, na.rm = TRUE)
      n_loo_obs <- nrow(loo_res$pointwise)
      loo_cache[[fit_model]] <- list(
        loo = loo_res,
        pointwise = loo_res$pointwise,
        LL = LL
      )
    }
    
    
    diag_res <- summarise_osbm_diagnostics(
      out_relab = mcmc_out_for_diag,
      regime    = regime_for_diag,
      K_max_hint = K_true,
      z_hat     = z_hat,
      n         = n,
      A         = A,
      m_items   = 500
    )
    
    res_row <- data.frame(
      gen_model  = gen_model,
      fit_model  = fit_model,
      K_target    = ifelse(is.finite(K_target), K_target, NA_integer_),
      K_true     = K_true,
      K_hat      = K_hat,
      K_p.mean   = K_p.mean,
      K_ci_lo    = K_ci[1L],
      K_ci_hi    = K_ci[2L],
      hierch     = hierch_label,           # weak / strong (via psi mean)
      kappa_mean = kappa_mean,
      kappa_var  = kappa_var,
      psi_mean   = psi_mean,
      psi_var    = psi_var,
      theta_ocrp = theta_ocrp,
      ocrp_attempts = z_draw$attempts,
      block_sizes_true = block_sizes_true,
      rate_kap   = kappa_mean,            # for backwards-compatible plots
      ari        = ari,
      vi         = vi,
      elpd       = elpd,
      elpd_se    = elpd_se,
      p_loo      = p_loo,
      p_loo_se   = p_loo_se,
      looic      = looic,
      looic_se   = looic_se,
      pk_bad     = pk_bad,
      pk_bad_count = pk_bad_count,
      pareto_k_max = pareto_k_max,
      pareto_k_mean = pareto_k_mean,
      n_loo_obs  = n_loo_obs,
      time       = time_taken,
      seed       = loop_seed,
      hyper_source = hyper_source,
      mu0        = fit_mu0_used,
      sigma0     = fit_sigma0_used,
      tau0       = fit_tau0_used,
      a_kappa_prior = fit_a_kappa_used,
      b_kappa_prior = fit_b_kappa_used,
      a_eta_prior = fit_a_eta_used,
      b_eta_prior = fit_b_eta_used,
      gamma_gn_prior = gamma_gn_tuned,
      dominance_target = demo_order_dominance_target,
      N0_psi     = demo_order_N0_psi,
      stringsAsFactors = FALSE
    )
    print(res_row)
    
    results_list[[length(results_list) + 1L]] <-
      cbind(res_row, as.data.frame(as.list(diag_res)))
  }
  
  result_df <- do.call(rbind, results_list)

  if (!is.null(result_df) && nrow(result_df) && length(loo_cache)) {
    model_names <- result_df$fit_model
    elpds <- setNames(result_df$elpd, model_names)
    best_model <- names(which.max(elpds))
    result_df$elpd_rank <- rank(-result_df$elpd, ties.method = "min")
    result_df$elpd_best_model <- best_model
    result_df$delta_elpd_best <- result_df$elpd - max(result_df$elpd, na.rm = TRUE)
    result_df$delta_elpd_best_se <- NA_real_
    result_df$delta_elpd_best_z <- NA_real_

    pointwise_elpd <- lapply(loo_cache, function(x) x$pointwise[, "elpd_loo"])
    delta_se <- function(lhs, rhs) {
      if (!lhs %in% names(pointwise_elpd) || !rhs %in% names(pointwise_elpd)) return(NA_real_)
      diff <- pointwise_elpd[[lhs]] - pointwise_elpd[[rhs]]
      if (length(diff) <= 1L) return(NA_real_)
      sqrt(length(diff) * stats::var(diff, na.rm = TRUE))
    }
    delta_est <- function(lhs, rhs) {
      if (!lhs %in% names(elpds) || !rhs %in% names(elpds)) return(NA_real_)
      elpds[[lhs]] - elpds[[rhs]]
    }

    for (m in model_names) {
      idx <- result_df$fit_model == m
      se <- if (m == best_model) 0 else delta_se(m, best_model)
      result_df$delta_elpd_best_se[idx] <- se
      result_df$delta_elpd_best_z[idx] <- if (is.finite(se) && se > 0) abs(result_df$delta_elpd_best[idx] / se) else NA_real_
    }

    result_df$delta_looic_best <- -2 * result_df$delta_elpd_best
    result_df$delta_looic_best[abs(result_df$delta_looic_best) < 1e-12] <- 0
    result_df$delta_looic_best_se <- 2 * result_df$delta_elpd_best_se
    result_df$delta_looic_best_se[abs(result_df$delta_looic_best_se) < 1e-12] <- 0
    result_df$delta_looic_best_z <- result_df$delta_elpd_best_z

    pair_specs <- list(
      wst_minus_sst = c("WST", "SST"),
      wst_minus_dc = c("WST", "DC-SBM"),
      sst_minus_dc = c("SST", "DC-SBM")
    )
    for (nm in names(pair_specs)) {
      lhs <- pair_specs[[nm]][1L]
      rhs <- pair_specs[[nm]][2L]
      est <- delta_est(lhs, rhs)
      se <- delta_se(lhs, rhs)
      result_df[[paste0("elpd_", nm)]] <- est
      result_df[[paste0("elpd_", nm, "_se")]] <- se
      result_df[[paste0("elpd_", nm, "_z")]] <- if (is.finite(se) && se > 0) abs(est / se) else NA_real_
      result_df[[paste0("looic_", nm)]] <- -2 * est
      result_df[[paste0("looic_", nm)]][abs(result_df[[paste0("looic_", nm)]]) < 1e-12] <- 0
      result_df[[paste0("looic_", nm, "_se")]] <- 2 * se
      result_df[[paste0("looic_", nm, "_se")]][abs(result_df[[paste0("looic_", nm, "_se")]]) < 1e-12] <- 0
      result_df[[paste0("looic_", nm, "_z")]] <- if (is.finite(se) && se > 0) abs(est / se) else NA_real_
    }

    if (isTRUE(demo_save_loo_pointwise) && length(pointwise_elpd)) {
      pointwise_dir <- file.path(run_dir, "loo_pointwise")
      dir.create(pointwise_dir, recursive = TRUE, showWarnings = FALSE)
      pw <- data.frame(
        gen_model = gen_model,
        K_true = K_true,
        density = density_label,
        hierch = hierch_label,
        rep_id = rep_id,
        scenario_id = scenario_id,
        theta_ocrp = theta_ocrp,
        dgp_partition = demo_dgp_partition,
        dyad_id = seq_len(nrow(DI)),
        i = DI[, 1],
        j = DI[, 2],
        Aij = A[cbind(DI[, 1], DI[, 2])],
        Aji = A[cbind(DI[, 2], DI[, 1])],
        Nij = A[cbind(DI[, 1], DI[, 2])] + A[cbind(DI[, 2], DI[, 1])],
        zi_true = z_true[DI[, 1]],
        zj_true = z_true[DI[, 2]],
        same_true = z_true[DI[, 1]] == z_true[DI[, 2]],
        stringsAsFactors = FALSE
      )
      for (m in names(loo_cache)) {
        key <- gsub("[^A-Za-z0-9]+", "_", m)
        pw[[paste0("elpd_", key)]] <- loo_cache[[m]]$pointwise[, "elpd_loo"]
        pw[[paste0("pareto_k_", key)]] <- loo_cache[[m]]$loo$diagnostics$pareto_k
      }
      if ("WST" %in% names(pointwise_elpd) && "SST" %in% names(pointwise_elpd)) {
        pw$delta_elpd_WST_minus_SST <- pointwise_elpd[["WST"]] - pointwise_elpd[["SST"]]
      }
      if ("WST" %in% names(pointwise_elpd) && "DC-SBM" %in% names(pointwise_elpd)) {
        pw$delta_elpd_WST_minus_DC <- pointwise_elpd[["WST"]] - pointwise_elpd[["DC-SBM"]]
      }
      if ("SST" %in% names(pointwise_elpd) && "DC-SBM" %in% names(pointwise_elpd)) {
        pw$delta_elpd_SST_minus_DC <- pointwise_elpd[["SST"]] - pointwise_elpd[["DC-SBM"]]
      }
      pointwise_file <- file.path(
        pointwise_dir,
        sprintf("loo_pointwise_rep%03d_%s_%s_K%d_%s_%s_seed%d.csv",
                ifelse(is.na(rep_id), 0L, as.integer(rep_id)),
                ifelse(nzchar(scenario_id), scenario_id, "scenario"),
                gen_model, K_true,
                ifelse(is.na(density_label), "density", density_label),
                hierch_label, seed)
      )
      write.csv(pw, pointwise_file, row.names = FALSE)
      result_df$loo_pointwise_file <- pointwise_file
    }
  }

  result_df
}

############################################################
## 4. Simple file lock for progressive writing
############################################################

# Very simple directory-based lock: no extra packages needed.
with_ilock <- function(lock_base, code) {
  lock_dir <- paste0(lock_base, ".lock")
  repeat {
    ok <- dir.create(lock_dir, showWarnings = FALSE)
    if (ok) break
    Sys.sleep(runif(1, 0.01, 0.05))
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  code()
}

append_results_with_lock <- function(res_df, out_file) {
  with_ilock(out_file, function() {
    if (file.exists(out_file)) {
      old <- read.csv(out_file, stringsAsFactors = FALSE)
      new <- dplyr::bind_rows(old, res_df)
    } else {
      new <- res_df
    }
    write.csv(new, out_file, row.names = FALSE)
  })
}

############################################################
## 5. Simulation grid (mean/var scenarios)
############################################################

# Reduced main-text OCRP-DGP grid:
# - the partition is sampled from OCRP(vartheta = 3) with a small minimum
#   block size, so K* is realized per dataset but singleton ranks are avoided
# - hard: moderate totals and identifiable but nontrivial directional signal
# - easy: dense totals and strong directional signal
# - both generative mechanisms (SST/WST)
K_grid <- c(NA_integer_)
kappa_mean_grid <- c(1.40, 4.00)
kappa_sd_grid <- c(0.45, 1.00)
kappa_labels <- c(`1.4` = "moderate", `4` = "dense")

hierarchy_specs <- data.frame(
  hierch = c("hard", "easy"),
  psi_mean = c(0.45, 1.40),
  psi_sd = c(0.15, 0.25),
  stringsAsFactors = FALSE
)

gen_models <- c("SST", "WST", "DC-SBM")
n_rep      <- 5
n_items    <- 60

K_grid <- as.integer(parse_numeric_vec_env("DEMOOCRPVAR_K_GRID", K_grid))
kappa_mean_grid <- parse_numeric_vec_env("DEMOOCRPVAR_KAPPA_MEAN_GRID", kappa_mean_grid)
kappa_sd_grid <- parse_numeric_vec_env("DEMOOCRPVAR_KAPPA_SD_GRID", kappa_sd_grid)
hierarchy_specs <- parse_hierarchy_specs_env(hierarchy_specs)
gen_models <- parse_char_vec_env("DEMOOCRPVAR_GEN_MODELS", gen_models, allowed = c("SST", "WST", "DC-SBM"))
scenario_specs <- parse_scenario_specs_env()
demo_factorial_grid <- tolower(Sys.getenv("DEMOOCRPVAR_FACTORIAL_GRID", unset = "false")) %in% c("1", "true", "yes")

run_start_time <- Sys.time()
run_stamp <- format(run_start_time, "%Y%m%d_%H%M%S")
run_id <- paste0("DemoOCRPvar_run_", run_stamp)
run_dir <- file.path("output", "simulation", "raw", "DemoOCRPvar_runs", run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
message("Run directory: ", run_dir)

if (tolower(Sys.getenv("DEMOOCRPVAR_SMOKE", unset = "false")) %in% c("1", "true", "yes")) {
  K_grid <- c(4)
  kappa_mean_grid <- c(1)
  kappa_sd_grid <- c(0.25)
  kappa_labels <- c(`1` = "smoke")
  hierarchy_specs <- data.frame(
    hierch = "weak",
    psi_mean = 0.25,
    psi_sd = 0.30,
    stringsAsFactors = FALSE
  )
  gen_models <- c("WST")
  scenario_specs <- data.frame(
    scenario_label = "smoke",
    gen_model = "WST",
    K_true = NA_integer_,
    kappa_mean = 1,
    kappa_sd = 0.25,
    hierch = "weak",
    psi_mean = 0.25,
    psi_sd = 0.30,
    stringsAsFactors = FALSE
  )
  n_rep <- 1
  n_items <- 24
}

n_rep_env <- suppressWarnings(as.integer(Sys.getenv("DEMOOCRPVAR_N_REP", unset = "")))
if (!is.na(n_rep_env) && n_rep_env >= 1L) n_rep <- n_rep_env

n_items_env <- suppressWarnings(as.integer(Sys.getenv("DEMOOCRPVAR_N_ITEMS", unset = "")))
if (!is.na(n_items_env) && n_items_env >= 4L) n_items <- n_items_env

demo_n_iter <- suppressWarnings(as.integer(Sys.getenv("DEMOOCRPVAR_N_ITER", unset = "10000")))
demo_burn   <- suppressWarnings(as.integer(Sys.getenv("DEMOOCRPVAR_BURN", unset = "1000")))
demo_thin   <- suppressWarnings(as.integer(Sys.getenv("DEMOOCRPVAR_THIN", unset = "2")))

if (is.na(demo_n_iter) || demo_n_iter < 1L) demo_n_iter <- 10000L
if (is.na(demo_burn) || demo_burn < 0L || demo_burn >= demo_n_iter) demo_burn <- max(0L, floor(demo_n_iter / 5))
if (is.na(demo_thin) || demo_thin < 1L) demo_thin <- 2L

if (!is.null(scenario_specs)) {
  grid_base <- scenario_specs %>%
    arrange(gen_model, K_true, kappa_mean, hierch) %>%
    mutate(
      psi_var = psi_sd^2,
      kappa_sd = ifelse(is.finite(kappa_sd), kappa_sd, 0.35),
      kappa_var = kappa_sd^2,
      density = unname(kappa_labels[as.character(kappa_mean)]),
      density = ifelse(is.na(density), paste0("kappa_", kappa_mean), density),
      theta_ocrp = if (is.finite(demo_theta_override)) {
        demo_theta_override
      } else {
        vapply(K_true, function(k) choose_theta_for_expected_K(n_items, k), numeric(1))
      },
      rate_kap = kappa_mean,
      scenario_id = ifelse(nzchar(scenario_label), scenario_label,
                           sprintf("scn_%03d", dplyr::row_number()))
    )
} else if (isTRUE(demo_factorial_grid)) {
  grid_base <- expand.grid(
    K_true     = K_grid,
    kappa_mean = kappa_mean_grid,
    kappa_sd   = kappa_sd_grid,
    gen_model  = gen_models,
    stringsAsFactors = FALSE
  ) %>%
    tidyr::crossing(hierarchy_specs) %>%
    arrange(gen_model, K_true, kappa_mean, kappa_sd, factor(hierch, levels = c("weak", "strong"))) %>%
    mutate(
      psi_var   = psi_sd^2,
      kappa_var = kappa_sd^2,
      density = unname(kappa_labels[as.character(kappa_mean)]),
      density = ifelse(is.na(density), paste0("kappa_", kappa_mean), density),
      theta_ocrp = if (is.finite(demo_theta_override)) {
        demo_theta_override
      } else {
        vapply(K_true, function(k) choose_theta_for_expected_K(n_items, k), numeric(1))
      },
      rate_kap  = kappa_mean,         # keep old name for compatibility
      scenario_id = sprintf("scn_%03d", dplyr::row_number())
    )
} else {
  grid_base <- data.frame(
    scenario_label = c("hard_sst_ocrp", "easy_sst_ocrp", "hard_wst_ocrp", "easy_wst_ocrp", "dcsbm_ocrp"),
    gen_model = c("SST", "SST", "WST", "WST", "DC-SBM"),
    K_true = NA_integer_,
    kappa_mean = c(1.40, 4.00, 1.40, 4.00, 1.60),
    kappa_sd = c(0.45, 1.00, 0.45, 1.00, 0.40),
    density = c("moderate", "dense", "moderate", "dense", "moderate"),
    hierch = c("hard", "easy", "hard", "easy", "misspecified"),
    psi_mean = c(0.45, 0.70, 0.45, 1.40, 0.70),
    psi_sd = c(0.12, 0.15, 0.15, 0.25, 0.20),
    stringsAsFactors = FALSE
  ) %>%
    filter(gen_model %in% gen_models) %>%
    arrange(gen_model, factor(hierch, levels = c("hard", "easy", "misspecified"))) %>%
    mutate(
      psi_var = psi_sd^2,
      kappa_var = kappa_sd^2,
      theta_ocrp = if (is.finite(demo_theta_override)) demo_theta_override else 3,
      rate_kap = kappa_mean,
      scenario_id = scenario_label
    )
}

grid_val <- grid_base %>%
  select(K_true, kappa_mean, kappa_sd, density, hierch, gen_model, psi_mean, psi_var,
         kappa_var, theta_ocrp, rate_kap, scenario_id)

message(
  sprintf(
    "Grid: %d scenarios x %d reps = %d tasks",
    nrow(grid_val), n_rep, nrow(grid_val) * n_rep
  )
)

# Task grid: replicate x scenario
task_grid <- expand.grid(
  rep_id = seq_len(n_rep),
  idx    = seq_len(nrow(grid_val)),
  stringsAsFactors = FALSE
) %>%
  arrange(rep_id, idx)

out_file <- file.path(run_dir, paste0("full_simulation_crossfit_progress_", run_id, ".csv"))

############################################################
## 6. Parallel execution with progressive writing
############################################################

slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")))
slurm_ntasks <- suppressWarnings(as.integer(Sys.getenv("SLURM_NTASKS", unset = "")))
env_cores <- suppressWarnings(as.integer(Sys.getenv("DEMOOCRPVAR_CORES", unset = "")))

if (!is.na(env_cores) && env_cores >= 1L) {
  n_cores <- env_cores
} else if (!is.na(slurm_cpus) && slurm_cpus >= 1L) {
  n_cores <- slurm_cpus
} else if (!is.na(slurm_ntasks) && slurm_ntasks >= 1L) {
  n_cores <- slurm_ntasks
} else {
  n_cores <- 4L
}

n_cores <- max(1L, min(n_cores, nrow(task_grid)))

message("Using ", n_cores, " cores.")

run_task <- function(tid) {
  rep_id <- task_grid$rep_id[tid]
  row_id <- task_grid$idx[tid]
  row    <- grid_val[row_id, ]
  
  seed <- 1000 + tid
  
  cat(sprintf(
    "Rep %d | Gen %s | K target %s | density %s | kappa_mean %.2f | kappa_sd %.2f | hierch %s | psi_mean %.2f | theta %.3f\n",
    rep_id, row$gen_model,
    ifelse(is.na(row$K_true), "random", as.character(row$K_true)),
    row$density, row$kappa_mean, row$kappa_sd, row$hierch,
    row$psi_mean, row$theta_ocrp
  ))
  
  res_df <- tryCatch(
    run_replicate_crossfit(
      gen_model   = row$gen_model,
      K_true      = row$K_true,
      n           = n_items,
      kappa_mean  = row$kappa_mean,
      kappa_var   = row$kappa_var,
      psi_mean    = row$psi_mean,
      psi_var     = row$psi_var,
      hierch_label = row$hierch,
      theta_ocrp  = row$theta_ocrp,
      seed        = seed,
      rep_id      = rep_id,
      scenario_id = row$scenario_id,
      density_label = row$density,
      n_iter      = demo_n_iter,
      burn        = demo_burn,
      thin        = demo_thin
    ),
    error = function(e) {
      message(
        sprintf(
          "Error in tid=%d, rep=%d, gen=%s, K=%d, kappa_mean=%.2f, hierch=%s:\n  %s",
          tid, rep_id, row$gen_model, row$K_true,
          row$kappa_mean, row$hierch,
          conditionMessage(e)
        )
      )
      NULL
    }
  )
  
  if (!is.null(res_df)) {
    res_df$rep_id <- rep_id
    res_df$kappa_mean <- row$kappa_mean
    res_df$kappa_sd   <- row$kappa_sd
    res_df$kappa_var  <- row$kappa_var
    res_df$psi_mean   <- row$psi_mean
    res_df$psi_var    <- row$psi_var
    res_df$hierch     <- row$hierch
    res_df$density    <- row$density
    res_df$theta_ocrp <- row$theta_ocrp
    res_df$rate_kap   <- row$kappa_mean
    res_df$eta_mean   <- demo_eta_mean
    res_df$eta_sd     <- demo_eta_sd
    res_df$eta_mode   <- demo_eta_mode
    
    append_results_with_lock(res_df, out_file)
  }
  
  res_df
}

run_task_batch <- function(tids) {
  if (!length(tids)) return(list())
  mclapply(
    X = tids,
    FUN = run_task,
    mc.cores = max(1L, min(n_cores, length(tids))),
    mc.preschedule = FALSE
  )
}

first_rep_tids <- which(task_grid$rep_id == 1L)
remaining_tids <- setdiff(seq_len(nrow(task_grid)), first_rep_tids)
snapshot_file <- NA_character_
timing_estimate_file <- NA_character_

if (demo_first_rep_snapshot) {
  message("Running first repetition snapshot batch: ", length(first_rep_tids), " tasks.")
  first_rep_start <- Sys.time()
  first_rep_list <- run_task_batch(first_rep_tids)
  first_rep_elapsed <- as.numeric(difftime(Sys.time(), first_rep_start, units = "secs"))

  snapshot_file <- file.path(run_dir, paste0("snapshot_after_rep1_", run_id, ".csv"))
  first_rep_table <- dplyr::bind_rows(first_rep_list)
  write.csv(first_rep_table, snapshot_file, row.names = FALSE)

  first_rep_waves <- ceiling(length(first_rep_tids) / n_cores)
  remaining_waves <- if (length(remaining_tids)) ceiling(length(remaining_tids) / n_cores) else 0L
  seconds_per_wave <- first_rep_elapsed / max(first_rep_waves, 1L)
  estimated_total_seconds <- seconds_per_wave * (first_rep_waves + remaining_waves)
  estimated_remaining_seconds <- seconds_per_wave * remaining_waves
  estimate_lines <- c(
    paste0("run_id: ", run_id),
    paste0("first_rep_tasks: ", length(first_rep_tids)),
    paste0("cores: ", n_cores),
    paste0("first_rep_parallel_waves: ", first_rep_waves),
    paste0("remaining_parallel_waves: ", remaining_waves),
    paste0("seconds_per_wave_estimate: ", round(seconds_per_wave, 3)),
    paste0("first_rep_elapsed_seconds: ", round(first_rep_elapsed, 3)),
    paste0("estimated_total_seconds: ", round(estimated_total_seconds, 3)),
    paste0("estimated_remaining_seconds_after_rep1: ", round(estimated_remaining_seconds, 3)),
    paste0("estimated_total_hours: ", round(estimated_total_seconds / 3600, 3)),
    paste0("estimated_remaining_hours_after_rep1: ", round(estimated_remaining_seconds / 3600, 3)),
    paste0("snapshot_csv: ", snapshot_file)
  )
  timing_estimate_file <- file.path(run_dir, "timing_estimate_after_rep1.txt")
  writeLines(estimate_lines, con = timing_estimate_file)
  message(sprintf(
    "First repetition took %.1f min; estimated total %.1f h, remaining %.1f h.",
    first_rep_elapsed / 60,
    estimated_total_seconds / 3600,
    estimated_remaining_seconds / 3600
  ))

  message("Running remaining repetitions: ", length(remaining_tids), " tasks.")
  remaining_list <- run_task_batch(remaining_tids)
  res_list <- c(first_rep_list, remaining_list)
} else {
  message("Running all tasks in a single batch: ", nrow(task_grid), " tasks.")
  res_list <- run_task_batch(seq_len(nrow(task_grid)))
}

final_table <- dplyr::bind_rows(res_list)

# Also write final combined table (in case you prefer a clean version)
final_file <- file.path(run_dir, paste0("full_simulation_crossfit_final_", run_id, ".csv"))
write.csv(final_table, final_file, row.names = FALSE)

git_hash <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)[1],
  error = function(e) NA_character_
)
if (!is.character(git_hash) || !length(git_hash) || is.na(git_hash) || !nzchar(git_hash)) {
  git_hash <- "NA"
}

format_K_target_grid <- function(x) {
  vals <- sort(unique(x))
  vals <- vals[!is.na(vals)]
  if (!length(vals)) return("random")
  paste(vals, collapse = ",")
}

format_scenario_specs_manifest <- function() {
  env_specs <- Sys.getenv("DEMOOCRPVAR_SCENARIO_SPECS", unset = "")
  if (nzchar(trimws(env_specs))) return(env_specs)
  if (!is.null(scenario_specs)) return(paste(unique(grid_val$scenario_id), collapse = ","))
  if (isTRUE(demo_factorial_grid)) return("factorial grid")
  "default OCRP main grid"
}

run_end_time <- Sys.time()
manifest_lines <- c(
  paste0("run_id: ", run_id),
  paste0("script: DemoOCRPvar.R"),
  paste0("dgp_partition: ", dplyr::case_when(
    demo_dgp_partition == "BALANCED" ~ "balanced block sizes",
    demo_dgp_partition == "OCRP" ~ "OCRP direct draw with random K_star",
    demo_dgp_partition == "OCRP_MIN" ~ paste0("OCRP direct draw with random K_star and min block size ", demo_min_block_size),
    demo_dgp_partition == "OCRP_FIXED_K_MIN" ~ paste0("OCRP conditioned on K_target with min block size ", demo_min_block_size),
    TRUE ~ "OCRP conditioned on K_target"
  )),
  paste0("min_block_size: ", demo_min_block_size),
  paste0("run_dir: ", run_dir),
  paste0("git_commit_short: ", git_hash),
  paste0("start_time: ", format(run_start_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("end_time: ", format(run_end_time, "%Y-%m-%d %H:%M:%S %Z")),
  paste0("elapsed_seconds: ", round(as.numeric(difftime(run_end_time, run_start_time, units = "secs")), 3)),
  paste0("smoke_mode: ", tolower(Sys.getenv("DEMOOCRPVAR_SMOKE", unset = "false"))),
  paste0("gen_models: ", paste(unique(grid_val$gen_model), collapse = ",")),
  paste0("scenario_specs: ", format_scenario_specs_manifest()),
  paste0("scenario_count: ", nrow(grid_val)),
  paste0("K_target_grid: ", format_K_target_grid(grid_val$K_true)),
  paste0("kappa_mean_grid: ", paste(sort(unique(grid_val$kappa_mean)), collapse = ",")),
  paste0("density_levels: ", paste(unique(grid_val$density), collapse = ",")),
  paste0("hierarchy_levels: ", paste(unique(grid_val$hierch), collapse = ",")),
  paste0("sst_psi_mode: ", demo_sst_psi_mode),
  paste0("n_rep: ", n_rep),
  paste0("n_items: ", n_items),
  paste0("eta_mode: ", demo_eta_mode),
  paste0("n_iter: ", demo_n_iter),
  paste0("burn: ", demo_burn),
  paste0("thin: ", demo_thin),
  paste0("K_init_ordered_models: n_items"),
  paste0("K_init_dcsbm: n_items"),
  paste0("K_expected_policy: realized K_star per replicate"),
  paste0("hyper_mode: ", if (demo_use_fixed_hypers) "fixed" else "calibrated"),
  paste0("first_rep_snapshot_enabled: ", demo_first_rep_snapshot),
  paste0("tau0_floor: ", demo_tau0_floor),
  paste0("tau0_override: ", if (is.finite(demo_tau0_override)) demo_tau0_override else "none"),
  paste0("sigma0_override: ", if (is.finite(demo_sigma0_override)) demo_sigma0_override else "none"),
  paste0("dominance_target: ", demo_order_dominance_target),
  paste0("N0_psi: ", demo_order_N0_psi),
  paste0("a_kappa_fixed: ", demo_a_kappa_fixed),
  paste0("b_kappa_fixed: ", demo_b_kappa_fixed),
  paste0("a_eta_fixed: ", demo_a_eta_fixed),
  paste0("b_eta_fixed: ", demo_b_eta_fixed),
  paste0("mu0_wst: ", demo_mu0_wst),
  paste0("sigma0_wst: ", demo_sigma0_wst),
  paste0("mu0_sst: ", demo_mu0_sst),
  paste0("sigma0_sst: ", demo_sigma0_sst),
  paste0("tau0_sst: ", demo_tau0_sst),
  paste0("gamma_gn_fixed: ", demo_gamma_gn_fixed),
  paste0("ordering_prior_mode: equivalence_class"),
  paste0("osbm_partition_prior: ", demo_osbm_partition_prior),
  paste0("theta_ocrp_policy: ", if (is.finite(demo_theta_override)) "fixed" else "per_K_expected"),
  paste0("theta_ocrp: ", paste(sort(unique(grid_val$theta_ocrp)), collapse = ",")),
  paste0("alpha_crp_policy: ", if (is.finite(demo_alpha_crp)) paste0("fixed=", demo_alpha_crp) else "match_theta_ocrp"),
  paste0("cores: ", n_cores),
  paste0("save_loo_pointwise: ", demo_save_loo_pointwise),
  paste0("first_rep_snapshot_csv: ", ifelse(is.na(snapshot_file), "disabled", snapshot_file)),
  paste0("timing_estimate_file: ", ifelse(is.na(timing_estimate_file), "disabled", timing_estimate_file)),
  paste0("progress_csv: ", out_file),
  paste0("final_csv: ", final_file),
  paste0("n_rows_final: ", nrow(final_table))
)
writeLines(manifest_lines, con = file.path(run_dir, "run_manifest.txt"))

############################################################
## 7. Post-processing (same spirit as before)
############################################################
# 
# # Clean types
# final_table <- final_table %>%
#   select(-starts_with("Unnamed")) %>%
#   mutate(
#     gen_model  = factor(gen_model,  levels = c("SST", "WST")),
#     fit_model  = factor(fit_model,  levels = c("SST", "WST", "DC-SBM")),
#     hierch     = factor(hierch,     levels = c("weak", "strong")),
#     K_true     = as.integer(K_true),
#     K_hat      = as.integer(K_hat),
#     rate_kap   = as.numeric(rate_kap),
#     kappa_mean = as.numeric(kappa_mean),
#     kappa_var  = as.numeric(kappa_var),
#     psi_mean   = as.numeric(psi_mean),
#     psi_var    = as.numeric(psi_var),
#     seed       = as.integer(seed)
#   )
# 
# # Main performance summary
# summary_main <- final_table %>%
#   group_by(gen_model, hierch, rate_kap, K_true, fit_model) %>%
#   summarise(
#     n_runs       = n(),
#     mean_K_hat   = mean(K_hat, na.rm = TRUE),
#     sd_K_hat     = sd(K_hat, na.rm = TRUE),
#     prop_K_true  = mean(K_hat == K_true, na.rm = TRUE),
#     prop_K_under = mean(K_hat < K_true, na.rm = TRUE),
#     prop_K_over  = mean(K_hat > K_true, na.rm = TRUE),
#     
#     mean_ari     = mean(ari, na.rm = TRUE),
#     sd_ari       = sd(ari, na.rm = TRUE),
#     mean_vi      = mean(vi, na.rm = TRUE),
#     
#     mean_elpd    = mean(elpd, na.rm = TRUE),
#     sd_elpd      = sd(elpd, na.rm = TRUE),
#     mean_looic   = mean(looic, na.rm = TRUE),
#     mean_pk_bad  = mean(pk_bad, na.rm = TRUE),
#     
#     .groups = "drop"
#   )
# 
# summary_main
# 
# # Diagnostics summary
# summary_diag <- final_table %>%
#   group_by(gen_model, hierch, rate_kap, K_true, fit_model) %>%
#   summarise(
#     mean_violation_rate = mean(violation_rate_mean, na.rm = TRUE),
#     sd_violation_rate   = sd(violation_rate_mean, na.rm = TRUE),
#     mean_viol_count     = mean(violation_count_mean, na.rm = TRUE),
#     mean_cross_mass     = mean(cross_mass_mean, na.rm = TRUE),
#     mean_cov_block_model = mean(coverage_block_model_avg, na.rm = TRUE),
#     mean_cov_block_emp   = mean(coverage_block_emp_avg,   na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# summary_diag
# 
# ############################################################
# ## 8. Example tables (same idea as your old ones)
# ############################################################
# 
# make_latex_table <- function(gen, K, hier) {
#   summary_main %>%
#     filter(gen_model == gen,
#            K_true    == K,
#            hierch    == hier) %>%
#     arrange(rate_kap, fit_model) %>%
#     select(rate_kap, fit_model,
#            mean_K_hat, prop_K_true,
#            mean_ari, mean_elpd, mean_pk_bad) %>%
#     kable(
#       format  = "latex",
#       booktabs = TRUE,
#       digits  = c(1, 1, 2, 2, 2, 1, 2),
#       caption = paste0(
#         gen, " generativo, $K^\\star=", K, "$, gerarchia ", hier,
#         ": recupero di $K$, ARI ed ELPD per modello di fit (prior mean/var)."
#       )
#     )
# }
# 
# # Examples:
# make_latex_table("SST", 3, "strong")
# make_latex_table("SST", 5, "strong")
# make_latex_table("WST", 3, "strong")
# make_latex_table("WST", 5, "strong")
# 
# tab_elpd_wide <- summary_main %>%
#   select(gen_model, hierch, K_true, rate_kap, fit_model, mean_elpd) %>%
#   mutate(
#     fit_model = as.character(fit_model),
#     rate_kap  = as.factor(rate_kap)
#   ) %>%
#   pivot_wider(
#     names_from  = fit_model,
#     values_from = mean_elpd
#   ) %>%
#   arrange(gen_model, hierch, K_true, rate_kap)
# 
# kable(
#   tab_elpd_wide,
#   format  = "latex",
#   booktabs = TRUE,
#   digits  = 1,
#   caption = "ELPD medio (prior mean/var) per combinazione di modello generativo, gerarchia, $K^\\star$ e densità (kappa_mean), a confronto tra modelli di fit."
# )
# 
# tab_violation <- summary_diag %>%
#   arrange(gen_model, hierch, K_true, rate_kap, fit_model) %>%
#   select(gen_model, hierch, K_true, rate_kap, fit_model,
#          mean_violation_rate, mean_viol_count,
#          mean_cross_mass, mean_cov_block_model)
# 
# kable(
#   tab_violation,
#   format  = "latex",
#   booktabs = TRUE,
#   digits  = c(0, 0, 0, 1, 0, 2, 1, 1, 2),
#   caption = "Statistiche di violazione dell'ordine e coverage (prior mean/var) per combinazione di scenario e modello di fit."
# )
# 
# ############################################################
# ## 9. Plots (same structure as before)
# ############################################################
# 
# ggplot(summary_main,
#        aes(x = rate_kap, y = mean_K_hat,
#            color = fit_model, group = interaction(fit_model, K_true))) +
#   geom_line() +
#   geom_point(aes(shape = factor(K_true)), size = 2) +
#   facet_grid(gen_model ~ hierch) +
#   labs(
#     x = "kappa_mean (densità blocchi)",
#     y = "K_hat medio",
#     color = "Modello di fit",
#     shape = "K_true",
#     title = "Recupero di K per scenario generativo (prior mean/var)"
#   ) +
#   theme_bw()
# 
# ggplot(summary_main,
#        aes(x = rate_kap, y = prop_K_true,
#            color = fit_model, group = interaction(fit_model, K_true))) +
#   geom_line() +
#   geom_point(aes(shape = factor(K_true)), size = 2) +
#   facet_grid(gen_model ~ hierch) +
#   scale_y_continuous(limits = c(0, 1)) +
#   labs(
#     x = "kappa_mean",
#     y = "Pr(K_hat = K_true)",
#     color = "Modello di fit",
#     shape = "K_true",
#     title = "Probabilità di recupero esatto di K (prior mean/var)"
#   ) +
#   theme_bw()
# 
# ggplot(summary_main,
#        aes(x = rate_kap, y = mean_elpd,
#            color = fit_model, group = interaction(fit_model, K_true))) +
#   geom_line() +
#   geom_point(aes(shape = factor(K_true)), size = 2) +
#   facet_grid(gen_model ~ hierch) +
#   labs(
#     x = "kappa_mean",
#     y = "ELPD medio (loo)",
#     color = "Modello di fit",
#     shape = "K_true",
#     title = "Confronto predittivo tra modelli di fit (prior mean/var)"
#   ) +
#   theme_bw()
# 
# ggplot(summary_diag,
#        aes(x = rate_kap, y = mean_violation_rate,
#            color = fit_model, group = fit_model)) +
#   geom_line() +
#   geom_point(size = 2) +
#   facet_grid(gen_model ~ hierch) +
#   labs(
#     x = "kappa_mean",
#     y = "Violazione media (psi)",
#     color = "Modello di fit",
#     title = "Grado di violazione dell'ordine tra blocchi (prior mean/var)"
#   ) +
#   theme_bw()
