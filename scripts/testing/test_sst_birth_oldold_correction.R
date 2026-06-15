

# =============================================================================
# PART B — sampler benchmark: corrected vs uncorrected birth move
# =============================================================================
cat("\n================================================================\n")
cat(" PART B — sampler benchmark on simulated SST data\n")
cat("================================================================\n\n")

suppressPackageStartupMessages(library(mcclust))
source("helper_folder/sim_study_helper.R")
source("helper_folder/helper.R")
source("helper_folder/mixing_moves.R")
source("helper_folder/estimating_K_helper.R")
source("helper_folder/SST_helpers.R")
source("helper_folder/WST_helpers.R")
source("helper_folder/Hyper_setup.R")
source("core/my_best_try_so_far.R")

run_one_chain <- function(A_sim, K_init, n_iter, burn, seed, hypers, score_mode) {
  modular_osbm_sampler(
    A = as.matrix(A_sim),
    K = K_init,
    n_iter = n_iter, burn = burn, thin = 1,
    psi_constraint = "SST" ,
    seed = seed,
    verbose = FALSE,
    partition_prior = "GN",
    gamma_gn = hypers$gamma,
    a_kappa = 1, b_kappa = 1,
    a_eta = 1, b_eta = 1,
    mu0 = hypers$mu0,
    sigma0 = hypers$sigma0,
    tau0 = max(hypers$tau0, 0.2),
    eta_identifiability = "block_sum",
    sst_birth_score_mode = score_mode,
    use_mixing_moves = TRUE
  )
}

bench_one <- function(rep, n_nodes = 40, K_star = 4, n_iter = 600, burn = 200) {
  set.seed(1000 + rep)
  z_star <- rep(seq_len(K_star), length.out = n_nodes)
  z_star <- z_star[sample.int(n_nodes)]
  eta_star <- rep(1, n_nodes)
  kappa_star <- matrix(3.0, K_star, K_star)
  psi_star <- cumsum(c(0.6, 0.7, 0.8))     # length K_star-1
  A_sim <- simulate_osbm(n_nodes, K_star, z_star, eta_star, kappa_star,
                         psi_star, regime = "SST")
  hypers <- get_principled_hypers_v2(
    A = as.matrix(A_sim),
    K_expected = K_star,
    c_kappa = 3,
    ordering_prior_mode = "equivalence_class"
  )

  res_exact <- run_one_chain(A_sim, K_init = n_nodes, n_iter, burn,
                             seed = 2000 + rep, hypers = hypers,
                             score_mode = "exact_nonlocal")
  res_local <- run_one_chain(A_sim, K_init = n_nodes, n_iter, burn,
                             seed = 3000 + rep, hypers = hypers,
                             score_mode = "local_approx")

  pull_z_chain <- function(r) {
    if (!is.null(r$z)) return(r$z)
    if (!is.null(r$z_chain)) return(r$z_chain)
    stop("Sampler result has neither z nor z_chain.")
  }
  pull_K  <- function(r) vapply(pull_z_chain(r), function(z) length(unique(z)), integer(1))
  pull_vi <- function(r) {
    z_chain <- pull_z_chain(r)
    n_keep <- length(z_chain)
    keep_idx <- seq(max(1, floor(n_keep * 0.5)), n_keep)
    vi_vals <- vapply(z_chain[keep_idx],
                      function(z) mcclust::vi.dist(z, z_star),
                      numeric(1))
    mean(vi_vals)
  }
  list(
    rep = rep,
    K_mode_exact = as.integer(names(sort(table(pull_K(res_exact)), decreasing=TRUE))[1]),
    K_mode_local = as.integer(names(sort(table(pull_K(res_local)), decreasing=TRUE))[1]),
    K_mean_exact = mean(pull_K(res_exact)),
    K_mean_local = mean(pull_K(res_local)),
    vi_exact     = pull_vi(res_exact),
    vi_local     = pull_vi(res_local),
    K_star     = K_star
  )
}

B <- as.integer(Sys.getenv("SST_OLDOLD_B", unset = "3"))
n_iter_default <- as.integer(Sys.getenv("SST_OLDOLD_N_ITER", unset = "600"))
burn_default <- as.integer(Sys.getenv("SST_OLDOLD_BURN", unset = "200"))
cat(sprintf("Running %d replicates (n=40, K*=4, n_iter=%d, burn=%d) ...\n\n",
            B, n_iter_default, burn_default))
out_list <- vector("list", B)
for (b in seq_len(B)) {
  cat(sprintf("  rep %d/%d ... ", b, B)); flush.console()
  t0 <- Sys.time()
  out_list[[b]] <- bench_one(b, n_iter = n_iter_default, burn = burn_default)
  cat(sprintf("done (%.1fs)\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

res_df <- do.call(rbind, lapply(out_list, function(x) {
  data.frame(rep = x$rep,
             K_star = x$K_star,
             K_mode_exact = x$K_mode_exact, K_mode_local = x$K_mode_local,
             K_mean_exact = x$K_mean_exact, K_mean_local = x$K_mean_local,
             vi_exact = x$vi_exact,         vi_local = x$vi_local)
}))

cat("\nPer-replicate results:\n")
print(res_df, row.names = FALSE, digits = 4)

cat("\nSummary across replicates:\n")
cat(sprintf("  K_mode  exact_nonlocal: %s\n", paste(res_df$K_mode_exact, collapse = " ")))
cat(sprintf("  K_mode  local_approx  : %s\n", paste(res_df$K_mode_local, collapse = " ")))
cat(sprintf("  mean VI exact_nonlocal: %.3f (sd %.3f)\n",
            mean(res_df$vi_exact),  sd(res_df$vi_exact)))
cat(sprintf("  mean VI local_approx  : %.3f (sd %.3f)\n",
            mean(res_df$vi_local), sd(res_df$vi_local)))

wins_K  <- sum(abs(res_df$K_mode_exact - res_df$K_star) <
               abs(res_df$K_mode_local - res_df$K_star))
ties_K  <- sum(abs(res_df$K_mode_exact - res_df$K_star) ==
               abs(res_df$K_mode_local - res_df$K_star))
wins_vi <- sum(res_df$vi_exact < res_df$vi_local)

cat(sprintf("\nExact nonlocal closer-to-truth K* mode in %d/%d replicates (ties: %d).\n",
            wins_K, B, ties_K))
cat(sprintf("Exact nonlocal lower VI to truth in %d/%d replicates.\n", wins_vi, B))

cat("\n=== ALL CHECKS COMPLETE ===\n")
