# ======================================================================
# Experiment driver: LOOIC + new order-conformity diagnostics (summarized)
# Author: Lapo (cleaned pipeline) — updated to drop per-draw storage
# ======================================================================

# OBSOLETE: this legacy simulation driver has been superseded by
# scripts/simulation/DemoKvar.R and is kept only for historical reference.

suppressPackageStartupMessages({
  library(Matrix)
  library(mcclust)
  library(mcclust.ext)   # minVI
  library(coda)
  library(dplyr)
  library(tidyr)
  library(BayesLogit)
  library(truncnorm)
  library(loo)
  library(fs)
  library(salso)
  library(readr)
  library(glue)
  library(future.apply)
  library(purrr)
})

plan(multisession, workers = 1)
RNGkind("L'Ecuyer-CMRG")
set.seed(42)



# ---------------- Main simulation driver -------------------------------

run_simulation_study <- function(
    K_vals = c(6),
    n_reps = 1, n = 100,
    n_iter = 10000, burn = 2000, thin = 1,
    gen_models = c("WST"),
    fit_models = c("WST"),
    verbose = TRUE, save_plots = FALSE, seed = 1234,
    m_item_triplets = 2000L       # per-draw item triple budget
){
  set.seed(seed)
  difficulties <- list(
    easy = list(rate_scale = 1, psi_scale = 1.0)
  )
  jobs <- expand.grid(
    K = K_vals,
    diff = names(difficulties),
    gen_model = gen_models,
    rep = seq_len(n_reps),
    stringsAsFactors = FALSE
  )
  
  job_results <- future.apply::future_lapply(seq_len(nrow(jobs)), function(j) {
    Sys.setenv("OMP_NUM_THREADS"="1","MKL_NUM_THREADS"="1",
               "OPENBLAS_NUM_THREADS"="1","VECLIB_MAXIMUM_THREADS"="1")
    options(mc.cores = 1)
    
    source("./core/modular mcmc.R")
    source("./helper_folder/helper.R")
    source("./helper_folder/transitivity_check_helper.R")
    source("./helper_folder/sim_study_helper.R")
    source("./core/DDCSBM.R")
    
    seed_job <- seed + j * 1009L
    set.seed(seed_job)
    
    K         <- jobs$K[j]
    diff      <- jobs$diff[j]
    gen_model <- jobs$gen_model[j]
    rep       <- jobs$rep[j]
    pars      <- difficulties[[diff]]
    
    # ----- Generate truth & data -----
    z_t   <- rep(1:K, length.out = n)
    eta_t <- runif(n, 0.8, 1.2); eta_t <- n * eta_t / sum(eta_t)
    kappa_t <- make_kappa(K, rate_scale = 1)
    
    if (gen_model == "WST") {
      psi_wst_t <- make_psi_wst(K, scale = pars$psi_scale)
      A <- simulate_osbm(n, K, z_t, eta_t, kappa_t, psi = psi_wst_t, regime = "WST")
      truth <- list(z = z_t, eta = eta_t, kappa = kappa_t, psi = psi_wst_t)
      A<- as.matrix(A)
      } else {
      psi_sst_t <- make_psi_sst(K, scale = pars$psi_scale)
      A <- simulate_osbm(n, K, z_t, eta_t, kappa_t, psi = psi_sst_t, regime = "SST")
      truth <- list(z = z_t, eta = eta_t, kappa = kappa_t, psi = psi_sst_t)
      A<- as.matrix(A)
      }
    
    rows <- list()
    
    nA  <- nrow(A)
    DI  <- build_dyad_index(nA)
    
    for (fit_model in fit_models) {
      if (verbose) cat(sprintf("[job %d/%d] GEN=%s | K=%d | %s | rep=%d\n",
                               j, nrow(jobs), gen_model, K, fit_model, rep))
      
      if (fit_model == "DC-SBM") {
        tm <- system.time({
          fit <- fit_dcsbm(as.matrix(A), K, iters = n_iter, burn = burn, thin = thin,
                           verbose = if (isTRUE(verbose)) 250 else FALSE,
                           seed = seed_job + rep)
        })
        
        # NEW: relabel by success (same as WST)
        fit_relab <- relabel_by_ECR(mcmc_out = fit,model = fit_model)
        
        # Use relabelled z for partition estimate
        z_hat <- estimate_partition(fit_relab$z)
        ari   <- fossil::adj.rand.index(z_hat, z_t)
        vi    <- mcclust::vi.dist(z_hat, z_t)
        
        # LOO on relabelled chains (mathematically identical, but keeps shapes aligned)
        LL <- loglik_matrix_dcsbm(as.matrix(A), fit, dyad_index = DI)
        loo_fit <- loo::loo(LL)
        elpd_loo <- loo_fit$estimates["elpd_loo","Estimate"]
        p_loo    <- loo_fit$estimates["p_loo","Estimate"]
        looic    <- loo_fit$estimates["looic","Estimate"]
        pareto_k_max   <- max(loo_fit$diagnostics$pareto_k)
        pareto_k_share <- mean(loo_fit$diagnostics$pareto_k > 0.7)
        
        # Build block probabilities for diagnostics consistent with relabel
        P_emp <- block_P_empirical(A, z_hat, K)
        
        # ---- summarized diagnostics (no per-draw storage)
        diag_summ <- summarise_dcsbm_diagnostics(fit_relab, z_hat, K, n, m_items = m_item_triplets,A=A)
        
        K_hat_VI     = length(unique(z_hat))
        z_hat_alt    = salso::salso(fit_relab$z,loss = salso::VI(a = 0.5))
        K_hat_VI_a05 =  length(unique(z_hat_alt))
        
        
        
        
        diag_df <- as.data.frame(as.list(unclass(diag_summ)), check.names = FALSE)
        rows[[length(rows)+1]] <- cbind(
          data.frame(
            gen_model, fit_model, K,
            K_hat_VI=K_hat_VI,
            K_hat_VI_a05 = K_hat_VI_a05,
            difficulty = diff, rep,
            ari, vi, elpd_loo, p_loo, looic, pareto_k_max, pareto_k_share,
            mae_eta = NA, bias_eta = NA, cov_eta = NA,
            mae_psi = NA, bias_psi = NA, cov_psi = NA,
            mae_kappa = NA, bias_kappa = NA, cov_kappa = NA,
            ess_eta = NA, ac_eta = NA, ess_psi = NA, ac_psi = NA, ess_kappa = NA, ac_kappa = NA,
            time_sec = unname(tm["elapsed"]),
            stringsAsFactors = FALSE
          ),
          diag_df,
          stringsAsFactors = FALSE
        )
        
      } else {
        regime <- if (fit_model == "WST") "WST" else "SST"
        tm <- system.time({
          out <- fit_osbm(as.matrix(A), K, truth, regime, n_iter, burn, thin, verbose)
        })
        z_hat <- estimate_partition(out$z)
        ari   <- fossil::adj.rand.index(z_hat, z_t)
        vi    <- mcclust::vi.dist(z_hat, z_t)
        # --- Relabel BEFORE getting z_hat
        out_relab <- if(fit_model=='SST'){
          relabel_sst_with_reversal(mcmc_out = out)
        }else if(fit_model == "WST") {
          relabel_by_ECR(mcmc_out = out,model = 'WST')
        }
        # parameter summaries (unchanged)
        eta_hat <- colMeans(out_relab$eta); mae_eta <- mean(abs(eta_hat - eta_t)); bias_eta <- mean(eta_hat - eta_t)
        ci_eta  <- apply(out_relab$eta, 2, quantile, c(0.025, 0.975)); cov_eta <- mean(eta_t >= ci_eta[1,] & eta_t <= ci_eta[2,])
        
        kappa_mean <- apply(out_relab$kappa, c(2,3), mean)
        mae_kappa  <- mean(abs(kappa_mean - kappa_t)); bias_kappa <- mean(kappa_mean - kappa_t)
        ci_k_low   <- apply(out_relab$kappa, c(2,3), quantile, 0.025)
        ci_k_high  <- apply(out_relab$kappa, c(2,3), quantile, 0.975)
        cov_kappa  <- mean(kappa_t >= ci_k_low & kappa_t <= ci_k_high)
        
        if (gen_model == "WST" && fit_model == "WST") {
          psi_hat <- apply(out_relab$psi, c(2,3), mean)
          ut <- upper.tri(matrix(0, K, K))
          mae_psi <- mean(abs(psi_hat[ut] - truth$psi[ut]))
          bias_psi <- mean(psi_hat[ut] - truth$psi[ut])
          ci_psi_lo <- apply(out_relab$psi, c(2,3), quantile, 0.025)
          ci_psi_hi <- apply(out_relab$psi, c(2,3), quantile, 0.975)
          cov_psi <- mean(truth$psi[ut] >= ci_psi_lo[ut] & truth$psi[ut] <= ci_psi_hi[ut])
        } else if (gen_model == "SST" && fit_model == "SST") {
          psi_hat  <- colMeans(out_relab$psi)
          mae_psi  <- mean(abs(psi_hat - truth$psi))
          bias_psi <- mean(psi_hat - truth$psi)
          ci_psi   <- apply(out_relab$psi, 2, quantile, c(0.025, 0.975))
          cov_psi  <- mean(truth$psi >= ci_psi[1,] & truth$psi <= ci_psi[2,])
        } else {
          mae_psi <- bias_psi <- cov_psi <- NA_real_
        }
        
        ess_eta <- mean(coda::effectiveSize(coda::mcmc(out_relab$eta)))
        ac_eta  <- mean(coda::autocorr.diag(coda::mcmc(out_relab$eta))[,2])
        
        kappamat <- matrix(out_relab$kappa, ncol = K*K)
        ess_kappa <- mean(coda::effectiveSize(coda::mcmc(kappamat)))
        ac_kappa  <- mean(coda::autocorr.diag(coda::mcmc(kappamat))[,2])
        
        if (!is.null(out_relab$psi)) {
          if (is.matrix(out_relab$psi)) {
            ess_psi <- mean(coda::effectiveSize(coda::mcmc(out_relab$psi)))
            ac_psi  <- mean(coda::autocorr.diag(coda::mcmc(out_relab$psi))[,2])
          } else {
            psimat <- if (length(dim(out_relab$psi)) == 3)
              matrix(out_relab$psi, ncol = K*K) else as.matrix(out_relab$psi)
            ess_psi <- mean(coda::effectiveSize(coda::mcmc(psimat)))
            ac_psi  <- mean(coda::autocorr.diag(coda::mcmc(psimat))[,2])
          }
        } else ess_psi <- ac_psi <- NA_real_
        
        # LOO
        LL <- if (regime == "WST")
          loglik_matrix_wst(A = as.matrix(A), mcmc_out = out, dyad_index = DI)
        else
          loglik_matrix_sst(A = as.matrix(A), mcmc_out = out, dyad_index = DI)
        
        loo_fit <- loo::loo(LL)
        elpd_loo <- loo_fit$estimates["elpd_loo","Estimate"]
        p_loo    <- loo_fit$estimates["p_loo","Estimate"]
        looic    <- loo_fit$estimates["looic","Estimate"]
        pareto_k_max   <- max(loo_fit$diagnostics$pareto_k)
        pareto_k_share <- mean(loo_fit$diagnostics$pareto_k > 0.7)
        
        # ---- summarized diagnostics (OSBM)
        diag_summ <- summarise_osbm_diagnostics(out_relab, regime, K, z_hat, n, m_items = m_item_triplets,A=A)
        diag_df <- as.data.frame(as.list(unclass(diag_summ)), check.names = FALSE)
        
        K_hat_VI     = length(unique(z_hat))
        z_hat_alt    = salso::salso(out_relab$z,loss = salso::VI(a = 0.5))
        K_hat_VI_a05 =  length(unique(z_hat_alt))
        
        
        
        rows[[length(rows)+1]] <- cbind(data.frame(
          gen_model, fit_model, 
          K = K,
          K_hat_VI = K_hat_VI,
          K_hat_VI_a05 = K_hat_VI_a05,
          difficulty = diff, rep,
          ari, vi, elpd_loo, p_loo, looic, pareto_k_max, pareto_k_share,
          mae_eta, bias_eta, cov_eta,
          mae_psi, bias_psi, cov_psi,
          mae_kappa, bias_kappa, cov_kappa,
          ess_eta, ac_eta, ess_psi, ac_psi, ess_kappa, ac_kappa,
          # new diag summaries (means, CIs, ESS, MCSE)
          time_sec = unname(tm["elapsed"])
        ),
        diag_df,
        stringsAsFactors = FALSE
        )
      }
    } # for fit_model
    
    dplyr::bind_rows(rows)
    
  }, future.seed = TRUE)
  
  res_summary <- dplyr::bind_rows(job_results)
  
  # Save only the summary
  write.csv(res_summary, "sim_results_summary_ext.csv", row.names = FALSE)
  saveRDS(res_summary, "sim_results_summary_ext.rds")
  
  res_summary
}

# ---------------- Run ---------------------------------------------------

sim_study <- run_simulation_study()
