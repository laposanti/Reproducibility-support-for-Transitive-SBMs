############################################################
## Mixing-moves comparison simulation study
##
## Generates data from known ground truth (SST & WST),
## runs the sampler WITH and WITHOUT mixing moves,
## compares posterior K, ARI, VI, and trace diagnostics.
############################################################

library(stats)
library(salso)
library(BayesLogit)
library(truncnorm)
library(fossil)
library(mcclust)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Matrix)
library(Rcpp)

source("./helper_folder/sim_study_helper.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./helper_folder/Hyper_setup.R")
source("./core/my_best_try_so_far.R")

sgn      <- function(x) ifelse(x > 0, 1L, ifelse(x < 0, -1L, 0L))
logistic <- function(x) 1 / (1 + exp(-x))

## ---- Log progress to file ------------------------------------------------
log_file <- "output/simulation/plots/mixing_moves_sim_log.txt"
dir.create("output/simulation/plots", recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)   # tee to console + file
sink(log_con, type = "message", append = TRUE)
cat(sprintf("=== Mixing-moves sim study started: %s ===\n", Sys.time()))

## ---- Data generation helpers (from DemoKvar.R) ---------------------------
sample_kappa_prior <- function(K, mean, var) {
  shape <- mean^2 / var; rate <- mean / var
  KAP <- matrix(0, K, K)
  for (k in 1:K) for (l in k:K) KAP[k, l] <- rgamma(1, shape = shape, rate = rate)
  KAP[lower.tri(KAP)] <- t(KAP)[lower.tri(KAP)]
  KAP
}

sample_psi_wst_prior <- function(K, mu_psi, var_psi) {
  sigma_psi <- sqrt(var_psi)
  PSI <- matrix(0, K, K)
  PSI[upper.tri(PSI)] <- truncnorm::rtruncnorm(
    n = K * (K - 1) / 2, a = 0, mean = mu_psi, sd = sigma_psi
  )
  PSI[lower.tri(PSI)] <- -t(PSI)[lower.tri(PSI)]
  PSI
}

sample_psi_sst_prior <- function(K, mu_psi, var_psi) {
  if (K <= 1) return(numeric(0))
  mu_delta  <- mu_psi  / (K - 1)
  tau_delta <- sqrt(var_psi) / (K - 1)
  deltas <- truncnorm::rtruncnorm(n = K - 1, a = 0, mean = mu_delta, sd = tau_delta)
  cumsum(deltas)
}

############################################################
## 1. Simulation settings
############################################################

# Scenarios: one SST and one WST, each with moderate signal
scenarios <- list(
  list(gen_model = "SST", K_true = 3, n = 40,
       kappa_mean = 2.0, kappa_var = 1.44,
       psi_mean = 1.0, psi_var = 0.09,
       hierch = "strong"),
  list(gen_model = "WST", K_true = 3, n = 40,
       kappa_mean = 2.0, kappa_var = 1.44,
       psi_mean = 1.0, psi_var = 0.09,
       hierch = "strong"),
  list(gen_model = "SST", K_true = 5, n = 60,
       kappa_mean = 2.0, kappa_var = 1.44,
       psi_mean = 0.70, psi_var = 0.09,
       hierch = "medium"),
  list(gen_model = "WST", K_true = 5, n = 60,
       kappa_mean = 2.0, kappa_var = 1.44,
       psi_mean = 0.70, psi_var = 0.09,
       hierch = "medium")
)

n_rep   <- 3L
n_iter  <- 6000L
burn    <- 1500L
thin    <- 2L
base_seed <- 42L

############################################################
## 2. Run all scenarios
############################################################

all_results <- list()
all_K_traces <- list()

for (sc_idx in seq_along(scenarios)) {
  sc <- scenarios[[sc_idx]]
  cat(sprintf("\n========== Scenario %d: %s K=%d n=%d %s ==========\n",
              sc_idx, sc$gen_model, sc$K_true, sc$n, sc$hierch))

  for (rep_id in seq_len(n_rep)) {
    seed <- base_seed + (sc_idx - 1) * 100 + rep_id
    set.seed(seed)

    # Generate data
    z_true     <- rep(seq_len(sc$K_true), length.out = sc$n)
    kappa_true <- sample_kappa_prior(sc$K_true, mean = sc$kappa_mean, var = sc$kappa_var)
    eta_true   <- runif(sc$n, 0.9, 1.1)

    if (sc$gen_model == "WST") {
      psi_true <- sample_psi_wst_prior(sc$K_true, mu_psi = sc$psi_mean, var_psi = sc$psi_var)
      A <- as.matrix(simulate_osbm(sc$n, sc$K_true, z_true, eta_true, kappa_true, psi_true, "WST"))
    } else {
      psi_true <- sample_psi_sst_prior(sc$K_true, mu_psi = sc$psi_mean, var_psi = sc$psi_var)
      A <- as.matrix(simulate_osbm(sc$n, sc$K_true, z_true, eta_true, kappa_true, psi_true, "SST"))
    }

    hypers <- get_corollary_calibrated_hypers(
      A = A, K_expected = 1,
      ordering_prior_mode = "equivalence_class",
      a_kappa = 2, a_eta = 2, mu0 = 1.0,
      gamma_bounds = c(0.3, 0.7)
    )

    # Run both configurations
    for (use_moves in c(FALSE, TRUE)) {
      label <- if (use_moves) "with_moves" else "baseline"
      cat(sprintf("  Rep %d, %s ... ", rep_id, label))

      set.seed(seed + if (use_moves) 500L else 0L)
      t0 <- proc.time()

      out <- tryCatch(
        modular_osbm_sampler(
          A = A, K = nrow(A), n_iter = n_iter, burn = burn, thin = thin,
          verbose = FALSE,
          psi_constraint = sc$gen_model ,
          partition_prior = "GN",
          gamma_gn        = hypers$gamma,
          a_kappa = hypers$a_kappa, b_kappa = hypers$b_kappa,
          a_eta   = hypers$a_eta,   b_eta   = hypers$b_eta,
          mu0 = hypers$mu0, sigma0 = hypers$sigma0, tau0 = hypers$tau0,
          use_mixing_moves = use_moves,
          seed = seed + if (use_moves) 500L else 0L
        ),
        error = function(e) { message("ERROR: ", e$message); NULL }
      )

      elapsed <- (proc.time() - t0)["elapsed"]

      if (is.null(out)) {
        cat("FAILED\n")
        next
      }

      # Extract K trace and z chain
      z_chain <- do.call(rbind, lapply(out$z, as.integer))
      K_count <- function(x) length(unique(x))
      K_trace <- apply(z_chain, 1, K_count)
      S <- nrow(z_chain)

      # Point estimate of partition
      z_hat <- salso::salso(z_chain, loss = salso::VI(), nRuns = 1L)
      ari   <- fossil::adj.rand.index(z_hat, z_true)
      vi    <- mcclust::vi.dist(z_hat, z_true)
      K_hat <- length(unique(z_hat))
      K_post_mean <- mean(K_trace)
      K_post_mode <- as.integer(names(sort(table(K_trace), decreasing = TRUE))[1])

      cat(sprintf("K_hat=%d ARI=%.3f VI=%.3f (%.1fs)\n",
                  K_hat, ari, vi, elapsed))

      row <- data.frame(
        scenario    = sc_idx,
        gen_model   = sc$gen_model,
        K_true      = sc$K_true,
        n           = sc$n,
        hierch      = sc$hierch,
        rep_id      = rep_id,
        method      = label,
        K_hat       = K_hat,
        K_post_mean = K_post_mean,
        K_post_mode = K_post_mode,
        ari         = ari,
        vi          = vi,
        time        = elapsed,
        stringsAsFactors = FALSE
      )
      all_results[[length(all_results) + 1L]] <- row

      # Store K trace for plotting
      all_K_traces[[length(all_K_traces) + 1L]] <- data.frame(
        scenario  = sc_idx,
        gen_model = sc$gen_model,
        K_true    = sc$K_true,
        hierch    = sc$hierch,
        rep_id    = rep_id,
        method    = label,
        iter      = seq_len(S),
        K         = K_trace,
        stringsAsFactors = FALSE
      )
    }
  }
}

############################################################
## 3. Aggregate results
############################################################

results_df <- dplyr::bind_rows(all_results)
K_traces_df <- dplyr::bind_rows(all_K_traces)

cat("\n\n========== SUMMARY TABLE ==========\n")
summary_table <- results_df %>%
  group_by(gen_model, K_true, hierch, method) %>%
  summarise(
    mean_ARI     = mean(ari, na.rm = TRUE),
    sd_ARI       = sd(ari, na.rm = TRUE),
    mean_VI      = mean(vi, na.rm = TRUE),
    sd_VI        = sd(vi, na.rm = TRUE),
    mean_K_hat   = mean(K_hat, na.rm = TRUE),
    mean_K_pmean = mean(K_post_mean, na.rm = TRUE),
    K_mode_correct = mean(K_post_mode == K_true, na.rm = TRUE),
    mean_time    = mean(time, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(summary_table))

# Save
dir.create("output/simulation/plots", recursive = TRUE, showWarnings = FALSE)
write.csv(results_df, "output/simulation/plots/mixing_moves_comparison_results.csv", row.names = FALSE)

############################################################
## 4. Plots
############################################################

# --- 4a. Posterior K distribution (bar plots) ---
K_posterior_df <- K_traces_df %>%
  group_by(scenario, gen_model, K_true, hierch, method, K) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(scenario, gen_model, K_true, hierch, method) %>%
  mutate(prob = count / sum(count)) %>%
  ungroup()

p_K_post <- ggplot(K_posterior_df,
                   aes(x = factor(K), y = prob, fill = method)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_vline(aes(xintercept = as.numeric(factor(K_true,
                  levels = levels(factor(K_posterior_df$K))))),
             linetype = "dashed", color = "red", linewidth = 0.6) +
  facet_wrap(~ gen_model + paste0("K=", K_true) + hierch,
             scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("baseline" = "steelblue", "with_moves" = "coral")) +
  labs(x = "K", y = "Posterior probability",
       title = "Posterior distribution of K: baseline vs mixing moves",
       fill = "Method") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("output/simulation/plots/mixing_moves_K_posterior.pdf", p_K_post, width = 10, height = 8)
ggsave("output/simulation/plots/mixing_moves_K_posterior.png", p_K_post, width = 10, height = 8, dpi = 150)
cat("Saved: output/simulation/plots/mixing_moves_K_posterior.pdf\n")

# --- 4b. K trace plots ---
p_K_trace <- ggplot(K_traces_df %>% filter(rep_id == 1),
                    aes(x = iter, y = K, color = method)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  geom_hline(aes(yintercept = K_true), linetype = "dashed", color = "black") +
  facet_wrap(~ gen_model + paste0("K=", K_true) + hierch,
             scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("baseline" = "steelblue", "with_moves" = "coral")) +
  labs(x = "Saved iteration", y = "K",
       title = "Trace of K: baseline vs mixing moves (rep 1)",
       color = "Method") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("output/simulation/plots/mixing_moves_K_trace.pdf", p_K_trace, width = 10, height = 8)
ggsave("output/simulation/plots/mixing_moves_K_trace.png", p_K_trace, width = 10, height = 8, dpi = 150)
cat("Saved: output/simulation/plots/mixing_moves_K_trace.pdf\n")

# --- 4c. ARI / VI comparison ---
metrics_long <- results_df %>%
  select(gen_model, K_true, hierch, rep_id, method, ari, vi) %>%
  pivot_longer(cols = c(ari, vi), names_to = "metric", values_to = "value") %>%
  mutate(metric = toupper(metric))

p_metrics <- ggplot(metrics_long,
                    aes(x = method, y = value, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  facet_wrap(~ gen_model + paste0("K=", K_true) + hierch + metric,
             scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("baseline" = "steelblue", "with_moves" = "coral")) +
  labs(x = "", y = "Value",
       title = "Partition recovery: ARI and VI distance",
       fill = "Method") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_blank())

ggsave("output/simulation/plots/mixing_moves_metrics.pdf", p_metrics, width = 12, height = 8)
ggsave("output/simulation/plots/mixing_moves_metrics.png", p_metrics, width = 12, height = 8, dpi = 150)
cat("Saved: output/simulation/plots/mixing_moves_metrics.pdf\n")

cat(sprintf("\n========== Done: %s ==========\n", Sys.time()))

## Close log sinks
sink(type = "message")
sink(type = "output")
close(log_con)
