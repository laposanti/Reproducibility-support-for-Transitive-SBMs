## -----------------------------------------------------------
## 0. Libraries and data
## -----------------------------------------------------------
get_bundle_root <- function() {
  env_root <- Sys.getenv("TRANSITIVE_SBM_BUNDLE_ROOT", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
  }
  wd <- getwd()
  if (dir.exists(file.path(wd, "scripts")) && dir.exists(file.path(wd, "helper_folder"))) {
    return(normalizePath(wd, winslash = "/", mustWork = TRUE))
  }
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grepl("^--file=", cmd)]
  if (!length(file_arg)) {
    return(normalizePath(wd, winslash = "/", mustWork = TRUE))
  }
  script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg[1L])),
                               winslash = "/", mustWork = TRUE)
  normalizePath(file.path(dirname(script_path), "../.."),
                winslash = "/", mustWork = TRUE)
}

setwd(get_bundle_root())
source("scripts/bundle_defaults.R", local = TRUE)

library(dplyr)
library(tidyr)
library(kableExtra)
library(ggplot2)
library(glue)
sim_file <- Sys.getenv("SIM_RESULTS_PATH", unset = "")
if (!nzchar(sim_file)) {
  sim_file <- bundle_resolve_simulation_results_csv(must_exist = TRUE)
}

sim_raw <- read.csv(sim_file, stringsAsFactors = FALSE)

## Clean basic types (mirrors what you already do at the end)
sim <- sim_raw %>%
  select(-dplyr::starts_with("Unnamed")) %>%
  mutate(
    gen_model  = factor(gen_model,  levels = c("SST", "WST")),
    fit_model  = factor(fit_model,  levels = c("SST", "WST", "DC-SBM")),
    hierch     = factor(hierch,     levels = c("weak", "strong")),
    K_true     = as.integer(K_true),
    K_hat      = as.integer(K_hat),
    kappa_mean = as.numeric(kappa_mean),
    kappa_var  = as.numeric(kappa_var),
    psi_mean   = as.numeric(psi_mean),
    psi_var    = as.numeric(psi_var),
    seed       = as.integer(seed)
  )



## -----------------------------------------------------------
## 1. Partition / K recovery summary
## -----------------------------------------------------------

summary_partition <- sim %>%
  group_by(gen_model, hierch, kappa_mean, K_true, fit_model) %>%
  summarise(
    n_runs       = n(),
    mean_K_hat   = mean(K_hat, na.rm = TRUE),
    sd_K_hat     = sd(K_hat, na.rm = TRUE),
    prop_K_true  = mean(K_hat == K_true, na.rm = TRUE),
    prop_K_under = mean(K_hat < K_true, na.rm = TRUE),
    prop_K_over  = mean(K_hat > K_true, na.rm = TRUE),
    
    mean_ari     = mean(ari, na.rm = TRUE),
    sd_ari       = sd(ari, na.rm = TRUE),
    mean_vi      = mean(vi,  na.rm = TRUE),
    
    mean_elpd    = mean(elpd,  na.rm = TRUE),
    sd_elpd      = sd(elpd,    na.rm = TRUE),
    mean_pk_bad  = mean(pk_bad, na.rm = TRUE),
    .groups = "drop"
  )

## Optional: save a compact CSV for me / for the appendix
write.csv(summary_partition,
          file = "./output/simulation/tables/sim_summary_partition.csv",
          row.names = FALSE)

## Small formatter for mean + 5–95% range across replicates
fmt_ci <- function(m, lo, hi, digits = 2) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
          m, lo, hi)
}
## -----------------------------------------------------------
## 2. kable table for partition/K recovery
## -----------------------------------------------------------

make_partition_table <- function(gen, K_star, hier,
                                 caption_prefix = "Simulation:") {
  
  tab <- summary_partition %>%
    filter(gen_model == gen,
           K_true    == K_star,
           hierch    == hier) %>%
    arrange(kappa_mean, fit_model) %>%
    transmute(
      `kappa_mean` = kappa_mean,
      Model        = as.character(fit_model),
      `E[K_hat]`   = round(mean_K_hat, 2),
      `Pr(K_hat = K^*)` = round(prop_K_true, 2),
      `ARI(z_hat, z^*)` = round(mean_ari, 2),
      `ELPD`            = round(mean_elpd, 1)   # drop if you don't want ELPD here
    )
  
  cap <- paste0(
    caption_prefix, " generative ",
    if (gen == "SST") "SST" else "WST",
    ", $K^*=", K_star, "$, hierarchy ",
    as.character(hier),
    ": partition recovery (mean block count, exact $K$ recovery, ARI, ELPD)by fitted model and density ($\\kappa_{\\text{mean}}$)."
  )
  
  kbl(tab,
      format   = "latex",
      booktabs = TRUE,
      caption  = cap,
      digits   = 2,
      escape   = FALSE) %>%
    kable_styling(full_width = FALSE)
}

## Example uses (and saving):
tab_part_SST_K3_strong <- make_partition_table("SST", K_star = 3, hier = "strong")
save_kable(tab_part_SST_K3_strong,
           "./output/simulation/tables/appendix/tab_sim_partition_SST_K3_strong.tex")

tab_part_WST_K5_weak <- make_partition_table("WST", K_star = 5, hier = "weak")
save_kable(tab_part_WST_K5_weak,
           "./output/simulation/tables/appendix/tab_sim_partition_WST_K5_weak.tex")



## -----------------------------------------------------------
## 3. Empirical block-level WST/SST summary
## -----------------------------------------------------------

summary_order_block <- sim %>%
  group_by(gen_model, hierch, kappa_mean, K_true, fit_model) %>%
  summarise(
    # WST, block level
    zetaW_blk_mean = mean(thetaW_block_emp_mean, na.rm = TRUE),
    zetaW_blk_lo   = quantile(thetaW_block_emp_mean, probs = 0.05, na.rm = TRUE),
    zetaW_blk_hi   = quantile(thetaW_block_emp_mean, probs = 0.95, na.rm = TRUE),
    
    # SST, block level
    zetaS_blk_mean = mean(thetaS_block_emp_mean, na.rm = TRUE),
    zetaS_blk_lo   = quantile(thetaS_block_emp_mean, probs = 0.05, na.rm = TRUE),
    zetaS_blk_hi   = quantile(thetaS_block_emp_mean, probs = 0.95, na.rm = TRUE),
    
    # (optional) coverage of transitivity-triggering triples at block level
    cov_block_emp_mean = mean(coverage_block_emp_avg, na.rm = TRUE),
    
    .groups = "drop"
  )

## Optional CSV for me / appendices
write.csv(summary_order_block,
          file = "./output/simulation/tables/sim_summary_order_block.csv",
          row.names = FALSE)
## -----------------------------------------------------------
## 5. Plots: K recovery
## -----------------------------------------------------------

ggplot(summary_partition,
       aes(x = kappa_mean, y = mean_K_hat,
           color = fit_model,
           group = interaction(fit_model, K_true))) +
  geom_line() +
  geom_point(aes(shape = factor(K_true)), size = 2) +
  facet_grid(gen_model ~ hierch) +
  labs(
    x = expression(kappa[mean]),
    y = expression(E[K[hat]]),
    color = "Fitted model",
    shape = "K*",
    title = "Recovery of the number of blocks across generative scenarios"
  ) +
  theme_bw()

ggplot(summary_partition,
       aes(x = kappa_mean, y = prop_K_true,
           color = fit_model,
           group = interaction(fit_model, K_true))) +
  geom_line() +
  geom_point(aes(shape = factor(K_true)), size = 2) +
  facet_grid(gen_model ~ hierch) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = expression(kappa[mean]),
    y = expression(Pr( K[hat] == K)),
    color = "Fitted model",
    shape = "K*",
    title = "Probability of exact K recovery across scenarios"
  ) +
  theme_bw()



## Long format for WST/SST block-level conformity
order_block_long <- summary_order_block %>%
  select(gen_model, hierch, kappa_mean, K_true, fit_model,
         zetaW_blk_mean, zetaS_blk_mean) %>%
  pivot_longer(
    cols      = c(zetaW_blk_mean, zetaS_blk_mean),
    names_to  = "type",
    values_to = "zeta"
  ) %>%
  mutate(
    type = recode(type,
                  zetaW_blk_mean = "WST",
                  zetaS_blk_mean = "SST")
  )

ggplot(order_block_long,
       aes(x = kappa_mean, y = zeta,
           color = fit_model,
           group = interaction(fit_model, type, K_true))) +
  geom_line() +
  geom_point(aes(shape = factor(K_true)), size = 2) +
  facet_grid(gen_model ~ interaction(hierch, type)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = expression(kappa[mean]),
    y = expression(widehat(zeta)),
    color = "Fitted model",
    shape = "K*",
    title = "Empirical block-level WST/SST conformity across scenarios"
  ) +
  theme_bw()


sim %>%
  filter(gen_model == 'WST')%>%
  ggplot(aes(x = K_true, y = ari,group = interaction(fit_model,K_true), color = fit_model))+
  geom_boxplot()+
  facet_wrap(~interaction(hierch,kappa_mean))

sim %>%
  filter(gen_model == 'SST')%>%
  ggplot(aes(x = K_true, y = ari,group = interaction(fit_model,K_true), color = fit_model))+
  geom_boxplot()+
  facet_wrap(~interaction(hierch,kappa_mean))



sim %>%
  filter(gen_model == 'WST')%>%
  group_by(fit_model,K_true,hierch,kappa_mean)%>%
  summarise(mean_ELPD = mean(elpd),
            mean_vi = mean(vi),
            mean_ari = mean(ari))

sim %>%
  filter(gen_model == 'SST')%>%
  ggplot(aes(x = K_true, y = elpd,group = interaction(fit_model,K_true), color = fit_model))+
  geom_boxplot()+
  facet_wrap(~interaction(hierch,kappa_mean))



# -------------------------------------------------------------------
# 2.1 Performance summary (ARI, K_hat, ELPD, PSIS)
# -------------------------------------------------------------------
summary_perf <- sim %>%
  group_by(gen_model, hierch, K_true, kappa_mean, fit_model) %>%
  summarise(
    n_runs       = n(),
    mean_K_hat   = mean(K_hat, na.rm = TRUE),
    sd_K_hat     = sd(K_hat,   na.rm = TRUE),
    prop_K_true  = mean(K_hat == K_true, na.rm = TRUE),
    prop_K_under = mean(K_hat <  K_true, na.rm = TRUE),
    prop_K_over  = mean(K_hat >  K_true, na.rm = TRUE),
    
    mean_ari     = mean(ari,  na.rm = TRUE),
    sd_ari       = sd(ari,    na.rm = TRUE),
    
    mean_elpd    = mean(elpd, na.rm = TRUE),
    sd_elpd      = sd(elpd,   na.rm = TRUE),
    
    mean_pk_bad  = mean(pk_bad, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(gen_model, hierch, K_true, kappa_mean) %>%
  mutate(
    best_elpd  = max(mean_elpd, na.rm = TRUE),
    delta_elpd = mean_elpd - best_elpd
  ) %>%
  ungroup()
# -------------------------------------------------------------------
# 2.2 Proportion of times each model has best ELPD
# -------------------------------------------------------------------
elpd_best <- sim %>%
  group_by(gen_model, hierch, K_true, kappa_mean, rep_id) %>%
  mutate(
    elpd_max = max(elpd, na.rm = TRUE),
    is_best  = (elpd == elpd_max)
  ) %>%
  ungroup() %>%
  group_by(gen_model, hierch, K_true, kappa_mean, fit_model) %>%
  summarise(
    prop_best_elpd = mean(is_best, na.rm = TRUE),
    .groups = "drop"
  )

summary_perf <- summary_perf %>%
  left_join(elpd_best,
            by = c("gen_model", "hierch", "K_true", "kappa_mean", "fit_model")) %>%
  mutate(
    prop_best_elpd = ifelse(is.na(prop_best_elpd), 0, prop_best_elpd)
  )

# -------------------------------------------------------------------
# 2.3 Formatting helpers for tables
# -------------------------------------------------------------------
summary_perf <- summary_perf %>%
  mutate(
    K_hat_str   = sprintf("%.1f (%.1f)", mean_K_hat, sd_K_hat),
    ARI_str     = sprintf("%.2f (%.2f)", mean_ari,   sd_ari),
    ELPD_str    = sprintf("%.1f (%.1f)", mean_elpd,  sd_elpd),
    dELPD_str   = sprintf("%.1f",         delta_elpd),
    pk_str      = sprintf("%.2f",         mean_pk_bad),
    propK_str   = sprintf("%.2f",         prop_K_true),
    propBest_str = sprintf("%.2f",        prop_best_elpd)
  )


# -------------------------------------------------------------------
# 3. Empirical block-level WST/SST + violation rate
# -------------------------------------------------------------------
summary_order <- sim %>%
  group_by(gen_model, hierch, K_true, kappa_mean, fit_model) %>%
  summarise(
    thetaW_blk_mean = mean(thetaW_block_emp_mean, na.rm = TRUE),
    thetaW_blk_lo   = as.numeric(quantile(thetaW_block_emp_mean, 0.10, na.rm = TRUE)),
    thetaW_blk_hi   = as.numeric(quantile(thetaW_block_emp_mean, 0.90, na.rm = TRUE)),
    
    thetaS_blk_mean = mean(thetaS_block_emp_mean, na.rm = TRUE),
    thetaS_blk_lo   = as.numeric(quantile(thetaS_block_emp_mean, 0.10, na.rm = TRUE)),
    thetaS_blk_hi   = as.numeric(quantile(thetaS_block_emp_mean, 0.90, na.rm = TRUE)),
    
    viol_mean       = mean(violation_rate_mean, na.rm = TRUE),
    viol_lo         = as.numeric(quantile(violation_rate_mean, 0.10, na.rm = TRUE)),
    viol_hi         = as.numeric(quantile(violation_rate_mean, 0.90, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    zetaW_blk_str = sprintf("%.2f [%.2f, %.2f]",
                            thetaW_blk_mean, thetaW_blk_lo, thetaW_blk_hi),
    zetaS_blk_str = sprintf("%.2f [%.2f, %.2f]",
                            thetaS_blk_mean, thetaS_blk_lo, thetaS_blk_hi),
    zeta_viol_str = sprintf("%.3f [%.3f, %.3f]",
                            viol_mean, viol_lo, viol_hi)
  )
# -------------------------------------------------------------------
# 4.1 Partition / K recovery table
# -------------------------------------------------------------------
make_sim_table_partition <- function(gen = "SST",
                                     hier = "strong",
                                     K    = 3) {
  summary_perf %>%
    filter(gen_model == gen,
           hierch    == hier,
           K_true    == K) %>%
    arrange(kappa_mean, fit_model) %>%
    transmute(
      `kappa_mean` = kappa_mean,
      Model        = fit_model,
      `$\\hat K$ (mean sd)` = K_hat_str,
      `Pr($\\hat K = K^\\star$)` = propK_str,
      `ARI` = ARI_str
    ) %>%
    kable(
      format   = "latex",
      booktabs = TRUE,
      escape   = FALSE,
      caption  = glue(
        "{gen} generative model, hierarchy {hier}, $K^\\star={K}$: ",
        "partition recovery and $K$ recovery across densities."
      )
    ) %>%
    kable_styling(full_width = FALSE, position = "center")
}
tab_part_SST_strong_K3 <- make_sim_table_partition("SST", "strong", 3)
save_kable(tab_part_SST_strong_K3,
           "./output/simulation/tables/appendix/tab_sim_partition_SST_strong_K3.tex")

# -------------------------------------------------------------------
# 4.2 ELPD-based model selection table
# -------------------------------------------------------------------
make_sim_table_elpd <- function(gen = "SST",
                                hier = "strong",
                                K    = 3) {
  summary_perf %>%
    filter(gen_model == gen,
           hierch    == hier,
           K_true    == K) %>%
    arrange(kappa_mean, fit_model) %>%
    transmute(
      `kappa_mean` = kappa_mean,
      Model        = fit_model,
      `ELPD`       = ELPD_str,
      `$\\Delta$ELPD` = dELPD_str,
      `PSIS tail`  = pk_str,
      `Pr(best ELPD)` = propBest_str
    ) %>%
    kable(
      format   = "latex",
      booktabs = TRUE,
      escape   = FALSE,
      caption  = glue(
        "{gen} generative model, hierarchy {hier}, $K^\\star={K}$: ",
        "ELPD comparison and PSIS diagnostics across densities."
      )
    ) %>%
    kable_styling(full_width = FALSE, position = "center")
}

tab_elpd_SST_strong_K3 <- make_sim_table_elpd("SST", "strong", 3)
save_kable(tab_elpd_SST_strong_K3,
           "./output/simulation/tables/appendix/tab_sim_elpd_SST_strong_K3.tex")



# -------------------------------------------------------------------
# 4.3 Block-level conformity + violations table
# -------------------------------------------------------------------
make_sim_table_order <- function(gen = "SST",
                                 hier = "strong",
                                 K    = 3) {
  summary_order %>%
    filter(gen_model == gen,
           hierch    == hier,
           K_true    == K) %>%
    arrange(kappa_mean, fit_model) %>%
    transmute(
      `kappa_mean` = kappa_mean,
      Model        = fit_model,
      `$\\widehat\\zeta_{\\mathrm{W}}^{\\mathrm{blk}}$` = zetaW_blk_str,
      `$\\widehat\\zeta_{\\mathrm{S}}^{\\mathrm{blk}}$` = zetaS_blk_str,
      `$\\zeta^{\\mathrm{viol}}_{\\hat z}$`             = zeta_viol_str
    ) %>%
    kable(
      format   = "latex",
      booktabs = TRUE,
      escape   = FALSE,
      caption  = glue(
        "{gen} generative model, hierarchy {hier}, $K^\\star={K}$: ",
        "empirical block-level WST/SST conformity and cross-block violation rates."
      )
    ) %>%
    kable_styling(full_width = FALSE, position = "center")
}
tab_order_SST_strong_K3 <- make_sim_table_order("SST", "strong", 3)
save_kable(tab_order_SST_strong_K3,
           "./output/simulation/tables/appendix/tab_sim_order_SST_strong_K3.tex")



for (gen in c("SST", "WST")) {
  for (hier in c("weak", "strong")) {
    for (K in c(3, 5)) {
      tab_part  <- make_sim_table_partition(gen, hier, K)
      tab_elpd  <- make_sim_table_elpd(gen, hier, K)
      tab_order <- make_sim_table_order(gen, hier, K)
      
      save_kable(tab_part,
                 glue("./output/simulation/tables/appendix/tab_sim_partition_{gen}_{hier}_K{K}.tex"))
      save_kable(tab_elpd,
                 glue("./output/simulation/tables/appendix/tab_sim_elpd_{gen}_{hier}_K{K}.tex"))
      save_kable(tab_order,
                 glue("./output/simulation/tables/appendix/tab_sim_order_{gen}_{hier}_K{K}.tex"))
    }
  }
}



# -------------------------------------------------------------------
# 5.1 ARI by K_true, facetted by (hierch, kappa_mean)
# -------------------------------------------------------------------

sst_ari_rows <- sim %>% filter(gen_model == "SST")
if (nrow(sst_ari_rows) > 0L) {
  Ari_boxplot <- ggplot(
    sst_ari_rows,
    aes(
      x = factor(K_true),
      y = ari,
      group = fit_model,
      colour = fit_model
    )
  ) +
    geom_boxplot(outlier.alpha = 0.3) +
    facet_wrap(~ interaction(hierch, kappa_mean)) +
    labs(
      x = expression(K^"*"),
      y = "ARI",
      colour = "Fitted model",
      title = "Partition recovery (ARI) under SST generative model"
    ) +
    theme_bw()
  ggsave("./output/simulation/plots/ari_sst.png", Ari_boxplot)
} else {
  message("Skipping ari_sst.png because the selected simulation results contain no SST-generated rows.")
}







#----


## -------------------------------------------------------------------
## 4. Compact paper tables (three tables only)
##    - Table 1: partition + K recovery
##    - Table 2: ELPD / model choice
##    - Table 3: order satisfaction (WST / SST)
##    Only lowest and highest kappa_mean are kept.
## -------------------------------------------------------------------

# pick only the lowest and highest kappa
kappa_vals <- sort(unique(summary_perf$kappa_mean))
kappa_keep <- c(min(kappa_vals), max(kappa_vals))

## --------------------------- 4.1 Partition / K recovery -------------
tab_partition <- summary_perf %>%
  filter(kappa_mean %in% kappa_keep) %>%
  arrange(gen_model, hierch, K_true, kappa_mean, fit_model) %>%
  transmute(
    `Generative` = gen_model,
    `Hierarchy`  = hierch,
    `K^*`        = K_true,
    `kappa_mean` = kappa_mean,
    `Model`      = fit_model,
    `E[K_hat] (sd)`        = K_hat_str,
    `Pr(K_hat = K^*)`      = propK_str,
    `ARI (mean sd)`        = ARI_str
  )

tab_partition_kbl <-
  kbl(tab_partition,
      format   = "latex",
      booktabs = TRUE,
      escape   = FALSE,
      caption  = "Partition recovery and K recovery by generative model, hierarchy, true K, density (kappa_mean) and fitted model (only lowest and highest densities).") %>%
  kable_styling(full_width = FALSE, position = "center")

save_kable(tab_partition_kbl,
           "./output/simulation/tables/tab_sim_partition_main.tex")

## --------------------------- 4.2 ELPD / model choice ----------------
tab_elpd <- summary_perf %>%
  filter(kappa_mean %in% kappa_keep) %>%
  arrange(gen_model, hierch, K_true, kappa_mean, fit_model) %>%
  transmute(
    `Generative` = gen_model,
    `Hierarchy`  = hierch,
    `K^*`        = K_true,
    `kappa_mean` = kappa_mean,
    `Model`      = fit_model,
    `ELPD (mean sd)` = ELPD_str,
    `ΔELPD`          = dELPD_str,
    `Pr(best ELPD)`  = propBest_str
  )

tab_elpd_kbl <-
  kbl(tab_elpd,
      format   = "latex",
      booktabs = TRUE,
      escape   = FALSE,
      caption  = "Predictive performance: ELPD, ΔELPD and probability of having the best ELPD by generative model, hierarchy, true K, density and fitted model (only lowest and highest densities).") %>%
  kable_styling(full_width = FALSE, position = "center")

save_kable(tab_elpd_kbl,
           "./output/simulation/tables/tab_sim_elpd_main.tex")

## --------------------------- 4.3 Order satisfaction -----------------
# use the summary_order already computed above
kappa_vals_order <- sort(unique(summary_order$kappa_mean))
kappa_keep_order <- c(min(kappa_vals_order), max(kappa_vals_order))

tab_order <- summary_order %>%
  filter(kappa_mean %in% kappa_keep_order) %>%
  arrange(gen_model, hierch, K_true, kappa_mean, fit_model) %>%
  transmute(
    `Generative` = gen_model,
    `Hierarchy`  = hierch,
    `K^*`        = K_true,
    `kappa_mean` = kappa_mean,
    `Model`      = fit_model,
    `$\\widehat{\\zeta}_W^{blk}$` = zetaW_blk_str,
    `$\\widehat{\\zeta}_S^{blk}$` = zetaS_blk_str,
    `$\\zeta^{viol}_{\\hat z}$`   = zeta_viol_str
  )

tab_order_kbl <-
  kbl(tab_order,
      format   = "latex",
      booktabs = TRUE,
      escape   = FALSE,
      caption  = "Empirical block-level WST/SST conformity and violation rate by generative model, hierarchy, true K, density and fitted model (only lowest and highest densities).") %>%
  kable_styling(full_width = FALSE, position = "center")

save_kable(tab_order_kbl,
           "./output/simulation/tables/tab_sim_order_main.tex")
