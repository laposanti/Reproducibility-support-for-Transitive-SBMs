#!/usr/bin/env Rscript

get_repo_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grepl("^--file=", cmd)]
  if (!length(file_arg)) return(getwd())
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = TRUE)
  normalizePath(file.path(dirname(script_path), "../.."), winslash = "/", mustWork = TRUE)
}

setwd(get_repo_root())

suppressPackageStartupMessages({
  library(coda)
  library(parallel)
  library(salso)
})

source("./helper_folder/sim_study_helper.R")
source("./helper_folder/helper.R")
source("./helper_folder/SST_helpers.R")
source("./helper_folder/WST_helpers.R")
source("./helper_folder/Hyper_setup.R")
source("./helper_folder/transitivity_check_helper.R")
source("./core/my_best_try_so_far.R")
source("./core/DCSBM_varK.R")

run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- Sys.getenv(
  "SUPP_STUDY_OUT_DIR",
  unset = file.path("output", "supplement_sensitivity", run_tag)
)
paper_table_dir <- file.path("output", "paper", "tables", "current")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paper_table_dir, recursive = TRUE, showWarnings = FALSE)

n_chains <- as.integer(Sys.getenv("SUPP_STUDY_CHAINS", unset = "3"))
n_iter <- as.integer(Sys.getenv("SUPP_STUDY_N_ITER", unset = "1200"))
burn <- as.integer(Sys.getenv("SUPP_STUDY_BURN", unset = "600"))
thin <- as.integer(Sys.getenv("SUPP_STUDY_THIN", unset = "3"))
mc_cores <- as.integer(Sys.getenv("SUPP_STUDY_CORES", unset = as.character(min(n_chains, 3L))))
seed_base <- as.integer(Sys.getenv("SUPP_STUDY_SEED", unset = "20260607"))

cat("=== Supplement convergence + sensitivity study ===\n")
cat("out_dir:", out_dir, "\n")
cat(sprintf("chains=%d, iter=%d, burn=%d, thin=%d, cores=%d\n\n",
            n_chains, n_iter, burn, thin, mc_cores))

setting_grid <- data.frame(
  setting_id = c("weak_reg", "default", "strong_reg"),
  setting_label = c("Weak prior", "Default prior", "Strong prior"),
  theta_ocrp = c(1.0, 0.5, 0.2),
  alpha_crp = c(1.0, 0.5, 0.2),
  b_scale = c(0.5, 1.0, 2.0),
  dir_scale = c(1.5, 1.0, 0.67),
  stringsAsFactors = FALSE
)
settings_filter <- trimws(Sys.getenv("SUPP_STUDY_SETTINGS", unset = ""))
if (nzchar(settings_filter)) {
  keep <- trimws(strsplit(settings_filter, ",", fixed = TRUE)[[1L]])
  setting_grid <- setting_grid[setting_grid$setting_id %in% keep, , drop = FALSE]
  if (!nrow(setting_grid)) {
    stop("SUPP_STUDY_SETTINGS filtered out every setting.", call. = FALSE)
  }
}

model_levels <- c("WST", "SST", "DCSBM")
model_labels <- c(WST = "WST--OSBM", SST = "Toeplitz SST--OSBM", DCSBM = "DC--SBM")

dataset_labels <- c(sim_sst = "Simulated SST dataset", moreno_sheep = "Bighorn sheep")

`%||%` <- function(x, y) if (is.null(x)) y else x

choose_dataset <- function(dataset = c("moreno_sheep")) {
  dataset <- match.arg(dataset)
  if (dataset == "moreno_sheep") {
    edges <- read.csv("./data/moreno_sheep/edges.csv", comment.char = "#", strip.white = TRUE)
    names(edges)[1:3] <- c("source", "target", "weight")
    edges$source <- as.integer(edges$source)
    edges$target <- as.integer(edges$target)
    edges$weight <- as.integer(edges$weight)
    if (min(edges$source, edges$target) == 0L) {
      edges$source <- edges$source + 1L
      edges$target <- edges$target + 1L
    }
    n_nodes <- max(c(edges$source, edges$target))
    A <- matrix(0L, n_nodes, n_nodes)
    for (r in seq_len(nrow(edges))) {
      i <- edges$source[r]
      j <- edges$target[r]
      w <- edges$weight[r]
      if (w > 0L) A[i, j] <- A[i, j] + w
    }
    diag(A) <- 0L
    return(A)
  }
  stop("Unknown dataset: ", dataset)
}

build_simulated_dataset <- function(seed = 99L) {
  set.seed(seed)
  n <- 60L
  K_true <- 5L
  z_true <- rep(seq_len(K_true), each = n / K_true)
  eta_true <- rep(1, n)
  kappa_true <- matrix(1.5, K_true, K_true)
  diag(kappa_true) <- 3.0
  psi_true <- (seq_len(K_true - 1L) * 1.3) / (K_true - 1L)
  A <- as.matrix(simulate_osbm(n, K_true, z_true, eta_true, kappa_true, psi_true, regime = "SST"))
  list(
    dataset_id = "sim_sst",
    dataset_label = dataset_labels[["sim_sst"]],
    A = A,
    truth = list(K_true = K_true, z_true = z_true, regime = "SST", psi_true = psi_true, kappa_true = kappa_true)
  )
}

build_real_dataset <- function() {
  list(
    dataset_id = "moreno_sheep",
    dataset_label = dataset_labels[["moreno_sheep"]],
    A = choose_dataset("moreno_sheep"),
    truth = NULL
  )
}

vi_between <- function(z1, z2) {
  stopifnot(length(z1) == length(z2))
  n <- length(z1)
  tab <- table(z1, z2)
  p <- tab / n
  p1 <- rowSums(p)
  p2 <- colSums(p)
  h1 <- -sum(p1[p1 > 0] * log(p1[p1 > 0]))
  h2 <- -sum(p2[p2 > 0] * log(p2[p2 > 0]))
  mi <- sum(p[p > 0] * (log(p[p > 0]) - log(outer(p1, p2)[p > 0])))
  h1 + h2 - 2 * mi
}

z_to_matrix_local <- function(z_obj) {
  if (is.matrix(z_obj)) return(z_obj)
  do.call(rbind, lapply(z_obj, as.integer))
}

extract_K_trace_local <- function(fit) {
  as.integer((fit$K_trace %||% fit$K))
}

summarise_K_trace_local <- function(K_trace) {
  K_trace <- K_trace[is.finite(K_trace)]
  tab <- table(K_trace)
  list(
    mode = as.integer(names(which.max(tab))),
    mean = mean(K_trace),
    lo = as.integer(stats::quantile(K_trace, 0.025, names = FALSE, type = 1)),
    hi = as.integer(stats::quantile(K_trace, 0.975, names = FALSE, type = 1))
  )
}

clamp_prob <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)

order_partition_by_block_score <- function(A, z_hat) {
  K <- length(unique(z_hat))
  N <- A + t(A)
  node_score <- numeric(length(z_hat))
  for (i in seq_along(z_hat)) {
    observed <- which(N[i, ] > 0 & seq_along(z_hat) != i)
    if (!length(observed)) {
      node_score[i] <- 0
    } else {
      node_score[i] <- mean(A[i, observed] / N[i, observed])
    }
  }
  block_score <- tapply(node_score, z_hat, mean)
  ord <- order(-block_score, seq_along(block_score))
  map <- integer(length(block_score))
  map[ord] <- seq_along(ord)
  unname(map[z_hat])
}

compute_violation_rate <- function(A, z_hat, reorder = FALSE) {
  z_use <- if (reorder) order_partition_by_block_score(A, z_hat) else z_hat
  idx <- which(A > 0, arr.ind = TRUE)
  if (!nrow(idx)) return(NA_real_)
  from_rank <- z_use[idx[, 1]]
  to_rank <- z_use[idx[, 2]]
  w <- A[idx]
  cross <- from_rank != to_rank
  cross_mass <- sum(w[cross])
  if (cross_mass <= 0) return(NA_real_)
  backward <- cross & from_rank > to_rank
  sum(w[backward]) / cross_mass
}

ordered_rho_mean <- function(lambda_mat) {
  rho <- lambda_mat / (lambda_mat + t(lambda_mat))
  rho <- clamp_prob(rho)
  diag(rho) <- 0.5
  score <- rowMeans(rho)
  rho <- rho[order(-score, seq_along(score)), order(-score, seq_along(score)), drop = FALSE]
  mean(rho[upper.tri(rho)])
}

osbm_scalar_traces <- function(fit, model) {
  data.frame(
    K = extract_K_trace_local(fit),
    eta_mean = vapply(fit$eta, mean, numeric(1)),
    volume_mean = vapply(fit$kappa, function(x) mean(x[upper.tri(x, diag = TRUE)]), numeric(1)),
    direction_mean = if (model == "WST") {
      vapply(fit$psi, function(x) {
        vals <- x[upper.tri(x)]
        mean(plogis(vals))
      }, numeric(1))
    } else {
      vapply(fit$psi, function(x) {
        if (!length(x)) return(0.5)
        mean(plogis(x))
      }, numeric(1))
    },
    stringsAsFactors = FALSE
  )
}

dcsbm_scalar_traces <- function(fit) {
  data.frame(
    K = as.integer(fit$K),
    eta_mean = rowMeans(fit$theta),
    volume_mean = vapply(fit$lambda, mean, numeric(1)),
    direction_mean = vapply(fit$lambda, ordered_rho_mean, numeric(1)),
    stringsAsFactors = FALSE
  )
}

safe_psrf <- function(chain_list) {
  chain_list <- lapply(chain_list, function(x) as.numeric(x[is.finite(x)]))
  min_len <- min(vapply(chain_list, length, integer(1)))
  if (!is.finite(min_len) || min_len < 20L) return(c(point = NA_real_, upper = NA_real_))
  chain_list <- lapply(chain_list, function(x) x[seq_len(min_len)])
  sds <- vapply(chain_list, stats::sd, numeric(1))
  if (all(sds == 0)) {
    vals <- vapply(chain_list, function(x) x[1], numeric(1))
    if (length(unique(round(vals, 10))) == 1L) return(c(point = 1, upper = 1))
    return(c(point = NA_real_, upper = NA_real_))
  }
  chain_list <- Map(function(x, sdx) {
    if (sdx == 0) x + stats::rnorm(length(x), sd = 1e-8) else x
  }, chain_list, sds)
  out <- tryCatch(
    coda::gelman.diag(coda::mcmc.list(lapply(chain_list, coda::mcmc)), autoburnin = FALSE, multivariate = FALSE)$psrf[1, ],
    error = function(e) c(NA_real_, NA_real_)
  )
  c(point = unname(out[1]), upper = unname(out[2]))
}

convergence_summary <- function(fit_list, model) {
  trace_list <- if (model == "DCSBM") {
    lapply(fit_list, dcsbm_scalar_traces)
  } else {
    lapply(fit_list, function(f) osbm_scalar_traces(f, model))
  }
  vars <- c("K", "eta_mean", "volume_mean", "direction_mean")
  stats <- lapply(vars, function(v) safe_psrf(lapply(trace_list, `[[`, v)))
  names(stats) <- vars
  data.frame(
    rhat_K = stats$K[["point"]],
    rhat_eta = stats$eta_mean[["point"]],
    rhat_volume = stats$volume_mean[["point"]],
    rhat_direction = stats$direction_mean[["point"]],
    rhat_upper_max = max(vapply(stats, function(x) x[["upper"]], numeric(1)), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

fit_one_chain <- function(A, dataset_id, model, setting_row, seed, K_expected) {
  cat(sprintf("[fit] dataset=%s setting=%s model=%s seed=%d\n",
              dataset_id, setting_row$setting_id, model, seed))
  if (model %in% c("WST", "SST")) {
    hypers <- get_principled_hypers_v2(A = A, K_expected = K_expected, c_kappa = 3)
    out <- modular_osbm_sampler(
      A = A,
      K = nrow(A),
      n_iter = n_iter,
      burn = burn,
      thin = thin,
      verbose = FALSE,
      psi_constraint = model,
      partition_prior = "OCRP",
      theta_ocrp = setting_row$theta_ocrp,
      eta_identifiability = "block_sum",
      a_kappa = 1,
      b_kappa = setting_row$b_scale,
      a_eta = 1,
      b_eta = 1,
      mu0 = hypers$mu0,
      sigma0 = hypers$sigma0 * setting_row$dir_scale,
      tau0 = max(hypers$tau0, 0.2) * setting_row$dir_scale,
      use_mixing_moves = TRUE,
      seed = seed
    )
    return(out)
  }
  fit_dcsbm_gibbs_gnedin(
    A = as.matrix(A),
    K_init = nrow(A),
    priors = list(
      a_eta = 1,
      b_eta = 1,
      a_lambda = 1,
      b_lambda = setting_row$b_scale,
      alpha_crp = setting_row$alpha_crp,
      partition_prior = "CRP"
    ),
    iters = n_iter,
    burn_in = burn,
    thin = thin,
    verbose = 0,
    seed = seed
  )
}

run_group <- function(dataset_obj, setting_row, model) {
  chain_seeds <- seed_base + seq_len(n_chains) + 1000L * match(model, model_levels) + 10000L * match(setting_row$setting_id, setting_grid$setting_id) + 100000L * match(dataset_obj$dataset_id, c("sim_sst", "moreno_sheep"))
  K_expected <- if (!is.null(dataset_obj$truth)) dataset_obj$truth$K_true else 5L
  fits <- if (mc_cores > 1L) {
    parallel::mclapply(chain_seeds, function(sd) {
      fit_one_chain(dataset_obj$A, dataset_obj$dataset_id, model, setting_row, sd, K_expected)
    }, mc.cores = mc_cores)
  } else {
    lapply(chain_seeds, function(sd) fit_one_chain(dataset_obj$A, dataset_obj$dataset_id, model, setting_row, sd, K_expected))
  }
  names(fits) <- paste0("chain", seq_along(fits))
  fits
}

format_k_summary <- function(K_trace) {
  ks <- summarise_K_trace_local(K_trace)
  sprintf("$%d\\;[%d,%d]$", ks$mode, ks$lo, ks$hi)
}

format_num <- function(x, digits = 2) {
  vapply(x, function(value) {
    if (!is.finite(value)) return("NA")
    formatC(value, format = "f", digits = digits)
  }, character(1))
}

format_pct <- function(x, digits = 1) {
  vapply(x, function(value) {
    if (!is.finite(value)) return("NA")
    paste0(formatC(100 * value, format = "f", digits = digits), "\\%")
  }, character(1))
}

write_tex_table <- function(lines, path) {
  writeLines(lines, con = path)
  cat("Wrote:", path, "\n")
}

build_convergence_table_tex <- function(df) {
  body <- apply(df, 1, function(row) {
    sprintf("%s & %s & %s & %s & %s & %s & %s %s%s",
            row[["Dataset"]], row[["Setting"]], row[["Model"]],
            row[["RhatK"]], row[["RhatEta"]], row[["RhatVol"]], row[["RhatMax"]],
            intToUtf8(92), intToUtf8(92))
  })
  c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\caption{Small convergence study based on \\citet{gelman1992}. For each dataset, prior setting, and fitted model we ran ", n_chains,
           " chains with ", n_iter, " iterations (burn-in ", burn, ", thin ", thin,
           "). Entries report the Gelman--Rubin point estimate for the occupied block count $K$, mean $\\eta$, and mean volume parameter, together with the largest upper 97.5\\% confidence limit across the monitored summaries. Values near 1 indicate good mixing.\\label{tab:supp-small-gr}}"),
    "\\begin{tabular}{lllrrrr}",
    "\\toprule",
    paste0("Dataset & Setting & Model & $\\hat R_K$ & $\\hat R_{\\eta}$ & $\\hat R_{\\mathrm{vol}}$ & max upper CI ", intToUtf8(92), intToUtf8(92)),
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

build_sim_table_tex <- function(df) {
  body <- apply(df, 1, function(row) {
    sprintf("%s & %s & %s & %s & %s %s%s",
            row[["Setting"]], row[["Model"]], row[["Khat"]], row[["VI"]], row[["RhatMax"]],
            intToUtf8(92), intToUtf8(92))
  })
  c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Sensitivity on one SST-generated dataset ($n=60$, $K^\\star=5$, moderate density and strong hierarchy). For each model and prior setting we pool all retained draws across chains, report the posterior block-count summary $\\hat K[2.5\\%,97.5\\%]$, the VI loss of the pooled \\texttt{minVI} partition against the truth, and the largest Gelman--Rubin upper confidence limit from Table~\\ref{tab:supp-small-gr}. Lower VI is better.\\label{tab:supp-small-sim-sensitivity}}",
    "\\begin{tabular}{llccc}",
    "\\toprule",
    paste0("Setting & Model & $\\hat K[2.5\\%,97.5\\%]$ & VI & max upper CI ", intToUtf8(92), intToUtf8(92)),
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

build_sheep_table_tex <- function(df) {
  body <- apply(df, 1, function(row) {
    sprintf("%s & %s & %s & %s & %s %s%s",
            row[["Setting"]], row[["Model"]], row[["Khat"]], row[["Violation"]], row[["RhatMax"]],
            intToUtf8(92), intToUtf8(92))
  })
  c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Sensitivity on the bighorn sheep network. We report the pooled posterior summary of $K$, the backward-flow rate $\\hat\\zeta_z$ at the pooled \\texttt{minVI} partition (with DC--SBM reordered by empirical block score), and the largest Gelman--Rubin upper confidence limit. Stable $K$ and small changes in $\\hat\\zeta_z$ indicate that the substantive ordered picture is robust to moderate prior perturbations.\\label{tab:supp-small-sheep-sensitivity}}",
    "\\begin{tabular}{llccc}",
    "\\toprule",
    paste0("Setting & Model & $\\hat K[2.5\\%,97.5\\%]$ & $\\hat\\zeta_z$ & max upper CI ", intToUtf8(92), intToUtf8(92)),
    "\\midrule",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

study_datasets <- list(build_simulated_dataset(), build_real_dataset())
all_fits <- list()
conv_rows <- list()
sim_rows <- list()
sheep_rows <- list()

for (dataset_obj in study_datasets) {
  for (setting_idx in seq_len(nrow(setting_grid))) {
    setting_row <- setting_grid[setting_idx, ]
    for (model in model_levels) {
      fit_key <- paste(dataset_obj$dataset_id, setting_row$setting_id, model, sep = "__")
      fits <- run_group(dataset_obj, setting_row, model)
      all_fits[[fit_key]] <- fits

      conv <- convergence_summary(fits, model)
      conv_rows[[length(conv_rows) + 1L]] <- data.frame(
        dataset = dataset_obj$dataset_id,
        setting_id = setting_row$setting_id,
        setting_label = setting_row$setting_label,
        model = model,
        model_label = model_labels[[model]],
        rhat_K = conv$rhat_K,
        rhat_eta = conv$rhat_eta,
        rhat_volume = conv$rhat_volume,
        rhat_direction = conv$rhat_direction,
        rhat_upper_max = conv$rhat_upper_max,
        stringsAsFactors = FALSE
      )

      K_trace_pool <- unlist(lapply(fits, extract_K_trace_local))
      z_pool <- do.call(rbind, lapply(fits, function(f) z_to_matrix_local(f$z)))
      z_hat <- estimate_partition(z_pool)
      rhat_max_fmt <- format_num(conv$rhat_upper_max, 3)

      if (!is.null(dataset_obj$truth)) {
        sim_rows[[length(sim_rows) + 1L]] <- data.frame(
          dataset = dataset_obj$dataset_id,
          setting_id = setting_row$setting_id,
          setting_label = setting_row$setting_label,
          model = model,
          model_label = model_labels[[model]],
          Khat = format_k_summary(K_trace_pool),
          VI = format_num(vi_between(z_hat, dataset_obj$truth$z_true), 3),
          RhatMax = rhat_max_fmt,
          stringsAsFactors = FALSE
        )
      } else {
        sheep_rows[[length(sheep_rows) + 1L]] <- data.frame(
          dataset = dataset_obj$dataset_id,
          setting_id = setting_row$setting_id,
          setting_label = setting_row$setting_label,
          model = model,
          model_label = model_labels[[model]],
          Khat = format_k_summary(K_trace_pool),
          Violation = format_pct(compute_violation_rate(dataset_obj$A, z_hat, reorder = (model == "DCSBM")), 1),
          RhatMax = rhat_max_fmt,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

conv_df <- do.call(rbind, conv_rows)
sim_df <- do.call(rbind, sim_rows)
sheep_df <- do.call(rbind, sheep_rows)

utils::write.csv(conv_df, file.path(out_dir, "convergence_summary.csv"), row.names = FALSE)
utils::write.csv(sim_df, file.path(out_dir, "simulation_sensitivity_summary.csv"), row.names = FALSE)
utils::write.csv(sheep_df, file.path(out_dir, "moreno_sheep_sensitivity_summary.csv"), row.names = FALSE)

conv_tex_df <- data.frame(
  Dataset = dataset_labels[conv_df$dataset],
  Setting = conv_df$setting_label,
  Model = unname(conv_df$model_label),
  RhatK = format_num(conv_df$rhat_K, 3),
  RhatEta = format_num(conv_df$rhat_eta, 3),
  RhatVol = format_num(conv_df$rhat_volume, 3),
  RhatMax = format_num(conv_df$rhat_upper_max, 3),
  stringsAsFactors = FALSE
)

sim_tex_df <- data.frame(
  Setting = sim_df$setting_label,
  Model = unname(sim_df$model_label),
  Khat = sim_df$Khat,
  VI = sim_df$VI,
  RhatMax = sim_df$RhatMax,
  stringsAsFactors = FALSE
)

sheep_tex_df <- data.frame(
  Setting = sheep_df$setting_label,
  Model = unname(sheep_df$model_label),
  Khat = sheep_df$Khat,
  Violation = sheep_df$Violation,
  RhatMax = sheep_df$RhatMax,
  stringsAsFactors = FALSE
)

conv_tex_path <- file.path(paper_table_dir, "supp_small_gr_table.tex")
sim_tex_path <- file.path(paper_table_dir, "supp_small_sim_sensitivity_table.tex")
sheep_tex_path <- file.path(paper_table_dir, "supp_small_sheep_sensitivity_table.tex")

write_tex_table(build_convergence_table_tex(conv_tex_df), conv_tex_path)
write_tex_table(build_sim_table_tex(sim_tex_df), sim_tex_path)
write_tex_table(build_sheep_table_tex(sheep_tex_df), sheep_tex_path)

manifest_path <- file.path(out_dir, "MANIFEST.txt")
writeLines(c(
  paste("Generated:", format(Sys.time())),
  paste("Convergence CSV:", file.path(out_dir, "convergence_summary.csv")),
  paste("Simulation sensitivity CSV:", file.path(out_dir, "simulation_sensitivity_summary.csv")),
  paste("Sheep sensitivity CSV:", file.path(out_dir, "moreno_sheep_sensitivity_summary.csv")),
  paste("Paper tables:", conv_tex_path, sim_tex_path, sheep_tex_path)
), con = manifest_path)
cat("Wrote:", manifest_path, "\n")

cat("\nStudy complete.\n")
cat("Convergence table:", conv_tex_path, "\n")
cat("Simulation table:", sim_tex_path, "\n")
cat("Sheep table:", sheep_tex_path, "\n")
