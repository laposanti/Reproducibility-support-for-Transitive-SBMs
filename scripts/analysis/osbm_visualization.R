# =============================================================================
# Visualization Script: OSBM (WST/SST) Heatmaps + Ordered Network Plot
# =============================================================================
# This script creates:
# 1. Annotated heatmap of posterior mean rho_kl for WST
# 2. Bar plot of posterior dominance parameters psi_d for SST
# 3. Combined visualization (rho heatmap + psi bars side by side)
# 4. Network plot with hierarchical layout: higher-ranked blocks on top,
#    forward edges in blue, backward (violation) edges in red
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(mcclust)
  library(mcclust.ext)
  library(igraph)
  library(ggraph)
  library(lpSolve)
  library(scales)
  library(viridis)
  library(coda)
})

# Source transitivity check helpers (needed for violation & ordering diagnostics)
source("helper_folder/transitivity_check_helper.R")

# --- Load adjacency matrix for a named dataset ---
choose_dataset_local <- function(dataset) {
  if (dataset == "mountain_goats") {
    matrix_files <- list.files("./data/ShizukaMcDonald_Data",
                               full.names = TRUE, pattern = "\\.csv$")
    n_each <- vapply(matrix_files, function(f) nrow(read.csv(f, row.names = 1)),
                     FUN.VALUE = integer(1))
    A <- as.matrix(read.csv(matrix_files[which.max(n_each)],
                            row.names = 1, check.names = FALSE))
  } else if (dataset == "citations_data") {
    A <- read.csv("./data/Citations_application/cross-citation-matrix.csv",
                  row.names = 1, header = TRUE, check.names = FALSE)
    diag(A) <- 0
  } else if (dataset == "macaques_data") {
    edge_list <- read.table("./data/macaques/out.moreno.txt")
    nodes <- sort(unique(c(edge_list[[1]], edge_list[[2]])))
    A <- matrix(0L, length(nodes), length(nodes), dimnames = list(nodes, nodes))
    for (i in seq_len(nrow(edge_list))) {
      A[edge_list[i, 1], edge_list[i, 2]] <- edge_list[i, "V3"]
    }
  } else if (dataset == "high_school") {
    edge_list <- read.csv("./data/high-school/edges.csv", comment.char = "#", header = FALSE)
    colnames(edge_list) <- c("source", "target", "weight")
    nodes <- sort(unique(c(edge_list$source, edge_list$target)))
    node_names <- as.character(nodes)
    A <- matrix(0, length(nodes), length(nodes), dimnames = list(node_names, node_names))
    for (i in seq_len(nrow(edge_list))) {
      A[as.character(edge_list$source[i]), as.character(edge_list$target[i])] <- edge_list$weight[i]
    }
  } else if (dataset == "moreno_sheep") {
    edge_list <- read.csv("./data/moreno_sheep/edges.csv", comment.char = "#", header = FALSE)
    colnames(edge_list) <- c("source", "target", "weight")
    nodes <- sort(unique(c(edge_list$source, edge_list$target)))
    node_names <- as.character(nodes)
    A <- matrix(0, length(nodes), length(nodes), dimnames = list(node_names, node_names))
    for (i in seq_len(nrow(edge_list))) {
      A[as.character(edge_list$source[i]), as.character(edge_list$target[i])] <- edge_list$weight[i]
    }
  } else if (dataset == "strauss_2019b") {
    edge_list <- read.csv("./data/Strauss_2019b/edges.csv", comment.char = "#", header = FALSE)

    colnames(edge_list) <- c("source", "target", "time", "date")
    edge_agg <- aggregate(time ~ source + target, data = edge_list, FUN = length)
    colnames(edge_agg) <- c("source", "target", "weight")
    nodes <- sort(unique(c(edge_agg$source, edge_agg$target)))
    node_names <- as.character(nodes)
    A <- matrix(0, length(nodes), length(nodes), dimnames = list(node_names, node_names))
    for (i in seq_len(nrow(edge_agg))) {
      A[as.character(edge_agg$source[i]), as.character(edge_agg$target[i])] <- edge_agg$weight[i]
    }
  } else {
    stop("choose_dataset_local: unknown dataset '", dataset, "'")
  }
  A <- as.matrix(A)
  stopifnot(nrow(A) == ncol(A))
  colnames(A) <- rownames(A)
  A
}

# --- Find fit file: run dir first (canonical, newer), ppc as fallback ---
find_fit_file_local <- function(dataset_name, model, ppc_dir) {
  fm <- if (model %in% c("DC-SBM", "DCSBM")) "DCSBM" else model
  # Prefer canonical run dir fit
  run_dir <- Sys.getenv("APP_RUN_DIR",
    unset = "./output/application/raw/application_run_20260411_163055")
  cand <- file.path(run_dir, paste0(dataset_name, "_", fm, "_fit.rds"))
  if (file.exists(cand)) return(cand)
  # Fallback: ppc directory (older run)
  ppc_path <- file.path(ppc_dir, dataset_name, fm,
                        paste0(dataset_name, "_", fm, "_fit.rds"))
  if (file.exists(ppc_path)) {
    message("osbm_visualization: using fallback fit from ppc dir: ", basename(ppc_path))
    return(ppc_path)
  }
  NULL  # not found
}

# --- Helper: Convert z list/matrix to matrix ---
z_to_matrix <- function(z) {
  if (is.list(z)) {
    do.call(rbind, z)
  } else if (is.matrix(z)) {
    z
  } else {
    matrix(z, nrow = 1)
  }
}

# --- Get point estimate of z using minVI ---
# Uses method="all" to match the canonical pipeline in app_var_analyze_new.R
get_z_hat <- function(z_mat) {
  psm <- mcclust::comp.psm(z_mat)
  res <- mcclust.ext::minVI(psm, cls.draw = z_mat, method = "all")
  z_hat <- if (is.matrix(res$cl)) res$cl[1, ] else res$cl
  as.integer(z_hat)
}

# --- Relabel blocks by decreasing per-node success score ---
# Matches the canonical pipeline in app_var_analyze_new.R: psm_and_zhat()
relabel_by_strength <- function(z_hat, A) {
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  K_hat <- max(z_hat, na.rm = TRUE)
  
  # Per-node success score (same as canonical pipeline)
  outdeg   <- rowSums(A)
  matches  <- A + t(A)
  success  <- outdeg / pmax(rowSums(matches), 1L)
  
  # Block means under z_hat
  means <- rep(-Inf, K_hat)
  for (k in seq_len(K_hat)) {
    idx <- which(z_hat == k)
    if (length(idx) > 0) means[k] <- mean(success[idx])
  }
  
  # Order blocks by decreasing strength (strongest = 1)
  ord_old   <- order(means, decreasing = TRUE, na.last = TRUE)
  p_old2new <- integer(K_hat)
  p_old2new[ord_old] <- seq_len(K_hat)
  
  z_new <- p_old2new[z_hat]
  block_map <- setNames(seq_len(K_hat), as.character(ord_old))
  
  list(z_new = as.integer(z_new), ord_old = ord_old, block_map = block_map)
}

# --- Pretty dataset/network title ---
# Canonical names matching Table~\ref{tab:datasets} (short forms used in
# tab:hierarchy-synopsis for plot titles).
pretty_dataset_label <- function(dataset_name) {
  key <- tolower(trimws(dataset_name))
  dplyr::case_when(
    key == "moreno_sheep"   ~ "Bighorn sheep",
    key == "strauss_2019b"  ~ "Spotted hyenas",
    key == "mountain_goats" ~ "Mountain goats",
    key == "citations_data" ~ "Stat. journals",
    key == "macaques_data"  ~ "Japanese macaques",
    key == "high_school"    ~ "High school",
    TRUE ~ tools::toTitleCase(gsub("_", " ", key))
  )
}

pretty_network_title <- function(dataset_name) {
  pretty_dataset_label(dataset_name)
}

pretty_model_label <- function(model_key) {
  dplyr::case_when(
    identical(model_key, "DCSBM") ~ "DC-SBM",
    TRUE ~ as.character(model_key)
  )
}

application_best_model <- function(dataset_name) {
  key <- tolower(trimws(dataset_name))
  dplyr::case_when(
    key == "mountain_goats" ~ "SST",
    key == "macaques_data" ~ "DCSBM",
    key %in% c("citations_data", "high_school") ~ "WST",
    TRUE ~ NA_character_
  )
}

save_plot_bundle <- function(plot_obj, out_dir, stem, width, height,
                             dpi = 360, bg = "white") {
  fs::dir_create(out_dir)
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))

  ggsave(
    filename = pdf_path,
    plot = plot_obj,
    width = width,
    height = height,
    device = grDevices::pdf,
    bg = bg,
    limitsize = FALSE
  )
  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )

  invisible(list(pdf = pdf_path, png = png_path))
}

save_plot_bundle_dirs <- function(plot_obj, dirs, stem, width, height,
                                  dpi = 360, bg = "white") {
  dirs <- unique(dirs[!is.na(dirs) & nzchar(dirs)])
  lapply(dirs, function(dir_path) {
    save_plot_bundle(
      plot_obj = plot_obj,
      out_dir = dir_path,
      stem = stem,
      width = width,
      height = height,
      dpi = dpi,
      bg = bg
    )
  })
}

# --- WST: Compute posterior mean rho matrix ---
compute_wst_rho_posterior <- function(fit, z_hat_new, block_map, modal_K = NULL) {
  psi_list <- fit$psi
  K_trace <- fit$K_trace
  S <- length(psi_list)
  K_out <- length(unique(z_hat_new))
  
  if (is.null(modal_K)) {
    modal_K <- K_out
  }
  if (is.null(K_trace)) {
    K_trace <- vapply(psi_list, function(x) if (is.matrix(x)) nrow(x) else NA_integer_, integer(1))
  }
  
  rho_sum <- matrix(0, K_out, K_out)
  count <- 0
  
  # Only use draws with the displayed K.
  modal_idx <- which(K_trace == modal_K)
  
  for (s in modal_idx) {
    psi_s <- psi_list[[s]]
    if (!is.matrix(psi_s)) next
    K_s <- nrow(psi_s)
    if (K_s != modal_K) next
    
    # Permute psi to match the new block ordering
    # block_map maps old -> new
    old_labels <- as.integer(names(block_map))
    
    # Only use blocks that exist in both
    valid_old <- old_labels[old_labels <= K_s]
    if (length(valid_old) < 2) next
    
    # Extract sub-matrix in the new order
    new_labels <- block_map[as.character(valid_old)]
    
    for (i in seq_along(valid_old)) {
      for (j in seq_along(valid_old)) {
        if (i == j) next
        oi <- valid_old[i]
        oj <- valid_old[j]
        ni <- new_labels[i]
        nj <- new_labels[j]
        
        val <- psi_s[oi, oj]
        rho_ij <- 1 / (1 + exp(-val))  # logistic transformation
        
        if (ni <= K_out && nj <= K_out) {
          rho_sum[ni, nj] <- rho_sum[ni, nj] + rho_ij
        }
      }
    }
    count <- count + 1
  }
  
  if (count > 0) {
    rho_mean <- rho_sum / count
  } else {
    rho_mean <- matrix(0.5, K_out, K_out)
  }
  
  diag(rho_mean) <- 0.5
  rho_mean
}

# --- SST: Compute posterior mean psi_d vector ---
compute_sst_psi_posterior <- function(fit, K_model = NULL, max_K = NULL) {
  psi_list <- fit$psi
  K_trace <- fit$K_trace
  S <- length(psi_list)

  if (is.null(K_trace)) {
    K_trace <- vapply(psi_list, function(x) length(x) + 1L, integer(1))
  }
  
  if (is.null(K_model)) {
    if (!is.null(max_K)) {
      K_model <- max_K + 1L
    } else {
      K_model <- as.integer(names(sort(table(K_trace), decreasing = TRUE))[1])
    }
  }
  max_K <- K_model - 1L
  
  # Accumulate psi for each distance d
  psi_sum <- rep(0, max_K)
  psi_count <- rep(0, max_K)
  
  for (s in which(K_trace == K_model)) {
    psi_s <- psi_list[[s]]
    if (is.null(psi_s) || length(psi_s) < max_K) next
    
    psi_sum <- psi_sum + psi_s[seq_len(max_K)]
    psi_count <- psi_count + 1
  }
  
  psi_mean <- ifelse(psi_count > 0, psi_sum / psi_count, NA)
  
  data.frame(
    distance = seq_len(max_K),
    psi_mean = psi_mean,
    rho_mean = 1 / (1 + exp(-psi_mean))  # dominance probability
  )
}

# --- DCSBM: Compute empirical block-level statistics from A using z_hat ---
# Note: We compute directly from A to avoid label-switching issues with MCMC lambda
compute_dcsbm_empirical <- function(A, z_hat_new) {
  blocks <- sort(unique(z_hat_new))
  K <- length(blocks)
  
  # Compute empirical block-to-block counts
  block_A <- matrix(0, K, K)
  block_N <- matrix(0, K, K)  # Number of possible edges
  
  for (i in seq_len(K)) {
    idx_i <- which(z_hat_new == blocks[i])
    n_i <- length(idx_i)
    for (j in seq_len(K)) {
      idx_j <- which(z_hat_new == blocks[j])
      n_j <- length(idx_j)
      block_A[i, j] <- sum(A[idx_i, idx_j, drop = FALSE])
      if (i == j) {
        block_N[i, j] <- n_i * (n_i - 1)  # Exclude diagonal
      } else {
        block_N[i, j] <- n_i * n_j
      }
    }
  }
  
  # Compute empirical lambda (rate per dyad)
  lambda_emp <- block_A / pmax(block_N, 1)
  
  # Compute rho: rho_kl = A_kl / (A_kl + A_lk)
  rho_emp <- matrix(0.5, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      if (i == j) {
        rho_emp[i, j] <- 0.5
      } else {
        total <- block_A[i, j] + block_A[j, i]
        if (total > 0) {
          rho_emp[i, j] <- block_A[i, j] / total
        } else {
          rho_emp[i, j] <- 0.5
        }
      }
    }
  }
  
  list(
    block_A = block_A,      # Raw edge counts
    block_N = block_N,      # Possible edges
    lambda_emp = lambda_emp, # Empirical rate
    rho_emp = rho_emp       # Empirical dominance probability
  )
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

# --- WST Rho Heatmap ---
plot_wst_rho_heatmap <- function(rho_mat, title = "WST: Block Interaction Matrix") {
  K <- nrow(rho_mat)
  
  df <- expand.grid(From = 1:K, To = 1:K) %>%
    mutate(rho = as.vector(rho_mat),
           label = ifelse(From == To, "", sprintf("%.2f", rho)))
  
  # Reverse y-axis so block 1 is on top
  ggplot(df, aes(x = factor(To), y = factor(From, levels = K:1), fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_viridis(option = "D", limits = c(0, 1), 
                       name = expression(hat(rho)[kl])) +
    labs(x = "Block l (weaker →)", 
         y = "Block k (stronger →)",
         title = title) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed()
}

# --- SST Psi Bar Plot ---
plot_sst_psi_bars <- function(psi_df, title = "SST: Distance-to-Dominance") {
  # Filter out NA and limit to reasonable range
  psi_df <- psi_df %>% filter(!is.na(psi_mean))
  
  ggplot(psi_df, aes(x = factor(distance), y = rho_mean)) +
    geom_col(fill = "#3B528B", width = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = sprintf("%.2f", rho_mean)), 
              vjust = -0.3, size = 3) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Distance d between blocks",
         y = expression("Dominance prob " ~ rho[d]),
         title = title) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major.x = element_blank()
    )
}

# --- DCSBM Lambda Heatmap ---
plot_dcsbm_lambda_heatmap <- function(lambda_mat, title = "DC-SBM: Block Interaction Rates") {
  K <- nrow(lambda_mat)
  
  df <- expand.grid(From = 1:K, To = 1:K) %>%
    mutate(lambda = as.vector(lambda_mat),
           label = sprintf("%.2f", lambda))
  
  # Reverse y-axis so block 1 is on top
  ggplot(df, aes(x = factor(To), y = factor(From, levels = K:1), fill = lambda)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_viridis(option = "C", 
                       name = expression(hat(lambda)[kl])) +
    labs(x = "Block l", 
         y = "Block k",
         title = title) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed()
}

# --- DCSBM Rho (Dominance) Heatmap ---
plot_dcsbm_rho_heatmap <- function(rho_mat, title = "DC-SBM: Dominance Probabilities") {
  K <- nrow(rho_mat)
  
  df <- expand.grid(From = 1:K, To = 1:K) %>%
    mutate(rho = as.vector(rho_mat),
           label = ifelse(From == To, "", sprintf("%.2f", rho)))
  
  # Reverse y-axis so block 1 is on top
  ggplot(df, aes(x = factor(To), y = factor(From, levels = K:1), fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_viridis(option = "D", limits = c(0, 1), 
                       name = expression(hat(rho)[kl])) +
    labs(x = "Block l (weaker →)", 
         y = "Block k (stronger →)",
         title = title) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed()
}

# --- SST Psi as Vertical Heatmap (Kx1) ---
plot_sst_psi_vertical <- function(psi_df, title = "SST", limits = c(0, 1)) {
  # Filter out NA and prepare for vertical heatmap
  psi_df <- psi_df %>% filter(!is.na(psi_mean))
  K <- max(psi_df$distance) + 1  # Distance goes from 1 to K-1, so K = max + 1
  
  # Create dataframe: distance d maps to probability of dominance at that distance
  df <- psi_df %>%
    mutate(
      block_from = 1,  # dummy x-axis
      block_to = distance,
      rho = rho_mean,
      label = sprintf("%.2f", rho_mean)
    )
  
  # Reverse y-axis so distance 1 is on top
  max_d <- max(df$distance)
  ggplot(df, aes(x = factor(1), y = factor(distance, levels = max_d:1), fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 3, color = "black") +
    scale_fill_viridis(option = "D", limits = limits,
                       name = expression(hat(rho)[d])) +
    labs(x = NULL, 
         y = "Distance d",
         title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "none"  # Will use shared legend
    )
}

# --- Combined WST Heatmap + SST Vertical + DCSBM Heatmap ---
plot_combined_rho_all <- function(rho_wst, psi_df, rho_dcsbm, dataset_name) {
  # Harmonize color scale: use [0, 1] for all
  limits <- c(0, 1)
  
  # WST heatmap (left)
  K_wst <- nrow(rho_wst)
  df_wst <- expand.grid(From = 1:K_wst, To = 1:K_wst) %>%
    mutate(rho = as.vector(rho_wst),
           label = ifelse(From == To, "", sprintf("%.2f", rho)))
  
  p_wst <- ggplot(df_wst, aes(x = factor(To), y = factor(From, levels = K_wst:1), fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.5, color = "black") +
    scale_fill_viridis(option = "D", limits = limits,
                       name = expression(hat(rho))) +
    labs(x = "Block l", y = "Block k", title = "WST") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    coord_fixed()
  
  # SST vertical heatmap (middle)
  psi_clean <- psi_df %>% filter(!is.na(psi_mean))
  max_d <- max(psi_clean$distance)
  df_sst <- psi_clean %>%
    mutate(rho = rho_mean, label = sprintf("%.2f", rho_mean))
  
  p_sst <- ggplot(df_sst, aes(x = factor(1), y = factor(distance, levels = max_d:1), fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.5, color = "black") +
    scale_fill_viridis(option = "D", limits = limits) +
    labs(x = NULL, y = "Distance d", title = "SST") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  # DCSBM heatmap (right)
  K_dc <- nrow(rho_dcsbm)
  df_dc <- expand.grid(From = 1:K_dc, To = 1:K_dc) %>%
    mutate(rho = as.vector(rho_dcsbm),
           label = ifelse(From == To, "", sprintf("%.2f", rho)))
  
  p_dc <- ggplot(df_dc, aes(x = factor(To), y = factor(From, levels = K_dc:1), fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.5, color = "black") +
    scale_fill_viridis(option = "D", limits = limits,
                       name = expression(hat(rho)[kl])) +
    labs(x = "Block l", y = "Block k", title = "DC-SBM") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_fixed()
  
  # Combine with patchwork - WST | SST | DCSBM
  p_wst + p_sst + p_dc + 
    plot_layout(widths = c(K_wst, 1.5, K_dc)) +
    plot_annotation(
      title = paste("Dominance Probabilities -", dataset_name),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
}

# --- Original Combined WST Heatmap + SST Bars (kept for backward compat) ---
plot_combined_heatmap_bars <- function(rho_mat, psi_df, dataset_name) {
  p1 <- plot_wst_rho_heatmap(rho_mat, title = paste("WST:", dataset_name))
  p2 <- plot_sst_psi_bars(psi_df, title = paste("SST:", dataset_name))
  
  p1 + p2 + plot_layout(widths = c(2, 1)) +
    plot_annotation(
      title = paste("Posterior Block Structure -", dataset_name),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
}

# =============================================================================
# VI DISTANCE TABLE BETWEEN PARTITIONS
# =============================================================================

#' Compute Variation of Information (VI) distance between two partitions
compute_vi_distance <- function(z1, z2) {
  # Ensure both have same length
  stopifnot(length(z1) == length(z2))
  n <- length(z1)
  
  # Get cluster labels
  labels1 <- sort(unique(z1))
  labels2 <- sort(unique(z2))
  K1 <- length(labels1)
  K2 <- length(labels2)
  
  # Compute contingency table
  cont <- matrix(0, K1, K2)
  for (i in seq_len(K1)) {
    for (j in seq_len(K2)) {
      cont[i, j] <- sum(z1 == labels1[i] & z2 == labels2[j])
    }
  }
  
  # Cluster sizes
  n1 <- rowSums(cont)
  n2 <- colSums(cont)
  
  # Entropies
  H1 <- -sum((n1/n) * log(n1/n + 1e-10))
  H2 <- -sum((n2/n) * log(n2/n + 1e-10))
  
  # Mutual information
  MI <- 0
  for (i in seq_len(K1)) {
    for (j in seq_len(K2)) {
      if (cont[i, j] > 0) {
        MI <- MI + (cont[i, j]/n) * log((n * cont[i, j]) / (n1[i] * n2[j]))
      }
    }
  }
  
  # VI = H(z1) + H(z2) - 2*MI
  VI <- H1 + H2 - 2 * MI
  
  # Normalized VI (divide by log(n) to get 0-1 range approximately)
  VI_norm <- VI / log(n)
  
  list(VI = VI, VI_norm = VI_norm, H1 = H1, H2 = H2, MI = MI, K1 = K1, K2 = K2)
}

#' Compute VI distance matrix between all model partitions
compute_vi_distance_table <- function(z_wst, z_sst, z_dcsbm) {
  models <- c("WST", "SST", "DCSBM")
  z_list <- list(WST = z_wst, SST = z_sst, DCSBM = z_dcsbm)
  
  # Remove NULL entries
  z_list <- z_list[!sapply(z_list, is.null)]
  models <- names(z_list)
  n_models <- length(models)
  
  if (n_models < 2) return(NULL)
  
  # Compute pairwise VI distances
  results <- list()
  for (i in seq_len(n_models - 1)) {
    for (j in (i + 1):n_models) {
      vi_res <- compute_vi_distance(z_list[[i]], z_list[[j]])
      results[[length(results) + 1]] <- data.frame(
        Model1 = models[i],
        Model2 = models[j],
        K1 = vi_res$K1,
        K2 = vi_res$K2,
        VI = round(vi_res$VI, 4),
        VI_norm = round(vi_res$VI_norm, 4)
      )
    }
  }
  
  do.call(rbind, results)
}

# =============================================================================
# AGONY: EXACT RANK PARTITION + MODEL COMPARISON
# =============================================================================

compact_partition_labels <- function(z) {
  z_int <- as.integer(z)
  lev <- sort(unique(z_int))
  if (!length(lev)) return(z_int)
  map <- setNames(seq_along(lev), lev)
  as.integer(map[as.character(z_int)])
}

midranks_from_partition <- function(z_ordered) {
  z_compact <- compact_partition_labels(z_ordered)
  if (!length(z_compact)) return(numeric(0))

  block_sizes <- tabulate(z_compact, nbins = max(z_compact))
  block_start <- cumsum(c(1, head(block_sizes, -1)))
  block_end <- cumsum(block_sizes)
  block_mid <- (block_start + block_end) / 2
  as.numeric(block_mid[z_compact])
}

kendall_tau_b_from_ranks <- function(rank1, rank2) {
  stopifnot(length(rank1) == length(rank2))
  n <- length(rank1)
  if (n < 2L) {
    return(list(
      tau_b = NA_real_,
      distance = NA_real_,
      concordant = 0L,
      discordant = 0L,
      tied_rank1 = 0L,
      tied_rank2 = 0L,
      tied_both = 0L
    ))
  }

  concordant <- 0L
  discordant <- 0L
  tied_rank1 <- 0L
  tied_rank2 <- 0L
  tied_both <- 0L

  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      s1 <- sign(rank1[i] - rank1[j])
      s2 <- sign(rank2[i] - rank2[j])

      if (s1 == 0 && s2 == 0) {
        tied_both <- tied_both + 1L
      } else if (s1 == 0) {
        tied_rank1 <- tied_rank1 + 1L
      } else if (s2 == 0) {
        tied_rank2 <- tied_rank2 + 1L
      } else if (s1 == s2) {
        concordant <- concordant + 1L
      } else {
        discordant <- discordant + 1L
      }
    }
  }

  denom <- sqrt((concordant + discordant + tied_rank1) *
                  (concordant + discordant + tied_rank2))
  tau_b <- if (denom > 0) (concordant - discordant) / denom else NA_real_

  list(
    tau_b = tau_b,
    distance = if (is.finite(tau_b)) (1 - tau_b) / 2 else NA_real_,
    concordant = concordant,
    discordant = discordant,
    tied_rank1 = tied_rank1,
    tied_rank2 = tied_rank2,
    tied_both = tied_both
  )
}

compute_agony_score <- function(A, ranks) {
  stopifnot(length(ranks) == nrow(A), nrow(A) == ncol(A))
  rows <- row(A)
  cols <- col(A)
  vals <- as.numeric(A)
  mask <- (rows != cols) & (vals > 0)
  if (!any(mask)) return(0)
  sum(vals[mask] * pmax(ranks[rows[mask]] - ranks[cols[mask]] + 1, 0))
}

solve_agony_partition <- function(A, timeout = 60L) {
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  stopifnot(is.matrix(A), nrow(A) == ncol(A))

  A <- as.matrix(A)
  diag(A) <- 0
  n <- nrow(A)

  edge_idx <- which(A > 0, arr.ind = TRUE)
  edge_idx <- edge_idx[edge_idx[, 1] != edge_idx[, 2], , drop = FALSE]
  m <- nrow(edge_idx)

  if (n == 0L) {
    return(list(
      ranks = integer(0),
      raw_ranks = integer(0),
      score = 0,
      objective = 0,
      status = 0L,
      total_mass = 0,
      K = 0L
    ))
  }

  if (m == 0L) {
    ranks <- rep(1L, n)
    return(list(
      ranks = ranks,
      raw_ranks = rep(0L, n),
      score = 0,
      objective = 0,
      status = 0L,
      total_mass = 0,
      K = 1L
    ))
  }

  n_vars <- n + m
  n_const <- 2L * m + 1L

  objective <- numeric(n_vars)
  objective[(n + 1L):n_vars] <- A[edge_idx]

  const.mat <- matrix(0, nrow = n_vars, ncol = n_const)
  const.dir <- c(rep(">=", m), rep(">=", m), "=")
  const.rhs <- c(rep(1, m), rep(0, m), 0)

  for (e in seq_len(m)) {
    i <- edge_idx[e, 1]
    j <- edge_idx[e, 2]
    slack_idx <- n + e

    # slack_e >= r_i - r_j + 1
    const.mat[slack_idx, e] <- 1
    const.mat[i, e] <- -1
    const.mat[j, e] <- 1

    # slack_e >= 0
    const.mat[slack_idx, m + e] <- 1
  }

  # Anchor the translation invariance: node 1 has raw rank 0.
  const.mat[1, n_const] <- 1

  lp_res <- lp(
    direction = "min",
    objective.in = objective,
    const.mat = const.mat,
    const.dir = const.dir,
    const.rhs = const.rhs,
    int.vec = seq_len(n),
    transpose.constraints = FALSE,
    timeout = as.integer(timeout)
  )

  if (lp_res$status != 0L) {
    stop("Agony optimization failed with lpSolve status ", lp_res$status)
  }

  raw_ranks <- as.integer(round(lp_res$solution[seq_len(n)]))
  ranks <- compact_partition_labels(raw_ranks)
  score <- compute_agony_score(A, ranks)

  list(
    ranks = ranks,
    raw_ranks = raw_ranks,
    score = score,
    objective = lp_res$objval,
    status = lp_res$status,
    total_mass = sum(A[row(A) != col(A)]),
    K = length(unique(ranks))
  )
}

build_block_weight_matrix <- function(A, z_vec) {
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  z_compact <- compact_partition_labels(z_vec)
  K <- max(z_compact, na.rm = TRUE)
  W <- matrix(0, nrow = K, ncol = K)

  for (i in seq_len(nrow(A))) {
    for (j in seq_len(ncol(A))) {
      if (i == j || A[i, j] <= 0) next
      W[z_compact[i], z_compact[j]] <- W[z_compact[i], z_compact[j]] + A[i, j]
    }
  }

  list(W = W, z_compact = z_compact, K = K)
}

score_partition_with_agony <- function(A, z_vec, timeout = 30L) {
  blk <- build_block_weight_matrix(A, z_vec)
  blk_fit <- solve_agony_partition(blk$W, timeout = timeout)
  node_ranks <- blk_fit$ranks[blk$z_compact]

  list(
    model_partition = blk$z_compact,
    node_ranks = node_ranks,
    K_partition = blk$K,
    K_levels = length(unique(node_ranks)),
    score = compute_agony_score(A, node_ranks),
    block_agony = blk_fit
  )
}

compute_agony_outputs <- function(A, z_hats, dataset_name, out_path,
                                  timeout_node = 60L, timeout_block = 30L) {
  fs::dir_create(out_path)

  node_names <- if (!is.null(rownames(A))) rownames(A) else as.character(seq_len(nrow(A)))
  agony_fit <- solve_agony_partition(A, timeout = timeout_node)
  agony_midranks <- midranks_from_partition(agony_fit$ranks)

  agony_partition_table <- data.frame(
    dataset = dataset_name,
    node = seq_len(nrow(A)),
    node_name = node_names,
    agony_rank = agony_fit$ranks,
    agony_midrank = agony_midranks,
    stringsAsFactors = FALSE
  )
  write.csv(
    agony_partition_table,
    file.path(out_path, paste0(dataset_name, "_agony_partition.csv")),
    row.names = FALSE
  )

  total_mass <- agony_fit$total_mass
  agony_rows <- list(
    Agony = data.frame(
      dataset = dataset_name,
      source = "Agony",
      K_partition = agony_fit$K,
      K_levels = agony_fit$K,
      agony_score = agony_fit$score,
      agony_score_norm = if (total_mass > 0) agony_fit$score / total_mass else NA_real_,
      agony_gap = 0,
      agony_ratio = 1,
      VI_to_agony = 0,
      VI_norm_to_agony = 0,
      tau_b_midrank_to_agony = 1,
      kendall_tau_midrank_distance_to_agony = 0,
      stringsAsFactors = FALSE
    )
  )

  agony_kendall_rows <- list()

  for (mod in names(z_hats)) {
    mod_fit <- score_partition_with_agony(A, z_hats[[mod]], timeout = timeout_block)
    vi_res <- compute_vi_distance(mod_fit$model_partition, agony_fit$ranks)
    mod_midranks <- midranks_from_partition(z_hats[[mod]])
    tau_res <- kendall_tau_b_from_ranks(agony_midranks, mod_midranks)

    agony_rows[[mod]] <- data.frame(
      dataset = dataset_name,
      source = mod,
      K_partition = mod_fit$K_partition,
      K_levels = mod_fit$K_levels,
      agony_score = mod_fit$score,
      agony_score_norm = if (total_mass > 0) mod_fit$score / total_mass else NA_real_,
      agony_gap = mod_fit$score - agony_fit$score,
      agony_ratio = if (agony_fit$score > 0) mod_fit$score / agony_fit$score else NA_real_,
      VI_to_agony = vi_res$VI,
      VI_norm_to_agony = vi_res$VI_norm,
      tau_b_midrank_to_agony = tau_res$tau_b,
      kendall_tau_midrank_distance_to_agony = tau_res$distance,
      stringsAsFactors = FALSE
    )

    agony_kendall_rows[[mod]] <- data.frame(
      dataset = dataset_name,
      reference = "Agony",
      model = mod,
      n_nodes = length(agony_midranks),
      agony_K = agony_fit$K,
      model_K = length(unique(compact_partition_labels(z_hats[[mod]]))),
      tau_b = tau_res$tau_b,
      kendall_tau_distance = tau_res$distance,
      concordant_pairs = tau_res$concordant,
      discordant_pairs = tau_res$discordant,
      ties_agony_only = tau_res$tied_rank1,
      ties_model_only = tau_res$tied_rank2,
      ties_both = tau_res$tied_both,
      stringsAsFactors = FALSE
    )
  }

  agony_table <- do.call(rbind, agony_rows)
  rownames(agony_table) <- NULL
  agony_table <- agony_table %>%
    arrange(match(source, c("Agony", "WST", "SST", "DCSBM")))

  agony_kendall_table <- if (length(agony_kendall_rows) > 0) {
    do.call(rbind, agony_kendall_rows)
  } else {
    NULL
  }

  model_only <- agony_table %>% filter(source != "Agony")
  agony_winner <- if (nrow(model_only) > 0) {
    best_score <- model_only %>% arrange(agony_gap, agony_score, source) %>% slice(1)
    best_vi <- model_only %>% arrange(VI_norm_to_agony, VI_to_agony, source) %>% slice(1)
    best_tau <- model_only %>%
      arrange(kendall_tau_midrank_distance_to_agony, desc(tau_b_midrank_to_agony), source) %>%
      slice(1)

    data.frame(
      dataset = dataset_name,
      agony_score = agony_fit$score,
      agony_K = agony_fit$K,
      closest_model_by_score = best_score$source,
      closest_model_score = best_score$agony_score,
      closest_model_agony_gap = best_score$agony_gap,
      closest_model_by_vi = best_vi$source,
      closest_model_vi_norm = best_vi$VI_norm_to_agony,
      closest_model_by_kendall_midrank = best_tau$source,
      closest_model_kendall_midrank_distance = best_tau$kendall_tau_midrank_distance_to_agony,
      stringsAsFactors = FALSE
    )
  } else {
    NULL
  }

  write.csv(
    agony_table,
    file.path(out_path, paste0(dataset_name, "_agony_comparison.csv")),
    row.names = FALSE
  )
  if (!is.null(agony_winner)) {
    write.csv(
      agony_winner,
      file.path(out_path, paste0(dataset_name, "_agony_winner.csv")),
      row.names = FALSE
    )
  }
  if (!is.null(agony_kendall_table)) {
    write.csv(
      agony_kendall_table,
      file.path(out_path, paste0(dataset_name, "_agony_kendall_tau_table.csv")),
      row.names = FALSE
    )
  }

  cat("    Saved agony partition and comparison tables\n")
  print(agony_table)
  if (!is.null(agony_kendall_table)) print(agony_kendall_table)
  if (!is.null(agony_winner)) print(agony_winner)

  list(
    agony_fit = agony_fit,
    agony_partition_table = agony_partition_table,
    agony_table = agony_table,
    agony_kendall_table = agony_kendall_table,
    agony_winner = agony_winner
  )
}

# =============================================================================
# NETWORK VISUALIZATION WITH ORDER LAYOUT
# =============================================================================

#' Plot network with ordered layout
#' @param A Adjacency matrix
#' @param z_hat Block assignments (1 = strongest)
#' @param edge_threshold Minimum edge weight to show (for clarity)
#' @param node_size_scale Scaling factor for node sizes
plot_ordered_network <- function(A, z_hat, edge_threshold = 0, 
                                  node_size_scale = 3,
                                  title = "Ordered Network",
                                  subtitle_text = NULL,
                                  alpha = 0.5,
                                  layout_style = c("tier_cluster", "fr_ordered", "stress_ordered", "kk_ordered", "tier_fan", "tier_line"),
                                  tier_spacing = 1.9,
                                  tier_line_spacing_scale = 1.0,
                                  precomputed_rate = NULL,
                                  precomputed_backward = NULL,
                                  precomputed_cross = NULL) {
  n <- nrow(A)
  layout_style <- match.arg(layout_style)

  # Canonical contract (post-processing cube): when precomputed violation
  # stats are supplied, the supplied z_hat IS the canonical partition and
  # MUST NOT be re-ranked here. Re-ranking inside the plot used to silently
  # produce a different backward-mass than the hierarchy diagnostics
  # (e.g. 31.7% vs 29.6% for citations_data / WST).
  if (!is.null(precomputed_rate)) {
    z_ranked <- as.integer(z_hat)
  } else {
    z_ranked <- .rank_partition_by_strength(
      A = A, z_vec = z_hat, alpha = alpha,
      order_direction = "strong_to_weak"
    )
  }
  K <- max(z_ranked, na.rm = TRUE)
  
  # Create igraph object
  g <- igraph::graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE)
  
  # Node attributes: use z_ranked for layout and direction
  V(g)$block <- z_ranked
  V(g)$tier <- factor(
    z_ranked,
    levels = seq_len(K),
    labels = paste("Tier", seq_len(K))
  )
  V(g)$block_rank <- z_ranked  # 1 = strongest
  
  # Edge attributes: forward vs backward
  edge_df <- igraph::as_data_frame(g, what = "edges")
  
  # Convert from/to to numeric indices
  node_ids <- V(g)$name
  if (is.null(node_ids)) node_ids <- as.character(1:n)
  
  from_idx <- match(edge_df$from, node_ids)
  to_idx <- match(edge_df$to, node_ids)
  
  edge_df$from_block <- z_ranked[from_idx]
  edge_df$to_block <- z_ranked[to_idx]
  
  # Direction: forward if from stronger (lower rank) to weaker (higher rank)
  # Remember: block 1 = strongest, higher number = weaker
  edge_df$is_forward <- edge_df$from_block < edge_df$to_block
  edge_df$is_backward <- edge_df$from_block > edge_df$to_block
  edge_df$is_within <- edge_df$from_block == edge_df$to_block

  # Node-level strength used for deterministic placement inside each tier.
  node_matches <- rowSums(A + t(A))
  node_success <- rowSums(A) / pmax(node_matches, 1)
  node_outdeg <- rowSums(A)
  
  # Compute backward edge statistics for annotation (weighted by edge weight)
  # Only cross-block edges count as potential violations (matching transitivity_check_helper)
  cross_block_edges <- edge_df %>% filter(!is_within)
  total_cross_block <- sum(cross_block_edges$weight)
  total_backward <- sum(cross_block_edges$weight[cross_block_edges$is_backward])
  pct_backward <- if (total_cross_block > 0) 100 * total_backward / total_cross_block else 0

  # Override with canonical cube values when supplied — these are the
  # numbers reported in hierarchy_diagnostics_overview.csv and the paper.
  if (!is.null(precomputed_rate)) {
    pct_backward <- 100 * as.numeric(precomputed_rate)
    if (!is.null(precomputed_backward)) total_backward    <- as.numeric(precomputed_backward)
    if (!is.null(precomputed_cross))    total_cross_block <- as.numeric(precomputed_cross)
  }
  
  E(g)$direction <- dplyr::case_when(
    edge_df$is_forward ~ "forward",
    edge_df$is_backward ~ "backward",
    TRUE ~ "within"
  )
  E(g)$color <- dplyr::case_when(
    edge_df$is_forward ~ "#2166AC",   # Blue for forward
    edge_df$is_backward ~ "#B2182B",  # Red for backward
    TRUE ~ "#808080"                  # Gray for within
  )
  
  # Tier centers used by all ordered layouts.
  tier_center_x <- numeric(K)
  tier_center_y <- numeric(K)
  for (r in seq_len(K)) {
    tier_center_y[r] <- -(r - 1) * tier_spacing
    x_slot <- c(0, 1, -1, -0.45, 0.45)[((r - 1) %% 5) + 1]
    tier_center_x[r] <- 0.95 * x_slot
  }

  # Node positions. By default, draw each tier as a compact group (cloud)
  # rather than as a single horizontal line, so intra-tier structure is visible.
  node_x <- numeric(n)
  node_y <- numeric(n)

  if (layout_style %in% c("fr_ordered", "stress_ordered", "kk_ordered")) {
    base_layout <- switch(
      layout_style,
      fr_ordered = igraph::layout_with_fr(g, weights = E(g)$weight),
      kk_ordered = igraph::layout_with_kk(g),
      stress_ordered = {
        if (requireNamespace("graphlayouts", quietly = TRUE)) {
          graphlayouts::layout_with_stress(g)
        } else {
          message("graphlayouts not installed; using FR fallback for stress_ordered")
          igraph::layout_with_fr(g, weights = E(g)$weight)
        }
      }
    )

    for (r in seq_len(K)) {
      idx <- which(z_ranked == r)
      nb <- length(idx)
      if (nb == 0) next

      y_center <- tier_center_y[r]
      x_center <- tier_center_x[r]

      if (nb == 1) {
        node_x[idx] <- x_center
        node_y[idx] <- y_center
        next
      }

      bx <- base_layout[idx, 1]
      by <- base_layout[idx, 2]
      bx_sd <- stats::sd(bx)
      by_sd <- stats::sd(by)
      bx_scaled <- if (is.finite(bx_sd) && bx_sd > 0) (bx - mean(bx)) / bx_sd else rep(0, nb)
      by_scaled <- if (is.finite(by_sd) && by_sd > 0) (by - mean(by)) / by_sd else rep(0, nb)

      rad_x <- min(0.98, 0.18 * sqrt(nb) + 0.10)
      rad_y <- min(0.56, 0.13 * sqrt(nb) + 0.06)
      node_x[idx] <- x_center + rad_x * tanh(0.8 * bx_scaled)
      node_y[idx] <- y_center + rad_y * tanh(0.8 * by_scaled)
    }
  } else {
  
  for (r in seq_len(K)) {
    idx <- which(z_ranked == r)
    nb <- length(idx)
    if (nb == 0) next

    # Stable within-tier ordering (stronger first, ties by out-degree then id).
    idx_ord <- idx[order(-node_success[idx], -node_outdeg[idx], idx)]
    y_center <- tier_center_y[r]
    x_center <- tier_center_x[r]

    if (layout_style == "tier_cluster") {
      # Compact cloud using a sunflower pattern (fills area, not just perimeter).
      if (nb == 1) {
        node_x[idx_ord] <- x_center
        node_y[idx_ord] <- y_center
      } else {
        phi <- pi * (3 - sqrt(5))
        seq_idx <- seq_len(nb)
        radial <- sqrt((seq_idx - 0.5) / nb)
        theta <- seq_idx * phi + 0.35 * r
        rad_x <- min(0.95, 0.19 * sqrt(nb) + 0.08)
        rad_y <- min(0.50, 0.13 * sqrt(nb) + 0.05)
        node_x[idx_ord] <- x_center + rad_x * radial * cos(theta)
        node_y[idx_ord] <- y_center + rad_y * radial * sin(theta)
      }
    } else if (layout_style == "tier_fan") {
      # Fan each tier into a shallow arc.
      x_seq <- if (nb == 1) 0 else seq(-1, 1, length.out = nb)
      node_x[idx_ord] <- x_center + x_seq * (0.85 + 0.03 * nb)
      node_y[idx_ord] <- y_center + 0.28 * (1 - x_seq^2)
    } else {
      # Legacy line-style tier layout.
      # Keep a tier-line look but cap horizontal spread so tiers do not get
      # visually squashed when using coord_equal().
      if (nb == 1) {
        node_x[idx_ord] <- x_center
        node_y[idx_ord] <- y_center
      } else {
        x_seq <- seq(-1, 1, length.out = nb)
        span_x <- (1.10 + 0.16 * sqrt(nb)) * tier_line_spacing_scale
        wave <- 0.040 * tier_line_spacing_scale * sin(seq_len(nb) * 2 * pi / nb)
        node_x[idx_ord] <- x_center + span_x * x_seq
        node_y[idx_ord] <- y_center + wave
      }
    }
  }
  }
  
  V(g)$x <- node_x
  V(g)$y <- node_y
  
  # Filter edges by threshold
  if (edge_threshold > 0) {
    g <- igraph::delete_edges(g, which(E(g)$weight < edge_threshold))
  }
  
  # Create layout matrix
  layout_mat <- cbind(node_x, node_y)
  
  x_rng <- range(node_x, na.rm = TRUE)
  y_rng <- range(node_y, na.rm = TRUE)
  x_span <- max(diff(x_rng), 1e-6)
  y_span <- max(diff(y_rng), 1e-6)

  # Y position for bottom annotation
  y_bottom <- y_rng[1] - 0.07 * y_span

  x_pad <- 0.025 * x_span
  y_top_pad <- 0.03 * y_span
  y_bottom_pad <- 0.04 * y_span

  x_limits <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
  y_limits <- c(min(y_bottom - y_bottom_pad, y_rng[1] - y_bottom_pad), y_rng[2] + y_top_pad)
  
  # Plot with ggraph
  edge_layer_main <- if (layout_style == "tier_line") {
    geom_edge_arc(aes(filter = direction != "backward", color = I(color), alpha = weight),
                  arrow = arrow(length = unit(2, "mm"), type = "closed"),
                  end_cap = circle(3, "mm"),
                  strength = 0.1)
  } else {
    geom_edge_link(aes(filter = direction != "backward", color = I(color), alpha = weight),
                   arrow = arrow(length = unit(2, "mm"), type = "closed"),
                   end_cap = circle(3, "mm"))
  }

  edge_layer_backward <- if (layout_style == "tier_line") {
    geom_edge_arc(aes(filter = direction == "backward", color = I(color)),
                  alpha = 0.9,
                  arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
                  end_cap = circle(3, "mm"),
                  strength = 0.1)
  } else {
    geom_edge_link(aes(filter = direction == "backward", color = I(color)),
                   alpha = 0.9,
                   arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
                   end_cap = circle(3, "mm"))
  }

  caption_text <- if (!is.null(subtitle_text) && nzchar(subtitle_text)) {
    subtitle_text
  } else if (!is.null(title) && nzchar(title)) {
    title
  } else {
    NULL
  }

  ggraph(g, layout = layout_mat) +
    # Edges
    edge_layer_main +
    edge_layer_backward +
    # Nodes (larger, bolder so the figure does not read zoomed-out)
    geom_node_point(aes(fill = tier),
                    shape = 21, size = node_size_scale * 1.6, stroke = 0.6) +
    geom_node_text(aes(label = name), size = 2.3, repel = TRUE,
             max.overlaps = 20, segment.color = NA) +
    # Backward edge statistics at bottom
    annotate("text", x = mean(range(node_x)), y = y_bottom,
         label = sprintf("Backward edges: %d (zeta[viol] = %.1f%%)",
                 total_backward, pct_backward),
             size = 4, color = "#B2182B", fontface = "bold") +
    scale_fill_viridis_d(
      option = "C",
      name = "Tiered block",
      drop = FALSE,
      guide = guide_legend(override.aes = list(size = 4.8, stroke = 0.45))
    ) +
    scale_edge_alpha(range = c(0.1, 0.8), guide = "none") +
    labs(caption = caption_text) +
    coord_equal(xlim = x_limits, ylim = y_limits, clip = "off") +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.caption = element_text(
        hjust = 0.5,
        size = 11,
        face = "bold",
        color = "gray35",
        margin = margin(t = 2)
      ),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.margin = margin(t = -2),
      legend.box.margin = margin(t = -2),
      plot.margin = margin(1, 1, 2, 1)
    )
}

compute_ordered_stress_layout <- function(
    A,
    z_ranked,
    layout_engine = c("stress", "fr"),
    tier_strength = 0.72,
    tier_spacing = 1.35,
    within_tier_jitter = 0.16,
    seed = 1
) {
  layout_engine <- match.arg(layout_engine)
  set.seed(seed)

  if (inherits(A, "Matrix")) A <- as.matrix(A)
  diag(A) <- 0

  n <- nrow(A)
  K <- max(z_ranked, na.rm = TRUE)

  g_und <- igraph::graph_from_adjacency_matrix(
    A + t(A),
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  if (igraph::ecount(g_und) > 0) {
    E(g_und)$weight <- pmax(E(g_und)$weight, 1e-6)
  }

  if (
    layout_engine == "stress" &&
    requireNamespace("graphlayouts", quietly = TRUE)
  ) {
    xy <- graphlayouts::layout_with_stress(g_und, weights = E(g_und)$weight)
  } else {
    xy <- igraph::layout_with_fr(
      g_und,
      weights = E(g_und)$weight,
      niter = 2000
    )
  }

  xy <- as.data.frame(xy)
  names(xy) <- c("x_raw", "y_raw")

  x_net <- as.numeric(scale(xy$x_raw))
  y_net <- as.numeric(scale(xy$y_raw))

  x_net[!is.finite(x_net)] <- 0
  y_net[!is.finite(y_net)] <- 0

  y_tier <- -tier_spacing * (z_ranked - 1)
  y <- tier_strength * y_tier + (1 - tier_strength) * y_net

  x <- x_net
  for (k in seq_len(K)) {
    idx <- which(z_ranked == k)
    if (length(idx) > 0) {
      x[idx] <- x[idx] - mean(x[idx], na.rm = TRUE)
      x[idx] <- x[idx] + rnorm(1, sd = 0.18)
      y[idx] <- y[idx] + rnorm(length(idx), sd = within_tier_jitter)
    }
  }

  x <- 1.45 * x

  data.frame(
    node = seq_len(n),
    x = x,
    y = y,
    block = z_ranked
  )
}

plot_ordered_network_stress <- function(
    A,
    z_hat,
    title = "Ordered network",
    subtitle_text = NULL,
    alpha = 0.5,
    layout_engine = c("stress", "fr"),
    tier_strength = 0.72,
    tier_spacing = 1.35,
    within_tier_jitter = 0.16,
    edge_threshold = 0,
    node_size_scale = 3.0,
    seed = 1,
    show_tier_guides = FALSE,
    show_node_labels = TRUE,
    precomputed_rate = NULL,
    precomputed_backward = NULL,
    precomputed_cross = NULL
) {
  layout_engine <- match.arg(layout_engine)

  if (inherits(A, "Matrix")) A <- as.matrix(A)
  diag(A) <- 0
  n <- nrow(A)

  if (!is.null(precomputed_rate)) {
    z_ranked <- as.integer(z_hat)
  } else {
    z_ranked <- .rank_partition_by_strength(
      A = A,
      z_vec = z_hat,
      alpha = alpha,
      order_direction = "strong_to_weak"
    )
  }

  K <- max(z_ranked, na.rm = TRUE)

  g <- igraph::graph_from_adjacency_matrix(
    A,
    mode = "directed",
    weighted = TRUE,
    diag = FALSE
  )

  V(g)$block <- z_ranked
  V(g)$tier <- factor(
    z_ranked,
    levels = seq_len(K),
    labels = paste("Tier", seq_len(K))
  )
  V(g)$name <- if (!is.null(rownames(A))) rownames(A) else as.character(seq_len(n))

  edge_df <- igraph::as_data_frame(g, what = "edges")
  node_ids <- V(g)$name

  from_idx <- match(edge_df$from, node_ids)
  to_idx   <- match(edge_df$to, node_ids)

  edge_df$from_block <- z_ranked[from_idx]
  edge_df$to_block   <- z_ranked[to_idx]

  edge_df$is_forward  <- edge_df$from_block < edge_df$to_block
  edge_df$is_backward <- edge_df$from_block > edge_df$to_block
  edge_df$is_within   <- edge_df$from_block == edge_df$to_block

  cross_block_edges <- edge_df %>% filter(!is_within)
  total_cross_block <- sum(cross_block_edges$weight)
  total_backward <- sum(cross_block_edges$weight[cross_block_edges$is_backward])
  pct_backward <- if (total_cross_block > 0) {
    100 * total_backward / total_cross_block
  } else {
    0
  }

  if (!is.null(precomputed_rate)) {
    pct_backward <- 100 * as.numeric(precomputed_rate)
    if (!is.null(precomputed_backward)) total_backward <- as.numeric(precomputed_backward)
    if (!is.null(precomputed_cross)) total_cross_block <- as.numeric(precomputed_cross)
  }

  E(g)$direction <- dplyr::case_when(
    edge_df$is_forward  ~ "forward",
    edge_df$is_backward ~ "backward",
    TRUE                ~ "within"
  )

  E(g)$edge_colour <- dplyr::case_when(
    E(g)$direction == "forward"  ~ "#2C7FB8",
    E(g)$direction == "backward" ~ "#C51B3A",
    TRUE                         ~ "#8C8C8C"
  )

  E(g)$edge_alpha <- dplyr::case_when(
    E(g)$direction == "forward"  ~ 0.22,
    E(g)$direction == "backward" ~ 0.78,
    TRUE                         ~ 0.28
  )

  ew <- E(g)$weight
  if (length(ew) > 0 && diff(range(ew, na.rm = TRUE)) > 0) {
    E(g)$edge_width <- scales::rescale(
      ew,
      to = c(0.18, 0.95),
      from = range(ew, na.rm = TRUE)
    )
  } else {
    E(g)$edge_width <- rep(0.55, igraph::ecount(g))
  }

  if (edge_threshold > 0) {
    g <- igraph::delete_edges(g, which(E(g)$weight < edge_threshold))
  }

  lay <- compute_ordered_stress_layout(
    A = A,
    z_ranked = z_ranked,
    layout_engine = layout_engine,
    tier_strength = tier_strength,
    tier_spacing = tier_spacing,
    within_tier_jitter = within_tier_jitter,
    seed = seed
  )

  layout_mat <- as.matrix(lay[, c("x", "y")])

  tier_df <- data.frame(
    block = seq_len(K),
    y = -tier_spacing * (seq_len(K) - 1)
  )

  x_rng <- range(lay$x, na.rm = TRUE)
  y_rng <- range(lay$y, na.rm = TRUE)

  y_span <- max(diff(y_rng), 1e-6)
  x_span <- max(diff(x_rng), 1e-6)
  y_bottom <- y_rng[1] - 0.07 * y_span

  x_pad <- 0.025 * x_span
  y_top_pad <- 0.03 * y_span
  y_bottom_pad <- 0.04 * y_span

  x_limits <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
  y_limits <- c(min(y_bottom - y_bottom_pad, y_rng[1] - y_bottom_pad), y_rng[2] + y_top_pad)

  caption_text <- if (!is.null(subtitle_text) && nzchar(subtitle_text)) {
    subtitle_text
  } else {
    title
  }

  ggraph(g, layout = layout_mat) +
    {
      if (show_tier_guides) {
        geom_hline(
          data = tier_df,
          aes(yintercept = y),
          inherit.aes = FALSE,
          linewidth = 0.28,
          linetype = "dashed",
          colour = "grey78"
        )
      }
    } +
    geom_edge_fan(
      aes(
        filter = direction != "backward",
        colour = I(edge_colour),
        alpha = I(edge_alpha),
        width = I(edge_width)
      ),
      arrow = arrow(length = unit(1.7, "mm"), type = "closed"),
      end_cap = circle(3.0, "mm"),
      strength = 0.28,
      lineend = "round"
    ) +
    geom_edge_fan(
      aes(
        filter = direction == "backward",
        colour = I(edge_colour),
        alpha = I(edge_alpha),
        width = I(edge_width)
      ),
      arrow = arrow(length = unit(2.0, "mm"), type = "closed"),
      end_cap = circle(3.2, "mm"),
      strength = 0.38,
      lineend = "round"
    ) +
    geom_node_point(
      aes(fill = tier),
      shape = 21,
      size = node_size_scale * 1.35,
      stroke = 0.55,
      colour = "grey12"
    ) +
    {
      if (show_node_labels) {
        geom_node_text(
          aes(label = name),
          size = 2.35,
          repel = TRUE,
          max.overlaps = 40,
          colour = "grey10",
          segment.color = NA
        )
      }
    } +
    annotate(
      "text",
      x = mean(x_rng),
      y = y_bottom,
      label = sprintf(
        "Backward edges: %d (zeta[viol] = %.1f%%)",
        round(total_backward),
        pct_backward
      ),
      size = 4.2,
      colour = "#C51B3A",
      fontface = "bold"
    ) +
    scale_fill_viridis_d(
      option = "C",
      name = "Tiered block",
      drop = FALSE,
      guide = guide_legend(override.aes = list(size = 4.8, stroke = 0.45))
    ) +
    scale_edge_width_identity() +
    scale_edge_alpha_identity() +
    coord_equal(xlim = x_limits, ylim = y_limits, clip = "off") +
    labs(caption = caption_text) +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.caption = element_text(
        hjust = 0.5,
        size = 11,
        face = "bold",
        colour = "grey35",
        margin = margin(t = 2)
      ),
      plot.margin = margin(1, 1, 2, 1),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.margin = margin(t = -2),
      legend.box.margin = margin(t = -2)
    )
}

#' Block-level flow matrix with two encodings per cell
#' 1) Inner square color: directional forward probability (red < 0.5, blue > 0.5)
#' 2) Inner square size: proportional to directed edge count in that cell
#'    scaled from smallest non-zero volume to full tile at maximum volume
#' Labels show count and direction-specific pairwise percentage
plot_block_network <- function(A, z_hat, title = "Block Flow Matrix", compact = FALSE,
                               show_cell_labels = TRUE, show_legend = !compact,
                               subtitle_text = NULL) {
  blocks <- sort(unique(z_hat))
  K <- length(blocks)
  
  # Aggregate edges between blocks
  block_A <- matrix(0, K, K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      idx_i <- which(z_hat == blocks[i])
      idx_j <- which(z_hat == blocks[j])
      block_A[i, j] <- sum(A[idx_i, idx_j])
    }
  }
  
  # Block sizes
  block_sizes <- as.numeric(table(z_hat)[as.character(blocks)])
  
  # Create data frame for heatmap with pairwise percentages
  df <- expand.grid(From = 1:K, To = 1:K) %>%
    mutate(
      count = as.vector(block_A),
      type = case_when(
        From < To ~ "forward",
        From > To ~ "backward",
        TRUE ~ "within"
      )
    )
  
  # Compute pairwise percentages for each (i,j) pair
  # For pair (i,j) with i < j: forward = A[i,j], backward = A[j,i]
  # pct_ij = A[i,j] / (A[i,j] + A[j,i])
  df <- df %>%
    rowwise() %>%
    mutate(
      # For off-diagonal: compute pairwise sum and percentage
      pair_min = min(From, To),
      pair_max = max(From, To),
      # Total for this pair (both directions)
      pair_total = block_A[pair_min, pair_max] + block_A[pair_max, pair_min],
      # Forward share for this pair (same value for both cells in the pair)
      forward_share_pair = if (type == "within" || pair_total == 0) {
        NA_real_
      } else {
        block_A[pair_min, pair_max] / pair_total
      },
      # Percentage for this direction
      pct = if (type == "within" || pair_total == 0) NA_real_ else 100 * count / pair_total,
      prob_dir = if (type == "within" || pair_total == 0) NA_real_ else count / pair_total,
      # Create label: "count (pct%)" for off-diagonal, just count for diagonal
      label = if (type == "within") {
        as.character(count)
      } else if (pair_total == 0) {
        "0"
      } else {
        sprintf("%d\n(%.0f%%)", count, pct)
      }
    ) %>%
    ungroup()
  
  # Compute overall statistics
  cross_block <- df %>% filter(type != "within")
  total_cross <- sum(cross_block$count)
  total_backward <- sum(cross_block$count[cross_block$type == "backward"])
  pct_backward <- if (total_cross > 0) 100 * total_backward / total_cross else 0
  
  max_count <- max(df$count, na.rm = TRUE)
  nonzero_counts <- df$count[df$count > 0]
  min_nonzero <- if (length(nonzero_counts)) min(nonzero_counts) else 0

  # Inner square side length in tile units [0,1]
  # 0 for zero-count cells; for non-zero cells scale min->small square and max->full tile
  min_side <- if (compact) 0.40 else 0.32
  df <- df %>%
    mutate(
      inner_side = dplyr::case_when(
        count <= 0 ~ 0,
        max_count <= min_nonzero ~ 1,
        TRUE ~ scales::rescale(count, to = c(min_side, 1), from = c(min_nonzero, max_count))
      )
    )

  # Label color: white on strong color, black otherwise
  df <- df %>%
    mutate(
      label_color = "white"
    )
  
  # Text size based on compact mode and K so manuscript panels stay legible.
  cell_text_size <- if (compact) {
    max(1.7, 3.1 - 0.11 * K)
  } else {
    max(2.3, 4.0 - 0.10 * K)
  }
  axis_text_size <- if (compact) {
    max(8.2, 10.8 - 0.15 * K)
  } else {
    max(9.2, 12.0 - 0.12 * K)
  }
  subtitle_full <- subtitle_text
  
  # Create the heatmap
  p <- ggplot(df, aes(x = factor(To), y = factor(From, levels = K:1))) +
    geom_tile(fill = "white", color = "grey75", linewidth = 0.8) +
    geom_tile(aes(width = inner_side, height = inner_side, fill = prob_dir),
              color = "grey20", linewidth = 0.2) +
    {
      if (show_cell_labels) {
        geom_text(aes(label = label), size = cell_text_size, fontface = "bold",
                  color = df$label_color,
                  lineheight = 0.85)
      }
    } +
    {
      if (show_cell_labels) {
        annotate("text", x = 1:K, y = 0.3,
                 label = paste0("n=", block_sizes),
                 size = if (compact) 2 else 2.5, color = "gray40")
      }
    } +
    scale_fill_gradientn(
      colours = c("#7F0000", "#B2182B", "#FFFFFF", "#2166AC", "#053061"),
      values = scales::rescale(c(0, 0.46, 0.50, 0.54, 1)),
      limits = c(0, 1), na.value = "#D9D9D9",
      name = "Forward share",
      guide = guide_colorbar(
        title.position = "top",
        barheight = grid::unit(15, "mm"),
        barwidth = grid::unit(2.2, "mm"),
        ticks = FALSE
      )
    ) +
    scale_x_discrete(position = "top") +
    labs(
      x = "To Block (weaker →)",
      y = "From Block (stronger →)",
      title = NULL,
      subtitle = subtitle_full
    ) +
    theme_minimal(base_size = if (compact) 10.5 else 11.5) +
    theme(
      axis.text = element_text(size = axis_text_size, face = "bold"),
      axis.title = element_text(size = axis_text_size + 0.2),
      axis.title.x.top = element_text(margin = margin(b = 5)),
      plot.title = element_blank(),
      plot.subtitle = if (is.null(subtitle_full)) element_blank() else element_text(
        hjust = 0.5,
        face = "bold",
        size = if (compact) max(8.6, 10.0 - 0.10 * K) else 10.0
      ),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = if (compact) 6.8 else 8.2, face = "bold"),
      legend.text = element_text(size = if (compact) 6.1 else 7.1),
      legend.key.height = grid::unit(if (compact) 12 else 16, "mm"),
      legend.key.width = grid::unit(if (compact) 2.2 else 2.8, "mm"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA)
    ) +
    coord_fixed(clip = "off")
  
  p
}

#' Combined block network plot for all three models (WST | SST | DCSBM)
plot_combined_block_networks <- function(A, z_wst, z_sst, z_dcsbm, dataset_name) {
  network_title <- pretty_network_title(dataset_name)
  # Create individual plots in compact mode
  p_wst <- plot_block_network(A, z_wst, title = "WST", compact = TRUE,
                              show_cell_labels = TRUE, show_legend = TRUE,
                              subtitle_text = "WST")
  p_sst <- plot_block_network(A, z_sst, title = "SST", compact = TRUE,
                              show_cell_labels = TRUE, show_legend = TRUE,
                              subtitle_text = "SST")
  p_dc <- plot_block_network(A, z_dcsbm, title = "DC-SBM", compact = TRUE,
                             show_cell_labels = TRUE, show_legend = TRUE,
                             subtitle_text = "DC-SBM")
  
  # Combine horizontally
  p_wst + p_sst + p_dc +
    plot_layout(ncol = 3, guides = "collect") +
    plot_annotation(
      title = paste("Block Networks -", network_title),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    ) &
    theme(legend.position = "right")
}

#' Combined clean block network plot (no cell annotations, one shared legend)
plot_combined_block_networks_clean <- function(A, z_wst, z_sst, z_dcsbm, dataset_name) {
  p_wst <- plot_block_network(A, z_wst, compact = TRUE,
                              show_cell_labels = FALSE, show_legend = TRUE,
                              subtitle_text = "WST") +
    labs(x = NULL, y = NULL)
  p_sst <- plot_block_network(A, z_sst, compact = TRUE,
                              show_cell_labels = FALSE, show_legend = TRUE,
                              subtitle_text = "SST") +
    labs(x = NULL, y = NULL)
  p_dc <- plot_block_network(A, z_dcsbm, compact = TRUE,
                             show_cell_labels = FALSE, show_legend = TRUE,
                             subtitle_text = "DC-SBM") +
    labs(x = NULL, y = NULL)

  p_wst + p_sst + p_dc +
    plot_layout(ncol = 3, guides = "collect") &
    theme(
      legend.position = "right",
      plot.margin = margin(3, 3, 3, 3),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

#' Combined ordered node-network plot (WST | SST | DCSBM)
plot_combined_ordered_networks <- function(A, z_wst, z_sst, z_dcsbm, dataset_name,
                                            viol_wst = NULL, viol_sst = NULL,
                                            viol_dcsbm = NULL) {
  network_title <- pretty_network_title(dataset_name)
  pr <- function(v) if (is.null(v)) NULL else v$rate
  pb <- function(v) if (is.null(v)) NULL else v$backward_mass
  pc <- function(v) if (is.null(v)) NULL else v$cross_mass
  p_wst <- plot_ordered_network(
    A, z_wst,
    node_size_scale = 2.4,
    title = NULL,
    subtitle_text = "WST",
    precomputed_rate = pr(viol_wst),
    precomputed_backward = pb(viol_wst),
    precomputed_cross = pc(viol_wst)
  ) + theme(legend.position = "none")

  p_sst <- plot_ordered_network(
    A, z_sst,
    node_size_scale = 2.4,
    title = NULL,
    subtitle_text = "SST",
    precomputed_rate = pr(viol_sst),
    precomputed_backward = pb(viol_sst),
    precomputed_cross = pc(viol_sst)
  ) + theme(legend.position = "none")

  p_dc <- plot_ordered_network(
    A, z_dcsbm,
    node_size_scale = 2.4,
    title = NULL,
    subtitle_text = "DC-SBM",
    precomputed_rate = pr(viol_dcsbm),
    precomputed_backward = pb(viol_dcsbm),
    precomputed_cross = pc(viol_dcsbm)
  ) + theme(legend.position = "none")

  p_wst + p_sst + p_dc +
    plot_layout(ncol = 3) +
    plot_annotation(
      title = paste("Ordered Node Networks -", network_title),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )
}

rank_nodes_within_blocks <- function(A, z_hat, alpha = 0.5) {
  z_ranked <- .rank_partition_by_strength(
    A = A,
    z_vec = z_hat,
    alpha = alpha,
    order_direction = "strong_to_weak"
  )

  matches <- A + t(A)
  success <- rowSums(A) / pmax(rowSums(matches), 1L)
  outdeg <- rowSums(A)

  node_order <- integer(0)
  block_rows <- list()
  block_id <- 0L

  for (block in seq_len(max(z_ranked, na.rm = TRUE))) {
    idx <- which(z_ranked == block)
    if (!length(idx)) next

    idx <- idx[order(-success[idx], -outdeg[idx], idx)]
    start_idx <- length(node_order) + 1L
    node_order <- c(node_order, idx)
    end_idx <- length(node_order)
    block_id <- block_id + 1L

    block_rows[[block_id]] <- data.frame(
      block = block,
      start = start_idx,
      end = end_idx,
      center = 0.5 * (start_idx + end_idx),
      size = length(idx)
    )
  }

  list(
    node_order = node_order,
    z_ranked = z_ranked,
    block_meta = bind_rows(block_rows)
  )
}

plot_ordered_success_matrix <- function(A, z_hat, dataset_name, model_key,
                                        alpha = 0.5) {
  order_info <- rank_nodes_within_blocks(A, z_hat, alpha = alpha)
  ord <- order_info$node_order
  block_meta <- order_info$block_meta
  z_ranked <- order_info$z_ranked

  A_ord <- A[ord, ord, drop = FALSE]
  matches_ord <- A_ord + t(A_ord)
  success_ord <- ifelse(matches_ord > 0, A_ord / matches_ord, NA_real_)
  diag(success_ord) <- NA_real_

  n <- nrow(A_ord)
  df <- expand.grid(sender = seq_len(n), receiver = seq_len(n)) %>%
    mutate(forward_share = as.vector(success_ord))

  boundary_pos <- block_meta$end[-nrow(block_meta)] + 0.5
  block_breaks <- block_meta$center
  block_labels_x <- paste0("B", block_meta$block)
  block_labels_y <- paste0("B", block_meta$block)

  ggplot(df, aes(x = receiver, y = sender, fill = forward_share)) +
    geom_tile(width = 1, height = 1) +
    geom_vline(xintercept = boundary_pos, linewidth = 0.45, color = "white") +
    geom_hline(yintercept = boundary_pos, linewidth = 0.45, color = "white") +
    scale_fill_gradient2(
      low = "#B2182B",
      mid = "white",
      high = "#2166AC",
      midpoint = 0.5,
      limits = c(0, 1),
      na.value = "#F2F2F2",
      name = "Forward share",
      guide = guide_colorbar(
        title.position = "top",
        barheight = grid::unit(15, "mm"),
        barwidth = grid::unit(2.1, "mm"),
        ticks = FALSE
      )
    ) +
    scale_x_continuous(
      position = "top",
      breaks = block_breaks,
      labels = block_labels_x,
      expand = c(0, 0)
    ) +
    scale_y_reverse(
      breaks = block_breaks,
      labels = block_labels_y,
      expand = c(0, 0)
    ) +
    coord_fixed() +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.ticks = element_blank(),
      axis.text.x.top = element_text(size = 6.0, face = "bold", lineheight = 0.9, margin = margin(b = 2)),
      axis.text.y = element_text(size = 6.3, face = "bold", margin = margin(r = 2)),
      axis.title = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 6.8, face = "bold"),
      legend.text = element_text(size = 6.1),
      legend.key.height = grid::unit(12, "mm"),
      legend.key.width = grid::unit(2.1, "mm"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(6, 4, 4, 4)
    )
}

compute_distance_profile_df <- function(A, z_hat, dataset_name, model_key,
                                        alpha = 0.5) {
  z_ranked <- .rank_partition_by_strength(
    A = A,
    z_vec = z_hat,
    alpha = alpha,
    order_direction = "strong_to_weak"
  )
  K <- max(z_ranked, na.rm = TRUE)
  rows <- list()
  row_id <- 0L

  for (k in seq_len(K - 1L)) {
    idx_k <- which(z_ranked == k)
    if (!length(idx_k)) next

    for (ell in seq((k + 1L), K)) {
      idx_ell <- which(z_ranked == ell)
      if (!length(idx_ell)) next

      forward <- sum(A[idx_k, idx_ell, drop = FALSE])
      backward <- sum(A[idx_ell, idx_k, drop = FALSE])
      total_mass <- forward + backward
      if (total_mass <= 0) next

      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        dataset = pretty_network_title(dataset_name),
        model = pretty_model_label(model_key),
        distance = ell - k,
        forward_share = forward / total_mass,
        total_mass = total_mass,
        pair = sprintf("B%d-B%d", k, ell)
      )
    }
  }

  bind_rows(rows)
}

plot_animal_hierarchy_profile <- function(results_list) {
  animal_specs <- data.frame(
    dataset_name = c("mountain_goats", "macaques_data"),
    model_key = c("SST", "DCSBM"),
    stringsAsFactors = FALSE
  )

  profile_rows <- list()
  for (i in seq_len(nrow(animal_specs))) {
    ds <- animal_specs$dataset_name[i]
    model_key <- animal_specs$model_key[i]
    res <- results_list[[ds]]
    if (is.null(res)) next

    z_hat <- switch(
      model_key,
      WST = res$z_hat_wst,
      SST = res$z_hat_sst,
      DCSBM = res$z_hat_dcsbm,
      NULL
    )
    if (is.null(z_hat) || is.null(res$A)) next

    profile_rows[[length(profile_rows) + 1L]] <- compute_distance_profile_df(
      A = res$A,
      z_hat = z_hat,
      dataset_name = ds,
      model_key = model_key
    )
  }

  profile_df <- bind_rows(profile_rows)
  if (!nrow(profile_df)) return(NULL)

  summary_df <- profile_df %>%
    group_by(dataset, model, distance) %>%
    summarise(
      mean_forward_share = weighted.mean(forward_share, w = total_mass),
      total_mass = sum(total_mass),
      .groups = "drop"
    )

  ggplot(profile_df, aes(x = distance, y = forward_share)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
    geom_point(
      aes(size = total_mass),
      color = "#4C78A8",
      alpha = 0.40,
      position = position_jitter(width = 0.08, height = 0)
    ) +
    geom_line(
      data = summary_df,
      aes(y = mean_forward_share, group = 1),
      color = "#B2182B",
      linewidth = 0.9
    ) +
    geom_point(
      data = summary_df,
      aes(y = mean_forward_share),
      color = "#B2182B",
      size = 2.2
    ) +
    facet_grid(
      . ~ dataset,
      scales = "free_x",
      space = "free_x",
      labeller = labeller(dataset = function(x) sub(" Network$", "", x))
    ) +
    scale_x_continuous(breaks = sort(unique(profile_df$distance))) +
    scale_y_continuous(
      limits = c(0.5, 1.0),
      breaks = seq(0.5, 1.0, by = 0.1),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_size_continuous(
      range = c(1.8, 5.8),
      name = "Block-pair mass",
      guide = guide_legend(
        title.position = "top",
        override.aes = list(alpha = 0.40, color = "#4C78A8")
      )
    ) +
    labs(
      x = "Block distance",
      y = "Empirical forward share"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      strip.background = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 6.8, face = "bold"),
      legend.text = element_text(size = 6.1),
      legend.key.height = grid::unit(10, "mm"),
      legend.key.width = grid::unit(2.5, "mm"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA)
    )
}

plot_best_model_ordered_matrices <- function(results_list) {
  dataset_order <- c("citations_data", "high_school")
  plot_list <- list()

  for (ds in dataset_order) {
    res <- results_list[[ds]]
    if (is.null(res) || is.null(res$A)) next

    model_key <- application_best_model(ds)
    z_hat <- switch(
      model_key,
      WST = res$z_hat_wst,
      SST = res$z_hat_sst,
      DCSBM = res$z_hat_dcsbm,
      NULL
    )
    if (is.null(z_hat)) next

    plot_list[[ds]] <- plot_ordered_success_matrix(
      A = res$A,
      z_hat = z_hat,
      dataset_name = ds,
      model_key = model_key
    )
  }

  if (!length(plot_list)) return(NULL)

  wrap_plots(plot_list, ncol = length(plot_list), guides = "collect") +
    plot_annotation(theme = theme(plot.title = element_blank(), plot.subtitle = element_blank())) &
    theme(
      legend.position = "right",
      plot.margin = margin(3, 3, 3, 3),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# =============================================================================
# SANITY CHECK FUNCTIONS: Violations, Orderings, Min-Feedback, Kendall's tau
# =============================================================================

# --- Compute empirical rho from A and z_hat (total block mass based) ---
compute_empirical_rho <- function(A, z_hat, alpha = 0.5) {
  blocks <- sort(unique(z_hat))
  K <- length(blocks)
  C_blk <- matrix(0, K, K)
  M_blk <- matrix(0, K, K)
  n <- nrow(A)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) {
        ki <- match(z_hat[i], blocks)
        kj <- match(z_hat[j], blocks)
        C_blk[ki, kj] <- C_blk[ki, kj] + A[i, j]
        M_blk[ki, kj] <- M_blk[ki, kj] + A[i, j] + A[j, i]
      }
    }
  }
  rho <- (C_blk + alpha) / (pmax(M_blk, 0) + 2 * alpha)
  diag(rho) <- 0.5
  list(rho = rho, C_blk = C_blk, M_blk = M_blk, blocks = blocks)
}

# --- Compute violation stats from draws (as in application.R diagnostics) ---
# WST/SST use intrinsic labels. DC-SBM ranks blocks by empirical strength.
compute_violation_table_from_draws <- function(fit, A, regime,
                                                alpha = 0.5,
                                                method_order = "bt") {
  # Extract z chain
  if (is.matrix(fit$z)) {
    z_chain <- fit$z
  } else if (is.list(fit$z)) {
    z_chain <- do.call(rbind, fit$z)
  } else {
    stop("fit$z must be a matrix or list")
  }
  S <- nrow(z_chain)
  n <- ncol(z_chain)

  violation_rate  <- rep(NA_real_, S)
  violation_count <- rep(NA_real_, S)
  cross_mass      <- rep(NA_real_, S)

  for (s in seq_len(S)) {
    z_s <- z_chain[s, ]
    z_ranked <- if (regime %in% c("WST", "SST")) {
      z_s
    } else {
      .rank_partition_by_strength(
        A = A, z_vec = z_s, alpha = alpha,
        order_direction = "strong_to_weak"
      )
    }
    vstats <- .violation_stats_from_ranked(A = A, z_ranked = z_ranked)
    violation_rate[s]  <- vstats$rate
    violation_count[s] <- vstats$count
    cross_mass[s]      <- vstats$cross_mass
  }

  # Summarize
  vr <- violation_rate[!is.na(violation_rate)]
  vc <- violation_count[!is.na(violation_count)]
  cm <- cross_mass[!is.na(cross_mass)]

  data.frame(
    violation_rate_mean  = mean(vr),
    violation_rate_lo    = if (length(vr) > 1) quantile(vr, 0.025) else NA_real_,
    violation_rate_hi    = if (length(vr) > 1) quantile(vr, 0.975) else NA_real_,
    violation_count_mean = mean(vc),
    cross_mass_mean      = mean(cm),
    n_draws              = length(vr),
    stringsAsFactors     = FALSE
  )
}

# --- Extract block ordering for a model from z_hat + A ---
# Returns the ordering (strongest-to-weakest) as block labels and node names
get_block_ordering <- function(z_hat, A, alpha = 0.5, method = "bt") {
  emp <- compute_empirical_rho(A, z_hat, alpha)
  ord_res <- order_from_rho(emp$rho, method = method,
                            order_direction = "strong_to_weak")
  blocks <- emp$blocks
  # ord_res$order gives positions: ord[1] is the index of the strongest block
  # Map back to original block labels
  ordered_blocks <- blocks[ord_res$order]
  list(
    ordered_blocks = ordered_blocks,
    scores         = ord_res$scores,
    block_map      = emp$blocks,
    rho            = emp$rho,
    order_index    = ord_res$order
  )
}

# --- Block-level Minimum Feedback Arc Set (MFAS) ordering ---
compute_min_feedback_ordering <- function(A, z_hat, alpha = 0.5, max_passes = 20L) {
  blocks <- sort(unique(z_hat))
  K <- length(blocks)
  if (K < 2) {
    return(list(ordered_blocks = blocks, backward_mass = 0, total_cross_mass = 0))
  }

  W <- matrix(0, K, K)
  n <- nrow(A)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j && z_hat[i] != z_hat[j]) {
        ki <- match(z_hat[i], blocks)
        kj <- match(z_hat[j], blocks)
        W[ki, kj] <- W[ki, kj] + A[i, j]
      }
    }
  }

  total_cross <- sum(W)
  net_out <- rowSums(W) - colSums(W)
  ord <- order(-net_out, seq_len(K))

  backward_mass_fn <- function(o) {
    bm <- 0
    for (a in seq_len(K - 1)) {
      for (b in (a + 1):K) {
        bm <- bm + W[o[b], o[a]]
      }
    }
    bm
  }

  improved <- TRUE
  passes <- 0L
  while (improved && passes < max_passes) {
    improved <- FALSE
    passes <- passes + 1L
    for (i in seq_len(K - 1)) {
      o_new <- ord
      o_new[c(i, i + 1)] <- o_new[c(i + 1, i)]
      if (backward_mass_fn(o_new) + 1e-12 < backward_mass_fn(ord)) {
        ord <- o_new
        improved <- TRUE
      }
    }
  }

  final_bm <- backward_mass_fn(ord)
  ordered_blocks <- blocks[ord]

  list(
    ordered_blocks = ordered_blocks,
    backward_mass  = final_bm,
    total_cross_mass = total_cross,
    backward_rate  = if (total_cross > 0) final_bm / total_cross else NA_real_,
    order_index    = ord,
    W              = W
  )
}

# =============================================================================
# NODE-LEVEL MFAS: ordering of individual nodes minimizing backward flow
# =============================================================================

#' Compute node-level MFAS ordering.
#' Finds the permutation sigma of {1,...,n} minimizing:
#'   B(sigma) = sum_{a < b} A[sigma(b), sigma(a)]
#' i.e. total weight of edges going from a node ranked weaker to one ranked stronger.
#'
#' Uses greedy sort by net-outflow + iterative adjacent-swap improvement.
#' For n up to ~100, this converges quickly.
compute_node_mfas <- function(A, max_passes = 30L) {
  n <- nrow(A)
  if (n < 2) return(list(ordering = seq_len(n), backward_mass = 0,
                          total_mass = sum(A), backward_rate = 0))

  # Total directed mass (excluding diagonal)
  diag(A) <- 0
  total_mass <- sum(A)

  # Greedy: sort nodes by net outflow  (rowSums - colSums)
  net_out <- rowSums(A) - colSums(A)
  ord <- order(-net_out, seq_len(n))

  # Backward mass: sum of A[ord[b], ord[a]] for a < b
  # Efficient incremental for adjacent swaps:
  # swapping positions i and i+1 only changes the contribution of that pair.
  # Before swap: backward from pair = A[ord[i+1], ord[i]]
  # After  swap: backward from pair = A[ord[i], ord[i+1]]
  # Everything else stays the same.

  # Full backward mass (needed once to initialize)
  compute_full_bm <- function(o) {
    bm <- 0
    for (a in seq_len(n - 1L)) {
      for (b in (a + 1L):n) {
        bm <- bm + A[o[b], o[a]]
      }
    }
    bm
  }

  current_bm <- compute_full_bm(ord)

  improved <- TRUE
  passes <- 0L
  while (improved && passes < max_passes) {
    improved <- FALSE
    passes <- passes + 1L
    for (i in seq_len(n - 1L)) {
      u <- ord[i]
      v <- ord[i + 1L]
      # Currently: u at position i, v at position i+1
      # Backward contribution from this pair: A[v, u] (v is below u, so edge v->u is backward)
      # After swap: u at position i+1, v at position i
      # Backward contribution from this pair: A[u, v]
      # Also need to account for interactions of u,v with all other nodes.

      # Full delta: for all other nodes w at position p:
      #   if p < i:   pair (p, i) and (p, i+1)
      #     before: bw += A[u, w_p] (for pos i) + A[v, w_p] (for pos i+1)  -- only backward part
      #     Hmm this is getting complicated. For correctness with O(n) per swap:

      # delta = sum over all other positions p != i, i+1:
      #   if p < i:  (A[v, ord[p]] - A[u, ord[p]]) + (A[u, ord[p]] - A[v, ord[p]])  = 0 ??
      # No, let me think carefully.
      # backward_mass = sum_{a<b} A[ord[b], ord[a]]
      # Swapping positions i and i+1 (values u <-> v):
      # Changes only pairs involving position i or i+1.

      # For pair (i, i+1) itself:
      #   before: A[v, u], after: A[u, v]  -> delta1 = A[u,v] - A[v,u]

      # For each position p < i, pair (p, i) and (p, i+1):
      #   before: A[u, ord[p]] + A[v, ord[p]]  (ord[i]=u is "below" ord[p], and ord[i+1]=v is below both)
      #   after:  A[v, ord[p]] + A[u, ord[p]]  -> same! delta = 0.

      # For each position p > i+1, pair (i, p) and (i+1, p):
      #   before: A[ord[p], u] + A[ord[p], v]
      #   after:  A[ord[p], v] + A[ord[p], u] -> same! delta = 0.

      # So the only change is from the pair (i, i+1):
      delta <- A[u, v] - A[v, u]
      if (delta < -1e-12) {
        # Swap improves (reduces backward mass)
        ord[i]     <- v
        ord[i + 1L] <- u
        current_bm <- current_bm + delta
        improved <- TRUE
      }
    }
  }

  node_names <- if (!is.null(rownames(A))) rownames(A) else as.character(seq_len(n))
  list(
    ordering      = ord,               # node indices, strongest first
    node_names    = node_names[ord],    # names in order
    backward_mass = current_bm,
    total_mass    = total_mass,
    backward_rate = if (total_mass > 0) current_bm / total_mass else NA_real_,
    n_passes      = passes
  )
}

# =============================================================================
# NODE-LEVEL VIOLATION RATE
# =============================================================================

#' Compute violation rate for a given node ordering.
#' ordering: integer vector of node indices (strongest first).
#' z_hat:   optional block assignments. If supplied, only cross-block edges
#'          are counted (matching .violation_stats_from_ranked in
#'          transitivity_check_helper.R). If NULL, all off-diagonal edges count.
#' Returns: backward_mass / total_mass plus counts.
node_violation_rate <- function(A, ordering, z_hat = NULL) {
  n <- length(ordering)
  # Build rank vector: rank[node] = position in ordering (1 = strongest)
  rank_vec <- integer(nrow(A))
  rank_vec[ordering] <- seq_len(n)

  total <- 0
  backward <- 0
  for (i in seq_len(nrow(A))) {
    for (j in seq_len(ncol(A))) {
      if (i == j) next
      if (A[i, j] <= 0) next
      # If z_hat supplied, skip within-block edges
      if (!is.null(z_hat) && z_hat[i] == z_hat[j]) next
      total <- total + A[i, j]
      # Edge i->j is backward if rank[i] > rank[j] (i is weaker than j)
      if (rank_vec[i] > rank_vec[j]) {
        backward <- backward + A[i, j]
      }
    }
  }

  list(
    backward_mass = backward,
    total_mass    = total,
    violation_rate = if (total > 0) backward / total else NA_real_
  )
}

# =============================================================================
# NODE-LEVEL WST / SST CONFORMITY
# =============================================================================

#' Compute WST and SST conformity at node level for a given ordering.
#' Uses smoothed pairwise dominance probabilities.
#' ordering: integer vector of node indices (strongest first).
#' T_triples: number of triples to sample (NULL = use all).
node_wst_sst_conformity <- function(A, ordering, alpha = 0.5,
                                     T_triples = 5000L, seed = 42L) {
  n <- length(ordering)
  if (n < 3) {
    return(list(thetaW = NA_real_, thetaS = NA_real_,
                prem = 0L, total = 0L, coverage = 0,
                thetaW_all = NA_real_, thetaS_all = NA_real_))
  }

  # Smoothed pairwise dominance: P[i,j] = (A[i,j] + alpha) / (A[i,j] + A[j,i] + 2*alpha)
  M <- A + t(A)
  P <- (A + alpha) / (M + 2 * alpha)
  diag(P) <- 0.5

  # Ordered triples: (ordering[a], ordering[b], ordering[c]) with a < b < c
  # Here a < b < c means a is strongest, c is weakest.
  all_n <- choose(n, 3)
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(T_triples) && T_triples < all_n) {
    # Sample triples of ordered positions
    triples <- matrix(NA_integer_, nrow = 3, ncol = T_triples)
    for (t in seq_len(T_triples)) {
      pos <- sort(sample.int(n, 3))
      triples[, t] <- ordering[pos]
    }
  } else {
    # Use all triples
    pos_triples <- utils::combn(n, 3)  # 3 x C(n,3), positions in ordering
    T_triples <- ncol(pos_triples)
    triples <- matrix(NA_integer_, nrow = 3, ncol = T_triples)
    for (t in seq_len(T_triples)) {
      triples[, t] <- ordering[pos_triples[, t]]
    }
  }

  # Score using the same logic as .score_triples
  # triples[1,t] = strongest, triples[2,t] = middle, triples[3,t] = weakest
  prem <- okW_prem <- okS_prem <- 0L
  okW_all <- okS_all <- 0L
  total <- ncol(triples)

  for (tt in seq_len(total)) {
    a <- triples[1, tt]; b <- triples[2, tt]; d <- triples[3, tt]
    rab <- P[a, b]; rbd <- P[b, d]; rad <- P[a, d]

    if (rad >= 0.5)             okW_all <- okW_all + 1L
    if (rad >= max(rab, rbd))   okS_all <- okS_all + 1L

    if (rab >= 0.5 && rbd >= 0.5) {
      prem <- prem + 1L
      if (rad >= 0.5)           okW_prem <- okW_prem + 1L
      if (rad >= max(rab, rbd)) okS_prem <- okS_prem + 1L
    }
  }

  list(
    thetaW     = if (prem > 0) okW_prem / prem else NA_real_,
    thetaS     = if (prem > 0) okS_prem / prem else NA_real_,
    prem       = prem,
    total      = total,
    coverage   = if (total > 0) prem / total else 0,
    thetaW_all = if (total > 0) okW_all / total else NA_real_,
    thetaS_all = if (total > 0) okS_all / total else NA_real_
  )
}

# =============================================================================
# DRAWS-BASED Z_HAT AND NODE ORDERING
# =============================================================================

#' Compact non-contiguous integer labels to contiguous 1-based.
#' E.g. {1,3,4,7} -> {1,2,3,4}.
.compact_labels_row <- function(z_row) {
  z_int <- as.integer(z_row)
  lev   <- sort(unique(z_int))
  if (length(lev) == 0L) return(z_int)
  # If already contiguous 1-based, return as-is (fast path)
  if (lev[1] == 1L && lev[length(lev)] == length(lev)) return(z_int)
  map <- setNames(seq_along(lev), lev)
  as.integer(map[as.character(z_int)])
}

#' Obtain z_hat from posterior draws.
#'
#' For WST/SST the draws are NOT relabelled (the ordered prior makes block
#' labels identifiable: block 1 = strongest).
#' For DCSBM every draw is first compacted to contiguous 1-based labels
#' (the MCMC can leave gaps when blocks are emptied), then relabelled so
#' that the block with the highest average per-node success probability
#' receives label 1, etc.
#'
#' z_hat is obtained via PSM / minVI(method="draws") on the (possibly
#' relabelled) chain.  No further relabelling of z_hat is performed.
#'
#' @param fit  Fitted model object (must contain \code{$z}).
#' @param A    Adjacency matrix.
#' @param model One of "WST", "SST", "DCSBM".
#' @return A list with components:
#'   \item{z_hat}{Integer vector length n – point-estimate partition.}
#'   \item{z_chain}{S x n integer matrix of (relabelled) draws used for
#'                  posterior probability computation.}
#'   \item{K}{Number of blocks in z_hat.}
get_z_hat_from_draws <- function(fit, A, model = c("WST", "SST", "DCSBM")) {
  model <- match.arg(model)
  if (inherits(A, "Matrix")) A <- as.matrix(A)

  z_chain <- z_to_matrix(fit$z)          # S x n
  S <- nrow(z_chain)
  n <- ncol(z_chain)

  # ---- DCSBM: compact + relabel every draw by block score ----
  if (model == "DCSBM") {
    N <- A + t(A)
    for (s in seq_len(S)) {
      # Compact to contiguous 1-based (DCSBM chains can have gaps)
      z_s <- .compact_labels_row(z_chain[s, ])
      K_s <- max(z_s, na.rm = TRUE)
      # Relabel by average per-node success: strongest block -> 1
      r_bar <- compute_block_scores(z_s, A, N)
      perm  <- order(r_bar, decreasing = TRUE, na.last = TRUE)
      inv_perm <- integer(K_s)
      inv_perm[perm] <- seq_len(K_s)
      z_chain[s, ] <- inv_perm[z_s]
    }
  }
  # WST / SST: draws kept as-is (ordered prior ensures label identity;
  #            MCMC already produces contiguous 1-based labels)

  # ---- Filter invalid rows ----
  valid_rows <- apply(z_chain, 1,
    function(r) all(is.finite(r)) && all(r >= 1))
  z_valid <- z_chain[valid_rows, , drop = FALSE]
  if (nrow(z_valid) == 0L) stop("No valid Z draws available.")

  # ---- PSM / minVI for z_hat ----
  psm   <- mcclust::comp.psm(z_valid)
  z_hat <- mcclust.ext::minVI(psm, cls.draw = z_valid, method = "draws")$cl
  z_hat <- as.integer(if (is.matrix(z_hat)) z_hat[1, ] else z_hat)

  list(
    z_hat   = z_hat,
    z_chain = z_valid,
    K       = max(z_hat, na.rm = TRUE)
  )
}

# =============================================================================
# POSTERIOR NODE ORDERING
# =============================================================================

#' Rank nodes using posterior block-membership probabilities.
#'
#' For each node i the posterior probability P(z_i = k) is estimated as the
#' fraction of MCMC draws in which node i is assigned to block k.
#' Nodes are then sorted lexicographically in descending order:
#'   first by P(z_i = 1)  (probability of being in the top block),
#'   breaking ties by P(z_i = 2),  then P(z_i = 3),  etc.
#'
#' @param z_chain S x n integer matrix of (relabelled) draws.
#'                Labels must satisfy: 1 = strongest block.
#' @return Integer vector of length n – node indices ordered strongest first.
posterior_node_ordering <- function(z_chain) {
  S     <- nrow(z_chain)
  n     <- ncol(z_chain)
  K_max <- max(z_chain, na.rm = TRUE)

  # P(z_i = k) for each node i and block k
  prob_mat <- matrix(0, nrow = n, ncol = K_max)
  for (k in seq_len(K_max)) {
    prob_mat[, k] <- colMeans(z_chain == k)
  }

  # Sort nodes lexicographically by (P(z_i=1), P(z_i=2), ...) descending.
  # Negate so that order() gives descending.
  sort_df <- as.data.frame(-prob_mat)
  names(sort_df) <- paste0("V", seq_len(K_max))
  ordering <- do.call(order, sort_df)
  ordering
}

# =============================================================================
# MODEL-INDUCED NODE ORDERING  (legacy wrapper – now delegates to
#                                posterior_node_ordering when z_chain supplied)
# =============================================================================

#' Extract the node ordering induced by a model's partition.
#'
#' If a z_chain is supplied (from get_z_hat_from_draws) the ordering is
#' computed via posterior block-membership probabilities.
#' Otherwise falls back to ranking blocks via .rank_partition_by_strength
#' and sorting nodes within blocks by smoothed mean success.
#'
#' @param z_hat   Point-estimate partition (integer vector).
#' @param A       Adjacency matrix.
#' @param alpha   Smoothing parameter (fallback mode only).
#' @param z_chain Optional S x n matrix of (relabelled) draws.
#' @return Integer vector of length n – node indices ordered strongest first.
model_node_ordering <- function(z_hat, A, alpha = 0.5, z_chain = NULL) {
  # ---- Posterior-probability mode ----
  if (!is.null(z_chain)) {
    return(posterior_node_ordering(z_chain))
  }

  # ---- Fallback: deterministic block-sort mode ----
  n <- nrow(A)
  z_ranked <- .rank_partition_by_strength(
    A = A, z_vec = z_hat, alpha = alpha,
    order_direction = "strong_to_weak"
  )
  K <- max(z_ranked, na.rm = TRUE)
  M <- A + t(A)
  P <- (A + alpha) / (M + 2 * alpha)
  diag(P) <- NA
  mean_success <- rowMeans(P, na.rm = TRUE)

  ordering <- integer(0)
  for (r in seq_len(K)) {
    idx <- which(z_ranked == r)
    ordering <- c(ordering, idx[order(-mean_success[idx])])
  }
  ordering
}

# --- Kendall's tau distance between two orderings ---
# Both orderings are vectors of labels (strongest to weakest).
# Kendall's tau distance = number of pairwise disagreements / total pairs.
kendall_tau_distance <- function(ord1, ord2) {
  # Both must be permutations of the same set
  if (length(ord1) != length(ord2)) {
    # If different K, we intersect
    common <- intersect(ord1, ord2)
    if (length(common) < 2) return(NA_real_)
    ord1 <- ord1[ord1 %in% common]
    ord2 <- ord2[ord2 %in% common]
  }
  K <- length(ord1)
  if (K < 2) return(0)

  # Convert ord2 to rank form relative to ord1
  rank1 <- match(ord1, ord1)  # identity: 1, 2, ..., K
  rank2 <- match(ord1, ord2)  # position of each element of ord1 in ord2

  # Count inversions in rank2
  n_disc <- 0L
  n_pairs <- choose(K, 2)
  for (i in seq_len(K - 1)) {
    for (j in (i + 1):K) {
      if (rank2[i] > rank2[j]) n_disc <- n_disc + 1L
    }
  }

  n_disc / n_pairs
}

# --- Get node-level ordering within each block ---
# Within each block, nodes are ordered by their empirical success rate
get_item_ordering <- function(z_hat, A, block_order) {
  n <- nrow(A)
  outdeg <- rowSums(A)
  matches <- A + t(A)
  success <- outdeg / pmax(rowSums(matches), 1)

  items_df <- data.frame(
    node = seq_len(n),
    block = z_hat,
    success = success,
    stringsAsFactors = FALSE
  )

  # If A has row/col names, use them
  if (!is.null(rownames(A))) {
    items_df$node_name <- rownames(A)
  } else {
    items_df$node_name <- as.character(seq_len(n))
  }

  # Assign block rank according to block_order
  items_df$block_rank <- match(items_df$block, block_order)

  # Sort: by block rank, then by descending success within block
  items_df <- items_df[order(items_df$block_rank, -items_df$success), ]
  items_df$overall_rank <- seq_len(nrow(items_df))

  items_df
}

# =============================================================================
# BUMP CHART: Node rankings across models + MFAS
# =============================================================================

#' Plot a bump chart showing how each node's rank changes across model orderings.
#'
#' @param A              adjacency matrix
#' @param z_hats         named list of z_hat vectors (e.g. WST, SST, DCSBM)
#' @param node_mfas      output of compute_node_mfas()
#' @param z_chains       named list of z_chain matrices (from get_z_hat_from_draws)
#' @param alpha          smoothing for mean-success tie-breaking (fallback)
#' @param dataset_name   used in the plot title
#' @return ggplot object
plot_bump_chart <- function(A, z_hats, node_mfas, z_chains = NULL,
                            alpha = 0.5, dataset_name = "") {
  n <- nrow(A)
  node_names <- if (!is.null(rownames(A))) rownames(A) else as.character(seq_len(n))

  # Gather orderings: MFAS + each model
  orderings <- list(MFAS = node_mfas$ordering)
  for (mod in names(z_hats)) {
    zc <- if (!is.null(z_chains)) z_chains[[mod]] else NULL
    orderings[[mod]] <- model_node_ordering(z_hats[[mod]], A, alpha = alpha,
                                            z_chain = zc)
  }

  # Build a data.frame: node x model -> rank
  rows <- list()
  for (mod_name in names(orderings)) {
    ord <- orderings[[mod_name]]
    rank_vec <- integer(n)
    rank_vec[ord] <- seq_len(n)
    rows[[mod_name]] <- data.frame(
      node      = seq_len(n),
      node_name = node_names,
      model     = mod_name,
      rank      = rank_vec,
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, rows)
  rownames(df) <- NULL

  # Determine model column order (MFAS first, then models alphabetically)
  model_levels <- c("MFAS", sort(setdiff(names(orderings), "MFAS")))
  df$model <- factor(df$model, levels = model_levels)

  # Colour each node by its MFAS rank (gradient from strong to weak)
  mfas_rank <- integer(n)
  mfas_rank[node_mfas$ordering] <- seq_len(n)
  df$mfas_rank <- mfas_rank[df$node]

  # For labelling: only label on the leftmost and rightmost columns
  df$label_left  <- ifelse(df$model == model_levels[1], df$node_name, NA_character_)
  df$label_right <- ifelse(df$model == model_levels[length(model_levels)], df$node_name, NA_character_)

  # Numeric x for geom_line
  df$x <- as.numeric(df$model)

  p <- ggplot(df, aes(x = x, y = rank, group = node, colour = mfas_rank)) +
    geom_line(alpha = 0.45, linewidth = 0.4) +
    geom_point(size = 1.2) +
    # Left labels
    geom_text(aes(label = label_left), hjust = 1, nudge_x = -0.08,
              size = 1.8, na.rm = TRUE, show.legend = FALSE) +
    # Right labels
    geom_text(aes(label = label_right), hjust = 0, nudge_x = 0.08,
              size = 1.8, na.rm = TRUE, show.legend = FALSE) +
    scale_x_continuous(
      breaks = seq_along(model_levels),
      labels = model_levels,
      expand = expansion(mult = c(0.20, 0.20))
    ) +
    scale_y_reverse(
      breaks = c(1, seq(5, n, by = 5), n),
      name   = "Rank (1 = strongest)"
    ) +
    scale_colour_viridis_c(option = "C", direction = -1,
                           name = "MFAS rank") +
    labs(
      title    = paste("Node Ranking Comparison \u2014", dataset_name),
      subtitle = "Lines connect the same node across orderings",
      x = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "gray40"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position = "right"
    )

  p
}

# =============================================================================
# MAIN FUNCTION: Generate All Visualizations for a Dataset
# =============================================================================

generate_osbm_visualizations <- function(dataset_name, 
                                          ppc_dir = "output/application/ppc",
                                          output_dir = NULL,
                                          dataset_plots_dir = "output/application/plots") {
  cat("Generating visualizations for:", dataset_name, "\n")

  # Canonical violation stats (same formula as build_post_processing.R cube
  # and helper_folder/transitivity_check_helper.R). Pre-computing them here
  # and passing the result into plot_ordered_network prevents that function
  # from re-ranking labels (which would silently produce a different
  # backward-mass than the hierarchy diagnostics).
  .canonical_viol <- function(A, z) {
    if (is.null(z)) return(NULL)
    if (inherits(A, "Matrix")) A <- as.matrix(A)
    rows <- row(A); cols <- col(A); vals <- A
    z_i  <- z[rows]; z_j <- z[cols]
    pos  <- vals > 0
    fwd  <- pos & (z_i < z_j); bwd <- pos & (z_i > z_j)
    cm   <- sum(vals[fwd | bwd]); bm <- sum(vals[bwd])
    list(rate = if (cm > 0) bm / cm else 0,
         backward_mass = bm, cross_mass = cm)
  }
  
  # Create output directories
  # New layout: output/application/plots/<dataset>/osbm_visualizations/
  # (The legacy duplication into output/application/plots/<dataset>/ at the
  # end of this function has been removed.)
  ds_main_dir <- file.path(dataset_plots_dir, dataset_name)
  if (is.null(output_dir) || !nzchar(output_dir)) {
    out_path <- file.path(ds_main_dir, "osbm_visualizations")
  } else {
    out_path <- file.path(output_dir, dataset_name)
  }
  fs::dir_create(out_path)
  fs::dir_create(ds_main_dir)
  
  # Load adjacency matrix: ppc rds first, then choose_dataset_local fallback
  a_obs_path <- file.path(ppc_dir, dataset_name, paste0(dataset_name, "_A_obs.rds"))
  if (file.exists(a_obs_path)) {
    A <- readRDS(a_obs_path)
  } else {
    cat("  A_obs.rds not found in ppc dir, loading from data sources...\n")
    A <- tryCatch(
      choose_dataset_local(dataset_name),
      error = function(e) {
        cat("  ERROR loading dataset:", conditionMessage(e), "\n")
        return(NULL)
      }
    )
    if (is.null(A)) return(invisible(NULL))
  }
  
  # Resolve fit file paths (ppc structure, then latest run dir fallback)
  wst_fit_path  <- find_fit_file_local(dataset_name, "WST",   ppc_dir)
  sst_fit_path  <- find_fit_file_local(dataset_name, "SST",   ppc_dir)
  
  # --- WST Analysis ---
  z_chain_wst <- NULL
  if (!is.null(wst_fit_path) && file.exists(wst_fit_path)) {
    cat("  Processing WST fit...\n")
    wst_fit <- readRDS(wst_fit_path)
    
    # Get z_hat from draws (no relabelling – ordered prior)
    cat("    Getting z_hat from draws...\n")
    wst_draws <- get_z_hat_from_draws(wst_fit, A, model = "WST")
    z_hat_wst    <- wst_draws$z_hat
    z_chain_wst  <- wst_draws$z_chain
    cat("    z_hat done (draws), K =", wst_draws$K, "\n")
    
    # Compute posterior mean rho  (identity block map – labels already ordered)
    cat("    Computing posterior rho...\n")
    K_wst     <- wst_draws$K
    id_map    <- setNames(seq_len(K_wst), as.character(seq_len(K_wst)))
    rho_mean  <- compute_wst_rho_posterior(wst_fit, z_hat_wst, id_map,
                         modal_K = K_wst)
    cat("    Posterior rho done, dim =", dim(rho_mean), "\n")
    
    # Plot WST rho heatmap
    cat("    Plotting rho heatmap...\n")
    p_rho <- plot_wst_rho_heatmap(rho_mean, 
                                  title = paste("WST:", dataset_name))
    ggsave(file.path(out_path, paste0(dataset_name, "_WST_rho_heatmap.png")),
           p_rho, width = 8, height = 7, dpi = 300)
    cat("    Rho heatmap saved\n")
    
    # Plot ordered network
    cat("    Plotting ordered network...\n")
    .vw <- .canonical_viol(A, z_hat_wst)
    p_net <- plot_ordered_network(A, z_hat_wst, 
                                   title = paste("WST Network -", dataset_name),
                                   precomputed_rate = .vw$rate,
                                   precomputed_backward = .vw$backward_mass,
                                   precomputed_cross = .vw$cross_mass)
    ggsave(file.path(out_path, paste0(dataset_name, "_WST_network.png")),
           p_net, width = 8, height = 9, dpi = 300)
    cat("    Network saved\n")
    
    # Plot block-level network
    cat("    Plotting block network...\n")
        p_block_net <- plot_block_network(A, z_hat_wst,
                  title = paste("WST Block Network -", dataset_name),
                  show_cell_labels = TRUE,
                  show_legend = FALSE)
    ggsave(file.path(out_path, paste0(dataset_name, "_WST_block_network.png")),
           p_block_net, width = 8, height = 10, dpi = 300)

        p_block_net_clean <- plot_block_network(A, z_hat_wst,
                  title = paste("WST Block Network -", dataset_name),
                  show_cell_labels = FALSE,
                  show_legend = TRUE)
        ggsave(file.path(out_path, paste0(dataset_name, "_WST_block_network_clean_legend.png")),
          p_block_net_clean, width = 8, height = 10, dpi = 300)
    
    cat("    Saved WST plots to:", out_path, "\n")
  } else {
    cat("  WST fit not found, skipping...\n")
    rho_mean <- NULL
    z_hat_wst <- NULL
  }
  
  # --- SST Analysis ---
  z_chain_sst <- NULL
  if (!is.null(sst_fit_path) && file.exists(sst_fit_path)) {
    cat("  Processing SST fit...\n")
    sst_fit <- readRDS(sst_fit_path)
    
    # Get z_hat from draws (no relabelling – ordered prior)
    cat("    Getting z_hat from draws (SST)...\n")
    sst_draws   <- get_z_hat_from_draws(sst_fit, A, model = "SST")
    z_hat_sst   <- sst_draws$z_hat
    z_chain_sst <- sst_draws$z_chain
    cat("    z_hat done (draws), K =", sst_draws$K, "\n")
    
    # Compute posterior mean psi_d using only draws with the displayed K
    psi_df <- compute_sst_psi_posterior(sst_fit, K_model = sst_draws$K)
    
    # Plot SST psi bars
    p_psi <- plot_sst_psi_bars(psi_df, 
                               title = paste("SST:", dataset_name))
    ggsave(file.path(out_path, paste0(dataset_name, "_SST_psi_bars.png")),
           p_psi, width = 8, height = 5, dpi = 300)
    
    # Plot ordered network for SST
    .vs <- .canonical_viol(A, z_hat_sst)
    p_net_sst <- plot_ordered_network(A, z_hat_sst, 
                                       title = paste("SST Network -", dataset_name),
                                       precomputed_rate = .vs$rate,
                                       precomputed_backward = .vs$backward_mass,
                                       precomputed_cross = .vs$cross_mass)
    ggsave(file.path(out_path, paste0(dataset_name, "_SST_network.png")),
           p_net_sst, width = 8, height = 9, dpi = 300)
    
    # Plot block-level network for SST
        p_block_net_sst <- plot_block_network(A, z_hat_sst,
                 title = paste("SST Block Network -", dataset_name),
                 show_cell_labels = TRUE,
                 show_legend = FALSE)
    ggsave(file.path(out_path, paste0(dataset_name, "_SST_block_network.png")),
           p_block_net_sst, width = 8, height = 10, dpi = 300)

        p_block_net_sst_clean <- plot_block_network(A, z_hat_sst,
                      title = paste("SST Block Network -", dataset_name),
                      show_cell_labels = FALSE,
                      show_legend = TRUE)
        ggsave(file.path(out_path, paste0(dataset_name, "_SST_block_network_clean_legend.png")),
          p_block_net_sst_clean, width = 8, height = 10, dpi = 300)
    
    cat("    Saved SST plots to:", out_path, "\n")
  } else {
    cat("  SST fit not found, skipping...\n")
    psi_df <- NULL
    z_hat_sst <- NULL
  }
  
  # --- DCSBM Analysis ---
  dcsbm_fit_path <- find_fit_file_local(dataset_name, "DCSBM", ppc_dir)
  dcsbm_res <- NULL
  z_hat_dcsbm <- NULL
  z_chain_dcsbm <- NULL
  
  if (!is.null(dcsbm_fit_path) && file.exists(dcsbm_fit_path)) {
    cat("  Processing DCSBM fit...\n")
    dcsbm_fit <- readRDS(dcsbm_fit_path)
    
    # Get z_hat from draws (relabel each draw by block score)
    cat("    Getting z_hat from draws (DCSBM, per-draw relabel)...\n")
    dc_draws      <- get_z_hat_from_draws(dcsbm_fit, A, model = "DCSBM")
    z_hat_dcsbm   <- dc_draws$z_hat
    z_chain_dcsbm <- dc_draws$z_chain
    cat("    z_hat done (draws), K =", dc_draws$K, "\n")
    
    # Compute empirical block-level statistics from A using relabeled z_hat
    cat("    Computing empirical block statistics...\n")
    dcsbm_res <- compute_dcsbm_empirical(A, z_hat_dcsbm)
    cat("    Empirical stats done, dim =", dim(dcsbm_res$lambda_emp), "\n")
    
    # Plot DCSBM lambda (empirical rate) heatmap
    cat("    Plotting lambda heatmap...\n")
    p_lambda <- plot_dcsbm_lambda_heatmap(dcsbm_res$lambda_emp, 
                                          title = paste("DC-SBM:", dataset_name))
    ggsave(file.path(out_path, paste0(dataset_name, "_DCSBM_lambda_heatmap.png")),
           p_lambda, width = 8, height = 7, dpi = 300)
    cat("    Lambda heatmap saved\n")
    
    # Plot DCSBM rho (dominance) heatmap
    cat("    Plotting rho heatmap...\n")
    p_rho_dc <- plot_dcsbm_rho_heatmap(dcsbm_res$rho_emp,
                                        title = paste("DC-SBM Dominance:", dataset_name))
    ggsave(file.path(out_path, paste0(dataset_name, "_DCSBM_rho_heatmap.png")),
           p_rho_dc, width = 8, height = 7, dpi = 300)
    cat("    Rho heatmap saved\n")
    
    # Plot ordered network for DCSBM
    cat("    Plotting ordered network...\n")
    .vd <- .canonical_viol(A, z_hat_dcsbm)
    p_net_dc <- plot_ordered_network(A, z_hat_dcsbm, 
                                      title = paste("DC-SBM Network -", dataset_name),
                                      precomputed_rate = .vd$rate,
                                      precomputed_backward = .vd$backward_mass,
                                      precomputed_cross = .vd$cross_mass)
    ggsave(file.path(out_path, paste0(dataset_name, "_DCSBM_network.png")),
           p_net_dc, width = 8, height = 9, dpi = 300)
    cat("    Network saved\n")
    
    # Plot block-level network for DCSBM
    cat("    Plotting block network...\n")
        p_block_net_dc <- plot_block_network(A, z_hat_dcsbm,
                title = paste("DC-SBM Block Network -", dataset_name),
                show_cell_labels = TRUE,
                show_legend = FALSE)
    ggsave(file.path(out_path, paste0(dataset_name, "_DCSBM_block_network.png")),
           p_block_net_dc, width = 8, height = 10, dpi = 300)

        p_block_net_dc_clean <- plot_block_network(A, z_hat_dcsbm,
                     title = paste("DC-SBM Block Network -", dataset_name),
                     show_cell_labels = FALSE,
                     show_legend = TRUE)
        ggsave(file.path(out_path, paste0(dataset_name, "_DCSBM_block_network_clean_legend.png")),
          p_block_net_dc_clean, width = 8, height = 10, dpi = 300)
    
    cat("    Saved DCSBM plots to:", out_path, "\n")
  } else {
    cat("  DCSBM fit not found, skipping...\n")
  }
  
  # --- Combined Plot (WST | SST | DCSBM with harmonized scale) ---
  if (!is.null(rho_mean) && !is.null(psi_df) && !is.null(dcsbm_res)) {
    # Use the new combined plot with all three models
    p_combined <- plot_combined_rho_all(rho_mean, psi_df, dcsbm_res$rho_emp, dataset_name)
    ggsave(file.path(out_path, paste0(dataset_name, "_combined_rho_all.png")),
           p_combined, width = 16, height = 7, dpi = 300)
    cat("    Saved combined rho plot (WST|SST|DCSBM)\n")
  } else if (!is.null(rho_mean) && !is.null(psi_df)) {
    # Fallback to old combined plot if DCSBM not available
    p_combined <- plot_combined_heatmap_bars(rho_mean, psi_df, dataset_name)
    ggsave(file.path(out_path, paste0(dataset_name, "_combined_heatmap_bars.png")),
           p_combined, width = 14, height = 7, dpi = 300)
    cat("    Saved combined plot (WST|SST)\n")
  }
  
  # --- VI Distance Table ---
  vi_table <- compute_vi_distance_table(z_hat_wst, z_hat_sst, z_hat_dcsbm)
  if (!is.null(vi_table)) {
    vi_table$dataset <- dataset_name
    write.csv(vi_table, file.path(out_path, paste0(dataset_name, "_VI_distance_table.csv")),
              row.names = FALSE)
    cat("    Saved VI distance table\n")
    print(vi_table)
  }
  
  # --- Combined Block Networks (WST | SST | DCSBM) ---
  if (!is.null(z_hat_wst) && !is.null(z_hat_sst) && !is.null(z_hat_dcsbm)) {
    cat("    Generating combined block networks...\n")
    p_block_combined <- plot_combined_block_networks(A, z_hat_wst, z_hat_sst, z_hat_dcsbm, dataset_name)
    ggsave(file.path(out_path, paste0(dataset_name, "_combined_block_networks.png")),
           p_block_combined, width = 15, height = 10, dpi = 300)

      p_block_combined_clean <- plot_combined_block_networks_clean(A, z_hat_wst, z_hat_sst, z_hat_dcsbm, dataset_name)
      ggsave(file.path(out_path, paste0(dataset_name, "_combined_block_networks_clean.png")),
        p_block_combined_clean, width = 15, height = 10, dpi = 300)

      p_node_combined <- plot_combined_ordered_networks(
        A, z_hat_wst, z_hat_sst, z_hat_dcsbm, dataset_name,
        viol_wst   = .canonical_viol(A, z_hat_wst),
        viol_sst   = .canonical_viol(A, z_hat_sst),
        viol_dcsbm = .canonical_viol(A, z_hat_dcsbm))
      ggsave(file.path(out_path, paste0(dataset_name, "_combined_ordered_networks.png")),
        p_node_combined, width = 18, height = 9, dpi = 300)

      cat("    Saved combined block networks (annotated + clean) and combined node networks\n")
  }
  
  # =========================================================================
  # SANITY CHECKS: Node-level MFAS, violations, WST/SST conformity,
  #                block-level comparison
  # =========================================================================
  cat("  --- Running sanity checks ---\n")

  # Collect fits, z_hats, and z_chains for all available models
  fits <- list()
  z_hats <- list()
  z_chains <- list()
  regimes <- character(0)

  if (file.exists(wst_fit_path) && !is.null(z_hat_wst)) {
    fits[["WST"]] <- readRDS(wst_fit_path)
    z_hats[["WST"]] <- z_hat_wst
    z_chains[["WST"]] <- z_chain_wst
    regimes <- c(regimes, "WST")
  }
  if (file.exists(sst_fit_path) && !is.null(z_hat_sst)) {
    fits[["SST"]] <- readRDS(sst_fit_path)
    z_hats[["SST"]] <- z_hat_sst
    z_chains[["SST"]] <- z_chain_sst
    regimes <- c(regimes, "SST")
  }
  dcsbm_fit_path2 <- file.path(ppc_dir, dataset_name, "DCSBM",
                                paste0(dataset_name, "_DCSBM_fit.rds"))
  if (file.exists(dcsbm_fit_path2) && !is.null(z_hat_dcsbm)) {
    fits[["DCSBM"]] <- readRDS(dcsbm_fit_path2)
    z_hats[["DCSBM"]] <- z_hat_dcsbm
    z_chains[["DCSBM"]] <- z_chain_dcsbm
    regimes <- c(regimes, "DCSBM")
  }

  n <- nrow(A)

  # ==== (A0) AGONY PARTITION + MODEL COMPARISON ====
  cat("    Computing agony partition and model comparison...\n")
  agony_outputs <- compute_agony_outputs(
    A = A,
    z_hats = z_hats,
    dataset_name = dataset_name,
    out_path = out_path
  )

  # ==== (A) NODE-LEVEL MFAS ORDERING ====
  cat("    Computing node-level MFAS ordering...\n")
  node_mfas <- compute_node_mfas(A, max_passes = 30L)
  cat(sprintf("      n=%d, backward_mass=%.0f, total_mass=%.0f, backward_rate=%.4f, passes=%d\n",
              n, node_mfas$backward_mass, node_mfas$total_mass,
              node_mfas$backward_rate, node_mfas$n_passes))

  # ==== (B) MODEL-INDUCED NODE ORDERINGS (posterior probabilities) ====
  cat("    Computing model-induced node orderings (posterior probabilities)...\n")
  model_node_ords <- list()
  for (mod in regimes) {
    model_node_ords[[mod]] <- model_node_ordering(z_hats[[mod]], A,
                                                   alpha = 0.5,
                                                   z_chain = z_chains[[mod]])
  }

  # ==== (C) NODE-LEVEL VIOLATION RATES ====
  cat("    Computing node-level violation rates...\n")

  # (C1) Violation rate for the MFAS ordering (all edges, no blocks)
  mfas_viol <- node_violation_rate(A, node_mfas$ordering, z_hat = NULL)

  # (C2) Violation rate for each model's induced node ordering
  #      Cross-block edges only — matching .violation_stats_from_ranked
  model_node_viols <- list()
  for (mod in regimes) {
    model_node_viols[[mod]] <- node_violation_rate(A, model_node_ords[[mod]],
                                                    z_hat = z_hats[[mod]])
  }

  # ==== (D) BLOCK-LEVEL VIOLATION RATES (point-estimate, matching reference) ====
  # Already computed inline in the comparison table below using
  # .violation_stats_from_ranked, which is the reference implementation.

  # ==== (E) NODE-LEVEL WST/SST CONFORMITY ====
  cat("    Computing node-level WST/SST conformity...\n")

  # (E1) For the MFAS node ordering
  mfas_conform <- node_wst_sst_conformity(A, node_mfas$ordering,
                                           alpha = 0.5, T_triples = 10000L, seed = 42L)
  cat(sprintf("      MFAS node: thetaW=%.4f, thetaS=%.4f (prem=%d, coverage=%.3f)\n",
              mfas_conform$thetaW, mfas_conform$thetaS,
              mfas_conform$prem, mfas_conform$coverage))

  # (E2) For each model's node ordering
  model_node_conform <- list()
  for (mod in regimes) {
    mc <- node_wst_sst_conformity(A, model_node_ords[[mod]],
                                   alpha = 0.5, T_triples = 10000L, seed = 42L)
    model_node_conform[[mod]] <- mc
    cat(sprintf("      %s node: thetaW=%.4f, thetaS=%.4f (prem=%d, coverage=%.3f)\n",
                mod, mc$thetaW, mc$thetaS, mc$prem, mc$coverage))
  }

  # ==== (F) BLOCK-LEVEL WST/SST CONFORMITY (from model rho) ====
  cat("    Computing block-level WST/SST conformity...\n")
  block_conform <- list()
  for (mod in regimes) {
    emp <- compute_empirical_rho(A, z_hats[[mod]], alpha = 0.5)
    diag_order <- if (mod %in% c("WST", "SST")) "identity" else "mean"
    bc <- block_diag_rates_ext(emp$rho, T_block = 1000L, seed = 123,
                                method_order = diag_order,
                                order_direction = "strong_to_weak")
    block_conform[[mod]] <- bc
    cat(sprintf("      %s block (K=%d): thetaW=%.4f, thetaS=%.4f (prem=%d, coverage=%.3f)\n",
                mod, length(unique(z_hats[[mod]])),
                ifelse(is.na(bc$thetaW), NA, bc$thetaW),
                ifelse(is.na(bc$thetaS), NA, bc$thetaS),
                bc$prem, bc$coverage))
  }

  # ==== (G) BUMP CHART: Node rankings across models ====
  cat("    Generating bump chart...\n")
  p_bump <- plot_bump_chart(A, z_hats, node_mfas, z_chains = z_chains,
                            alpha = 0.5, dataset_name = dataset_name)
  ggsave(file.path(out_path, paste0(dataset_name, "_bump_chart.png")),
         p_bump, width = 10, height = max(8, n * 0.18), dpi = 300)
  cat("    Saved bump chart\n")

  # ==== (H) KENDALL'S TAU BETWEEN NODE ORDERINGS ====
  cat("    Computing Kendall's tau between node orderings...\n")
  tau_rows <- list()
  all_node_ords <- c(list(MFAS = node_mfas$ordering), model_node_ords)
  ord_names <- names(all_node_ords)

  if (length(ord_names) >= 2) {
    for (i in seq_len(length(ord_names) - 1)) {
      for (j in (i + 1):length(ord_names)) {
        m1 <- ord_names[i]
        m2 <- ord_names[j]
        tau_val <- kendall_tau_distance(all_node_ords[[m1]], all_node_ords[[m2]])
        tau_rows[[length(tau_rows) + 1]] <- data.frame(
          dataset = dataset_name,
          ordering1 = m1,
          ordering2 = m2,
          kendall_tau = round(tau_val, 4),
          n_nodes = n,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  kendall_table <- if (length(tau_rows) > 0) do.call(rbind, tau_rows) else NULL

  # ==== ASSEMBLE COMPARISON TABLE ====
  # One row per "source": MFAS (node), each model (node), each model (block)
  cat("    Assembling comparison table...\n")
  comp_rows <- list()

  # Row for MFAS node ordering
  comp_rows[["MFAS_node"]] <- data.frame(
    dataset          = dataset_name,
    source           = "MFAS",
    level            = "node",
    K                = n,
    violation_rate   = round(mfas_viol$violation_rate, 6),
    backward_mass    = mfas_viol$backward_mass,
    total_mass       = mfas_viol$total_mass,
    thetaW           = round(mfas_conform$thetaW, 4),
    thetaS           = round(mfas_conform$thetaS, 4),
    thetaW_all       = round(mfas_conform$thetaW_all, 4),
    thetaS_all       = round(mfas_conform$thetaS_all, 4),
    premise_coverage = round(mfas_conform$coverage, 4),
    n_triples        = mfas_conform$total,
    stringsAsFactors = FALSE
  )

  for (mod in regimes) {
    K_mod <- length(unique(z_hats[[mod]]))

    # Model node-level row
    mnv <- model_node_viols[[mod]]
    mnc <- model_node_conform[[mod]]
    comp_rows[[paste0(mod, "_node")]] <- data.frame(
      dataset          = dataset_name,
      source           = mod,
      level            = "node",
      K                = K_mod,
      violation_rate   = round(mnv$violation_rate, 6),
      backward_mass    = mnv$backward_mass,
      total_mass       = mnv$total_mass,
      thetaW           = round(mnc$thetaW, 4),
      thetaS           = round(mnc$thetaS, 4),
      thetaW_all       = round(mnc$thetaW_all, 4),
      thetaS_all       = round(mnc$thetaS_all, 4),
      premise_coverage = round(mnc$coverage, 4),
      n_triples        = mnc$total,
      stringsAsFactors = FALSE
    )

    # Model block-level row — use .violation_stats_from_ranked for consistency
    #   with transitivity_check_helper.R reference
    z_ranked <- .rank_partition_by_strength(
      A = A, z_vec = z_hats[[mod]], alpha = 0.5,
      order_direction = "strong_to_weak"
    )
    vstats <- .violation_stats_from_ranked(A = A, z_ranked = z_ranked)
    bc  <- block_conform[[mod]]
    comp_rows[[paste0(mod, "_block")]] <- data.frame(
      dataset          = dataset_name,
      source           = mod,
      level            = "block",
      K                = K_mod,
      violation_rate   = round(vstats$rate, 6),
      backward_mass    = vstats$count,
      total_mass       = vstats$cross_mass,
      thetaW           = round(ifelse(is.na(bc$thetaW), NA, bc$thetaW), 4),
      thetaS           = round(ifelse(is.na(bc$thetaS), NA, bc$thetaS), 4),
      thetaW_all       = round(ifelse(is.na(bc$thetaW_all), NA, bc$thetaW_all), 4),
      thetaS_all       = round(ifelse(is.na(bc$thetaS_all), NA, bc$thetaS_all), 4),
      premise_coverage = round(bc$coverage, 4),
      n_triples        = if (bc$coverage > 0) as.integer(round(bc$prem / bc$coverage)) else bc$prem,
      stringsAsFactors = FALSE
    )
  }

  comparison_table <- do.call(rbind, comp_rows)
  rownames(comparison_table) <- NULL

  # Save tables
  write.csv(comparison_table,
            file.path(out_path, paste0(dataset_name, "_node_vs_block_comparison.csv")),
            row.names = FALSE)
  cat("    Saved node-vs-block comparison table\n")
  cat("\n    === Node-level vs Block-level Comparison ===\n")
  print(comparison_table)

  if (!is.null(kendall_table)) {
    write.csv(kendall_table,
              file.path(out_path, paste0(dataset_name, "_node_kendall_tau_table.csv")),
              row.names = FALSE)
    cat("    Saved node-level Kendall's tau table\n")
    print(kendall_table)
  }

  # Save item-level orderings (MFAS + model-induced)
  cat("    Saving item-level orderings...\n")
  item_rows <- list()
  node_names <- if (!is.null(rownames(A))) rownames(A) else as.character(seq_len(n))

  # MFAS ordering
  mfas_df <- data.frame(
    rank        = seq_len(n),
    node        = node_mfas$ordering,
    node_name   = node_names[node_mfas$ordering],
    model       = "MFAS",
    dataset     = dataset_name,
    stringsAsFactors = FALSE
  )
  item_rows[["MFAS"]] <- mfas_df

  # Model node orderings
  for (mod in regimes) {
    ord_mod <- model_node_ords[[mod]]
    mod_df <- data.frame(
      rank      = seq_len(n),
      node      = ord_mod,
      node_name = node_names[ord_mod],
      model     = mod,
      dataset   = dataset_name,
      stringsAsFactors = FALSE
    )
    item_rows[[mod]] <- mod_df
  }
  item_ordering_table <- do.call(rbind, item_rows)
  rownames(item_ordering_table) <- NULL
  write.csv(item_ordering_table,
            file.path(out_path, paste0(dataset_name, "_item_ordering_table.csv")),
            row.names = FALSE)
  cat("    Saved item-level ordering table\n")

  cat("  --- Sanity checks complete ---\n")

  # NOTE: legacy file.copy duplicating osbm_visualizations contents into the
  # parent <ds>/ folder has been removed. Plots now live only in
  # output/application/plots/<dataset>/osbm_visualizations/.

  cat("Done with", dataset_name, "\n\n")
  
  invisible(list(
    A = A,
    rho_mean_wst = rho_mean,
    psi_df = psi_df,
    dcsbm_res = dcsbm_res,
    z_hat_wst = z_hat_wst,
    z_hat_sst = z_hat_sst,
    z_hat_dcsbm = z_hat_dcsbm,
    z_chains = z_chains,
    agony_fit = agony_outputs$agony_fit,
    agony_partition_table = agony_outputs$agony_partition_table,
    agony_table = agony_outputs$agony_table,
    agony_kendall_table = agony_outputs$agony_kendall_table,
    agony_winner = agony_outputs$agony_winner,
    comparison_table = comparison_table,
    kendall_table = kendall_table,
    node_mfas = node_mfas
  ))
}

export_application_manuscript_figures <- function(results_list,
                                                  figure_dir = file.path("tex file", "Figures"),
                                                  manuscript_dir = file.path("output", "application", "plots", "_cross_dataset", "manuscript"),
                                                  run_id = NULL) {
  # Resolve run_id from env var if not provided
  if (is.null(run_id)) {
    rd <- Sys.getenv("APP_RUN_DIR", unset = "")
    if (nzchar(rd)) run_id <- basename(rd)
  }
  paper_fig_dir <- if (!is.null(run_id) && nzchar(run_id))
    file.path("output", "paper", "figures", run_id) else NULL

  fs::dir_create(figure_dir)
  fs::dir_create(manuscript_dir)
  if (!is.null(paper_fig_dir)) fs::dir_create(paper_fig_dir)

  for (ds in names(results_list)) {
    res <- results_list[[ds]]
    if (is.null(res) || is.null(res$A)) next
    if (is.null(res$z_hat_wst) || is.null(res$z_hat_sst) || is.null(res$z_hat_dcsbm)) next

    p_clean <- plot_combined_block_networks_clean(
      A = res$A,
      z_wst = res$z_hat_wst,
      z_sst = res$z_hat_sst,
      z_dcsbm = res$z_hat_dcsbm,
      dataset_name = ds
    )
    save_plot_bundle_dirs(
      plot_obj = p_clean,
      dirs = c(figure_dir, manuscript_dir, paper_fig_dir),
      stem = paste0(ds, "_combined_block_networks_clean"),
      width = 7.2,
      height = 2.95
    )
  }

  p_human <- plot_best_model_ordered_matrices(results_list)
  if (!is.null(p_human)) {
    save_plot_bundle_dirs(
      plot_obj = p_human,
      dirs = c(figure_dir, manuscript_dir, paper_fig_dir),
      stem = "human_best_model_ordered_matrices",
      width = 7.0,
      height = 3.55
    )
  }

  p_animal <- plot_animal_hierarchy_profile(results_list)
  if (!is.null(p_animal)) {
    save_plot_bundle_dirs(
      plot_obj = p_animal,
      dirs = c(figure_dir, manuscript_dir, paper_fig_dir),
      stem = "animal_hierarchy_profile",
      width = 6.7,
      height = 3.0
    )
  }

  # Provenance manifest in the versioned paper-figures dir
  if (!is.null(paper_fig_dir)) {
    manifest_path <- file.path(paper_fig_dir, "MANIFEST.txt")
    writeLines(c(
      paste0("run_id: ", run_id),
      paste0("generated_by: scripts/analysis/osbm_visualization.R"),
      paste0("generated_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
    ), manifest_path)
  }
}

# =============================================================================
# RUN FOR ALL DATASETS
# =============================================================================
run_all_osbm_visualizations <- function() {
  # Get available datasets from ppc dir
  ppc_dir <- "output/application/ppc"
  datasets <- list.dirs(ppc_dir, recursive = FALSE, full.names = FALSE)
  datasets <- datasets[!datasets %in% c(".DS_Store")]
  
  # Also pick up datasets present in the hardcoded run dir but not yet in ppc
  run_dir <- Sys.getenv("APP_RUN_DIR",
    unset = "./output/application/raw/application_run_20260411_163055")
  if (dir.exists(run_dir)) {
    run_wst_files <- list.files(run_dir, pattern = "_WST_fit\\.rds$", full.names = FALSE)
    run_datasets  <- sub("_WST_fit\\.rds$", "", run_wst_files)
    datasets <- unique(c(datasets, run_datasets))
  }
  
  cat("Found datasets:", paste(datasets, collapse = ", "), "\n\n")
  
  # Accumulate cross-dataset tables
  all_results <- list()
  all_agony_tables <- list()
  all_agony_kendall_tables <- list()
  all_agony_winners <- list()
  all_comparison_tables <- list()
  all_kendall_tables    <- list()

  for (ds in datasets) {
    tryCatch({
      res <- generate_osbm_visualizations(ds, ppc_dir = ppc_dir)
      all_results[[ds]] <- res
      all_agony_tables[[ds]] <- res$agony_table
      all_agony_kendall_tables[[ds]] <- res$agony_kendall_table
      all_agony_winners[[ds]] <- res$agony_winner
      all_comparison_tables[[ds]] <- res$comparison_table
      all_kendall_tables[[ds]]    <- res$kendall_table
    }, error = function(e) {
      cat("Error processing", ds, ":", conditionMessage(e), "\n")
    })
  }

  if (length(all_results) > 0) {
    export_application_manuscript_figures(all_results)
    cat("Saved manuscript-oriented application figures\n")
  }

  # Save combined cross-dataset summary tables
  out_root <- "output/application/plots/_cross_dataset"
  fs::dir_create(out_root)
  if (length(all_agony_tables) > 0) {
    combined_agony <- do.call(rbind, Filter(Negate(is.null), all_agony_tables))
    write.csv(combined_agony,
              file.path(out_root, "all_datasets_agony_comparison.csv"),
              row.names = FALSE)
    cat("\nSaved combined agony comparison:",
        file.path(out_root, "all_datasets_agony_comparison.csv"), "\n")
    print(combined_agony)
  }
  if (length(all_agony_kendall_tables) > 0) {
    combined_agony_kendall <- do.call(rbind, Filter(Negate(is.null), all_agony_kendall_tables))
    if (nrow(combined_agony_kendall) > 0) {
      write.csv(combined_agony_kendall,
                file.path(out_root, "all_datasets_agony_kendall_tau.csv"),
                row.names = FALSE)
      cat("Saved combined agony Kendall tau table:",
          file.path(out_root, "all_datasets_agony_kendall_tau.csv"), "\n")
      print(combined_agony_kendall)
    }
  }
  if (length(all_agony_winners) > 0) {
    combined_agony_winners <- do.call(rbind, Filter(Negate(is.null), all_agony_winners))
    if (nrow(combined_agony_winners) > 0) {
      write.csv(combined_agony_winners,
                file.path(out_root, "all_datasets_agony_winners.csv"),
                row.names = FALSE)
      cat("Saved combined agony winners:",
          file.path(out_root, "all_datasets_agony_winners.csv"), "\n")
      print(combined_agony_winners)
    }
  }
  if (length(all_comparison_tables) > 0) {
    combined_comp <- do.call(rbind, Filter(Negate(is.null), all_comparison_tables))
    write.csv(combined_comp,
              file.path(out_root, "all_datasets_node_vs_block_comparison.csv"),
              row.names = FALSE)
    cat("\nSaved combined node-vs-block comparison:",
        file.path(out_root, "all_datasets_node_vs_block_comparison.csv"), "\n")
    cat("\n=== CROSS-DATASET NODE-vs-BLOCK COMPARISON ===\n")
    print(combined_comp)
  }
  if (length(all_kendall_tables) > 0) {
    combined_kt <- do.call(rbind, Filter(Negate(is.null), all_kendall_tables))
    if (nrow(combined_kt) > 0) {
      write.csv(combined_kt,
                file.path(out_root, "all_datasets_node_kendall_tau.csv"),
                row.names = FALSE)
      cat("Saved combined Kendall tau table:",
          file.path(out_root, "all_datasets_node_kendall_tau.csv"), "\n")
    }
  }
  
  cat("\nAll visualizations and sanity checks complete!\n")
}

# Auto-run when script is executed directly
if (sys.nframe() == 0) {
  run_all_osbm_visualizations()
}
