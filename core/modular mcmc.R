############################################################################
# Modular OSBM Gibbs – clean main (WST -> pair ψ, SST -> distance ψ)
#
# Returns a list with z, psi, kappa, eta. Any parameter not in `free`
# is fixed to the value in `truth`.
#
# External functions expected in the environment:
#   Pair (WST):
#     - z_update_osbm_pair(i, A, z, eta, kappa, psi, alpha_vec, log1p_tab, log_kappa,
#                          i_idx, j_idx, N_edge, edge_by_node)
#     - aggregate_by_pair(), pair_totals(), draw_omega_pair(), update_psi_pair(),
#       make_log1p_tables_pair()
#   Distance (SST):
#     - z_update_osbm_distance(i, A, z, eta, kappa, psi, alpha_vec, log1p_tab, log_kappa,
#                              i_idx, j_idx, N_edge, edge_by_node)
#     - aggregate_by_distance(), distance_totals(), draw_omega_bar(),
#       make_log1p_tables(), update_psi_sst()
############################################################################

modular_osbm_sampler <- function(
    A, K, truth = NA,
    free   = c("psi","kappa","eta","z"),
    n_iter = 4000, burn = 500, thin = 1,
    verbose = FALSE,
    psi_constraint = c("WST","SST"),
    seed = NULL,
    hyper = NULL
){
  ## ---- hygiene -------------------------------------------------------------
  op <- options(
    warnPartialMatchArgs   = TRUE,
    warnPartialMatchDollar = TRUE,
    warnPartialMatchAttr   = TRUE
  )
  on.exit(options(op), add = TRUE)
  
  if (!is.null(seed)) set.seed(as.integer(seed))
  psi_constraint <- match.arg(psi_constraint)
  
  if (!all(free %in% c("psi","kappa","eta","z")))
    stop("`free` must be subset of {'psi','kappa','eta','z'}.", call. = FALSE)

  needs <- setdiff(c("psi","kappa","eta","z"), free)
  if (length(needs) > 0) {
    if (exists(".assert_truth", mode = "function", inherits = TRUE)) {
      .assert_truth(truth, needs)
    } else {
      if (is.null(truth) || is.na(truth)) {
        stop("`truth` must be a named list when parameters are fixed.", call. = FALSE)
      }
      miss <- setdiff(needs, names(truth))
      if (length(miss)) {
        stop("`truth` is missing: ", paste(miss, collapse = ", "), call. = FALSE)
      }
    }
  }
  
  .assert_scalar_int(K, "K")
  .assert_scalar_int(n_iter, "n_iter")
  .assert_scalar_int(burn, "burn")
  .assert_scalar_int(thin, "thin")
  if (burn >= n_iter) stop("`burn` must be < `n_iter`.", call. = FALSE)
  if (thin < 1L)      stop("`thin` must be >= 1.", call. = FALSE)
  
  n <- nrow(A)
  .assert_scalar_int(n, "n (nrow(A))")
  .assert_matrix(A, n, "A")
  if (!inherits(A, "Matrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  
  if (!requireNamespace("BayesLogit", quietly = TRUE))
    stop("Package 'BayesLogit' is required.", call. = FALSE)
  if (psi_constraint == "SST" &&
      !exists("update_psi_sst", mode = "function", inherits = TRUE))
    stop("`psi_constraint='SST'` requested but `update_psi_sst()` is not visible.", call. = FALSE)
  
  ## ---- hyper-parameters ----------------------------------------------------
  default_hyper <- list(
    a_kappa = 1, b_kappa = 1,
    a_eta   = 1, b_eta   = 1,
    mu0     = 1, sigma0  = 10,   # prior for ψ
    tau0    = 5                 # extra prior scale (SST)
  )
  if (!is.null(hyper)) {
    if (!is.list(hyper)) stop("`hyper` must be a list when provided.", call. = FALSE)
    unknown <- setdiff(names(hyper), names(default_hyper))
    if (length(unknown)) {
      stop("Unknown entries in `hyper`: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    default_hyper[names(hyper)] <- hyper
  }
  hyper <- default_hyper
  
  scalar_numeric <- function(x, name, positive = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
      stop("`hyper$", name, "` must be a finite numeric scalar.", call. = FALSE)
    }
    if (positive && x <= 0) {
      stop("`hyper$", name, "` must be > 0.", call. = FALSE)
    }
  }
  scalar_numeric(hyper$a_kappa, "a_kappa", positive = TRUE)
  scalar_numeric(hyper$b_kappa, "b_kappa", positive = TRUE)
  scalar_numeric(hyper$a_eta,   "a_eta",   positive = TRUE)
  scalar_numeric(hyper$b_eta,   "b_eta",   positive = TRUE)
  scalar_numeric(hyper$mu0,     "mu0")
  scalar_numeric(hyper$sigma0,  "sigma0",  positive = TRUE)
  scalar_numeric(hyper$tau0,    "tau0",    positive = TRUE)
  
  ## ---- flags ---------------------------------------------------------------
  psi_est   <- "psi"   %in% free
  z_est     <- "z"     %in% free
  eta_est   <- "eta"   %in% free
  kappa_est <- "kappa" %in% free
  
  ## ---- precompute edge list ------------------------------------------------
  A_sym <- A + t(A)
  idx   <- which(A_sym > 0 & upper.tri(A_sym), arr.ind = TRUE)
  i_idx <- idx[,1];  j_idx <- idx[,2]
  N_edge <- A_sym[idx]
  
  edge_by_node <- vector("list", n)
  for (e in seq_len(nrow(idx))) {
    edge_by_node[[i_idx[e]]] <- c(edge_by_node[[i_idx[e]]], e)
    edge_by_node[[j_idx[e]]] <- c(edge_by_node[[j_idx[e]]], e)
  }
  
  ## ---- initial values ------------------------------------------------------
  # z
  z <- if (z_est) {
    sample.int(K, n, replace = TRUE)
  } else {
    .assert_len(truth$z, n, "truth$z")
    as.integer(truth$z)
  }
  alpha_vec <- rep(5, K)
  
  # ψ representation decided by psi_constraint:
  #   WST -> pair      (KxK, upper-tri positive)
  #   SST -> distance  (length K-1)
  psi_mode <- if (psi_constraint == "WST") "pair" else "distance"
  
  # psi
  if (psi_est) {
    if (psi_mode == "pair") {
      psi <- matrix(0, K, K)
      psi[upper.tri(psi)] <- seq_len(K*(K-1)/2) / (K*(K-1)/2)
    } else { # distance
      psi <- seq(1, K, length.out = K-1)
    }
  } else {
    if (psi_mode == "pair") {
      if (!is.matrix(truth$psi) || any(dim(truth$psi) != c(K,K)))
        stop("With `psi_constraint='WST'`, `truth$psi` must be a KxK matrix (upper-tri used).",
             call. = FALSE)
      psi <- truth$psi
    } else {
      .assert_len(truth$psi, K-1, "truth$psi (SST distance mode)")
      psi <- as.numeric(truth$psi)
    }
  }
  
  # eta
  if (eta_est) {
    eta <- rgamma(n, hyper$a_eta, hyper$b_eta)
    n_k <- as.integer(tabulate(z, nbins = K))
    for (k in seq_len(K)) {
      idx_k <- which(z == k)
      if (!length(idx_k)) next
      s_k <- sum(eta[idx_k])
      if (s_k > 0) eta[idx_k] <- n_k[k] * eta[idx_k] / s_k
      else         eta[idx_k] <- n_k[k] / length(idx_k)
    }
  } else {
    .assert_len(truth$eta, n, "truth$eta")
    eta <- truth$eta
  }
  
  # kappa
  if (kappa_est) {
    p <- pmin(z[i_idx], z[j_idx]); q <- pmax(z[i_idx], z[j_idx])
    p <- factor(p, levels = 1:K);  q <- factor(q, levels = 1:K)
    Rkl <- as.matrix(xtabs(N_edge ~ p + q, drop.unused.levels = FALSE))
    
    E_k   <- tapply(eta,   factor(z, levels = 1:K), sum);   E_k[is.na(E_k)] <- 0
    eta2k <- tapply(eta^2, factor(z, levels = 1:K), sum); eta2k[is.na(eta2k)] <- 0
    
    Tkl <- outer(E_k, E_k, `*`)
    diag(Tkl) <- pmax((E_k^2 - eta2k) / 2, 0)
    
    kappa <- (Rkl + hyper$a_kappa) / (Tkl + hyper$b_kappa)
    kappa[lower.tri(kappa)] <- t(kappa)[lower.tri(kappa)]
  } else {
    .assert_matrix(truth$kappa, K, "truth$kappa")
    kappa <- truth$kappa
  }
  log_kappa <- log(pmax(kappa, 1e-15))
  
  # log(1+exp) tables for ψ (needed by z-updater)
  log1p_tab <- if (psi_mode == "pair") make_log1p_tables_pair(psi) else make_log1p_tables(psi)
  
  ## ---- storage -------------------------------------------------------------
  keep_seq    <- seq(burn + 1L, n_iter, by = thin)
  n_keep      <- length(keep_seq)
  
  draws_z     <- matrix(0L, n_keep, n)
  draws_kappa <- array(0, dim = c(n_keep, K, K))
  draws_eta   <- matrix(0,  n_keep, n)
  if (psi_mode == "pair") {
    draws_psi <- array(0, dim = c(n_keep, K, K))  # upper-tri meaningful
  } else {
    draws_psi <- matrix(0, n_keep, K-1)
  }
  
  keep <- 0L
  t0 <- proc.time()[3]
  
  ## ---- dispatchers ---------------------------------------------------------
  z_update_fun <- switch(
    psi_mode,
    pair     = get("z_update_osbm_pair",     mode = "function"),
    distance = get("z_update_osbm_distance", mode = "function")
  )
  
  ## ---- main Gibbs loop -----------------------------------------------------
  for (it in seq_len(n_iter)) {
    
    ## ---- z block -----------------------------------------------------------
    if (z_est) {
      for (i in seq_len(n)) {
        z[i] <- z_update_fun(
          i,
          A = A, z = z, eta = eta, kappa = kappa,
          psi = psi,
          alpha_vec = alpha_vec,
          log1p_tab = log1p_tab,
          log_kappa = log_kappa,
          i_idx = i_idx, j_idx = j_idx, N_edge = N_edge,
          edge_by_node = edge_by_node
        )
      }
    }
    
    ## ---- psi block ---------------------------------------------------------
    if (psi_est) {
      if (psi_mode == "pair") {                  # WST
        bar_y_mat <- aggregate_by_pair(K,
                                       z_i = z[i_idx],
                                       z_j = z[j_idx],
                                       A_ij = A[cbind(i_idx, j_idx)],
                                       N_edge = N_edge)
        B_mat     <- pair_totals(K,
                                 z_i = z[i_idx],
                                 z_j = z[j_idx],
                                 N_edge = N_edge)
        bar_omega <- draw_omega_pair(B_mat, psi_mat = psi)
        
        psi <- update_psi_pair(bar_y = bar_y_mat, bar_omega = bar_omega,
                               mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2,
                               trunc = TRUE)  # half-Normal => ψ_{kl} >= 0
        
        log1p_tab <- make_log1p_tables_pair(psi)
        
      } else {                                   # SST (distance)
        agg      <- aggregate_by_distance(K,
                                          z_i = z[i_idx],
                                          z_j = z[j_idx],
                                          A_ij = A[cbind(i_idx, j_idx)],
                                          N_edge = N_edge)
        bar_y    <- agg$bar_y
        if (exists("assert_sst_invariants", mode = "function", inherits = TRUE)) {
          assert_sst_invariants(z, kappa, psi, where = sprintf("it=%d (psi update)", it))
        } else if (length(psi) != max(K - 1L, 0L)) {
          stop(sprintf("SST invariant failed: length(psi)=%d != K-1=%d",
                       length(psi), max(K - 1L, 0L)), call. = FALSE)
        }
        B_d      <- distance_totals(K, z_i = z[i_idx], z_j = z[j_idx], N_edge = N_edge)
        bar_omega <- draw_omega_bar(B_d = B_d, psi = psi)
        
        psi <- update_psi_sst(K = K, bar_y = bar_y, bar_omega = bar_omega,
                              psi_curr = psi,
                              mu0 = hyper$mu0, sig2_0 = hyper$sigma0^2,
                              tau2_0 = hyper$tau0^2,
                              n_inner_sweeps = 1)
        
        log1p_tab <- make_log1p_tables(psi)
      }
    }
    
    ## ---- eta block ---------------------------------------------------------
    if (eta_est) {
      G_i <- Matrix::rowSums(A) + Matrix::colSums(A)
      E_k <- tapply(eta, factor(z, levels = 1:K), sum); E_k[is.na(E_k)] <- 0
      
      for (i in seq_len(n)) {
        k      <- z[i]
        E_k[k] <- E_k[k] - eta[i]
        rate_i <- hyper$b_eta + sum(kappa[k, ] * E_k)
        eta[i] <- rgamma(1,
                         shape = hyper$a_eta + G_i[i],
                         rate  = max(rate_i, 1e-10))
        E_k[k] <- E_k[k] + eta[i]
      }
      
      # blockwise normalization Σ_{i:z_i=k} η_i = n_k
      n_k <- as.integer(tabulate(z, nbins = K))
      for (k in seq_len(K)) {
        idx_k <- which(z == k); if (!length(idx_k)) next
        s_k <- sum(eta[idx_k])
        if (s_k > 0) {
          sc <- n_k[k] / s_k
          eta[idx_k] <- eta[idx_k] * sc
        } else {
          eta[idx_k] <- n_k[k] / length(idx_k)
        }
      }
    }
    
    ## ---- kappa block -------------------------------------------------------
    if (kappa_est) {
      p <- pmin(z[i_idx], z[j_idx]); q <- pmax(z[i_idx], z[j_idx])
      p <- factor(p, levels = 1:K);  q <- factor(q, levels = 1:K)
      
      Rkl <- as.matrix(xtabs(N_edge ~ p + q, drop.unused.levels = FALSE))
      
      E_k   <- tapply(eta,   factor(z, levels = 1:K), sum);   E_k[is.na(E_k)] <- 0
      eta2k <- tapply(eta^2, factor(z, levels = 1:K), sum); eta2k[is.na(eta2k)] <- 0
      
      Tkl <- outer(E_k, E_k, `*`)
      diag(Tkl) <- pmax((E_k^2 - eta2k) / 2, 0)
      
      kappa <- matrix(0, K, K)
      for (k in seq_len(K)) for (l in k:K) {
        val <- rgamma(1,
                      shape = hyper$a_kappa + Rkl[k, l],
                      rate  = hyper$b_kappa + Tkl[k, l])
        kappa[k, l] <- kappa[l, k] <- val
      }
      log_kappa <- log(pmax(kappa, 1e-15))
    }
    
    ## ---- save draws --------------------------------------------------------
    if (it %in% keep_seq) {
      keep <- keep + 1L
      draws_z[keep, ]      <- z
      draws_kappa[keep,, ] <- kappa
      draws_eta[keep, ]    <- eta
      if (psi_mode == "pair") {
        draws_psi[keep,, ] <- psi
      } else {
        draws_psi[keep, ]  <- psi
      }
    }
    
    if (verbose && it %% 500 == 0)
      cat("iter", it, "/", n_iter, "saved:", keep, " ",
          round(proc.time()[3]-t0,2), "s \r")
  }
  
  invisible(list(
    z     = draws_z,
    psi   = draws_psi,
    kappa = draws_kappa,
    eta   = draws_eta,
    keep  = keep_seq,
    meta  = list(psi_mode = psi_mode, psi_constraint = psi_constraint)
  ))
}



set.seed(432)
library(Matrix)

DEMO <- F
if(DEMO){
  ## --- ground truth (easy) ---
  n    <- 100;  K <- 3
  
  z_t  <- rep(1:K, length.out = n)
  psi_t= seq(1, 5, length.out = K-1)
  kappa_t <- matrix(K*K,nrow = K,ncol = K)
  kappa_t[lower.tri(kappa_t)] = t(kappa_t)[lower.tri(kappa_t)]
  
  eta_t <- runif(n,0.8,1.2)
  eta_t <- n*eta_t/sum(eta_t)
  ## --- simulate network ---
  A   <- Matrix(0, n, n, sparse = TRUE)
  for(i in 1:(n-1)) for(j in (i+1):n){
    lam <- eta_t[i]*eta_t[j]*kappa_t[z_t[i], z_t[j]]
    N   <- rpois(1, lam)
    if(N){
      d <- abs(z_t[i]-z_t[j])
      pr<- if(d==0) 0.5 else
        logistic(if(z_t[i] < z_t[j]) psi_t[d] else -psi_t[d])
      fwd <- rbinom(1, N, pr)
      A[i, j] <- fwd;  A[j, i] <- N - fwd
    }
  }
  
  ## --- run MCMC that infers the chosen parameters ---
  alpha_vec <- rep(5, K)     # same α you use inside the sampler
  
  ## --- 1. Sanity check -------------------------------------------------------
  run_osbm_checks(A         = A,
                  z         = z_t,
                  eta       = eta_t,
                  kappa     = kappa_t,
                  psi       = psi_t,
                  alpha_vec = alpha_vec,
                  verbose   = TRUE)
  #choice of the parameters to be inferred
  
  free   = c('z','kappa','psi','eta')
  
  n_iter = 500
  burn = 100
  thin = 1
  truth=list(kappa = kappa_t,
             psi   = psi_t,
             eta   = eta_t,
             z     = z_t)
  #---------------------
  #Running the sampler
  #---------------------
  out <- modular_osbm_sampler(
    A      = A,
    K      = K,
    truth  = truth,
    free   = free,
    n_iter = n_iter, 
    burn = burn, thin = thin,verbose = T)
  
  
  #-------------------
  # Inference: point estimate and uncertainty quantification
  #-------------------
  
  #adjusting for label switching
  out_relab = relabel_osbm(out,z_ref = z_t)
  
  #point estimate for z and Adj.Rand Index
  psm  <- comp.psm(out$z)             
  hat  <- minVI(psm)$cl
  cat("Adjusted Rand =", fossil::adj.rand.index(hat, z_t), "\n")
  
  
  colMeans(out_relab$psi) 
  apply(out_relab$kappa, c(2:3), mean)[, ]   
  mean(abs(colMeans(out$eta) - eta_t))
  
  #----- Plotting similarity and adjacency matrices-------------
  # Create row and column indices
  indices <- expand.grid(row = 1:n, col = 1:n)
  
  z_df <- data.frame(items = 1:n, 
                     z = as.vector(z_t),
                     unique_identifier = runif(n, -0.001,0.001))%>%
    mutate(unique_identifier = unique_identifier+z)
  
  
  # Convert the matrix to a data frame
  z_df_complete <- data.frame(
    row = indices$row,
    col = indices$col,
    similarity_value = NA,
    Y = NA
  )
  
  for (i in seq_len(nrow(z_df_complete))) {
    z_df_complete$Y[i] <- A[z_df_complete$col[i], z_df_complete$row[i]]
  }
  for (i in seq_len(nrow(z_df_complete))) {
    z_df_complete$similarity_value[i] <- psm[z_df_complete$col[i], z_df_complete$row[i]]
  }
  for(i in seq_len(nrow(z_df_complete))){
    z_df_complete$marginal_victories_row[i] <- sum(A[z_df_complete$row[i],]/z_df$unique_identifier)/sum(A[,z_df_complete$row[i]]*z_df$unique_identifier)
  }
  for(i in seq_len(nrow(z_df_complete))){
    z_df_complete$marginal_victories_col[i] <- sum(A[z_df_complete$col[i],]/z_df$unique_identifier)/sum(A[,z_df_complete$col[i]]*z_df$unique_identifier)
  }
  
  plot_df = z_df_complete%>%
    inner_join(z_df, by = c("row" = "items")) %>%
    dplyr::rename(row_z = z) %>%
    inner_join(z_df, by = c("col" = "items")) %>%
    dplyr::rename(col_z = z) %>%
    mutate(row = factor(row, levels = unique(row[order(row_z, -marginal_victories_row)])),
           col = factor(col, levels = unique(col[order(col_z, -marginal_victories_col, decreasing = TRUE)]))) 
  
  v_lines_list = list()
  for (k in 1:(K-1)) {
    lines_df = plot_df %>% filter(row_z == k, col_z == k)
    v_lines_list[[k]] <- lines_df$row[which.min(lines_df$marginal_victories_row)]
  }
  
  
  adjacency_m<- ggplot(plot_df, aes(x = row, y = col)) +
    geom_tile(aes(fill = Y), color = "gray", show.legend = FALSE) +
    scale_fill_gradient(low = "white", high = "black") +
    geom_ysidetile(aes(color = factor(col_z)), show.legend = FALSE, width = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  adjacency_m #actual data
  
  similarity_m <- ggplot(plot_df, aes(x = row, y = col)) +
    geom_tile(aes(fill = similarity_value), color = "gray", show.legend = FALSE) +
    scale_fill_gradient(low = "white", high = "black") +
    geom_vline(xintercept = as.numeric(unlist(v_lines_list)) + 0.5, color = 'red3') +
    geom_hline(yintercept = as.numeric(unlist(v_lines_list)) - 0.5, color = 'red3') +
    geom_ysidetile(aes(color = factor(col_z)), show.legend = FALSE, width = 0.5) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  similarity_m #posterior similarity matrix
  
  
}
