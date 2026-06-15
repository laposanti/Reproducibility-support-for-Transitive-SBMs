# Deterministic per-cell seed helper
seed_from <- function(...) {
  key <- paste(..., collapse = "|")
  if ("digest2int" %in% getNamespaceExports("digest")) {
    s <- digest::digest2int(key)
  } else {
    h  <- digest::digest(key, algo = "xxhash32", serialize = FALSE)
    hi <- strtoi(substr(h, 1, 4), base = 16L)
    lo <- strtoi(substr(h, 5, 8), base = 16L)
    s  <- as.integer((hi * 65536 + lo) %% .Machine$integer.max)
  }
  if (is.na(s) || s <= 0L) s <- 1L
  if (s >= .Machine$integer.max) s <- .Machine$integer.max - 1L
  s
}


# --- Plugin: parallel streaming save with file locking -----------------------
fit_osbm_model <- function(A, K,
                           ordering = c("NONE","WST","SST"),
                           n_iter = 5000, burn = 2000, thin = 2,
                           verbose = FALSE) {
  ordering <- match.arg(ordering)
  
  out <- modular_osbm_sampler(
    A = A, K = K,
    free = c("z","kappa","psi","eta"),
    n_iter = n_iter, burn = burn, thin = thin,
    verbose = verbose, psi_constraint = ordering
  )
  
  # Relabel with orientation logic (assumed available from your helpers)
  out_relab <- relabel_osbm(out, ordering = ordering)
  
  # Posterior means (plug-in)
  z_hat     <- mcclust.ext::minVI(mcclust::comp.psm(out$z))$cl
  psi_hat   <- colMeans(out_relab$psi)
  kappa_hat <- apply(out_relab$kappa, c(2,3), mean)
  eta_hat   <- colMeans(out_relab$eta)
  
  list(out = out_relab,
       z_hat = z_hat,
       psi_hat = psi_hat,
       kappa_hat = kappa_hat,
       eta_hat = eta_hat)
}
# Create header once if file doesn't exist, otherwise append only rows.
.append_rows_locked <- function(path, df, lock_timeout_ms = 60000L) {
  lockfile <- paste0(path, ".lock")
  lk <- filelock::lock(lockfile, timeout = lock_timeout_ms)
  on.exit(filelock::unlock(lk), add = TRUE)
  
  if (!file.exists(path)) {
    data.table::fwrite(df, path, append = FALSE)
  } else {
    # Defensive: order columns as in existing file header
    con <- file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    hdr <- strsplit(readLines(con, n = 1L), ",", fixed = TRUE)[[1]]
    dt  <- as.data.table(df)
    setcolorder(dt, hdr)
    data.table::fwrite(dt, path, append = TRUE)
  }
}

# A single cell of your design
.run_one_cell <- function(K, n, gen_model, diff, rep,
                          n_iter, burn, thin, verbose) {
  # Deterministic per-cell seed
  set.seed(seed_from(K, n, gen_model, diff, rep))
  
  dat <- generate_osbm_data(n = n, K = K,
                            gen_model = gen_model,
                            difficulty = diff)
  A     <- dat$A
  truth <- dat$truth
  
  # Fit three orderings
  ords <- c("NONE","WST","SST")
  fits <- lapply(ords, function(ord) {
    fit_osbm_model(A = A, K = K,
                   ordering = ord,
                   n_iter = n_iter, burn = burn, thin = thin,
                   verbose = verbose)
  })
  names(fits) <- ords
  
  # Metrics per fit (using LOOIC/WAIC; PPL dropped)
  scores <- lapply(fits, function(fit) {
    z_hat <- fit$z_hat
    ari <- fossil::adj.rand.index(z_hat, truth$z)
    vi  <- mcclust::vi.dist(z_hat, truth$z)
    
    mae_eta <- mean(abs(fit$eta_hat - truth$eta))
    mae_kap <- mean(abs(fit$kappa_hat - truth$kappa))
    mae_psi <- mean(abs(fit$psi_hat - truth$psi))
    
    ll_mat <- make_loglik_matrix_osbm(A, fit$out)   # assumed in helpers
    ic <- compute_loo_waic(ll_mat)                  # assumed in helpers
    fit$loo  <- ic$loo_obj
    fit$waic <- ic$waic_obj
    
    psi_draws <- fit$out$psi
    pos_rate  <- mean(apply(psi_draws, 1, function(v) all(v > 0)))
    mono_rate <- mean(apply(psi_draws, 1, function(v) all(diff(v) >= 0)))
    
    c(ari=ari, vi=vi,
      mae_eta=mae_eta, mae_kappa=mae_kap, mae_psi=mae_psi,
      pos_rate=pos_rate, mono_rate=mono_rate,
      looic = ic$looic, waic = ic$waic)
  })
  score_mat <- do.call(rbind, scores)
  
  # Choose best by *lowest* LOOIC; tie-break by WAIC, then deterministic order
  loo_min <- min(score_mat[, "looic"], na.rm = TRUE)
  cand <- which(score_mat[, "looic"] == loo_min)
  if (length(cand) > 1) {
    waic_min <- min(score_mat[cand, "waic"], na.rm = TRUE)
    cand <- cand[score_mat[cand, "waic"] == waic_min]
  }
  chosen_name <- rownames(score_mat)[min(cand)]
  
  df <- data.frame(
    n           = n,
    K           = K,
    gen_model   = gen_model,
    difficulty  = diff,
    rep         = rep,
    fit_model   = rownames(score_mat),
    ari         = score_mat[,"ari"],
    vi          = score_mat[,"vi"],
    mae_eta     = score_mat[,"mae_eta"],
    mae_kappa   = score_mat[,"mae_kappa"],
    mae_psi     = score_mat[,"mae_psi"],
    pos_rate    = score_mat[,"pos_rate"],
    mono_rate   = score_mat[,"mono_rate"],
    looic       = score_mat[,"looic"],
    waic        = score_mat[,"waic"],
    chosen      = (rownames(score_mat) == chosen_name),
    stringsAsFactors = FALSE
  )
  df
}

run_simulation_study_streaming <- function(
    K_vals      = c(3,5,7),
    n_vals      = c(80),
    gen_model   = c('SST','WST','NONE'),
    difficulties= c("dense","medium","sparse"),
    n_reps      = 3,
    n_iter = 5000, burn = 2000, thin = 2,
    seed = 123, verbose = TRUE,
    out_path = "sim_results_stream.csv",
    workers = max(1L, parallel::detectCores() - 1L),
    plan_strategy = c("multisession","multicore"),
    use_progress = TRUE 
) {
  # ---- Guard: output file exists? ----
  if (file.exists(out_path)) {
    ans <- readline(sprintf('File "%s" already exists. Delete and continue? (y/n): ', out_path))
    if (tolower(ans) %in% c("y", "yes")) {
      message("Deleting existing file and proceeding...")
      file.remove(out_path)
      if (file.exists(paste0(out_path, ".lock"))) file.remove(paste0(out_path, ".lock"))
    } else {
      stop("Simulation aborted: output file already exists.")
    }
  }
  set.seed(seed)
  
  plan_strategy <- match.arg(plan_strategy)
  if (plan_strategy == "multicore" && .Platform$OS.type == "windows") {
    warning("multicore not available on Windows; falling back to multisession")
    plan_strategy <- "multisession"
  }
  if (plan_strategy == "multicore") {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::multisession, workers = workers)
  }
  on.exit(future::plan(future::sequential), add = TRUE)
  
  # Build design grid
  grid <- expand.grid(
    K          = K_vals,
    n          = n_vals,
    gen_model  = gen_model,
    diff       = difficulties,
    rep        = seq_len(n_reps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Lay down header once (schema from a tiny dry-run)
  proto <- .run_one_cell(
    K = grid$K[1], n = grid$n[1],
    gen_model = grid$gen_model[1],
    diff = grid$diff[1], rep = grid$rep[1],
    n_iter = 2, burn = 0, thin = 1, verbose = FALSE
  )[0, ]
  .append_rows_locked(out_path, proto)
  
  # Progress handler
  if (use_progress) {
    if (length(progressr::handlers()) == 0L) {
      progressr::handlers(global = TRUE)
      progressr::handlers(progressr::handler_txtprogressbar(clear = FALSE))
    }
  }
  
  # Parallel loop
  runner <- function() {
    p <- if (use_progress) progressr::progressor(steps = nrow(grid)) else NULL
    
    future.apply::future_lapply(
      seq_len(nrow(grid)),
      function(i) {
        g <- grid[i, ]
        message(sprintf("[task %d/%d] K=%d n=%d gen=%s diff=%s rep=%d",
                        i, nrow(grid), g$K, g$n, g$gen_model, g$diff, g$rep))
        
        df_block <- .run_one_cell(
          K = g$K, n = g$n, gen_model = g$gen_model, diff = g$diff, rep = g$rep,
          n_iter = n_iter, burn = burn, thin = thin, verbose = verbose
        )
        
        .append_rows_locked(out_path, df_block)
        
        if (!is.null(p)) {
          p(sprintf("K=%d n=%d gen=%s diff=%s rep=%d",
                    g$K, g$n, g$gen_model, g$diff, g$rep))
        }
        
        best_idx <- which(df_block$chosen)[1]
        list(
          i = i,
          key = g,
          best = df_block$fit_model[best_idx],
          looic_best = df_block$looic[best_idx],
          waic_best  = df_block$waic[best_idx]
        )
      },
      future.seed = TRUE
    )
  }
  
  if (use_progress) {
    progressr::with_progress(runner())
  } else {
    runner()
  }
}

options(progressr.enable = TRUE)
progressr::handlers(progressr::handler_txtprogressbar(clear = FALSE))

# --- Helpers to read and summarise mid/run or after completion ---------------
read_stream_results <- function(path = "sim_results_stream.csv") {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- data.table::fread(path)
  df[, gen_model := factor(gen_model, levels = c("NONE","WST","SST"))]
  df[, fit_model := factor(fit_model, levels = c("NONE","WST","SST"))]
  as.data.frame(df)
}

confusion_from_stream <- function(df) {
  df |>
    dplyr::filter(chosen) |>
    dplyr::count(gen_model, fit_model, name = "count") |>
    tidyr::complete(gen_model, fit_model, fill = list(count = 0))
}
