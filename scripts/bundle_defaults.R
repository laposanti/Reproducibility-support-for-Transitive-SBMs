# Shared path helpers for the minimal paper-reproduction bundle.

bundle_get_bundle_root <- function() {
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

  script_path <- normalizePath(
    gsub("~\\+~", " ", sub("^--file=", "", file_arg[1L])),
    winslash = "/",
    mustWork = TRUE
  )
  candidate <- dirname(script_path)
  for (root_candidate in unique(c(
    candidate,
    dirname(candidate),
    dirname(dirname(candidate))
  ))) {
    if (dir.exists(file.path(root_candidate, "scripts")) &&
        dir.exists(file.path(root_candidate, "helper_folder"))) {
      return(normalizePath(root_candidate, winslash = "/", mustWork = TRUE))
    }
  }
  normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
}

bundle_latest_subdir <- function(base_dir, prefix = "") {
  if (!dir.exists(base_dir)) return("")

  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[dir.exists(dirs)]
  if (nzchar(prefix)) {
    dirs <- dirs[startsWith(basename(dirs), prefix)]
  }
  if (!length(dirs)) return("")

  dirs <- dirs[order(basename(dirs))]
  dirs[[length(dirs)]]
}

bundle_resolve_application_run_dir <- function(
  run_dir = Sys.getenv("APP_RUN_DIR", unset = ""),
  must_exist = TRUE
) {
  if (!nzchar(run_dir)) {
    run_dir <- bundle_latest_subdir(
      file.path("output", "application", "raw"),
      prefix = "application_run_"
    )
  }

  if (!nzchar(run_dir)) {
    if (must_exist) {
      stop(
        "No application run found. Run `Rscript scripts/01_run_application_mcmc.R` first, ",
        "or set APP_RUN_DIR=output/application/raw/<run_id>."
      )
    }
    return("")
  }

  if (must_exist && !dir.exists(run_dir)) {
    stop("Application run directory does not exist: ", run_dir)
  }

  run_dir
}

bundle_resolve_simulation_results_csv <- function(
  results_path = Sys.getenv("SIM_RESULTS_PATH", unset = ""),
  must_exist = TRUE
) {
  if (!nzchar(results_path)) {
    latest_run <- bundle_latest_subdir(
      file.path("output", "simulation", "raw", "DemoKvar_runs"),
      prefix = "DemoKvar_run_"
    )

    if (nzchar(latest_run)) {
      csvs <- list.files(
        latest_run,
        pattern = "^full_simulation_crossfit_final_.*\\.csv$",
        full.names = TRUE
      )
      if (length(csvs)) {
        csvs <- csvs[order(basename(csvs))]
        results_path <- csvs[[length(csvs)]]
      }
    }
  }

  if (!nzchar(results_path)) {
    if (must_exist) {
      stop(
        "No simulation results CSV found. Run `Rscript scripts/02_run_main_simulation_study.R` first, ",
        "or set SIM_RESULTS_PATH=output/simulation/raw/DemoKvar_runs/<run_id>/full_simulation_crossfit_final_<run_id>.csv."
      )
    }
    return("")
  }

  if (must_exist && !file.exists(results_path)) {
    stop("Simulation results CSV does not exist: ", results_path)
  }

  results_path
}

bundle_defaults <- local({
  app_run_dir <- bundle_resolve_application_run_dir(must_exist = FALSE)
  sim_csv <- bundle_resolve_simulation_results_csv(must_exist = FALSE)

  list(
    bundle_name = "Transitive SBM",
    manuscript_title = "Ordering Stochastic Block Models via prior transitivity",
    canonical_simulation_results_csv = sim_csv,
    canonical_simulation_driver = file.path(
      "scripts", "simulation", "run_paper_main_simulation_grid.R"
    ),
    canonical_application_run_dir = app_run_dir,
    canonical_application_driver = file.path(
      "scripts", "application", "run_application_model_fits.R"
    ),
    canonical_postprocessing_dir = if (nzchar(app_run_dir)) {
      file.path("output", "posterior_post_processing", basename(app_run_dir))
    } else {
      ""
    }
  )
})
