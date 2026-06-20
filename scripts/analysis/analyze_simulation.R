cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]),
                               winslash = "/", mustWork = TRUE)
  setwd(normalizePath(file.path(dirname(script_path), "../.."),
                      winslash = "/", mustWork = TRUE))
}

Sys.setenv(TRANSITIVE_SBM_BUNDLE_ROOT = getwd())
source("scripts/analysis/build_simulation_tables_and_figures.R", chdir = FALSE)
