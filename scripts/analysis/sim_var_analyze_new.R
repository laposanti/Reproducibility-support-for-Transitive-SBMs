cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]),
                               winslash = "/", mustWork = TRUE)
  setwd(normalizePath(file.path(dirname(script_path), "../.."),
                      winslash = "/", mustWork = TRUE))
}

source("scripts/analysis/build_simulation_crossfit_tables.R", chdir = FALSE)
