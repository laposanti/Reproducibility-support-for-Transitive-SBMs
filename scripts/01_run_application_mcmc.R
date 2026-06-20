#!/usr/bin/env Rscript
# Public entry point for bundled application raw fits.

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
if (length(file_arg)) {
  script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg[1L])),
                               winslash = "/", mustWork = TRUE)
  setwd(normalizePath(file.path(dirname(script_path), ".."),
                      winslash = "/", mustWork = TRUE))
}

Sys.setenv(TRANSITIVE_SBM_BUNDLE_ROOT = getwd())
source("scripts/application/run_application_model_fits.R", chdir = FALSE)
