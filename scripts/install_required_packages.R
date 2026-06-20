#!/usr/bin/env Rscript
# Install the R packages required by the minimal paper-reproduction workflow.

required_packages <- c(
  "BayesLogit",
  "Matrix",
  "Rcpp",
  "truncnorm",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "ggplot2",
  "patchwork",
  "viridis",
  "igraph",
  "ggraph",
  "graphlayouts",
  "mcclust",
  "mcclust.ext",
  "salso",
  "loo",
  "coda",
  "fossil",
  "scales",
  "fs",
  "glue",
  "stringr",
  "kableExtra",
  "knitr",
  "purrr",
  "lpSolve",
  "reshape2",
  "ggside",
  "plotly",
  "htmlwidgets",
  "fs",
  "rmarkdown"
)
missing <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (!length(missing)) {
  cat("All requested packages are already installed.\n")
  quit(status = 0)
}

cat("Installing", length(missing), "packages:\n")
cat(paste0(" - ", missing, collapse = "\n"), "\n")
install.packages(missing)
