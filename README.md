# Reproducing The Paper Results

This folder is a self-contained reproduction bundle for the manuscript figures,
tables, and cached numerical results used by `ArXiv Submission.tex`.

Run commands from this folder:

```sh
cd "tex file/reproduce_paper_results"
```

The bundle includes cached MCMC results, post-processing outputs, CSV summaries,
data files, and regenerated plots. This lets a reader rebuild every plot without
rerunning the expensive MCMC. Full rerun entry points are also included.

## Requirements

Use R with the packages used by the project scripts, in particular:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "readr", "patchwork", "viridis",
  "igraph", "ggraph", "graphlayouts", "mcclust", "mcclust.ext",
  "salso", "loo", "coda", "truncnorm", "fossil", "Matrix",
  "BayesLogit", "Rcpp", "lpSolve", "scales", "fs"
))
```

Some scripts run without optional packages by using fallbacks, but the list above
matches the environment used to create the bundled results.

## Folder Contents

- `data/`: all empirical network data used by the application analysis.
- `core/`, `helper_folder/`: model fitting and helper functions.
- `scripts/`: renamed entry-point scripts plus copied helper script folders.
- `output/application/raw/application_run_20260529_110306/`: cached application MCMC fits and CSV summaries.
- `output/posterior_post_processing/application_run_20260529_110306/`: canonical posterior partitions, diagnostics, and K summaries used by plots.
- `output/simulation/raw/`: simulation result CSVs used for ARI/VI figures.
- `output/simulation/plots/2026-06-15_sst_labels/`: regenerated simulation plots with `SST` labels.
- `output/paper/figures/application_run_20260529_110306_sst_labels/`: regenerated application figures with `SST` labels.
- `paper_figures/all_figures/`: exact figure files copied into the manuscript figure folder.
- `paper_figures/updated plots/`: collected updated application and simulation plots.

## Renamed Entry Scripts

- `scripts/01_run_application_mcmc.R`: full application MCMC run.
- `scripts/02_run_main_simulation_study.R`: full main simulation study.
- `scripts/03_build_application_postprocessing_cube.R`: builds canonical z-hat/K/diagnostic cube from application fits.
- `scripts/04_build_paper_tables.R`: rebuilds paper tables from cached application results.
- `scripts/05_plot_paper_application_figures.R`: rebuilds application network figures and BT diagnostic plot.
- `scripts/06_plot_simulation_recovery_figures.R`: rebuilds ARI/VI simulation figures.
- `scripts/07_plot_support_geometry.R`: rebuilds support-geometry diagnostics.
- `scripts/08_plot_mirrored_ocrp_diagnostics.R`: rebuilds mirrored OCRP prior diagnostics.
- `scripts/09_build_bradley_terry_delta_plot.R`: rebuilds the BT additivity diagnostic only.
- `scripts/10_plot_prior_satisfaction_rate.R`: archived diagnostic script for prior satisfaction checks.

The copied `scripts/analysis`, `scripts/application`, `scripts/simulation`,
`scripts/diagnostics`, and `scripts/testing` folders are dependencies for these
entry points.

## Fast Plot Reproduction

### Application Network Figures

These commands regenerate the application plots used in the manuscript without
refitting models. The output folder name is intentionally separate from the
cached manuscript copy.

```sh
APP_RUN_DIR=output/application/raw/application_run_20260529_110306 \
APP_PAPER_FIGURES_DIR=output/paper/figures/rebuilt_application_sst_labels \
APP_PAPER_TABLES_DIR=output/paper/tables/rebuilt_application_sst_labels \
Rscript scripts/05_plot_paper_application_figures.R
```

This reproduces:

- `moreno_sheep_SST_network_tier_line.png`
- `moreno_sheep_DCSBM_network_tier_line.png`
- `<dataset>_combined_block_networks_clean.pdf/png` for all six datasets
- `bt_delta_wst_applications.pdf/png`

The manuscript currently uses the copies in `paper_figures/all_figures/`.

### Simulation ARI/VI Figures

```sh
SIM_RESULTS_PATH=output/simulation/raw/full_simulation_crossfit_final_DemoKvar_run_20260302_153429.csv \
SIM_PLOTS_OUTPUT_DIR=output/simulation/plots/rebuilt_sst_labels \
Rscript scripts/06_plot_simulation_recovery_figures.R
```

This reproduces:

- `vi_boxplot_WST_gen.pdf/png`
- `vi_boxplot_SST_gen.pdf/png`
- `vi_boxplot_combined.pdf/png`
- `ari_boxplot_combined.pdf/png`
- `ari_boxplot_all.*`, `vi_boxplot_all.*`, and summary CSVs

### Support Geometry Figure

```sh
Rscript scripts/07_plot_support_geometry.R
```

The paper uses `support_3d_shaded_geometry.png`, cached in
`paper_figures/all_figures/` and `output/diagnostics/support_geometry/`.

### Mirrored OCRP Diagnostic Figures

```sh
Rscript scripts/08_plot_mirrored_ocrp_diagnostics.R
```

This rebuilds the OCRP prior figures under `output/simulation/ocrp_tests_mirrored/`,
including the K prior, position profile, end-block, max-block, and equal-size
diagnostics used by the supplement.

## Rebuilding Tables

From cached application fits:

```sh
APP_RUN_DIR=output/application/raw/application_run_20260529_110306 \
Rscript scripts/04_build_paper_tables.R
```

This rebuilds model-selection, hierarchy-diagnostic, violation-rate, and
application supplement tables under `output/paper/tables/application_run_20260529_110306/`.

To rebuild only the Bradley-Terry additivity diagnostic:

```sh
APP_RUN_DIR=output/application/raw/application_run_20260529_110306 \
APP_PAPER_TABLES_DIR=output/paper/tables/rebuilt_bt_delta \
APP_PAPER_FIGURES_DIR=output/paper/figures/rebuilt_bt_delta \
Rscript scripts/09_build_bradley_terry_delta_plot.R
```

## Full Reruns

Full reruns are computationally expensive and may not reproduce byte-identical
MCMC output unless the same R version, package versions, and random-number
streams are used.

Application MCMC:

```sh
APP_N_ITER=10000 APP_BURN=3000 APP_THIN=2 APP_SEED=42 \
Rscript scripts/01_run_application_mcmc.R
```

Main simulation study:

```sh
Rscript scripts/02_run_main_simulation_study.R
```

After a full application rerun, use the new run directory in:

```sh
APP_RUN_DIR=output/application/raw/<new_application_run_id> \
Rscript scripts/03_build_application_postprocessing_cube.R
```

Then rerun the table and plot scripts above with the same `APP_RUN_DIR`.

## Figure-To-File Map

- Main simulation recovery figure: `scripts/06_plot_simulation_recovery_figures.R`
  from `output/simulation/raw/full_simulation_crossfit_final_DemoKvar_run_20260302_153429.csv`.
- Main bighorn sheep tier-line figure: `scripts/05_plot_paper_application_figures.R`
  from `output/application/raw/application_run_20260529_110306/` and
  `output/posterior_post_processing/application_run_20260529_110306/`.
- Application block-flow figures: same application plot script and cached run.
- Support geometry figure: `scripts/07_plot_support_geometry.R`.
- BT additivity diagnostic: `scripts/09_build_bradley_terry_delta_plot.R`, also
  called by `scripts/05_plot_paper_application_figures.R`.
- Mirrored OCRP prior figures: `scripts/08_plot_mirrored_ocrp_diagnostics.R`.
- Prior satisfaction-rate figure: cached as
  `paper_figures/all_figures/satisfaction_rate.png`.

## Label Update Note

The application and simulation plot scripts in this bundle use `SST` in figure
labels after the geometric section of the paper. The old `Toeplitz SST` labels
are not used in the regenerated application plots.
