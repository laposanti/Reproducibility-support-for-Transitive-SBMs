# Ordered SBMs Reproduction Bundle

Minimal reproduction repository for the paper on ordered stochastic block models with transitivity-constrained priors.

The repository now contains only:

- the samplers and helper code needed for the paper-facing workflow
- the six application datasets used in the paper
- the public scripts needed to regenerate simulations, application fits, tables, and figures
- a small preview set for this README

No cached raw runs are bundled. The only committed rendered assets are the four preview images in `docs/previews/`. All real results are regenerated locally under `output/`, which is ignored by Git.

## Preview

| Support Geometry | Simulation Recovery | Application Figure | Bradley-Terry Delta |
| --- | --- | --- | --- |
| ![Support geometry](docs/previews/support_geometry.png) | ![Simulation VI](docs/previews/simulation_vi.png) | ![Application network](docs/previews/application_network.png) | ![BT delta](docs/previews/bt_delta.png) |

## Repository Layout

```text
core/
  transitive_sbm_sampler.R
  DCSBM_varK.R
  ppc_checks.R

helper_folder/
  load_sampler_helpers.R
  config/
  diagnostics/
  io/
  models/
  simulation/

data/
  mountain_goats/adjacency_matrix.csv
  citations_data/adjacency_matrix.csv
  macaques_data/edge_list.tsv
  high_school/edges.csv
  moreno_sheep/edges.csv
  Strauss_2019b/edges.csv

scripts/
  01_run_application_mcmc.R
  02_run_main_simulation_study.R
  03_build_application_postprocessing_cube.R
  04_build_paper_tables.R
  05_plot_paper_application_figures.R
  06_build_simulation_tables_and_figures.R
  07_plot_support_geometry.R
  09_build_bradley_terry_delta_plot.R
```

## Installation

Install the required R packages:

```sh
Rscript scripts/install_required_packages.R
```

## Quick Validation

Run the two lightweight checks first:

```sh
Rscript scripts/testing/test_hierarchy_metrics.R
Rscript scripts/testing/quick_smoke_test.R
```

Then use the fast smoke paths to confirm that the cleaned repository runs from a blank `output/` tree:

```sh
DEMOKVAR_SMOKE=1 Rscript scripts/02_run_main_simulation_study.R
Rscript scripts/06_build_simulation_tables_and_figures.R
env APP_DATASETS=moreno_sheep APP_N_ITER=800 APP_BURN=200 APP_THIN=2 APP_SEED=1 Rscript scripts/01_run_application_mcmc.R
Rscript scripts/03_build_application_postprocessing_cube.R
Rscript scripts/04_build_paper_tables.R
Rscript scripts/05_plot_paper_application_figures.R
Rscript scripts/07_plot_support_geometry.R
Rscript scripts/09_build_bradley_terry_delta_plot.R
```

The `moreno_sheep` application smoke run is only a validation shortcut. For the full paper-facing application rebuild, run `scripts/01_run_application_mcmc.R` with its default six-dataset configuration.

## Public Workflow

### 1. Application Study

Run the application fits:

```sh
Rscript scripts/01_run_application_mcmc.R
```

This creates a new directory:

```text
output/application/raw/application_run_<timestamp>/
```

The post-processing and paper builders automatically use the latest application run when `APP_RUN_DIR` is not set. You can still force a specific run with `APP_RUN_DIR=...`.

To run a smaller subset during validation or debugging, pass `APP_DATASETS=<comma-separated names>`.

Build the canonical post-processing cube:

```sh
Rscript scripts/03_build_application_postprocessing_cube.R
```

Build the paper tables:

```sh
Rscript scripts/04_build_paper_tables.R
```

Build the paper figures:

```sh
Rscript scripts/05_plot_paper_application_figures.R
```

Rebuild only the Bradley-Terry delta summary:

```sh
Rscript scripts/09_build_bradley_terry_delta_plot.R
```

### 2. Simulation Study

Run the simulation grid:

```sh
Rscript scripts/02_run_main_simulation_study.R
```

Notes:

- the default reproduction setting now uses `3` replicates
- override with `DEMOKVAR_N_REP=<n>` if you want a different replicate count
- `DEMOKVAR_SMOKE=1` runs the tiny smoke grid

This creates a new directory:

```text
output/simulation/raw/DemoKvar_runs/DemoKvar_run_<timestamp>/
```

Build the simulation tables and figures from the latest simulation run:

```sh
Rscript scripts/06_build_simulation_tables_and_figures.R
```

Or target a specific results CSV:

```sh
SIM_RESULTS_PATH=output/simulation/raw/DemoKvar_runs/<run_id>/full_simulation_crossfit_final_<run_id>.csv \
Rscript scripts/06_build_simulation_tables_and_figures.R
```

### 3. Support Geometry Figure

```sh
Rscript scripts/07_plot_support_geometry.R
```

## Outputs

Generated files are written under `output/`:

- `output/application/raw/`
- `output/posterior_post_processing/`
- `output/paper/tables/`
- `output/paper/figures/`
- `output/simulation/raw/`
- `output/simulation/tables/`
- `output/simulation/plots/`
- `output/diagnostics/support_geometry/`

These outputs are not committed. A fresh clone of the repository should have an effectively empty `output/` directory except for `.gitkeep`.

## Verified Commands

The following commands were rerun successfully after the cleanup:

```sh
Rscript scripts/testing/test_hierarchy_metrics.R
Rscript scripts/testing/quick_smoke_test.R
DEMOKVAR_SMOKE=1 Rscript scripts/02_run_main_simulation_study.R
Rscript scripts/06_build_simulation_tables_and_figures.R
env APP_DATASETS=moreno_sheep APP_N_ITER=800 APP_BURN=200 APP_THIN=2 APP_SEED=1 Rscript scripts/01_run_application_mcmc.R
Rscript scripts/03_build_application_postprocessing_cube.R
Rscript scripts/04_build_paper_tables.R
Rscript scripts/05_plot_paper_application_figures.R
Rscript scripts/07_plot_support_geometry.R
Rscript scripts/09_build_bradley_terry_delta_plot.R
```

The public entry-point scripts resolve the repository root from their own path, so they can be launched either from the repository root or via an absolute script path from another working directory.
