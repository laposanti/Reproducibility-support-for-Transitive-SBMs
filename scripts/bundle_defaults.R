bundle_defaults <- local({
  list(
    bundle_name = "Transitive SBM",
    manuscript_title = "Ordered Stochastic Block Models via Polya-Gamma data augmentation",
    canonical_simulation_results_csv = file.path(
      "output", "simulation", "raw",
      "full_simulation_crossfit_final_DemoKvar_run_20260302_153429.csv"
    ),
    canonical_simulation_driver = file.path(
      "scripts", "simulation", "run_paper_main_simulation_grid.R"
    ),
    canonical_application_run_dir = file.path(
      "output", "application", "raw", "application_run_20260529_110306"
    ),
    canonical_application_driver = file.path(
      "scripts", "application", "run_application_model_fits.R"
    ),
    canonical_postprocessing_dir = file.path(
      "output", "posterior_post_processing", "application_run_20260529_110306"
    )
  )
})
