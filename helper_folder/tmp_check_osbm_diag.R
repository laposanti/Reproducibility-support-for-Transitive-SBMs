source('./helper_folder/transitivity_check_helper.R')
A <- readRDS('posterior_predictive_checks/high_school/high_school_A_obs.rds')
sst <- readRDS('posterior_predictive_checks/high_school/SST/high_school_SST_fit.rds')
wst <- readRDS('posterior_predictive_checks/high_school/WST/high_school_WST_fit.rds')

z_sst <- if (is.matrix(sst$z)) sst$z else do.call(rbind, sst$z)
z_wst <- if (is.matrix(wst$z)) wst$z else do.call(rbind, wst$z)

zhat_sst <- salso::salso(z_sst, loss = salso::VI(), nRuns = 1L, nCores = 1L)
zhat_wst <- salso::salso(z_wst, loss = salso::VI(), nRuns = 1L, nCores = 1L)

d1 <- summarise_osbm_diagnostics(
  sst,
  regime = 'SST',
  K_max_hint = length(unique(zhat_sst)),
  z_hat = zhat_sst,
  n = nrow(A),
  A = A,
  method_order = 'bt',
  T_block = NULL
)

d2 <- summarise_osbm_diagnostics(
  wst,
  regime = 'WST',
  K_max_hint = length(unique(zhat_wst)),
  z_hat = zhat_wst,
  n = nrow(A),
  A = A,
  method_order = 'bt',
  T_block = NULL
)

cat('SST violation_rate_zhat=', d1['violation_rate_zhat'], '\n', sep = '')
cat('WST violation_rate_zhat=', d2['violation_rate_zhat'], '\n', sep = '')
