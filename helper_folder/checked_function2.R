## ------------------------------------------------------------
## Utilities: truncated Normal > 0 without packages
## ------------------------------------------------------------

rtruncnorm_pos <- function(n, mean, sd) {
  a <- pnorm(0, mean, sd)             # P(Theta <= 0)
  u <- runif(n, min = a, max = 1)     # U ~ Unif(Phi(0), 1)
  qnorm(u, mean, sd)
}

## numerically stable log-mean-exp
logmeanexp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

## ------------------------------------------------------------
## WST half-line integral: analytic and Monte Carlo
## Integral: I = ∫_{0}^{∞} exp(y*psi - 0.5*omega*psi^2) * TN(psi|mu0, sigma0^2) dpsi
## where TN is Normal(mu0,sigma0^2) truncated to (0,∞) with proper normalizer.
## Closed-form:
##   a = omega + 1/sigma0^2
##   b = y + mu0/sigma0^2
##   log I = -log(sigma0) - 0.5*log(a) + b^2/(2a) - mu0^2/(2*sigma0^2)
##           + log Phi(b/sqrt(a)) - log Phi(mu0/sigma0)
## ------------------------------------------------------------

wst_integral_analytic <- function(y, omega, mu0, sigma0, return_log = FALSE) {
  a <- omega + 1/sigma0^2
  b <- y + mu0/sigma0^2
  logI <- -log(sigma0) - 0.5*log(a) + (b^2)/(2*a) - (mu0^2)/(2*sigma0^2) +
    pnorm(b/sqrt(a), log.p = TRUE) - pnorm(mu0/sigma0, log.p = TRUE)
  if (return_log) return(logI) else return(exp(logI))
}

wst_integral_mc <- function(y, omega, mu0, sigma0, nsim = 2e5, seed = 1L, return_log = FALSE) {
  set.seed(seed)
  psi <- rtruncnorm_pos(nsim, mu0, sigma0)
  # Work in log space: log f = y*psi - 0.5*omega*psi^2
  lvals <- y*psi - 0.5*omega*psi^2
  log_est <- logmeanexp(lvals)  # since E[f(psi)] under TN prior (proper) is the integral
  if (return_log) return(log_est) else return(exp(log_est))
}

## ------------------------------------------------------------
## SST increment integral (extreme slots):
## I = ∫_{0}^{∞} exp( (yK - omegaK*psiKm1)*delta - 0.5*omegaK*delta^2 ) * TN(delta|m, v) d delta
## Closed form is identical to WST with substitutions:
##   y <- yK - omegaK*psiKm1, mu0 <- m, sigma0^2 <- v, omega <- omegaK
## ------------------------------------------------------------

sst_increment_analytic <- function(yK, omegaK, psiKm1, m, v, return_log = FALSE) {
  yshift <- yK - omegaK * psiKm1
  wst_integral_analytic(y = yshift, omega = omegaK, mu0 = m, sigma0 = sqrt(v),
                        return_log = return_log)
}

sst_increment_mc <- function(yK, omegaK, psiKm1, m, v, nsim = 2e5, seed = 2L, return_log = FALSE) {
  set.seed(seed)
  delta <- rtruncnorm_pos(nsim, m, sqrt(v))
  lvals <- (yK - omegaK*psiKm1)*delta - 0.5*omegaK*delta^2
  log_est <- logmeanexp(lvals)
  if (return_log) return(log_est) else return(exp(log_est))
}

## ------------------------------------------------------------
## Unconstrained (sanity check): full-line Normal prior (no truncation)
## I = ∫_{R} exp(y*theta - 0.5*omega*theta^2) * N(theta|mu0, sigma0^2) d theta
## log I = -log(sigma0) - 0.5*log(a) + b^2/(2a) - mu0^2/(2*sigma0^2), with a,b as above
## ------------------------------------------------------------

unconstrained_analytic <- function(y, omega, mu0, sigma0, return_log = FALSE) {
  a <- omega + 1/sigma0^2
  b <- y + mu0/sigma0^2
  logI <- -log(sigma0) - 0.5*log(a) + (b^2)/(2*a) - (mu0^2)/(2*sigma0^2)
  if (return_log) return(logI) else return(exp(logI))
}

unconstrained_mc <- function(y, omega, mu0, sigma0, nsim = 2e5, seed = 3L, return_log = FALSE) {
  set.seed(seed)
  th <- rnorm(nsim, mu0, sigma0)
  lvals <- y*th - 0.5*omega*th^2
  # Importance weight is exactly the Normal density; since we sample from the prior,
  # E[f(theta)] is again the integral.
  log_est <- logmeanexp(lvals)
  if (return_log) return(log_est) else return(exp(log_est))
}

## ------------------------------------------------------------
## Quick unit tests / demonstrations
## ------------------------------------------------------------

run_tests <- function() {
  cat("=== WST tests ===\n")
  pars <- list(
    list(y=0.0, omega=0.0, mu0=0.2, s0=1.0),        # should be 1
    list(y=1.3, omega=0.5, mu0=0.0, s0=1.0),
    list(y=-0.7, omega=3.0, mu0=0.5, s0=0.8),
    list(y=2.0, omega=10.0, mu0=0.1, s0=0.5)
  )
  for (p in pars) {
    la <- wst_integral_analytic(p$y, p$omega, p$mu0, p$s0, TRUE)
    lm <- wst_integral_mc(p$y, p$omega, p$mu0, p$s0, nsim = 3e5, seed = 42, return_log = TRUE)
    cat(sprintf("WST: y=%.2f, om=%.2f, mu0=%.2f, s0=%.2f  |  log diff = %.3e  (rel diff ~ %.3e)\n",
                p$y, p$omega, p$mu0, p$s0, la - lm, exp(la - lm) - 1))
  }
  
  cat("\n=== SST increment tests (extreme slots) ===\n")
  sst <- list(
    list(yK=0.0, omegaK=0.0, psiKm1=0.7, m=0.2, v=1.0), # should be 1
    list(yK=1.1, omegaK=0.8, psiKm1=0.4, m=0.1, v=0.7),
    list(yK=-0.3, omegaK=4.0, psiKm1=0.2, m=0.4, v=0.3),
    list(yK=1.8, omegaK=9.0, psiKm1=0.9, m=0.3, v=0.2)
  )
  for (p in sst) {
    la <- sst_increment_analytic(p$yK, p$omegaK, p$psiKm1, p$m, p$v, TRUE)
    lm <- sst_increment_mc(p$yK, p$omegaK, p$psiKm1, p$m, p$v, nsim = 3e5, seed = 777, return_log = TRUE)
    cat(sprintf("SST: yK=%.2f, omK=%.2f, psiK-1=%.2f, m=%.2f, v=%.2f  |  log diff = %.3e  (rel diff ~ %.3e)\n",
                p$yK, p$omegaK, p$psiKm1, p$m, p$v, la - lm, exp(la - lm) - 1))
  }
  
  cat("\n=== Unconstrained sanity checks ===\n")
  upars <- list(
    list(y=0, omega=0, mu0=0.3, s0=1.2),               # should be 1
    list(y=0.7, omega=0.3, mu0=-0.2, s0=0.9),
    list(y=-1.2, omega=6.0, mu0=0.5, s0=0.7)
  )
  for (p in upars) {
    la <- unconstrained_analytic(p$y, p$omega, p$mu0, p$s0, TRUE)
    lm <- unconstrained_mc(p$y, p$omega, p$mu0, p$s0, nsim = 2e5, seed = 99, return_log = TRUE)
    cat(sprintf("UNC: y=%.2f, om=%.2f, mu0=%.2f, s0=%.2f  |  log diff = %.3e  (rel diff ~ %.3e)\n",
                p$y, p$omega, p$mu0, p$s0, la - lm, exp(la - lm) - 1))
  }
}
## ------------------------------------------------------------
## Utilities: truncated Normal > 0 without packages
## ------------------------------------------------------------

rtruncnorm_pos <- function(n, mean, sd) {
  a <- pnorm(0, mean, sd)             # P(Theta <= 0)
  u <- runif(n, min = a, max = 1)     # U ~ Unif(Phi(0), 1)
  qnorm(u, mean, sd)
}

## numerically stable log-mean-exp
logmeanexp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

## =========================
##   FAST integral helpers
## =========================

## Cache for WST half-line integral constants (mu0, sigma0)
make_wst_cache <- function(mu0, sigma0) {
  inv_sig2   <- 1 / (sigma0^2)
  mu_over_s2 <- mu0 * inv_sig2
  const0     <- -log(sigma0) - (mu0^2) / (2 * sigma0^2) - pnorm(mu0 / sigma0, log.p = TRUE)
  list(inv_sig2 = inv_sig2, mu_over_s2 = mu_over_s2, const0 = const0,
       mu0 = mu0, sigma0 = sigma0)
}

## Cache for SST-increment constants (m, v)
make_sst_cache <- function(m, v) {
  inv_v     <- 1 / v
  m_over_v  <- m * inv_v
  const0    <- -0.5 * log(v) - (m^2) / (2 * v) - pnorm(m / sqrt(v), log.p = TRUE)
  list(inv_v = inv_v, m_over_v = m_over_v, const0 = const0,
       m = m, v = v)
}

## Cache for unconstrained Normal prior constants (mu0, sigma0)
make_unconstrained_cache <- function(mu0, sigma0) {
  inv_sig2   <- 1 / (sigma0^2)
  mu_over_s2 <- mu0 * inv_sig2
  const0     <- -log(sigma0) - (mu0^2) / (2 * sigma0^2)
  list(inv_sig2 = inv_sig2, mu_over_s2 = mu_over_s2, const0 = const0,
       mu0 = mu0, sigma0 = sigma0)
}

## ------------------------------------------------------------
## WST half-line integral: analytic, fast, and Monte Carlo
## I = ∫_{0}^{∞} exp(y*psi - 0.5*omega*psi^2) * TN(psi|mu0, sigma0^2) dpsi
## Closed-form:
##   a = omega + 1/sigma0^2
##   b = y + mu0/sigma0^2
##   log I = -log(sigma0) - 0.5*log(a) + b^2/(2a) - mu0^2/(2*sigma0^2)
##           + log Phi(b/sqrt(a)) - log Phi(mu0/sigma0)
## ------------------------------------------------------------

wst_integral_analytic <- function(y, omega, mu0, sigma0, return_log = FALSE) {
  a <- omega + 1/sigma0^2
  b <- y + mu0/sigma0^2
  logI <- -log(sigma0) - 0.5*log(a) + (b^2)/(2*a) - (mu0^2)/(2*sigma0^2) +
    pnorm(b/sqrt(a), log.p = TRUE) - pnorm(mu0/sigma0, log.p = TRUE)
  if (return_log) return(logI) else return(exp(logI))
}

## FAST: vectorized, reusing cache to avoid recomputing constants in mu0,sigma0
wst_integral_fast <- function(y, omega, cache = NULL, mu0 = NULL, sigma0 = NULL, return_log = FALSE) {
  if (is.null(cache)) {
    if (is.null(mu0) || is.null(sigma0)) stop("Provide either `cache` or (`mu0`,`sigma0`).")
    cache <- make_wst_cache(mu0, sigma0)
  }
  a    <- omega + cache$inv_sig2
  b    <- y + cache$mu_over_s2
  logI <- cache$const0 - 0.5 * log(a) + (b^2) / (2 * a) + pnorm(b / sqrt(a), log.p = TRUE)
  if (return_log) logI else exp(logI)
}

wst_integral_mc <- function(y, omega, mu0, sigma0, nsim = 2e5, seed = 1L, return_log = FALSE) {
  set.seed(seed)
  psi <- rtruncnorm_pos(nsim, mu0, sigma0)
  lvals <- y*psi - 0.5*omega*psi^2
  log_est <- logmeanexp(lvals)
  if (return_log) return(log_est) else return(exp(log_est))
}

## ------------------------------------------------------------
## SST increment integral (extreme slots):
## I = ∫_{0}^{∞} exp( (yK - omegaK*psiKm1)*delta - 0.5*omegaK*delta^2 ) * TN(delta|m, v) d delta
## Closed form is identical to WST with substitutions y <- yK - omegaK*psiKm1, mu0 <- m, sigma0^2 <- v, omega <- omegaK
## ------------------------------------------------------------

sst_increment_analytic <- function(yK, omegaK, psiKm1, m, v, return_log = FALSE) {
  yshift <- yK - omegaK * psiKm1
  wst_integral_analytic(y = yshift, omega = omegaK, mu0 = m, sigma0 = sqrt(v),
                        return_log = return_log)
}

## FAST: vectorized, reuse cache(m,v)
sst_increment_fast <- function(yK, omegaK, psiKm1, cache = NULL, m = NULL, v = NULL, return_log = FALSE) {
  if (is.null(cache)) {
    if (is.null(m) || is.null(v)) stop("Provide either `cache` or (`m`,`v`).")
    cache <- make_sst_cache(m, v)
  }
  yshift <- yK - omegaK * psiKm1
  a      <- omegaK + cache$inv_v
  b      <- yshift + cache$m_over_v
  logI   <- cache$const0 - 0.5 * log(a) + (b^2) / (2 * a) + pnorm(b / sqrt(a), log.p = TRUE)
  if (return_log) logI else exp(logI)
}

sst_increment_mc <- function(yK, omegaK, psiKm1, m, v, nsim = 2e5, seed = 2L, return_log = FALSE) {
  set.seed(seed)
  delta <- rtruncnorm_pos(nsim, m, sqrt(v))
  lvals <- (yK - omegaK*psiKm1)*delta - 0.5*omegaK*delta^2
  log_est <- logmeanexp(lvals)
  if (return_log) return(log_est) else return(exp(log_est))
}

## ------------------------------------------------------------
## Unconstrained (sanity check): full-line Normal prior (no truncation)
## I = ∫_{R} exp(y*theta - 0.5*omega*theta^2) * N(theta|mu0, sigma0^2) d theta
## log I = -log(sigma0) - 0.5*log(a) + b^2/(2a) - mu0^2/(2*sigma0^2)
## ------------------------------------------------------------

unconstrained_analytic <- function(y, omega, mu0, sigma0, return_log = FALSE) {
  a <- omega + 1/sigma0^2
  b <- y + mu0/sigma0^2
  logI <- -log(sigma0) - 0.5*log(a) + (b^2)/(2*a) - (mu0^2)/(2*sigma0^2)
  if (return_log) return(logI) else return(exp(logI))
}

## FAST: vectorized with cache(mu0, sigma0)
unconstrained_fast <- function(y, omega, cache = NULL, mu0 = NULL, sigma0 = NULL, return_log = FALSE) {
  if (is.null(cache)) {
    if (is.null(mu0) || is.null(sigma0)) stop("Provide either `cache` or (`mu0`,`sigma0`).")
    cache <- make_unconstrained_cache(mu0, sigma0)
  }
  a    <- omega + cache$inv_sig2
  b    <- y + cache$mu_over_s2
  logI <- cache$const0 - 0.5 * log(a) + (b^2) / (2 * a)
  if (return_log) logI else exp(logI)
}

unconstrained_mc <- function(y, omega, mu0, sigma0, nsim = 2e5, seed = 3L, return_log = FALSE) {
  set.seed(seed)
  th <- rnorm(nsim, mu0, sigma0)
  lvals <- y*th - 0.5*omega*th^2
  log_est <- logmeanexp(lvals)
  if (return_log) return(log_est) else return(exp(log_est))
}

## ------------------------------------------------------------
## Three-way unit tests / demonstrations (analytic vs fast vs MC)
## ------------------------------------------------------------

run_tests <- function() {
  cat("=== WST tests (analytic vs FAST vs MC) ===\n")
  pars <- list(
    list(y=0.0, omega=0.0, mu0=0.2, s0=1.0),        # should be log I = 0
    list(y=1.3, omega=0.5, mu0=0.0, s0=1.0),
    list(y=-0.7, omega=3.0, mu0=0.5, s0=0.8),
    list(y=2.0, omega=10.0, mu0=0.1, s0=0.5)
  )
  for (p in pars) {
    cache <- make_wst_cache(p$mu0, p$s0)
    la <- wst_integral_analytic(p$y, p$omega, p$mu0, p$s0, TRUE)
    lf <- wst_integral_fast     (p$y, p$omega, cache = cache, return_log = TRUE)
    lm <- wst_integral_mc       (p$y, p$omega, p$mu0, p$s0, nsim = 3e5, seed = 42, return_log = TRUE)
    cat(sprintf("WST: y=%.2f, om=%.2f, mu0=%.2f, s0=%.2f  |  fast-analytic=%.2e,  analytic-MC=%.2e\n",
                p$y, p$omega, p$mu0, p$s0, lf - la, la - lm))
  }
  
  cat("\n=== SST increment tests (analytic vs FAST vs MC) ===\n")
  sst <- list(
    list(yK=0.0, omegaK=0.0, psiKm1=0.7, m=0.2, v=1.0), # should be log I = 0
    list(yK=1.1, omegaK=0.8, psiKm1=0.4, m=0.1, v=0.7),
    list(yK=-0.3, omegaK=4.0, psiKm1=0.2, m=0.4, v=0.3),
    list(yK=1.8, omegaK=9.0, psiKm1=0.9, m=0.3, v=0.2)
  )
  for (p in sst) {
    cache <- make_sst_cache(p$m, p$v)
    la <- sst_increment_analytic(p$yK, p$omegaK, p$psiKm1, p$m, p$v, TRUE)
    lf <- sst_increment_fast     (p$yK, p$omegaK, p$psiKm1, cache = cache, return_log = TRUE)
    lm <- sst_increment_mc       (p$yK, p$omegaK, p$psiKm1, p$m, p$v, nsim = 3e5, seed = 777, return_log = TRUE)
    cat(sprintf("SST: yK=%.2f, omK=%.2f, psiK-1=%.2f, m=%.2f, v=%.2f  |  fast-analytic=%.2e,  analytic-MC=%.2e\n",
                p$yK, p$omegaK, p$psiKm1, p$m, p$v, lf - la, la - lm))
  }
  
  cat("\n=== Unconstrained sanity checks (analytic vs FAST vs MC) ===\n")
  upars <- list(
    list(y=0, omega=0, mu0=0.3, s0=1.2),               # should be log I = 0
    list(y=0.7, omega=0.3, mu0=-0.2, s0=0.9),
    list(y=-1.2, omega=6.0, mu0=0.5, s0=0.7)
  )
  for (p in upars) {
    cache <- make_unconstrained_cache(p$mu0, p$s0)
    la <- unconstrained_analytic(p$y, p$omega, p$mu0, p$s0, TRUE)
    lf <- unconstrained_fast     (p$y, p$omega, cache = cache, return_log = TRUE)
    lm <- unconstrained_mc       (p$y, p$omega, p$mu0, p$s0, nsim = 2e5, seed = 99, return_log = TRUE)
    cat(sprintf("UNC: y=%.2f, om=%.2f, mu0=%.2f, s0=%.2f  |  fast-analytic=%.2e,  analytic-MC=%.2e\n",
                p$y, p$omega, p$mu0, p$s0, lf - la, la - lm))
  }
}

## Run the tests
run_tests()

## Run the tests
run_tests()


