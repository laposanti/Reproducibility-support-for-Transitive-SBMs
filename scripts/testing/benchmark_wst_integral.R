# ======================================================================
# Benchmark: three methods to compute
#   log ∫₀^∞ Binom(A | n, σ(ψ)) · N⁺(ψ; μ₀, σ₀²) dψ
#
# Methods:
#   1. integrate()  — R adaptive quadrature (QUADPACK), treated as ground truth
#   2. Laplace      — mode + Hessian (O(1) arithmetic)
#   3. Gauss–Hermite in C++ via Rcpp (compiled at source time)
#
# Reports: max |error|, mean |error|, speed ratios
# ======================================================================

setwd("/Users/lapo_santi/Desktop/Nial/polya-transitive-sbm/")
library(Rcpp)

# ------------------------------------------------------------------
# 1. Original integrate() version (ground truth)
# ------------------------------------------------------------------
log1pexp <- function(x) ifelse(x > 0, x + log1p(exp(-x)), log1p(exp(x)))
logsig   <- function(x) -log1pexp(-x)
log1m_sig <- function(x) -log1pexp(x)

log_int_binom_truncnorm_integrate <- function(A, n, mu0, sigma0, rel.tol = 1e-12) {
  if (n == 0) return(0)
  logZ <- pnorm(mu0 / sigma0, log.p = TRUE)

  logf <- function(psi) {
    A * logsig(psi) + (n - A) * log1m_sig(psi) +
      dnorm(psi, mean = mu0, sd = sigma0, log = TRUE) - logZ
  }

  ub  <- 40
  opt <- optimize(function(x) -logf(x), interval = c(0, ub))
  m0  <- max(logf(0), logf(opt$minimum), logf(ub))
  g   <- function(psi) exp(logf(psi) - m0)
  out <- integrate(g, lower = 0, upper = Inf, rel.tol = rel.tol)

  if (!is.finite(out$value) || out$value <= 0)
    return(NA_real_)
  log(out$value) + m0
}

# ------------------------------------------------------------------
# 2. Laplace approximation (improved: analytic completion-of-square
#    for the prior part, Laplace only for the Binomial-logistic part)
# ------------------------------------------------------------------
log_int_binom_truncnorm_laplace <- function(A, n, mu0, sigma0) {
  if (n == 0) return(0)
  logZ <- pnorm(mu0 / sigma0, log.p = TRUE)

  logf <- function(psi) {
    A * logsig(psi) + (n - A) * log1m_sig(psi) +
      dnorm(psi, mean = mu0, sd = sigma0, log = TRUE) - logZ
  }

  ub  <- 40
  opt <- optimize(function(x) -logf(x), interval = c(0, ub))
  psi_star <- opt$minimum

  p_star <- 1 / (1 + exp(-psi_star))
  H      <- n * p_star * (1 - p_star) + 1 / sigma0^2
  s_post <- 1 / sqrt(H)    # posterior std dev
  m_post <- psi_star        # posterior mode ≈ posterior mean

  logf_star <- logf(psi_star)

  # Full Laplace on the half-line: Gaussian(m,s²) truncated to [0,∞)
  # log ∫₀^∞ ≈ logf(ψ*) + ½ log(2π/H) + log Φ(m/s) 
  logf_star + 0.5 * log(2 * pi / H) + pnorm(m_post / s_post, log.p = TRUE)
}

# ------------------------------------------------------------------
# 3. Gauss–Hermite quadrature in C++ (Rcpp, compiled inline)
# ------------------------------------------------------------------
sourceCpp(code = '
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// 30-point Gauss–Hermite nodes and weights (physicist convention: ∫ f(x) exp(-x²) dx)
static const int GH_N = 30;
static const double gh_nodes[GH_N] = {
  -6.86334529352990, -6.13827922012394, -5.53314715156751, -4.98891896858994,
  -4.48305535709252, -4.00390860386123, -3.54444387315535, -3.09997052958645,
  -2.66713212453561, -2.24339146776389, -1.82674114360368, -1.41552780019819,
  -1.00833827104672, -0.60392106236773, -0.20112857654803,
   0.20112857654803,  0.60392106236773,  1.00833827104672,  1.41552780019819,
   1.82674114360368,  2.24339146776389,  2.66713212453561,  3.09997052958645,
   3.54444387315535,  4.00390860386123,  4.48305535709252,  4.98891896858994,
   5.53314715156751,  6.13827922012394,  6.86334529352990
};
static const double gh_weights[GH_N] = {
  2.90825470013122e-21, 2.81033360223599e-17, 2.87860708054870e-14, 8.10618629746304e-12,
  9.17858042437852e-10, 5.10852245077594e-08, 1.57909488732608e-06, 2.93872522892298e-05,
  3.48310124318686e-04, 2.73792247488994e-03, 1.47038297048267e-02, 5.51441768702342e-02,
  1.46735847540890e-01, 2.80130930839212e-01, 3.86394889541814e-01,
  3.86394889541814e-01, 2.80130930839212e-01, 1.46735847540890e-01, 5.51441768702342e-02,
  1.47038297048267e-02, 2.73792247488994e-03, 3.48310124318686e-04, 2.93872522892298e-05,
  1.57909488732608e-06, 5.10852245077594e-08, 9.17858042437852e-10, 8.10618629746304e-12,
  2.87860708054870e-14, 2.81033360223599e-17, 2.90825470013122e-21
};

static inline double logsig(double x) {
  return (x > 0.0) ? -std::log1p(std::exp(-x)) : (x - std::log1p(std::exp(x)));
}

// [[Rcpp::export]]
double log_int_binom_truncnorm_gh_cpp(double A, double n, double mu0, double sigma0) {
  if (n == 0.0) return 0.0;

  double logZ = R::pnorm(mu0 / sigma0, 0.0, 1.0, 1, 1);

  // log-integrand on the ORIGINAL scale ψ ∈ [0,∞)
  auto logf = [&](double psi) -> double {
    return A * logsig(psi) + (n - A) * logsig(-psi) +
           R::dnorm(psi, mu0, sigma0, 1) - logZ;
  };

  // --- Find mode on [0, 40] by golden-section search ---
  double lo_s = 0.0, hi_s = 40.0;
  const double gr = 0.6180339887498949;
  double c = hi_s - gr * (hi_s - lo_s);
  double d = lo_s + gr * (hi_s - lo_s);
  for (int iter = 0; iter < 80; ++iter) {
    if (logf(c) > logf(d)) { hi_s = d; d = c; c = hi_s - gr * (hi_s - lo_s); }
    else                    { lo_s = c; c = d; d = lo_s + gr * (hi_s - lo_s); }
    if (hi_s - lo_s < 1e-12) break;
  }
  double psi_star = 0.5 * (lo_s + hi_s);

  // --- Hessian at mode ---
  double p_star = 1.0 / (1.0 + std::exp(-psi_star));
  double H = n * p_star * (1.0 - p_star) + 1.0 / (sigma0 * sigma0);
  double s = 1.0 / std::sqrt(H);

  // --- Change of variable: ψ = exp(u), dψ = exp(u) du ---
  // ∫₀^∞ f(ψ) dψ = ∫_{-∞}^{∞} f(exp(u)) exp(u) du
  // Define g(u) = logf(exp(u)) + u  (the log-transformed integrand)
  // Find mode of g and its curvature, then apply GH centered there.

  auto logg = [&](double u) -> double {
    double psi = std::exp(u);
    if (psi > 1e6) return -1e300;
    return logf(psi) + u;
  };

  // Mode of g(u): start from log(max(psi_star, 0.01))
  double u_init = std::log(std::max(psi_star, 0.01));

  // Golden-section on g in [u_init - 15, u_init + 15]
  double u_lo = u_init - 15.0, u_hi = u_init + 15.0;
  c = u_hi - gr * (u_hi - u_lo);
  d = u_lo + gr * (u_hi - u_lo);
  for (int iter = 0; iter < 80; ++iter) {
    if (logg(c) > logg(d)) { u_hi = d; d = c; c = u_hi - gr * (u_hi - u_lo); }
    else                    { u_lo = c; c = d; d = u_lo + gr * (u_hi - u_lo); }
    if (u_hi - u_lo < 1e-12) break;
  }
  double u_star = 0.5 * (u_lo + u_hi);

  // Numerical second derivative of g at u_star
  double h = 1e-4;
  double g_plus  = logg(u_star + h);
  double g_0     = logg(u_star);
  double g_minus = logg(u_star - h);
  double d2g = (g_plus - 2.0 * g_0 + g_minus) / (h * h);
  double H_u = std::max(-d2g, 1e-8);  // ensure positive curvature
  double s_u = 1.0 / std::sqrt(H_u);

  // GH quadrature: ∫ exp(g(u)) du ≈ (s_u√2) Σ w_i exp(t_i²) exp(g(u_star + √2 s_u t_i))
  double sqrt2 = std::sqrt(2.0);
  double max_log = -1e300;
  double log_contribs[GH_N];

  for (int i = 0; i < GH_N; ++i) {
    double t = gh_nodes[i];
    double u = u_star + sqrt2 * s_u * t;
    log_contribs[i] = std::log(gh_weights[i]) + t * t + logg(u);
    if (log_contribs[i] > max_log) max_log = log_contribs[i];
  }

  double sum_exp = 0.0;
  for (int i = 0; i < GH_N; ++i) {
    if (log_contribs[i] > max_log - 60.0)
      sum_exp += std::exp(log_contribs[i] - max_log);
  }

  return std::log(sqrt2 * s_u) + max_log + std::log(sum_exp);
}
')

# ------------------------------------------------------------------
# 4. Test grid: coverage of the (A, n, mu0, sigma0) space
# ------------------------------------------------------------------
cat("=== Benchmark: log_int_binom_truncnorm ===\n\n")

set.seed(42)
n_tests <- 500

# Structured grid covering edge cases and typical values
test_cases <- data.frame(
  A = integer(0), n = integer(0), mu0 = numeric(0), sigma0 = numeric(0)
)

# (a) Uniform draws from realistic ranges
for (i in 1:200) {
  nn <- sample(1:300, 1)
  AA <- sample(0:nn, 1)
  mm <- runif(1, -1, 5)
  ss <- runif(1, 0.3, 6)
  test_cases <- rbind(test_cases, data.frame(A = AA, n = nn, mu0 = mm, sigma0 = ss))
}

# (b) Small counts (birth moves with sparse data)
for (i in 1:100) {
  nn <- sample(1:10, 1)
  AA <- sample(0:nn, 1)
  mm <- runif(1, 0, 3)
  ss <- runif(1, 0.5, 4)
  test_cases <- rbind(test_cases, data.frame(A = AA, n = nn, mu0 = mm, sigma0 = ss))
}

# (c) Extreme imbalance (A ≈ 0 or A ≈ n)
for (i in 1:100) {
  nn <- sample(10:200, 1)
  AA <- if (runif(1) < 0.5) sample(0:2, 1) else nn - sample(0:2, 1)
  mm <- runif(1, 0, 4)
  ss <- runif(1, 0.5, 5)
  test_cases <- rbind(test_cases, data.frame(A = AA, n = nn, mu0 = mm, sigma0 = ss))
}

# (d) Large sigma0 (vague prior) — calibrate_full often produces these
for (i in 1:100) {
  nn <- sample(1:100, 1)
  AA <- sample(0:nn, 1)
  mm <- runif(1, 0, 2)
  ss <- runif(1, 2, 10)
  test_cases <- rbind(test_cases, data.frame(A = AA, n = nn, mu0 = mm, sigma0 = ss))
}

n_tests <- nrow(test_cases)
cat(sprintf("Testing %d cases\n\n", n_tests))

# ------------------------------------------------------------------
# 5. Accuracy comparison
# ------------------------------------------------------------------
res <- data.frame(
  A = test_cases$A, n = test_cases$n, mu0 = test_cases$mu0, sigma0 = test_cases$sigma0,
  truth = NA_real_, laplace = NA_real_, gh_cpp = NA_real_
)

for (i in seq_len(n_tests)) {
  tc <- test_cases[i, ]
  res$truth[i]   <- log_int_binom_truncnorm_integrate(tc$A, tc$n, tc$mu0, tc$sigma0)
  res$laplace[i] <- log_int_binom_truncnorm_laplace(tc$A, tc$n, tc$mu0, tc$sigma0)
  res$gh_cpp[i]  <- log_int_binom_truncnorm_gh_cpp(tc$A, tc$n, tc$mu0, tc$sigma0)
}

res$err_laplace <- res$laplace - res$truth
res$err_gh      <- res$gh_cpp  - res$truth

cat("--- Accuracy (error = method - integrate) ---\n")
cat(sprintf("Laplace:  max|err|=%.2e  mean|err|=%.2e  median|err|=%.2e  bias=%.2e\n",
            max(abs(res$err_laplace), na.rm = TRUE),
            mean(abs(res$err_laplace), na.rm = TRUE),
            median(abs(res$err_laplace), na.rm = TRUE),
            mean(res$err_laplace, na.rm = TRUE)))
cat(sprintf("GH(C++): max|err|=%.2e  mean|err|=%.2e  median|err|=%.2e  bias=%.2e\n",
            max(abs(res$err_gh), na.rm = TRUE),
            mean(abs(res$err_gh), na.rm = TRUE),
            median(abs(res$err_gh), na.rm = TRUE),
            mean(res$err_gh, na.rm = TRUE)))

# Worst cases
cat("\n--- 5 worst Laplace errors ---\n")
worst_lap <- head(res[order(-abs(res$err_laplace)), ], 5)
print(worst_lap[, c("A", "n", "mu0", "sigma0", "truth", "laplace", "err_laplace")])

cat("\n--- 5 worst GH errors ---\n")
worst_gh <- head(res[order(-abs(res$err_gh)), ], 5)
print(worst_gh[, c("A", "n", "mu0", "sigma0", "truth", "gh_cpp", "err_gh")])

# ------------------------------------------------------------------
# 6. Speed comparison
# ------------------------------------------------------------------
cat("\n--- Speed benchmark (microbenchmark, 100 reps each) ---\n")

if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  install.packages("microbenchmark", repos = "https://cloud.r-project.org")
}
library(microbenchmark)

# Pick a representative medium-difficulty case
A_bm <- 45; n_bm <- 120; mu0_bm <- 1.5; s0_bm <- 3.0

bm <- microbenchmark(
  integrate  = log_int_binom_truncnorm_integrate(A_bm, n_bm, mu0_bm, s0_bm),
  laplace    = log_int_binom_truncnorm_laplace(A_bm, n_bm, mu0_bm, s0_bm),
  gh_cpp     = log_int_binom_truncnorm_gh_cpp(A_bm, n_bm, mu0_bm, s0_bm),
  times = 500L
)
print(bm)

cat("\nSpeedup factors (median time relative to integrate):\n")
med_times <- tapply(bm$time, bm$expr, median)
cat(sprintf("  Laplace / integrate:  %.1fx faster\n", med_times["integrate"] / med_times["laplace"]))
cat(sprintf("  GH(C++) / integrate:  %.1fx faster\n", med_times["integrate"] / med_times["gh_cpp"]))
cat(sprintf("  GH(C++) / Laplace:    %.1fx faster\n", med_times["laplace"]   / med_times["gh_cpp"]))

# Repeat for a batch call (vectorised over all test cases) to measure iteration-level impact
cat("\n--- Batch speed (all test cases × 1 rep) ---\n")
t_int <- system.time({
  for (i in seq_len(n_tests))
    log_int_binom_truncnorm_integrate(test_cases$A[i], test_cases$n[i],
                                      test_cases$mu0[i], test_cases$sigma0[i])
})[["elapsed"]]

t_lap <- system.time({
  for (i in seq_len(n_tests))
    log_int_binom_truncnorm_laplace(test_cases$A[i], test_cases$n[i],
                                    test_cases$mu0[i], test_cases$sigma0[i])
})[["elapsed"]]

t_gh <- system.time({
  for (i in seq_len(n_tests))
    log_int_binom_truncnorm_gh_cpp(test_cases$A[i], test_cases$n[i],
                                   test_cases$mu0[i], test_cases$sigma0[i])
})[["elapsed"]]

cat(sprintf("  integrate:  %.4fs  (%.1f µs/call)\n", t_int, 1e6 * t_int / n_tests))
cat(sprintf("  Laplace:    %.4fs  (%.1f µs/call)\n", t_lap, 1e6 * t_lap / n_tests))
cat(sprintf("  GH(C++):   %.4fs  (%.1f µs/call)\n", t_gh,  1e6 * t_gh  / n_tests))

cat("\n=== Done ===\n")
