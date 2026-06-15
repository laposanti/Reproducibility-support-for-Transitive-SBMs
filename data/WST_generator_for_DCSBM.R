logistic <- function(x) 1/(1+exp(-x))

simulate_wst_network <- function(z, eta, kappa, psi) {
  n <- length(z)
  A <- matrix(0L, n, n)
  N <- matrix(0L, n, n)
  
  for (i in 1:(n-1)) for (j in (i+1):n) {
    ki <- z[i]; kj <- z[j]
    rate <- eta[i] * eta[j] * kappa[ki, kj]
    Nij <- rpois(1L, rate)
    if (Nij == 0L) next
    
    p <- logistic(psi[ki, kj])
    wij <- rbinom(1L, Nij, p)
    A[i,j] <- wij
    A[j,i] <- Nij - wij
    N[i,j] <- Nij
    N[j,i] <- Nij
  }
  list(A=A, N=N)
}

make_wst_dataset <- function(
    K_true = 10,
    block_sizes = c(18,14,12,10,9,9,8,7,7,6),
    games_per_node_target = 25,   # <-- key knob
    within_mult = 3,              # more assortativity helps K recovery
    delta = 0.25,                 # step in psi by rank distance
    psi_cap = 2.0,                # don’t cap too early
    eta_sd = 0.05,
    seed = 1
) {
  set.seed(seed)
  stopifnot(length(block_sizes)==K_true)
  n <- sum(block_sizes)
  z_true <- rep(seq_len(K_true), times=block_sizes)
  
  eta_true <- pmax(0.3, rnorm(n, 1, eta_sd))
  
  # Choose m0_between so that the *largest* block gets about games_per_node_target
  nk_max <- max(block_sizes)
  # rough expected games per node ≈ (K-1 + 2*within_mult) * m0_between / n_k
  m0_between <- games_per_node_target * nk_max / ( (K_true-1) + 2*within_mult )
  
  kappa <- matrix(0, K_true, K_true)
  for (k in 1:K_true) for (l in 1:K_true) {
    nk <- block_sizes[k]; nl <- block_sizes[l]
    if (k == l) {
      m_within <- within_mult * m0_between
      kappa[k,l] <- 2 * m_within / max(1, nk*(nk-1))
    } else {
      kappa[k,l] <- m0_between / (nk*nl)
    }
  }
  kappa <- (kappa + t(kappa))/2
  
  # WST psi: ordered, increasing with distance, with mild jitter
  psi_true <- matrix(0, K_true, K_true)
  for (k in 1:K_true) for (l in (k+1):K_true) {
    base <- pmin(delta*(l-k), psi_cap)
    val <- base + rnorm(1, 0, 0.03)
    psi_true[k,l] <- val
    psi_true[l,k] <- -val
  }
  diag(psi_true) <- 0
  
  sim <- simulate_wst_network(z_true, eta_true, kappa, psi_true)
  
  cat("n =", n, "K_true =", K_true, "\n")
  cat("Median games/node:", median(rowSums(sim$N)), "\n")
  cat("5%-95% games/node:", quantile(rowSums(sim$N), c(.05,.95)), "\n")
  
  list(
    A = sim$A, N = sim$N,
    z_true = z_true,
    eta_true = eta_true,
    kappa_true = kappa,
    psi_true = psi_true,
    block_sizes = block_sizes,
    m0_between = m0_between
  )
}

toy <- make_wst_dataset(seed = 123)
A <- toy$A
