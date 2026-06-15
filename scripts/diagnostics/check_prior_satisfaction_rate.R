## =========================================================
##  WST / MST (global) and SST (via permutation) checks
##  Lapo — corrected version with witnesses + plot
## =========================================================

set.seed(123)

# ---------- Matrix generator ----------
sample_M <- function(K) {
  M <- matrix(0.5, K, K)
  for (i in 1:(K-1)) for (j in (i+1):K) {
    u <- runif(1)
    M[i,j] <- u; M[j,i] <- 1 - u
  }
  M
}

# ---------- Global WST/MST checks (order-invariant) ----------
# WST: for all ordered triples (i,j,k), if M[i,j]>=1/2 and M[j,k]>=1/2 then M[i,k]>=1/2
check_WST_global <- function(M) {
  K <- nrow(M)
  for (i in 1:K) for (j in 1:K) for (k in 1:K) {
    if (length(unique(c(i,j,k))) < 3) next
    if (M[i,j] >= 0.5 && M[j,k] >= 0.5 && M[i,k] < 0.5) {
      return(list(ok = FALSE, witness = list(i=i,j=j,k=k,
                                             rij=M[i,j], rjk=M[j,k], rik=M[i,k])))
    }
  }
  list(ok = TRUE, witness = NULL)
}

# MST: if M[i,j]>=1/2 and M[j,k]>=1/2 then M[i,k] >= min(M[i,j], M[j,k])
check_MST_global <- function(M) {
  K <- nrow(M)
  for (i in 1:K) for (j in 1:K) for (k in 1:K) {
    if (length(unique(c(i,j,k))) < 3) next
    if (M[i,j] >= 0.5 && M[j,k] >= 0.5) {
      thr <- min(M[i,j], M[j,k])
      if (M[i,k] < thr) {
        return(list(ok = FALSE, witness = list(i=i,j=j,k=k,
                                               rij=M[i,j], rjk=M[j,k], rik=M[i,k], thr=thr)))
      }
    }
  }
  list(ok = TRUE, witness = NULL)
}

# ---------- SST via permutation (bivariate isotonic) ----------
# After permuting by pi, require: for all i<j, for all k:
# row monotonicity: P[i,k] >= P[j,k];  column monotonicity: P[k,i] <= P[k,j]
satisfies_SST_perm <- function(M, pi) {
  P <- M[pi, pi, drop=FALSE]
  K <- nrow(P)
  if (K < 2) return(TRUE)
  for (i in 1:(K-1)) for (j in (i+1):K) {
    # rows nonincreasing
    if (any(P[i, ] < P[j, ] - 1e-12)) return(FALSE)
    # columns nondecreasing
    if (any(P[, i] > P[, j] + 1e-12)) return(FALSE)
  }
  TRUE
}

# Generate all permutations (factorial; OK up to ~K=7)
all_perms <- function(K) {
  perms <- list(1L)
  if (K == 1L) return(perms)
  for (m in 2:K) {
    new_perms <- vector("list", length(perms) * m)
    idx <- 1L
    for (p in perms) for (pos in 0:(m-1)) {
      new_perms[[idx]] <- append(p, m, after = pos)
      idx <- idx + 1L
    }
    perms <- new_perms
  }
  lapply(perms, as.integer)
}

exists_SST <- function(M) {
  for (pi in all_perms(nrow(M))) {
    if (satisfies_SST_perm(M, pi)) return(list(ok=TRUE, pi=pi))
  }
  list(ok=FALSE, pi=NULL)
}

# ---------- Experiment + plotting ----------
run_experiment <- function(K_vals = 3:6, nrep = 200) {
  summary <- data.frame()
  witnesses <- list()
  for (K in K_vals) {
    cat(sprintf("\n=== K = %d ===\n", K))
    cnt <- c(WST=0L, MST=0L, SST=0L)
    wit <- list(WST=NULL, MST=NULL, SST=NULL)
    for (r in 1:nrep) {
      M <- sample_M(K)
      wst <- check_WST_global(M)
      mst <- check_MST_global(M)
      sst <- exists_SST(M)
      
      if (wst$ok) cnt["WST"] <- cnt["WST"] + 1L else if (is.null(wit$WST)) wit$WST <- list(M=M, witness=wst$witness)
      if (mst$ok) cnt["MST"] <- cnt["MST"] + 1L else if (is.null(wit$MST)) wit$MST <- list(M=M, witness=mst$witness)
      if (sst$ok) cnt["SST"] <- cnt["SST"] + 1L else if (is.null(wit$SST)) wit$SST <- list(M=M, pi=NULL)
      
      if (r %% 50 == 0) cat(sprintf("  ...rep %d/%d\r", r, nrep))
    }
    rat <- 100 * cnt / nrep
    cat(sprintf("WST: %.1f%%   MST: %.1f%%   SST: %.1f%%\n", rat["WST"], rat["MST"], rat["SST"]))
    summary <- rbind(summary,
                     data.frame(K=K, regime=factor(names(cnt), levels=c("WST","MST","SST")),
                                percent=as.numeric(rat)))
    witnesses[[as.character(K)]] <- wit
  }
  list(summary=summary, witnesses=witnesses)
}

# ---- run
res <- run_experiment(K_vals = 3:6, nrep = 150)

# ---- plot
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
ggplot(res$summary, aes(x=K, y=percent, color=regime)) +
  geom_point(size=3) + geom_line() +
  labs(x="K", y="Satisfaction rate (%)", color="Regime",
       title="") +
  theme_minimal(base_size = 12)

# ---- print one witness if any (for K=6)
wit6 <- res$witnesses[["6"]]
if (!is.null(wit6$WST)) {
  cat("\n[WST witness for K=6]\n")
  print(round(wit6$WST$M, 3))
  print(wit6$WST$witness)
}
if (!is.null(wit6$MST)) {
  cat("\n[MST witness for K=6]\n")
  print(round(wit6$MST$M, 3))
  print(wit6$MST$witness)
}
