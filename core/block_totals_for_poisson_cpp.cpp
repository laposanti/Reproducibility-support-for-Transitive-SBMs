#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List block_totals_for_poisson_cpp(
    Rcpp::IntegerVector z,      // length n, labels in {1,...,K} or NA
    Rcpp::NumericVector eta,    // length n
    Rcpp::IntegerVector i_idx,  // length m, 1-based
    Rcpp::IntegerVector j_idx,  // length m, 1-based
    Rcpp::NumericVector N_edge, // length m, dyad totals (R)
    int K
) {
  int n = z.size();
  int m = N_edge.size();
  
  // 1) E_k and eta2k
  Rcpp::NumericVector E_k(K);
  Rcpp::NumericVector eta2k(K);
  
  for (int i = 0; i < n; ++i) {
    int zi = z[i];
    if (zi == NA_INTEGER) continue;
    if (zi < 1 || zi > K) continue;
    int k = zi - 1;
    double ei = eta[i];
    E_k[k]   += ei;
    eta2k[k] += ei * ei;
  }
  
  // 2) Rkl from observed dyads (zeros contribute 0 anyway)
  Rcpp::NumericMatrix Rkl(K, K);
  for (int e = 0; e < m; ++e) {
    int u = i_idx[e] - 1;
    int v = j_idx[e] - 1;
    if (u < 0 || u >= n || v < 0 || v >= n) continue;
    
    int zu = z[u];
    int zv = z[v];
    if (zu == NA_INTEGER || zv == NA_INTEGER) continue;
    if (zu < 1 || zu > K || zv < 1 || zv > K) continue;
    
    int p = (zu <= zv) ? zu : zv;
    int q = (zu <= zv) ? zv : zu;
    int row = p - 1;
    int col = q - 1;
    
    Rkl(row, col) += N_edge[e];
    // no Tkl accumulation here
  }
  
  // 3) Tkl from ALL pairs: outer-product exposure + within-block correction
  Rcpp::NumericMatrix Tkl(K, K);
  
  for (int k = 0; k < K; ++k) {
    for (int l = 0; l < K; ++l) {
      Tkl(k, l) = E_k[k] * E_k[l];
    }
    double diag = 0.5 * (E_k[k] * E_k[k] - eta2k[k]);
    if (diag < 0.0) diag = 0.0; // numeric safety
    Tkl(k, k) = diag;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("Rkl")   = Rkl,
    Rcpp::Named("Tkl")   = Tkl,
    Rcpp::Named("E_k")   = E_k,
    Rcpp::Named("eta2k") = eta2k
  );
}
