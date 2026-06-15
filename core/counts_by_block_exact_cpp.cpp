#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List counts_by_block_exact_cpp(
    int i,
    const NumericMatrix& A,
    const IntegerVector& z,
    const IntegerVector& i_idx,
    const IntegerVector& j_idx,
    const NumericVector& N_edge,   // <- accept numeric (matches your R)
    const List& edge_by_node,
    int K
) {
  int n = z.size();
  if (i < 1 || i > n) stop("counts_by_block_exact_cpp: i out of range.");
  if (A.nrow() != n || A.ncol() != n) stop("counts_by_block_exact_cpp: A dims != length(z).");
  if (K < 1) stop("counts_by_block_exact_cpp: K < 1.");
  if (edge_by_node.size() != n) stop("counts_by_block_exact_cpp: edge_by_node wrong length.");
  
  NumericVector c_plus(K);
  NumericVector N_tot(K);
  
  IntegerVector e_list = edge_by_node[i - 1];
  int L = e_list.size();
  if (L == 0) {
    return List::create(_["c_plus"] = c_plus, _["N_tot"] = N_tot);
  }
  
  int m = N_edge.size();
  for (int idx = 0; idx < L; ++idx) {
    int e1 = e_list[idx];
    if (e1 < 1 || e1 > m) stop("counts_by_block_exact_cpp: edge index out of range.");
    int e = e1 - 1;
    
    int u = i_idx[e];
    int v = j_idx[e];
    if (u < 1 || u > n || v < 1 || v > n) stop("counts_by_block_exact_cpp: i_idx/j_idx out of range.");
    
    int j = (i == u) ? v : u;
    int ell = z[j - 1];
    if (ell < 1 || ell > K) stop("counts_by_block_exact_cpp: z has label outside 1..K.");
    
    double n_e = N_edge[e];
    double a_uv = A(u - 1, v - 1);
    double a_i_to_j = (i == u) ? a_uv : (n_e - a_uv);
    
    c_plus[ell - 1] += a_i_to_j;
    N_tot [ell - 1] += n_e;
  }
  
  return List::create(_["c_plus"] = c_plus, _["N_tot"] = N_tot);
}
