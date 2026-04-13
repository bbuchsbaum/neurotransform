#include <Rcpp.h>
#include <limits>
#include <map>
#include <vector>

//' Fast triplet to dgCMatrix assembly with duplicate aggregation
//'
//' Converts triplet format (i, j, x) to sparse column matrix.
//'
//' @param i 1-based row indices
//' @param j 1-based col indices
//' @param x values
//' @param nrow number of rows
//' @param ncol number of cols
//' @param threads number of threads (unused currently)
//' @return dgCMatrix sparse matrix
//' @keywords internal
// [[Rcpp::export]]
Rcpp::S4 cpp_triplets_to_dgC(const Rcpp::IntegerVector& i,
                             const Rcpp::IntegerVector& j,
                             const Rcpp::NumericVector& x,
                             int nrow, int ncol,
                             int threads = 1) {
  const R_xlen_t nnz = x.size();
  if (i.size() != nnz || j.size() != nnz) {
    Rcpp::stop("i, j, x must have same length");
  }
  if (nrow < 0 || ncol < 0) {
    Rcpp::stop("nrow and ncol must be non-negative");
  }

  std::vector<std::map<int, double> > cols(static_cast<std::size_t>(ncol));

  for (R_xlen_t k = 0; k < nnz; ++k) {
    const int row = i[k] - 1;
    const int col = j[k] - 1;
    if (row < 0 || row >= nrow || col < 0 || col >= ncol) {
      Rcpp::stop("triplet index out of bounds");
    }
    cols[static_cast<std::size_t>(col)][row] += x[k];
  }

  R_xlen_t nnz_out = 0;
  for (int col = 0; col < ncol; ++col) {
    for (std::map<int, double>::const_iterator it = cols[static_cast<std::size_t>(col)].begin();
         it != cols[static_cast<std::size_t>(col)].end(); ++it) {
      if (it->second != 0.0) {
        ++nnz_out;
      }
    }
  }

  if (nnz_out > static_cast<R_xlen_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("sparse matrix has too many non-zero entries for dgCMatrix");
  }

  Rcpp::IntegerVector p(ncol + 1);
  Rcpp::IntegerVector i_out(static_cast<int>(nnz_out));
  Rcpp::NumericVector x_out(static_cast<int>(nnz_out));

  R_xlen_t offset = 0;
  p[0] = 0;
  for (int col = 0; col < ncol; ++col) {
    for (std::map<int, double>::const_iterator it = cols[static_cast<std::size_t>(col)].begin();
         it != cols[static_cast<std::size_t>(col)].end(); ++it) {
      if (it->second == 0.0) {
        continue;
      }
      i_out[static_cast<int>(offset)] = it->first;
      x_out[static_cast<int>(offset)] = it->second;
      ++offset;
    }
    p[col + 1] = static_cast<int>(offset);
  }

  Rcpp::S4 out("dgCMatrix");
  out.slot("Dim") = Rcpp::IntegerVector::create(nrow, ncol);
  out.slot("p") = p;
  out.slot("i") = i_out;
  out.slot("x") = x_out;
  out.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
  return out;
}
