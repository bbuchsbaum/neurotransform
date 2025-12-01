#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
  const std::size_t nnz = x.size();
  if (i.size() != nnz || j.size() != nnz) {
    Rcpp::stop("i, j, x must have same length");
  }
  if (nnz == 0) {
    Rcpp::S4 empty("dgCMatrix");
    empty.slot("Dim") = Rcpp::IntegerVector::create(nrow, ncol);
    empty.slot("p") = Rcpp::IntegerVector(ncol + 1, 0);
    empty.slot("i") = Rcpp::IntegerVector(0);
    empty.slot("x") = Rcpp::NumericVector(0);
    empty.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
    return empty;
  }

  arma::umat locations(2, nnz);
  arma::vec values(nnz);

#ifdef _OPENMP
  #pragma omp parallel for num_threads(threads)
#endif
  for (std::size_t k = 0; k < nnz; ++k) {
    int row = i[k] - 1;
    int col = j[k] - 1;
    if (row < 0 || row >= nrow || col < 0 || col >= ncol) {
      // Rcpp::stop is not thread-safe; mark sentinel and stop later
      locations(0, k) = std::numeric_limits<arma::uword>::max();
      continue;
    }
    locations(0, k) = static_cast<arma::uword>(row);
    locations(1, k) = static_cast<arma::uword>(col);
    values(k) = x[k];
  }

  if (locations.has_nan() ||
      locations.min() == std::numeric_limits<arma::uword>::max()) {
    Rcpp::stop("triplet index out of bounds");
  }

  arma::sp_mat mat(locations, values, nrow, ncol, /*sort*/ true, /*check_zeros*/ true);
  return Rcpp::wrap(mat); // RcppArmadillo returns dgCMatrix
}
