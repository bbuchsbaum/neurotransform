#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// Convert dgCMatrix S4 to arma::sp_mat
arma::sp_mat as_spmat(const Rcpp::S4& mat) {
  Rcpp::IntegerVector dims = mat.slot("Dim");
  arma::uword n_rows = dims[0];
  arma::uword n_cols = dims[1];
  Rcpp::IntegerVector p = mat.slot("p");    // column pointers
  Rcpp::IntegerVector i = mat.slot("i");    // row indices
  Rcpp::NumericVector x = mat.slot("x");    // values
  arma::uword nnz = x.size();
  arma::umat locations(2, nnz);
  arma::vec values(x.begin(), nnz, false, true);

  arma::uword idx = 0;
  for (arma::uword col = 0; col < n_cols; ++col) {
    for (arma::uword k = p[col]; k < static_cast<arma::uword>(p[col + 1]); ++k) {
      locations(0, idx) = static_cast<arma::uword>(i[k]);   // row
      locations(1, idx) = col;                              // col
      ++idx;
    }
  }
  return arma::sp_mat(locations, values, n_rows, n_cols, true, false);
}

//' Apply sparse projector matrix to data
//'
//' Fast C++ sparse matrix multiplication.
//'
//' @param projMat dgCMatrix sparse projector
//' @param data Numeric matrix of data
//' @param threads Number of OpenMP threads
//' @return Projected data matrix
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_apply_projector(const Rcpp::S4& projMat,
                                        const Rcpp::NumericMatrix& data,
                                        int threads = 1) {
  arma::sp_mat P = as_spmat(projMat);
  arma::mat X = Rcpp::as<arma::mat>(data); // copy

#ifdef _OPENMP
  if (threads > 0) omp_set_num_threads(threads);
#endif

  arma::mat Y = P * X;

  Rcpp::NumericMatrix out(Y.n_rows, Y.n_cols);
  std::copy(Y.begin(), Y.end(), out.begin());
  return out;
}

//' Apply chain of affine transformations to coordinates
//'
//' Efficiently applies multiple 4x4 affine matrices to coordinates.
//'
//' @param coords Numeric matrix (N x 3) of coordinates
//' @param matrices List of 4x4 affine matrices
//' @return Transformed coordinates
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_apply_affine_chain(const Rcpp::NumericMatrix& coords,
                                           const Rcpp::List& matrices) {
  if (matrices.size() == 0) return coords;
  arma::mat C = Rcpp::as<arma::mat>(coords); // N x 3 copy
  arma::mat A = arma::eye(4,4);
  for (int k = 0; k < matrices.size(); ++k) {
    Rcpp::NumericMatrix m = matrices[k];
    if (m.nrow() != 4 || m.ncol() != 4) Rcpp::stop("Affine matrices must be 4x4");
    arma::mat Mk(m.begin(), 4, 4, false, true);
    A = Mk * A;
  }
  arma::mat hom(C.n_rows, 4, arma::fill::ones);
  hom.cols(0,2) = C;
  arma::mat tmp = (A * hom.t()).t();
  arma::mat res = tmp.cols(0,2);
  Rcpp::NumericMatrix out(res.n_rows, res.n_cols);
  std::copy(res.begin(), res.end(), out.begin());
  return out;
}
