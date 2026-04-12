#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

//' Apply sparse projector matrix to data
//'
//' Fast C++ sparse matrix multiplication.
//'
//' @param projMat dgCMatrix sparse projector
//' @param data Numeric matrix of data
//' @param threads Number of OpenMP threads
//' @return Projected data matrix
//' @keywords internal
//
// =============================================================================
// cpp_apply_projector - direct dgCMatrix SpMV
// =============================================================================
//
// This is a deliberate specialization that bypasses arma::sp_mat. Do not
// "simplify" to arma::sp_mat * arma::mat.
//
// RATIONALE
// ---------
// The previous arma::sp_mat implementation had two fatal problems at
// realistic neuroimaging scales:
//
// 1. ARMA 32-BIT OVERFLOW.
//    Without ARMA_64BIT_WORD (the default build), arma::sp_mat stores
//    element counts in 32-bit uwords. arma::SpMat::init() validates that
//    n_rows * n_cols fits in a uword. A surface-to-volume projector of the
//    bread-and-butter shape 300k target voxels x 32k source vertices
//    produces n_rows * n_cols ~= 9.6e9, which overflows and hard-crashes
//    with "SpMat::init(): requested size is too large; suggest to enable
//    ARMA_64BIT_WORD". All four of the following shapes crash the
//    arma::sp_mat path:
//      - 300k x 32k  (fsLR-32k -> MNI-2mm)
//      - 32k x 300k  (MNI-2mm  -> fsLR-32k)
//      - 1M  x 32k   (extreme-tall)
//      - 32k x 1M    (extreme-wide)
//
// 2. COPY OVERHEAD AND SERIAL SpMV.
//    arma::sp_mat * arma::mat allocates a full dense output (n_target x T)
//    and copies input and output once each. arma's SpMV does not
//    parallelize over the time axis T. The direct-slot path below reads
//    dgCMatrix slots in-place, writes the output Rcpp::NumericMatrix
//    directly, and parallelizes the outer loop over T.
//
// CORRECTNESS
// -----------
// Outputs match arma::sp_mat * arma::mat to 1e-10 on all shapes where the
// arma path does not crash. See tools/bench_apply_projector.R for the
// reproducer and measured numbers.
//
// WHEN TO RECONSIDER
// ------------------
// This specialization becomes unnecessary if BOTH of the following hold:
//   (a) ARMA_64BIT_WORD is enabled in the Armadillo build, AND
//   (b) a benchmark reruns on the same shapes and shows arma::sp_mat is at
//       least competitive across single- and multi-threaded runs.
// Until then, do not change this to arma::sp_mat.
// =============================================================================
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_apply_projector(const Rcpp::S4& projMat,
                                        const Rcpp::NumericMatrix& data,
                                        int threads = 1) {
  if (threads < 1) threads = 1;

  Rcpp::IntegerVector dims = projMat.slot("Dim");
  const int n_rows = dims[0];
  const int n_cols = dims[1];

  if (data.nrow() != n_cols) {
    Rcpp::stop("Data has %d rows but projector expects %d source elements",
               data.nrow(), n_cols);
  }

  Rcpp::IntegerVector p = projMat.slot("p");
  Rcpp::IntegerVector i = projMat.slot("i");
  Rcpp::NumericVector x = projMat.slot("x");

  if (p.size() != n_cols + 1) {
    Rcpp::stop("Invalid dgCMatrix: p slot has length %d but expected %d",
               p.size(), n_cols + 1);
  }
  if (i.size() != x.size()) {
    Rcpp::stop("Invalid dgCMatrix: i slot has length %d but x slot has length %d",
               i.size(), x.size());
  }

  const int n_time = data.ncol();
  Rcpp::NumericMatrix out(n_rows, n_time);

  const double* X = data.begin();
  double* Y = out.begin();

#ifdef _OPENMP
  if (threads > 1) {
    omp_set_num_threads(threads);
    #pragma omp parallel for
    for (int c = 0; c < n_time; ++c) {
      for (int col = 0; col < n_cols; ++col) {
        const int k0 = p[col];
        const int k1 = p[col + 1];
        if (k0 == k1) continue;
        const double xcol = X[col + c * n_cols];
        if (xcol == 0.0) continue;
        for (int k = k0; k < k1; ++k) {
          const int row = i[k];
          Y[row + c * n_rows] += x[k] * xcol;
        }
      }
    }
  } else {
    for (int c = 0; c < n_time; ++c) {
      for (int col = 0; col < n_cols; ++col) {
        const int k0 = p[col];
        const int k1 = p[col + 1];
        if (k0 == k1) continue;
        const double xcol = X[col + c * n_cols];
        if (xcol == 0.0) continue;
        for (int k = k0; k < k1; ++k) {
          const int row = i[k];
          Y[row + c * n_rows] += x[k] * xcol;
        }
      }
    }
  }
#else
  for (int c = 0; c < n_time; ++c) {
    for (int col = 0; col < n_cols; ++col) {
      const int k0 = p[col];
      const int k1 = p[col + 1];
      if (k0 == k1) continue;
      const double xcol = X[col + c * n_cols];
      if (xcol == 0.0) continue;
      for (int k = k0; k < k1; ++k) {
        const int row = i[k];
        Y[row + c * n_rows] += x[k] * xcol;
      }
    }
  }
#endif

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
