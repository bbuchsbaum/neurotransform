#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Generate world coordinates for full volume grid
//'
//' @param dim integer vector length 3
//' @param affine 4x4 voxel-to-world
//' @return Numeric matrix (N x 3) of world coordinates
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_volume_world_coords(const Rcpp::IntegerVector& dim,
                                            const Rcpp::NumericMatrix& affine) {
  if (dim.size() != 3) Rcpp::stop("dim must be length 3");
  if (affine.nrow() != 4 || affine.ncol() != 4) Rcpp::stop("affine must be 4x4");
  int nx = dim[0], ny = dim[1], nz = dim[2];
  const std::size_t n = static_cast<std::size_t>(nx) * ny * nz;
  Rcpp::NumericMatrix out(n, 3);

  std::size_t idx = 0;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i, ++idx) {
        double x = affine(0,0)*i + affine(0,1)*j + affine(0,2)*k + affine(0,3);
        double y = affine(1,0)*i + affine(1,1)*j + affine(1,2)*k + affine(1,3);
        double z = affine(2,0)*i + affine(2,1)*j + affine(2,2)*k + affine(2,3);
        out(idx, 0) = x;
        out(idx, 1) = y;
        out(idx, 2) = z;
      }
    }
  }
  return out;
}
