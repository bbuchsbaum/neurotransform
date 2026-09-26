#ifndef NEUROTRANSFORM_SURFACE_PROJECTION_H
#define NEUROTRANSFORM_SURFACE_PROJECTION_H
#include <RcppArmadillo.h>
#include <algorithm>
#include <array>
namespace neurotransform {
inline void validate_mesh(const arma::mat& verts, const arma::imat& faces) {
  if (verts.n_cols != 3 || faces.n_cols != 3)
    Rcpp::stop("vertices and faces must have three columns");
  if (!verts.is_finite()) Rcpp::stop("vertices must be finite");
  for (arma::uword j = 0; j < faces.n_elem; ++j) {
    if (faces[j] < 0 || faces[j] >= static_cast<arma::sword>(verts.n_rows))
      Rcpp::stop("faces must contain valid zero-based vertex indices");
  }
}

}
#endif
