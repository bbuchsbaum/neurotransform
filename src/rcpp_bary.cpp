#include <RcppArmadillo.h>
#include "surface_projection.h"
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute barycentric interpolation weights for surface mesh
//'
//' Finds the closest point on the mesh, including edge and vertex projections.
//' With `closest=FALSE`, only interior orthogonal projections are considered.
//' Degenerate faces are ignored; unsupported queries have no triplets.
//' This is Euclidean projection, not radial ray intersection.
//'
//' @param coords Numeric matrix (N x 3) of query points
//' @param vertices Numeric matrix (V x 3) of mesh vertices
//' @param faces Integer matrix (F x 3) of face indices (0-based)
//' @param closest Include closest edge/vertex points when the plane projection
//'   is outside a triangle (used for spherical resampling).
//' @return List with rows, cols, vals for sparse matrix construction
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_barycentric_weights(const Rcpp::NumericMatrix& coords,
                                   const Rcpp::NumericMatrix& vertices,
                                   const Rcpp::IntegerMatrix& faces,
                                   bool closest = true) {
  if (coords.ncol() != 3 || vertices.ncol() != 3 || faces.ncol() != 3)
    Rcpp::stop("coords, vertices and faces must have three columns");
  arma::mat pts = Rcpp::as<arma::mat>(coords);
  arma::mat verts = Rcpp::as<arma::mat>(vertices);
  // copy into Armadillo integer matrix to avoid strides/ownership surprises
  arma::imat f = Rcpp::as<arma::imat>(faces);
  neurotransform::validate_mesh(verts, f);

  std::vector<int> rows;
  std::vector<int> cols;
  std::vector<double> vals;

#ifdef _OPENMP
  #pragma omp parallel
  {
    std::vector<int> lrows; lrows.reserve(1024);
    std::vector<int> lcols; lcols.reserve(1024);
    std::vector<double> lvals; lvals.reserve(1024);

    #pragma omp for nowait
    for (int i = 0; i < pts.n_rows; ++i) {
      arma::uword face_idx;
      arma::rowvec bary(3);
      if (!neurotransform::bary_point(pts.row(i), verts, f, closest, face_idx, bary)) continue;
      arma::irowvec vids = f.row(face_idx);
      for (int k = 0; k < 3; ++k) {
        if (bary[k] > 0) {
          lrows.push_back(i+1);
          lcols.push_back(vids[k] + 1); // R 1-based
          lvals.push_back(bary[k]);
        }
      }
    }
    #pragma omp critical
    {
      rows.insert(rows.end(), lrows.begin(), lrows.end());
      cols.insert(cols.end(), lcols.begin(), lcols.end());
      vals.insert(vals.end(), lvals.begin(), lvals.end());
    }
  }
#else
  for (int i = 0; i < pts.n_rows; ++i) {
    arma::uword face_idx;
    arma::rowvec bary(3);
    if (!neurotransform::bary_point(pts.row(i), verts, f, closest, face_idx, bary)) continue;
    arma::irowvec vids = f.row(face_idx);
    for (int k = 0; k < 3; ++k) {
      if (bary[k] > 0) {
        rows.push_back(i+1);
        cols.push_back(vids[k] + 1);
        vals.push_back(bary[k]);
      }
    }
  }
#endif

  return Rcpp::List::create(
    Rcpp::Named("rows") = rows,
    Rcpp::Named("cols") = cols,
    Rcpp::Named("vals") = vals
  );
}
