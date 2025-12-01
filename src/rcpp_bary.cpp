#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

// helper to find containing triangle and barycentric coords via brute-force nearest
static bool bary_point(const arma::rowvec& pt, const arma::mat& verts, const arma::imat& faces, arma::uword& face_idx, arma::rowvec& bary) {
  double best = std::numeric_limits<double>::infinity();
  arma::rowvec best_bary(3);
  arma::uword best_face = 0;
  for (arma::uword f = 0; f < faces.n_rows; ++f) {
    arma::rowvec v1 = verts.row(faces(f,0));
    arma::rowvec v2 = verts.row(faces(f,1));
    arma::rowvec v3 = verts.row(faces(f,2));
    arma::rowvec v0 = v2 - v1;
    arma::rowvec v1p = v3 - v1;
    arma::rowvec vp = pt - v1;
    double d00 = arma::dot(v0, v0);
    double d01 = arma::dot(v0, v1p);
    double d11 = arma::dot(v1p, v1p);
    double d20 = arma::dot(vp, v0);
    double d21 = arma::dot(vp, v1p);
    double denom = d00 * d11 - d01 * d01 + 1e-15;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;
    if (u >= -1e-4 && v >= -1e-4 && w >= -1e-4) {
      // inside or near edge
      double dist = std::abs(u) + std::abs(v) + std::abs(w);
      if (dist < best) {
        best = dist;
        best_bary = {u, v, w};
        best_face = f;
      }
    }
  }
  if (best == std::numeric_limits<double>::infinity()) return false;
  face_idx = best_face;
  // clamp and renormalize
  best_bary = arma::clamp(best_bary, 0.0, 1.0);
  best_bary = best_bary / arma::accu(best_bary);
  bary = best_bary;
  return true;
}

//' Compute barycentric interpolation weights for surface mesh
//'
//' Projects points onto a triangular mesh and computes barycentric weights.
//'
//' @param coords Numeric matrix (N x 3) of query points
//' @param vertices Numeric matrix (V x 3) of mesh vertices
//' @param faces Integer matrix (F x 3) of face indices (0-based)
//' @return List with rows, cols, vals for sparse matrix construction
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_barycentric_weights(const Rcpp::NumericMatrix& coords,
                                   const Rcpp::NumericMatrix& vertices,
                                   const Rcpp::IntegerMatrix& faces) {
  arma::mat pts = Rcpp::as<arma::mat>(coords);
  arma::mat verts = Rcpp::as<arma::mat>(vertices);
  // copy into Armadillo integer matrix to avoid strides/ownership surprises
  arma::imat f = Rcpp::as<arma::imat>(faces);

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
      if (!bary_point(pts.row(i), verts, f, face_idx, bary)) continue;
      arma::irowvec vids = f.row(face_idx);
      for (int k = 0; k < 3; ++k) {
        if (bary[k] > 1e-10) {
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
    if (!bary_point(pts.row(i), verts, f, face_idx, bary)) continue;
    arma::irowvec vids = f.row(face_idx);
    for (int k = 0; k < 3; ++k) {
      if (bary[k] > 1e-10) {
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
