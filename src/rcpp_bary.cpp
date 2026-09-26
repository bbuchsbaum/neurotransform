#include <RcppArmadillo.h>
#include "surface_index.h"
#include <chrono>
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
//' @param indexed Use the spatial index; FALSE retains exhaustive search for validation.
//' @return List with sparse triplets, query diagnostics and build/query timings
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_barycentric_weights(const Rcpp::NumericMatrix& coords,
                                   const Rcpp::NumericMatrix& vertices,
                                   const Rcpp::IntegerMatrix& faces,
                                   bool closest = true, bool indexed = true) {
  if (coords.ncol()!=3) Rcpp::stop("coords must have three columns");
  auto start=std::chrono::steady_clock::now();
  neurotransform::SurfaceIndex mesh(Rcpp::as<arma::mat>(vertices), Rcpp::as<arma::imat>(faces));
  auto built=std::chrono::steady_clock::now();
  int n=coords.nrow();
  std::vector<neurotransform::SurfaceProjection> projections(n);
  // Interrupt checks stay on the R thread, outside OpenMP workers.
  for (int begin=0;begin<n;begin+=4096) {
    Rcpp::checkUserInterrupt();
    int end=std::min(n,begin+4096);
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i=begin;i<end;++i) {
      neurotransform::Point3 p={{coords(i,0),coords(i,1),coords(i,2)}};
      projections[i]=mesh.query(p,closest,indexed);
    }
  }
  auto queried=std::chrono::steady_clock::now();
  std::vector<int> rows,cols;
  std::vector<double> vals;
  rows.reserve(3*n);cols.reserve(3*n);vals.reserve(3*n);
  Rcpp::NumericVector distances(n,NA_REAL);
  Rcpp::IntegerVector face(n,NA_INTEGER),tested(n);
  for (int i=0;i<n;++i) {
    const auto& p=projections[i];tested[i]=p.tested;
    if (p.face<0) continue;
    face[i]=p.face+1;distances[i]=std::sqrt(p.distance2);
    for (int k=0;k<3;++k) if (p.weights[k]>0) {
      rows.push_back(i+1);cols.push_back(p.ids[k]+1);vals.push_back(p.weights[k]);
    }
  }
  return Rcpp::List::create(Rcpp::Named("rows")=rows,Rcpp::Named("cols")=cols,Rcpp::Named("vals")=vals,
    Rcpp::Named("distance")=distances,Rcpp::Named("face")=face,Rcpp::Named("triangle_tests")=tested,
    Rcpp::Named("timing")=Rcpp::List::create(
      Rcpp::Named("index_seconds")=std::chrono::duration<double>(built-start).count(),
      Rcpp::Named("query_seconds")=std::chrono::duration<double>(queried-built).count()));
}
