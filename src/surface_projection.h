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

// Find the nearest orthogonal projection supported by a triangle. In closest
// mode, also consider triangle edges/vertices (Euclidean mesh projection).
inline bool bary_point(const arma::rowvec& pt, const arma::mat& verts,
                       const arma::imat& faces, bool closest,
                       arma::uword& face_idx, arma::rowvec& bary) {
  if (!pt.is_finite()) return false;
  double best = std::numeric_limits<double>::infinity();
  arma::rowvec best_bary(3);
  arma::uword best_face = 0;
  std::array<arma::sword, 3> best_ids = {{0, 0, 0}};
  const double tol = 64 * std::numeric_limits<double>::epsilon();
  for (arma::uword f = 0; f < faces.n_rows; ++f) {
    // Canonical order makes winding and face-row permutations immaterial,
    // including the tie-break for exactly equidistant, disjoint triangles.
    std::array<arma::sword, 3> ids = {{faces(f,0), faces(f,1), faces(f,2)}};
    std::sort(ids.begin(), ids.end());
    arma::rowvec v1 = verts.row(ids[0]);
    arma::rowvec v2 = verts.row(ids[1]);
    arma::rowvec v3 = verts.row(ids[2]);
    arma::rowvec v0 = v2 - v1;
    arma::rowvec v1p = v3 - v1;
    arma::rowvec vp = pt - v1;
    double d00 = arma::dot(v0, v0);
    double d01 = arma::dot(v0, v1p);
    double d11 = arma::dot(v1p, v1p);
    double d20 = arma::dot(vp, v0);
    double d21 = arma::dot(vp, v1p);
    arma::rowvec normal = arma::cross(v0, v1p);
    double denom = arma::dot(normal, normal);
    double scale2 = std::max(d00, std::max(d11, arma::dot(v3 - v2, v3 - v2)));
    // A relative degeneracy check, not a fixed additive perturbation of the
    // denominator: valid weights must not depend on mesh units or radius.
    if (!std::isfinite(denom) || denom <= tol * scale2 * scale2) continue;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;
    arma::rowvec candidate = {u, v, w};
    if (u >= -tol && v >= -tol && w >= -tol) {
      candidate = arma::clamp(candidate, 0.0, 1.0);
      candidate /= arma::accu(candidate);
    } else if (closest) {
      // If the plane projection is outside, the closest point lies on an
      // edge. Project onto each segment; endpoints are included by clamping.
      double edge_best = std::numeric_limits<double>::infinity();
      arma::rowvec corners[3] = {v1, v2, v3};
      for (int a = 0; a < 3; ++a) {
        int b = (a + 1) % 3;
        arma::rowvec edge = corners[b] - corners[a];
        arma::rowvec delta = pt - corners[a];
        double t = arma::dot(delta, edge) / arma::dot(edge, edge);
        t = std::max(0.0, std::min(1.0, t));
        arma::rowvec residual = delta - t * edge;
        double dist = arma::dot(residual, residual);
        if (dist < edge_best) {
          edge_best = dist;
          candidate.zeros();
          candidate[a] = 1 - t;
          candidate[b] = t;
        }
      }
    } else {
      continue;
    }
    arma::rowvec residual = vp - candidate[1] * v0 - candidate[2] * v1p;
    double dist = arma::dot(residual, residual);
    if (!std::isfinite(dist)) continue;
    if (dist < best || (dist == best && ids < best_ids)) {
      best = dist;
      best_ids = ids;
      for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
          if (faces(f,k) == ids[j]) best_bary[k] = candidate[j];
        }
      }
      best_face = f;
    }
  }
  if (best == std::numeric_limits<double>::infinity()) return false;
  face_idx = best_face;
  bary = best_bary;
  return true;
}

}
#endif
