#include <RcppArmadillo.h>
#include <array>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include "surface_projection.h"

namespace {
using Point = std::array<double, 3>;
// Filtered determinant: ambiguous signs are rejected, never guessed. The
// conservative bound dominates the rounding error of the six triple products.
std::pair<double, double> determinant(const Point& a, const Point& b, const Point& c) {
  double terms[6] = {a[0]*b[1]*c[2], a[1]*b[2]*c[0], a[2]*b[0]*c[1],
                    a[2]*b[1]*c[0], a[1]*b[0]*c[2], a[0]*b[2]*c[1]};
  double value = terms[0]+terms[1]+terms[2]-terms[3]-terms[4]-terms[5];
  double permanent = 0;
  for (double x : terms) permanent += std::abs(x);
  return {value, 64 * std::numeric_limits<double>::epsilon() * permanent +
                 64 * std::numeric_limits<double>::min()};
}
double dot(const Point& a, const Point& b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
struct Edge {
  int first, slot, direction, count = 1;
  Edge(int f, int k, int d): first(f), slot(k), direction(d) {}
};
}

// [[Rcpp::export]]
Rcpp::List cpp_validate_surface(const Rcpp::NumericMatrix& vertices,
                                const Rcpp::IntegerMatrix& faces,
                                bool spherical = true, double radius_tolerance = 1.001) {
  if (!std::isfinite(radius_tolerance) || radius_tolerance <= 1)
    Rcpp::stop("radius_tolerance must be finite and greater than one");
  arma::mat v = Rcpp::as<arma::mat>(vertices);
  arma::imat f = Rcpp::as<arma::imat>(faces);
  neurotransform::validate_mesh(v, f);
  if (v.n_rows == 0) Rcpp::stop("vertices must be non-empty");
  const int nv = v.n_rows, nf = f.n_rows;
  std::map<std::string, std::set<int>> issues;
  std::vector<Point> unit(nv);
  double min_radius = R_PosInf, max_radius = 0;
  for (int i = 0; i < nv; ++i) {
    double r = std::hypot(v(i,0), std::hypot(v(i,1), v(i,2)));
    min_radius = std::min(min_radius, r); max_radius = std::max(max_radius, r);
    if (!(r > 0) || !std::isfinite(r)) {
      if (spherical) issues["invalid_radius_vertices"].insert(i+1);
    } else unit[i] = {{v(i,0)/r, v(i,1)/r, v(i,2)/r}};
  }
  if (spherical && !(min_radius * radius_tolerance > max_radius))
    issues["radius_spread"].insert(1);
  // Test shape after safe scaling. Spherical projection uses unit directions,
  // so admission checks the geometry that will actually be rescaled.
  arma::mat quality = v;
  if (spherical) {
    for (int i=0; i<nv; ++i) for (int k=0; k<3; ++k) quality(i,k)=unit[i][k];
  } else {
    double scale=arma::abs(v).max();
    if (scale>0) quality/=scale;
  }
  if (nf == 0) issues["no_faces"].insert(1);
  std::vector<std::vector<int>> incident(nv);
  std::vector<std::array<int,3>> adj(nf, {{-1,-1,-1}}), relation(nf, {{0,0,0}});
  std::unordered_map<uint64_t, Edge> edges;
  edges.reserve(3 * static_cast<size_t>(nf));
  std::map<std::array<int,3>, int> unique_faces;
  const double epsilon = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < nf; ++i) {
    if (i % 4096 == 0) Rcpp::checkUserInterrupt();
    std::array<int,3> ids = {{static_cast<int>(f(i,0)), static_cast<int>(f(i,1)), static_cast<int>(f(i,2))}};
    auto sorted = ids; std::sort(sorted.begin(), sorted.end());
    if (sorted[0] == sorted[1] || sorted[1] == sorted[2]) {
      issues["repeated_vertex_faces"].insert(i+1); continue;
    }
    auto duplicate = unique_faces.emplace(sorted, i);
    if (!duplicate.second) {
      issues["duplicate_faces"].insert(i+1);
      issues["duplicate_faces"].insert(duplicate.first->second+1);
    }
    arma::rowvec e1 = quality.row(ids[1])-quality.row(ids[0]), e2 = quality.row(ids[2])-quality.row(ids[0]);
    arma::rowvec normal = arma::cross(e1,e2), e3 = e2-e1;
    double scale2 = std::max(arma::dot(e1,e1),std::max(arma::dot(e2,e2),arma::dot(e3,e3)));
    if (!(arma::dot(normal,normal) > 64*epsilon*scale2*scale2))
      issues["degenerate_faces"].insert(i+1);
    for (int k = 0; k < 3; ++k) {
      incident[ids[k]].push_back(i);
      int a = ids[k], b = ids[(k+1)%3], direction = a < b ? 1 : -1;
      uint64_t key = (static_cast<uint64_t>(std::min(a,b)) << 32) | static_cast<uint32_t>(std::max(a,b));
      auto entry = edges.emplace(key, Edge(i,k,direction));
      if (!entry.second) {
        Edge& edge = entry.first->second; ++edge.count;
        if (edge.count == 2) {
          adj[i][k] = edge.first; adj[edge.first][edge.slot] = i;
          relation[i][k] = relation[edge.first][edge.slot] = -direction*edge.direction;
        } else {
          issues["nonmanifold_edge_faces"].insert(i+1);
          issues["nonmanifold_edge_faces"].insert(edge.first+1);
        }
      }
    }
  }
  int boundary_edges = 0;
  for (const auto& entry : edges) {
    if (entry.second.count == 1) {
      ++boundary_edges;
      if (spherical) issues["boundary_edge_faces"].insert(entry.second.first+1);
    }
  }
  for (int i = 0; i < nv; ++i) {
    if (i % 4096 == 0) Rcpp::checkUserInterrupt();
    if (incident[i].empty()) {
      if (spherical) issues["unused_vertices"].insert(i+1);
      continue;
    }
    // With closed two-face edges, a connected incident-face link is a cycle.
    if (!spherical || boundary_edges) continue;
    std::unordered_set<int> reached;
    std::vector<int> pending(1, incident[i][0]);
    while (!pending.empty()) {
      int face = pending.back(); pending.pop_back();
      if (!reached.insert(face).second) continue;
      for (int k = 0; k < 3; ++k)
        if ((f(face,k)==i || f(face,(k+1)%3)==i) && adj[face][k]>=0)
          pending.push_back(adj[face][k]);
    }
    if (reached.size()!=incident[i].size()) issues["nonmanifold_vertices"].insert(i+1);
  }
  std::vector<int> orientation(nf,0);
  int components = 0;
  for (int start = 0; start < nf; ++start) {
    if (orientation[start]) continue;
    ++components; orientation[start] = 1;
    std::vector<int> pending(1,start);
    while (!pending.empty()) {
      int face = pending.back(); pending.pop_back();
      for (int k = 0; k < 3; ++k) {
        int other = adj[face][k];
        if (other < 0) continue;
        int expected = orientation[face]*relation[face][k];
        if (!orientation[other]) { orientation[other] = expected; pending.push_back(other); }
        else if (orientation[other]!=expected) issues["nonorientable_faces"].insert(other+1);
      }
    }
  }
  int euler = nv - static_cast<int>(edges.size()) + nf;
  if (spherical && components!=1) issues["disconnected_faces"].insert(components);
  if (spherical && euler!=2) issues["not_sphere_topology"].insert(euler);
  double area = 0, compensation = 0, min_margin = R_PosInf;
  int global_sign = 0, degree = -1, ray_attempt = 0;
  const double area_tolerance = std::max(1e-7, 128*epsilon*nf);
  if (spherical && issues.empty()) {
    for (int i = 0; i < nf; ++i) {
      const Point &a = unit[f(i,0)], &b = unit[f(i,1)], &c = unit[f(i,2)];
      auto d = determinant(a,b,c);
      if (std::abs(d.first) <= d.second) { issues["indeterminate_radial_faces"].insert(i+1); continue; }
      min_margin = std::min(min_margin, std::abs(d.first)/std::max(d.second,std::numeric_limits<double>::min()));
      int sign = d.first*orientation[i] > 0 ? 1 : -1;
      if (!global_sign) global_sign = sign;
      else if (sign != global_sign) issues["folded_faces"].insert(i+1);
      double value = 2*std::atan2(std::abs(d.first),1+dot(a,b)+dot(b,c)+dot(c,a));
      double y = value-compensation, next = area+y;
      compensation = (next-area)-y; area=next;
    }
    // Positively oriented spherical triangles form a branched cover. Count
    // interiors at a generic ray to determine its integer degree without
    // treating a transcendental area tolerance as a topological certificate.
    if (issues.empty()) for (int attempt=1; attempt<=64; ++attempt) {
      Rcpp::checkUserInterrupt();
      Point ray = {{std::sin(attempt*1.234567+.13), std::cos(attempt*2.345678+.27), std::sin(attempt*3.456789+.39)}};
      int hits=0; bool ambiguous=false;
      for (int i=0; i<nf; ++i) {
        bool outside=false, boundary=false;
        for (int k=0; k<3; ++k) {
          auto d=determinant(unit[f(i,k)], unit[f(i,(k+1)%3)],ray);
          double oriented=d.first*orientation[i]*global_sign;
          if (oriented < -d.second) { outside=true; break; }
          if (oriented <= d.second) boundary=true;
        }
        if (!outside) {
          if (boundary) { ambiguous=true; break; }
          ++hits;
        }
      }
      if (!ambiguous) { degree=hits; ray_attempt=attempt; break; }
    }
    if (issues.empty() && degree!=1) issues[degree<0 ? "indeterminate_degree" : "radial_degree_not_one"].insert(degree);
    if (std::abs(area-4*std::acos(-1.0)) > area_tolerance) issues["spherical_area_mismatch"].insert(1);
  }
  Rcpp::List diagnostics;
  for (const auto& entry:issues) diagnostics.push_back(Rcpp::wrap(std::vector<int>(entry.second.begin(),entry.second.end())),entry.first);
  return Rcpp::List::create(Rcpp::Named("valid")=issues.empty(), Rcpp::Named("issues")=diagnostics,
    Rcpp::Named("vertices")=nv, Rcpp::Named("faces")=nf, Rcpp::Named("edges")=edges.size(),
    Rcpp::Named("components")=components, Rcpp::Named("boundary_edges")=boundary_edges,
    Rcpp::Named("euler_characteristic")=euler, Rcpp::Named("radius_range")=Rcpp::NumericVector::create(min_radius,max_radius),
    Rcpp::Named("radial_degree")=degree, Rcpp::Named("degree_ray_attempt")=ray_attempt,
    Rcpp::Named("spherical_area")=area, Rcpp::Named("area_tolerance")=area_tolerance,
    Rcpp::Named("minimum_orientation_margin")=min_margin);
}
