#ifndef NEUROTRANSFORM_SURFACE_INDEX_H
#define NEUROTRANSFORM_SURFACE_INDEX_H
#include "surface_projection.h"
#include <numeric>
#include <vector>

namespace neurotransform {
using Point3 = std::array<double,3>;
inline Point3 subtract(const Point3& a, const Point3& b) {
  return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};
}
inline double dot3(const Point3& a, const Point3& b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
inline Point3 cross3(const Point3& a, const Point3& b) {
  return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
}
inline double convex_distance2(const Point3& point, const std::array<Point3,3>& vertices,
                                const Point3& weights) {
  Point3 delta;
  for (int k=0;k<3;++k) {
    double coordinate=0;
    for (int a=0;a<3;++a) coordinate+=weights[a]*vertices[a][k];
    double lower=std::min(vertices[0][k],std::min(vertices[1][k],vertices[2][k]));
    double upper=std::max(vertices[0][k],std::max(vertices[1][k],vertices[2][k]));
    // Reconstruct from original corners (not rounded edge differences). The
    // clamp also makes the projection distance consistent with the AABB bound
    // in floating-point arithmetic, including highly translated triangles.
    coordinate=std::max(lower,std::min(upper,coordinate));
    delta[k]=point[k]-coordinate;
  }
  return dot3(delta,delta);
}
struct SurfaceProjection {
  int face = -1;
  std::array<int,3> ids = {{0,0,0}};
  Point3 weights = {{0,0,0}};
  double distance2 = std::numeric_limits<double>::infinity();
  int tested = 0;
};
struct SurfaceTriangle {
  std::array<int,3> ids;
  std::array<Point3,3> vertices;
  Point3 e1,e2,lower,upper,center;
  double d00,d01,d11,denom;
  int face;
};
struct SurfaceNode {
  Point3 lower,upper;
  int begin,end,left=-1,right=-1;
};

class SurfaceIndex {
  std::vector<SurfaceTriangle> triangles;
  std::vector<int> order;
  std::vector<SurfaceNode> nodes;
  int vertex_count;

  int build(int begin,int end) {
    SurfaceNode node;
    node.begin=begin; node.end=end;
    node.lower.fill(std::numeric_limits<double>::infinity());
    node.upper.fill(-std::numeric_limits<double>::infinity());
    for (int p=begin;p<end;++p) for (int k=0;k<3;++k) {
      const auto& t=triangles[order[p]];
      node.lower[k]=std::min(node.lower[k],t.lower[k]);
      node.upper[k]=std::max(node.upper[k],t.upper[k]);
    }
    int id=nodes.size(); nodes.push_back(node);
    if (end-begin<=8) return id;
    int axis=0;
    for (int k=1;k<3;++k) if (node.upper[k]-node.lower[k]>node.upper[axis]-node.lower[axis]) axis=k;
    int mid=begin+(end-begin)/2;
    std::nth_element(order.begin()+begin,order.begin()+mid,order.begin()+end,[&](int a,int b) {
      double x=triangles[a].center[axis],y=triangles[b].center[axis];
      return x<y || (x==y && triangles[a].ids<triangles[b].ids);
    });
    int left=build(begin,mid),right=build(mid,end);
    nodes[id].left=left; nodes[id].right=right;
    return id;
  }

  double bound(const Point3& p,int index) const {
    const auto& node=nodes[index];
    double result=0;
    for (int k=0;k<3;++k) {
      double d=std::max(0.0,std::max(node.lower[k]-p[k],p[k]-node.upper[k]));
      result+=d*d;
    }
    // Outward-rounded boxes plus a downward margin on squared-distance
    // arithmetic ensure near-tie nodes are searched rather than falsely pruned.
    return std::nextafter(result*(1-64*std::numeric_limits<double>::epsilon()),0.0);
  }

  void project(const Point3& pt,const SurfaceTriangle& t,bool closest,SurfaceProjection& best) const {
    ++best.tested;
    Point3 delta=subtract(pt,t.vertices[0]);
    double d20=dot3(delta,t.e1),d21=dot3(delta,t.e2);
    double v=(t.d11*d20-t.d01*d21)/t.denom;
    double w=(t.d00*d21-t.d01*d20)/t.denom;
    Point3 bary={{1-v-w,v,w}};
    // An exactly coincident query has an exact one-vertex representation. This
    // preserves structural zeros (important for categorical and adaptive support)
    // instead of introducing roundoff contributors through the Gram solve.
    for (int a=0;a<3;++a) if (pt==t.vertices[a]) {
      bary={{0,0,0}};bary[a]=1;break;
    }
    const double tol=64*std::numeric_limits<double>::epsilon();
    if (bary[0]>=-tol && bary[1]>=-tol && bary[2]>=-tol) {
      double sum=0;
      for (double& x:bary) { x=std::max(0.0,std::min(1.0,x));sum+=x; }
      for (double& x:bary) x/=sum;
    } else if (closest) {
      double edge_best=std::numeric_limits<double>::infinity();
      for (int a=0;a<3;++a) {
        int b=(a+1)%3;
        Point3 edge=subtract(t.vertices[b],t.vertices[a]),d=subtract(pt,t.vertices[a]);
        double u=std::max(0.0,std::min(1.0,dot3(d,edge)/dot3(edge,edge)));
        Point3 candidate={{0,0,0}};candidate[a]=1-u;candidate[b]=u;
        double distance=convex_distance2(pt,t.vertices,candidate);
        if (distance<edge_best) { edge_best=distance;bary=candidate; }
      }
    } else return;
    double distance=convex_distance2(pt,t.vertices,bary);
    if (!std::isfinite(distance)) return;
    if (distance<best.distance2 || (distance==best.distance2 &&
        (t.ids<best.ids || (t.ids==best.ids && (best.face<0 || t.face<best.face))))) {
      best.face=t.face;best.ids=t.ids;best.weights=bary;best.distance2=distance;
    }
  }

public:
  SurfaceIndex(const arma::mat& vertices,const arma::imat& faces):vertex_count(vertices.n_rows) {
    validate_mesh(vertices,faces);
    triangles.reserve(faces.n_rows);
    for (arma::uword f=0;f<faces.n_rows;++f) {
      if (f%4096==0) Rcpp::checkUserInterrupt();
      SurfaceTriangle t;
      t.face=f;
      for (int k=0;k<3;++k) t.ids[k]=faces(f,k);
      std::sort(t.ids.begin(),t.ids.end());
      for (int a=0;a<3;++a) for (int k=0;k<3;++k) t.vertices[a][k]=vertices(t.ids[a],k);
      t.e1=subtract(t.vertices[1],t.vertices[0]);t.e2=subtract(t.vertices[2],t.vertices[0]);
      t.d00=dot3(t.e1,t.e1);t.d01=dot3(t.e1,t.e2);t.d11=dot3(t.e2,t.e2);
      auto normal=cross3(t.e1,t.e2),third=subtract(t.vertices[2],t.vertices[1]);
      t.denom=dot3(normal,normal);
      double scale2=std::max(t.d00,std::max(t.d11,dot3(third,third)));
      if (!std::isfinite(t.denom) || t.denom<=64*std::numeric_limits<double>::epsilon()*scale2*scale2) continue;
      for (int k=0;k<3;++k) {
        t.lower[k]=std::nextafter(std::min(t.vertices[0][k],std::min(t.vertices[1][k],t.vertices[2][k])), -std::numeric_limits<double>::infinity());
        t.upper[k]=std::nextafter(std::max(t.vertices[0][k],std::max(t.vertices[1][k],t.vertices[2][k])), std::numeric_limits<double>::infinity());
        t.center[k]=t.vertices[0][k]/3+t.vertices[1][k]/3+t.vertices[2][k]/3;
      }
      triangles.push_back(t);
    }
    order.resize(triangles.size());std::iota(order.begin(),order.end(),0);
    nodes.reserve(2*triangles.size());
    if (!triangles.empty()) build(0,triangles.size());
  }
  int nvertices() const { return vertex_count; }
  SurfaceProjection query(const Point3& point,bool closest,bool indexed=true) const {
    SurfaceProjection result;
    if (nodes.empty() || !std::isfinite(point[0]) || !std::isfinite(point[1]) || !std::isfinite(point[2])) return result;
    if (!indexed) { for (const auto& t:triangles) project(point,t,closest,result);return result; }
    std::vector<std::pair<int,double>> stack;stack.reserve(64);stack.emplace_back(0,bound(point,0));
    while (!stack.empty()) {
      auto item=stack.back();stack.pop_back();
      if (item.second>result.distance2) continue;
      const auto& node=nodes[item.first];
      if (node.left<0) {
        for (int p=node.begin;p<node.end;++p) project(point,triangles[order[p]],closest,result);
      } else {
        double left=bound(point,node.left),right=bound(point,node.right);
        if (left<=right) {
          if (right<=result.distance2) stack.emplace_back(node.right,right);
          if (left<=result.distance2) stack.emplace_back(node.left,left);
        } else {
          if (left<=result.distance2) stack.emplace_back(node.left,left);
          if (right<=result.distance2) stack.emplace_back(node.right,right);
        }
      }
    }
    return result;
  }
};

inline Rcpp::NumericVector sample_surface(const SurfaceIndex& mesh,
                                         const Rcpp::NumericMatrix& coords,
                                         const Rcpp::NumericVector& data,bool closest) {
  if (coords.ncol()!=3) Rcpp::stop("coords must have three columns");
  int n=coords.nrow(),nv=mesh.nvertices(),k=1;
  bool matrix_data=data.hasAttribute("dim");
  if (matrix_data) {
    Rcpp::IntegerVector dims=data.attr("dim");
    if (dims.size()!=2 || dims[0]!=nv) Rcpp::stop("data must have one row per vertex");
    k=dims[1];
  }
  if (data.size()!=static_cast<R_xlen_t>(nv)*k) Rcpp::stop("data must have one value or row per vertex");
  Rcpp::NumericVector out(static_cast<R_xlen_t>(n)*k,NA_REAL);
  if (matrix_data) out.attr("dim")=Rcpp::IntegerVector::create(n,k);
  for (int begin=0;begin<n;begin+=4096) {
    Rcpp::checkUserInterrupt();
    int end=std::min(n,begin+4096);
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i=begin;i<end;++i) {
      Point3 point={{coords(i,0),coords(i,1),coords(i,2)}};
      auto p=mesh.query(point,closest);
      if (p.face<0) continue;
      for (int col=0;col<k;++col) {
        double value=0;
        for (int j=0;j<3;++j) if (p.weights[j]>0) value+=p.weights[j]*data[p.ids[j]+nv*col];
        out[i+n*col]=value;
      }
    }
  }
  return out;
}
}
#endif
