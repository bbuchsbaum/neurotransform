#include <RcppArmadillo.h>
#include <limits>
#include <unordered_map>
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]

//' Compute ribbon sampling weights for volume-to-surface projection
//'
//' Samples along cortical ribbon (inner to outer surface) and accumulates
//' trilinear interpolation weights.
//'
//' @param inner Inner surface coordinates (N x 3)
//' @param outer Outer surface coordinates (N x 3)
//' @param world_to_vox 4x4 world-to-voxel transform
//' @param dims Volume dimensions (3)
//' @param n_steps Number of samples along ribbon
//' @param mask Optional logical mask
//' @return List with rows, cols, vals for sparse matrix construction
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ribbon_weights(const Rcpp::NumericMatrix& inner,
                              const Rcpp::NumericMatrix& outer,
                              const Rcpp::NumericMatrix& world_to_vox,
                              const Rcpp::IntegerVector& dims,
                              int n_steps,
                              const Rcpp::LogicalVector& mask) {
  const int n = inner.nrow();
  const int nx = dims[0], ny = dims[1], nz = dims[2];
  const R_xlen_t nvox = static_cast<R_xlen_t>(nx) * ny * nz;
  if (nvox > std::numeric_limits<int>::max()) {
    Rcpp::stop("volume too large for cpp_ribbon_weights index representation");
  }
  const bool use_mask = mask.size() == nx * ny * nz;

  std::vector<int> rows; rows.reserve(n * 40); // heuristic
  std::vector<int> cols; cols.reserve(n * 40);
  std::vector<double> vals; vals.reserve(n * 40);

#ifdef _OPENMP
  #pragma omp parallel
  {
    std::vector<int> lrows; lrows.reserve(1024);
    std::vector<int> lcols; lcols.reserve(1024);
    std::vector<double> lvals; lvals.reserve(1024);

    #pragma omp for nowait
    for (int i = 0; i < n; ++i) {
      std::unordered_map<int,double> accum;
      accum.reserve(32);
      double x0 = inner(i, 0), y0 = inner(i, 1), z0 = inner(i, 2);
      double x1 = outer(i, 0), y1 = outer(i, 1), z1 = outer(i, 2);

      for (int s = 0; s <= n_steps; ++s) {
        double alpha = static_cast<double>(s) / static_cast<double>(n_steps);
        double wx = x0 + (x1 - x0) * alpha;
        double wy = y0 + (y1 - y0) * alpha;
        double wz = z0 + (z1 - z0) * alpha;

        // homogeneous multiply world_to_vox (row-major from R)
        double vx = world_to_vox(0,0)*wx + world_to_vox(0,1)*wy + world_to_vox(0,2)*wz + world_to_vox(0,3);
        double vy = world_to_vox(1,0)*wx + world_to_vox(1,1)*wy + world_to_vox(1,2)*wz + world_to_vox(1,3);
        double vz = world_to_vox(2,0)*wx + world_to_vox(2,1)*wy + world_to_vox(2,2)*wz + world_to_vox(2,3);

        if (vx < -0.5 || vy < -0.5 || vz < -0.5 || vx > nx - 0.5 || vy > ny - 0.5 || vz > nz - 0.5) {
          continue;
        }

        int xfloor = static_cast<int>(std::floor(vx));
        int yfloor = static_cast<int>(std::floor(vy));
        int zfloor = static_cast<int>(std::floor(vz));
        double fx = vx - xfloor, fy = vy - yfloor, fz = vz - zfloor;
        double wxs[2] = {1.0 - fx, fx};
        double wys[2] = {1.0 - fy, fy};
        double wzs[2] = {1.0 - fz, fz};

        for (int iz = 0; iz < 2; ++iz) {
          for (int iy = 0; iy < 2; ++iy) {
            for (int ix = 0; ix < 2; ++ix) {
              int xi = std::min(std::max(xfloor + ix, 0), nx - 1);
              int yi = std::min(std::max(yfloor + iy, 0), ny - 1);
              int zi = std::min(std::max(zfloor + iz, 0), nz - 1);
              const int lin = static_cast<int>(
                static_cast<R_xlen_t>(xi) +
                static_cast<R_xlen_t>(yi) * nx +
                static_cast<R_xlen_t>(zi) * nx * ny
              ); // 0-based
              if (use_mask && (mask[lin] == FALSE || mask[lin] == NA_LOGICAL)) {
                continue;
              }
              double w = wxs[ix] * wys[iy] * wzs[iz];
              accum[lin] += w;
            }
          }
        }
      } // s

      double total = 0.0;
      for (const auto& kv : accum) total += kv.second;
      if (total <= 0.0) continue;
      for (const auto& kv : accum) {
        lrows.push_back(i + 1);          // 1-based row
        lcols.push_back(kv.first + 1);   // 1-based col
        lvals.push_back(kv.second / total);
      }
    } // i

    #pragma omp critical
    {
      rows.insert(rows.end(), lrows.begin(), lrows.end());
      cols.insert(cols.end(), lcols.begin(), lcols.end());
      vals.insert(vals.end(), lvals.begin(), lvals.end());
    }
  }
#else
  for (int i = 0; i < n; ++i) {
    std::unordered_map<int,double> accum;
    accum.reserve(32);
    double x0 = inner(i, 0), y0 = inner(i, 1), z0 = inner(i, 2);
    double x1 = outer(i, 0), y1 = outer(i, 1), z1 = outer(i, 2);

    for (int s = 0; s <= n_steps; ++s) {
      double alpha = static_cast<double>(s) / static_cast<double>(n_steps);
      double wx = x0 + (x1 - x0) * alpha;
      double wy = y0 + (y1 - y0) * alpha;
      double wz = z0 + (z1 - z0) * alpha;

      double vx = world_to_vox(0,0)*wx + world_to_vox(0,1)*wy + world_to_vox(0,2)*wz + world_to_vox(0,3);
      double vy = world_to_vox(1,0)*wx + world_to_vox(1,1)*wy + world_to_vox(1,2)*wz + world_to_vox(1,3);
      double vz = world_to_vox(2,0)*wx + world_to_vox(2,1)*wy + world_to_vox(2,2)*wz + world_to_vox(2,3);

      if (vx < -0.5 || vy < -0.5 || vz < -0.5 || vx > nx - 0.5 || vy > ny - 0.5 || vz > nz - 0.5) {
        continue;
      }

      int xfloor = static_cast<int>(std::floor(vx));
      int yfloor = static_cast<int>(std::floor(vy));
      int zfloor = static_cast<int>(std::floor(vz));
      double fx = vx - xfloor, fy = vy - yfloor, fz = vz - zfloor;
      double wxs[2] = {1.0 - fx, fx};
      double wys[2] = {1.0 - fy, fy};
      double wzs[2] = {1.0 - fz, fz};

      for (int iz = 0; iz < 2; ++iz) {
        for (int iy = 0; iy < 2; ++iy) {
          for (int ix = 0; ix < 2; ++ix) {
            int xi = std::min(std::max(xfloor + ix, 0), nx - 1);
            int yi = std::min(std::max(yfloor + iy, 0), ny - 1);
            int zi = std::min(std::max(zfloor + iz, 0), nz - 1);
            const int lin = static_cast<int>(
              static_cast<R_xlen_t>(xi) +
              static_cast<R_xlen_t>(yi) * nx +
              static_cast<R_xlen_t>(zi) * nx * ny
            );
            if (use_mask && (mask[lin] == FALSE || mask[lin] == NA_LOGICAL)) {
              continue;
            }
            double w = wxs[ix] * wys[iy] * wzs[iz];
            accum[lin] += w;
          }
        }
      }
    } // s

    double total = 0.0;
    for (const auto& kv : accum) total += kv.second;
    if (total <= 0.0) continue;
    for (const auto& kv : accum) {
      rows.push_back(i + 1);
      cols.push_back(kv.first + 1);
      vals.push_back(kv.second / total);
    }
  }
#endif

  return Rcpp::List::create(
    Rcpp::Named("rows") = rows,
    Rcpp::Named("cols") = cols,
    Rcpp::Named("vals") = vals
  );
}
