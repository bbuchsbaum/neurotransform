# neurotransform

`neurotransform` is a lightweight R package for geometric transforms in neuroimaging.
It provides a small, consistent transform kernel for:

- affine transforms,
- nonlinear warp fields,
- volume-to-surface sampling,
- surface-to-surface mappings,
- coordinate-convention utilities for common neuroimaging tools.

The package is designed around a single idea: spatial transforms should be first-class objects with clear composition, inversion, and Jacobian semantics.

## Status

`neurotransform` is under active development.
The current API already covers the main transform types and resampling paths, but the package should still be treated as early-stage software.

## Installation

From GitHub:

```r
install.packages("remotes")
remotes::install_github("bbuchsbaum/neurotransform")
```

Core dependencies are intentionally small.
Some file formats and workflows use suggested packages such as `RNifti`, `hdf5r`, `Matrix`, `knitr`, and `rmarkdown`.

## What The Package Covers

- `IdentityMorphism`, `Affine3DMorphism`, `Warp3DMorphism`, `VolToSurfMorphism`, and `SurfToSurfMorphism`
- pullback-style coordinate transforms via `transform()`
- transform composition with `compose()`
- inversion and adjoint access where mathematically available
- Jacobian matrices and determinants for affine and warp morphisms
- repeated volume resampling through cached `ResamplingPlan`s
- import/export helpers for ANTs, FSL, AFNI, dense fields, and related affine utilities

## Quick Start

```r
library(neurotransform)

# An affine mapping from "native" to "template" space.
# With pullback semantics, it maps template coordinates back into native space.
mat <- diag(4)
mat[1:3, 4] <- c(10, 20, 30)

aff <- Affine3DMorphism(
  source = "native",
  target = "template",
  matrix = mat
)

coords_template <- matrix(c(
  0, 0, 0,
  5, 5, 5
), ncol = 3, byrow = TRUE)

coords_native <- transform(aff, coords_template)
coords_native
```

You can compose transforms directly:

```r
aff1 <- Affine3DMorphism("native", "intermediate", diag(4))
aff2 <- Affine3DMorphism("intermediate", "template", diag(4))

path <- compose(aff1, aff2)
transform(path, coords_template)
```

And resample a volume into a target grid:

```r
vol <- array(runif(27), dim = c(3, 3, 3))
grid <- grid_spec(dims = c(3, 3, 3), affine = diag(4), domain = "template")

out <- resample_volume(
  data = vol,
  morphism = aff,
  target = grid,
  method = "linear"
)
```

## Coordinate Semantics

The package uses pullback semantics throughout.
A morphism from domain A to domain B answers the question:

“Given coordinates in B, where should I sample in A?”

That convention makes composition and resampling consistent, especially for image interpolation.

The package also includes helpers for the common convention mismatches across:

- ANTs,
- FSL,
- AFNI,
- FreeSurfer / tkRAS,
- RAS/LPS conversions.

## Documentation

- Intro vignette: [vignettes/introduction.Rmd](./vignettes/introduction.Rmd)
- Developer notes: [DEVELOPER_GUIDE.md](./DEVELOPER_GUIDE.md)
- Design notes: [Vision.md](./Vision.md)
- Pkgdown site: https://bbuchsbaum.github.io/neurotransform/

## Why This Exists

Most neuroimaging projects need correct coordinate transforms, but not an entire processing framework.
`neurotransform` aims to be the reusable geometry layer: small enough to embed, explicit enough to reason about, and tested closely enough to trust for affine/warp math and resampling.

## License

MIT. See [LICENSE](./LICENSE).
