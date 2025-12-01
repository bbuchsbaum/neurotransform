# neurotransform API Design

## Goal

The best neuroimaging transform library on earth. Most elegant, most beautiful, with a category-theoretic flair that illuminates rather than intimidates.

---

## Design Philosophy

**Morphisms are the stars.** A transform is a morphism between coordinate spaces. You compose them. You invert them (when you can). You apply them to coordinates. That's it.

**Verbs are few and powerful.** We don't need twenty functions when five will do.

**Nouns are clear.** `Morphism`, `Affine3D`, `Warp3D`, `VolToSurf`, `SurfToSurf`. You know what they are.

**Convention over configuration.** RAS millimeters internally. Loaders handle the mess of FSL/AFNI/ANTs/FreeSurfer conventions so you don't have to.

**Domain-agnostic.** If you want to know *what* a domain is (MNI152 2mm, native T1, fsaverage, etc.), use higher-level packages like neurofunctor. Here we only promise consistent behavior given consistent hashes.

---

## Domain Hashes: What We Assume

neurotransform doesn't create or interpret domain hashes—it just threads them through morphisms for consistency checking.

**Contract:**
- Domain hashes are **opaque strings** from our point of view
- Two domains are identical **iff their hashes compare equal**
- Hashes should be **stable** under re-loading a dataset
- Hashes may differ between datasets (that's fine—different data, different domains)
- We never look inside a hash or try to parse it

This keeps us cleanly separated from the Domain/Geometry layer in neurofunctor.

---

## The Verbs (Functions)

### Core Transform Operations

| Verb | Signature | What it does |
|------|-----------|--------------|
| **`transform`** | `transform(morphism, coords)` | Apply a single morphism to coordinates (pullback) |
| **`transform_path`** | `transform_path(path, coords)` | Apply a sequence of morphisms efficiently |
| **`compose`** | `compose(f, g)` | Compose two morphisms into a `MorphismPath` |
| **`invert`** | `invert(morphism)` | Get the inverse morphism (errors if non-invertible) |
| **`adjoint`** | `adjoint(morphism)` | Get the adjoint (generalized inverse) for non-invertible morphisms |

### Paths and Composition

| Verb | Signature | What it does |
|------|-----------|--------------|
| **`is_valid_path`** | `is_valid_path(path)` | TRUE if path is composable (domains chain correctly) |
| **`validate_path`** | `validate_path(path)` | Errors on invalid path, returns `invisible(path)` otherwise |

### Introspection

| Verb | Signature | What it does |
|------|-----------|--------------|
| **`morphism_kind`** | `morphism_kind(m)` | Returns: "identity", "affine3d", "warp3d", "vol2surf", "surf2surf" |
| **`is_invertible`** | `is_invertible(m)` | TRUE only if `inverse_type == "exact"` |
| **`has_adjoint`** | `has_adjoint(m)` | TRUE if `adjoint()` will succeed |
| **`source_of`** | `source_of(m)` | Source domain hash |
| **`target_of`** | `target_of(m)` | Target domain hash |

### Loader Registry

| Verb | Signature | What it does |
|------|-----------|--------------|
| **`register_loader`** | `register_loader(name, fn)` | Register a warp loader function |
| **`get_loader`** | `get_loader(name)` | Retrieve a registered loader |
| **`list_loaders`** | `list_loaders()` | List all registered loader names |

#### Loader Function Contract

A loader function must have this signature:

```r
my_loader <- function(path) {

  # Load warp field from disk and return:
 list(
    array = numeric_vector,    # Displacement values, flattened (3 × X × Y × Z)
    dim = integer_vector,      # c(X, Y, Z) - spatial dimensions
    world_to_vox = matrix_4x4  # Maps world coords (mm) to voxel indices
  )
}
register_loader("my_format", my_loader)
```

The default "rnifti" loader uses RNifti to load NIfTI displacement fields.

### Coordinate Utilities

| Verb | Signature | What it does |
|------|-----------|--------------|
| **`ras_to_lps`** | `ras_to_lps(coords)` | RAS → LPS (ANTs/ITK convention) |
| **`lps_to_ras`** | `lps_to_ras(coords)` | LPS → RAS |
| **`tkras_to_ras`** | `tkras_to_ras(coords, c_ras)` | FreeSurfer tkRAS → scanner RAS |
| **`ras_to_tkras`** | `ras_to_tkras(coords, c_ras)` | Scanner RAS → FreeSurfer tkRAS |
| **`apply_affine`** | `apply_affine(coords, mat)` | Apply 4×4 affine to coordinates |
| **`invert_affine`** | `invert_affine(mat)` | Invert a 4×4 affine matrix |
| **`compose_affines`** | `compose_affines(A, B)` | Compose affines: B(A(x)) = BA |

---

## Composition, Paths, and the Pullback Direction

### The Pullback Model

All transforms in neurotransform follow the **pullback** convention:

> Given a morphism `f : A → B`, calling `transform(f, coords_in_B)` returns `coords_in_A`.

This matches the standard neuroimaging convention where a "forward" warp (e.g., native→MNI) is applied by giving it MNI coordinates and getting back native coordinates. The warp tells you "where in the source to sample."

### Path Application Order

For a path of morphisms `[f : A→B, g : B→C]`:

```
transform_path(list(f, g), coords_in_C)  →  coords_in_A
```

Internally, this applies **last to first**: `f(g(coords_in_C))`.

```
     f           g
A ←――――― B ←――――― C
         ↑
    pullback direction
```

**Rule:** Given a path `f : A→B, g : B→C`, `transform_path(list(f, g), coords_in_C)` returns `coords_in_A` by applying g then f (last to first).

### What `compose()` Returns

`compose(f, g)` returns a `MorphismPath` object—a lightweight wrapper around `list(f, g)` with class semantics.

```r
path <- compose(f, g)      # Returns MorphismPath
class(path)                # "MorphismPath"
transform(path, coords)    # Dispatches to transform_path internally
```

This keeps the category flair ("morphisms compose") while reusing `transform_path` under the hood.

### Path Validation

Before applying a path, you may want to check that domains chain correctly:

```r
f <- Affine3DMorphism(source = "A", target = "B", matrix = mat1)
g <- Warp3DMorphism(source = "B", target = "C", warp_path = "warp.nii.gz", warp_type = "ants")
h <- Affine3DMorphism(source = "X", target = "Y", matrix = mat2)  # Doesn't chain!

is_valid_path(list(f, g))  # TRUE: A→B→C
is_valid_path(list(f, h))  # FALSE: B ≠ X

validate_path(list(f, g))  # Returns invisible(path)
validate_path(list(f, h))  # Error: target of morphism 1 ("B") != source of morphism 2 ("X")
```

---

## Inverses vs. Adjoints

### The Distinction

| Concept | Meaning | When Available | Use Case |
|---------|---------|----------------|----------|
| **Inverse** | Exact geometric inverse: `f⁻¹(f(x)) = x` | `inverse_type == "exact"` | Affines, warps with inverse_path |
| **Adjoint** | "Best effort" reverse mapping, not a true inverse | `inverse_type == "adjoint"` | VolToSurf backprojection |

### `is_invertible()` vs `has_adjoint()`

```r
is_invertible(affine)       # TRUE  - exact inverse exists
is_invertible(warp)         # TRUE if inverse_path provided, FALSE otherwise
is_invertible(vol2surf)     # FALSE - no geometric inverse

has_adjoint(affine)         # TRUE  - inverse is also an adjoint
has_adjoint(vol2surf)       # TRUE  - adjoint exists for backprojection
has_adjoint(warp_no_inv)    # FALSE - neither inverse nor adjoint
```

### Using `invert()` and `adjoint()`

```r
# Exact inverse
aff_inv <- invert(affine)           # Works
warp_inv <- invert(warp_with_inv)   # Works

# Adjoint for non-invertible morphisms
v2s_adj <- adjoint(vol2surf)        # Works - returns backprojection morphism

# Errors
invert(vol2surf)                    # Error: not invertible (use adjoint())
adjoint(warp_no_inverse)            # Error: no adjoint available
```

### The `inverse_type` Slot

| Value | `is_invertible()` | `has_adjoint()` | `invert()` | `adjoint()` |
|-------|-------------------|-----------------|------------|-------------|
| `"exact"` | TRUE | TRUE | Works | Works (same as invert) |
| `"approximate"` | FALSE | TRUE | Error | Works |
| `"adjoint"` | FALSE | TRUE | Error | Works |
| `"none"` | FALSE | FALSE | Error | Error |

---

## The Nouns (Classes)

### Base Class

```
Morphism (virtual)
├── id: character
├── source: character (domain hash)
├── target: character (domain hash)
├── kind: character
├── inverse_type: character ("exact", "approximate", "adjoint", "none")
├── inverse_quality: numeric (0-1)
└── cost: numeric
```

### Concrete Morphism Types

```
IdentityMorphism < Morphism
└── (no additional slots)

Affine3DMorphism < Morphism
└── matrix: 4×4 numeric

Warp3DMorphism < Morphism
├── warp_path: character
├── warp_type: character ("ants", "fsl", "afni", "freesurfer")
├── def_type: character ("relative", "absolute")
└── inverse_path: character

VolToSurfMorphism < Morphism
├── method: character ("trilinear", "ribbon", "nearest")
├── ribbon_inner: character (path)
├── ribbon_outer: character (path)
└── n_ribbon_samples: integer

SurfToSurfMorphism < Morphism
├── method: character ("sphere", "area", "sulc")
├── source_sphere: character (path)
└── target_sphere: character (path)
```

---

## Constructors

Clean, predictable, no surprises.

```r
# Identity: domain stays the same
IdentityMorphism(domain)

# Affine: 4×4 matrix maps target coords → source coords
Affine3DMorphism(source, target, matrix, cost = 1, method_tag = "anatomical")

# Warp: displacement field on disk
Warp3DMorphism(source, target, warp_path,
               warp_type = c("ants", "fsl", "afni", "freesurfer"),
               def_type = c("relative", "absolute"),
               inverse_path = NULL, cost = 1.5)

# Volume to surface sampling
VolToSurfMorphism(source, target,
                  method = c("trilinear", "ribbon", "nearest"),
                  ribbon_inner = NULL, ribbon_outer = NULL,
                  n_ribbon_samples = 6L, cost = 2)

# Surface to surface (sphere-based)
SurfToSurfMorphism(source, target,
                   method = c("sphere", "area", "sulc"),
                   source_sphere = NULL, target_sphere = NULL, cost = 1)
```

---

## Compatibility with neurofunctor

We keep these **identical** to neurofunctor's current API:

| neurofunctor | neurotransform | Status |
|--------------|----------------|--------|
| `IdentityMorphism()` | `IdentityMorphism()` | Same |
| `Affine3DMorphism()` | `Affine3DMorphism()` | Same |
| `Warp3DMorphism()` | `Warp3DMorphism()` | Same |
| `VolToSurfMorphism()` | `VolToSurfMorphism()` | Same |
| `SurfToSurfMorphism()` | `SurfToSurfMorphism()` | Same |
| `transform_path()` | `transform_path()` | Same |
| `transform_coords()` | `transform()` | **Renamed** (generic method stays for compat) |
| `compose()` | `compose()` | Same |
| `invert()` | `invert()` | Same |
| `morphism_kind()` | `morphism_kind()` | Same |
| `source_domain()` | `source_of()` | **Renamed** (alias provided) |
| `target_domain()` | `target_of()` | **Renamed** (alias provided) |
| `register_warp_loader()` | `register_loader()` | **Renamed** (alias provided) |
| `lps_to_ras()` | `lps_to_ras()` | Same |
| `ras_to_lps()` | `ras_to_lps()` | Same |
| `tkras_to_scanner_ras()` | `tkras_to_ras()` | **Simplified name** (alias provided) |
| `scanner_ras_to_tkras()` | `ras_to_tkras()` | **Simplified name** (alias provided) |

### Compatibility Aliases

For seamless migration, we provide aliases:

```r
# These work but emit no warnings (permanent aliases)
transform_coords <- transform
source_domain <- source_of
target_domain <- target_of
register_warp_loader <- register_loader
tkras_to_scanner_ras <- tkras_to_ras
scanner_ras_to_tkras <- ras_to_tkras
```

---

## What's New (Beyond neurofunctor)

### `is_invertible()`

```r
is_invertible(affine)      # TRUE - exact inverse
is_invertible(warp)        # TRUE if inverse_path provided
is_invertible(vol2surf)    # FALSE - adjoint only
```

### Cleaner names

- `transform()` instead of `transform_coords()` as the primary verb
- `source_of()` / `target_of()` - reads better: "source of this morphism"
- `register_loader()` - we're not just loading warps anymore

### Better error messages

```r
invert(vol2surf)
# Error: VolToSurfMorphism is not invertible (has adjoint only).
# Use adjoint() at the projector level for backprojection.
```

---

## Tool Convention Handling

This is where the real complexity lives. Each neuroimaging tool has its own coordinate conventions, warp representations, and affine semantics. neurotransform handles this so you don't have to.

### The Problem

| Tool | Coordinate System | Warp Type | Affine Format | Quirks |
|------|------------------|-----------|---------------|--------|
| **ANTs** | LPS | Relative displacement | ITK .txt or .mat | Displacement in LPS coords |
| **FSL** | Scaled voxels | Relative OR Absolute | FLIRT .mat (4×4) | Needs source/ref geometry to interpret |
| **AFNI** | RAI | Relative displacement | .aff12.1D (3×4) | Z-axis flipped vs RAS |
| **FreeSurfer** | tkRAS | N/A (surface) | N/A | Offset by c_ras from scanner RAS |

### Our Solution: Convention-Aware Warp Application

When you create a `Warp3DMorphism` with `warp_type`, we handle the convention internally:

```r
# ANTs warp (LPS displacement field)
warp_ants <- Warp3DMorphism(
  source = "native", target = "mni",
  warp_path = "ants_warp.nii.gz",
  warp_type = "ants"        # We flip LPS↔RAS internally
)

# FSL FNIRT warp (may be absolute coordinates!)
warp_fsl <- Warp3DMorphism(
  source = "native", target = "mni",
  warp_path = "fnirt_warp.nii.gz",
  warp_type = "fsl",
  def_type = "absolute"     # We convert absolute→relative internally
)

# AFNI 3dQwarp (RAI displacement)
warp_afni <- Warp3DMorphism(
  source = "native", target = "mni",
  warp_path = "afni_warp.nii.gz",
  warp_type = "afni"        # We flip RAI↔RAS internally
)
```

### Internal Processing by Tool

#### ANTs (LPS)
```
Input coords (RAS) → flip to LPS → sample displacement → add to coords → flip back to RAS
```

#### FSL (Scaled Voxels + Absolute/Relative)
```
If absolute: convert field to displacement (subtract voxel coords)
Then: standard displacement sampling
```
Plus: `detect_warp_type()` heuristic for FNIRT fields that don't declare their type.

#### AFNI (RAI)
```
Input coords (RAS) → flip Z to RAI → sample displacement → add to coords → flip Z back to RAS
```

### Coordinate Conversion Utilities

For working with tool-specific coordinates directly:

```r
# ANTs/ITK ↔ Internal
lps_to_ras(coords)
ras_to_lps(coords)

# AFNI ↔ Internal
rai_to_ras(coords)
ras_to_rai(coords)

# FreeSurfer ↔ Internal
tkras_to_ras(coords, c_ras)
ras_to_tkras(coords, c_ras)

# FSL scaled coords ↔ World mm
fsl_to_world(coords, src_affine, ref_affine)
world_to_fsl(coords, src_affine, ref_affine)
```

### Affine Loading Helpers

Tool-specific affine matrix readers that return internal (RAS) convention:

```r
# FSL FLIRT matrix → internal affine
read_flirt(mat_path, source_affine, ref_affine)

# AFNI .aff12.1D → internal affine
read_afni_aff12(path)

# ANTs/ITK .txt → internal affine
read_itk_affine(path)
```

### Warp Type Detection

For FSL FNIRT warps that might be absolute or relative:

```r
# Heuristic detection
detect_warp_type(warp_path)
# Returns: "relative" or "absolute"

# How it works:
# - Sample random voxels
# - If field values correlate highly with world coords → absolute
# - If field values are small displacements → relative
```

---

## What We Explicitly Don't Have

- No `Domain` class
- No `Geometry` class
- No `BrainGraph`
- No `Projector`
- No `Field`
- No pathfinding
- No ingestion pipelines

Those belong in neurofunctor. We're the kernel.

---

## Usage Examples

### Basic transform

```r
library(neurotransform)

# Create an affine morphism
aff <- Affine3DMorphism(
  source = "native_hash",
  target = "mni_hash",
  matrix = my_4x4_matrix
)

# Transform some coordinates
coords_mni <- matrix(c(0, 0, 0, 10, 20, 30), ncol = 3, byrow = TRUE)
coords_native <- transform(aff, coords_mni)
```

### Compose transforms

```r
# Native → MNI (affine) then MNI → template (warp)
aff <- Affine3DMorphism(source = "native", target = "mni", matrix = mat1)
warp <- Warp3DMorphism(source = "mni", target = "template",
                        warp_path = "warp.nii.gz", warp_type = "ants")

# Compose into a path
path <- list(aff, warp)
template_coords <- matrix(c(0, 0, 0), ncol = 3)
native_coords <- transform_path(path, template_coords)
```
### Inversion

```r
# Affines invert exactly
aff_inv <- invert(aff)

# Warps invert if inverse_path was provided
warp_with_inv <- Warp3DMorphism(
  source = "mni", target = "template",
  warp_path = "forward.nii.gz",
  inverse_path = "inverse.nii.gz",
  warp_type = "ants"
)
warp_inv <- invert(warp_with_inv)  # Works

# Warps without inverse_path error
warp_no_inv <- Warp3DMorphism(source = "a", target = "b",
                               warp_path = "warp.nii.gz", warp_type = "fsl")
invert(warp_no_inv)  # Error: no inverse available
```

### Coordinate conventions

```r
# Coming from ANTs (LPS)?
ras_coords <- lps_to_ras(ants_coords)

# FreeSurfer surface?
ras_coords <- tkras_to_ras(surface_coords, c_ras = c(0, 0, 0))
```

---

## Summary: The Complete Public API

### Constructors (5)
- `IdentityMorphism()`
- `Affine3DMorphism()`
- `Warp3DMorphism()`
- `VolToSurfMorphism()`
- `SurfToSurfMorphism()`

### Transform Verbs (5)
- `transform()` — apply single morphism (pullback)
- `transform_path()` — apply sequence efficiently
- `compose()` — compose morphisms into `MorphismPath`
- `invert()` — get exact inverse (errors if unavailable)
- `adjoint()` — get generalized inverse for non-invertible morphisms

### Paths (2)
- `is_valid_path()` — check if path domains chain correctly
- `validate_path()` — error on invalid, return invisible otherwise

### Introspection (5)
- `morphism_kind()` — "identity", "affine3d", "warp3d", "vol2surf", "surf2surf"
- `is_invertible()` — TRUE only if `inverse_type == "exact"`
- `has_adjoint()` — TRUE if `adjoint()` will succeed
- `source_of()` — source domain hash
- `target_of()` — target domain hash

### Loader Registry (3)
- `register_loader()` — register a warp loader function
- `get_loader()` — retrieve a registered loader
- `list_loaders()` — list all registered loader names

### Coordinate Conventions (11)
- `ras_to_lps()` / `lps_to_ras()` — ANTs/ITK
- `ras_to_rai()` / `rai_to_ras()` — AFNI
- `tkras_to_ras()` / `ras_to_tkras()` — FreeSurfer
- `fsl_to_world()` / `world_to_fsl()` — FSL scaled coords
- `apply_affine()` / `invert_affine()` / `compose_affines()` — generic affine math

### Tool-Specific Readers (4)
- `read_flirt()` — FSL FLIRT .mat → internal affine
- `read_afni_aff12()` — AFNI .aff12.1D → internal affine
- `read_itk_affine()` — ANTs/ITK .txt → internal affine
- `detect_warp_type()` — heuristic for FSL FNIRT absolute vs relative

### Classes (2)
- `Morphism` (and subclasses) — the core transform abstraction
- `MorphismPath` — lightweight wrapper for composed morphism sequences

### Compatibility Aliases (6)
- `transform_coords` → `transform`
- `source_domain` → `source_of`
- `target_domain` → `target_of`
- `register_warp_loader` → `register_loader`
- `tkras_to_scanner_ras` → `tkras_to_ras`
- `scanner_ras_to_tkras` → `ras_to_tkras`

**Total: 35 exports + 6 aliases = 41 public symbols**

Still minimal. Every function earns its place. Handles every major neuroimaging tool's coordinate mess. Beautiful.
