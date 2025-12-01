# neurotransform Developer Guide

## Our Identity & Goal

**neurotransform aims to be the most elegant, most beautiful neuroimaging transform library on earth.**

We embrace a category-theoretic flair—not to intimidate, but because morphisms, composition, and invertibility are *exactly* the right mental model for spatial transforms. When you compose two transforms, you're composing morphisms. When you ask "can I go back?", you're asking about invertibility. The math isn't scary decoration; it's the cleanest way to think about what we're actually doing.

**Our principles:**

- **Elegant by design**: Clean abstractions that map naturally to the problem domain
- **Useful above all**: Category theory serves the user, never the other way around
- **Correct first, fast second**: Get the math right, then optimize
- **Minimal and complete**: Everything you need, nothing you don't

---

## Letter from neurofunctor Maintainer

> As the engineer shepherding neurofunctor, here's how I'd love to see neurotransform shaped so it's a delight to depend on.

### What "Beautiful" Looks Like

1. **Single mandate**: Be the geometry/transform kernel. Own morphisms (identity/affine/warp/vol→surf/surf→surf), coordinate conventions, transform application, and loader seams. No graph, projectors, domains, or ingestion pipelines.

2. **Stable, slim API**:
   - Constructors + generics: `IdentityMorphism`, `Affine3DMorphism`, `Warp3DMorphism`, `VolToSurfMorphism`, `SurfToSurfMorphism`; `transform_path`, `transform_coords`, `compose`, `invert`, `morphism_kind`
   - Loader registry: `register_warp_loader`, default RNifti loader; morphism-level caching
   - Coordinate helpers: RAS/LPS/tkRAS/FSL/AFNI flips; affine utils

3. **Minimal deps**: `methods`, `Rcpp`, `digest`; Suggests: `RNifti`, `testthat`. Keep neurosurf/neuroim2 out.

4. **Fast & correct**: C++ warp sampling/composition, affine chains; tests that pin relative/absolute warps, AFNI/FSL conventions, vol→surf sampling basics.

### Concrete Instructions for Painless Integration

| Guideline | Details |
|-----------|---------|
| **Keep domain-agnostic** | source/target remain character hashes; no Domain/Geometry classes |
| **Exports** | Export morphism constructors, `transform_path`, `morphism_kind`, `register_warp_loader`. Tag internals with `@keywords internal` |
| **Versioning/semver** | Bump minor for new features, major for slot/API changes; add NEWS entries for user-visible changes |
| **Compatibility shim** | Keep function names identical to neurofunctor; avoid gratuitous renames |
| **Tests we rely on** | Parity tests for `transform_path` vs sequential (affine, relative warp, absolute warp), AFNI/FSL convention checks, loader registry custom-loader test |
| **Build friendliness** | No external headers beyond Rcpp/Armadillo; `compileAttributes` committed; R CMD check clean on mac/win/linux |
| **Docs** | Short README section on coordinate conventions and loader registry; note RNifti is optional but default |
| **Performance sanity** | Keep warp load/apply sub-millisecond per-call for small warps; document any loader layout changes |

### How neurofunctor Will Wire In

1. Add `neurotransform` as an Import, re-export key constructors/`transform_path` for one release, then remove old shims
2. Projector/ingest will call `transform_path` and morphism constructors from neurotransform namespace
3. Loader stays default RNifti unless user overrides

---

## Bootstrap Plan

### Philosophy: Copy First, Then Diverge

We **copy** code from neurofunctor rather than move it. This is safer and lets both packages work independently during the transition. neurotransform grows its own identity; neurofunctor keeps working. Eventually neurofunctor will depend on neurotransform and shed its copies, but that's a later step.

### Scope: What Gets Copied from neurofunctor

#### R Code
- `morphism.R` — classes/constructors/generics
- `coordinates.R` — coordinate convention utilities
- `warp_transform.R` — warp transform logic
- `warp_chain.R` — chained warp composition
- `warp_loader.R` — loader registry
- `transform_path.R` — main path transform function
- Small utilities: `compute_hash`, `new_cache_env` if needed

#### C++ Code
- `rcpp_warp.cpp`
- `rcpp_apply.cpp`
- `rcpp_worldcoords.cpp`
- `rcpp_bary.cpp`
- `rcpp_ribbon.cpp`
- Generated `RcppExports.*`

#### Tests
- `transform_path` parity
- morphism kind helpers
- warp loader registry
- affine/warp transform correctness (relative & absolute)
- AFNI/FSL/ANTs minimal mocks

---

### Package Scaffold

#### DESCRIPTION
```
Package: neurotransform
Title: Geometric Transforms for Neuroimaging Data
Version: 0.1.0
Description: A self-contained, dependency-light library for geometric
    transforms in neuroimaging. Provides affine and nonlinear (warp) mappings
    between coordinate systems, plus volume-to-surface mapping primitives.
    Emphasizes correctness, elegance, and composability.
Imports:
    methods,
    Rcpp,
    digest
Suggests:
    RNifti,
    testthat (>= 3.0.0)
LinkingTo: Rcpp, RcppArmadillo
```

#### NAMESPACE
```r
useDynLib(neurotransform, .registration = TRUE)

# Morphism constructors
export(IdentityMorphism)
export(Affine3DMorphism)
export(Warp3DMorphism)
export(VolToSurfMorphism)
export(SurfToSurfMorphism)

# Core functions
export(transform_path)
export(transform_coords)
export(compose)
export(invert)
export(morphism_kind)

# Loader registry
export(register_warp_loader)

# Coordinate utilities
export(ras_to_lps)
export(lps_to_ras)
# ... etc
```

---

### API Alignment

| Component | Specification |
|-----------|---------------|
| **Morphism slots** | Keep identical semantics; source/target = domain hash strings |
| **Default loader** | RNifti via registry: `register_warp_loader("rnifti", load_warp_rnifti)` |
| **Public path function** | `transform_path(path, coords)` |
| **Generics** | `transform_coords`, `compose`, `invert`, `morphism_kind`, etc. |

---

### Testing Strategy

Unit tests in `tests/testthat/`:

```
test-morphism-constructors.R    # morphism constructors & kind helpers
test-transform-path.R           # transform_path vs sequential
                                # (affine-only, affine+warp, absolute warp)
test-loader-registry.R          # loader registry custom loader
test-warp-conventions.R         # AFNI/FSL absolute/relative warp mock sampling
```

**Critical**: Tests must have no dependency on neurofunctor objects.

---

### Future: neurofunctor Adaptation (Later Phase)

Once neurotransform is stable and tested:

1. **Add dependency**: `Imports: neurotransform (>= 0.1.0)`

2. **Re-export for compatibility**: One-release window with re-exports
   ```r
   # compat.R (temporary)
   #' @export
   #' @importFrom neurotransform transform_path
   transform_path <- neurotransform::transform_path
   ```

3. **Remove duplicated code**: Delete transform R/src files from neurofunctor

4. **Update projector/ingest**: Use imported `transform_path` and morphism constructors

---

### CI/Build Requirements

- Add CI for neurotransform (R CMD check + mac/win/linux builds)
- Ensure compiled code is standalone (no neurosurf/neuroim2 headers)
- `Rcpp::compileAttributes()` output committed

---

### Timeline

| Phase | Tasks |
|-------|-------|
| **Phase 1** | Scaffold package, copy R files, adapt namespaces, run unit tests |
| **Phase 2** | Copy C++ + compileAttributes, fix includes, R CMD check clean |
| **Phase 3** | Documentation, README, coordinate convention guide |
| **Phase 4** | (Later) neurofunctor wiring, re-exports, deprecate old code |

---

## Summary Checklist

- [ ] Package scaffold with correct DESCRIPTION/NAMESPACE
- [ ] All morphism classes copied with identical slots
- [ ] `transform_path` and generics exported
- [ ] Loader registry with default RNifti
- [ ] C++ code copied with standalone headers
- [ ] Unit tests passing without neurofunctor deps
- [ ] R CMD check clean on all platforms
- [ ] README with coordinate conventions docs
- [ ] NEWS file documenting v0.1.0
- [ ] (Later) neurofunctor updated with imports and re-exports
