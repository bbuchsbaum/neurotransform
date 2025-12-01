# Design Review Feedback

## Overall Assessment

> You're actually in really good shape. The 22+6 surface basically does what the Vision + Dev Guide promise: morphisms as first-class things, a tiny verb set, RAS-centric geometry, warp loader registry, and no graph / projector creep.

No gaping holes, but places where:
- The docs and API don't fully line up
- Some design decisions are implied but not nailed down

---

## 1. Composition & Paths: Story Almost There, Needs Nailing

### Current State
- `transform(morphism, coords)` – pullback for a single morphism
- `transform_path(path, coords)` – efficient application of a sequence
- `compose(f, g)` – described as returning g ∘ f
- Paths represented as `list(morphism, ...)` (e.g., `list(aff, warp)`)

### Missing Bits

#### 1a. What does `compose(f, g)` actually return?

Docs say "Compose two morphisms: g ∘ f (f first, then g)" but don't say:
- Is the result itself a Morphism?
- Or just a `list(f, g)` (i.e., a path)?
- Or a specialized `CompositeMorphism` / `MorphismPath` class?

**Suggestion (minimalist path):**
- Define `compose(f, g)` to return a `MorphismPath` object, which is just a list of morphisms with a light S3/S7 class
- Have `transform()` dispatch on `MorphismPath` by delegating to `transform_path()`
- `transform_path()` can accept either a bare list or a `MorphismPath` (but always returns coords, never a morphism)

This keeps the category flair ("morphisms compose") while reusing the one efficient implementation (`transform_path`) without introducing deep new machinery.

#### 1b. Direction & Order in `transform_path()`

Example given:
```r
# Native → MNI (aff) then MNI → template (warp)
aff  <- Affine3DMorphism(source = "native", target = "mni", matrix = mat1)
warp <- Warp3DMorphism(source = "mni", target = "template", warp_path = "warp.nii.gz")

path <- list(aff, warp)
template_coords <- matrix(c(0, 0, 0), ncol = 3)
native_coords   <- transform_path(path, template_coords)
```

Given:
- Base morphism semantics: matrix maps target → source
- `transform(m, coords)` is documented as a pullback: target coords → source coords

Then for path `[aff: native→mni, warp: mni→template]`, the correct pullback from template → native means:
```
native <- aff( warp( template ) )
```

i.e., `transform_path()` must apply from **last to first**.

**Action items:**
- Explicitly write in docs:
  > Given a path f : A→B, g : B→C, `transform_path(list(f, g), coords_in_C)` returns coords_in_A by applying g then f (last to first).
- Add one diagram with arrows and pullbacks

#### 1c. Path Validation Helper (tiny but high-value)

We already have `source_of()` / `target_of()` and domain hashes. Almost everyone will want to check that a path is composable before calling `transform_path`.

**Suggestion:** Add one helper:
```r
validate_path(path)  # errors on mismatched source/target, returns invisible(path) otherwise
```
or
```r
is_valid_path(path)  # TRUE/FALSE
```

Keeps core verbs at 4 while centralizing domain-consistency logic.

---

## 2. Inverses, Adjoints, and Non-Invertible Morphisms

### Current State

Base Morphism slot layout anticipates a richer story:
```
Morphism (virtual)
├── inverse_type: character ("exact", "approximate", "adjoint", "none")
├── inverse_quality: numeric (0-1)
└── ...
```

Vision doc says:
> apply transforms … and expose adjoint/metadata for non-invertible cases

But public verbs only give:
- `invert(morphism)` – errors if non-invertible
- `is_invertible(morphism)` – TRUE if `invert()` will succeed

### Missing: Adjoint vs. Non-Invertible

For `VolToSurfMorphism` and parts of `SurfToSurfMorphism`, an adjoint (e.g., back-projection) makes sense but isn't a true inverse.

Right now the only public hook is `inverse_type` metadata; there's no way to actually *use* an adjoint.

**Option A (add one verb):**
```r
adjoint(morphism)
```
- Returns a morphism when `inverse_type == "adjoint"`, errors otherwise
- `is_invertible()` stays about exact inverses
- `adjoint()` is explicitly "not quite inverse" but useful

**Option B (extend `invert()`):**
```r
invert(morphism, mode = c("exact", "approximate", "adjoint", "any"))
```
- Default `mode = "exact"` keeps current semantics
- `mode = "adjoint"` uses adjoint if present; errors otherwise

**Recommendation:** Option A keeps verb set conceptually cleaner:
- `invert` = honest isomorphism
- `adjoint` = "generalized inverse"

### What does `is_invertible()` mean exactly?

With the `inverse_type` slot, decide:
- Is `is_invertible()` TRUE only when `inverse_type == "exact"`?
- Or also when `inverse_type == "approximate"`?

**Strong recommendation:** Only for "exact". For approximate inverses/adjoints, rely on `adjoint()` / extra parameters to `invert()`.

Make this explicit in docs to avoid surprises.

---

## 3. Domain Hashes and "What is a Domain?"

### Current State

Docs consistently say:
- `source_of(m)` / `target_of(m)` returns a domain hash
- Vision & Dev Guide emphasize: no Domain/Geometry classes here; that lives in neurofunctor

### Action Items

#### 3a. Document the contract for domain hashes

Even if neurotransform never creates them, spell out assumptions:
- They're **opaque strings** from neurotransform's POV
- Two domains are considered identical iff their hashes compare equal
- Hashes should be stable under re-loading a dataset but may be different between datasets

Add one short subsection in API_DESIGN.md: **"Domain hashes: what we assume"**

#### 3b. Note explicitly that geometry/orientation lives outside

Already say "No Domain/Geometry/BrainGraph/Projector/Field" in API_DESIGN. Add one small note near intro:

> If you want to know *what* a domain is (MNI152 2mm, native T1, fsaverage, etc.), you use higher-level packages (e.g., neurofunctor). Here we only promise consistent behavior given consistent hashes.

Keeps people from smuggling domain semantics back into neurotransform.

---

## 4. Coordinate Utilities: Coverage and Expectations

### Current Exported Utilities
- `ras_to_lps()`, `lps_to_ras()`
- `tkras_to_ras()`, `ras_to_tkras()`
- `apply_affine()`, `invert_affine()`, `compose_affines()`

### Mismatch

Vision / Dev Guide hint at broader set ("FSL/AFNI flips"). That's the only real mismatch noticed.

**Option 1: Clarify FSL/AFNI conventions handled inside loaders**

If plan is:
- FSL/AFNI warps loaded via `register_loader()` implementations
- Those loaders internally normalize everything to RAS mm

…then say that explicitly and remove "FSL/AFNI flips" from coordinate-utility marketing bullets.

**Option 2: Add one or two general helpers**

If users will call them directly, minimal additions:
```r
# classic NIfTI-style helpers
ijk_to_ras(coords, affine)  # apply qform/sform
ras_to_ijk(coords, affine)
```

Keeps surface area small but covers most common "real world" conversion task. Already have `apply_affine()`, so these could just be examples built on it.

**Either way:** Make Vision / Dev Guide text match actual exported list.

---

## 5. Warp Loader Registry Ergonomics

### Current Public API
- `register_loader(name, fn)`
- `get_loader(name)`

### Missing in Docs

**What shape must a loader function have?**

Example contract worth writing down explicitly:
```r
# Loader function contract:
my_loader <- function(path) {

  # Returns: list(array, dim, world_to_vox)
  # - array: numeric vector of displacement values (flattened 4D: 3 × X × Y × Z)
  # - dim: integer vector c(X, Y, Z)
  # - world_to_vox: 4×4 matrix mapping world coords to voxel indices
}
register_loader("my_format", my_loader)
```

Even a short bullet list in API_DESIGN.md ("Loader function contract") prevents guessing.

**Optional tiny extra:** `list_loaders()` could be nice but not essential. Registry is already small and focused.

---

## 6. Tiny Nits & Consistency Things

### Naming Drift Between Docs

Vision/Dev Guide still talk about `transform_coords` and `register_warp_loader` as primaries; API_DESIGN has `transform()` + aliases.

**Action:** Update prose to use new names; describe old ones as compatibility aliases only.

### Cost & inverse_quality Slots

Nicely defined but no dedicated accessor. That's fine; people can inspect slots directly.

If meant to be used by external packages (e.g., path planners in neurofunctor), consider:
- Documenting in a "Morphism metadata" section
- Maybe adding `morphism_meta(m)` that returns list with cost, inverse_type, inverse_quality

---

## What NOT to Add (Yet)

Explicitly not for v0, even though tempting:
- Jacobian / gradient transforms (`transform_jacobian`, etc.)
- Volume/surface data resampling (`warp_volume`, `sample_surface`)
- Any Domain/Geometry/Graph/Projector abstraction
- On-disk warp composition / writing new warps

Those belong either in higher-level packages or in "vNext" once core kernel has baked.

---

## TL;DR

The 22 exports + 6 aliases are basically sufficient for the "geometry kernel" described.

**The only real "holes" are:**
1. `compose()` semantics & relation to `transform_path` / path representation
2. The adjoint/approximate inverse story, given `inverse_type` slot and Vision promises
3. Documentation around domain hashes, loader contracts, and coordinate helper scope

> If you tighten those up, I think you'll have exactly what you're aiming for: a small, compositional, category-flavored transform kernel that feels inevitable rather than "clever."
