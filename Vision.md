# neurotransform Vision

## What It Is

A self-contained, dependency-light library for geometric transforms in neuroimaging—affine and nonlinear (warp) mappings between coordinate systems—plus volume↔surface mapping primitives. It is **not** a projector/graph/runtime; it's the reusable transform kernel.

## Core Responsibilities

1. **Define morphism types** (identity, affine3d, warp3d, vol→surf, surf→surf) with stable S4/S7-like interfaces.

2. **Apply transforms** (`transform_path`, `transform_coords`), compose/invert where mathematically valid, and expose adjoint/metadata for non-invertible cases.

3. **Handle coordinate conventions** (RAS/LPS, tkRAS, FSL/AFNI flips), affine utilities, and warp sampling/composition with C++ acceleration.

4. **Provide a pluggable warp loader registry** (default RNifti) with morphism-level caching.

5. **Ship minimal, well-covered tests** for transform correctness; no reliance on neurofunctor graphs/domains.

## What It Explicitly Leaves Out

- No Domain/Geometry abstractions, projectors, fields, or graph/pathfinding.
- No dataset ingestion pipelines (fMRIPrep/FS/ANTs/FSL/AFNI); only generic loaders for warps/affines.
- No heavy neuroimaging I/O; NIfTI access only via optional RNifti loader.

## Niche in the Neuroimaging Ecosystem

- A **lightweight, embeddable transform kernel** for R-based neuroimaging tools needing accurate coordinate mappings without pulling in full processing stacks.

- **Bridges common tool outputs** (ANTs/FSL/AFNI/FreeSurfer) into a unified transform API, emphasizing correct conventions and invertibility metadata.

- Positioned as a **"geometry math" layer**: you bring your geometry/data, it gives you fast, correct coordinate transforms (including volume→surface).

- **Suitable as a dependency** for higher-level packages (e.g., neurofunctor, surface/ROI toolkits, QC utilities) and for scripting small transform tasks without larger frameworks.

## Domain

Neuroanatomical spatial transforms on MRI/CT-derived data:

- 3D affine and nonlinear warps in millimeter RAS space
- Surface coordinate mappings on standard meshes
- Inputs/outputs are numeric matrices (coords) and displacement fields; no dependence on specific domain classes
