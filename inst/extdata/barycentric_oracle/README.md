# Ordinary barycentric projection oracle

`oracle.json` contains complete impulse-response weights generated independently
with Connectome Workbench 2.2.1, `wb_command -metric-resample ... BARYCENTRIC`.
The fixture records the Workbench version/commit, executable and generator
SHA-256 hashes, exact commands, and GIFTI file hashes. No neurotransform code is
used to generate expected values.

The three synthetic cases each have 66 target vertices: an octahedron with
shared edges/vertices, an octahedron with asymmetric rotated queries, and an
irregular 18-vertex sphere. Each input metric column is one vertex impulse, so
the output identifies all interpolation weights, not just constant-field
preservation. Inputs are stored at GIFTI float32 precision. The absolute weight
tolerance is 2e-6, allowing float32 arithmetic and radius normalization in
Workbench; analytic double-precision tests use 1e-12.

Regenerate with Python/numpy and an official Workbench installation:

```sh
python3 tools/generate_barycentric_oracle.py /path/to/wb_command /tmp/new-oracle
```

The generator retains GIFTI inputs and outputs in that new directory. Copy its
`oracle.json` here only after reviewing version and numerical differences.

Spherical plans normalize both meshes to a common radius and use Euclidean
closest-point projection onto triangles, including edges/vertices. This is
ordinary barycentric interpolation, not radial intersection or ADAP_BARY_AREA.
The generic `SurfToSurfMorphism` explicitly retains interior-only orthogonal
projection and returns `NA` for unsupported queries. Degenerate faces are not
valid support. Spherical plans error if there is no valid triangle support.

Workbench's implementation can be inspected in
[SurfaceResamplingHelper](https://github.com/Washington-University/workbench/blob/01164ffa47f2778088bd6ff472ec9cd9a57f5b42/src/Files/SurfaceResamplingHelper.cxx)
and [SignedDistanceHelper](https://github.com/Washington-University/workbench/blob/01164ffa47f2778088bd6ff472ec9cd9a57f5b42/src/Files/SignedDistanceHelper.cxx).
Its [metric-resample documentation](https://www.humanconnectome.org/software/workbench-command/-metric-resample)
distinguishes BARYCENTRIC from ADAP_BARY_AREA.

`tools/compare_barycentric_templates.py` additionally compares six nonconstant
metrics at 64 seeded target vertices per direction/hemisphere on local
fsaverage 164k and registered fsLR 32k spheres. It records geometry and binary
hashes, selected indices, commands, errors, and R session information. This is
sampled validation, not full-density qualification. The brute-force kernel
still scans all faces for each query; no acceleration or performance claim is
made. Mask, label, missing-data, area-correction and route-admission behavior
are not qualified by these fixtures.
