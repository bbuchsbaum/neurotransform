# Native transform interpretation QA

Build an offline HTML viewer, a PDF, and per-case PNG panels from actual native
ANTs, AFNI, and FSL resampling commands. The report compares neurotransform's
public coordinate and image-resampling paths with those independent outputs.
Failures are retained in the report.

## Run

Requirements: Python 3.12 or newer, Docker with Linux/amd64 emulation where
needed, R, the package's dependencies, and `pkgload`, `jsonlite`, and `RNifti`.
From the repository root:

```sh
python3 -m venv /private/tmp/neurotransform-qa-venv
/private/tmp/neurotransform-qa-venv/bin/pip install -r tools/visual_qa/requirements.txt
/private/tmp/neurotransform-qa-venv/bin/python tools/visual_qa/run.py --strict
open output/visual-qa/index.html
```

Docker uses immutable image digests declared in `run.py`. Native processes run
without network access. The images may need downloading before the first run.
The tested producers are ANTs 2.6.5.dev1-gfdce4d2, AFNI 26.1.04, and FSL 5.0.9.
These results do not qualify other producer versions.

`--strict` writes the complete report and then exits **2** if any case fails or
is unavailable. Without it, a successfully generated diagnostic report exits 0
even when comparisons fail. Producer or infrastructure errors always stop the
run with an error. Use a separate `--output PATH` to retain multiple runs.

```sh
# Verify native input/output hashes, then reevaluate the current package.
python tools/visual_qa/run.py --reuse-native --strict
# Verify recorded data hashes and rebuild presentation from recorded results.
python tools/visual_qa/run.py --report-only --strict
# Run just the ten convention phantoms.
python tools/visual_qa/run.py --synthetic-only --output /private/tmp/qa-synthetic --strict
# Check that scoring cannot hide invalid samples or coordinate-only failures.
python -m unittest discover -s tools/visual_qa
```

The generated directory is gitignored and excluded from R package builds.
Keep `index.html` and `data.js` together. No server or network connection is
required to view the report. The PDF and PNGs are standalone exports.

## Evidence and gates

Every case uses the same source image, target grid, transform, and linear
interpolation in the native tool and neurotransform. Native tools also resample
three source RAS coordinate ramps. Those ramps independently reveal where each
target voxel samples the source, even where image texture is uninformative.

The phantom is asymmetric, with unequal voxel spacing, distinct landmarks,
nonzero transform centers, oblique grids, both affine handednesses, and
spatially varying fields. H5 cases exercise both noncommuting component orders.
AFNI's `-nwarp` list evaluates left to right; the package's `MorphismPath` list
evaluates right to left, so the harness explicitly reverses that list.

The gates were fixed before inspecting candidate results:

| Quantity | Required value |
| --- | --- |
| Maximum 3D source-coordinate error | ≤ 0.02 mm |
| Image RMSE / native interior intensity range | ≤ 0.0002 |
| Image absolute residual p99 / that range | ≤ 0.001 |
| Finite candidate interior samples | 100% |
| Native-defined interior size | ≥ 1,000 samples |

The interior is native-resampled source support > 0.999, with two voxels removed
at both source and target edges. The candidate cannot select the mask. Invalid
candidate values remain failures; finite-only error summaries never override
that requirement. Full images remain visible, including excluded padding.
Magenta means an invalid value. This is a numerical interpolation interior,
not a brain segmentation mask; background is included.

Deliberate X/Y sign errors, inverse affine direction, and swapped component
order must fail the same gates and produce ≥ 0.1 mm coordinate error or invalid
samples. A case cannot pass if an applicable negative control is ineffective.
These controls probe sensitivity; they are not exhaustive fault coverage.

The viewer has synchronized orthogonal slice controls, checkerboard/overlay
comparison, fixed-scale residual maps, coordinate errors, a deformation
lattice, and seven fixed target landmarks. Negative-control selection updates
the comparison, lattice, and landmark errors. Jacobians are diagnostics of the
correct candidate only, not independently qualified by this report. Display
volumes are quantized for size; all scores use original floating-point data.
Voxel planes retain their recorded orientation, with explicit axis labels;
pixel aspect ratios account for voxel spacing.

## Current findings and limits

The initial native comparison exposes failures in the **relative and absolute
FSL dense-field loaders** on the recorded left-handed grid. Maximum coordinate
errors are approximately 3.60 mm and 101 mm, respectively; the absolute case
also has invalid interior samples. These failures are retained in the pre-fix snapshot under `before-fsl-fix/`.
The corrected candidate uses explicit source/reference geometry and is scored
against the same native files and unchanged gates. An
independent scaled-voxel-to-RAS calculation explains the discrepancy for this
matched-grid example, but is only a diagnostic: it does not replace the scored
package output or establish a general fix for arbitrary source/target grids.

The available FSL producer is 5.0.9. FNIRT coefficient decoding, newer FSL
versions, and optimizer accuracy are outside this report's qualification.
The eleven pre-existing FNIRT tests with missing reference fixtures remain
unresolved; dense applywarp examples do not substitute for them.

Acquired T1 images are replayed at 4 mm. Existing AFNI and ANTs transforms are
retained, and fresh native ANTs affine and FSL rigid registrations provide
known input/transform pairing. The legacy ANTs composite covers only 4.5% of
the target under the native interior rule. Its pairing is unverified and the
report labels it as limited coverage even if the numerical comparison passes.
The fresh affine-only ANTs H5 composite was initially unavailable. Explicit H5
loading now routes it through the ITK affine reader. The same deterministic native
registration is also written as an affine file and evaluated through the
linear-transform API. These are separate cases; one cannot qualify the other.
Visual agreement with a native replay establishes interpretation agreement for
those inputs, not successful anatomical registration or clinical suitability.

`manifest.json` records producer image digests, exact argument arrays, original
anatomical input hashes, native input/output hashes, candidate source hashes,
candidate output hashes, and renderer hashes. Native shell scripts and logs
are retained per case; `candidate-results.json` records R session information.
The candidate is the working tree, including uncommitted changes. Reused data
must match recorded hashes. `--report-only` displays the recorded candidate;
use `--reuse-native` to evaluate subsequent package changes.

## Visual verification

`browser_qa.cjs` uses an installed Playwright headless browser and an ephemeral
owned profile. Follow the repository browser policy: audit before and after,
never attach to other sessions, and do not use system Chrome. It checks every
case, slice changes, all controls, JavaScript errors, and desktop/mobile
overflow, then closes its browser and records ownership and closure.

```sh
node ~/.local/share/agent-policy/browser-automation-guard.mjs --audit
node tools/visual_qa/browser_qa.cjs output/visual-qa
node ~/.local/share/agent-policy/browser-automation-guard.mjs --audit
pdftoppm -png -r 80 output/visual-qa/visual-qa.pdf /private/tmp/qa-page
```

Inspect the rendered PDF pages and browser screenshots after presentation
changes. A passing browser check says nothing about numerical correctness.

Primary convention references:
[ANTs transforms](https://github.com/ANTsX/ANTsPy/wiki/ANTs-transform-concepts-and-file-formats),
[AFNI composition](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dNwarpApply.html),
[AFNI matrices](https://afni.nimh.nih.gov/pub/dist/doc/misc/Registration_Matrix_Save.html),
[FSL FNIRT/applywarp](https://fsl.fmrib.ox.ac.uk/fsl/docs/registration/fnirt/user_guide.html).
