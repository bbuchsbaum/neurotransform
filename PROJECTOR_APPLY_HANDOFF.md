# Handoff: `cpp_apply_projector` is broken on real projector shapes

**Audience:** neurotransform engineers
**Status:** neurotransform HEAD (`3c68a8b`, the "Remove Grid@domain architectural leak" commit) ships a `cpp_apply_projector` that **hard-crashes** on the shapes of projector a real user would hand it. This is a bug, not an optimization question. The fix is drop-in and already validated. This document has everything you need to take it from here without talking to anyone.

**Why you are hearing about this now:** during a deduplication pass between neurotransform and neurofunctor, we ran a benchmark of neurotransform's `cpp_apply_projector` against neurofunctor's independently-developed specialized version (`rcpp_projector_apply_direct.cpp`). neurotransform's version failed with `SpMat::init(): requested size is too large` on every realistic shape we tried. neurofunctor's version handled them all and was 1.66×–5.61× faster on shapes where both ran. Details below.

---

## 1. TL;DR

- **Bug:** `neurotransform::cpp_apply_projector` calls `arma::sp_mat(locations, values, n_rows, n_cols, ...)` where `locations` is an `arma::umat`. On the default 32-bit Armadillo build (no `ARMA_64BIT_WORD`), `arma::sp_mat::init()` validates that `n_rows * n_cols` fits in a 32-bit `uword`. For a 300k-voxel target × 32k-vertex source projector — **the standard fsLR-32k → MNI-2mm resampler** — the product is ≈ 9.6 × 10⁹, which overflows and throws. The function is unusable for production neuroimaging shapes.
- **Fix:** Replace the function body with a direct-`dgCMatrix`-slot SpMV that never constructs an `arma::sp_mat`. ~60 lines of C++. Full source below. This implementation already exists, is already tested, and is already live in neurofunctor as `src/rcpp_projector_apply_direct.cpp`.
- **Bonus:** After porting, neurofunctor can retire its duplicate copy and route all `cpp_apply_projector` calls through `neurotransform:::cpp_apply_projector` via its existing `R/cpp_shims.R` pattern, eliminating the last remaining `.cpp` file unique to neurofunctor and closing out the cross-package C++ deduplication.
- **Effort:** ~30 minutes of engineering work plus one benchmark run. No algorithm design, no new dependencies.

---

## 2. Reproducer

From a checkout of neurotransform and a checkout of neurofunctor side-by-side, run this one-liner. It builds a 300k × 32k sparse projector (4 nnz/row) and calls the current `cpp_apply_projector`:

```r
library(Matrix)
set.seed(1)
n_rows <- 300000L
n_cols <- 32000L
nnz_per_row <- 4L
i <- rep(seq_len(n_rows), each = nnz_per_row)
j <- sample.int(n_cols, n_rows * nnz_per_row, replace = TRUE)
x <- runif(length(i))
P <- sparseMatrix(i = i, j = j, x = x, dims = c(n_rows, n_cols), repr = "C")
D <- matrix(runif(n_cols * 50), n_cols, 50)

neurotransform:::cpp_apply_projector(P, D, threads = 1L)
# Error: SpMat::init(): requested size is too large;
#        suggest to enable ARMA_64BIT_WORD
```

This is not a synthetic edge case. This is "surface-to-volume resampling of a standard subject," the very thing a transform kernel package exists to do. A dense `50M voxels × 200 timepoint` fMRI workload is within an order of magnitude of this, and also crashes.

For completeness, the same projector + data combination running through neurofunctor's specialized version:

```r
neurofunctor:::cpp_apply_projector(P, D, threads = 1L)
# Returns a 300000 x 50 numeric matrix in ~2 seconds. No crash.
```

---

## 3. Benchmark evidence

Benchmark script: `neurofunctor/tools/bench_apply_projector.R`
Raw output preserved at: `neurofunctor/tools/bench_results/apply_projector_baseline.txt`

Six shapes were run single-threaded and with 4 threads:

| label            | shape                       | nf_ms (direct) | nt_ms (arma) | nt/nf |
|------------------|-----------------------------|---------------:|-------------:|------:|
| small-balanced   | 50k × 50k,   4 nnz, T=50    | 112            | 115          | 1.03× |
| fsLR-vol-small   | 100k × 32k,  8 nnz, T=100   | 628            | 1678         | 2.67× |
| fsLR-vol-medium  | 300k × 32k,  8 nnz, T=100   | 1 run          | **crash**    |   —   |
| vol-fsLR-medium  | 32k × 300k,  6 nnz, T=100   | 1 run          | **crash**    |   —   |
| tall-extreme     | 1M × 32k,    8 nnz, T=50    | 1 run          | **crash**    |   —   |
| wide-extreme     | 32k × 1M,    6 nnz, T=50    | 1 run          | **crash**    |   —   |

Single-thread: 1.66× geomean speedup (direct vs arma) on shapes where arma runs at all.
Four-thread:   5.61× geomean speedup — because arma's `sp_mat * mat` does not parallelize over the right-hand-side column dimension and the direct-slot loop does.

On every shape where both ran, outputs matched the arma version **to 1e-10 tolerance**. Correctness is not in dispute.

**Four out of six realistic shapes crash neurotransform.** That is the headline.

---

## 4. Root cause

Current `neurotransform/src/rcpp_apply.cpp`:

```cpp
arma::sp_mat as_spmat(const Rcpp::S4& mat) {
  Rcpp::IntegerVector dims = mat.slot("Dim");
  arma::uword n_rows = dims[0];
  arma::uword n_cols = dims[1];
  ...
  arma::umat locations(2, nnz);
  arma::vec  values(x.begin(), nnz, false, true);
  ...
  return arma::sp_mat(locations, values, n_rows, n_cols, true, false);
}

Rcpp::NumericMatrix cpp_apply_projector(const Rcpp::S4& projMat,
                                        const Rcpp::NumericMatrix& data,
                                        int threads = 1) {
  arma::sp_mat P = as_spmat(projMat);
  arma::mat    X = Rcpp::as<arma::mat>(data);     // full input copy
  arma::mat    Y = P * X;                         // full output allocation
  Rcpp::NumericMatrix out(Y.n_rows, Y.n_cols);
  std::copy(Y.begin(), Y.end(), out.begin());     // full output copy
  return out;
}
```

Three independent problems, in order of severity:

1. **The crash.** The internal `arma::SpMat::init()` validator computes `n_rows * n_cols` as an `arma::uword`. On a default Armadillo build `arma::uword` is `unsigned int` (32-bit), which caps the element count at ≈ 4.3 × 10⁹. Our bread-and-butter shape of 300k × 32k yields 9.6 × 10⁹ and overflows. The validator throws. This happens *regardless* of `nnz` — it's a dense-element-count check, not a storage check. A 5-nnz-per-row projector with only 1.5 M nonzeros still dies because `n_rows * n_cols` is too big for the validator.

2. **The copies.** The function copies the input data (`Rcpp::as<arma::mat>(data)`), allocates a full output (`arma::mat Y`), then copies the output back into an Rcpp matrix (`std::copy`). On a realistic 1.2 GB fMRI volume, that's 3.6 GB of transient memory traffic per call. This is why the function also happens to be slow on shapes where it doesn't crash.

3. **Serial SpMV.** `arma::sp_mat * arma::mat` does not parallelize across the columns of the dense right-hand side, even with OpenMP available. Setting `threads > 1` on an fMRI-shaped workload does nothing. This is why the 4-thread benchmark shows a 5.61× gap — the direct loop parallelizes across the `T` axis and the arma version doesn't.

---

## 5. The fix: port neurofunctor's direct-`dgCMatrix` SpMV

Do not "simplify" to `arma::sp_mat`. The benchmark is the proof.

Drop-in replacement for the body of `neurotransform/src/rcpp_apply.cpp`. This is byte-identical (minus comments and whitespace) to `neurofunctor/src/rcpp_projector_apply_direct.cpp`, which is currently in production in neurofunctor and passing 358 tests:

```cpp
// =============================================================================
// cpp_apply_projector — direct dgCMatrix SpMV
// =============================================================================
//
// This is a deliberate specialization that bypasses arma::sp_mat. Do not
// "simplify" to arma::sp_mat * arma::mat.
//
// RATIONALE
// ---------
// The previous arma::sp_mat implementation had two fatal problems at
// realistic neuroimaging scales:
//
// 1. ARMA 32-BIT OVERFLOW.
//    Without ARMA_64BIT_WORD (the default build), arma::sp_mat stores
//    element counts in 32-bit uwords. arma::SpMat::init() validates that
//    n_rows * n_cols fits in a uword. A surface-to-volume projector of the
//    bread-and-butter shape 300k target voxels x 32k source vertices
//    produces n_rows * n_cols ~= 9.6e9, which overflows and hard-crashes
//    with "SpMat::init(): requested size is too large; suggest to enable
//    ARMA_64BIT_WORD". All four of the following shapes crash the
//    arma::sp_mat path:
//      - 300k x 32k  (fsLR-32k -> MNI-2mm)
//      - 32k x 300k  (MNI-2mm  -> fsLR-32k)
//      - 1M  x 32k   (extreme-tall)
//      - 32k x 1M    (extreme-wide)
//
// 2. COPY OVERHEAD AND SERIAL SpMV.
//    arma::sp_mat * arma::mat allocates a full dense output (n_target x T)
//    and copies input and output once each. arma's SpMV does not
//    parallelize over the time axis T. The direct-slot path below reads
//    dgCMatrix slots in-place, writes the output Rcpp::NumericMatrix
//    directly, and parallelizes the outer loop over T.
//
// CORRECTNESS
// -----------
// Outputs match arma::sp_mat * arma::mat to 1e-10 on all shapes where the
// arma path does not crash. See tools/bench_apply_projector.R (copied from
// neurofunctor) for the reproducer and measured numbers.
//
// WHEN TO RECONSIDER
// ------------------
// This specialization becomes unnecessary if BOTH of the following hold:
//   (a) ARMA_64BIT_WORD is enabled in the Armadillo build, AND
//   (b) a benchmark reruns on the same shapes and shows arma::sp_mat is at
//       least competitive across single- and multi-threaded runs.
// Until then, do not change this to arma::sp_mat.
// =============================================================================

#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_apply_projector(const Rcpp::S4& projMat,
                                        const Rcpp::NumericMatrix& data,
                                        int threads = 1) {
  // Read dgCMatrix slots directly; no arma::sp_mat conversion.
  if (threads < 1) threads = 1;

  Rcpp::IntegerVector dims = projMat.slot("Dim");
  const int n_rows = dims[0];
  const int n_cols = dims[1];

  if (data.nrow() != n_cols) {
    Rcpp::stop("Data has %d rows but projector expects %d source elements",
               data.nrow(), n_cols);
  }

  Rcpp::IntegerVector p = projMat.slot("p");    // column pointers (length n_cols + 1)
  Rcpp::IntegerVector i = projMat.slot("i");    // 0-based row indices (length nnz)
  Rcpp::NumericVector x = projMat.slot("x");    // values (length nnz)

  if (p.size() != n_cols + 1) {
    Rcpp::stop("Invalid dgCMatrix: p slot has length %d but expected %d",
               p.size(), n_cols + 1);
  }
  if (i.size() != x.size()) {
    Rcpp::stop("Invalid dgCMatrix: i slot has length %d but x slot has length %d",
               i.size(), x.size());
  }

  const int n_time = data.ncol();
  Rcpp::NumericMatrix out(n_rows, n_time); // initialized to zeros

  const double* X = data.begin();
  double* Y = out.begin();

#ifdef _OPENMP
  if (threads > 1) {
    omp_set_num_threads(threads);
    #pragma omp parallel for
    for (int c = 0; c < n_time; ++c) {
      for (int col = 0; col < n_cols; ++col) {
        const int k0 = p[col];
        const int k1 = p[col + 1];
        if (k0 == k1) continue;
        const double xcol = X[col + c * n_cols];
        if (xcol == 0.0) continue;
        for (int k = k0; k < k1; ++k) {
          const int row = i[k];
          Y[row + c * n_rows] += x[k] * xcol;
        }
      }
    }
  } else {
    for (int c = 0; c < n_time; ++c) {
      for (int col = 0; col < n_cols; ++col) {
        const int k0 = p[col];
        const int k1 = p[col + 1];
        if (k0 == k1) continue;
        const double xcol = X[col + c * n_cols];
        if (xcol == 0.0) continue;
        for (int k = k0; k < k1; ++k) {
          const int row = i[k];
          Y[row + c * n_rows] += x[k] * xcol;
        }
      }
    }
  }
#else
  for (int c = 0; c < n_time; ++c) {
    for (int col = 0; col < n_cols; ++col) {
      const int k0 = p[col];
      const int k1 = p[col + 1];
      if (k0 == k1) continue;
      const double xcol = X[col + c * n_cols];
      if (xcol == 0.0) continue;
      for (int k = k0; k < k1; ++k) {
        const int row = i[k];
        Y[row + c * n_rows] += x[k] * xcol;
      }
    }
  }
#endif
  return out;
}
```

### Notes on the replacement

- It is pure Rcpp. It still depends on RcppArmadillo (left via the `[[Rcpp::depends(RcppArmadillo)]]` attribute) because the rest of the file — specifically `cpp_apply_affine_chain` and `as_spmat` — uses arma. If you are removing `as_spmat` because nothing else in the file uses it after this change, go ahead and remove it. `cpp_apply_affine_chain` stays: it is still called from `neurotransform/R/warp_chain.R:89,99`.
- The `n_rows` and `n_cols` are `int`, matching R's integer size. R's integers are 32-bit but the arithmetic `row + c * n_rows` can overflow at ~2B elements in the output matrix. For pure safety margin on extreme outputs you can widen these to `std::size_t`; we left them as `int` in neurofunctor because an output bigger than 2B doubles is already 16 GB and beyond what the test shapes touch, but it is an easy follow-up.
- Nothing in the function constructs a dense intermediate matrix. Peak extra memory is the `Rcpp::NumericMatrix out(n_rows, n_time)` allocation — exactly one output-sized buffer, initialized to zero.
- Correctness against the old arma version is independently verified by `tools/bench_apply_projector.R` in neurofunctor (1e-10 tolerance on every shape where the arma version does not crash).

---

## 6. Execution plan

### Phase A — Fix neurotransform (this is the whole bug fix)

1. Replace the `cpp_apply_projector` function body in `neurotransform/src/rcpp_apply.cpp` with the source block from §5.
2. Keep the header comment block. Future readers need to understand why this is not `arma::sp_mat`.
3. Regenerate exports:
   ```r
   Rcpp::compileAttributes("path/to/neurotransform")
   devtools::document("path/to/neurotransform")
   ```
4. Rebuild and run the existing test suite:
   ```r
   devtools::test("path/to/neurotransform")
   ```
   Expected result: **497 tests / 0 failed / 11 skipped**, identical to current HEAD.
5. Copy `neurofunctor/tools/bench_apply_projector.R` into `neurotransform/tools/` (adjust the two `neurofunctor:::` / `neurotransform:::` references — you only need the `nt_apply` side now; drop the direct comparison, or keep it pointed at neurofunctor until you retire that copy in Phase B). Run it and confirm:
   - all shapes complete without error (including the four that used to crash)
   - correctness check still passes on the small shapes
   - performance is at least as good as neurofunctor's numbers in `tools/bench_results/apply_projector_baseline.txt`
6. Commit. Suggested message:
   ```
   Replace arma::sp_mat SpMV with direct dgCMatrix implementation

   The previous cpp_apply_projector crashed with
   "SpMat::init(): requested size is too large" on projector shapes
   as small as 300k x 32k (the standard fsLR-32k -> MNI-2mm resampler)
   because arma::sp_mat's internal size validator uses a 32-bit uword
   on default Armadillo builds.

   The new implementation reads dgCMatrix slots directly, writes the
   output Rcpp::NumericMatrix in place, and parallelizes over the
   time-axis columns. Measured 1.66x single-thread and 5.61x four-thread
   speedup vs the old arma path (geomean), and correct output to 1e-10
   on every shape where the old path runs.

   Full rationale and benchmark evidence: PROJECTOR_APPLY_HANDOFF.md.
   ```

### Phase B — Retire the duplicate in neurofunctor (optional but recommended)

Once Phase A is merged, neurofunctor's local `src/rcpp_projector_apply_direct.cpp` is redundant — its sole reason to exist was "neurotransform's version is broken." That reason is now gone. Neurofunctor can route `cpp_apply_projector` through `neurotransform:::cpp_apply_projector` using the same shim pattern it already uses for the other 5 cross-package cpp kernels.

Steps in the neurofunctor repo:

1. `git rm src/rcpp_projector_apply_direct.cpp`
2. In `R/cpp_shims.R`, add `cpp_apply_projector <- neurotransform:::cpp_apply_projector` alongside the existing five shims and remove the "explicitly NOT shimmed" comment about `cpp_apply_projector`.
3. `Rcpp::compileAttributes()` + `devtools::document()`. This will drop the `_neurofunctor_cpp_apply_projector` symbol from `src/RcppExports.cpp` and `R/RcppExports.R`. The neurofunctor `.so` is now effectively empty (only Rcpp's infrastructure stub remains).
4. `devtools::test()` — expect **358 tests / 0 failed / 5 skipped**, identical to current HEAD.
5. Verify call-site resolution:
   ```r
   f <- get("cpp_apply_projector", envir = asNamespace("neurofunctor"))
   environmentName(environment(f))
   # expect "neurotransform"
   ```
6. Commit in neurofunctor. Suggested message:
   ```
   Retire rcpp_projector_apply_direct.cpp; route via neurotransform

   neurotransform now ships the direct-dgCMatrix cpp_apply_projector
   (see neurotransform commit <SHA>), so neurofunctor no longer needs
   its own specialized copy. The R-level call sites in R/projector.R
   resolve unchanged via R/cpp_shims.R.
   ```

After Phase B, neurofunctor has **zero** `.cpp` files of its own (aside from the autogenerated `RcppExports.cpp`). All transform math lives in neurotransform. The package boundary is clean.

---

## 7. Verification protocol

Two checks. Both should pass before Phase A is considered done.

### Check 1 — Correctness

Run the neurofunctor benchmark script (it cross-checks both versions at 1e-10 tolerance on every shape that doesn't crash). After Phase A, the "crashed" entries become timings and the correctness column should be `TRUE` everywhere.

### Check 2 — Performance parity

The new neurotransform implementation should match neurofunctor's published numbers for the shapes that work today and should simply not crash on the other four. Acceptable regression budget: ±5% on the two shapes that currently run. Anything worse suggests a copy-paste error.

### Check 3 (optional but smart) — No memory regression

Run each benchmark shape through `gc()` before and after. The direct-slot implementation should allocate at most one output-sized buffer per call. If you see multiple GB of transient allocation, something was copied that didn't need to be.

---

## 8. Why not just enable `ARMA_64BIT_WORD`?

This is the obvious "smaller" fix. It is not the right fix.

- **It only addresses the crash.** You still pay the full `as_spmat()` conversion cost, the `Rcpp::as<arma::mat>(data)` input copy, the dense `arma::mat Y` output allocation, and the output-copy-back. The benchmark measured those costs on shapes where the arma path runs: 1.66× single-thread, 5.61× 4-thread disadvantage. Those costs do not depend on 32- vs 64-bit indices; they are structural to the arma path.
- **`ARMA_64BIT_WORD` is a compile-time flag you don't control.** You can set it yourself via `-DARMA_64BIT_WORD`, but downstream packages that also include `RcppArmadillo.h` get whatever their own build sets, and if those builds don't match at the header level you get ODR problems. Setting the flag package-wide across the neurotransform + neurofunctor + future-consumer dependency graph is a much harder ask than porting a 60-line function.
- **Doubling the index width doubles SpMat storage.** For a 50M-nonzero projector you are paying 400 MB of index memory instead of 200 MB. The direct-slot path doesn't store any extra indices at all — it reads from R's `dgCMatrix` slots in place.
- **Serial SpMV is still serial.** Even with 64-bit words, `arma::sp_mat * arma::mat` does not parallelize across the RHS column dimension. The multi-threaded gap stays 5×+.
- **You still pay the `Rcpp::as<arma::mat>(data)` input copy.** For a 1.2 GB fMRI volume that is a full 1.2 GB of avoidable memory traffic per call.

Enable `ARMA_64BIT_WORD` if you want, but do it in addition to this fix, not in place of it. And re-run the benchmark before you decide the arma path is good enough again — the header comment in the replacement spells out that exact reevaluation criterion.

---

## 9. Questions you may reasonably ask

**Q: Is the direct-slot loop correct for all `dgCMatrix` variants?**
Yes within the `dgCMatrix` contract: column-compressed, 0-based row indices in the `i` slot, column pointers in `p` with `length(p) == n_cols + 1`, and `length(i) == length(x) == nnz`. The function validates the slot invariants at the top and `Rcpp::stop`s on violation. Logical sparse variants (`lgCMatrix`), symmetric (`dsCMatrix`), and triangular (`dtCMatrix`) are not `dgCMatrix` and should be coerced at the caller before this function sees them — same as before.

**Q: What about numerical drift vs arma?**
The direct loop is a straight textbook `y[i] += A[i,j] * x[j]` accumulation in double. Its order-of-summation is deterministic (column-major outer loop, nnz inner loop). The measured tolerance vs arma is 1e-10 across all shapes in the benchmark. Any drift is in the last few bits.

**Q: Is the OpenMP `parallel for` over `c` safe?**
Yes. Each thread writes to a disjoint column slice of `out` (`Y[row + c * n_rows]`) with fixed `c`. There is no overlap between iterations of the `c` loop. No atomics, no locks, no races.

**Q: Does this need `Rcpp::NumericMatrix` to be zero-initialized?**
Yes, and it is. Rcpp's `NumericMatrix(n_rows, n_time)` zero-initializes by default. The accumulator pattern `Y[...] += ...` depends on that.

**Q: Should I keep the RcppArmadillo dependency?**
As long as `cpp_apply_affine_chain` (and anything else in the same compilation unit or elsewhere in `src/`) uses `arma::mat`, yes. Removing it is a separate cleanup, not blocked by this fix.

**Q: Is there any risk to downstream packages that currently call `neurotransform::cpp_apply_projector` and pass unusual inputs?**
The public signature is unchanged: `(const Rcpp::S4& projMat, const Rcpp::NumericMatrix& data, int threads)`. Anything that worked before still works. Anything that was crashing before now works. There is no semantic change for callers.

---

## 10. File references

**In neurofunctor (for context and the source of the replacement):**

- `neurofunctor/src/rcpp_projector_apply_direct.cpp` — the production-tested version you are porting. Header comment spells out the rationale.
- `neurofunctor/tools/bench_apply_projector.R` — the benchmark. Self-contained, uses `neurofunctor:::cpp_apply_projector` and `neurotransform:::cpp_apply_projector` side by side.
- `neurofunctor/tools/bench_results/apply_projector_baseline.txt` — raw run output from 2026-04-12. The measured numbers in this document come from here.
- `neurofunctor/R/cpp_shims.R` — the pattern for routing cross-package C++ kernels. Phase B will add `cpp_apply_projector` to this file.
- `neurofunctor/planning_docs/cpp-dedup-investigation.md` — the original investigation plan that surfaced this issue. Broader than just this kernel but §5 specifically discusses `rcpp_apply.cpp`.

**In neurotransform (where you will be working):**

- `neurotransform/src/rcpp_apply.cpp` — the file you are editing.
- `neurotransform/R/RcppExports.R:13-15` — the R wrapper that dispatches `.Call("_neurotransform_cpp_apply_projector", ...)`. No change needed; `compileAttributes()` regenerates it untouched.
- `neurotransform/R/warp_chain.R:89,99` — the live callers of `cpp_apply_affine_chain` in the same compilation unit. Leave them alone, they're fine.
- `neurotransform/tests/testthat/` — all 497 tests should keep passing. There are no tests exercising extreme projector shapes today, so the existing suite is unaffected. Consider adding one. A minimal version could reuse the `make_projector()` helper from the benchmark script with a shape like 100k × 32k + T=20.

**Relevant commits (the audit trail for why this was found):**

In neurofunctor:
- `a0eccdd` — "Deduplicate R-level overlap with neurotransform"
- `72a7d73` — "Delete dead C++ kernels and document specialized SpMV rationale" (this is where the header rationale was first written)
- `3b40641` — "Delete 5 byte-identical C++ kernels; route to neurotransform via shims"

In neurotransform:
- `3c68a8b` — "Remove Grid@domain architectural leak" (current HEAD at the time of this report)

---

## 11. One thing you should know that isn't in the code

The deduplication pass that surfaced this did not start as a performance investigation. It started as "both packages have byte-identical `rcpp_bary.cpp`, `rcpp_ribbon.cpp`, etc. — let's delete neurofunctor's copies and route through neurotransform." When we got to `rcpp_apply.cpp` we noticed it was the only *divergent* file in the set, dug in, and the crash appeared on the first realistic shape we tried. We were not hunting for a bug. The bug fell out of routine housekeeping.

Which means **this bug is almost certainly reachable by anyone who tries to use `neurotransform::cpp_apply_projector` on their first real projector.** Fix it first, then ship.
