"""Compare both directions and hemispheres against Workbench BARYCENTRIC.

Usage: python3 tools/compare_barycentric_templates.py WB_COMMAND INPUTS OUTPUT [--full]
INPUTS contains fsaverage 164k and fsLR 32k (space-fsaverage) spheres. OUTPUT
must be new. Set R_LIBS to the rebuilt neurotransform library being tested.
Without --full this samples 64 targets; --full compares every target. Neither
mode establishes area/mask qualification.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import platform
import re

import numpy as np
from generate_barycentric_oracle import read_gifti, write_gifti, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("workbench", "inputs", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--full", action="store_true", help="Compare every target, including reordered source faces")
    args = parser.parse_args()
    wb, inputs, out = args.workbench.resolve(), args.inputs.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    # Fixed in advance: float32 coordinates at radius 100 and sub-mm template
    # triangles amplify coordinate rounding in weights; all fields are in [-1,1].
    tolerance, seed, sample_size = 5e-5, 20926, 64
    receipt = dict(method="BARYCENTRIC", tolerance=tolerance, seed=seed, full=args.full,
                   target_sample_size=None if args.full else sample_size,
                   platform=platform.platform(), threads=os.environ.get("OMP_NUM_THREADS"),
                   revision=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
                   version=subprocess.check_output([str(wb), "-version"], text=True),
                   binary_sha256=sha(wb), generator_sha256=sha(__file__), cases=[])
    manifest_path = os.environ.get("NEUROTRANSFORM_BUILD_MANIFEST")
    if manifest_path: receipt["build_manifest"] = json.loads(Path(manifest_path).read_text())
    script = out / "compare.R"
    script.write_text('''
library(neurotransform)
library(jsonlite)
args <- commandArgs(TRUE)
folder <- args[1]
read <- function(name) as.matrix(read.csv(file.path(folder, name), header=FALSE))
vertices <- read("vertices.csv")
faces <- read("faces.csv") + 1L
query <- read("query.csv")
data <- read("data.csv")
expected <- read("expected.csv")
elapsed <- system.time({
  plan <- surface_resampling_plan(surface_mesh(query), surface_mesh(vertices, faces))
})[["elapsed"]]
application_elapsed <- system.time(got <- apply_surface_resampling(plan, data, normalize="none"))[["elapsed"]]
row_mass <- numeric(nrow(query))
mass <- rowsum(plan$vals, plan$rows, reorder=FALSE)
row_mass[as.integer(rownames(mass))] <- mass[,1]
geometry_valid <- all(is.finite(plan$vals)) && all(plan$vals >= 0) &&
  all(plan$rows >= 1 & plan$rows <= nrow(query)) && all(plan$cols >= 1 & plan$cols <= nrow(vertices)) &&
  max(abs(row_mass-1)) < 1e-12
reordered <- surface_resampling_plan(surface_mesh(query), surface_mesh(vertices, faces[nrow(faces):1,,drop=FALSE]))
permutation_error <- max(abs(apply_surface_resampling(reordered,data,normalize="none")-got))
error <- max(abs(got - expected))
worst <- arrayInd(order(abs(got-expected),decreasing=TRUE)[seq_len(min(10,length(got)))],dim(got))
receipt <- list(max_absolute_error=error, elapsed_seconds=elapsed,
                rms_error=sqrt(mean((got-expected)^2)),
                error_quantiles=quantile(abs(got-expected),c(.5,.95,.99,1)),
                worst=worst, worst_query_coordinates=query[worst[,1],,drop=FALSE],
                field_max_errors=apply(abs(got-expected),2,max),
                application_seconds=application_elapsed, native_timing=plan$timing,
                max_row_sum_error=max(abs(row_mass-1)), geometry_valid=geometry_valid,
                permutation_error=permutation_error, targets=nrow(query), nnz=length(plan$vals),
                passed=geometry_valid && permutation_error < 1e-12 && is.finite(error) && error < as.numeric(args[2]),
                dll_sha256=digest::digest(file=getLoadedDLLs()[["neurotransform"]][["path"]], algo="sha256"),
                session=capture.output(sessionInfo()))
write_json(receipt, file.path(folder, "result.json"), auto_unbox=TRUE, pretty=TRUE, digits=NA)
quit(status=if (receipt$passed) 0L else 1L)
''')
    for hemi in ("L", "R"):
        paths = [inputs / f"tpl-fsaverage_hemi-{hemi}_den-164k_sphere.surf.gii",
                 inputs / f"tpl-fsLR_space-fsaverage_hemi-{hemi}_den-32k_sphere.surf.gii"]
        for source, target in (paths, paths[::-1]):
            name = hemi + "-" + ("164k-to-32k" if source == paths[0] else "32k-to-164k")
            folder = out / name
            folder.mkdir()
            vertices, faces = read_gifti(source)
            query = read_gifti(target)[0]
            p = vertices / np.linalg.norm(vertices, axis=1)[:, None]
            data = np.column_stack([p, p[:, 0] * p[:, 1],
                                    np.sin(np.arange(len(p)) * .017), np.cos(np.arange(len(p)) * .031)]).astype(np.float32)
            rng=np.random.default_rng(seed)
            impulses=np.zeros((len(p),4),dtype=np.float32)
            impulses[rng.choice(len(p),4,replace=False),np.arange(4)]=1
            data=np.column_stack([data,rng.uniform(-1,1,len(p)).astype(np.float32),impulses])
            write_gifti(folder / "input.func.gii", [("NIFTI_INTENT_SHAPE", col) for col in data.T])
            command = [str(wb), "-metric-resample", str(folder / "input.func.gii"), str(source),
                       str(target), "BARYCENTRIC", str(folder / "output.func.gii")]
            run = subprocess.run(command, capture_output=True, text=True, check=True)
            expected = np.column_stack(read_gifti(folder / "output.func.gii"))
            indices = np.arange(len(query)) if args.full else np.sort(np.random.default_rng(seed).choice(len(query), sample_size, replace=False))
            for filename, values in [("vertices", vertices), ("faces", faces), ("data", data),
                                      ("query", query[indices]), ("expected", expected[indices])]:
                np.savetxt(folder / (filename + ".csv"), values, delimiter=",", fmt="%.17g")
            comparison_argv=["/usr/bin/time", "-l" if platform.system()=="Darwin" else "-v",
                             "Rscript", str(script), str(folder), str(tolerance)]
            compare = subprocess.run(comparison_argv,
                                     capture_output=True, text=True, env=os.environ.copy())
            (folder / "comparison.log").write_text(compare.stdout + compare.stderr)
            result = json.loads((folder / "result.json").read_text()) if (folder / "result.json").exists() else {}
            peak=re.search(r"(\d+)\s+maximum resident set size",compare.stderr)
            if peak: result["peak_resident_bytes"]=int(peak.group(1))
            peak_linux=re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)",compare.stderr)
            if peak_linux: result["peak_resident_bytes"]=1024*int(peak_linux.group(1))
            receipt["cases"].append(dict(name=name, source=str(source), target=str(target),
                source_sha256=sha(source), target_sha256=sha(target), target_vertices=len(query),
                selected_target_indices_zero_based="all" if args.full else indices.tolist(), argv=command,
                comparison_argv=comparison_argv,
                workbench_exit=run.returncode, workbench_stderr=run.stderr,
                comparison_exit=compare.returncode, result=result))
            (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
            print(name, result.get("max_absolute_error"), "exit", compare.returncode, flush=True)
    if any(case["comparison_exit"] != 0 for case in receipt["cases"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
