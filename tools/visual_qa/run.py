"""Build native-reference visual QA. Run from any directory; see README.md."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

import nibabel as nib
import numpy as np
from nibabel.processing import resample_to_output
import SimpleITK as sitk

ROOT = Path(__file__).resolve().parents[2]
IMAGES = {
    "ANTs": "antsx/ants@sha256:59c45f54a1f1dc69134f63bec91a726e41c71c64a16cc21cda0b54526910a3c3",
    "AFNI": "afni/afni_make_build@sha256:6bc2b04e92e0874d7cf252006be35f2671438904c233f844633c40a8fc7cf9bf",
    "FSL": "brainlife/fsl@sha256:fbd262c385e9de22aa58bf7b6311cbd5cd96c7b4eaff151f879191e869bf224e",
}
# Frozen before observing candidate results. Scale is the native reference range.
GATES = {"coordinate_max_mm": .02, "image_nrmse": .0002,
         "image_p99_relative": .001, "finite_fraction": 1.,
         "negative_control_min_mm": .1}


def save(path, data, affine, vector=False):
    im = nib.Nifti1Image(np.asarray(data, dtype=np.float32), affine)
    im.set_qform(affine, code=1)
    im.set_sform(affine, code=1)
    im.header.set_xyzt_units("mm")
    if vector:
        im.header.set_intent("vector")
    nib.save(im, path)


def world_grid(shape, affine):
    ijk = np.stack(np.meshgrid(*[np.arange(n) for n in shape], indexing="ij"), axis=-1)
    return ijk @ affine[:3, :3].T + affine[:3, 3]


def phantom(shape, affine):
    p = world_grid(shape, affine)
    x, y, z = np.moveaxis(p, -1, 0)
    inside = (x/24)**2 + (y/31)**2 + (z/30)**2 < 1
    a = inside * (.25 + .12*np.cos(x/5)*np.cos(y/7)*np.cos(z/6))
    a += .22 * (((x+10)/8)**2 + ((y-6)/17)**2 + (z/19)**2 < 1)
    a += .13 * (((x-11)/9)**2 + ((y+4)/20)**2 + ((z+2)/21)**2 < 1)
    for k, (cx, cy, cz, radius) in enumerate([(-15,-12,-12,4), (13,11,13,5),
                                              (-8,14,4,3), (9,-17,0,3)]):
        a[(x-cx)**2+(y-cy)**2+(z-cz)**2 < radius**2] = .7 + k*.09
    # Three asymmetric rectangular notches remain visible in orthogonal views.
    a[(x < -12) & (abs(y) < 2) & (abs(z) < 15)] = .95
    a[(x > 8) & (abs(y+8) < 2) & (abs(z) < 8)] = .08
    return a


def supporting_images(folder, source):
    im = nib.load(source)
    shape, aff = im.shape[:3], im.affine
    xyz = world_grid(shape, aff)
    for axis in range(3):
        save(folder / f"coord{axis}.nii.gz", xyz[..., axis], aff)
    support = np.zeros(shape)
    support[2:-2, 2:-2, 2:-2] = 1
    save(folder / "support.nii.gz", support, aff)


def native(out, family, case_id, commands, records):
    folder = out / case_id
    script = folder / "native.sh"
    script.write_text("#!/bin/bash\nset -euo pipefail\nexport OMP_NUM_THREADS=1 ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1\ncd /work/" + shlex.quote(case_id) + "\n" +
                      "\n".join(shlex.join(c) for c in commands) + "\n")
    args = ["docker", "run", "--rm", "--network", "none", "--platform", "linux/amd64",
            "--entrypoint", "/bin/bash", "-v", f"{out}:/work", IMAGES[family],
            f"/work/{case_id}/native.sh"]
    print(f"Native {family}: {case_id}", flush=True)
    with (folder / "native.log").open("w") as log:
        result = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT)
    records.append({"case": case_id, "image": IMAGES[family], "commands": commands,
                    "returncode": result.returncode})
    if result.returncode:
        raise RuntimeError(f"Native producer failed; see {folder / 'native.log'}")


def apply_commands(family, transform, kind, extra=None):
    commands = []
    for name in ["source", "coord0", "coord1", "coord2", "support"]:
        inp, output = name + ".nii.gz", "native_" + name + ".nii.gz"
        if family == "ANTs":
            c = ["antsApplyTransforms", "-d", "3", "-i", inp, "-r", "target.nii.gz",
                 "-o", output, "-n", "Linear", "-t", transform]
        elif family == "AFNI":
            if kind == "affine":
                c = ["3dAllineate", "-overwrite", "-input", inp, "-master", "target.nii.gz",
                     "-1Dmatrix_apply", transform, "-final", "linear", "-prefix", output]
            else:
                c = ["3dNwarpApply", "-overwrite", "-source", inp, "-master", "target.nii.gz",
                     "-nwarp", transform, "-interp", "linear", "-ainterp", "linear", "-prefix", output]
        elif kind == "affine":
            c = ["flirt", "-in", inp, "-ref", "target.nii.gz", "-applyxfm", "-init", transform,
                 "-interp", "trilinear", "-out", output, "-paddingsize", "0"]
        else:
            c = ["applywarp", "--in=" + inp, "--ref=target.nii.gz", "--warp=" + transform,
                 "--out=" + output, "--interp=trilinear", "--" + (extra or "rel")]
        commands.append(c)
    return commands


def prepare(out, records, include_anatomy):
    cases = []
    shape = (41, 45, 39)
    cardinal = np.diag([1.5, 1.8, 2.1, 1.])
    cardinal[:3, 3] = -np.array(shape)//2 * np.array([1.5, 1.8, 2.1])
    theta = .17
    rot = np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0], [0,0,1.]])
    oblique = cardinal.copy()
    oblique[:3, :3] = rot @ cardinal[:3, :3]
    oblique[:3, 3] = -oblique[:3, :3] @ ((np.array(shape)-1)/2)

    def new_case(cid, family, title, kind, aff, target_aff=None, scope="synthetic"):
        folder = out / cid
        folder.mkdir(parents=True, exist_ok=True)
        save(folder / "source.nii.gz", phantom(shape, aff), aff)
        save(folder / "target.nii.gz", phantom(shape, target_aff if target_aff is not None else aff),
             target_aff if target_aff is not None else aff)
        case = {"id": cid, "family": family, "title": title, "kind": kind, "scope": scope,
                "source": "source.nii.gz", "target": "target.nii.gz", "notes": []}
        cases.append(case)
        return folder, case

    def affine_lps(folder):
        tx = sitk.AffineTransform(3)
        tx.SetMatrix((1.02, -.12, .02, .11, .98, .01, -.01, .02, 1.01))
        tx.SetCenter((4., -6., 3.))
        tx.SetTranslation((2.4, -1.8, 1.2))
        sitk.WriteTransform(tx, str(folder / "affine.mat"))
        return tx

    for order in ["affine", "affine_warp", "warp_affine"]:
        folder, case = new_case("ants_"+order, "ANTs", {"affine":"Centered affine",
                     "affine_warp":"H5: affine + warp", "warp_affine":"H5: warp + affine"}[order],
                     "affine" if order == "affine" else "h5", oblique)
        tx = affine_lps(folder)
        if order != "affine":
            p = world_grid(shape, oblique)
            d = np.stack([1.8*np.sin(p[...,1]/18)*np.cos(p[...,2]/25),
                          1.4*np.sin(p[...,0]/17), .9*np.cos(p[...,1]/23)], axis=-1)
            save(folder / "warp.nii.gz", (d*[-1,-1,1])[..., None, :], oblique, vector=True)
            warp = sitk.DisplacementFieldTransform(sitk.ReadImage(str(folder / "warp.nii.gz"), sitk.sitkVectorFloat64))
            composite = sitk.CompositeTransform([tx,warp] if order == "affine_warp" else [warp,tx])
            sitk.WriteTransform(composite, str(folder / "composite.h5"))
        case["transform"] = "affine.mat" if order == "affine" else "composite.h5"
        case["notes"] = ["Oblique grid; unequal spacing; nonzero rotation center."]

    for kind in ["affine", "warp", "composite"]:
        folder, case = new_case("afni_"+kind, "AFNI", {"affine":"DICOM affine", "warp":"DICOM displacement",
                    "composite":"DICOM affine + warp"}[kind], kind, cardinal)
        tx = affine_lps(folder)
        A = np.array(tx.GetMatrix()).reshape(3,3)
        center = np.array(tx.GetCenter())
        M = np.column_stack([A, np.array(tx.GetTranslation()) + center - A@center])
        np.savetxt(folder / "affine.aff12.1D", M.reshape(1,12), fmt="%.12g")
        p = world_grid(shape, cardinal)
        d = np.stack([1.8*np.sin(p[...,1]/18)*np.cos(p[...,2]/25),
                      1.4*np.sin(p[...,0]/17), .9*np.cos(p[...,1]/23)], axis=-1)
        save(folder / "warp.nii.gz", d*[-1,-1,1], cardinal)
        case["transform"] = {"affine":"affine.aff12.1D", "warp":"warp.nii.gz",
                              "composite":"affine.aff12.1D warp.nii.gz"}[kind]
        case["notes"] = ["Native 3dAllineate or 3dNwarpApply, explicitly linear interpolation."]

    for hand in ["right", "left"]:
        source_aff = oblique.copy()
        if hand == "left":
            source_aff[:3, 3] += source_aff[:3, 0]*(shape[0]-1)
            source_aff[:3, 0] *= -1
        target_aff = oblique.copy()
        target_aff[:3, 3] += target_aff[:3, 0]*(shape[0]-1)
        target_aff[:3, 0] *= -1
        folder, case = new_case("fsl_affine_"+hand, "FSL", "FLIRT: "+hand+"-handed source",
                                 "affine", source_aff, target_aff)
        M = np.array([[1.01,.05,.01,2.],[-.04,.99,.02,-1.5],[0,-.02,1.01,.8],[0,0,0,1.]])
        np.savetxt(folder / "affine.mat", M, fmt="%.12g")
        case["transform"] = "affine.mat"
        case["notes"] = ["Oblique grids; target left-handed; native FLIRT matrix in scaled voxel coordinates."]

    for representation in ["relative", "absolute"]:
        aff = cardinal.copy()
        aff[0,0] *= -1
        aff[0,3] = -aff[0,3]
        folder, case = new_case("fsl_"+representation, "FSL", "Dense warp: "+representation,
                                 "warp", aff)
        p = world_grid(shape, aff)
        d = np.stack([1.8*np.sin(p[...,1]/18)*np.cos(p[...,2]/25),
                      1.4*np.sin(p[...,0]/17), .9*np.cos(p[...,1]/23)], axis=-1)
        save(folder / "seed.nii.gz", d, aff)
        case["transform"] = "warp.nii.gz"
        case["representation"] = representation
        case["pre_commands"] = [["convertwarp", "--ref=target.nii.gz", "--warp1=seed.nii.gz",
                "--rel", "--"+("relout" if representation == "relative" else "absout"), "--out=warp.nii.gz"]]
        case["notes"] = ["Native convertwarp/applywarp replay of a known smooth field; not FNIRT coefficient or optimizer qualification."]

    if include_anatomy:
        real = [
          ("ants_anatomy", "ANTs", "h5", "chris/chris_t1.nii.gz", "chris/ants/chris_in_mni.nii.gz", "chris/ants/chris_to_mni_Composite.h5"),
          ("ants_anatomy_registered", "ANTs", "h5", "chris/chris_t1.nii.gz", "chris/ants/chris_in_mni.nii.gz", None),
          ("ants_anatomy_affine", "ANTs", "affine", "chris/chris_t1.nii.gz", "chris/ants/chris_in_mni.nii.gz", None),
          ("afni_anatomy", "AFNI", "affine", "afni/ss_sub-1001_T1w.nii.gz", "afni/mni.nii", "afni/sub-1001_T1w_in_mni_Allin.aff12.1D"),
          ("fsl_anatomy", "FSL", "affine", "chris/chris_t1.nii.gz", "chris/ants/chris_in_mni.nii.gz", None)]
        for cid, family, kind, src, ref, transform in real:
            folder, case = new_case(cid, family, "Anatomical replay", kind, cardinal, scope="anatomical")
            inputs = [ROOT/"inst/extdata"/src, ROOT/"inst/extdata"/ref]
            if not all(p.exists() for p in inputs):
                case["missing"] = "Acquired source or target fixture unavailable"
                continue
            for path, name in zip(inputs, ["source.nii.gz", "target.nii.gz"]):
                im = resample_to_output(nib.load(path), voxel_sizes=4., order=1)
                save(folder/name, im.get_fdata(), im.affine)
            case["input_provenance"] = [{"path":str(p.relative_to(ROOT)), "sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs]
            if transform:
                original = ROOT/"inst/extdata"/transform
                if not original.exists():
                    case["missing"] = "Acquired transform fixture unavailable"
                    continue
                shutil.copyfile(original, folder/original.name)
                case["transform"] = original.name
                case["input_provenance"].append({"path":str(original.relative_to(ROOT)), "sha256":hashlib.sha256(original.read_bytes()).hexdigest()})
            elif family == "ANTs":
                case["transform"] = "registration_Composite.h5" if kind == "h5" else "registration_0GenericAffine.mat"
                case["title"] = "Anatomy: affine-only H5" if kind == "h5" else "Anatomy: native affine file"
                case["pre_commands"] = [["antsRegistration", "-d", "3", "--float", "1",
                    "--output", "[registration_,registration.nii.gz]", "--write-composite-transform", "1" if kind == "h5" else "0",
                    "--initial-moving-transform", "[target.nii.gz,source.nii.gz,1]",
                    "--transform", "Affine[0.1]", "--metric", "MI[target.nii.gz,source.nii.gz,1,32,Regular,0.25]",
                    "--convergence", "[100x50x20,1e-6,10]", "--shrink-factors", "4x2x1",
                    "--smoothing-sigmas", "2x1x0vox", "--random-seed", "1729"]]
            else:
                case["transform"] = "affine.mat"
                case["pre_commands"] = [["flirt", "-in", "source.nii.gz", "-ref", "target.nii.gz",
                    "-omat", "affine.mat", "-out", "registration.nii.gz", "-dof", "6", "-cost", "normmi",
                    "-searchrx", "-10", "10", "-searchry", "-10", "10", "-searchrz", "-10", "10"]]
            case["notes"] = ["Acquired T1 anatomy resampled to 4 mm for this replay. Native output regenerated with the exact same inputs and transform.",
                             "This compares transform interpretation, not registration quality."]
            if cid == "ants_anatomy":
                case["title"] = "Legacy anatomy: limited coverage"
                case["notes"].append("Legacy transform pairing is unverified. Only 4.5% of the target has native interior support; numerical agreement does not qualify anatomical alignment.")

    for case in cases:
        if case.get("missing"):
            continue
        folder = out/case["id"]
        supporting_images(folder, folder/"source.nii.gz")
        native(out, case["family"], case["id"], case.get("pre_commands", []) +
               apply_commands(case["family"], case["transform"], case["kind"],
                              "abs" if case.get("representation") == "absolute" else "rel"), records)
    return cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(ROOT/"output/visual-qa"))
    parser.add_argument("--reuse-native", action="store_true")
    parser.add_argument("--synthetic-only", action="store_true")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--strict", action="store_true", help="Exit 2 if any case fails or is unavailable")
    args = parser.parse_args()
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    records = []
    if args.reuse_native or args.report_only:
        manifest = json.loads((out/"manifest.json").read_text())
        if manifest["gates"] != GATES:
            raise RuntimeError("Stored gates differ from frozen gates")
        for name, expected in manifest["native_hashes"].items():
            if hashlib.sha256((out/name).read_bytes()).hexdigest() != expected:
                raise RuntimeError(f"Native input/output changed: {name}")
    else:
        manifest = {"schema": 1, "gates": GATES, "images": IMAGES, "commands": records,
                    "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        manifest["cases"] = prepare(out, records, not args.synthetic_only)
        manifest["native_hashes"] = {str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest()
            for case in manifest["cases"] for p in (out/case["id"]).iterdir()
            if p.is_file() and not p.name.startswith(("ours", "wrong_")) and p.suffix not in (".log", ".sh")}
        (out/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    if not args.report_only:
        env = {**os.environ, "LC_ALL":"en_US.UTF-8", "OMP_NUM_THREADS":"1"}
        with (out/"neurotransform.log").open("w") as log:
            result = subprocess.run(["Rscript", str(ROOT/"tools/visual_qa/compare.R"), str(out)],
                                    cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Candidate evaluation failed; see {out/'neurotransform.log'}")
        manifest["candidate"] = {
            "head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
            "dirty":bool(subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True)),
            "source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                for d in ["R", "src", "tools/visual_qa"] for p in (ROOT/d).rglob("*")
                if p.is_file() and p.suffix in (".R", ".cpp", ".h")},
            "result_sha256":hashlib.sha256((out/"candidate-results.json").read_bytes()).hexdigest(),
            "output_sha256":{str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest()
                for case in manifest["cases"] for p in (out/case["id"]).iterdir()
                if p.is_file() and p.name.startswith(("ours", "wrong_"))}}
        (out/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    for name, expected in manifest["candidate"].get("output_sha256", {}).items():
        if hashlib.sha256((out/name).read_bytes()).hexdigest() != expected:
            raise RuntimeError(f"Candidate output changed: {name}")
    if hashlib.sha256((out/"candidate-results.json").read_bytes()).hexdigest() != manifest["candidate"]["result_sha256"]:
        raise RuntimeError("Candidate results changed")
    manifest["renderer_sha256"] = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (ROOT/"tools/visual_qa").iterdir() if p.suffix in (".py", ".html", ".txt", ".cjs")}
    (out/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    subprocess.run([sys.executable,str(ROOT/"tools/visual_qa/report.py"), str(out)],check=True)
    print(f"Report: {out/'index.html'}", flush=True)
    if args.strict:
        scores=json.loads((out/"summary.json").read_text())
        if scores["failed"] or scores["blocked"]:
            sys.exit(2)


if __name__ == "__main__":
    main()
