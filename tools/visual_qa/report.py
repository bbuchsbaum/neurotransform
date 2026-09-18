"""Score immutable native references and render an offline HTML/PNG/PDF report."""
from pathlib import Path
import base64
import json
import os
import shutil
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "neurotransform-qa-matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import nibabel as nib
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.utils import ImageReader

HERE = Path(__file__).resolve().parent
COLORS = {"ink":"#122a3c", "muted":"#506677", "teal":"#087c82", "red":"#b74239", "amber":"#996419"}


def load(folder, name):
    return nib.load(folder/name).get_fdata(dtype=np.float32)


def finite_number(x):
    return float(x) if np.isfinite(x) else None


def json_safe(x):
    if isinstance(x, dict):
        return {k:json_safe(v) for k,v in x.items()}
    if isinstance(x, (list, tuple)):
        return [json_safe(v) for v in x]
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def score(native, ours, native_coords, ours_coords, mask, gates):
    count = int(mask.sum())
    finite = np.isfinite(ours) & np.isfinite(ours_coords).all(axis=-1)
    scale = float(np.ptp(native[mask])) if count else 0.
    all_valid = mask & finite
    error = np.linalg.norm(ours_coords-native_coords, axis=-1)
    residual = np.abs(ours-native)
    fraction = float(finite[mask].mean()) if count else 0.
    # Finite-only diagnostics do not hide invalid points: finite_fraction is
    # independently required to equal 1, and its denominator is the fixed mask.
    maximum = float(error[all_valid].max()) if all_valid.any() else np.inf
    rms_mm = float(np.sqrt(np.mean(error[all_valid]**2))) if all_valid.any() else np.inf
    nrmse = float(np.sqrt(np.mean(residual[all_valid]**2))/scale) if all_valid.any() and scale else np.inf
    p99 = float(np.quantile(residual[all_valid], .99)/scale) if all_valid.any() and scale else np.inf
    passed = (count >= 1000 and fraction == gates["finite_fraction"] and
              maximum <= gates["coordinate_max_mm"] and nrmse <= gates["image_nrmse"] and
              p99 <= gates["image_p99_relative"])
    return {"pass": bool(passed), "samples": count, "mask_fraction": float(mask.mean()),
            "finite_fraction": fraction, "invalid_samples": int((mask & ~finite).sum()),
            "coordinate_max_mm": finite_number(maximum), "coordinate_rms_mm": finite_number(rms_mm),
            "image_nrmse": finite_number(nrmse), "image_p99_relative": finite_number(p99),
            "intensity_range": scale}, error, residual


def pack(a, low, high):
    good = np.isfinite(a)
    v = np.full(a.shape, 65535, dtype="<u2")
    v[good] = np.round(np.clip((a[good]-low)/(high-low), 0, 1)*65534).astype(np.uint16)
    return {"low":float(low), "high":float(high),
            "bytes":base64.b64encode(v.ravel(order="F").tobytes()).decode()}


def slice2(a, axis, k):
    return np.rot90(np.take(a, k, axis=axis))


def num(value, digits=3):
    if value is None:
        return "unavailable"
    return f"{value:.{digits}g}"


def panel(path, case, volumes, metrics, native_coords, ours_coords, mask, window, spacing, axis_codes):
    fig = plt.figure(figsize=(14, 9), facecolor="#f4f7f8")
    gs = fig.add_gridspec(4, 4, height_ratios=[1,1,1,.72],
                          left=.045, right=.96, bottom=.055, top=.87, hspace=.27, wspace=.12)
    state = case["status"].upper()
    state_color = COLORS["teal"] if case["status"] == "pass" else COLORS["red"]
    fig.text(.045,.955,f"{case['family']} / {case['title']}", fontsize=20, weight="bold", color=COLORS["ink"])
    fig.text(.045,.919,case["scope"].upper()+"  ·  linear interpolation  ·  native tool reference",fontsize=10,color=COLORS["muted"])
    fig.text(.955,.947,state, ha="right",fontsize=18,weight="bold",color=state_color)
    headings = ["Native reference", "neurotransform", "Checkerboard", "Absolute difference"]
    vmax = window[1]
    errorscale = max(metrics["intensity_range"]*.01, 1e-6)
    cmap = plt.get_cmap("magma").copy()
    cmap.set_bad("#ee46d0")
    grayscale = plt.get_cmap("gray").copy()
    grayscale.set_bad("#ee46d0")
    for row, (axis,label) in enumerate([(2,"Axial"),(1,"Coronal"),(0,"Sagittal")]):
        k = volumes["native"].shape[axis]//2
        a,b = slice2(volumes["native"],axis,k), slice2(volumes["ours"],axis,k)
        yy,xx=np.indices(a.shape)
        board=(xx//5+yy//5)%2
        images=[a,b,np.where(board,a,b),slice2(volumes["error"],axis,k)]
        for col,img in enumerate(images):
            ax=fig.add_subplot(gs[row,col])
            displayed_axes = [i for i in range(3) if i != axis]
            ax.imshow(img, cmap=cmap if col==3 else grayscale, vmin=0,
                      vmax=errorscale if col==3 else vmax, interpolation="nearest",
                      aspect=spacing[displayed_axes[1]]/spacing[displayed_axes[0]])
            ax.set_xticks([]); ax.set_yticks([])
            opposite = dict(R="L", L="R", A="P", P="A", S="I", I="S")
            horizontal, vertical = [axis_codes[i] for i in displayed_axes]
            for x,y,label_code in [(.03,.5,opposite[horizontal]),(.97,.5,horizontal),
                                   (.5,.97,vertical),(.5,.03,opposite[vertical])]:
                ax.text(x,y,label_code,transform=ax.transAxes,ha="center",va="center",
                        fontsize=7,color="#7d9fad")
            for spine in ax.spines.values(): spine.set_visible(False)
            if row==0: ax.set_title(headings[col],fontsize=11,color=COLORS["ink"],pad=8)
            if col==0: ax.set_ylabel(f"{label}\nindex {k}",fontsize=9,color=COLORS["muted"])
    ax=fig.add_subplot(gs[3,:2]); ax.axis("off")
    text=(f"Max coordinate error   {num(metrics['coordinate_max_mm'])} mm\n"
          f"RMS coordinate error   {num(metrics['coordinate_rms_mm'])} mm\n"
          f"Image NRMSE   {num(metrics['image_nrmse'])}     /     gate 0.0002\n"
          f"Native-defined interior   {metrics['mask_fraction']:.1%} of target\n"
          f"Invalid candidate samples   {metrics['invalid_samples']:,}")
    ax.text(0,.95,text,va="top",fontsize=10,linespacing=1.6,color=COLORS["ink"])
    ax=fig.add_subplot(gs[3,2]); ax.axis("off")
    variants=case.get("variant_metrics",[])
    lines=["DELIBERATE ERRORS"]
    for v in variants:
        lines += [v["label"],f"{num(v['metrics']['coordinate_max_mm'])} mm  /  {'REJECTED' if v['rejected'] else 'NOT REJECTED'}"]
    ax.text(0,.95,"\n".join(lines),va="top",fontsize=8.5,linespacing=1.5,color=COLORS["muted"])
    ax=fig.add_subplot(gs[3,3]); ax.axis("off")
    ax.text(0,.95,"DIFFERENCE SCALE\n0 to 1% of native range\n\nMagenta = invalid sample\n\nFull volume shown; gates use\nthe fixed native interior.",
            va="top",fontsize=8.5,linespacing=1.5,color=COLORS["muted"])
    fig.savefig(path,dpi=160,facecolor=fig.get_facecolor())
    plt.close(fig)


def pdf_report(out, summary):
    pdf = canvas.Canvas(str(out/"visual-qa.pdf"), pagesize=(1008,648))
    pdf.setTitle("neurotransform - Native transform interpretation QA")
    pdf.setAuthor("neurotransform visual QA")
    pdf.setFillColor(HexColor("#f4f7f8")); pdf.rect(0,0,1008,648,fill=1,stroke=0)
    pdf.setFillColor(HexColor(COLORS["ink"]))
    pdf.setFont("Helvetica-Bold",26); pdf.drawString(45,595,"Transform interpretation / visual QA")
    pdf.setFont("Helvetica",12); pdf.drawString(45,568,"Native references, fixed gates, and deliberately incorrect interpretations")
    pdf.setFont("Helvetica-Bold",14)
    pdf.drawString(45,523,f"{summary['passed']} passed    {summary['failed']} failed    {summary['blocked']} unavailable")
    pdf.setFont("Helvetica",10)
    y=490
    for case in summary["cases"]:
        color=COLORS["teal"] if case["status"]=="pass" else COLORS["red"] if case["status"]=="fail" else COLORS["amber"]
        pdf.setFillColor(HexColor(color)); pdf.drawString(45,y,case["status"].upper())
        pdf.setFillColor(HexColor(COLORS["ink"])); pdf.drawString(120,y,f"{case['family']} / {case['title']}")
        if case.get("metrics"):
            pdf.drawRightString(955,y,f"max coordinate error {num(case['metrics']['coordinate_max_mm'])} mm")
        y-=18
    pdf.setFont("Helvetica",10); pdf.setFillColor(HexColor(COLORS["muted"]))
    for line in ["Coordinate gate: <= 0.02 mm. Image NRMSE <= 0.0002; relative p99 residual <= 0.001.",
                 "Interior mask: native-resampled two-voxel source margin > 0.999, plus two target voxels.",
                 "Every interior sample must be finite. Native support determines the mask; the candidate does not.",
                 "FSL 5.0.9 is the locally available pinned producer. Newer FSL versions and FNIRT coefficients are unqualified.",
                 "Anatomical replays assess interpretation of saved transforms, not registration accuracy.",
                 "See index.html for synchronized slices, negative controls, landmarks, grids, and downloadable evidence."]:
        pdf.drawString(45,y-12,line); y-=17
    pdf.showPage()
    for i,case in enumerate(summary["cases"]):
        panel_file=out/"panels"/(case["id"]+".png")
        if panel_file.exists():
            pdf.drawImage(ImageReader(str(panel_file)),0,0,width=1008,height=648)
            pdf.setFont("Helvetica",8); pdf.setFillColor(HexColor(COLORS["muted"]))
            pdf.drawRightString(985,10,f"{i+2} / {len(summary['cases'])+1}")
            pdf.showPage()
        else:
            import textwrap
            pdf.setFillColor(HexColor(COLORS["ink"])); pdf.setFont("Helvetica-Bold",24)
            pdf.drawString(45,585,f"{case['family']} / {case['title']}")
            pdf.setFillColor(HexColor(COLORS["amber"])); pdf.setFont("Helvetica-Bold",16)
            pdf.drawString(45,530,"UNAVAILABLE - not counted as passing evidence")
            pdf.setFont("Helvetica",12); pdf.setFillColor(HexColor(COLORS["ink"]))
            for j,line in enumerate(textwrap.wrap(case.get("reason","No output"),width=110)):
                pdf.drawString(45,480-j*20,line)
            pdf.drawString(45,350,"Native outputs and commands remain in the case directory for inspection.")
            pdf.showPage()
    pdf.save()


def main():
    out=Path(sys.argv[1]).resolve()
    (out/"panels").mkdir(exist_ok=True)
    manifest=json.loads((out/"manifest.json").read_text())
    results={c["id"]:c for c in json.loads((out/"candidate-results.json").read_text())["cases"]}
    summary={"gates":manifest["gates"],"passed":0,"failed":0,"blocked":0,"cases":[],
             "candidate":manifest.get("candidate"),"images":manifest["images"],
             "mask":"Native support > 0.999 with two-voxel source and target margins; >= 1000 samples.",
             "versions":{"ANTs":"2.6.5.dev1-gfdce4d2", "AFNI":"26.1.04", "FSL":"5.0.9"},
             "limitations":["FSL 5.0.9 only; newer versions have not been checked.",
                "FNIRT coefficient fields and optimizer accuracy are not admitted by these dense-field replays.",
                "Image padding outside the native-defined interior is displayed but is not part of the numerical gates.",
                "Anatomical input images use a 4 mm replay grid; registration quality is not scored."]}
    visual=[]
    for config in manifest["cases"]:
        case=dict(config); folder=out/case["id"]
        result=results.get(case["id"],{})
        if case.get("missing") or result.get("status") != "evaluated":
            case.update(status="unavailable",reason=case.get("missing") or result.get("message","Candidate not evaluated"))
            summary["blocked"]+=1; summary["cases"].append(case); visual.append(case); continue
        im=nib.load(folder/"native_source.nii.gz")
        ref=im.get_fdata(dtype=np.float32); ours=load(folder,"ours_source.nii.gz")
        native_coords=np.stack([load(folder,f"native_coord{k}.nii.gz") for k in range(3)],axis=-1)
        ours_coords=load(folder,"ours_coords.nii.gz")
        mask=load(folder,"native_support.nii.gz")>.999
        interior=np.zeros(ref.shape,dtype=bool); interior[2:-2,2:-2,2:-2]=True
        mask &= interior
        metrics,error,residual=score(ref,ours,native_coords,ours_coords,mask,manifest["gates"])
        window=[0.,float(np.quantile(load(folder,"source.nii.gz"),.998))]
        if window[1]<=0: window[1]=1.
        error_scale=max(metrics["intensity_range"]*.01,1e-6)
        data={"native":pack(ref,*window),"ours":pack(ours,*window),"error":pack(residual,0,error_scale),
              "coord_error":pack(np.where(mask,error,np.nan),0,.1),
              "jacobian":pack(np.where(mask,load(folder,"ours_jacobian.nii.gz"),np.nan),.5,1.5)}
        case["variant_metrics"]=[]
        variant_coordinates={}
        for v in result["variants"]:
            image=load(folder,v["id"]+".nii.gz"); coords=load(folder,v["id"]+"_coords.nii.gz")
            variant_coordinates[v["id"]]=coords
            vm,ve,vr=score(ref,image,native_coords,coords,mask,manifest["gates"])
            rejected=(not vm["pass"] and (vm["invalid_samples"]>0 or
                       vm["coordinate_max_mm"] is None or vm["coordinate_max_mm"]>=manifest["gates"]["negative_control_min_mm"]))
            case["variant_metrics"].append({**v,"metrics":vm,"rejected":rejected})
            data[v["id"]]=pack(image,*window)
            data[v["id"]+"_error"]=pack(vr,0,error_scale)
            data[v["id"]+"_coord_error"]=pack(np.where(mask,ve,np.nan),0,.1)
        case["metrics"]=metrics
        if case["id"] in ("fsl_relative", "fsl_absolute"):
            # Diagnostic only: for this matched source/target left-handed grid,
            # independently express scaled-voxel components in source RAS.
            # This does NOT replace the package output or turn its failure green.
            src=nib.load(folder/"source.nii.gz")
            aff=src.affine
            d=load(folder,"warp.nii.gz")
            basis=aff[:3,:3]/np.linalg.norm(aff[:3,:3],axis=0)
            if case["representation"]=="relative":
                ijk=np.stack(np.meshgrid(*[np.arange(n) for n in src.shape],indexing="ij"),axis=-1)
                diagnostic=d@basis.T + ijk@aff[:3,:3].T + aff[:3,3]
            else:
                diagnostic=d@basis.T + aff[:3,3]
            error_mm=float(np.linalg.norm(diagnostic-native_coords,axis=-1)[mask].max())
            case["finding"]={"summary":"The dense FSL loader treats scaled voxel vectors/coordinates as RAS. This case requires a basis conversion"+
                            (" and the source-grid origin." if case["representation"]=="absolute" else "."),
                            "diagnostic_basis_error_mm":error_mm,
                            "qualification":"Diagnostic calculation only. The actual package output remains the scored candidate."}
            if metrics["pass"]:
                case["finding"]["summary"] = "Regression fixed: FSL scaled-voxel values are converted using source and reference geometry before evaluating the RAS displacement."
                case["finding"]["resolved"] = True
        # A numerically passing interpretation cannot receive PASS if its
        # deliberate-error controls are ineffective.
        case["status"]="pass" if metrics["pass"] and all(v["rejected"] for v in case["variant_metrics"]) else "fail"
        summary["passed" if case["status"]=="pass" else "failed"]+=1
        panel(out/"panels"/(case["id"]+".png"),case,
              {"native":ref,"ours":ours,"error":residual},metrics,native_coords,ours_coords,mask,window,
              np.linalg.norm(im.affine[:3,:3],axis=0), nib.aff2axcodes(im.affine))
        summary["cases"].append(case)
        # Fixed target landmarks; invalid or excluded points remain marked.
        n=np.array(ref.shape)
        landmark_indices=np.round((n-1)*np.array([[.5,.5,.5],[.33,.5,.5],[.67,.5,.5],
                        [.5,.33,.5],[.5,.67,.5],[.5,.5,.33],[.5,.5,.67]])).astype(int)
        landmarks=[]
        for j,ijk in enumerate(landmark_indices):
            idx=tuple(ijk)
            landmarks.append({"id":j+1,"ijk":ijk.tolist(),"native":native_coords[idx].tolist(),
                               "ours":ours_coords[idx].tolist(),"inside":bool(mask[idx]),
                               "alternatives":{key:arr[idx].tolist() for key,arr in variant_coordinates.items()}})
        # Middle axial deformation lattice expressed as source RAS x/y.
        grid=[]
        for axis in [0,1]:
            for fixed in range(5,int(n[axis])-4,5):
                inds=[]
                for varying in range(3,int(n[1-axis])-3):
                    ijk=[0,0,int(n[2]//2)]; ijk[axis]=fixed; ijk[1-axis]=varying; inds.append(tuple(ijk))
                grid.append({"native":[native_coords[k][:2].tolist() if mask[k] else None for k in inds],
                             "ours":[ours_coords[k][:2].tolist() if mask[k] else None for k in inds],
                             "alternatives":{key:[arr[k][:2].tolist() if mask[k] else None for k in inds]
                                             for key,arr in variant_coordinates.items()}})
        visual.append({**case,"dims":list(ref.shape),"affine":im.affine.tolist(),
                       "axis_codes":list(nib.aff2axcodes(im.affine)),"window":window,"error_scale":error_scale,
                       "volumes":data,"landmarks":landmarks,"grid":grid})
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    payload={**summary,"cases":visual}
    # Non-finite diagnostics are explicit nulls in the display payload.
    raw=json.dumps(json_safe(payload),allow_nan=False,separators=(",",":"))
    (out/"data.js").write_text("window.QA="+raw+";\n")
    shutil.copyfile(HERE/"report.html",out/"index.html")
    pdf_report(out,summary)
    print(json.dumps({k:summary[k] for k in ["passed","failed","blocked"]}),flush=True)


if __name__=="__main__":
    main()
