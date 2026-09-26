"""Independent local least-squares oracle for retained full-template failures.
Run after /tmp/diagnose-barycentric.R writes per-case diagnose.json. This is a
failure investigation, not an alternative threshold or a qualification pass.
"""
import json, sys, subprocess
from pathlib import Path
import numpy as np
from generate_barycentric_oracle import read_gifti, write_gifti, sha

base, wb = Path(sys.argv[1]), Path(sys.argv[2])
receipt=json.loads((base/'receipt.json').read_text())
reports=[]
for case in receipt['cases']:
 folder=base/case['name']; diag=json.loads((folder/'diagnose.json').read_text())
 vertices,faces=read_gifti(case['source']); query=read_gifti(case['target'])[0][np.array(diag['indices'])-1]
 unit=vertices.astype(float); unit=100*unit/np.linalg.norm(unit,axis=1)[:,None]
 qp=query.astype(float); qp=100*qp/np.linalg.norm(qp,axis=1)[:,None]
 triangles=unit[faces]; lo=triangles.min(axis=1); hi=triangles.max(axis=1)
 records=[]; all_nodes=set()
 for row,p in zip(diag['indices'],qp):
  upper=np.min(np.sum((unit-p)**2,axis=1))
  lower=np.sum(np.maximum(0,np.maximum(lo-p,p-hi))**2,axis=1)
  candidates=np.where(lower<=upper)[0]
  best=None
  for fi in candidates:
   t=triangles[fi]; edge=(t[1:]-t[0]).T
   vw=np.linalg.lstsq(edge,p-t[0],rcond=None)[0]
   bary=np.r_[1-vw.sum(),vw]
   trials=[bary] if min(bary)>=0 else []
   for a,b in ((0,1),(1,2),(2,0)):
    e=t[b]-t[a]; u=np.clip(np.dot(e,p-t[a])/np.dot(e,e),0,1)
    w=np.zeros(3);w[a]=1-u;w[b]=u;trials.append(w)
   for w in trials:
    d=np.sum((p-w@t)**2)
    if best is None or d<best[0]:best=(d,fi,w)
  native={col-1:w for rr,col,w in zip(diag['rows'],diag['cols'],diag['weights']) if rr==row}
  oracle={int(c):float(w) for c,w in zip(faces[best[1]],best[2]) if w>0}
  all_nodes.update(faces[candidates].ravel().tolist())
  records.append(dict(target_one_based=row,native=native,oracle=oracle,distance=float(np.sqrt(best[0])),candidates=len(candidates),
                      native_oracle_max_error=max(abs(native.get(c,0)-oracle.get(c,0)) for c in set(native)|set(oracle))))
 nodes=sorted(all_nodes)
 # Add three filler vertices solely to give Workbench a triangle in the target file.
 target=np.vstack((query,read_gifti(case['target'])[0][:3]))
 write_gifti(folder/'diagnostic-target.surf.gii', [('NIFTI_INTENT_POINTSET',target),('NIFTI_INTENT_TRIANGLE',np.array([[len(query),len(query)+1,len(query)+2]]))])
 impulse=np.zeros((len(vertices),len(nodes)),dtype=np.float32);impulse[nodes,np.arange(len(nodes))]=1
 write_gifti(folder/'diagnostic-impulses.func.gii',[('NIFTI_INTENT_SHAPE',c) for c in impulse.T])
 argv=[str(wb),'-metric-resample',str(folder/'diagnostic-impulses.func.gii'),case['source'],str(folder/'diagnostic-target.surf.gii'),'BARYCENTRIC',str(folder/'diagnostic-output.func.gii')]
 run=subprocess.run(argv,capture_output=True,text=True,check=True)
 expected=np.column_stack(read_gifti(folder/'diagnostic-output.func.gii'))[:len(query)]
 # Reproduce Workbench's float32 radius normalization, then independently
 # compare its selected edge with the native face in double precision.
 vf=vertices*(np.float32(100)/np.sqrt(np.sum(vertices*vertices,axis=1)))[:,None]
 qf=query*(np.float32(100)/np.sqrt(np.sum(query*query,axis=1)))[:,None]
 for r,e,pf in zip(records,expected,qf.astype(float)):
  r['workbench']={int(c):float(w) for c,w in zip(nodes,e) if w>0}
  r['workbench_row_sum']=float(sum(e)); r['weight_l1_difference']=sum(abs(r['native'].get(c,0)-r['workbench'].get(c,0)) for c in set(r['native'])|set(r['workbench']))
  if len(r['workbench'])==2 and len(r['native'])==3:
   corners=vf[list(r['native'])].astype(float)
   vw=np.linalg.lstsq((corners[1:]-corners[0]).T,pf-corners[0],rcond=None)[0]
   bw=np.r_[1-vw.sum(),vw]
   ends=vf[list(r['workbench'])].astype(float); edge=ends[1]-ends[0]
   t=np.clip(np.dot(pf-ends[0],edge)/np.dot(edge,edge),0,1)
   r['float_geometry_face_weights']=bw.tolist()
   r['float_geometry_edge_minus_face_distance2']=float(np.sum((pf-ends[0]-t*edge)**2)-np.sum((pf-bw@corners)**2))
 report=dict(name=case['name'],queries=records,argv=argv,exit=run.returncode,stderr=run.stderr)
 reports.append(report)
 (folder/'diagnostic-result.json').write_text(json.dumps(report,indent=2)+'\n')
 print(case['name'],len(records),'queries; native/oracle max',max(r['native_oracle_max_error'] for r in records),flush=True)
(base/'diagnostic-receipt.json').write_text(json.dumps(dict(generator_sha256=sha(__file__),cases=reports),indent=2)+'\n')
