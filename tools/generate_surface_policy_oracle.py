"""Generate Workbench mask/label fixtures independent of neurotransform.
Usage: python3 tools/generate_surface_policy_oracle.py WB NEW_OUTPUT
"""
import json,sys,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from generate_barycentric_oracle import sphere,read_gifti,write_gifti,sha

wb,out=Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve()
out.mkdir(parents=True,exist_ok=False)
v,f=sphere(1); q,qf=sphere(2)
v=v@np.array([[1,.17,.09],[.04,1.1,.12],[.02,-.07,.9]])
q=q@np.array([[1,.231,.173],[-.219,1,.13],[-.15,-.127,1]])
v=(100*v/np.linalg.norm(v,axis=1)[:,None]).astype(np.float32)
q=(100*q/np.linalg.norm(q,axis=1)[:,None]).astype(np.float32)
write_gifti(out/'source.surf.gii',[('NIFTI_INTENT_POINTSET',v),('NIFTI_INTENT_TRIANGLE',f)])
write_gifti(out/'target.surf.gii',[('NIFTI_INTENT_POINTSET',q),('NIFTI_INTENT_TRIANGLE',qf)])
labels=np.array([10,10,20,40,20,40,10,10,20,20,40,40,10,20,40,10,20,40],dtype=np.int32)
write_gifti(out/'labels.label.gii',[('NIFTI_INTENT_LABEL',labels)])
tree=ET.parse(out/'labels.label.gii');root=tree.getroot()
root.find('DataArray').set('DataType','NIFTI_TYPE_INT32')
for key,name,color in [(0,'unassigned',(0,0,0)),(10,'A',(1,0,0)),(20,'B',(0,1,0)),(40,'C',(0,0,1))]:
 ET.SubElement(root.find('LabelTable'),'Label',Key=str(key),Red=str(color[0]),Green=str(color[1]),Blue=str(color[2]),Alpha='1').text=name
tree.write(out/'labels.label.gii',encoding='UTF-8',xml_declaration=True)
# Basis columns prove every nonzero interpolation weight, plus a nonconstant field.
data=np.column_stack((np.eye(len(v)),np.sin(np.arange(len(v))*.37)))
write_gifti(out/'input.func.gii',[('NIFTI_INTENT_SHAPE',col) for col in data.T])
cases=[];commands=[]
for name,roi in [('all',np.ones(len(v))),('hemisphere',np.where(v[:,2]>-20,.2,-1)),('isolated',np.arange(len(v))==7),('disconnected',np.isin(np.arange(len(v)),[0,1,12,16])),('empty',np.zeros(len(v)))]:
 write_gifti(out/f'{name}-roi.func.gii',[('NIFTI_INTENT_SHAPE',roi)])
 case=dict(name=name,source_mask=roi.astype(float).tolist())
 for kind,extra in [('metric',[]),('aggregate',[]),('largest',['-largest'])]:
  metric=kind=='metric';output=out/f'{name}-{kind}.gii';valid=out/f'{name}-{kind}-valid.func.gii'
  argv=[str(wb),'-metric-resample' if metric else '-label-resample',str(out/('input.func.gii' if metric else 'labels.label.gii')),str(out/'source.surf.gii'),str(out/'target.surf.gii'),'BARYCENTRIC',str(output),'-current-roi',str(out/f'{name}-roi.func.gii'),'-valid-roi-out',str(valid)]+extra
  run=subprocess.run(argv,capture_output=True,text=True,check=True)
  commands.append(dict(argv=argv,exit=run.returncode,stderr=run.stderr))
  case[kind]=np.column_stack(read_gifti(output)).tolist() if metric else read_gifti(output)[0].tolist()
  case['valid']=read_gifti(valid)[0].tolist()
 cases.append(case)
fixture=dict(method='BARYCENTRIC',tolerance=2e-6,vertices=v.tolist(),faces=f.tolist(),query=q.tolist(),query_faces=qf.tolist(),data=data.tolist(),labels=labels.tolist(),cases=cases,commands=commands,
 version=subprocess.check_output([str(wb),'-version'],text=True),binary_sha256=sha(wb),generator_sha256=sha(__file__),files={p.name:sha(p) for p in out.glob('*.gii')})
(out/'oracle.json').write_text(json.dumps(fixture,indent=2)+'\n')
print(out/'oracle.json')
