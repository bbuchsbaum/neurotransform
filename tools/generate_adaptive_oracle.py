"""Generate unequal-area ADAP_BARY_AREA basis, ROI and label oracles.
Usage: python3 tools/generate_adaptive_oracle.py WB NEW_OUTPUT
"""
import json,sys,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from generate_barycentric_oracle import sphere,read_gifti,write_gifti,sha

wb,out=Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve()
out.mkdir(parents=True,exist_ok=False)
coarse=np.array([[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]],dtype=float)
coarse/=np.linalg.norm(coarse,axis=1)[:,None]
cf=np.array([[0,1,2],[0,3,1],[0,2,3],[1,3,2]],dtype=np.int32)
center=coarse[:3].sum(axis=0);center/=np.linalg.norm(center)
fine=np.vstack((coarse,center));ff=np.vstack((cf[1:],[[0,1,4],[1,2,4],[2,0,4]]))
small,sf=sphere(1);large,lf=sphere(2)
small=small@np.array([[1,.17,.09],[.04,1.1,.12],[.02,-.07,.9]])
large=large@np.array([[1,.231,.173],[-.219,1,.13],[-.15,-.127,1]])
geometries=[('tetra-up',coarse,cf,fine,ff),('tetra-down',fine,ff,coarse,cf),('irregular-up',small,sf,large,lf),('irregular-down',large,lf,small,sf)]
cases=[];commands=[]
for name,v,f,q,qf in geometries:
 v=(100*v/np.linalg.norm(v,axis=1)[:,None]).astype(np.float32)
 q=(100*q/np.linalg.norm(q,axis=1)[:,None]).astype(np.float32)
 a=(1+np.arange(len(v))*.31+np.sin(np.arange(len(v)))**2).astype(np.float32)
 b=(2+np.arange(len(q))*.27+np.cos(np.arange(len(q)))**2).astype(np.float32)
 data=np.column_stack((np.eye(len(v)),np.sin(np.arange(len(v))*.37),np.ones(len(v))))
 labels=(10*(1+(np.arange(len(v))//2)%3)).astype(np.int32)
 folder=out/name;folder.mkdir()
 for filename,verts,faces in [('source',v,f),('target',q,qf)]:
  write_gifti(folder/f'{filename}.surf.gii',[('NIFTI_INTENT_POINTSET',verts),('NIFTI_INTENT_TRIANGLE',faces)])
 for filename,values in [('source-area',a),('target-area',b)]:
  write_gifti(folder/f'{filename}.func.gii',[('NIFTI_INTENT_SHAPE',values)])
 write_gifti(folder/'input.func.gii',[('NIFTI_INTENT_SHAPE',col) for col in data.T])
 write_gifti(folder/'labels.label.gii',[('NIFTI_INTENT_LABEL',labels)])
 tree=ET.parse(folder/'labels.label.gii');root=tree.getroot();root.find('DataArray').set('DataType','NIFTI_TYPE_INT32')
 for key,color in [(0,(0,0,0)),(10,(1,0,0)),(20,(0,1,0)),(30,(0,0,1))]:
  ET.SubElement(root.find('LabelTable'),'Label',Key=str(key),Red=str(color[0]),Green=str(color[1]),Blue=str(color[2]),Alpha='1').text='label'+str(key)
 tree.write(folder/'labels.label.gii',encoding='UTF-8',xml_declaration=True)
 for maskname,roi in [('all',np.ones(len(v))),('partial',np.where(np.arange(len(v))%3,.2,0)),('empty',np.zeros(len(v)))]:
  write_gifti(folder/f'{maskname}-roi.func.gii',[('NIFTI_INTENT_SHAPE',roi)])
  case=dict(name=name+'-'+maskname,vertices=v.tolist(),faces=f.tolist(),query=q.tolist(),query_faces=qf.tolist(),source_areas=a.tolist(),target_areas=b.tolist(),data=data.tolist(),labels=labels.tolist(),source_mask=roi.tolist())
  for kind,extra in [('metric',[]),('aggregate',[]),('largest',['-largest'])]:
   metric=kind=='metric';output=folder/f'{maskname}-{kind}.gii';valid=folder/f'{maskname}-{kind}-valid.func.gii'
   argv=[str(wb),'-metric-resample' if metric else '-label-resample',str(folder/('input.func.gii' if metric else 'labels.label.gii')),str(folder/'source.surf.gii'),str(folder/'target.surf.gii'),'ADAP_BARY_AREA',str(output),'-area-metrics',str(folder/'source-area.func.gii'),str(folder/'target-area.func.gii'),'-current-roi',str(folder/f'{maskname}-roi.func.gii'),'-valid-roi-out',str(valid)]+extra
   run=subprocess.run(argv,capture_output=True,text=True,check=True)
   commands.append(dict(argv=argv,exit=run.returncode,stderr=run.stderr))
   case[kind]=np.column_stack(read_gifti(output)).tolist() if metric else read_gifti(output)[0].tolist()
   case['valid']=read_gifti(valid)[0].tolist()
  cases.append(case)
fixture=dict(method='ADAP_BARY_AREA',tolerance=2e-6,cases=cases,commands=commands,
 version=subprocess.check_output([str(wb),'-version'],text=True),binary_sha256=sha(wb),generator_sha256=sha(__file__),files={str(p.relative_to(out)):sha(p) for p in out.glob('*/*.gii')})
(out/'oracle.json').write_text(json.dumps(fixture,indent=2)+'\n')
print(out/'oracle.json')
