"""Full-template ADAP_BARY_AREA comparison with anatomical areas, ROIs, labels.
Usage: python3 tools/compare_adaptive_templates.py WB INPUTS INPUT_LOCK NEW_OUTPUT
R_LIBS must point to the frozen rebuilt engine. All input hashes must match lock.
"""
import json,os,sys,subprocess,platform,re,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from generate_barycentric_oracle import read_gifti,write_gifti,sha
wb,inputs,lockpath,out=map(lambda p:Path(p).resolve(),sys.argv[1:])
out.mkdir(parents=True,exist_ok=False)
lock=json.loads(lockpath.read_text());locked={a['path']:a for a in lock['assets']}
def check(path):
 assert sha(path)==locked[path.name]['sha256'], 'Input hash mismatch: '+str(path)
 return sha(path)
script=Path(__file__).with_suffix('.R'); (out/'compare.R').write_bytes(script.read_bytes())
receipt=dict(method='ADAP_BARY_AREA',tolerance=5e-5,seed=20926,full=True,
  platform=platform.platform(),threads=os.environ.get('OMP_NUM_THREADS'),
  revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
  version=subprocess.check_output([str(wb),'-version'],text=True),binary_sha256=sha(wb),
  generator_sha256=sha(__file__),comparison_sha256=sha(script),input_lock_sha256=sha(lockpath),cases=[])
manifest_path=os.environ.get('NEUROTRANSFORM_BUILD_MANIFEST')
if manifest_path:receipt['build_manifest']=json.loads(Path(manifest_path).read_text())
for hemi in ('L','R'):
 spheres=[inputs/f'tpl-fsaverage_hemi-{hemi}_den-164k_sphere.surf.gii',inputs/f'tpl-fsLR_space-fsaverage_hemi-{hemi}_den-32k_sphere.surf.gii']
 areas=[inputs/f'tpl-fsaverage_hemi-{hemi}_den-164k_desc-vaavg_midthickness.shape.gii',inputs/f'tpl-fsLR_hemi-{hemi}_den-32k_desc-vaavg_midthickness.shape.gii']
 masks=[inputs/f'tpl-fsaverage_den-164k_hemi-{hemi}_desc-nomedialwall_dparc.label.gii',inputs/f'tpl-fsLR_hemi-{hemi}_den-32k_desc-nomedialwall_dparc.label.gii']
 for direction in (0,1):
  source,target=spheres[direction],spheres[1-direction]
  sa,ta=areas[direction],areas[1-direction];mask=masks[direction]
  hashes={str(p):check(p) for p in (source,target,sa,ta,mask)}
  v,f=read_gifti(source);q,qf=read_gifti(target)
  a=read_gifti(sa)[0].reshape(-1); b=read_gifti(ta)[0].reshape(-1)
  assert len(a)==len(v) and len(b)==len(q) and np.all(a>0) and np.all(b>0)
  p=v/np.linalg.norm(v,axis=1)[:,None];rng=np.random.default_rng(receipt['seed'])
  impulses=np.zeros((len(v),4),dtype=np.float32); impulses[rng.choice(len(v),4,replace=False),np.arange(4)]=1
  data=np.column_stack([p,p[:,0]*p[:,1],np.sin(np.arange(len(v))*.017),np.cos(np.arange(len(v))*.031),rng.uniform(-1,1,len(v)),impulses,np.ones(len(v))]).astype(np.float32)
  labels=(10*rng.integers(1,5,len(v))).astype(np.int32)
  for roi_name in ('all','medial-wall'):
   name=hemi+('-164k-to-32k-' if direction==0 else '-32k-to-164k-')+roi_name
   folder=out/name;folder.mkdir()
   roi=np.ones(len(v),dtype=np.float32) if roi_name=='all' else read_gifti(mask)[0].reshape(-1).astype(np.float32)
   write_gifti(folder/'roi.func.gii',[('NIFTI_INTENT_SHAPE',roi)])
   write_gifti(folder/'input.func.gii',[('NIFTI_INTENT_SHAPE',col) for col in data.T])
   write_gifti(folder/'input.label.gii',[('NIFTI_INTENT_LABEL',labels)])
   tree=ET.parse(folder/'input.label.gii');root=tree.getroot();root.find('DataArray').set('DataType','NIFTI_TYPE_INT32')
   for key in (0,10,20,30,40):
    ET.SubElement(root.find('LabelTable'),'Label',Key=str(key),Red=str(key/40),Green='0',Blue='0',Alpha='1').text='label'+str(key)
   tree.write(folder/'input.label.gii',encoding='UTF-8',xml_declaration=True)
   commands=[];expected_labels=[]
   for kind,extra in [('metric',[]),('aggregate',[]),('largest',['-largest'])]:
    metric=kind=='metric';output=folder/(kind+'.gii');valid=folder/(kind+'-valid.func.gii')
    argv=[str(wb),'-metric-resample' if metric else '-label-resample',str(folder/('input.func.gii' if metric else 'input.label.gii')),str(source),str(target),'ADAP_BARY_AREA',str(output),'-area-metrics',str(sa),str(ta),'-current-roi',str(folder/'roi.func.gii'),'-valid-roi-out',str(valid)]+extra
    run=subprocess.run(argv,capture_output=True,text=True,check=True)
    commands.append(dict(argv=argv,exit=run.returncode,stderr=run.stderr))
    if metric:
     expected=np.column_stack(read_gifti(output));valid_values=read_gifti(valid)[0].reshape(-1)
    else:expected_labels.append(read_gifti(output)[0].reshape(-1))
   for filename,values in [('vertices',v),('faces',f),('query',q),('query_faces',qf),('data',data),('expected',expected),('source_area',a),('target_area',b),('roi',roi),('labels',labels),('expected_labels',np.column_stack(expected_labels)),('valid',valid_values)]:
    np.savetxt(folder/(filename+'.csv'),values,delimiter=',',fmt='%.17g')
   (folder/'config.json').write_text(json.dumps(dict(source_area_sha256=hashes[str(sa)],target_area_sha256=hashes[str(ta)])))
   argv=['/usr/bin/time','-l' if platform.system()=='Darwin' else '-v','Rscript',str(out/'compare.R'),str(folder),str(receipt['tolerance'])]
   run=subprocess.run(argv,capture_output=True,text=True,env=os.environ.copy())
   (folder/'comparison.log').write_text(run.stdout+run.stderr)
   result=json.loads((folder/'result.json').read_text()) if (folder/'result.json').exists() else {}
   rss=re.search(r'(\d+)\s+maximum resident set size',run.stderr)
   if rss:result['peak_resident_bytes']=int(rss.group(1))
   rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',run.stderr)
   if rss:result['peak_resident_bytes']=1024*int(rss.group(1))
   receipt['cases'].append(dict(name=name,input_hashes=hashes,commands=commands,comparison_argv=argv,comparison_exit=run.returncode,result=result))
   (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
   print(name,result.get('max_absolute_error'),result.get('label_mismatches'),'exit',run.returncode,flush=True)
if any(c['comparison_exit']!=0 for c in receipt['cases']):raise SystemExit(1)
