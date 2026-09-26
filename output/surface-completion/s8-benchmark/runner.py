import json,pathlib,subprocess,os,re,hashlib
out=pathlib.Path('output/surface-completion/s8-benchmark');out.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
receipt=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 script_sha256=sha('tools/benchmark_surface_plan.R'),threads=os.environ.get('OMP_NUM_THREADS'),
 budget=dict(construction_seconds=60,additional_resident_bytes=1073741824),cases=[])
for file in sorted(pathlib.Path('/tmp').glob('neurotransform-benchmark-*.rds')):
 name=file.stem.replace('neurotransform-benchmark-','')
 argv=['/usr/bin/time','-l','Rscript','tools/benchmark_surface_plan.R',str(file),str(out/(name+'.json'))]
 run=subprocess.run(argv,capture_output=True,text=True,env=os.environ.copy())
 (out/(name+'.log')).write_text(run.stdout+run.stderr)
 result=json.loads((out/(name+'.json')).read_text()) if (out/(name+'.json')).exists() else {}
 peak=re.search(r'(\d+)\s+maximum resident set size',run.stderr)
 if peak:
  result['peak_resident_bytes']=int(peak.group(1))
  result['additional_peak_resident_bytes']=max(0,result['peak_resident_bytes']-result['baseline_resident_bytes'])
 result['passed']=run.returncode==0 and result.get('passed',False) and result.get('additional_peak_resident_bytes',float('inf'))<1073741824
 receipt['cases'].append(dict(name=name,input_sha256=sha(file),argv=argv,exit=run.returncode,result=result))
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(name,result.get('construction_seconds'),result.get('additional_peak_resident_bytes'),result['passed'],flush=True)
if not all(c['result']['passed'] for c in receipt['cases']):raise SystemExit(1)
