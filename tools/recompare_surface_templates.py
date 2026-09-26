"""Reapply a frozen engine to retained Workbench outputs without regenerating them.
Usage: python3 tools/recompare_surface_templates.py RECEIPT NEW_OUTPUT [R_SCRIPT]
Omit R_SCRIPT to reuse the ordinary comparison script from the original attempt.
NEUROTRANSFORM_BUILD_MANIFEST records the tested engine's source and DLL hashes.
"""
import json,os,sys,subprocess,platform,re
from pathlib import Path
from generate_barycentric_oracle import sha
source,out=Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve()
original=json.loads(source.read_text());out.mkdir(parents=True,exist_ok=False)
script=Path(sys.argv[3]).resolve() if len(sys.argv)>3 else source.parent/'compare.R'
code=script.read_text()
# Older ordinary harnesses wrote into the input directory. Route their output
# to the new attempt while preserving the original comparison and its failure.
code=code.replace('file.path(folder, "result.json")','if(length(args)>2L) args[3] else file.path(folder, "result.json")')
(out/'compare.R').write_text(code)
receipt={k:v for k,v in original.items() if k not in ('cases','build_manifest','revision')}
receipt.update(input_receipt=str(source),input_receipt_sha256=sha(source),comparison_sha256=sha(out/'compare.R'),
 generator_sha256=sha(__file__),revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
 build_manifest=json.loads(Path(os.environ['NEUROTRANSFORM_BUILD_MANIFEST']).read_text()),cases=[])
for case in original['cases']:
 folder=Path(case['comparison_argv'][-2]);target=out/case['name'];target.mkdir()
 # Bind the retained tabular inputs to this reapplication attempt.
 input_hashes={p.name:sha(p) for p in folder.glob('*.csv')}
 argv=['/usr/bin/time','-l' if platform.system()=='Darwin' else '-v','Rscript',str(out/'compare.R'),str(folder),str(original['tolerance']),str(target/'result.json')]
 run=subprocess.run(argv,capture_output=True,text=True,env=os.environ.copy())
 (target/'comparison.log').write_text(run.stdout+run.stderr)
 result=json.loads((target/'result.json').read_text()) if (target/'result.json').exists() else {}
 rss=re.search(r'(\d+)\s+maximum resident set size',run.stderr)
 if rss:result['peak_resident_bytes']=int(rss.group(1))
 rss=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',run.stderr)
 if rss:result['peak_resident_bytes']=1024*int(rss.group(1))
 entry={k:v for k,v in case.items() if k not in ('result','comparison_argv','comparison_exit')}
 entry.update(comparison_argv=argv,comparison_exit=run.returncode,input_csv_sha256=input_hashes,result=result)
 receipt['cases'].append(entry)
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(case['name'],result.get('max_absolute_error'),result.get('label_mismatches'),'exit',run.returncode,flush=True)
if any(c['comparison_exit']!=0 for c in receipt['cases']):raise SystemExit(1)
