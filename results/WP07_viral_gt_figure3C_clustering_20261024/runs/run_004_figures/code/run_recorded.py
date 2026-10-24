"""Record code/config/runtime/input identity BEFORE starting any scientific subprocess."""
import argparse,datetime,hashlib,importlib.metadata,importlib,inspect,json,os,platform,shutil,subprocess,sys,time
from pathlib import Path
import common as c
ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['science','figures','validation']);ap.add_argument('script');a=ap.parse_args()
for d in ('runs','tables','graphs','assignments','stability','provenance','figures','runtime/matplotlib_cache'):(c.OUT/d).mkdir(parents=True,exist_ok=True)
existing=list((c.OUT/'runs').glob('run_*'));run_id=f'run_{len(existing)+1:03d}_{a.stage}';run=c.OUT/'runs'/run_id;run.mkdir();(run/'code').mkdir()
inputs=c.read(c.OUT/'INPUT_MANIFEST.csv')
for r in inputs:
 assert c.sha(Path(r['path']))==r['sha256'],r['path']
sourcefiles=list((c.OUT/'scripts').glob('*.py'))+[c.OUT/'config.json',c.OUT/'INPUT_MANIFEST.csv']
sources=[]
for p in sourcefiles:
 dest=run/'code'/p.name;shutil.copyfile(p,dest);sources.append(dict(path=str(p),snapshot=dest.relative_to(c.OUT).as_posix(),sha256=c.sha(p)))
versions={p:importlib.metadata.version(p) for p in ['numpy','scipy','scikit-learn','igraph','leidenalg','texttable','matplotlib','pillow']}
# Capture distribution/module file identities at execution time, not retrospectively.
package_hashes=[]
for package in ['numpy','scipy','scikit-learn','igraph','leidenalg','matplotlib']:
 dist=importlib.metadata.distribution(package)
 for rel in dist.files or []:
  p=Path(dist.locate_file(rel))
  if p.is_file() and (p.suffix.lower() in ('.py','.pyd','.dll') or p.name in ('METADATA','RECORD')):
   package_hashes.append(dict(distribution=package,path=str(p),size_bytes=p.stat().st_size,sha256=c.sha(p)))
c.write(f'runs/{run_id}/dependency_file_hashes.csv',package_hashes)
modules={name:importlib.import_module(name) for name in ['numpy','scipy','sklearn','igraph','leidenalg','matplotlib']}
from sklearn.cluster import HDBSCAN
api=dict(HDBSCAN_signature=str(inspect.signature(HDBSCAN)),Leiden_optimise_signature=str(inspect.signature(modules['leidenalg'].Optimiser.optimise_partition)),Leiden_optimizer_has_max_comm_size=hasattr(modules['leidenalg'].Optimiser(),'max_comm_size'))
command=[sys.executable,'-B',str(c.OUT/'scripts'/a.script)]
record=dict(run_id=run_id,stage=a.stage,status='PREPARED',command=command,cwd=str(c.ROOT),python=sys.version,executable=sys.executable,platform=platform.platform(),started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),configuration=c.CONFIG,configuration_sha256=c.sha(c.OUT/'config.json'),source_files=sources,input_manifest_sha256=c.sha(c.OUT/'INPUT_MANIFEST.csv'),software_versions=versions,distribution_inventory=sorted({d.metadata['Name']:d.version for d in importlib.metadata.distributions() if d.metadata.get('Name')}.items()),loaded_module_paths={n:m.__file__ for n,m in modules.items()},explicit_environment={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','MPLCONFIGDIR']},api=api,dependency_file_manifest=f'runs/{run_id}/dependency_file_hashes.csv')
(run/'before.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
print('PRE-EXECUTION RECORD SAVED',run_id,'source',c.sha(c.OUT/'scripts'/a.script),'versions',versions,flush=True)
start=time.monotonic()
with (run/'stdout.txt').open('wb') as stdout,(run/'stderr.txt').open('wb') as stderr:
 p=subprocess.run(command,cwd=c.ROOT,stdout=stdout,stderr=stderr,env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1'})
unchanged=[]
for r in inputs:unchanged.append(dict(path=r['path'],before=r['sha256'],after=c.sha(Path(r['path']))))
source_unchanged=all(c.sha(Path(r['path']))==r['sha256'] for r in sources)
record.update(status='SUCCEEDED' if p.returncode==0 else 'FAILED',exit_status=p.returncode,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start,inputs_unchanged=all(r['before']==r['after'] for r in unchanged),executable_sources_unchanged=source_unchanged)
c.write(f'runs/{run_id}/input_integrity_after.csv',unchanged)
outputs=[dict(path=p.relative_to(c.OUT).as_posix(),size_bytes=p.stat().st_size,sha256=c.sha(p)) for p in sorted(c.OUT.rglob('*')) if p.is_file() and not any(part in ('runtime','runs') for part in p.relative_to(c.OUT).parts)]
c.write(f'runs/{run_id}/output_manifest.csv',outputs)
record['output_manifest']=f'runs/{run_id}/output_manifest.csv';(run/'after.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
print('POST-EXECUTION',run_id,'exit',p.returncode,'inputs unchanged',record['inputs_unchanged'],'source unchanged',source_unchanged,flush=True)
print((run/'stdout.txt').read_text(encoding='utf-8',errors='replace')[-6000:],flush=True)
print((run/'stderr.txt').read_text(encoding='utf-8',errors='replace')[-3000:],flush=True)
if p.returncode!=0:raise SystemExit(p.returncode)
assert record['inputs_unchanged'] and source_unchanged
