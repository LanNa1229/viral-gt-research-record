"""Record a file-based supplementary stability calculation without any fitting."""
import datetime,hashlib,json,shutil,subprocess,sys
from pathlib import Path
OUT=Path(__file__).resolve().parents[1];run=OUT/'runs'/f'run_{len(list((OUT/"runs").glob("run_*")))+1:03d}_supplementary_stability';run.mkdir();(run/'code').mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
sources=[OUT/'scripts/cross_grid_audit.py',OUT/'scripts/common.py',Path(__file__),OUT/'config.json'];inputs=list((OUT/'assignments').glob('*_leiden_grid.csv'))+[OUT/'tables/duplicate_context_registry.csv']
for p in sources:shutil.copyfile(p,run/'code'/p.name)
command=[sys.executable,'-B',str(OUT/'scripts/cross_grid_audit.py')];r={'run_id':run.name,'stage':'supplementary_stability','command':command,'cwd':str(OUT),'python':sys.version,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_hashes':{str(p):sha(p) for p in sources},'input_hashes':{str(p):sha(p) for p in inputs},'new_fitting':False};(run/'before.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
p=subprocess.run(command,cwd=OUT,capture_output=True);(run/'stdout.txt').write_bytes(p.stdout);(run/'stderr.txt').write_bytes(p.stderr);r.update(exit_status=p.returncode,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources_unchanged=all(sha(p)==r['source_hashes'][str(p)] for p in sources),inputs_unchanged=all(sha(p)==r['input_hashes'][str(p)] for p in inputs));(run/'after.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8');print(p.stdout.decode(errors='replace'));print(p.stderr.decode(errors='replace'));raise SystemExit(p.returncode)
