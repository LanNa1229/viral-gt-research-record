"""Record standard-library documentation generation; no fitting or plotting."""
import datetime,hashlib,json,platform,shutil,subprocess,sys
from pathlib import Path
OUT=Path(__file__).resolve().parents[1];run=OUT/'runs'/f'run_{len(list((OUT/"runs").glob("run_*")))+1:03d}_documentation';run.mkdir();(run/'code').mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
sources=[OUT/'scripts/write_records_final.py',OUT/'scripts/common.py',OUT/'scripts/record_documentation.py',OUT/'config.json']
for p in sources:shutil.copyfile(p,run/'code'/p.name)
inputs=[OUT/'RESULTS.json',OUT/'VALIDATION.json']+list((OUT/'tables').glob('*.csv'))+list((OUT/'stability').glob('*.csv'))+list((OUT/'runs').glob('*/before.json'))+list((OUT/'runs').glob('*/after.json'))
command=[sys.executable,'-B',str(OUT/'scripts/write_records_final.py')];record={'run_id':run.name,'stage':'documentation','command':command,'cwd':str(OUT),'python':sys.version,'platform':platform.platform(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_hashes':{str(p):sha(p) for p in sources},'input_hashes':{str(p):sha(p) for p in inputs},'scope':'Author Methods/readout from recorded results; standard-library file I/O and arithmetic only. No fitting or annotation.'}
(run/'before.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
r=subprocess.run(command,cwd=OUT,capture_output=True);(run/'stdout.txt').write_bytes(r.stdout);(run/'stderr.txt').write_bytes(r.stderr);record.update(exit_status=r.returncode,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_unchanged=all(sha(p)==record['source_hashes'][str(p)] for p in sources),inputs_unchanged=all(sha(p)==record['input_hashes'][str(p)] for p in inputs));(run/'after.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');print(r.stdout.decode(errors='replace'));print(r.stderr.decode(errors='replace'));print('DOCUMENTATION',r.returncode,record['source_unchanged'],record['inputs_unchanged']);raise SystemExit(r.returncode)
