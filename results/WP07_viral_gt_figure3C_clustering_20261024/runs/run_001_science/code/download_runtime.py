import datetime,hashlib,json,subprocess,sys
from pathlib import Path
OUT=Path(__file__).resolve().parents[1];log=OUT/'logs/download_recorded';log.mkdir(exist_ok=True)
command=[sys.executable,'-m','pip','download','--only-binary=:all:','--no-deps','--disable-pip-version-check','--dest',str(OUT/'runtime/wheels'),'igraph==1.0.0','leidenalg==0.12.0','texttable==1.7.2']
record={'command':command,'cwd':str(OUT),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(log/'before.json').write_text(json.dumps(record,indent=2)+'\n')
r=subprocess.run(command,cwd=OUT,capture_output=True);(log/'stdout.txt').write_bytes(r.stdout);(log/'stderr.txt').write_bytes(r.stderr);record.update(exit_status=r.returncode,files=[p.name for p in (OUT/'runtime/wheels').iterdir()]);(log/'after.json').write_text(json.dumps(record,indent=2)+'\n');print(r.stdout.decode(errors='replace'));print(r.stderr.decode(errors='replace'));print(record['files']);raise SystemExit(r.returncode)
