"""Install only recorded local Leiden wheels; retain setup stdout/stderr and exit code."""
import datetime,hashlib,json,subprocess,sys
from pathlib import Path
OUT=Path(__file__).resolve().parents[1];w=Path(r'C:\Users\rahnn\.codex\visualizations\2026\10\24\01a16c5b-2688-7373-bfe8-092a74217b20\clustering_runtime\wheels');log=OUT/'logs/setup_install_04';log.mkdir(exist_ok=True)
command=[sys.executable,'-m','pip','install','--no-index','--no-deps','--no-cache-dir','--disable-pip-version-check','--target',str(w.parent/'deps'),'--report',str(log/'pip_report.json')]+[str(p) for p in sorted(w.glob('*.whl'))]
record={'command':command,'cwd':str(OUT),'python':sys.version,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'wheel_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in w.glob('*.whl')}}
(log/'before.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
r=subprocess.run(command,cwd=OUT,capture_output=True);(log/'stdout.txt').write_bytes(r.stdout);(log/'stderr.txt').write_bytes(r.stderr)
record.update(exit_status=r.returncode,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(log/'after.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
print(r.stdout.decode(errors='replace'));print(r.stderr.decode(errors='replace'));raise SystemExit(r.returncode)


