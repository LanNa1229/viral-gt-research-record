"""Execute only local source curation, table construction and QC; never Git."""
import ast,datetime,hashlib,json,platform,shutil,subprocess,sys
from pathlib import Path
OUT=Path(__file__).resolve().parent
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main():
 run=OUT/'executions'/('run_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
 run.mkdir(parents=True,exist_ok=False);code=run/'executed_code';code.mkdir()
 scripts=['curation_common.py','curate_priority_records.py','build_reference_database.py','validate_reference_database.py','run_curation.py']
 for name in scripts:
  ast.parse((OUT/name).read_text(encoding='utf-8'));shutil.copyfile(OUT/name,code/name)
 version=dict(executable=sys.executable,python=sys.version,platform=platform.platform(),packages='Python standard library only',
  argv=sys.argv,working_directory=str(Path.cwd()),executed_code_sha256={n:sha(code/n) for n in scripts},
  no_model_dependencies=True,no_commits_or_pushes=True,started_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
 (run/'before.json').write_text(json.dumps(version,indent=2)+'\n',encoding='utf-8')
 (OUT/'SOFTWARE_VERSIONS.json').write_text(json.dumps(version,indent=2)+'\n',encoding='utf-8')
 log=(OUT/'EXECUTION_LOG.md').open('w',encoding='utf-8')
 log.write('# Executed curation pipeline\n\nRun evidence directory: '+str(run)+'\n\n');log.flush()
 for name in scripts[1:4]:
  cmd=[sys.executable,'-B',str(OUT/name)]
  print('RUN',name,flush=True)
  start=datetime.datetime.now(datetime.timezone.utc).isoformat()
  p=subprocess.run(cmd,cwd=OUT,capture_output=True,encoding='utf-8')
  (run/(name+'.stdout.txt')).write_text(p.stdout,encoding='utf-8')
  (run/(name+'.stderr.txt')).write_text(p.stderr,encoding='utf-8')
  record=dict(command=cmd,started_UTC=start,finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),exit_code=p.returncode,
   source_sha256=sha(OUT/name),stdout_path=name+'.stdout.txt',stderr_path=name+'.stderr.txt')
  (run/(name+'.execution.json')).write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
  log.write('## '+name+'\n\n'+json.dumps(record,indent=2)+'\n\n');log.flush()
  print(p.stdout[-16000:],flush=True)
  if p.returncode:
   print(p.stderr,flush=True);log.close();raise SystemExit(p.returncode)
 log.write('All stages completed successfully. No scientific model calculations, commits or pushes.\n');log.close()
 after=dict(status='PASS',finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
  outputs=[p.name for p in OUT.glob('*.tsv')],validation_sha256=sha(OUT/'validation_results.json'))
 (run/'after.json').write_text(json.dumps(after,indent=2)+'\n',encoding='utf-8')
 print('PASS local curation pipeline. Output-manifest finalization is a separate logged step.',flush=True)
if __name__=='__main__':
 sys.stdout.reconfigure(encoding='utf-8');main()
