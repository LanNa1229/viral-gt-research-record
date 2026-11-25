"""Run a local plotting/QC stage with contemporaneous code and launch records."""
import ast,datetime,hashlib,json,os,platform,shutil,subprocess,sys
from pathlib import Path
OUT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 stage=sys.argv[1];assert stage in ['prepare_figure_inputs.py','plot_integrated_umap.py','validate_figure.py']
 run=OUT/'executions'/(stage.removesuffix('.py')+'_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
 run.mkdir(parents=True,exist_ok=False)
 files=[OUT/stage,Path(__file__)]
 if (OUT/'config.json').exists():files.append(OUT/'config.json')
 for p in files:
  if p.suffix=='.py':ast.parse(p.read_text(encoding='utf-8'))
  shutil.copyfile(p,run/p.name)
 record=dict(command=[sys.executable,'-B',str(OUT/stage)],cwd=str(OUT),started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,executable=sys.executable,platform=platform.platform(),code_configuration_sha256={p.name:sha(p) for p in run.iterdir() if p.is_file()},input_manifest_sha256=sha(OUT/'INPUT_MANIFEST_SHA256.tsv') if (OUT/'INPUT_MANIFEST_SHA256.tsv').exists() else 'Created by preparation before joins/counts',model_fitting=False,Git_mutations=False)
 (run/'before.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
 env=dict(os.environ);env['PYTHONIOENCODING']='utf-8';env['MPLBACKEND']='Agg';env['MPLCONFIGDIR']=str(OUT/'runtime/mplconfig')
 p=subprocess.run(record['command'],cwd=OUT,env=env,capture_output=True,encoding='utf-8')
 (run/'stdout.txt').write_text(p.stdout,encoding='utf-8');(run/'stderr.txt').write_text(p.stderr,encoding='utf-8')
 (run/'after.json').write_text(json.dumps(dict(exit_code=p.returncode,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2)+'\n',encoding='utf-8')
 print(p.stdout,flush=True);print(p.stderr,flush=True)
 raise SystemExit(p.returncode)
if __name__=='__main__':main()
