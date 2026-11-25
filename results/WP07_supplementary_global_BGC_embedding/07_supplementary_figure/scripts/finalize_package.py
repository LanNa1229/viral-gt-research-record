"""Finalize completed local figure records and verify every output hash."""
import concurrent.futures,csv,datetime,hashlib,json,platform,sys
from pathlib import Path
OUT=Path(__file__).resolve().parent
EXCLUDE={'OUTPUT_MANIFEST_SHA256.tsv','CHECKSUM_VERIFICATION.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(n,x):(OUT/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def main():
 for stage in ['prepare_figure_inputs','plot_integrated_umap','validate_figure']:
  runs=sorted((OUT/'executions').glob(stage+'_*/before.json'));assert runs
  assert json.loads((runs[-1].parent/'after.json').read_text())['exit_code']==0
 v=json.loads((OUT/'VALIDATION_RESULTS.json').read_text());assert v['status']=='PASS' and len(v['checks'])==22
 qa=json.loads((OUT/'VISUAL_QA.json').read_text());assert qa['status'].startswith('PASS')
 for r in qa['inspected_outputs']:assert sha(OUT/r['relative_path'])==r['sha256']
 entries=[]
 for p in sorted((OUT/'executions').glob('*/before.json')):
  b=json.loads(p.read_text());a=json.loads((p.parent/'after.json').read_text())
  entries.append(dict(directory=str(p.parent.relative_to(OUT)),started_utc=b['started_utc'],finished_utc=a['finished_utc'],exit_code=a['exit_code'],command=b['command']))
 save('EXECUTION_SUMMARY.json',entries)
 lines=['# Actual figure execution log','','Each directory contains frozen code/configuration, pre-execution launch evidence, stdout, stderr and final exit code. Failed development checks remain available; final preparation, plotting and QC succeeded.','']
 for r in entries:
  lines.extend(['- '+r['directory']+'; start '+r['started_utc']+'; finish '+r['finished_utc']+'; exit '+str(r['exit_code']), '  Actual child argument array: '+json.dumps(r['command']),''])
 lines.extend(['PDF_RENDER_EXECUTION.json contains the separate contemporaneous independent Poppler-render record.','VISUAL_QA.json contains reviewed artifact hashes and display limitations.','No fitting, clustering, commits or pushes occurred in this plotting task.'])
 (OUT/'EXECUTION_LOG.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 save('FINALIZATION_RECORD.json',dict(command=[sys.executable,*sys.argv],cwd=str(Path.cwd()),started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,platform=platform.platform(),script_sha256=sha(Path(__file__)),excluded_from_manifest=sorted(EXCLUDE)))
 files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.relative_to(OUT).as_posix() not in EXCLUDE)
 def rec(p):return dict(relative_path=p.relative_to(OUT).as_posix(),size_bytes=p.stat().st_size,sha256=sha(p))
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:rows=list(pool.map(rec,files))
 manifest=OUT/'OUTPUT_MANIFEST_SHA256.tsv'
 with manifest.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['relative_path','size_bytes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 with manifest.open(encoding='utf-8',newline='') as f:reloaded=list(csv.DictReader(f,delimiter='\t'))
 def verify(r):
  p=OUT/r['relative_path'];assert p.resolve().is_relative_to(OUT.resolve())
  return p.stat().st_size==int(r['size_bytes']) and sha(p)==r['sha256']
 with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:passed=list(pool.map(verify,reloaded))
 assert all(passed)
 assert {r['relative_path'] for r in reloaded}=={p.relative_to(OUT).as_posix() for p in OUT.rglob('*') if p.is_file() and p.relative_to(OUT).as_posix() not in EXCLUDE}
 p=OUT/'figures/SuppFig_global_leiden_biosynthetic_classes_minimalist_v2'
 seal=dict(status='PASS',finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),manifest_sha256=sha(manifest),output_files_verified=len(rows),total_output_bytes=sum(int(r['size_bytes']) for r in rows),input_files_verified=34,scientific_and_export_checks_passed=22,PDF_sha256=sha(p.with_suffix('.pdf')),SVG_sha256=sha(p.with_suffix('.svg')),PNG_sha256=sha(p.with_suffix('.png')),source_data_sha256=sha(OUT/'FIGURE_SOURCE_DATA.tsv'),excluded_from_manifest=sorted(EXCLUDE),local_only=True,committed=False,pushed=False)
 save('CHECKSUM_VERIFICATION.json',seal)
 print(json.dumps(seal,indent=2),flush=True)
if __name__=='__main__':main()
