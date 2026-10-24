"""Prepare a NEW reproduction workspace from archived inputs. Never runs analysis."""
import argparse,csv,hashlib,json,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workspace-root',type=Path,required=True);p.add_argument('--previous-archive',type=Path,required=True);p.add_argument('--graph-runtime',type=Path,required=True);p.add_argument('--analysis-dependencies',default='analysis/bgcprophet_reference/.deps');a=p.parse_args()
assert not a.workspace_root.exists(),'Choose a new nonexistent workspace; never overwrite originals.'
assert (a.previous_archive/'original').is_dir() and a.graph_runtime.is_dir()
a.workspace_root.mkdir(parents=True)
old=a.workspace_root/'analysis/viral_gt_umap_20261024';shutil.copytree(a.previous_archive/'original',old)
with (HERE/'ARCHIVED_INPUT_MANIFEST.tsv').open(encoding='utf-8',newline='') as f:inputs=list(csv.DictReader(f,delimiter='\t'))
for row in inputs:
 dest=a.workspace_root/row['workspace_relative_path']
 if row['archive_path']:
  src=HERE/row['archive_path'];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest)
 assert dest.exists() and hashlib.sha256(dest.read_bytes()).hexdigest()==row['sha256'],row['workspace_relative_path']
out=a.workspace_root/'analysis/viral_gt_clustering_20261024';out.mkdir();shutil.copytree(HERE/'scripts',out/'scripts')
cfg=json.loads((HERE/'config.json').read_text(encoding='utf-8-sig'));cfg.update(workspace_root=str(a.workspace_root.resolve()),graph_runtime=str(a.graph_runtime.resolve()),analysis_dependencies=a.analysis_dependencies);(out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n',encoding='utf-8')
for row in inputs:row['path']=str(a.workspace_root/row['workspace_relative_path'])
fields=['path','workspace_relative_path','role','size_bytes','sha256']
with (out/'INPUT_MANIFEST.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fields,lineterminator='\n');w.writeheader();w.writerows({k:row[k] for k in fields} for row in inputs)
print('Restored and hash-verified',len(inputs),'input/protected identities. No scientific code executed.')
print('New output directory:',out)
print('No previous runs/results copied. A later fit will receive fresh code/configuration/input/runtime/exit records.')
