"""Restore archived files to a NEW workspace. Does not fit UMAP or run analysis."""
import argparse,csv,hashlib,shutil
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--workspace',type=Path,required=True);p.add_argument('--code-stage',choices=['fit','render'],default='fit');a=p.parse_args()
assert not a.workspace.exists(),'Choose a new, nonexistent directory; originals must not be overwritten.'
a.workspace.mkdir(parents=True)
with (HERE/'INPUT_MANIFEST.tsv').open(encoding='utf-8',newline='') as f:inputs=list(csv.DictReader(f,delimiter='\t'))
missing=[]
for r in inputs:
 dest=a.workspace/r['workspace_relative_path']
 if r['availability']!='ARCHIVED':missing.append(r);continue
 src=HERE/r['archive_path'];assert hashlib.sha256(src.read_bytes()).hexdigest()==r['execution_sha256']
 dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest)
out=a.workspace/'analysis/viral_gt_umap_20261024'
shutil.copytree(HERE/'original',out)
if a.code_stage=='fit':shutil.copyfile(HERE/'executed_code/exploratory_umap.fit_executed.py',out/'scripts/exploratory_umap.py')
with (a.workspace/'EXTERNAL_INPUTS_TO_ACQUIRE.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,list(inputs[0]),delimiter='\t');w.writeheader();w.writerows(missing)
print('Restored files only; no scientific code executed. Acquire and SHA256-verify',len(missing),'external files. See EXTERNAL_INPUTS_TO_ACQUIRE.tsv.')
print('Recreate the recorded Python/package environment separately; no historical wheel/binary lock is available.')
