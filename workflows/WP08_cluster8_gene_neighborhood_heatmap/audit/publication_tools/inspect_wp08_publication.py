"""Read publication inputs and git history; do not change scientific data."""
import csv,json,os,subprocess
from collections import Counter
from pathlib import Path
BASE=Path('/mnt/d/Projects/Viral GT cluster 8')
REPO=Path('/home/rahnn/projects/viral-gt-research-record')
PKG=BASE/'WP08_release_candidate/workflows/WP08_cluster8_gene_neighborhood_heatmap'
OUT=BASE/'WP08_github_publication_20261102'
def run(args,cwd=REPO):
 r=subprocess.run(args,cwd=cwd,capture_output=True,text=True,env=dict(os.environ,GIT_OPTIONAL_LOCKS='0'))
 return dict(command=args,returncode=r.returncode,stdout=r.stdout,stderr=r.stderr)
def main():
 OUT.mkdir(exist_ok=True)
 commands=[['git','branch','--show-current'],['git','status','--porcelain=v1'],['git','log','-5','--oneline'],['git','ls-remote','origin','refs/heads/main'],['git','diff','--cached','--name-status'],['git','config','--get','credential.helper']]
 history=[run(a) for a in commands]
 print(json.dumps(history,indent=2))
 (OUT/'git_history_before.json').write_text(json.dumps(history,indent=2))
 audit=[]
 for rec in csv.DictReader((PKG/'originals/260630_Temp/cluster8_analysis_manifest.fixed.tsv').open(),delimiter='\t'):
  raw=PKG/'originals/260630_Temp'/rec['candidate_id']/'final/foldseek_all_hits.tsv'
  with raw.open() as f:
   records=list(csv.DictReader(f,delimiter='\t'))
  fields={n:sorted(set(r.get(n,'') for r in records)) for n in ['run_id','structure_file','raw_json_file','submitted_at','databases_submitted']}
  audit.append(dict(candidate=rec['candidate_id'],raw=str(raw),fields=fields))
 (OUT/'raw_input_references.json').write_text(json.dumps(audit,indent=2))
 for r in audit:
  print('CANDIDATE',r['candidate'],'RUNS',r['fields']['run_id'],'STRUCTURE',r['fields']['structure_file'][:2],'RAW',r['fields']['raw_json_file'][:2])
 directories=[]
 for p in sorted(Path('/home/rahnn').glob('2606*Foldseek*')):
  counts=Counter();sizes=Counter();total=0
  for q in p.rglob('*'):
   if q.is_file():
    try: st=q.stat()
    except FileNotFoundError:continue
    counts[q.suffix]+=1;sizes[q.suffix]+=st.st_size;total+=st.st_size
  directories.append(dict(path=str(p),bytes=total,extensions=dict(counts),extension_bytes=dict(sizes),top_names=[q.name for q in p.iterdir()]))
 (OUT/'related_directory_inventory.json').write_text(json.dumps(directories,indent=2))
 print('DIRECTORIES',json.dumps(directories,indent=2))
 for name in ['LICENSE-MIT.txt','LICENSE-CC-BY-4.0.md','THIRD_PARTY_NOTICES.md']:
  p=REPO/name
  if p.exists():print('NOTICE',name,p.read_text())
if __name__=='__main__': main()
