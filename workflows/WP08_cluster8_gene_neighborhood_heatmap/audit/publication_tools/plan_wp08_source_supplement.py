"""Inventory existing source files explicitly referenced by the archived hit tables."""
import csv,hashlib,json
from pathlib import Path
BASE=Path('/mnt/d/Projects/Viral GT cluster 8');OUT=BASE/'WP08_github_publication_20261102'
PKG=BASE/'WP08_release_candidate/workflows/WP08_cluster8_gene_neighborhood_heatmap'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 records=[];loci=[];missing=[]
 for cand in csv.DictReader((PKG/'originals/260630_Temp/cluster8_analysis_manifest.fixed.tsv').open(),delimiter='\t'):
  c=cand['candidate_id'];raw=PKG/'originals/260630_Temp'/c/'final/foldseek_all_hits.tsv'
  with raw.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
  refs={}
  for r in rows:refs[r['query_id']]=(r['structure_file'],r['raw_json_file'])
  roots=set();structroots=set()
  for query,(structure,rawjson) in sorted(refs.items()):
   p=Path(rawjson)
   if not p.exists():
    old=p; matches=list(Path('/home/rahnn').glob('2606*Foldseek*/'+p.parent.parent.name+'/raw_results/'+p.name))
    if len(matches)==1:p=matches[0]
    else:missing.append(dict(candidate=c,query=query,kind='raw_json',reference=str(old),matches=[str(m) for m in matches]));continue
   root=p.parent.parent;batch=root.parent;roots.add(root);structroots.add(batch/'structures')
   struct=batch/structure
   if not struct.is_file():
    original_struct=Path(rawjson).parent.parent.parent/structure
    if original_struct.is_file():struct=original_struct;structroots.add(struct.parent)
   if not struct.is_file():missing.append(dict(candidate=c,query=query,kind='query_structure',reference=str(struct)))
   else:records.append(dict(candidate=c,kind='query_structure',query=query,source=str(struct),bytes=struct.stat().st_size))
   records.append(dict(candidate=c,kind='raw_api_json',query=query,source=str(p),historical_reference=rawjson,bytes=p.stat().st_size))
  # Existing per-query tables and execution ledger/logs belong to the exact run roots.
  extra=[]
  for root in roots:
   for p in root.rglob('*'):
    if not p.is_file() or ':Zone.Identifier' in p.name or p.suffix=='.pyc': continue
    # All result files from this referenced run, including historical final/annotation variants.
    extra.append(dict(candidate=c,kind='run_record',source=str(p),bytes=p.stat().st_size))
  records.extend(extra)
  for structroot in structroots:
   # Include successful and no-hit query models, not just rows appearing in hit tables.
   if structroot.is_dir() and (c!='B2ZYD6' or structroot.parent.name=='260617_Foldseek_batch'):
    records.extend(dict(candidate=c,kind='query_structure',source=str(p),bytes=p.stat().st_size) for p in structroot.rglob('*') if p.is_file() and ':Zone.Identifier' not in p.name)
  batchdirs={p.parent for p in structroots}|{p.parent for p in roots}
  for batch in batchdirs:
   for p in batch.iterdir():
    if p.is_file() and ':Zone.Identifier' not in p.name and (p.suffix in ['.py','.tsv','.bak'] or p.name=='how_to_run.txt'):
     records.append(dict(candidate=c,kind='source_script_or_manifest',source=str(p),bytes=p.stat().st_size))
  loci.append(dict(candidate=c,queries=len(refs),run_roots=[str(p) for p in roots],query_structure_roots=[str(p) for p in structroots],historical_raw_sha256=sha(raw),
    run_raw_tables=[dict(path=str(p),sha256=sha(p),matches_historical=sha(p)==sha(raw)) for root in roots for p in root.glob('final/foldseek_all_hits.tsv')]))
 unique={r['source']:r for r in records};records=list(unique.values())
 report=dict(files=records,loci=loci,missing=missing,bytes=sum(r['bytes'] for r in records),largest=sorted(records,key=lambda r:r['bytes'],reverse=True)[:10],over_100_MiB=[r for r in records if r['bytes']>=100*1024*1024])
 (OUT/'source_supplement_plan.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(dict(file_count=len(records),total_bytes=report['bytes'],missing=missing,largest=report['largest'],over_100_MiB=report['over_100_MiB'],loci=loci),indent=2))
if __name__=='__main__':main()
