import csv,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis/bgcprophet_reference/.deps'))
import numpy as np
ARCH=Path(r'\\wsl.localhost\Ubuntu\home\rahnn\projects\viral-gt-research-record\results\WP07_viral_gt_figure3C_umap_20261024');OLD=ROOT/'analysis/viral_gt_umap_20261024'
def read(p,sep=','):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=sep))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest={r['path'].replace('\\','/'):r['sha256'] for r in read(ARCH/'original/output_checksums.csv')}
files=['tables/representative_15gene_embeddings.npy','tables/included_locus_representatives.csv','joint_umap_samples.csv','tables/locus_registry.csv','tables/window_locus_anchor_crosswalk.csv','tables/representative_neighbors_top20.csv']
checks=[]
for name in files:
 assert sha(OLD/name)==manifest[name]==sha(ARCH/'original'/name),name
 checks.append({'path':str(OLD/name),'sha256':manifest[name],'archive_evidence':'original/output_checksums.csv '+name})
x=np.load(OLD/files[0]);meta=read(OLD/files[1]);sel=[i for i,r in enumerate(meta) if r['candidate_status']=='candidate']
assert x.shape==(556,320) and len(meta)==556
if len(sel)!=287:raise SystemExit('STOP: selected cohort discrepancy: '+str(len(sel)))
assert len({r['viral_locus_id'] for r in meta})==556
samples=read(OLD/'joint_umap_samples.csv');v=[r for r in samples if r['sample_type']=='viral_locus']
assert [r['viral_locus_id'] for r in meta]==[r['viral_locus_id'] for r in v]
assert [r['representative_window_id'] for r in meta]==[r['representative_window_id'] for r in v]
refpath=ROOT/'analysis/viral_glycan_benchmark_20261016/mibig_CORE_ONLY_embeddings.npy';metapath=refpath.with_name('mibig_matched_reference_metadata.tsv')
im={r['workspace_relative_path']:r for r in read(ARCH/'INPUT_MANIFEST.tsv','\t')}
for p in [refpath,metapath]:
 key=p.relative_to(ROOT).as_posix();assert sha(p)==im[key]['execution_sha256'];checks.append({'path':str(p),'sha256':sha(p),'archive_evidence':'INPUT_MANIFEST.tsv '+key})
ref=np.load(refpath);rm=read(metapath,'\t');sac=sum('Saccharide' in r['label'].split() for r in rm)
assert ref.shape==(2502,320) and sac==195
result={'status':'PASS','primary_selected_shape':[len(sel),320],'all_viral_shape':list(x.shape),'reference_shape':list(ref.shape),'saccharide_references':sac,'global_joint_count':len(sel)+len(ref),'saccharide_joint_count':len(sel)+sac,'selected_full15_count':sum(meta[i]['representative_n_genes']=='15' for i in sel),'verified_files':checks,'row_alignment':'Saved representative metadata matches original viral UMAP rows and IDs; scientific runner will also verify representative rows against frozen source vectors.'}
(OUT/'PREFLIGHT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))
