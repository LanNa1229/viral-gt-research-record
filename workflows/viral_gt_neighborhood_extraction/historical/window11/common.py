"""Shared paths and frozen source imports for the paired length experiment."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parents[1]
PRIMARY=ROOT/'analysis/viral_gt_mining_20261017'
BASE=PRIMARY/'runs/baseline_594_v2'
RUN=PRIMARY/'runs/window11_594_v1'
CTRL=PRIMARY/'runs/controls11_existing_v1'
PRIOR=ROOT/'analysis/viral_gt_interpretation_20261017'
RESCUE=ROOT/'analysis/viral_gt_retrieval_rescue_20261017'
REF=ROOT/'analysis/viral_glycan_benchmark_20261016'
sys.path.insert(0,str(BASE/'implementation_snapshot/src'))
from viral_gt_mining import extract as extraction,inference,analysis as baseline_analysis
from viral_gt_mining.extract import read_tsv as read,write_tsv as write,json_write as save,sha_file as sha
import numpy as np
import json,collections,hashlib,re
SEED=20261017
LABELS=inference.LABELS
GT=re.compile(r'glycosyl\s*transferase|glycosyltransferase|galactosyltransferase|glucosyltransferase|mannosyltransferase|fucosyltransferase|rhamnosyltransferase|xylosyltransferase|sialyltransferase|glucuronosyltransferase',re.I)
CARB=re.compile(r'glycosyl|glycan|polysacchar|saccharide|oligosacchar|carbohydrat|sugar|glucos|galactos|mannos|fucos|rhamnos|xylos|sialic|neuramin|glucuron|glucan|chitin|cellulos|dTDP|UDP[- ]|GDP[- ]mann|lipopolysacchar|capsular',re.I)
UNKNOWN=re.compile(r'hypothetical|uncharacterized|unknown|\bDUF\d*\b|^ORF\s*\d*$',re.I)
def yes(x):return str(x).lower() in {'true','1','yes'}
def verify_prior():
    counts={}
    for directory,manifest in [(BASE,PRIOR/'freeze/baseline_files.tsv'),(PRIOR,PRIOR/'artifact_checksums.tsv'),(RESCUE,RESCUE/'artifact_checksums.tsv')]:
        rows=read(manifest)
        for r in rows:
            path=directory/r['path'].replace('\\','/')
            assert sha(path)==r['sha256'],str(path)
        counts[directory.name]=len(rows)
    return counts
def grouped_genes(directory):
    result=collections.defaultdict(list)
    for g in read(directory/'viral_gt_window_genes.tsv'):result[g['window_id']].append(g)
    for v in result.values():v.sort(key=lambda g:int(g['relative_position']))
    return result
def props(genes,gt_hashes):
    return dict(n_genes=len(genes),annotated_gt_count=sum(bool(GT.search(g.get('product',''))) for g in genes),
        gt_proxy_count=sum(bool(GT.search(g.get('product',''))) or g['sequence_sha256'] in gt_hashes for g in genes),
        carbohydrate_annotation_count=sum(bool(CARB.search(g.get('product',''))) for g in genes),
        unknown_annotation_count=sum(not g.get('product','').strip() or bool(UNKNOWN.search(g['product'])) for g in genes),
        total_aa=sum(len(g['sequence']) for g in genes),n_long=sum(len(g['sequence'])>1022 for g in genes),
        n_reverse=sum(int(g['strand'])==-1 for g in genes))
def status(ids,pred):
    if any(pred.get(w,False) for w in ids):return 'candidate'
    return 'not_selected' if all(w in pred for w in ids) else 'unresolved'
def bootstrap(values,seed=SEED,n=5000):
    values=np.array(values,float);rng=np.random.default_rng(seed)
    if not len(values):return dict(n=0,mean=None,low=None,high=None)
    # Genome means are the explicit resampling units; vectorized in small batches.
    samples=np.concatenate([values[rng.integers(len(values),size=(min(250,n-i),len(values)))].mean(1) for i in range(0,n,250)])
    return dict(n=len(values),mean=float(values.mean()),low=float(np.quantile(samples,.025)),high=float(np.quantile(samples,.975)))
def load_results(directory):
    p={r['window_id']:r for r in read(directory/'viral_bgc_predictions.tsv')}
    loc={r['window_id']:r for r in read(directory/'saccharide_localization.tsv')}
    m=read(directory/'embedding_manifest.tsv');v=np.load(directory/'viral_bgc_embeddings.npy',allow_pickle=False)
    assert len(p)==len(loc)==len(m)==len(v)
    for i,r in enumerate(m):assert int(r['embedding_row'])==i and r['window_id'] in p and r['window_id'] in loc
    return p,loc,{r['window_id']:v[i] for i,r in enumerate(m)}
