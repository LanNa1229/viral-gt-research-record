"""Non-core contextual mean: saved-vector algebra only. No model code imported."""
from pathlib import Path
import os,sys,json,hashlib,time,platform,importlib.metadata,itertools
ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parents[1]
OLD=PROJECT/"analysis/core_context_representation_geometry_20260910"
CONS=PROJECT/"analysis/core_context_consensus_step3_20260911"
sys.path.insert(0,str(OLD/".deps"));sys.dont_write_bytecode=True
for key in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"]:os.environ[key]="1"
os.environ["MPLCONFIGDIR"]=str(ROOT/".cache/matplotlib")
import numpy as np,pandas as pd
from scipy import stats
SEED=1729;BOOT=10000;KS=[5,10,20,50];DEFS=["any_saccharide","pure_saccharide"]
def read(p):return pd.read_csv(p,sep="\t")
def tab(p,data):pd.DataFrame(data).to_csv(ROOT/p,sep="\t",index=False,float_format="%.17g")
def js(p,data):(ROOT/p).write_text(json.dumps(data,indent=2,default=str)+"\n")
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()
def integrity():
 before=read(ROOT/"provenance/input_sha256_before.tsv");out=[]
 for r in before.itertuples():
  after=sha(r.path);out.append(dict(path=r.path,sha256_before=r.sha256,sha256_after=after,unchanged=after==r.sha256))
 tab("provenance/input_integrity_after.tsv",out)
 assert all(x["unchanged"] for x in out),"Protected input modified"
 return len(out)
def describe(v):
 v=np.asarray(v);v=v[np.isfinite(v)]
 return dict(n=len(v),mean=np.mean(v),median=np.median(v),sd=np.std(v,ddof=1),minimum=np.min(v),maximum=np.max(v),
  q25=np.quantile(v,.25),q75=np.quantile(v,.75),p05=np.quantile(v,.05),p95=np.quantile(v,.95))
def raw():
 I=json.loads((ROOT/"provenance/inputs.json").read_text())
 return {r:np.load(I[k]).astype(np.float64) for r,k in [("A","published"),("B","core"),("N","noncore")]}
def cosine(a,b):
 den=np.linalg.norm(a,axis=-1)*np.linalg.norm(b,axis=-1)
 return np.divide(np.sum(a*b,axis=-1),den,out=np.full(den.shape,np.nan),where=den>0)
def bootdiff(d,seed):
 rng=np.random.default_rng(seed);means=[];medians=[]
 for start in range(0,BOOT,100):
  b=d[rng.integers(len(d),size=(min(100,BOOT-start),len(d)))]
  means.extend(b.mean(1));medians.extend(np.median(b,axis=1))
 mc=np.quantile(means,[.025,.975]);dc=np.quantile(medians,[.025,.975])
 return dict(mean_difference=np.mean(d),median_difference=np.median(d),mean_ci95_low=mc[0],mean_ci95_high=mc[1],
  median_ci95_low=dc[0],median_ci95_high=dc[1],fraction_improved=np.mean(d>1e-12),fraction_equal=np.mean(abs(d)<=1e-12),fraction_worsened=np.mean(d< -1e-12))
