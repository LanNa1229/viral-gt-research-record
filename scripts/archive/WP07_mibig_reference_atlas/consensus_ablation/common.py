"""Shared utilities for corrected consensus Step 3; never imports/runs Step 2."""
from pathlib import Path
import os,sys,json,hashlib,time,platform,importlib.metadata
ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT.parents[1]
PREV=PROJECT/"analysis/core_context_representation_geometry_20260910"
sys.path.insert(0,str(PREV/".deps"))
sys.dont_write_bytecode=True
for k in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"]:os.environ[k]="1"
os.environ["MPLCONFIGDIR"]=str(ROOT/".cache/matplotlib")
os.environ["NUMBA_CACHE_DIR"]=str(ROOT/".cache/numba")
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":10,"pdf.fonttype":42,"axes.spines.top":False,"axes.spines.right":False})
SEED=1729;N_PERM=10000;N_BOOT=10000;KS=[5,10,20,50]
DEFS=["any_saccharide","pure_saccharide"]
def tab(path,data):
 p=ROOT/path;p.parent.mkdir(parents=True,exist_ok=True)
 pd.DataFrame(data).to_csv(p,sep="\t",index=False,float_format="%.17g")
def js(path,data): (ROOT/path).write_text(json.dumps(data,indent=2,default=str)+"\n")
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()
def read(p):return pd.read_csv(p,sep="\t")
def load():
 m=read(ROOT/"results/consensus_metadata.tsv")
 a={r:np.load(ROOT/"results"/f"{n}_consensus_embeddings.npy") for r,n in [("A","published"),("B","core_context")]}
 return m,a
def save(fig,name,source):
 fig.savefig(ROOT/"figures"/(name+".png"),dpi=400,bbox_inches="tight")
 fig.savefig(ROOT/"figures"/(name+".pdf"),bbox_inches="tight")
 plt.close(fig)
 with open(ROOT/"provenance/figure_sources.tsv","a") as f:f.write(name+"\t"+source+"\n")
def versions():
 return {p:importlib.metadata.version(p) for p in ["numpy","pandas","scipy","scikit-learn","matplotlib","umap-learn","numba","pynndescent"]}
def integrity():
 before=read(ROOT/"provenance/input_sha256_before.tsv")
 out=[]
 for p,h in zip(before.path,before.sha256):
  actual=sha(p);out.append(dict(path=p,sha256_before=h,sha256_after=actual,unchanged=actual==h))
 tab("provenance/input_integrity_after.tsv",out)
 assert all(x["unchanged"] for x in out)
 return len(out)
