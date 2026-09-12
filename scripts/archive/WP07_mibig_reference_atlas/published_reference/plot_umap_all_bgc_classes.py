"""Add class colors to saved UMAP coordinates; never fits or imports UMAP.
Run: python -B scripts/plot_umap_all_bgc_classes.py
Uses project-local analysis dependencies. Existing outputs are never overwritten.
"""
import os,sys,csv,json,hashlib
from pathlib import Path
from collections import Counter
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/".deps"))
os.environ["MPLCONFIGDIR"]=str(ROOT/".cache/matplotlib")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=ROOT/"results"
CANONICAL=("Alkaloid","Terpene","NRP","Polyketide","RiPP","Saccharide","Other")
CLASSES=CANONICAL+("Hybrid",)
COLORS=("#CC79A7","#009E73","#0072B2","#D55E00","#E69F00","#56B4E9","#777777","#7B3294")
def read(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()
targets=[R/"umap_all_bgc_classes.tsv",R/"figures/umap_all_bgc_classes.png",
         R/"figures/umap_all_bgc_classes.pdf",R/"umap_all_bgc_classes_provenance.json"]
if any(p.exists() for p in targets):raise FileExistsError("New outputs already exist; refusing to overwrite.")
before={p:sha(p) for p in R.rglob("*") if p.is_file()}
coords=read(R/"umap_coordinates.tsv");meta=read(R/"consensus_metadata.tsv")
assert len(coords)==len(meta)==2502
lookup={r["BGC"]:r for r in meta}
assert len(lookup)==2502 and len({r["BGC"] for r in coords})==2502
assert set(lookup)=={r["BGC"] for r in coords}
rows=[]
for r in coords:
    label=lookup[r["BGC"]]["label"]
    tokens=label.split()
    assert tokens and set(tokens)<=set(CANONICAL),label
    assert len(tokens)==len(set(tokens)),label
    cls=tokens[0] if len(tokens)==1 else "Hybrid"
    rows.append(dict(BGC=r["BGC"],original_label=label,plot_class=cls,
                     UMAP001=r["UMAP001"],UMAP002=r["UMAP002"]))
# The original UMAP result was float32. Casting saved round-trip decimals back
# to float32 recovers the exact coordinates used by the original figure.
xy=np.array([[r["UMAP001"],r["UMAP002"]] for r in rows],dtype=np.float32)
assert xy.shape==(2502,2) and np.isfinite(xy).all()
counts=Counter(r["plot_class"] for r in rows)
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "pdf.fonttype":42,"ps.fonttype":42,"savefig.dpi":600})
# Reproduce the original umap_saccharide scatter autoscaling, not its UMAP fit.
# Both panels include all points; verify that their limits agree.
limits=[]
for pure in (False,True):
    positive=np.array([(r["original_label"].strip()=="Saccharide") if pure else
                       ("Saccharide" in r["original_label"]) for r in rows])
    ref,ax=plt.subplots()
    ax.scatter(xy[~positive,0],xy[~positive,1],s=6,alpha=.5)
    ax.scatter(xy[positive,0],xy[positive,1],s=13,alpha=.8)
    limits.append((ax.get_xlim(),ax.get_ylim()))
    plt.close(ref)
assert limits[0]==limits[1]
fig,ax=plt.subplots(figsize=(6.8,4.2))
classes=np.array([r["plot_class"] for r in rows])
for cls,color in zip(CLASSES,COLORS):
    mask=classes==cls
    ax.scatter(xy[mask,0],xy[mask,1],s=10,c=color,alpha=.8,
               linewidths=0,rasterized=True,label=f"{cls} (n={counts[cls]:,})")
ax.set(xlim=limits[0][0],ylim=limits[0][1],xlabel="UMAP 1",ylabel="UMAP 2",
       title="Biosynthetic classes of 2,502 consensus BGCs")
ax.legend(loc="center left",bbox_to_anchor=(1.02,.5),frameon=False,
          title="Biosynthetic class",markerscale=1.6,fontsize=9)
fig.tight_layout()
with targets[0].open("x",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
for p in targets[1:3]:fig.savefig(p,bbox_inches="tight",dpi=600)
plt.close(fig)
assert all(p.exists() and sha(p)==h for p,h in before.items()),"An existing result changed"
provenance=dict(script_sha256=sha(Path(__file__)),command=sys.argv,
    numpy_version=np.__version__,matplotlib_version=matplotlib.__version__,
    coordinate_source="results/umap_coordinates.tsv",
    coordinate_sha256=before[R/"umap_coordinates.tsv"],
    metadata_sha256=before[R/"consensus_metadata.tsv"],
    canonical_classes=CANONICAL,hybrid_rule="More than one whitespace-separated canonical label",
    counts={c:counts[c] for c in CLASSES},xlim=limits[0][0],ylim=limits[0][1],
    limits_method="Original float32 coordinates, original scatter groups and Matplotlib autoscaling",
    umap_recomputed=False,existing_results_verified_unchanged=len(before),
    colors=dict(zip(CLASSES,COLORS)),figure_inches=[6.8,4.2],dpi=600)
with targets[3].open("x",encoding="utf-8") as f:json.dump(provenance,f,indent=2)
print(json.dumps(provenance,indent=2))

