from common import *
import re
def main():
 t=time.time()
 I=json.loads((PREV/"provenance/inputs.json").read_text())
 protected=[]
 for d in [PROJECT/"analysis/bgcprophet_reference",PREV,Path(I["canonical_directory"])]:
  for parent,dirs,files in os.walk(d):
   dirs[:]=[v for v in dirs if v not in [".deps",".cache",".tmp","__pycache__",".git"]]
   protected.extend(Path(parent)/v for v in files)
 protected.append(Path(I["csv"]))
 protected=sorted(set(protected))
 tab("provenance/input_sha256_before.tsv",[dict(path=str(p),sha256=sha(p)) for p in protected])
 js("provenance/inputs.json",I)
 proof=read(Path(I["canonical_directory"])/"input_integrity_before.tsv")
 assert sha(I["published"]) in proof.astype(str).values
 summary=json.loads((Path(I["canonical_directory"])/"full_summary.json").read_text())
 assert summary["passed"] and summary["unmasked_reference_exact_rows"]==12510 and summary["unmasked_reference_max_abs_difference"]==0
 A=np.load(I["published"]);B=np.load(I["core"])
 am=read(I["published_metadata"]);bm=read(I["core_metadata"])
 assert am.equals(bm) and len(am)==12510 and am.ID.is_unique
 assert np.array_equal(am.embedding_row,np.arange(12510))
 csv=pd.read_csv(I["csv"])
 assert np.array_equal(csv.ID,am.ID) and np.array_equal(csv.labels,am.labels)
 for a in [A,B]:assert a.shape==(12510,320) and np.isfinite(a).all() and a.dtype==np.float32
 parts=am.ID.str.extract(r"^(BGC[0-9]{7})_(-1|0|1|2|3)$")
 assert not parts.isna().any().any()
 am["base_BGC"]=parts[0];am["suffix"]=parts[1].astype(int)
 classes=json.loads(Path(I["palette"]).read_text())["colors"]
 rows=[];groups=[]
 for bg,g in am.groupby("base_BGC",sort=True):
  g=g.sort_values("suffix")
  assert len(g)==5 and set(g.suffix)=={-1,0,1,2,3} and g.labels.nunique()==1
  label=g.labels.iloc[0].strip();tokens=label.split()
  assert tokens and len(tokens)==len(set(tokens)) and set(tokens)<=set(classes)-{"Hybrid"}
  ix=g.embedding_row.to_numpy();groups.append(ix)
  rows.append(dict(base_BGC=bg,label=label,plot_class=tokens[0] if len(tokens)==1 else "Hybrid",
    any_saccharide="Saccharide" in tokens,pure_saccharide=label=="Saccharide",
    source_row_IDs=json.dumps(g.ID.tolist()),source_row_indices=json.dumps(ix.tolist()),
    augmentation_suffixes=json.dumps(g.suffix.tolist()),n_augmentations=5))
 assert len(rows)==2502
 ix=np.array(groups)
 for name,array in [("published",A),("core_context",B)]:
  consensus=array[ix].mean(axis=1,dtype=np.float64)
  assert consensus.shape==(2502,320) and consensus.dtype==np.float64 and np.isfinite(consensus).all()
  assert np.all(np.linalg.norm(consensus,axis=1)>0)
  np.save(ROOT/"results"/f"{name}_consensus_embeddings.npy",consensus)
 tab("results/consensus_metadata.tsv",rows)
 tab("tables/class_counts.tsv",pd.DataFrame(rows).groupby("plot_class").size().rename("n_BGC").reset_index())
 js("provenance/parameters.json",dict(seed=SEED,bootstrap_count=N_BOOT,permutation_count=N_PERM,ks=KS,
  consensus="Raw arithmetic mean of five source vectors, float64 accumulation and float64 storage; suffix order -1,0,1,2,3",
  knn="Original 320-D cosine; float64 transient L2-normalized copy; diagonal excluded; stable argsort; lexicographic base_BGC breaks exact ties",
  bootstrap="Paired positive-query BGC resampling with replacement; percentile 95% CI for mean and median B-minus-A, fixed neighbor graphs",
  permutation="Same permutation of the two-column Saccharide membership tuple used for A and B; all 2502 labels permuted; positive queries recomputed",
  pca=dict(svd_solver="full",n_components=320,whiten=False,standardize=False),
  umap=dict(metric="cosine",n_neighbors=[15,30,50],min_dist=.1,n_components=2,random_state=SEED,transform_seed=SEED,n_jobs=1,n_epochs=500,init="spectral"),
  python=sys.version,platform=platform.platform(),packages=versions(),dependency_path=str(PREV/".deps")))
 js("results/consensus_validation.json",dict(passed=True,n_rows=12510,n_base=2502,suffixes=[-1,0,1,2,3],
  row_alignment_exact=True,labels_constant=True,unmasked_A_exact_extraction_rows=12510,
  any_saccharide=sum(r["any_saccharide"] for r in rows),pure_saccharide=sum(r["pure_saccharide"] for r in rows),
  protected_files=len(protected),runtime_seconds=time.time()-t))
 print(json.dumps(json.loads((ROOT/"results/consensus_validation.json").read_text())),flush=True)
if __name__=="__main__":main()
