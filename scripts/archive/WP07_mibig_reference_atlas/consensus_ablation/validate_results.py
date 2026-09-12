from common import *
from scipy.spatial.distance import cdist
from PIL import Image
def main():
 m,arrays=load();I=json.loads((ROOT/"provenance/inputs.json").read_text())
 checks=[]
 def check(name,value,detail=""):
  checks.append(dict(check=name,passed=bool(value),detail=detail))
  if not value:
   tab("provenance/independent_validation.tsv",checks);raise AssertionError(name+": "+str(detail))
 check("2502 unique metadata rows",len(m)==2502 and m.base_BGC.is_unique)
 rawmeta=read(I["core_metadata"])
 for rep,key in [("A","published"),("B","core")]:
  raw=np.load(I[key]);recomputed=[]
  for r in m.itertuples():
   ix=json.loads(r.source_row_indices);ids=json.loads(r.source_row_IDs)
   check(f"{rep} {r.base_BGC} five-source mapping",len(ix)==5 and rawmeta.ID.iloc[ix].tolist()==ids and {int(i.rsplit("_",1)[1]) for i in ids}=={-1,0,1,2,3})
   total=np.zeros(320,np.float64)
   for i in ix:total+=raw[i].astype(np.float64)
   recomputed.append(total/5)
  check(rep+" independently summed raw consensus",np.array_equal(recomputed,arrays[rep]))
  check(rep+" finite float64",arrays[rep].dtype==np.float64 and np.isfinite(arrays[rep]).all())
  d=cdist(arrays[rep],arrays[rep],metric="cosine");np.fill_diagonal(d,np.inf)
  order=np.argsort(d,axis=1,kind="stable")[:,:50]
  saved=np.load(ROOT/"results"/f"{rep}_neighbor_indices_top50.npy")
  check(rep+" all neighbors independent scipy cosine",np.array_equal(order,saved))
  check(rep+" all queries excluded",not np.any(saved==np.arange(2502)[:,None]))
  frame=read(ROOT/f"tables/{rep}_nearest_neighbors_top50.tsv")
  check(rep+" neighbor TSV identity",np.array_equal(frame.neighbor_BGC.to_numpy(),m.base_BGC.to_numpy()[saved].ravel()))
  check(rep+" neighbor TSV distances",np.allclose(frame.cosine_distance,np.take_along_axis(d,saved,axis=1).ravel(),atol=1e-12,rtol=1e-10))
  p=read(ROOT/f"results/pca_{rep}_coordinates.tsv")
  co=read(ROOT/f"results/pca_{rep}_components.tsv").drop(columns="PC").to_numpy()
  center=read(ROOT/f"results/pca_{rep}_center.tsv")["mean"].to_numpy()
  scores=p[[f"PC{j}" for j in range(1,321)]].to_numpy()
  check(rep+" PCA full reconstruction",np.allclose(scores@co+center,arrays[rep],atol=1e-10,rtol=1e-10))
 for p in sorted((ROOT/"results").glob("*coordinates.tsv")):
  f=read(p)
  check(p.name+" exact metadata join",np.array_equal(f.base_BGC,m.base_BGC) and np.array_equal(f.label,m.label))
  xyz=f.filter(regex=r"^(PC|UMAP)[0-9]+$").to_numpy()
  check(p.name+" all finite",np.isfinite(xyz).all())
 s=read(ROOT/"tables/A_vs_B_saccharide_knn_summary.tsv")
 detail=read(ROOT/"tables/saccharide_knn_per_query.tsv")
 nn={rep:np.load(ROOT/f"results/{rep}_neighbor_indices_top50.npy") for rep in ["A","B"]}
 rng=np.random.default_rng(SEED+1)
 for definition in DEFS:
  mask=m[definition].to_numpy(bool);count=int(mask.sum());draws=rng.integers(count,size=(N_BOOT,count))
  for k in KS:
   row=s[(s.definition==definition)&(s.k==k)].iloc[0]
   av=mask[nn["A"][mask,:k]].sum(1)/k;bv=mask[nn["B"][mask,:k]].sum(1)/k
   for rep,v in [("A",av),("B",bv)]:
    check(f"{definition} {k} {rep} precision",abs(v.mean()-row[rep+"_precision"])<1e-14)
   delta=bv-av;ci=np.percentile(delta[draws].mean(1),[2.5,97.5])
   check(f"{definition} {k} mean delta CI",np.allclose(ci,[row.paired_bootstrap_CI_low,row.paired_bootstrap_CI_high],atol=1e-14,rtol=0))
   check(f"{definition} {k} paired fraction partition",abs(row.fraction_queries_B_better+row.fraction_queries_equal+row.fraction_queries_B_worse-1)<1e-14)
 perms=np.load(ROOT/"results/shared_label_permutations.npy",mmap_mode="r")
 null=np.load(ROOT/"results/permutation_null_precision.npy",mmap_mode="r")
 rng=np.random.default_rng(SEED)
 check("shared permutation shape",perms.shape==(10000,2502))
 for i in range(N_PERM):
  check(f"shared permutation {i} seed and bijection",np.array_equal(perms[i],rng.permutation(2502)))
 for i in list(range(50))+list(range(N_PERM-50,N_PERM)):
  for j,definition in enumerate(DEFS):
   mask=m[definition].to_numpy(bool)[perms[i]]
   for r,rep in enumerate(["A","B"]):
    for h,k in enumerate(KS):
     actual=np.mean([np.count_nonzero(mask[neighbors[:k]])/k for neighbors in nn[rep][mask]])
     check(f"null independent {i} {rep} {definition} {k}",abs(actual-null[i,r,j,h])<1e-14)
 for row in s.itertuples():
  j=DEFS.index(row.definition);h=KS.index(row.k)
  for r,rep in enumerate(["A","B"]):
   p=(1+np.count_nonzero(null[:,r,j,h]>=getattr(row,rep+"_precision")))/(N_PERM+1)
   check(f"empirical p {rep} {j} {h}",abs(p-getattr(row,rep+"_empirical_p"))<1e-14)
 old=read(PROJECT/"analysis/bgcprophet_reference/results/cosine_neighbors_top20.tsv")
 prevmap=old.groupby("BGC",sort=True).neighbor.agg(list).to_dict()
 ov=np.array([len(set(prevmap[bg])&set(m.base_BGC.to_numpy()[neighbors[:20]]))/20 for bg,neighbors in zip(m.base_BGC,nn["A"])])
 tab("tables/A_previous_consensus_neighbor_reproduction.tsv",dict(base_BGC=m.base_BGC,overlap_at_20=ov))
 oldmat=np.load(PROJECT/"analysis/bgcprophet_reference/results/consensus_embeddings.npy")
 oldm=read(PROJECT/"analysis/bgcprophet_reference/results/consensus_metadata.tsv")
 oldmat=oldmat[pd.Index(oldm.BGC).get_indexer(m.base_BGC)]
 js("results/A_reference_geometry_check.json",dict(mean_top20_overlap=ov.mean(),min_top20_overlap=ov.min(),consensus_max_abs_difference=float(np.max(abs(oldmat-arrays["A"]))),old_dtype=str(oldmat.dtype),new_dtype=str(arrays["A"].dtype)))
 for p in (ROOT/"figures").glob("*.png"):
  with Image.open(p) as im:im.verify()
  check(p.name+" readable PNG",True)
  check(p.stem+" PDF exists",(p.with_suffix(".pdf")).stat().st_size>1000)
 count=integrity();check("protected inputs unchanged",True,str(count)+" files")
 tab("provenance/independent_validation.tsv",checks)
 js("results/final_validation.json",dict(passed=True,checks=len(checks),protected_files_unchanged=count,previous_A_neighbor_mean_overlap=float(ov.mean())))
 print(json.dumps(json.loads((ROOT/"results/final_validation.json").read_text())),flush=True)
if __name__=="__main__":main()
