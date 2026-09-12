from common import *
def main():
 t=time.time();I=json.loads((CONS/"provenance/inputs.json").read_text())
 paths=[]
 for directory in [PROJECT/"analysis/bgcprophet_reference",OLD,CONS,Path(I["canonical_directory"])]:
  for parent,dirs,files in os.walk(directory):
   dirs[:]=[d for d in dirs if d not in [".deps",".cache",".tmp","__pycache__",".git"]]
   paths.extend(Path(parent)/name for name in files)
 paths.append(Path(I["csv"]))
 paths=sorted(set(paths))
 tab("provenance/input_sha256_before.tsv",[dict(path=str(p),sha256=sha(p)) for p in paths])
 prior=read(CONS/"provenance/primary_inputs.tsv")
 for key in ["published","published_metadata","core","core_metadata","csv"]:
  q=prior[prior.path==I[key]];assert len(q)==1 and q.sha256.iloc[0]==sha(I[key])
 A=np.load(I["published"]);B=np.load(I["core"])
 am=read(I["published_metadata"]);bm=read(I["core_metadata"]);csv=pd.read_csv(I["csv"])
 assert A.shape==B.shape==(12510,320) and A.dtype==B.dtype==np.float32
 assert np.isfinite(A).all() and np.isfinite(B).all()
 assert am.equals(bm) and len(am)==12510 and am.ID.is_unique
 assert np.array_equal(am.ID,csv.ID) and np.array_equal(am.labels,csv.labels)
 assert np.array_equal(am.embedding_row,np.arange(12510))
 valid=read(Path(I["canonical_directory"])/"core_context_full_validation.tsv")
 assert np.array_equal(valid.ID,am.ID) and np.array_equal(valid.embedding_row,am.embedding_row)
 n=csv.sentence.str.split().str.len().to_numpy()
 assert np.array_equal(n,valid.core_length) and np.array_equal(n,valid.selected_positions)
 assert valid.selected_count_matches_core.all() and valid.unmasked_reference_exact.all()
 assert np.all((n>0)&(n<128))
 parts=am.ID.str.extract(r"^(BGC[0-9]{7})_(-1|0|1|2|3)$");assert not parts.isna().any().any()
 meta=pd.DataFrame(dict(embedding_row=am.embedding_row,ID=am.ID,base_BGC=parts[0],augmentation=parts[1].astype(int),
  label=am.labels,core_gene_count=n,noncore_position_count=128-n,interrupted_core_status=valid.core_interrupted))
 groups=[];bgs=[]
 for bg,g in meta.groupby("base_BGC",sort=True):
  g=g.sort_values("augmentation")
  assert len(g)==5 and set(g.augmentation)=={-1,0,1,2,3} and g.label.nunique()==1 and g.core_gene_count.nunique()==1
  groups.append(g.embedding_row.to_numpy());bgs.append(bg)
 assert len(bgs)==2502
 ix=np.array(groups);np.save(ROOT/"results/group_row_indices.npy",ix)
 A=A.astype(np.float64);B=B.astype(np.float64);N=(128*A-n[:,None]*B)/(128-n[:,None])
 assert N.shape==(12510,320) and np.isfinite(N).all()
 np.save(ROOT/"results/noncore_context_embeddings.npy",N)
 tab("results/noncore_context_metadata.tsv",meta)
 recon=n[:,None]/128*B+(128-n[:,None])/128*N;err=abs(recon-A)
 tol=1e-12+1e-12*abs(A);passed=np.all(err<=tol,axis=1);assert passed.all()
 relnorm=np.linalg.norm(recon-A,axis=1)/np.linalg.norm(A,axis=1)
 component=abs(A)>1e-8
 tab("tables/noncore_reconstruction_validation.tsv",meta.assign(max_absolute_error=err.max(1),mean_absolute_error=err.mean(1),relative_L2_error=relnorm,passed=passed))
 js("results/reconstruction_summary.json",dict(shape=list(N.shape),dtype=str(N.dtype),finite=True,
  min_core_count=int(n.min()),max_core_count=int(n.max()),min_noncore_count=int((128-n).min()),rows_passed=int(passed.sum()),
  max_absolute_error=float(err.max()),mean_absolute_error=float(err.mean()),max_row_relative_L2_error=float(relnorm.max()),
  max_element_relative_error_where_abs_A_gt_1e_8=float((err[component]/abs(A[component])).max()),
  tolerance="abs(error) <= 1e-12 + 1e-12*abs(A); no bitwise equality required",source_precision="Saved A/B are float32; algebra and N storage are float64",
  n_matches_sentence_and_extraction=True,interrupted_rows=int(valid.core_interrupted.sum()),runtime_seconds=time.time()-t))
 I["noncore"]=str(ROOT/"results/noncore_context_embeddings.npy")
 js("provenance/inputs.json",I)
 js("provenance/parameters.json",dict(python=sys.version,platform=platform.platform(),bootstrap_count=BOOT,bootstrap_seed=SEED,
  comparison_seeds="1729 + paired-test index (0..5)",equality_tolerance=1e-12,ks=KS,
  reconstruction="N=(128*A-n*B)/(128-n); raw float64 algebra",
  centroid="Mean of five raw vectors, float64; no individual normalization",
  knn="Exact corrected Step 3 implementation: transient float64 cosine normalization, 1-Z@Z.T, diagonal infinity, stable argsort, lexicographic BGC ties",
  permutations="Reuse exact shared_label_permutations.npy from corrected consensus Step 3; verify seed/count/group order",
  no_inference=True,no_PCA_or_UMAP=True))
 print((ROOT/"results/reconstruction_summary.json").read_text(),flush=True)
if __name__=="__main__":main()
