from common import *
def main():
 t=time.time();arrays=raw();ix=np.load(ROOT/"results/group_row_indices.npy")
 meta=read(ROOT/"results/noncore_context_metadata.tsv");cm=read(CONS/"results/consensus_metadata.tsv")
 assert np.array_equal(meta.base_BGC.to_numpy()[ix[:,0]],cm.base_BGC)
 assert np.array_equal(meta.label.to_numpy()[ix[:,0]],cm.label)
 label=cm.label.str.strip()
 assert np.array_equal(cm.any_saccharide,label.str.split().apply(lambda t:"Saccharide" in t))
 assert np.array_equal(cm.pure_saccharide,label=="Saccharide")
 cons={r:a[ix].mean(1,dtype=np.float64) for r,a in arrays.items()}
 reproduction=[]
 for r,name in [("A","published"),("B","core_context")]:
  old=np.load(CONS/"results"/f"{name}_consensus_embeddings.npy")
  assert np.array_equal(cons[r],old),"Consensus mismatch: STOP"
  reproduction.append(dict(representation=r,consensus_exact=True,max_absolute_difference=float(abs(cons[r]-old).max())))
 for r,a in cons.items():
  assert a.shape==(2502,320) and np.isfinite(a).all() and np.all(np.linalg.norm(a,axis=1)>0)
  np.save(ROOT/f"results/{r}_consensus_embeddings.npy",a)
 cm["core_gene_count"]=meta.core_gene_count.to_numpy()[ix[:,0]]
 cm["noncore_position_count"]=128-cm.core_gene_count
 tab("results/consensus_metadata.tsv",cm)
 nn={};prec={};perquery=[];ties=[]
 for r,a in cons.items():
  z=a/np.linalg.norm(a,axis=1,keepdims=True);d=1-z@z.T;np.fill_diagonal(d,np.inf)
  order=np.argsort(d,axis=1,kind="stable");nn[r]=order[:,:50]
  if r in ["A","B"]:
   assert np.array_equal(nn[r],np.load(CONS/f"results/{r}_neighbor_indices_top50.npy")),"A/B neighbor mismatch: STOP"
  for k in KS:ties.append(dict(representation=r,k=k,n_boundary_ties=int(np.sum(d[np.arange(len(cm)),order[:,k-1]]==d[np.arange(len(cm)),order[:,k]]))))
  np.save(ROOT/f"results/{r}_neighbor_indices_top50.npy",nn[r])
  tab(f"tables/{r}_nearest_neighbors_top50.tsv",dict(query_BGC=np.repeat(cm.base_BGC.to_numpy(),50),
    rank=np.tile(np.arange(1,51),len(cm)),neighbor_BGC=cm.base_BGC.to_numpy()[nn[r]].ravel(),
    cosine_distance=np.take_along_axis(d,nn[r],axis=1).ravel()))
  for definition in DEFS:
   mask=cm[definition].to_numpy(bool);hits=mask[nn[r][mask]].cumsum(1)[:,np.array(KS)-1]/np.array(KS)
   prec[r,definition]=hits.mean(0)
   for h,k in enumerate(KS):perquery.extend(dict(representation=r,definition=definition,k=k,query_BGC=bg,positive_neighbor_fraction=v) for bg,v in zip(cm.base_BGC[mask],hits[:,h]))
 prior=read(CONS/"tables/A_vs_B_saccharide_knn_summary.tsv")
 for row in prior.itertuples():
  for r in ["A","B"]:
   assert abs(prec[r,row.definition][KS.index(row.k)]-getattr(row,r+"_precision"))<1e-12,"Prior precision mismatch: STOP"
 tab("tables/consensus_reproduction.tsv",reproduction);tab("tables/neighbor_boundary_ties.tsv",ties)
 tab("tables/A_B_N_saccharide_per_query.tsv",perquery)
 print("A/B consensus, top-50 neighbors, and all corrected precisions reproduced.",flush=True)
 path=CONS/"results/shared_label_permutations.npy";perm=np.load(path,mmap_mode="r")
 assert perm.shape==(10000,2502)
 rng=np.random.default_rng(SEED)
 for i in range(10000):assert np.array_equal(perm[i],rng.permutation(2502)),"Shared permutations incompatible"
 old_axes=json.loads((CONS/"provenance/permutation_array_axes.json").read_text())
 assert old_axes["definitions"]==DEFS and old_axes["k"]==KS
 membership=cm[DEFS].to_numpy(bool);null=np.empty((10000,2,4))
 for i in range(10000):
  lab=membership[perm[i]]
  for j in range(2):
   mask=lab[:,j];null[i,j]=(mask[nn["N"][mask]].cumsum(1)[:,np.array(KS)-1]/np.array(KS)).mean(0)
 np.save(ROOT/"results/N_permutation_null_precision.npy",null)
 js("provenance/permutation_reuse.json",dict(path=str(path),sha256=sha(path),shape=list(perm.shape),exact_seed_regeneration_verified=True,
  seed=SEED,assignments_reused_not_replaced=True,metadata_order_exact=True,definitions=DEFS,ks=KS))
 tab("tables/N_permutation_null_precision.tsv",[dict(permutation=i,definition=definition,k=k,null_precision=null[i,j,h])
   for i in range(10000) for j,definition in enumerate(DEFS) for h,k in enumerate(KS)])
 rows=[]
 for r in ["A","B","N"]:
  for j,definition in enumerate(DEFS):
   count=int(cm[definition].sum());prev=count/len(cm)
   for h,k in enumerate(KS):
    obs=prec[r,definition][h]
    row=dict(representation=r,definition=definition,k=k,positive_BGCs=count,total_BGCs=len(cm),precision_at_k=obs,
      prevalence=prev,self_excluded_null_expectation=(count-1)/(len(cm)-1),fold_enrichment=obs/prev,enrichment_percentage_points=100*(obs-prev))
    if r=="N":
     nv=null[:,j,h];row.update(null_mean=nv.mean(),empirical_p=(1+np.count_nonzero(nv>=obs))/10001,
       z_score=(obs-nv.mean())/nv.std(ddof=1),permutations=10000)
    else:
     old=prior[(prior.definition==definition)&(prior.k==k)].iloc[0]
     row.update(null_mean=old[r+"_null_mean"],empirical_p=old[r+"_empirical_p"],z_score=old[r+"_z_score"],permutations=10000)
    rows.append(row)
 tab("tables/A_B_N_saccharide_knn_summary.tsv",rows)
 overlaps=[];summary=[]
 for a,b in [("A","B"),("A","N"),("B","N")]:
  for k in KS:
   ov=np.array([len(set(x[:k])&set(y[:k]))/k for x,y in zip(nn[a],nn[b])])
   overlaps.extend(dict(comparison=a+" vs "+b,k=k,base_BGC=bg,overlap=v) for bg,v in zip(cm.base_BGC,ov))
   for scope,mask in [("all",np.ones(len(cm),bool)),("any_saccharide",cm.any_saccharide.to_numpy(bool)),("pure_saccharide",cm.pure_saccharide.to_numpy(bool)),("non_saccharide",~cm.any_saccharide.to_numpy(bool))]:
    summary.append(dict(comparison=a+" vs "+b,k=k,group=scope,**describe(ov[mask])))
 tab("tables/A_B_N_neighbor_overlap.tsv",overlaps);tab("tables/neighbor_overlap_summary.tsv",summary)
 prevov=read(CONS/"tables/neighbor_overlap_per_BGC.tsv")
 newov=pd.DataFrame(overlaps).query("comparison == 'A vs B'")
 joined=prevov.merge(newov,on=["base_BGC","k"],validate="one_to_one")
 assert np.allclose(joined.overlap_x,joined.overlap_y,atol=1e-15,rtol=0)
 js("results/geometry_validation.json",dict(passed=True,A_consensus_exact=True,B_consensus_exact=True,AB_neighbors_exact=True,
  AB_precision_reproduced=True,AB_overlap_reproduced=True,any_saccharide=int(cm.any_saccharide.sum()),pure_saccharide=int(cm.pure_saccharide.sum()),runtime_seconds=time.time()-t))
 print(pd.DataFrame(rows).pivot(index=["definition","k"],columns="representation",values="precision_at_k").to_string(),flush=True)
 print(pd.DataFrame(summary).query("k==20 and group=='all'")[["comparison","mean","median"]].to_string(index=False),flush=True)
if __name__=="__main__":main()
