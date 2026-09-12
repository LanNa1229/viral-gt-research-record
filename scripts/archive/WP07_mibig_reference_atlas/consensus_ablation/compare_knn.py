from common import *
def main():
 t=time.time();m,arrays=load();n=len(m);nn={};values={};detail=[];ties=[]
 for rep,a in arrays.items():
  z=a/np.linalg.norm(a,axis=1,keepdims=True);d=1-z@z.T;np.fill_diagonal(d,np.inf)
  nn[rep]=np.argsort(d,axis=1,kind="stable")[:,:50]
  np.save(ROOT/"results"/f"{rep}_neighbor_indices_top50.npy",nn[rep])
  tab(f"tables/{rep}_nearest_neighbors_top50.tsv",dict(query_BGC=np.repeat(m.base_BGC.to_numpy(),50),
   rank=np.tile(np.arange(1,51),n),neighbor_BGC=m.base_BGC.to_numpy()[nn[rep]].ravel(),
   cosine_distance=np.take_along_axis(d,nn[rep],axis=1).ravel()))
  for k in KS:
   order=np.argsort(d,axis=1,kind="stable")
   boundary=d[np.arange(n),order[:,k-1]]==d[np.arange(n),order[:,k]]
   ties.append(dict(representation=rep,k=k,queries_with_boundary_ties=int(boundary.sum())))
  for definition in DEFS:
   mask=m[definition].to_numpy(bool);values[rep,definition]=mask[nn[rep]].cumsum(axis=1)[:,np.array(KS)-1]/np.array(KS)
   for j,k in enumerate(KS):
    detail.extend(dict(representation=rep,definition=definition,k=k,query_BGC=bg,positive_neighbor_fraction=float(v))
     for bg,v in zip(m.base_BGC[mask],values[rep,definition][mask,j]))
 tab("tables/knn_boundary_ties.tsv",ties);tab("tables/saccharide_knn_per_query.tsv",detail)
 ref=read(PROJECT/"analysis/bgcprophet_reference/results/saccharide_knn_summary.tsv")
 check=[]
 for r in ref.itertuples():
  actual=float(values["A",r.definition][m[r.definition].to_numpy(bool),KS.index(r.k)].mean())
  check.append(dict(definition=r.definition,k=r.k,reference=r.precision_at_k,recomputed=actual,
    absolute_difference=abs(actual-r.precision_at_k),passed=abs(actual-r.precision_at_k)<1e-12))
 tab("tables/A_reference_reproduction.tsv",check)
 if not all(r["passed"] for r in check):raise RuntimeError("A consensus reference discrepancy: STOP before A/B interpretation")
 print("A reference gate PASS: all six stored consensus Precision@k values reproduced.",flush=True)
 memberships=m[DEFS].to_numpy(bool);rng=np.random.default_rng(SEED)
 null=np.empty((N_PERM,2,2,4),dtype=np.float64)
 permutations=np.empty((N_PERM,n),dtype=np.uint16)
 for b in range(N_PERM):
  p=rng.permutation(n);permutations[b]=p
  labels=memberships[p]
  for j in range(2):
   mask=labels[:,j]
   for r,rep in enumerate(["A","B"]):
    hits=mask[nn[rep][mask]].cumsum(1)[:,np.array(KS)-1]/np.array(KS)
    null[b,r,j]=hits.mean(0)
  if (b+1)%2000==0:print(f"Shared-label permutations {b+1}/{N_PERM}",flush=True)
 np.save(ROOT/"results/shared_label_permutations.npy",permutations)
 np.save(ROOT/"results/permutation_null_precision.npy",null)
 js("provenance/permutation_array_axes.json",dict(axes=["permutation","representation","definition","k"],representations=["A","B"],definitions=DEFS,k=KS))
 nullrows=[];summary=[];overlap=[]
 rng=np.random.default_rng(SEED+1)
 for j,definition in enumerate(DEFS):
  mask=memberships[:,j];count=int(mask.sum());prev=count/n
  draws=rng.integers(count,size=(N_BOOT,count))
  for h,k in enumerate(KS):
   av=values["A",definition][mask,h];bv=values["B",definition][mask,h];delta=bv-av
   boots=delta[draws];ci=np.quantile(boots.mean(1),[.025,.975]);mic=np.quantile(np.median(boots,axis=1),[.025,.975])
   row=dict(definition=definition,k=k,positive_BGCs=count,prevalence=prev,self_excluded_null_expectation=(count-1)/(n-1),
    A_precision=av.mean(),B_precision=bv.mean(),A_fold_enrichment=av.mean()/prev,B_fold_enrichment=bv.mean()/prev,
    A_enrichment_percentage_points=100*(av.mean()-prev),B_enrichment_percentage_points=100*(bv.mean()-prev),
    B_minus_A_precision=delta.mean(),B_minus_A_percentage_points=100*delta.mean(),median_paired_difference=np.median(delta),
    paired_bootstrap_CI_low=ci[0],paired_bootstrap_CI_high=ci[1],median_bootstrap_CI_low=mic[0],median_bootstrap_CI_high=mic[1],
    fraction_queries_B_better=(delta>0).mean(),fraction_queries_equal=(delta==0).mean(),fraction_queries_B_worse=(delta<0).mean())
   for r,rep in enumerate(["A","B"]):
    observed=row[rep+"_precision"];nv=null[:,r,j,h]
    row[rep+"_null_mean"]=nv.mean();row[rep+"_null_sd"]=nv.std(ddof=1)
    row[rep+"_empirical_p"]=(1+np.count_nonzero(nv>=observed))/(N_PERM+1)
    row[rep+"_z_score"]=(observed-nv.mean())/nv.std(ddof=1) if nv.std(ddof=1)>0 else np.nan
    nullrows.extend(dict(permutation=b,representation=rep,definition=definition,k=k,null_precision=float(x)) for b,x in enumerate(nv))
   summary.append(row)
 tab("tables/A_vs_B_saccharide_knn_summary.tsv",summary)
 tab("tables/permutation_null_precision.tsv",nullrows)
 del nullrows
 per=[]
 for k in KS:
  ov=np.array([len(set(a[:k])&set(b[:k]))/k for a,b in zip(nn["A"],nn["B"])])
  per.extend(dict(base_BGC=bg,k=k,overlap=float(v)) for bg,v in zip(m.base_BGC,ov))
  for name,mask in [("all",np.ones(n,bool)),("any_saccharide",memberships[:,0]),("non_any_saccharide",~memberships[:,0]),("pure_saccharide",memberships[:,1]),("non_pure_saccharide",~memberships[:,1])]:
   v=ov[mask];overlap.append(dict(group=name,k=k,n_BGC=len(v),mean_overlap=v.mean(),median_overlap=np.median(v),q25=np.quantile(v,.25),q75=np.quantile(v,.75)))
 tab("tables/neighbor_overlap_per_BGC.tsv",per);tab("tables/neighbor_overlap_summary.tsv",overlap)
 earlier=read(PREV/"tables/knn_class_enrichment.tsv");comparisons=[]
 for r in pd.DataFrame(summary).itertuples():
  for rep,oldrep in [("A","published_style"),("B","core_context")]:
   group="Any Saccharide" if r.definition=="any_saccharide" else "Saccharide"
   old=earlier[(earlier.representation==oldrep)&(earlier.group==group)&(earlier.k==r.k)]
   assert len(old)==1
   comparisons.append(dict(representation=rep,definition=r.definition,k=r.k,
    superseded_zero_precision=old.iloc[0].observed_same_class_neighbor_fraction,
    consensus_precision=getattr(r,rep+"_precision")))
 tab("tables/consensus_vs_superseded_zero.tsv",comparisons)
 js("results/knn_completion.json",dict(passed=True,runtime_seconds=time.time()-t,permutations=N_PERM))
 print(pd.DataFrame(summary)[["definition","k","A_precision","B_precision","B_minus_A_precision","paired_bootstrap_CI_low","paired_bootstrap_CI_high"]].to_string(index=False),flush=True)
if __name__=="__main__":main()
