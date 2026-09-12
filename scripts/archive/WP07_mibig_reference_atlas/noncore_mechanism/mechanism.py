from common import *
def main():
 t=time.time();arrays=raw();ix=np.load(ROOT/"results/group_row_indices.npy")
 meta=read(ROOT/"results/noncore_context_metadata.tsv");base=meta.iloc[ix[:,0]].reset_index(drop=True)
 wide=base[["base_BGC","label","core_gene_count","noncore_position_count"]].copy()
 pairs=[];centers=[];globalrows=[];garrays={r:a[ix] for r,a in arrays.items()}
 metrics=["mean_pairwise_cosine","median_pairwise_cosine","min_pairwise_cosine","std_pairwise_cosine","mean_centroid_cosine_distance","max_centroid_cosine_distance"]
 for rep,g in garrays.items():
  assert np.all(np.linalg.norm(g,axis=2)>0)
  pp=np.stack([cosine(g[:,i],g[:,j]) for i,j in itertools.combinations(range(5),2)],axis=1)
  cent=g.mean(1,dtype=np.float64);assert np.all(np.linalg.norm(cent,axis=1)>0)
  cd=1-cosine(g,cent[:,None,:])
  values=[pp.mean(1),np.median(pp,axis=1),pp.min(1),pp.std(1,ddof=1),cd.mean(1),cd.max(1)]
  for key,v in zip(metrics,values):
   wide[rep+"_"+key]=v;globalrows.append(dict(representation=rep,metric=key,unit="base_BGC",**describe(v)))
  globalrows.append(dict(representation=rep,metric="all_pairwise_cosines",unit="augmentation_pair_descriptive",**describe(pp.ravel())))
  for col,(i,j) in enumerate(itertools.combinations(range(5),2)):
   pairs.extend(dict(representation=rep,base_BGC=bg,ID1=a,ID2=b,cosine=float(v)) for bg,a,b,v in zip(base.base_BGC,meta.ID.to_numpy()[ix[:,i]],meta.ID.to_numpy()[ix[:,j]],pp[:,col]))
  centers.extend(dict(representation=rep,base_BGC=bg,ID=id,centroid_cosine_distance=float(v))
    for bg,id,v in zip(np.repeat(base.base_BGC.to_numpy(),5),meta.ID.to_numpy()[ix.ravel()],cd.ravel()))
 tab("tables/A_B_N_augmentation_stability.tsv",wide)
 tab("tables/A_B_N_stability_summary.tsv",globalrows)
 tab("tables/all_augmentation_pairwise_cosines.tsv",pairs)
 tab("tables/augmentation_centroid_distances.tsv",centers)
 tests=[]
 for left,right in [("B","A"),("N","A"),("N","B")]:
  for metric in ["mean_pairwise_cosine","min_pairwise_cosine"]:
   d=(wide[left+"_"+metric]-wide[right+"_"+metric]).to_numpy()
   tests.append(dict(comparison=left+" - "+right,metric=metric,n_BGC=len(d),**bootdiff(d,SEED+len(tests))))
 tab("tables/paired_stability_differences.tsv",tests)
 stronger=(wide.A_mean_pairwise_cosine>wide.B_mean_pairwise_cosine+1e-12)&(wide.A_mean_pairwise_cosine>wide.N_mean_pairwise_cosine+1e-12)
 js("results/stability_comparison.json",dict(fraction_A_more_stable_than_B=float(np.mean(wide.A_mean_pairwise_cosine>wide.B_mean_pairwise_cosine+1e-12)),
  fraction_A_more_stable_than_N=float(np.mean(wide.A_mean_pairwise_cosine>wide.N_mean_pairwise_cosine+1e-12)),
  fraction_A_more_stable_than_both=float(stronger.mean()),count_A_more_stable_than_both=int(stronger.sum()),equality_tolerance=1e-12))
 deviations={rep:g-g.mean(1)[:,None,:] for rep,g in garrays.items()}
 wB=base.core_gene_count.to_numpy()/128;wN=1-wB
 core=deviations["B"]*wB[:,None,None];noncore=deviations["N"]*wN[:,None,None]
 da=deviations["A"];vectorerr=np.max(abs(da-core-noncore),axis=2)
 eb=np.sum(core**2,axis=2);en=np.sum(noncore**2,axis=2);ea=np.sum(da**2,axis=2)
 cross=2*np.sum(core*noncore,axis=2);dotcos=cosine(core,noncore)
 energyerr=abs(ea-eb-en-cross)
 assert np.all(vectorerr<=1e-12) and np.all(energyerr<=1e-12+1e-12*abs(ea))
 frame=meta.iloc[ix.ravel()].reset_index(drop=True).copy()
 for key,v in dict(wB=np.repeat(wB,5),wN=np.repeat(wN,5),core_term_squared_norm=eb.ravel(),noncore_term_squared_norm=en.ravel(),
  A_deviation_squared_norm=ea.ravel(),cross_term=cross.ravel(),deviation_cosine=dotcos.ravel(),
  deviation_cosine_defined=np.isfinite(dotcos).ravel(),deviation_identity_max_abs_error=vectorerr.ravel(),variance_identity_abs_error=energyerr.ravel()).items():frame[key]=v
 tab("tables/core_noncore_cancellation.tsv",frame)
 aggregate=wide.copy()
 aggregate["wB"]=wB;aggregate["wN"]=wN
 aggregate["mean_cross_term"]=cross.mean(1);aggregate["median_cross_term"]=np.median(cross,axis=1)
 aggregate["fraction_augmentations_negative_cross"]=(cross<0).mean(1)
 aggregate["mean_deviation_cosine"]=np.nanmean(dotcos,axis=1)
 aggregate["mean_core_energy"]=eb.mean(1);aggregate["mean_noncore_energy"]=en.mean(1)
 aggregate["mean_A_deviation_energy"]=ea.mean(1);denom=(eb+en).mean(1)
 aggregate["variance_reduction_fraction_from_cross"]=np.divide(-cross.mean(1),denom,out=np.full(len(denom),np.nan),where=denom>0)
 aggregate["noncore_fraction_of_uncoupled_energy"]=np.divide(en.mean(1),denom,out=np.full(len(denom),np.nan),where=denom>0)
 aggregate["A_minus_B_stability"]=wide.A_mean_pairwise_cosine-wide.B_mean_pairwise_cosine
 aggregate["core_gene_count_bin"]=pd.cut(base.core_gene_count,[0,5,10,20,40,80,128],labels=["1-5","6-10","11-20","21-40","41-80","81-127"]).astype(str)
 tab("tables/per_BGC_cancellation.tsv",aggregate)
 summary=[]
 for scope,values in [("augmentation_cross_term",cross.ravel()),("augmentation_deviation_cosine",dotcos.ravel()),
  ("per_BGC_mean_cross_term",cross.mean(1)),("per_BGC_variance_reduction_fraction",aggregate.variance_reduction_fraction_from_cross)]:
  summary.append(dict(metric=scope,**describe(values)))
 tab("tables/cancellation_summary.tsv",summary)
 canceldict=dict(fraction_augmentation_negative_cross=float((cross<0).mean()),fraction_BGC_negative_mean_cross=float((cross.mean(1)<0).mean()),
  median_cross_term=float(np.median(cross)),median_deviation_cosine=float(np.nanmedian(dotcos)),undefined_cosines=int((~np.isfinite(dotcos)).sum()),
  median_BGC_variance_reduction_fraction=float(np.nanmedian(aggregate.variance_reduction_fraction_from_cross)),
  pooled_variance_reduction_fraction=float(-cross.sum()/(eb+en).sum()),
  max_deviation_identity_error=float(vectorerr.max()),max_variance_identity_error=float(energyerr.max()))
 js("results/cancellation_overview.json",canceldict)
 correlations=[]
 for x,y in [("core_gene_count",r+"_mean_pairwise_cosine") for r in ["A","B","N"]]+[
  ("core_gene_count","A_minus_B_stability"),("noncore_position_count","A_minus_B_stability"),
  ("core_gene_count","mean_cross_term"),("core_gene_count","variance_reduction_fraction_from_cross"),
  ("core_gene_count","noncore_fraction_of_uncoupled_energy"),
  ("mean_cross_term","A_minus_B_stability"),("variance_reduction_fraction_from_cross","A_minus_B_stability")]:
  a=aggregate[x].to_numpy();b=aggregate[y].to_numpy();mask=np.isfinite(a)&np.isfinite(b)
  rho,p=stats.spearmanr(a[mask],b[mask]);correlations.append(dict(x=x,y=y,n_BGC=int(mask.sum()),spearman_rho=rho,p_value=p,absolute_rho=abs(rho)))
 tab("tables/gene_count_mechanism_correlations.tsv",correlations)
 strata=[]
 for name,g in aggregate.groupby("core_gene_count_bin",sort=False):
  row=dict(core_gene_count_bin=name,n_BGC=len(g),median_core_count=g.core_gene_count.median())
  for col in ["A_mean_pairwise_cosine","B_mean_pairwise_cosine","N_mean_pairwise_cosine","mean_cross_term","variance_reduction_fraction_from_cross","noncore_fraction_of_uncoupled_energy","A_minus_B_stability"]:
   row["median_"+col]=g[col].median()
  row["fraction_BGC_negative_mean_cross"]=(g.mean_cross_term<0).mean()
  strata.append(row)
 tab("tables/core_count_stratified_mechanism.tsv",strata)
 old=read(OLD/"tables/base_bgc_augmentation_stability.tsv")
 js("results/mechanism_completion.json",dict(passed=True,runtime_seconds=time.time()-t))
 print(pd.DataFrame(globalrows).query("metric == 'mean_pairwise_cosine'")[["representation","mean","median"]].to_string(index=False),flush=True)
 print(json.dumps(canceldict,indent=2),flush=True)
 print(pd.DataFrame(correlations).to_string(index=False),flush=True)
if __name__=="__main__":main()
