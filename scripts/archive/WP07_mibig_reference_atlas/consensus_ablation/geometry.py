from common import *
from sklearn.decomposition import PCA
import umap
I=json.loads((ROOT/"provenance/inputs.json").read_text())
COLORS=json.loads(Path(I["palette"]).read_text())["colors"]
CLASSES=list(COLORS)
TITLES={"A":"A: published-style consensus","B":"B: core-context consensus"}
def classes(ax,xy,m,title,xlabel,ylabel):
 for c in CLASSES:
  mask=m.plot_class.to_numpy()==c
  ax.scatter(xy[mask,0],xy[mask,1],s=10,alpha=.7,c=COLORS[c],linewidths=0,rasterized=True,label=f"{c} ({mask.sum()})")
 ax.set(title=title,xlabel=xlabel,ylabel=ylabel)
 ax.legend(bbox_to_anchor=(1.02,1),loc="upper left",fontsize=8)
def highlight(ax,xy,mask,title,xlabel,ylabel):
 ax.scatter(xy[~mask,0],xy[~mask,1],s=8,c="0.78",alpha=.45,linewidths=0,rasterized=True,label=f"Negative ({(~mask).sum()})")
 ax.scatter(xy[mask,0],xy[mask,1],s=15,c=COLORS["Saccharide"],alpha=.85,edgecolors="0.2",linewidths=.2,rasterized=True,label=f"Positive ({mask.sum()})")
 ax.set(title=title,xlabel=xlabel,ylabel=ylabel);ax.legend(fontsize=8)
def coords(name,xy,m,prefix):
 tab("results/"+name,pd.concat([m.reset_index(drop=True),pd.DataFrame(xy,columns=[f"{prefix}{j+1}" for j in range(xy.shape[1])])],axis=1))
def main():
 t=time.time();(ROOT/"provenance/figure_sources.tsv").write_text("")
 assert json.loads((ROOT/"results/knn_completion.json").read_text())["passed"]
 m,arrays=load();pcs={};models={};vr=[];thresholds=[];diag=[]
 for rep,a in arrays.items():
  model=PCA(n_components=320,svd_solver="full",whiten=False)
  pcs[rep]=model.fit_transform(a);models[rep]=model
  coords(f"pca_{rep}_coordinates.tsv",pcs[rep],m,"PC")
  tab(f"results/pca_{rep}_components.tsv",pd.DataFrame(model.components_,columns=[f"dim_{j+1}" for j in range(320)]).assign(PC=range(1,321)))
  tab(f"results/pca_{rep}_center.tsv",dict(dimension=range(1,321),mean=model.mean_))
  for j in range(320):
   vr.append(dict(representation=rep,PC=j+1,explained_variance=model.explained_variance_[j],explained_ratio=model.explained_variance_ratio_[j],cumulative_ratio=model.explained_variance_ratio_[:j+1].sum()))
   diag.append(dict(representation=rep,dimension=j+1,mean=a[:,j].mean(),sd=a[:,j].std(ddof=1),minimum=a[:,j].min(),maximum=a[:,j].max()))
  for q in [.5,.75,.8,.9,.95]:
   thresholds.append(dict(representation=rep,variance_fraction=q,n_PCs=int(np.searchsorted(model.explained_variance_ratio_.cumsum(),q)+1)))
 tab("tables/pca_explained_variance.tsv",vr);tab("tables/pca_variance_thresholds.tsv",thresholds);tab("tables/dimension_scale_diagnostics.tsv",diag)
 both=np.concatenate([pcs[r][:,:2] for r in ["A","B"]]);lo=both.min(0);hi=both.max(0);pad=(hi-lo)*.05
 for rep in arrays:
  model=models[rep];labels=[f"PC{j+1} ({model.explained_variance_ratio_[j]:.1%})" for j in range(2)]
  source=f"results/pca_{rep}_coordinates.tsv"
  fig,ax=plt.subplots(figsize=(7.6,5))
  classes(ax,pcs[rep][:,:2],m,TITLES[rep],*labels);ax.set(xlim=(lo[0]-pad[0],hi[0]+pad[0]),ylim=(lo[1]-pad[1],hi[1]+pad[1]))
  save(fig,f"pca_{rep}_pc1_pc2_classes",source)
  fig,axs=plt.subplots(1,2,figsize=(11,4.5))
  for ax,definition in zip(axs,DEFS):
   highlight(ax,pcs[rep][:,:2],m[definition].to_numpy(bool),definition.replace("_"," ").title(),*labels)
   ax.set(xlim=(lo[0]-pad[0],hi[0]+pad[0]),ylim=(lo[1]-pad[1],hi[1]+pad[1]))
  fig.suptitle(TITLES[rep]);fig.tight_layout();save(fig,f"pca_{rep}_pc1_pc2_saccharide",source)
  fig,axs=plt.subplots(1,2,figsize=(9,3.8))
  axs[0].plot(range(1,21),model.explained_variance_ratio_[:20],"o-",color="#0072B2")
  axs[0].set(xlabel="PC",ylabel="Explained variance ratio",ylim=(0,.5),xticks=[1,5,10,15,20])
  axs[1].plot(range(1,321),model.explained_variance_ratio_.cumsum(),color="#0072B2")
  axs[1].set(xlabel="Number of PCs",ylabel="Cumulative variance",ylim=(0,1.01))
  fig.suptitle(TITLES[rep]);fig.tight_layout();save(fig,f"pca_{rep}_scree","tables/pca_explained_variance.tsv")
 params=dict(metric="cosine",min_dist=.1,n_components=2,random_state=SEED,transform_seed=SEED,n_jobs=1,n_epochs=500,init="spectral")
 runs=[];xy={}
 for rep,a in arrays.items():
  for k in [15,30,50]:
   start=time.time()
   xy[rep,k]=umap.UMAP(n_neighbors=k,**params).fit_transform(a)
   assert xy[rep,k].shape==(2502,2) and np.isfinite(xy[rep,k]).all()
   coords(f"umap_{rep}_n{k}_coordinates.tsv",xy[rep,k],m,"UMAP")
   runs.append(dict(representation=rep,n_neighbors=k,seconds=time.time()-start,**params))
   print(f"UMAP {rep} n_neighbors={k} finished in {time.time()-start:.1f}s",flush=True)
 tab("provenance/umap_runs.tsv",runs)
 for rep in arrays:
  source=f"results/umap_{rep}_n30_coordinates.tsv"
  fig,ax=plt.subplots(figsize=(7.6,5));classes(ax,xy[rep,30],m,TITLES[rep],"UMAP 1","UMAP 2")
  save(fig,f"umap_{rep}_classes",source)
  fig,axs=plt.subplots(1,2,figsize=(11,4.5))
  for ax,definition in zip(axs,DEFS):highlight(ax,xy[rep,30],m[definition].to_numpy(bool),definition.replace("_"," ").title(),"UMAP 1","UMAP 2")
  fig.suptitle(TITLES[rep]);fig.tight_layout();save(fig,f"umap_{rep}_saccharide",source)
 fig,axs=plt.subplots(2,3,figsize=(14,8))
 for i,rep in enumerate(["A","B"]):
  for j,k in enumerate([15,30,50]):
   highlight(axs[i,j],xy[rep,k],m.any_saccharide.to_numpy(bool),f"{rep}: neighbors={k}","UMAP 1","UMAP 2")
 fig.suptitle("Consensus UMAP robustness: Any Saccharide; all other settings fixed")
 fig.tight_layout();save(fig,"umap_robustness",";".join(f"results/umap_{r}_n{k}_coordinates.tsv" for r in arrays for k in [15,30,50]))
 s=read(ROOT/"tables/A_vs_B_saccharide_knn_summary.tsv")
 fig,axs=plt.subplots(1,2,figsize=(11,4.5))
 for ax,definition in zip(axs,DEFS):
  g=s[s.definition==definition]
  for rep,col,marker in [("A","#0072B2","o"),("B","#D55E00","s")]:
   ax.plot(g.k,g[rep+"_precision"]*100,marker=marker,color=col,label=TITLES[rep])
  ax.axhline(g.prevalence.iloc[0]*100,color=".4",ls="--",label="Global prevalence")
  ax.set(title=definition.replace("_"," ").title(),xlabel="k",ylabel="Precision@k (%)",ylim=(0,100),xticks=KS)
  ax.legend(fontsize=8)
 fig.tight_layout();save(fig,"saccharide_knn_comparison","tables/A_vs_B_saccharide_knn_summary.tsv")
 fig,axs=plt.subplots(1,2,figsize=(10,4))
 for ax,definition in zip(axs,DEFS):
  g=s[s.definition==definition];y=g.B_minus_A_precision*100
  ax.errorbar(g.k,y,yerr=np.array([y-g.paired_bootstrap_CI_low*100,g.paired_bootstrap_CI_high*100-y]),fmt="o",capsize=4,color="#0072B2")
  ax.axhline(0,c=".4",ls="--");ax.set(title=definition.replace("_"," ").title(),xlabel="k",ylabel="B minus A (percentage points)",xticks=KS,ylim=(-25,2))
 fig.tight_layout();save(fig,"saccharide_paired_difference","tables/A_vs_B_saccharide_knn_summary.tsv")
 ov=read(ROOT/"tables/neighbor_overlap_summary.tsv")
 fig,ax=plt.subplots(figsize=(7,4))
 for group,c,marker in [("all","#777777","o"),("any_saccharide","#0072B2","s"),("non_any_saccharide","#D55E00","^")]:
  g=ov[ov.group==group];ax.plot(g.k,g.mean_overlap,marker=marker,c=c,label=group.replace("_"," "))
 ax.set(xlabel="k",ylabel="Mean neighbor overlap",ylim=(0,1),xticks=KS);ax.legend()
 save(fig,"neighbor_overlap","tables/neighbor_overlap_summary.tsv")
 js("results/geometry_completion.json",dict(passed=True,runtime_seconds=time.time()-t,figures=len(list((ROOT/"figures").glob("*.png")))))
if __name__=="__main__":main()
