"""Reproducible analysis of frozen HPC BGC-Prophet outputs; no inference."""
import os, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".deps"))
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[key] = "1"
for key, folder in (("MPLCONFIGDIR", ".cache/matplotlib"), ("NUMBA_CACHE_DIR", ".cache/numba")):
    path = ROOT / folder
    path.mkdir(parents=True, exist_ok=True)
    os.environ[key] = str(path)
import argparse, csv, hashlib, itertools, json, platform, re, datetime
import importlib.metadata as metadata
import numpy as np

VERSION = "1.0.1"
SEED = 1729
def sha_array(array):
    return hashlib.sha256(array.tobytes()).hexdigest()
def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()
def read(path, delimiter="\t"):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter=delimiter))
def table(name, rows, fields=None):
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else [])
    with open(OUT / name, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
def matrix(name, x, ids, prefix):
    table(name, (dict(BGC=b, **{f"{prefix}{j+1:03}": v for j,v in enumerate(row)}) for b,row in zip(ids,x)))
def quantiles(values, metric):
    return [{"metric":metric, "statistic":name, "value":float(value)} for name,value in
            [("n",len(values)),("mean",np.mean(values)),("std",np.std(values,ddof=1)),
             *[(f"q{q:g}",np.percentile(values,q)) for q in (0,1,5,25,50,75,95,99,100)]]]
def validate():
    paths = [DATA / f for f in ("mibig_12510_embeddings.npy","mibig_12510_metadata.tsv",
                               "mibig_12510_probabilities.npy","BGC_train_dataset_classify.csv")]
    paths.append(ROOT.parents[1] / "BGC_train_dataset_classify.csv")
    hashes = {str(p):sha(p) for p in paths}
    x = np.load(paths[0], allow_pickle=False)
    m = read(paths[1])
    p = np.load(paths[2], allow_pickle=False)
    c = read(paths[3], ",")
    checks=[]
    def check(name, condition, detail):
        checks.append(dict(check=name,passed=bool(condition),detail=detail))
        table("validation.tsv",checks)
        if not condition: raise ValueError(f"{name}: {detail}")
    check("embedding_shape",x.shape==(12510,320),str(x.shape))
    check("embedding_finite",np.isfinite(x).all(),str(np.isfinite(x).sum()))
    check("nonzero_embedding_norms",np.all(np.linalg.norm(x,axis=1)>0),"Required for cosine")
    check("metadata_csv_rows",len(m)==len(c)==12510,f"{len(m)}, {len(c)}")
    idx=[int(r["embedding_row"]) for r in m]
    check("embedding_row_bijection",sorted(idx)==list(range(12510)),"Unique indices 0..12509")
    m=sorted(m,key=lambda r:int(r["embedding_row"]))
    check("unique_sample_ids",len({r["ID"] for r in m})==12510,"12,510 unique IDs")
    check("metadata_csv_alignment",all(all(a[k]==b[k] for k in ("ID","labels","isBGC")) for a,b in zip(m,c)),
          "Exact ID, labels, isBGC match in embedding_row order")
    check("source_csv_identity",hashes[str(paths[3])]==hashes[str(paths[4])],"Analysis copy SHA256 equals root CSV")
    check("probability_shape_finite",p.shape==(12510,7) and np.isfinite(p).all(),str(p.shape))
    check("probability_range",np.all((p>=0)&(p<=1)),"Multilabel probabilities: no row-sum constraint")
    groups={}
    for i,r in enumerate(m):
        match=re.fullmatch(r"(BGC[0-9]{7})_(-1|0|1|2|3)",r["ID"])
        check_id=match is not None
        if not check_id: raise ValueError("Unexpected sample ID "+r["ID"])
        groups.setdefault(match[1],[]).append((int(match[2]),i))
    check("bgc_count",len(groups)==2502,str(len(groups)))
    check("five_augmentations",all(sorted(s for s,i in g)==[-1,0,1,2,3] for g in groups.values()),"Exact suffix set for each BGC")
    check("constant_bgc_metadata",all(len({(c[i]["labels"],c[i]["sentence"]) for s,i in g})==1 for g in groups.values()),"Label and original sentence constant")
    ids=sorted(groups)
    ix=np.array([[i for s,i in sorted(groups[b])] for b in ids])
    labels=[m[i[0]]["labels"] for i in ix]
    proteins=[len(c[i[0]]["sentence"].split()) for i in ix]
    table("input_checksums.tsv",[dict(path=k,sha256=v,bytes=Path(k).stat().st_size) for k,v in hashes.items()])
    return x.astype(np.float64),m,ids,ix,labels,proteins,hashes

def qc(x,m,ids,ix,labels,proteins):
    z=x/np.linalg.norm(x,axis=1,keepdims=True)
    pairs=list(itertools.combinations(range(5),2))
    cos=np.stack([np.sum(z[ix[:,a]]*z[ix[:,b]],axis=1) for a,b in pairs],axis=1).clip(-1,1)
    table("augmentation_pairwise_cosines.tsv",(
        dict(BGC=bgc,label=labels[i],sample_a=m[ix[i,a]]["ID"],sample_b=m[ix[i,b]]["ID"],
             embedding_row_a=ix[i,a],embedding_row_b=ix[i,b],cosine=cos[i,j])
        for i,bgc in enumerate(ids) for j,(a,b) in enumerate(pairs)))
    means=cos.mean(axis=1); mins=cos.min(axis=1)
    q1,q3=np.quantile(means,[.25,.75]); fence=q1-1.5*(q3-q1)
    rows=[dict(BGC=b,label=labels[i],protein_count=proteins[i],n_pairs=10,
               mean_cosine=means[i],min_cosine=mins[i],median_cosine=np.median(cos[i]),
               max_cosine=cos[i].max(),std_cosine=cos[i].std(ddof=1),
               min_raw_norm=np.linalg.norm(x[ix[i]],axis=1).min(),
               max_raw_norm=np.linalg.norm(x[ix[i]],axis=1).max(),
               below_095=bool(mins[i]<.95),relative_low_tail=bool(means[i]<fence))
          for i,b in enumerate(ids)]
    table("augmentation_stability.tsv",rows)
    ranked=sorted(rows,key=lambda r:(r["min_cosine"],r["mean_cosine"],r["BGC"]))
    table("augmentation_stability_ranked.tsv",[dict(rank=j+1,**r) for j,r in enumerate(ranked)])
    table("unstable_bgcs.tsv",[r for r in ranked if r["below_095"]],list(rows[0]))
    stats=quantiles(cos.ravel(),"all_pairwise_cosines")+quantiles(means,"bgc_mean_cosine")+quantiles(mins,"bgc_min_cosine")
    table("stability_global_distribution.tsv",stats)
    summary=dict(script_version=VERSION,n_bgc=len(ids),n_pairs=cos.size,minimum_pairwise=float(mins.min()),
                 median_pairwise=float(np.median(cos)),mean_pairwise=float(cos.mean()),
                 unstable_threshold=.95,n_below_095=int((mins<.95).sum()),
                 relative_fence=float(fence),n_relative_low_tail=int((means<fence).sum()))
    (OUT/"qc_summary.json").write_text(json.dumps(summary,indent=2))
    print("QC",json.dumps(summary),flush=True)
    return cos,rows,summary

def figures(cos,rows,coords=None,scores=None,pca=None,knn=None,any_s=None,pure_s=None,hs=None,assign=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                         "pdf.fonttype":42,"ps.fonttype":42,"savefig.dpi":600})
    F=OUT/"figures";F.mkdir(exist_ok=True)
    def save(fig,name):
        fig.tight_layout()
        for ext in ("pdf","png"): fig.savefig(F/(name+"."+ext),bbox_inches="tight")
        plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(9,3.5))
    ax[0].hist(cos.ravel(),bins=60,color="#3574A8")
    ax[0].set(xlabel="Pairwise cosine similarity",ylabel="Augmentation pairs",title="25,020 pairs within 2,502 BGCs")
    for key,label,color in (("min_cosine","BGC minimum","#C87522"),("mean_cosine","BGC mean","#3574A8")):
        vals=np.sort([r[key] for r in rows])
        ax[1].plot(vals,np.arange(1,len(vals)+1)/len(vals),label=label,color=color)
    ax[1].axvline(.95,ls="--",color="0.4",lw=1,label="QC cutoff 0.95")
    ax[1].set(xlabel="Cosine similarity",ylabel="Cumulative fraction of BGCs",title="Stability across five augmentations")
    ax[1].legend(fontsize=8)
    save(fig,"augmentation_stability")
    worst=sorted(rows,key=lambda r:r["min_cosine"])[:20][::-1]
    fig,ax=plt.subplots(figsize=(7,6))
    ax.scatter([r["min_cosine"] for r in worst],range(20),color="#C87522")
    ax.set_yticks(range(20),[r["BGC"] for r in worst])
    ax.axvline(.95,ls="--",c="0.4")
    ax.set(xlabel="Minimum of 10 pairwise cosine similarities",title="20 BGCs with the lowest padding stability")
    save(fig,"least_stable_bgcs")
    if coords is None:return
    def scatter(ax,xy,mask,title):
        ax.scatter(xy[~mask,0],xy[~mask,1],s=6,c="#C5C5C5",alpha=.5,rasterized=True,label=f"Other (n={(~mask).sum():,})")
        ax.scatter(xy[mask,0],xy[mask,1],s=13,c="#176BA0",alpha=.8,rasterized=True,label=f"Positive (n={mask.sum():,})")
        ax.set_title(title);ax.legend(fontsize=7,loc="best")
    fig,ax=plt.subplots(1,2,figsize=(9,3.6))
    scatter(ax[0],scores,any_s,"PCA: any Saccharide label")
    ax[0].set(xlabel=f"PC1 ({pca.explained_variance_ratio_[0]:.1%})",ylabel=f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    ax[1].plot(np.arange(1,321),np.cumsum(pca.explained_variance_ratio_),c="#3574A8")
    ax[1].set(xlabel="Number of principal components",ylabel="Cumulative explained variance",ylim=(0,1.01),title="Centered raw consensus embeddings")
    save(fig,"pca_qc")
    fig,ax=plt.subplots(1,2,figsize=(9,3.8))
    for a,mask,title in zip(ax,[any_s,pure_s],["Any Saccharide label","Pure Saccharide only"]):
        scatter(a,coords,mask,title);a.set(xlabel="UMAP 1",ylabel="UMAP 2")
    fig.suptitle("Consensus BGC visualization; separation is not quantitative evidence",fontsize=10)
    save(fig,"umap_saccharide")
    fig,ax=plt.subplots(1,2,figsize=(8,3.6))
    for a,definition,title in zip(ax,["any_saccharide","pure_saccharide"],["Any Saccharide label","Pure Saccharide only"]):
        vals=[r for r in knn if r["definition"]==definition]
        a.plot([r["k"] for r in vals],[r["precision_at_k"] for r in vals],"o-",c="#176BA0",label="Positive-query Precision@k")
        a.axhline(vals[0]["global_prevalence"],color="#C87522",ls="--",label="Global prevalence")
        for r in vals:a.annotate(f'{r["precision_at_k"]:.1%}',(r["k"],r["precision_at_k"]),xytext=(0,6),textcoords="offset points",ha="center",fontsize=8)
        a.set(title=title,xlabel="k (self excluded)",ylabel="Fraction positive",ylim=(0,1.05),xticks=[5,10,20]);a.legend(fontsize=7)
    save(fig,"saccharide_knn_enrichment")
    fig,ax=plt.subplots(1,2,figsize=(9,3.7))
    labs,counts=np.unique(assign,return_counts=True)
    ax[0].bar(range(len(labs)),counts,color=["#AAAAAA" if l<0 else "#3574A8" for l in labs])
    ax[0].set_xticks(range(len(labs)),["Noise" if l<0 else str(l) for l in labs],rotation=90)
    ax[0].set(ylabel="BGC count",xlabel="Exploratory cluster",title="HDBSCAN: size 30, samples 10")
    grid=np.array([r["noise_fraction"] for r in hs]).reshape(3,3)
    im=ax[1].imshow(grid,vmin=0,vmax=1,cmap="Blues")
    ax[1].set(xticks=range(3),xticklabels=[5,10,20],yticks=range(3),yticklabels=[15,30,60],xlabel="min_samples",ylabel="min_cluster_size",title="Noise fraction / cluster count")
    for i in range(3):
        for j in range(3):
            r=hs[i*3+j];ax[1].text(j,i,f'{r["noise_fraction"]:.0%}\n{r["n_clusters"]} clusters',ha="center",va="center",color="white" if grid[i,j]>.5 else "black",fontsize=8)
    fig.colorbar(im,ax=ax[1],label="Noise fraction")
    save(fig,"hdbscan_exploration")

def downstream(x,ids,ix,labels,proteins,cos,rows):
    from sklearn.decomposition import PCA
    from sklearn.cluster import HDBSCAN
    import umap
    consensus=x[ix].mean(axis=1)
    if consensus.shape!=(2502,320) or not np.isfinite(consensus).all() or np.any(np.linalg.norm(consensus,axis=1)==0):
        raise ValueError("Invalid consensus")
    matrix("consensus_embeddings.tsv",consensus,ids,"dim_")
    np.save(OUT/"consensus_embeddings.npy",consensus,allow_pickle=False)
    any_s=np.array(["Saccharide" in label for label in labels])
    pure_s=np.array([label.strip()=="Saccharide" for label in labels])
    table("consensus_metadata.tsv",[dict(BGC=b,label=labels[i],any_saccharide=bool(any_s[i]),pure_saccharide=bool(pure_s[i]),n_augmentations=5,protein_count=proteins[i],padding_qc_flag=rows[i]["below_095"]) for i,b in enumerate(ids)])
    pca=PCA(n_components=320,svd_solver="full",whiten=False)
    scores=pca.fit_transform(consensus)
    matrix("pca_scores.tsv",scores,ids,"PC")
    matrix("pca_components.tsv",pca.components_,[f"PC{i+1}" for i in range(320)],"dim_")
    table("pca_feature_means.tsv",[dict(dimension=i+1,mean=v) for i,v in enumerate(pca.mean_)])
    table("pca_explained_variance.tsv",[dict(PC=i+1,variance=v,ratio=r,cumulative=c) for i,(v,r,c) in enumerate(zip(pca.explained_variance_,pca.explained_variance_ratio_,np.cumsum(pca.explained_variance_ratio_)))])
    print("Consensus and PCA saved; fitting UMAP",flush=True)
    coords=umap.UMAP(n_components=2,n_neighbors=15,min_dist=.1,metric="cosine",random_state=SEED,transform_seed=SEED,n_jobs=1,n_epochs=500).fit_transform(consensus)
    matrix("umap_coordinates.tsv",coords,ids,"UMAP")
    z=consensus/np.linalg.norm(consensus,axis=1,keepdims=True)
    distance=np.clip(1-z@z.T,0,2);distance=(distance+distance.T)/2;np.fill_diagonal(distance,0)
    search=distance.copy();np.fill_diagonal(search,np.inf)
    # Stable sorting of lexicographically ordered BGCs gives deterministic distance ties.
    neighbors=np.argsort(search,axis=1,kind="stable")[:,:20]
    table("cosine_neighbors_top20.tsv",(dict(BGC=ids[i],rank=j+1,neighbor=ids[n],cosine_distance=distance[i,n],neighbor_label=labels[n],any_saccharide=bool(any_s[n]),pure_saccharide=bool(pure_s[n])) for i in range(len(ids)) for j,n in enumerate(neighbors[i])))
    knn=[];per_query=[]
    for definition,mask in (("any_saccharide",any_s),("pure_saccharide",pure_s)):
        count=int(mask.sum());prev=count/len(ids)
        for k in (5,10,20):
            hit=mask[neighbors[mask,:k]].sum(axis=1);prec=float(np.mean(hit/k))
            knn.append(dict(definition=definition,k=k,n_positive=count,n_total=len(ids),precision_at_k=prec,global_prevalence=prev,fold_over_prevalence=prec/prev,excess_percentage_points=100*(prec-prev),self_excluded_null_prevalence=(count-1)/(len(ids)-1)))
            per_query.extend(dict(definition=definition,BGC=ids[i],k=k,positive_neighbors=int(h),precision_at_k=h/k) for i,h in zip(np.flatnonzero(mask),hit))
    table("saccharide_knn_summary.tsv",knn);table("saccharide_knn_per_bgc.tsv",per_query)
    print("kNN",json.dumps(knn),flush=True)
    hs=[];allassign=[];primary=None;probs=None
    for size in (15,30,60):
        for samples in (5,10,20):
            model=HDBSCAN(copy=True,min_cluster_size=size,min_samples=samples,metric="precomputed",algorithm="brute",n_jobs=1,cluster_selection_method="eom",allow_single_cluster=False)
            before_distance=sha_array(distance)
            assignment=model.fit_predict(distance)
            if sha_array(distance)!=before_distance: raise RuntimeError("HDBSCAN mutated original distances")
            hs.append(dict(min_cluster_size=size,min_samples=samples,n_clusters=len(set(assignment)-{-1}),n_noise=int((assignment==-1).sum()),noise_fraction=float((assignment==-1).mean())))
            allassign.extend(dict(BGC=b,min_cluster_size=size,min_samples=samples,cluster=int(assignment[i]),membership_strength=float(model.probabilities_[i])) for i,b in enumerate(ids))
            if size==30 and samples==10:primary=assignment.copy();probs=model.probabilities_.copy()
    table("hdbscan_sensitivity.tsv",hs);table("hdbscan_all_assignments.tsv",allassign)
    table("hdbscan_primary_assignments.tsv",[dict(BGC=b,label=labels[i],cluster=int(primary[i]),membership_strength=probs[i]) for i,b in enumerate(ids)])
    table("hdbscan_primary_composition.tsv",[dict(cluster=int(cluster),label=label,n=int(sum((primary==cluster)&(np.array(labels)==label)))) for cluster in sorted(set(primary)) for label in sorted(set(labels)) if sum((primary==cluster)&(np.array(labels)==label))])
    figures(cos,rows,coords,scores,pca,knn,any_s,pure_s,hs,primary)
    summary=dict(consensus_shape=list(consensus.shape),pca_pc1=float(pca.explained_variance_ratio_[0]),pca_pc2=float(pca.explained_variance_ratio_[1]),
                 pca_components_90=int(np.searchsorted(np.cumsum(pca.explained_variance_ratio_),.9)+1),
                 hdbscan_primary=next(r for r in hs if r["min_cluster_size"]==30 and r["min_samples"]==10),knn=knn)
    (OUT/"analysis_summary.json").write_text(json.dumps(summary,indent=2))
    print("DOWNSTREAM",json.dumps(summary),flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage",choices=["qc","downstream","all"],default="all")
    parser.add_argument("--data-dir",type=Path,default=ROOT/"data")
    parser.add_argument("--output-dir",type=Path,default=ROOT/"results")
    args=parser.parse_args();DATA=args.data_dir;OUT=args.output_dir;OUT.mkdir(parents=True,exist_ok=True)
    x,m,ids,ix,labels,proteins,hashes=validate()
    cos,rows,summary=qc(x,m,ids,ix,labels,proteins)
    if args.stage=="qc":figures(cos,rows)
    else:downstream(x,ids,ix,labels,proteins,cos,rows)
    after={p:sha(Path(p)) for p in hashes}
    table("input_integrity_after.tsv",[dict(path=p,sha256_before=h,sha256_after=after[p],unchanged=h==after[p]) for p,h in hashes.items()])
    if after!=hashes:raise RuntimeError("Input changed during analysis")
    versions={d.metadata["Name"]:d.version for d in metadata.distributions(path=[str(ROOT/".deps")])}
    table("package_versions.tsv",[dict(package=k,version=v) for k,v in sorted(versions.items())])
    (ROOT/"requirements.lock.txt").write_text("\n".join(f"{k}=={v}" for k,v in sorted(versions.items()))+"\n")
    provenance=dict(script_version=VERSION,script_sha256=sha(Path(__file__)),utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    python=sys.version,platform=platform.platform(),command=sys.argv,seed=SEED,stage=args.stage,
                    umap=dict(n_neighbors=15,min_dist=.1,metric="cosine",n_epochs=500,n_jobs=1,init="spectral"),
                    consensus="Arithmetic mean of five raw float64-converted embeddings",knn="Original 320-D cosine; self excluded; stable lexicographic ties",
                    hdbscan="sklearn; precomputed 320-D cosine; primary size30/samples10; EOM",input_hashes=hashes,versions=versions)
    (OUT/("provenance_"+args.stage+".json")).write_text(json.dumps(provenance,indent=2))
    print("DONE; original input checksums unchanged",flush=True)

