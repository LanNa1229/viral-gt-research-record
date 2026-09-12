"""Independently check saved numerical outputs and complete the analysis report."""
import sys, json, csv, hashlib, datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/".deps"))
import numpy as np
from scipy.spatial.distance import pdist, cdist
R=ROOT/"results"
def read(name):
    with open(R/name,encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write(name,rows):
    with open(R/name,"w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
def matrix(name):
    rows=read(name);return [r["BGC"] for r in rows],np.array([[float(v) for k,v in r.items() if k!="BGC"] for r in rows])
checks=[]
def check(name,ok,detail=""):
    checks.append(dict(check=name,passed=bool(ok),detail=detail));write("independent_verification.tsv",checks)
    if not ok:raise ValueError(name)
x=np.load(ROOT/"data/mibig_12510_embeddings.npy",allow_pickle=False).astype(float)
with open(ROOT/"data/mibig_12510_metadata.tsv",encoding="utf-8-sig",newline="") as f:metadata=list(csv.DictReader(f,delimiter="\t"))
groups={}
for r in metadata:groups.setdefault(r["ID"].rsplit("_",1)[0],[]).append((int(r["ID"].rsplit("_",1)[1]),int(r["embedding_row"])))
ids,y=matrix("consensus_embeddings.tsv")
expected=np.array([np.mean(x[[i for _,i in sorted(groups[b])]],axis=0) for b in ids])
check("raw_arithmetic_means",np.array_equal(y,expected),"All 2502 x 320 saved values equal recomputed float64 means")
check("npy_tsv_identical",np.array_equal(y,np.load(R/"consensus_embeddings.npy",allow_pickle=False)))
pairs=read("augmentation_pairwise_cosines.tsv");observed=np.array([float(r["cosine"]) for r in pairs])
independent=np.concatenate([1-pdist(x[[i for _,i in sorted(groups[b])]],metric="cosine") for b in ids])
check("all_pairwise_scipy_validation",len(pairs)==25020 and np.allclose(observed,independent,rtol=0,atol=1e-12),f"Max absolute difference {abs(observed-independent).max():.3g}")
neighbors=read("cosine_neighbors_top20.tsv");d=cdist(y,y,metric="cosine");np.fill_diagonal(d,np.inf)
expected_neighbors=np.argsort(d,axis=1,kind="stable")[:,:20]
lookup={b:i for i,b in enumerate(ids)}
actual=np.array([lookup[r["neighbor"]] for r in neighbors]).reshape(2502,20)
check("all_neighbors_independent_cosine",np.array_equal(actual,expected_neighbors),"SciPy cdist, self excluded, stable ties")
check("neighbor_distances",np.allclose(np.array([float(r["cosine_distance"]) for r in neighbors]).reshape(2502,20),np.take_along_axis(d,actual,axis=1),atol=1e-12,rtol=0))
cm=read("consensus_metadata.tsv")
for definition in ("any_saccharide","pure_saccharide"):
    mask=np.array([r[definition]=="True" for r in cm])
    for row in read("saccharide_knn_summary.tsv"):
        if row["definition"]==definition:
            k=int(row["k"]);v=mask[expected_neighbors[mask,:k]].mean()
            check(f"{definition}_precision_{k}",abs(v-float(row["precision_at_k"]))<1e-14)
pcids,scores=matrix("pca_scores.tsv");_,components=matrix("pca_components.tsv")
mean=np.array([float(r["mean"]) for r in read("pca_feature_means.tsv")])
check("pca_reconstruction",pcids==ids and np.allclose(scores@components+mean,y,atol=1e-10,rtol=0))
uid,coords=matrix("umap_coordinates.tsv")
check("umap_consensus_rows",uid==ids and coords.shape==(2502,2) and np.isfinite(coords).all())
check("hdbscan_consensus_rows",[r["BGC"] for r in read("hdbscan_primary_assignments.tsv")]==ids)
check("hdbscan_sensitivity_rows",len(read("hdbscan_all_assignments.tsv"))==9*2502)
for r in read("input_integrity_after.tsv"):
    h=hashlib.sha256(Path(r["path"]).read_bytes()).hexdigest()
    check("input_unchanged:"+Path(r["path"]).name,h==r["sha256_before"]==r["sha256_after"])
from PIL import Image
for p in sorted((R/"figures").glob("*.png")):
    with Image.open(p) as im:im.verify()
check("figure_pairs",len(list((R/"figures").glob("*.png")))==6 and len(list((R/"figures").glob("*.pdf")))==6)
qc=json.loads((R/"qc_summary.json").read_text());analysis=json.loads((R/"analysis_summary.json").read_text())
write("qc_summary.tsv",[dict(metric=k,value=v) for k,v in qc.items()])
write("analysis_summary.tsv",[dict(metric=k,value=v) for k,v in analysis.items() if not isinstance(v,(dict,list))]+[dict(metric="hdbscan_"+k,value=v) for k,v in analysis["hdbscan_primary"].items()])
cohorts=[]
for definition,mask in [("all",np.ones(2502,dtype=bool)),("any_saccharide",np.array([r["any_saccharide"]=="True" for r in cm])),("pure_saccharide",np.array([r["pure_saccharide"]=="True" for r in cm]))]:
    flags=np.array([r["padding_qc_flag"]=="True" for r in cm])
    cohorts.append(dict(cohort=definition,n=int(mask.sum()),n_min_cosine_below_095=int((mask&flags).sum()),fraction_flagged=float(flags[mask].mean())))
write("stability_cohort_summary.tsv",cohorts)
h=analysis["hdbscan_primary"]
paragraph=f"""
All 13 structural checks passed: 12,510 finite 320-D embeddings, exact metadata/CSV alignment, and 2,502 BGCs with five augmentations each. Original input SHA256 values remained unchanged. Independent numerical verification passed for all within-BGC cosines, all cosine neighbors, consensus arithmetic, PCA reconstruction, and saved coordinate alignment.

Among 25,020 pairwise cosines, the median was {qc['median_pairwise']:.6f}, the mean {qc['mean_pairwise']:.6f}, and the minimum {qc['minimum_pairwise']:.6f}. The minimum-pair cutoff of 0.95 flagged {qc['n_below_095']:,}/{qc['n_bgc']:,} BGCs ({qc['n_below_095']/qc['n_bgc']:.1%}). The relative mean-cosine cutoff was {qc['relative_fence']:.6f}, flagging {qc['n_relative_low_tail']} BGCs. Padding sensitivity is substantial despite a high overall median.

The five lowest minimum pairwise cosines were:

| BGC | Label | Minimum cosine |
|---|---|---:|
"""
for r in read("augmentation_stability_ranked.tsv")[:5]:paragraph+=f"| {r['BGC']} | {r['label']} | {float(r['min_cosine']):.6f} |\n"
paragraph+=f"""
All 2,502 BGCs were retained in the (2502,320) consensus. PC1 explained {analysis['pca_pc1']:.1%} and PC2 {analysis['pca_pc2']:.1%} of variance; {analysis['pca_components_90']} components reached 90% cumulative variance.

| Saccharide definition | Positive / total | Global prevalence | k | Precision@k | Fold / prevalence |
|---|---:|---:|---:|---:|---:|
"""
for r in analysis["knn"]:paragraph+=f"| {r['definition']} | {r['n_positive']}/{r['n_total']} | {r['global_prevalence']:.2%} | {r['k']} | {r['precision_at_k']:.2%} | {r['fold_over_prevalence']:.2f} |\n"
paragraph+=f"""
Primary HDBSCAN found {h['n_clusters']} clusters and classified {h['n_noise']} BGCs ({h['noise_fraction']:.1%}) as noise. The nine-setting sensitivity table accompanies the exploratory assignments. These results do not establish biological classes.

| Padding QC cohort | BGCs | Flagged | Fraction |
|---|---:|---:|---:|
"""
for r in cohorts:paragraph+=f"| {r['cohort']} | {r['n']} | {r['n_min_cosine_below_095']} | {r['fraction_flagged']:.1%} |\n"
paragraph+="""
This is a descriptive reference analysis. No p-values or held-out performance estimates are claimed. High Saccharide neighborhood enrichment must be considered together with supervised training context and observed augmentation sensitivity.

Run the following after the main pipeline to reproduce independent checks, cohort QC summary, and this results section:

    python -B scripts/verify_and_summarize.py

The full execution arguments and versions are in results/provenance_downstream.json. Six figures were exported as PDF and PNG.
"""
readme=ROOT/"README.md"
readme.write_text(readme.read_text(encoding="utf-8-sig").split("## Results")[0]+"## Results\n"+paragraph,encoding="utf-8")
(R/"verification_provenance.json").write_text(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,command=sys.argv),indent=2))
print(json.dumps(dict(checks=len(checks),passed=True,cohorts=cohorts),indent=2))

