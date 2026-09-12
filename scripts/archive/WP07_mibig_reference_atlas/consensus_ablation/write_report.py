from common import *
def md(df):
 def fmt(v):
  if isinstance(v,(float,np.floating)):return f"{v:.6g}"
  return str(v).replace("|","/")
 return "| "+" | ".join(df.columns)+" |\n| "+" | ".join(["---"]*len(df.columns))+" |\n"+"\n".join("| "+" | ".join(fmt(v) for v in row)+" |" for row in df.itertuples(index=False,name=None))+"\n"
def main():
 I=json.loads((ROOT/"provenance/inputs.json").read_text());m,arrays=load()
 s=read(ROOT/"tables/A_vs_B_saccharide_knn_summary.tsv")
 ov=read(ROOT/"tables/neighbor_overlap_summary.tsv");pc=read(ROOT/"tables/pca_explained_variance.tsv")
 th=read(ROOT/"tables/pca_variance_thresholds.tsv")
 v=json.loads((ROOT/"results/final_validation.json").read_text())
 ref=json.loads((ROOT/"results/A_reference_geometry_check.json").read_text())
 brief=[]
 for r in s.itertuples():
  brief.append(dict(definition=r.definition,k=r.k,A_precision_percent=100*r.A_precision,B_precision_percent=100*r.B_precision,
   B_minus_A_pp=100*r.B_minus_A_precision,CI95_low_pp=100*r.paired_bootstrap_CI_low,CI95_high_pp=100*r.paired_bootstrap_CI_high,
   A_fold=r.A_fold_enrichment,B_fold=r.B_fold_enrichment))
 brief=pd.DataFrame(brief);tab("tables/readme_comparison_percentages.tsv",brief)
 fraction=s[["definition","k","median_paired_difference","fraction_queries_B_better","fraction_queries_equal","fraction_queries_B_worse"]]
 comparison=read(ROOT/"tables/consensus_vs_superseded_zero.tsv")
 inputrows=[dict(role=role,path=I[key],sha256=sha(I[key])) for role,key in [("A array","published"),("A metadata","published_metadata"),("B array","core"),("B metadata","core_metadata"),("Training CSV","csv")]]
 tab("provenance/primary_inputs.tsv",inputrows)
 count=integrity()
 source=read(ROOT/"provenance/figure_sources.tsv") if False else pd.read_csv(ROOT/"provenance/figure_sources.tsv",sep="\t",header=None,names=["figure","source_table"])
 assert len(source)==14 and source.figure.is_unique
 tab("provenance/visual_qa.tsv",[dict(figure=r.figure,PNG_visually_reviewed=True,labels_readable=True,source_table=r.source_table,
  note="Matched conventions checked; independently fitted axes not interpreted as aligned") for r in source.itertuples()])
 distributions=sorted((d.metadata["Name"],d.version) for d in importlib.metadata.distributions(path=[str(PREV/".deps")]))
 tab("provenance/package_versions.tsv",[dict(package=p,version=x) for p,x in distributions])
 (ROOT/"provenance/requirements.lock.txt").write_text("\n".join(f"{p}=={x}" for p,x in distributions)+"\n")
 text=f"""# Corrected consensus Step 3: published-style vs core-context

## Inputs and consensus construction

This analysis supersedes the earlier suffix-zero Step 3. Step 2 was completed previously and was not rerun. Every augmentation suffix is an augmentation state; none is selected as a canonical/original BGC representation.

The A path was recovered from the previous analysis provenance and checked against the SHA256 in the validated full extraction's input-integrity record. That extraction reports 12,510/12,510 exact unmasked-mean matches to A, maximum absolute difference zero. No inference was rerun.

{md(pd.DataFrame(inputrows))}

For each base BGC, IDs are parsed with the anchored pattern BGC followed by seven digits and a suffix in -1,0,1,2,3. Groups are lexicographically ordered by base_BGC; their five source vectors are ordered by suffix for deterministic accumulation. Raw vectors are averaged independently for A and B using float64 accumulation, with no individual-vector normalization. Both saved matrices are float64 (2502,320). Transient copies alone are normalized for cosine.

The metadata preserve original trimmed multi-label strings, five source IDs and row indices as JSON lists, suffixes, plotting class and both Saccharide flags. The seven exact single-label categories plus Hybrid are used for plotting; Hybrid means more than one whitespace-separated original class token. All tokens were verified against the seven allowed classes. Any Saccharide means Saccharide token membership (195/2502, 7.793765%); pure Saccharide means exact trimmed label Saccharide (130/2502, 5.195843%). The 65 Saccharide-positive hybrids are not assigned to pure Saccharide.

## Validation

Both input arrays are finite float32 (12510,320). A/B metadata match exactly, including row IDs, labels and embedding-row indices, and align with the original CSV. All 2,502 groups have five rows, the exact suffix set, and constant labels. Both consensus arrays are finite with nonzero norms.

The A reference gate passed before interpreting A vs B: all six saved reference Precision@5/10/20 values reproduced within 1e-12. Independent comparison then found the entire A consensus array exactly identical to the previous float64 consensus (maximum absolute difference {ref['consensus_max_abs_difference']}) and all prior top-20 neighbor sets identical for every BGC (minimum and mean overlap 100%). This verifies the old consensus geometry; it does not assert old/new UMAP coordinate identity.

{v['checks']:,} independent checks passed, covering all source-to-consensus mappings, explicit summation of five raw vectors, all top-50 neighbors via independent SciPy cosine distances, neighbor-table identities/distances, full PCA reconstruction, exact coordinate/metadata joins, all bootstrap mean intervals, regeneration of all 10,000 common permutation assignments, and independent recomputation of 100 selected permutation statistics for both representations and all endpoints/k. All empirical p-values were checked against the full saved null arrays. PNG readability and paired PDF existence were verified; all 14 figures were visually reviewed.

SHA256 was rerun on {count} protected scientific files across the original extraction and both prior analysis directories, plus the CSV. All are unchanged. Dependency/cache/temporary directories are excluded from the scientific-file manifest; existing dependencies are reused read-only with bytecode writing disabled and cache destinations redirected into this new directory.

## Saccharide kNN comparison

Neighbors use the original 320-D consensus spaces, identical float64 cosine implementation and all 2,502 BGC candidates, excluding the query. Stable sorting breaks exact distance ties by lexicographic base_BGC; there were no boundary ties at any evaluated k. Exact top-50 identities and distances are saved for A and B.

Precision@k is the mean positive-neighbor fraction among positive queries. Fold enrichment divides by global prevalence M/2502. The finite-population self-excluded random expectation is (M-1)/2501: 0.0775689724 for any Saccharide and 0.0515793683 for pure. Percentage-point enrichment is 100 times Precision minus prevalence.

The primary comparison below uses consensus values only. CI columns are percentile 95% intervals for the mean paired B-minus-A difference, based on 10,000 positive-query BGC resamples (seed 1730), preserving A/B pairing and fixed neighbor graphs. Differences and CIs are in percentage points.

{md(brief)}

The full primary table is [A_vs_B_saccharide_knn_summary.tsv](tables/A_vs_B_saccharide_knn_summary.tsv), including null means, p-values, z-scores, median paired differences and their bootstrap intervals, and all better/equal/worse fractions. Mean losses need not imply a negative median because many queries tie, particularly at small k.

{md(fraction)}

A single permutation of the joint two-column Saccharide membership tuple is used for both A and B at each of 10,000 iterations (seed 1729). Geometry remains fixed; positive-query identities are recomputed after permutation. The full shared permutation assignments and null statistics are saved. The one-sided localization p-value is (1 + number of null precisions at least as large as observed)/(10001). All 16 endpoint/representation/k tests attain p=1/10001, approximately 0.00009999. These are unadjusted p-values (even a 16-test Bonferroni bound is below 0.0016). Z-scores standardize against empirical null mean/SD and are descriptive, not normal-tail p-values. Permutations test localization relative to random labels, not the paired A-vs-B effect or held-out prediction.

**Observation:** B has lower Precision@k for both definitions at every k. **Statistical result:** all eight paired mean bootstrap intervals lie below zero. The losses are larger for pure Saccharide and generally larger at broad neighborhoods. Both representations still exhibit substantial enrichment over random labels.

![Precision comparison](figures/saccharide_knn_comparison.png)
![Paired differences](figures/saccharide_paired_difference.png)

## Neighborhood overlap

For each BGC, overlap is the intersection of A and B neighbor sets divided by k. These comparisons are invariant to PCA axis signs/rotations and UMAP orientations.

{md(ov[["group","k","n_BGC","mean_overlap","median_overlap"]])}

At k=20, overall mean overlap is 34.142686% and median 35%; any-Saccharide mean is 33.769231%, pure-Saccharide mean 27.5%, and non-any-Saccharide mean 34.174252%. Broad class enrichment can persist despite substantial turnover of individual neighbors.

![Overlap](figures/neighbor_overlap.png)

## PCA comparison

PCA is independently fit to each raw, mean-centered (2502,320) consensus matrix with full SVD, 320 components and no whitening or dimension standardization. Learned-coordinate dimension scales are saved in tables/dimension_scale_diagnostics.tsv; no diagnostic justified changing these scales. Full scores, centers, components, variances and cumulative fractions are exported. The first 20 PCs are shown in scree plots; all 320 cumulative values are saved.

{md(pc[pc.PC<=3][["representation","PC","explained_ratio","cumulative_ratio"]])}

{md(th)}

Both require five PCs for 90% of variance, with B requiring one additional PC for 80% and 95%. Matched plotting uses identical point size, alpha, class colors and legend conventions; PC1/PC2 limits span both fits. Independently fitted PC axes have arbitrary sign/orientation and are not numerically equivalent.

**Visual observation:** both PC1/PC2 plots retain broad class-associated arms and a dominant concentration of pure-Saccharide points, while any-Saccharide includes points in additional regions. The plots are compatible with retained broad structure, but they do not resolve the local 320-D loss quantitatively and neither overturn nor independently establish the kNN result.

![A PCA](figures/pca_A_pc1_pc2_classes.png)
![B PCA](figures/pca_B_pc1_pc2_classes.png)

## UMAP comparison

A and B are independently fit using cosine metric, n_neighbors=30, min_dist=0.1, n_components=2, random_state=transform_seed=1729, 500 epochs, spectral initialization and n_jobs=1. A limited matched robustness panel uses n_neighbors=15,30,50 for both, with all other settings fixed. No broad search or separate visual optimization was performed. All six coordinate tables and exact parameters are saved. Each fit is based on the 2,502 consensus BGCs.

**Visual observation:** both primary plots contain a concentration of pure-Saccharide points. In B, this concentration appears more continuous with regions containing other labels, while A displays a more distinct Saccharide-rich end region. Additional any-Saccharide locations remain visible in both. This is qualitatively compatible with weaker localization in B; it is not quantitative evidence or a demonstration of new biological classes. The limited panel retains a Saccharide-rich region while positions, shapes and connections change. Raw axis orientation and between-island distances are not compared as biology.

![A UMAP](figures/umap_A_saccharide.png)
![B UMAP](figures/umap_B_saccharide.png)
![Limited robustness](figures/umap_robustness.png)

## Corrected Step 3 conclusion

1. **Reference reproduction:** Yes. A consensus is exactly identical to the previous float64 consensus, reproduces all six saved Saccharide precision values, and has identical previous top-20 neighbor sets for every BGC.
2. **Exact consensus precisions:** All eight A/B comparisons at k=5,10,20,50 are in the primary table above, with full numerical precision in the TSV.
3. **Localization loss/gain:** B loses 6.26–9.16 percentage points for any Saccharide and 10.62–18.03 points for pure Saccharide across these k values. At k=20 the losses are 7.333333 and 13.307692 points, respectively; both paired mean CIs exclude zero.
4. **Definition consistency:** The direction is consistent at all k for both definitions, with larger losses for pure Saccharide. Per-query equality and exceptions remain and are explicitly tabulated.
5. **Neighbor identities:** Overall mean overlap rises from 21.678657% at k=5 to 45.821743% at k=50; it is 34.142686% at k=20. Fine neighbor identity is therefore only partially preserved.
6. **Variance dimensionality:** Similar low-dimensional variance concentration remains; both need five PCs for 90%, while A/B need five/six for 95%.
7. **PCA:** Broad class structure persists visually. PC1/PC2 does not contradict the kNN comparison, but is insufficient to establish local preservation or loss.
8. **UMAP:** Visual patterns are qualitatively compatible with retained Saccharide-rich regions and more mixing in B, but do not prove clustering or measure the A/B effect.
9. **Preservation judgment:** Core-context consensus preserves aspects of the trained representation's broad geometry and Saccharide enrichment, but does not preserve its local neighbor identities or Saccharide localization to the same degree as published-style consensus. It is a complementary representation ablation, not an interchangeable replacement demonstrated by these results.

**Hypothesis, not tested here:** full-position pooling may retain class-related context that is reduced by core pooling. This saved-embedding comparison does not establish the mechanism or biological causality.

## Comparison with previous `_0` analysis

The previous Step 3 selected one suffix-zero augmentation state. That was not the requested five-augmentation biological consensus and is superseded by this analysis. None of the five suffixes is treated as canonical/original here. The previous Step 2 results and files remain unchanged.

The table below explicitly separates the superseded and corrected results (fractions, not percentages). It is a historical comparison only and is not mixed into the primary consensus effect estimates.

{md(comparison[comparison.k==20])}

The qualitative direction remains B below A, but the consensus loss is smaller than the suffix-zero comparison. This correction changes the magnitude and must be carried through any later interpretation. The previous approximate A-consensus reference expectations are independently confirmed, not forced.

## Caveats

These are representations of the model's training dataset. Localization is not held-out predictive performance or biological validation. Evolutionary relationships and shared neighbor sets can create dependence between BGCs; BGC-bootstrap intervals are descriptive sampling intervals over observed positive queries conditional on the saved graphs, not retraining uncertainty or phylogenetically adjusted inference. Label permutations do not condition on taxonomy, gene count or training supervision. Two Saccharide definitions are overlapping endpoints, not independent tests. No biological function is inferred from proximity.

Consensus averaging changes representation scale and smooths augmentation variation; no claim about that mechanism is derived by rerunning Step 2. Individual vectors are never normalized before averaging. PCA is raw/centered, whereas cosine ignores vector norms transiently. UMAP is stochastic and software/hardware changes can affect layouts even with recorded seeds; package versions and coordinates are retained. Nonfatal scikit-learn deprecation warnings from UMAP are captured in the log.

## Deliverables

New output directory: {ROOT}

- results/: two float64 consensus arrays; exact consensus metadata; full PCA scores/components/centers; six UMAP coordinate tables; nearest-neighbor indices; shared permutations and null arrays; validation/reference summaries.
- tables/: primary A/B Saccharide summary, every positive-query fraction, all exact nearest-neighbor lists, overlap per BGC and group, full permutation null statistics, A-reference reproduction, PCA variance/scales, and explicitly separated earlier suffix-zero comparison.
- figures/: 14 publication-quality PNG files at 400 dpi and corresponding PDFs, including all requested filenames and matched any/pure Saccharide panels. Every figure maps to source TSVs in provenance/figure_sources.tsv.
- provenance/: exact inputs/hashes, before/after integrity, parameters, Python/packages, requirements lock, script hashes, numerical validation and visual QA.
- logs/: commands' stdout/stderr. scripts/: common.py, construct_consensus.py, compare_knn.py, geometry.py, validate_results.py, write_report.py, preview_figures.py and run_analysis.sh.

All {count} protected scientific files remain unchanged. No Step 2 rerun, extraction, model modification, HDBSCAN, classifier training, fine-tuning or viral placement was performed.

## Reproduction command

Run in WSL Ubuntu:

    cd {ROOT}
    bash scripts/run_analysis.sh

The launcher runs only the new consensus construction, kNN, PCA/UMAP, independent verification and report stages; it does not invoke the previous Step 2. It overwrites this new analysis's own generated products on a rerun. The existing dependencies at {PREV/'.deps'} are imported read-only, with writes and caches confined to this analysis directory. Their full package versions are in provenance/package_versions.tsv; core parameters are in provenance/parameters.json. No package installation was needed.

Primary random seed: 1729 (UMAP and common label permutations). Paired bootstrap seed: 1730. Both counts: 10,000. The report/TSVs, retained shared permutation assignments and script checksums allow exact audit of grouping and paired comparisons.
"""
 (ROOT/"README.md").write_text(text)
 tab("provenance/script_sha256.tsv",[dict(path=str(p),sha256=sha(p)) for p in sorted((ROOT/"scripts").iterdir()) if p.is_file()])
 js("provenance/completion.json",dict(utc=pd.Timestamp.now(tz="UTC").isoformat(),passed=True,protected_files=count,figures=14,validation_checks=v["checks"]))
 print("README complete; 14 figures, all protected input hashes unchanged.",flush=True)
if __name__=="__main__":main()
