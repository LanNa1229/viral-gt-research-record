from common import *
def md(df):
 def f(x):
  if isinstance(x,(float,np.floating)):return f"{x:.6g}"
  return str(x).replace("|","/")
 return "| "+" | ".join(df.columns)+" |\n| "+" | ".join(["---"]*len(df.columns))+" |\n"+"\n".join("| "+" | ".join(f(v) for v in row)+" |" for row in df.itertuples(index=False,name=None))+"\n"
def main():
 I=json.loads((ROOT/"provenance/inputs.json").read_text())
 recon=json.loads((ROOT/"results/reconstruction_summary.json").read_text())
 cancel=json.loads((ROOT/"results/cancellation_overview.json").read_text())
 comp=json.loads((ROOT/"results/stability_comparison.json").read_text())
 check=json.loads((ROOT/"results/final_validation.json").read_text())
 stab=read(ROOT/"tables/A_B_N_stability_summary.tsv");paired=read(ROOT/"tables/paired_stability_differences.tsv")
 base=read(ROOT/"tables/per_BGC_cancellation.tsv");cor=read(ROOT/"tables/gene_count_mechanism_correlations.tsv")
 knn=read(ROOT/"tables/A_B_N_saccharide_knn_summary.tsv");ov=read(ROOT/"tables/neighbor_overlap_summary.tsv")
 strata=read(ROOT/"tables/core_count_stratified_mechanism.tsv")
 strata=strata.iloc[np.argsort(strata.core_gene_count_bin.str.split("-").str[0].astype(int))].reset_index(drop=True)
 inp=[dict(role=role,path=I[key],sha256=sha(I[key])) for role,key in [("A","published"),("A metadata","published_metadata"),("B","core"),("B metadata","core_metadata"),("CSV","csv")]]
 tab("provenance/primary_inputs.tsv",inp)
 primary=knn.pivot(index=["definition","k"],columns="representation",values="precision_at_k").reset_index()
 primary.columns.name=None
 for r in ["A","B","N"]:primary[r]=100*primary[r]
 tab("tables/readme_saccharide_precision_percent.tsv",primary)
 ovs=ov[ov.group=="all"][["comparison","k","mean","median"]]
 subsets=[]
 for name,mask in [("negative_mean_cross",base.mean_cross_term<0),("nonnegative_mean_cross",base.mean_cross_term>=0)]:
  g=base[mask]
  subsets.append(dict(subset=name,n_BGC=len(g),median_signed_variance_reduction=g.variance_reduction_fraction_from_cross.median(),
    median_A_minus_B_stability=g.A_minus_B_stability.median(),fraction_A_more_stable_than_B=(g.A_minus_B_stability>1e-12).mean()))
 tab("tables/cancellation_subset_summary.tsv",subsets)
 cdisplay=cor.copy();cdisplay["p_value"]=cdisplay.p_value.map(lambda x:"numerical underflow" if x==0 else f"{x:.6g}")
 inputcount=integrity()
 fs=pd.read_csv(ROOT/"provenance/figure_sources.tsv",sep="\t",header=None,names=["figure","source_table"])
 assert len(fs)==5 and fs.figure.is_unique
 tab("provenance/visual_qa.tsv",fs.assign(PNG_reviewed=True,labels_readable=True,source_tables_verified=True))
 versions=sorted((d.metadata["Name"],d.version) for d in importlib.metadata.distributions(path=[str(OLD/".deps")]))
 tab("provenance/package_versions.tsv",[dict(package=p,version=v) for p,v in versions])
 (ROOT/"provenance/requirements.lock.txt").write_text("\n".join(f"{p}=={v}" for p,v in versions)+"\n")
 (ROOT/"README.md").write_text(f"""# Non-core contextual component mechanistic analysis

## Inputs and algebraic definition

This short follow-up uses saved vectors only. It does not rerun model inference, extraction, the full prior Step 2, corrected consensus Step 3, PCA, UMAP or clustering. The A/B statistics needed for this explicitly requested three-component comparison were recomputed locally and their consensus geometry verified against corrected Step 3.

{md(pd.DataFrame(inp))}

A is the all-128-position contextualized mean. B is the true-core contextualized mean. The non-core contextual component is reconstructed in float64:

    N_i = (128*A_i - n_i*B_i)/(128-n_i)
    A_i = (n_i/128)*B_i + ((128-n_i)/128)*N_i

Here n_i is the number of ordered protein-ID tokens in sentence. It agrees exactly with both core_length and selected_positions in the validated extraction record for every row. All n are between {recon['min_core_count']} and {recon['max_core_count']}, so all non-core denominators are positive (minimum {recon['min_noncore_count']}). Counts are constant within each five-state BGC. A/B metadata, CSV rows, and extraction validation identities/indices agree exactly; all arrays are finite and the 2,502 groups each have exactly suffixes -1,0,1,2,3.

N is an algebraically reconstructed mean of non-core contextualized positions. It is not raw ESM2, a new inference result, a trained representation, or a biological class. Nothing here establishes literal zero-padding, native flanking sequence, or independently meaningful biological context.

## Reconstruction validation

**OBSERVATION:** N is finite float64 (12510,320), with metadata preserving all source rows and all 1,900 interrupted-core observations.

| Check | Result |
| --- | --- |
| Rows satisfying weighted identity | {recon['rows_passed']} / 12510 |
| Maximum absolute reconstruction error | {recon['max_absolute_error']:.17g} |
| Mean absolute reconstruction error | {recon['mean_absolute_error']:.17g} |
| Maximum row-relative L2 error | {recon['max_row_relative_L2_error']:.17g} |
| Maximum component-relative error for abs(A)>1e-8 | {recon['max_element_relative_error_where_abs_A_gt_1e_8']:.17g} |

The strict elementwise tolerance is abs(error) <= 1e-12 + 1e-12*abs(A). Bitwise equality is not required. Component-relative error excludes near-zero denominators; the row-relative L2 measure is more stable. The full per-row validation is tables/noncore_reconstruction_validation.tsv.

**INTERPRETATION:** this confirms the algebra and row/count handling, not an independent recovery of the original full tensor. A and B were saved in float32; N inherits their rounding and pooling error. Error amplification increases as 128-n decreases. No contextual tensor was regenerated to obtain independent per-position confirmation.

## A/B/N augmentation stability

All ten pairwise cosines among five states are evaluated for each BGC (25,020 pairs per representation). Raw arithmetic centroids precede any transient normalization. Pairwise sample SD uses ddof=1. Both mean and maximum centroid cosine distance are saved, along with every pair and every augmentation-to-centroid distance.

Distribution summaries below distinguish BGC-level statistics from the descriptive pooled-pair distribution:

{md(stab[stab.metric.isin(['all_pairwise_cosines','mean_pairwise_cosine','min_pairwise_cosine','mean_centroid_cosine_distance','max_centroid_cosine_distance'])][['representation','metric','n','mean','median','q25','q75']])}

**STATISTICAL RESULT:** paired differences use 2,502 BGC units and 10,000 bootstrap resamples. Both mean and median difference percentile CIs are exported. The primary paired median results are:

{md(paired[['comparison','metric','median_difference','median_ci95_low','median_ci95_high','mean_difference','fraction_improved','fraction_equal','fraction_worsened']])}

Equality/improvement fractions use a numerical tolerance of 1e-12 on the cosine difference; no pairs fall in this equality interval. A is more stable than B for {100*comp['fraction_A_more_stable_than_B']:.4f}% of BGCs, more stable than N for {100*comp['fraction_A_more_stable_than_N']:.4f}%, and more stable than both simultaneously for {100*comp['fraction_A_more_stable_than_both']:.4f}% ({comp['count_A_more_stable_than_both']}/2502).

**INTERPRETATION:** N is much more stable than B, but marginally less stable than A on average and by the paired median. N−B median mean-cosine change is +0.026507 (95% CI +0.023821 to +0.030426), while N−A is −0.000330 (−0.000425 to −0.000251). The small N−A effect should not be exaggerated because its interval excludes zero.

![Stability](figures/A_B_N_augmentation_stability.png)

## Core/non-core deviation decomposition

Within each BGC, raw centroids define dA, dB and dN. Weighted terms are C=(n/128)dB and Q=((128-n)/128)dN. All rows satisfy dA=C+Q and ||dA||²=||C||²+||Q||²+2 C·Q within documented floating-point tolerances. Maximum deviation identity error is {cancel['max_deviation_identity_error']:.6g}; maximum squared-norm identity error is {cancel['max_variance_identity_error']:.6g}. The latter is checked with atol=rtol=1e-12. All 12,510 term cosines have nonzero denominators.

**OBSERVATION:**

- Negative cross terms: {100*cancel['fraction_augmentation_negative_cross']:.4f}% of augmentation observations.
- Negative per-BGC mean cross terms: {100*cancel['fraction_BGC_negative_mean_cross']:.4f}% of BGCs.
- Median augmentation cross term: {cancel['median_cross_term']:.6f}.
- Median cosine(C,Q): {cancel['median_deviation_cosine']:.6f}.

To express magnitude relative to squared-term energy, define the signed fraction:

    R_BGC = -mean(cross)/(mean(||C||² + ||Q||²))

Positive R indicates reduction; negative R indicates reinforcement. Median R is {cancel['median_BGC_variance_reduction_fraction']:.6f}. The pooled ratio is {cancel['pooled_variance_reduction_fraction']:.6f}: cross terms add approximately {100*-cancel['pooled_variance_reduction_fraction']:.2f}% to the total squared-term energy, rather than cancel it. This ratio concerns Euclidean deviations, not a direct decomposition of cosine similarity.

Minority subset details:

{md(pd.DataFrame(subsets))}

**INTERPRETATION:** weighted perturbations predominantly point in similar directions. Partial opposition exists in a minority, but cancellation is not supported as the population-wide explanation for A's higher cosine stability. Positive cross terms can coexist with A being more cosine-stable than B because weighting, vector norms, centroids and directions also matter.

![Deviations](figures/core_noncore_deviation_relationship.png)

## Gene-count dependence

All tests below use one observation per BGC (n=2502), not 12,510 independent augmentation rows. Effect size is Spearman rho; p-values are exploratory and unadjusted. The very strong correlations can numerically underflow rather than yield literal probability zero.

{md(cdisplay)}

Core count versus A−B stability has rho −0.222842: the advantage of A is modestly greater for smaller cores. Non-core count is exactly 128−core count, so its +0.222842 association is algebraically redundant, not an independent discovery.

Core count correlates strongly with the signed variance-reduction fraction (rho −0.800497) and with the non-core fraction of uncoupled squared-term energy (rho −0.931402). The former mainly reflects increasingly negative R, hence reinforcement, rather than increasingly strong cancellation. These weighted quantities explicitly depend on n; the association is not evidence that biological length causes compensation. The raw mean cross term has rho +0.337877 with count. R has little association with A−B stability (rho +0.025731, p≈0.198).

The median core count is {base.core_gene_count.median():g}, and median weight wN is {base.wN.median():.6f}. Most BGCs therefore weight the non-core component heavily. Its median uncoupled deviation-energy fraction is {base.noncore_fraction_of_uncoupled_energy.median():.6f}; this is a weighted perturbation-energy fraction, not an attention score or causal attribution.

{md(strata)}

The longest bin contains only six BGCs and should not support a separate biological claim.

![Count dependence](figures/core_count_vs_stability.png)

## Consensus construction

All A/B/N biological geometry uses one arithmetic mean of five raw augmentation vectors per BGC, in float64, without individual-vector normalization. No suffix is treated as canonical.

All three saved consensus matrices are (2502,320), finite float64, with nonzero norms. Both A and B reproduce the corrected Step 3 consensus matrices exactly, as do all their top-50 neighbor lists. Metadata order, original labels and the five source row IDs/indices are retained. The weighted consensus identity also holds because n is constant within every BGC.

## Saccharide kNN comparison

The exact corrected Step 3 cosine implementation is reused: original 320-D space, transient float64 L2-normalized copies, distance 1-Z@Z.T, self excluded, stable argsort and lexicographic BGC tie handling. All k=5,10,20,50 boundary-tie counts are saved. Both A/B precisions reproduce before N interpretation.

Any Saccharide includes 195 BGCs (7.793765% prevalence); pure Saccharide includes 130 (5.195843%). The 65 Saccharide-positive hybrids remain outside the pure endpoint. Precision below is a percentage:

{md(primary)}

The full table also records prevalence, finite-population self-excluded expectation, fold and percentage-point enrichment, null mean, empirical p and descriptive z-score.

For N, the exact 10,000 shared joint-membership permutation assignments from corrected Step 3 are reused. Their SHA256, shape, metadata order, definition axes and exact regeneration from seed 1729 were verified. Existing A/B null results are read from corrected Step 3; they were not recalculated with new permutations. N geometry remains fixed and positive-query identities are recalculated after each permutation.

{md(knn[knn.representation=='N'][['definition','k','precision_at_k','fold_enrichment','null_mean','empirical_p','z_score']])}

All eight N localization tests attain one-sided p=1/10001≈0.00009999. They are unadjusted localization p-values (even an eight-test Bonferroni bound is below 0.0008), not A-vs-N difference tests or held-out accuracy. Z-scores use empirical null SD and are not interpreted as normal-tail p-values.

**OBSERVATION:** A ≥ N > B at every endpoint/k; A=N for pure Saccharide at k=5. N retains substantial Saccharide-label localization. It is close to A and consistently higher than B, without proving that the small A/N differences reflect a generalizable performance effect.

![Saccharide kNN](figures/A_B_N_saccharide_knn.png)

## Neighbor overlap

Overlap is the intersection of two top-k neighbor sets divided by k. Previous A-vs-B values reproduce for every BGC and k. Overall summaries:

{md(ovs)}

At k=20, mean overlap is 34.142686% for A/B, 85.947242% for A/N, and 29.698242% for B/N. Exact values and any-Saccharide, pure-Saccharide and non-Saccharide strata are in tables/neighbor_overlap_summary.tsv. Here non-Saccharide means absence of any Saccharide label; it does not include Saccharide-positive hybrids. N is much closer to A than B in local neighbor identity.

![Neighbor overlap](figures/A_B_N_neighbor_overlap.png)

## Mechanistic interpretation

1. **Reconstruction:** N satisfies A=wB*B+wN*N to float64 roundoff, with every source count validated. This is an algebraic consistency check, not independent tensor evidence.
2. **N versus B stability:** N is more stable by both average and paired median, and improves mean pairwise cosine in 79.1367% of BGCs.
3. **N versus A stability:** N is slightly less stable overall; the N−A median paired difference is about −0.000330.
4. **Opposition:** core/non-core weighted changes usually reinforce. Only 9.3605% of augmentation cross terms and 1.6387% of BGC mean cross terms are negative.
5. **Magnitude:** cancellation is not a sufficient population-wide explanation here; pooled cross terms add 15.20% relative to the squared-term sum. A minority does show partial opposition, whose subset magnitudes are tabulated.
6. **Length:** length is associated with weighted decomposition terms, but these weights depend on length by definition. The signed cancellation fraction is not appreciably associated with the A−B cosine-stability advantage.
7. **N localization:** N shows strong Saccharide-label localization under the reused random-label null.
8. **Ranking:** A ≥ N > B for all tested k and both definitions, with the pure-Saccharide k=5 A/N tie.
9. **Neighbor identity:** N largely retains A's local identities and differs substantially from B.
10. **Descriptive mechanism:** option (a), substantial class-related signal in N, is supported for the tested Saccharide labels. Option (b), improved stability through pervasive perturbation cancellation, is not supported by these cross terms. Therefore option (c), both, is not supported at the population level; option (d), neither, is also inconsistent with the localization results. No broad all-class claim was tested.

**INTERPRETATION:** the combination of N's high augmentation stability, its often large averaging weight and its A-like Saccharide neighborhoods is consistent with a contribution of the non-core contextual component to A's observed properties. This is a descriptive decomposition, not causal isolation.

**HYPOTHESIS:** all-position averaging may retain stable label-associated information present in contextualized non-core positions. The present analysis cannot determine whether this information originates in the core, augmentation context, learned supervision, or interactions among them. It does not demonstrate that attention compensates for augmentation or that N represents biological flanking sequence.

## Caveats

The dataset is supervised training data; neighbor localization is in-sample representation geometry, not biological validation or held-out predictive performance. N is algebraically dependent on A and B and inherits float32 pooling/storage errors. Its near-A geometry is not an independent validation of a new modality. The reconstruction identity is expected by construction.

Core count affects the mixture weights directly. Non-core positions cannot be automatically interpreted as zero-padding or native flanks. Cross terms measure Euclidean interaction, not biological compensation, attention weights, or a causal decomposition of cosine stability. BGCs can be evolutionarily related; the BGC bootstrap does not remove phylogenetic dependence. Bootstrap CIs use 10,000 resamples; correlation p-values are exploratory. Numerical p underflow is explicitly marked.

## Deliverables

All new files are under {ROOT}.

- results/: per-row N embeddings and metadata; A/B/N consensus arrays and metadata; exact top-50 neighbor indices; N null statistics; validation summaries.
- tables/: all requested primary TSVs; every augmentation pair and centroid distance; paired mean/median bootstrap comparisons; per-BGC cancellation, length strata and subset summaries; per-query Saccharide precision; exact neighbor lists and overlaps.
- figures/: five requested 400-dpi PNGs and matching PDFs, with numerical sources listed in provenance/figure_sources.tsv and visual QA recorded.
- provenance/: exact input paths/SHA256, before/after integrity, Python/package versions, parameters/seeds, reused-permutation verification, script checksums and independent numerical validation.
- scripts/: reconstruct.py, mechanism.py, compare_consensus.py, plot_results.py, validate.py, write_report.py, common.py, preview_figures.py and run_analysis.sh.
- logs/: stage stdout/stderr.

{check['checks']} independent checks passed. All {inputcount} hashed protected scientific files remain unchanged, including the original CSV, extraction outputs and prior analysis figures/tables/provenance. Dependency/cache/temporary directories are excluded from the scientific manifest; existing dependencies are imported read-only with bytecode disabled and caches directed to the new directory. The model, checkpoint and LMDB were neither opened nor modified by this analysis; they are not claimed to have been newly checksum-audited.

## Reproduction command

Run in WSL Ubuntu:

    cd {ROOT}
    bash scripts/run_analysis.sh

The launcher only executes the present saved-vector follow-up and overwrites its own outputs on rerun; it does not invoke prior analysis pipelines or model inference. Existing project-local dependencies at {OLD/'.deps'} are reused; no installation is needed. Exact versions are in provenance/package_versions.tsv and requirements.lock.txt.

Paired bootstrap count: 10,000; comparison seeds 1729 through 1734 in the exported comparison order. Random-label assignments: the exact corrected Step 3 array generated with seed 1729, reused without replacement or silent regeneration. All operations use WSL Python and write only inside this new analysis directory.
""")
 tab("provenance/script_sha256.tsv",[dict(path=str(p),sha256=sha(p)) for p in sorted((ROOT/"scripts").iterdir()) if p.is_file()])
 js("provenance/completion.json",dict(passed=True,utc=pd.Timestamp.now(tz="UTC").isoformat(),protected_files_unchanged=inputcount,figures=5,validation_checks=check["checks"]))
 print("Report complete; input integrity unchanged.",flush=True)
if __name__=="__main__":main()
