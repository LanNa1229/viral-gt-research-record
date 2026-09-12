# Non-core contextual component mechanistic analysis

## Inputs and algebraic definition

This short follow-up uses saved vectors only. It does not rerun model inference, extraction, the full prior Step 2, corrected consensus Step 3, PCA, UMAP or clustering. The A/B statistics needed for this explicitly requested three-component comparison were recomputed locally and their consensus geometry verified against corrected Step 3.

| role | path | sha256 |
| --- | --- | --- |
| A | /home/rahnn/260908_MiBG/analysis/bgcprophet_reference/data/mibig_12510_embeddings.npy | 9a747f4a73e4d76074c95e6d2e8fbeb329e688e7ed27bda11d8a5ed8efb7d038 |
| A metadata | /home/rahnn/260908_MiBG/analysis/bgcprophet_reference/data/mibig_12510_metadata.tsv | faaee7f876a610b6b0c55e7dd982c8556d32e2cd4ff179c03855234c5e7706f9 |
| B | /home/rahnn/260908_MiBG/outputs/core_context_full_20260910/results_validatednode/mibig_12510_core_context_embeddings.npy | 96982f0eb361b628ba0f6845234bc84ac592d7b9784f650b929efc76b4ffbc1d |
| B metadata | /home/rahnn/260908_MiBG/outputs/core_context_full_20260910/results_validatednode/mibig_12510_core_context_metadata.tsv | faaee7f876a610b6b0c55e7dd982c8556d32e2cd4ff179c03855234c5e7706f9 |
| CSV | /home/rahnn/260908_MiBG/BGC_train_dataset_classify.csv | 78dab471e8ef7dcdec16f6b4eb82ef15ef01b8bf900743c14238d26186a1d7f9 |


A is the all-128-position contextualized mean. B is the true-core contextualized mean. The non-core contextual component is reconstructed in float64:

    N_i = (128*A_i - n_i*B_i)/(128-n_i)
    A_i = (n_i/128)*B_i + ((128-n_i)/128)*N_i

Here n_i is the number of ordered protein-ID tokens in sentence. It agrees exactly with both core_length and selected_positions in the validated extraction record for every row. All n are between 1 and 115, so all non-core denominators are positive (minimum 13). Counts are constant within each five-state BGC. A/B metadata, CSV rows, and extraction validation identities/indices agree exactly; all arrays are finite and the 2,502 groups each have exactly suffixes -1,0,1,2,3.

N is an algebraically reconstructed mean of non-core contextualized positions. It is not raw ESM2, a new inference result, a trained representation, or a biological class. Nothing here establishes literal zero-padding, native flanking sequence, or independently meaningful biological context.

## Reconstruction validation

**OBSERVATION:** N is finite float64 (12510,320), with metadata preserving all source rows and all 1,900 interrupted-core observations.

| Check | Result |
| --- | --- |
| Rows satisfying weighted identity | 12510 / 12510 |
| Maximum absolute reconstruction error | 2.2204460492503131e-16 |
| Mean absolute reconstruction error | 3.6392125317071396e-19 |
| Maximum row-relative L2 error | 5.1787536764868139e-17 |
| Maximum component-relative error for abs(A)>1e-8 | 8.25640580332871e-13 |

The strict elementwise tolerance is abs(error) <= 1e-12 + 1e-12*abs(A). Bitwise equality is not required. Component-relative error excludes near-zero denominators; the row-relative L2 measure is more stable. The full per-row validation is tables/noncore_reconstruction_validation.tsv.

**INTERPRETATION:** this confirms the algebra and row/count handling, not an independent recovery of the original full tensor. A and B were saved in float32; N inherits their rounding and pooling error. Error amplification increases as 128-n decreases. No contextual tensor was regenerated to obtain independent per-position confirmation.

## A/B/N augmentation stability

All ten pairwise cosines among five states are evaluated for each BGC (25,020 pairs per representation). Raw arithmetic centroids precede any transient normalization. Pairwise sample SD uses ddof=1. Both mean and maximum centroid cosine distance are saved, along with every pair and every augmentation-to-centroid distance.

Distribution summaries below distinguish BGC-level statistics from the descriptive pooled-pair distribution:

| representation | metric | n | mean | median | q25 | q75 |
| --- | --- | --- | --- | --- | --- | --- |
| A | mean_pairwise_cosine | 2502 | 0.932363 | 0.96268 | 0.91325 | 0.982026 |
| A | min_pairwise_cosine | 2502 | 0.860863 | 0.920402 | 0.818988 | 0.963216 |
| A | mean_centroid_cosine_distance | 2502 | 0.0282588 | 0.0151219 | 0.00722623 | 0.0356102 |
| A | max_centroid_cosine_distance | 2502 | 0.0615008 | 0.0317995 | 0.0136858 | 0.0769411 |
| A | all_pairwise_cosines | 25020 | 0.932363 | 0.969034 | 0.920936 | 0.98749 |
| B | mean_pairwise_cosine | 2502 | 0.867491 | 0.924871 | 0.807946 | 0.969423 |
| B | min_pairwise_cosine | 2502 | 0.734747 | 0.837668 | 0.603514 | 0.937315 |
| B | mean_centroid_cosine_distance | 2502 | 0.0574754 | 0.0307486 | 0.0123296 | 0.0808626 |
| B | max_centroid_cosine_distance | 2502 | 0.122818 | 0.0666319 | 0.0244152 | 0.174975 |
| B | all_pairwise_cosines | 25020 | 0.867491 | 0.940758 | 0.835599 | 0.978985 |
| N | mean_pairwise_cosine | 2502 | 0.931796 | 0.961598 | 0.912447 | 0.981258 |
| N | min_pairwise_cosine | 2502 | 0.859035 | 0.918006 | 0.817153 | 0.961668 |
| N | mean_centroid_cosine_distance | 2502 | 0.0285044 | 0.0155516 | 0.00752645 | 0.0361621 |
| N | max_centroid_cosine_distance | 2502 | 0.0624281 | 0.0325621 | 0.0143413 | 0.0777392 |
| N | all_pairwise_cosines | 25020 | 0.931796 | 0.968752 | 0.920374 | 0.987245 |


**STATISTICAL RESULT:** paired differences use 2,502 BGC units and 10,000 bootstrap resamples. Both mean and median difference percentile CIs are exported. The primary paired median results are:

| comparison | metric | median_difference | median_ci95_low | median_ci95_high | mean_difference | fraction_improved | fraction_equal | fraction_worsened |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B - A | mean_pairwise_cosine | -0.0273862 | -0.03097 | -0.0251267 | -0.0648718 | 0.166267 | 0 | 0.833733 |
| B - A | min_pairwise_cosine | -0.0574793 | -0.0636155 | -0.0508427 | -0.126116 | 0.188249 | 0 | 0.811751 |
| N - A | mean_pairwise_cosine | -0.000330257 | -0.000424971 | -0.000250733 | -0.000566797 | 0.417266 | 0 | 0.582734 |
| N - A | min_pairwise_cosine | -0.000778873 | -0.000984021 | -0.000582222 | -0.00182861 | 0.41247 | 0 | 0.58753 |
| N - B | mean_pairwise_cosine | 0.0265065 | 0.0238205 | 0.0304264 | 0.064305 | 0.791367 | 0 | 0.208633 |
| N - B | min_pairwise_cosine | 0.0554639 | 0.0493816 | 0.0621704 | 0.124287 | 0.770184 | 0 | 0.229816 |


Equality/improvement fractions use a numerical tolerance of 1e-12 on the cosine difference; no pairs fall in this equality interval. A is more stable than B for 83.3733% of BGCs, more stable than N for 58.2734%, and more stable than both simultaneously for 41.6467% (1042/2502).

**INTERPRETATION:** N is much more stable than B, but marginally less stable than A on average and by the paired median. N−B median mean-cosine change is +0.026507 (95% CI +0.023821 to +0.030426), while N−A is −0.000330 (−0.000425 to −0.000251). The small N−A effect should not be exaggerated because its interval excludes zero.

![Stability](figures/A_B_N_augmentation_stability.png)

## Core/non-core deviation decomposition

Within each BGC, raw centroids define dA, dB and dN. Weighted terms are C=(n/128)dB and Q=((128-n)/128)dN. All rows satisfy dA=C+Q and ||dA||²=||C||²+||Q||²+2 C·Q within documented floating-point tolerances. Maximum deviation identity error is 9.71445e-16; maximum squared-norm identity error is 1.42109e-14. The latter is checked with atol=rtol=1e-12. All 12,510 term cosines have nonzero denominators.

**OBSERVATION:**

- Negative cross terms: 9.3605% of augmentation observations.
- Negative per-BGC mean cross terms: 1.6387% of BGCs.
- Median augmentation cross term: 0.207681.
- Median cosine(C,Q): 0.662973.

To express magnitude relative to squared-term energy, define the signed fraction:

    R_BGC = -mean(cross)/(mean(||C||² + ||Q||²))

Positive R indicates reduction; negative R indicates reinforcement. Median R is -0.151636. The pooled ratio is -0.152015: cross terms add approximately 15.20% to the total squared-term energy, rather than cancel it. This ratio concerns Euclidean deviations, not a direct decomposition of cosine similarity.

Minority subset details:

| subset | n_BGC | median_signed_variance_reduction | median_A_minus_B_stability | fraction_A_more_stable_than_B |
| --- | --- | --- | --- | --- |
| negative_mean_cross | 41 | 0.0209179 | 0.0720008 | 0.902439 |
| nonnegative_mean_cross | 2461 | -0.155493 | 0.026784 | 0.832588 |


**INTERPRETATION:** weighted perturbations predominantly point in similar directions. Partial opposition exists in a minority, but cancellation is not supported as the population-wide explanation for A's higher cosine stability. Positive cross terms can coexist with A being more cosine-stable than B because weighting, vector norms, centroids and directions also matter.

![Deviations](figures/core_noncore_deviation_relationship.png)

## Gene-count dependence

All tests below use one observation per BGC (n=2502), not 12,510 independent augmentation rows. Effect size is Spearman rho; p-values are exploratory and unadjusted. The very strong correlations can numerically underflow rather than yield literal probability zero.

| x | y | n_BGC | spearman_rho | p_value | absolute_rho |
| --- | --- | --- | --- | --- | --- |
| core_gene_count | A_mean_pairwise_cosine | 2502 | 0.23395 | 1.8804e-32 | 0.23395 |
| core_gene_count | B_mean_pairwise_cosine | 2502 | 0.292053 | 2.17167e-50 | 0.292053 |
| core_gene_count | N_mean_pairwise_cosine | 2502 | 0.225213 | 3.88727e-30 | 0.225213 |
| core_gene_count | A_minus_B_stability | 2502 | -0.222842 | 1.58955e-29 | 0.222842 |
| noncore_position_count | A_minus_B_stability | 2502 | 0.222842 | 1.58955e-29 | 0.222842 |
| core_gene_count | mean_cross_term | 2502 | 0.337877 | 7.35397e-68 | 0.337877 |
| core_gene_count | variance_reduction_fraction_from_cross | 2502 | -0.800497 | numerical underflow | 0.800497 |
| core_gene_count | noncore_fraction_of_uncoupled_energy | 2502 | -0.931402 | numerical underflow | 0.931402 |
| mean_cross_term | A_minus_B_stability | 2502 | 0.245622 | 1.06921e-35 | 0.245622 |
| variance_reduction_fraction_from_cross | A_minus_B_stability | 2502 | 0.0257306 | 0.198229 | 0.0257306 |


Core count versus A−B stability has rho −0.222842: the advantage of A is modestly greater for smaller cores. Non-core count is exactly 128−core count, so its +0.222842 association is algebraically redundant, not an independent discovery.

Core count correlates strongly with the signed variance-reduction fraction (rho −0.800497) and with the non-core fraction of uncoupled squared-term energy (rho −0.931402). The former mainly reflects increasingly negative R, hence reinforcement, rather than increasingly strong cancellation. These weighted quantities explicitly depend on n; the association is not evidence that biological length causes compensation. The raw mean cross term has rho +0.337877 with count. R has little association with A−B stability (rho +0.025731, p≈0.198).

The median core count is 13, and median weight wN is 0.898438. Most BGCs therefore weight the non-core component heavily. Its median uncoupled deviation-energy fraction is 0.980164; this is a weighted perturbation-energy fraction, not an attention score or causal attribution.

| core_gene_count_bin | n_BGC | median_core_count | median_A_mean_pairwise_cosine | median_B_mean_pairwise_cosine | median_N_mean_pairwise_cosine | median_mean_cross_term | median_variance_reduction_fraction_from_cross | median_noncore_fraction_of_uncoupled_energy | median_A_minus_B_stability | fraction_BGC_negative_mean_cross |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1-5 | 443 | 4 | 0.933437 | 0.8235 | 0.933528 | 0.125683 | -0.0431113 | 0.998524 | 0.0784188 | 0.0248307 |
| 6-10 | 611 | 8 | 0.958133 | 0.907043 | 0.956521 | 0.214093 | -0.102899 | 0.992855 | 0.0322814 | 0.00981997 |
| 11-20 | 732 | 14 | 0.966831 | 0.935039 | 0.965569 | 0.321755 | -0.191855 | 0.973686 | 0.0216879 | 0.0191257 |
| 21-40 | 573 | 27 | 0.968075 | 0.942079 | 0.967941 | 0.412311 | -0.34774 | 0.886702 | 0.0206213 | 0.0157068 |
| 41-80 | 137 | 48 | 0.974905 | 0.958067 | 0.97015 | 0.361144 | -0.45984 | 0.598542 | 0.0159079 | 0 |
| 81-127 | 6 | 105 | 0.99121 | 0.990053 | 0.982733 | 0.0591088 | -0.21152 | 0.101464 | 0.000626846 | 0.166667 |


The longest bin contains only six BGCs and should not support a separate biological claim.

![Count dependence](figures/core_count_vs_stability.png)

## Consensus construction

All A/B/N biological geometry uses one arithmetic mean of five raw augmentation vectors per BGC, in float64, without individual-vector normalization. No suffix is treated as canonical.

All three saved consensus matrices are (2502,320), finite float64, with nonzero norms. Both A and B reproduce the corrected Step 3 consensus matrices exactly, as do all their top-50 neighbor lists. Metadata order, original labels and the five source row IDs/indices are retained. The weighted consensus identity also holds because n is constant within every BGC.

## Saccharide kNN comparison

The exact corrected Step 3 cosine implementation is reused: original 320-D space, transient float64 L2-normalized copies, distance 1-Z@Z.T, self excluded, stable argsort and lexicographic BGC tie handling. All k=5,10,20,50 boundary-tie counts are saved. Both A/B precisions reproduce before N interpretation.

Any Saccharide includes 195 BGCs (7.793765% prevalence); pure Saccharide includes 130 (5.195843%). The 65 Saccharide-positive hybrids remain outside the pure endpoint. Precision below is a percentage:

| definition | k | A | B | N |
| --- | --- | --- | --- | --- |
| any_saccharide | 5 | 76.9231 | 70.6667 | 76.2051 |
| any_saccharide | 10 | 76.6154 | 69.0256 | 75.2821 |
| any_saccharide | 20 | 74.9487 | 67.6154 | 74.0256 |
| any_saccharide | 50 | 72.2462 | 63.0872 | 71.4359 |
| pure_saccharide | 5 | 92.9231 | 82.3077 | 92.9231 |
| pure_saccharide | 10 | 93 | 81.1538 | 92.3077 |
| pure_saccharide | 20 | 92.5769 | 79.2692 | 92.4231 |
| pure_saccharide | 50 | 92.5077 | 74.4769 | 92.2154 |


The full table also records prevalence, finite-population self-excluded expectation, fold and percentage-point enrichment, null mean, empirical p and descriptive z-score.

For N, the exact 10,000 shared joint-membership permutation assignments from corrected Step 3 are reused. Their SHA256, shape, metadata order, definition axes and exact regeneration from seed 1729 were verified. Existing A/B null results are read from corrected Step 3; they were not recalculated with new permutations. N geometry remains fixed and positive-query identities are recalculated after each permutation.

| definition | k | precision_at_k | fold_enrichment | null_mean | empirical_p | z_score |
| --- | --- | --- | --- | --- | --- | --- |
| any_saccharide | 5 | 0.762051 | 9.7777 | 0.0775373 | 9.999e-05 | 63.8465 |
| any_saccharide | 10 | 0.752821 | 9.65927 | 0.0775473 | 9.999e-05 | 86.4416 |
| any_saccharide | 20 | 0.740256 | 9.49806 | 0.0775741 | 9.999e-05 | 112.878 |
| any_saccharide | 50 | 0.714359 | 9.16578 | 0.077549 | 9.999e-05 | 147.682 |
| pure_saccharide | 5 | 0.929231 | 17.8841 | 0.051598 | 9.999e-05 | 80.5957 |
| pure_saccharide | 10 | 0.923077 | 17.7657 | 0.0516318 | 9.999e-05 | 111.008 |
| pure_saccharide | 20 | 0.924231 | 17.7879 | 0.0515705 | 9.999e-05 | 151.653 |
| pure_saccharide | 50 | 0.922154 | 17.7479 | 0.0515735 | 9.999e-05 | 215.236 |


All eight N localization tests attain one-sided p=1/10001≈0.00009999. They are unadjusted localization p-values (even an eight-test Bonferroni bound is below 0.0008), not A-vs-N difference tests or held-out accuracy. Z-scores use empirical null SD and are not interpreted as normal-tail p-values.

**OBSERVATION:** A ≥ N > B at every endpoint/k; A=N for pure Saccharide at k=5. N retains substantial Saccharide-label localization. It is close to A and consistently higher than B, without proving that the small A/N differences reflect a generalizable performance effect.

![Saccharide kNN](figures/A_B_N_saccharide_knn.png)

## Neighbor overlap

Overlap is the intersection of two top-k neighbor sets divided by k. Previous A-vs-B values reproduce for every BGC and k. Overall summaries:

| comparison | k | mean | median |
| --- | --- | --- | --- |
| A vs B | 5 | 0.216787 | 0.2 |
| A vs B | 10 | 0.271343 | 0.2 |
| A vs B | 20 | 0.341427 | 0.35 |
| A vs B | 50 | 0.458217 | 0.46 |
| A vs N | 5 | 0.794404 | 0.8 |
| A vs N | 10 | 0.829416 | 0.9 |
| A vs N | 20 | 0.859472 | 0.9 |
| A vs N | 50 | 0.898122 | 0.92 |
| B vs N | 5 | 0.178657 | 0.2 |
| B vs N | 10 | 0.229337 | 0.2 |
| B vs N | 20 | 0.296982 | 0.3 |
| B vs N | 50 | 0.412902 | 0.42 |


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

All new files are under /home/rahnn/260908_MiBG/analysis/noncore_context_mechanism_20260911.

- results/: per-row N embeddings and metadata; A/B/N consensus arrays and metadata; exact top-50 neighbor indices; N null statistics; validation summaries.
- tables/: all requested primary TSVs; every augmentation pair and centroid distance; paired mean/median bootstrap comparisons; per-BGC cancellation, length strata and subset summaries; per-query Saccharide precision; exact neighbor lists and overlaps.
- figures/: five requested 400-dpi PNGs and matching PDFs, with numerical sources listed in provenance/figure_sources.tsv and visual QA recorded.
- provenance/: exact input paths/SHA256, before/after integrity, Python/package versions, parameters/seeds, reused-permutation verification, script checksums and independent numerical validation.
- scripts/: reconstruct.py, mechanism.py, compare_consensus.py, plot_results.py, validate.py, write_report.py, common.py, preview_figures.py and run_analysis.sh.
- logs/: stage stdout/stderr.

481 independent checks passed. All 278 hashed protected scientific files remain unchanged, including the original CSV, extraction outputs and prior analysis figures/tables/provenance. Dependency/cache/temporary directories are excluded from the scientific manifest; existing dependencies are imported read-only with bytecode disabled and caches directed to the new directory. The model, checkpoint and LMDB were neither opened nor modified by this analysis; they are not claimed to have been newly checksum-audited.

## Reproduction command

Run in WSL Ubuntu:

    cd /home/rahnn/260908_MiBG/analysis/noncore_context_mechanism_20260911
    bash scripts/run_analysis.sh

The launcher only executes the present saved-vector follow-up and overwrites its own outputs on rerun; it does not invoke prior analysis pipelines or model inference. Existing project-local dependencies at /home/rahnn/260908_MiBG/analysis/core_context_representation_geometry_20260910/.deps are reused; no installation is needed. Exact versions are in provenance/package_versions.tsv and requirements.lock.txt.

Paired bootstrap count: 10,000; comparison seeds 1729 through 1734 in the exported comparison order. Random-label assignments: the exact corrected Step 3 array generated with seed 1729, reused without replacement or silent regeneration. All operations use WSL Python and write only inside this new analysis directory.
