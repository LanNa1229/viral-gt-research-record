# BGC-Prophet MIBiG consensus embedding analysis

This pipeline analyzes the final 320-dimensional Transformer representations extracted on the validated HPC environment. It performs no Transformer inference and does not modify the source files.

## Reproduce

Run with Python 3.12 and the packages in requirements.lock.txt. The completed run used the Codex bundled Windows Python and analysis-only dependencies installed into this directory's .deps; no BGC-Prophet package or inference environment was installed or changed.

From this directory, with that Python interpreter:

    python -B scripts/analyze_reference.py --stage qc
    python -B scripts/analyze_reference.py --stage downstream

Or reproduce both stages in order:

    python -B scripts/analyze_reference.py --stage all

The script detects this directory relative to its own location; the working directory is immaterial. Optional --data-dir and --output-dir parameters select inputs and results. Dependencies can be reproduced with:

    python -m pip install --target .deps --no-cache-dir -r requirements.lock.txt

The lock records the exact installed versions. Platform-specific wheels will differ across operating systems. Exact numerical reproduction is expected within the recorded runtime, versions, and thread settings; cross-platform bitwise identity is not promised. Existing analysis results may be overwritten on rerun; data files are opened read-only.

## Validation and order of operations

The pipeline requires a finite, nonzero (12510,320) embedding array; unique metadata indices covering 0..12509; exact ID, label, and isBGC agreement with the CSV in embedding-row order; 2502 base accessions each with exactly suffixes -1,0,1,2,3; and constant labels and original sentence within each BGC. The analysis CSV must be byte-identical to the root CSV. The (12510,7) probability array is checked for finite values in [0,1], without assuming its column-to-class mapping or requiring rows to sum to one. Probabilities are not used in this analysis.

Alignment checks validate the supplied metadata contract. Without independent HPC extraction logs, array contents alone cannot prove that the producer assigned the correct biological sample to each row.

All 10 unordered pairwise cosine similarities are computed from the five raw embeddings for each BGC. A minimum cosine below 0.95 is an operational instability flag selected before looking at stability values, not a validated biological threshold. A separate relative flag marks mean cosines below Q1 - 1.5*IQR across BGCs. Complete rankings and distributions are retained. The 25,020 pairs are dependent observations nested within 2,502 BGCs.

QC files are written before any consensus is constructed. All BGCs passing structural checks, including instability flags, remain in the downstream analysis. Consensus is the arithmetic mean of the five raw 320-D vectors, accumulated in float64. Individual embeddings are not normalized before averaging.

## Representations and neighborhoods

PCA centers the raw consensus vectors, uses full SVD, and performs neither feature standardization nor whitening. It is a QC representation.

UMAP receives the 2502 consensus vectors: cosine metric, 15 neighbors, min_dist 0.1, 2 components, spectral initialization, 500 epochs, random_state and transform_seed 1729, one thread. UMAP is only a visualization; its apparent separation is not quantitative evidence for classes.

kNN uses cosine distance directly in the original 320-D consensus space, excluding the query itself. All 2502 BGCs are neighbor candidates. Ties use lexicographic accession order deterministically. Query sets and positive neighbor definitions are:
- Any Saccharide: the source label contains Saccharide, including hybrid labels.
- Pure Saccharide: the entire trimmed source label equals Saccharide. Hybrid BGCs remain in the candidate pool but are negative for this definition.

Precision@k is the mean positive-neighbor fraction among positive query BGCs, for k=5,10,20. It is compared with global prevalence M/N; fold enrichment and percentage-point differences are provided. The additional finite-population, self-excluded random baseline is (M-1)/(N-1). These are descriptive results, not significance tests.

HDBSCAN uses sklearn.cluster.HDBSCAN on the precomputed cosine distance matrix from original 320-D consensus vectors, not UMAP. Primary parameters are min_cluster_size=30, min_samples=10, EOM selection, allow_single_cluster=False, algorithm=brute, n_jobs=1, copy=True. A checksum assertion confirms that every fit leaves the original cosine matrix unchanged. Sensitivity covers min_cluster_size=15,30,60 crossed with min_samples=5,10,20. sklearn's min_samples includes the point itself. Cluster -1 denotes noise; membership strength is not a probability of a biological class. Clusters and label compositions are exploratory.

These are training-set representations from a supervised model. Label enrichment is in-sample representation structure, not held-out predictive performance or independent confirmation of biological classes. Consensus averaging reduces augmentation multiplicity but does not establish padding invariance.

## Outputs

All numerical results are TSV, with additional NPY/JSON files for convenient reuse and provenance:
- validation.tsv, input_checksums.tsv, input_integrity_after.tsv
- augmentation_pairwise_cosines.tsv, augmentation_stability.tsv, augmentation_stability_ranked.tsv, unstable_bgcs.tsv, stability_global_distribution.tsv
- consensus_embeddings.tsv and .npy, consensus_metadata.tsv
- pca_scores.tsv, pca_components.tsv, pca_feature_means.tsv, pca_explained_variance.tsv
- umap_coordinates.tsv
- cosine_neighbors_top20.tsv, saccharide_knn_per_bgc.tsv, saccharide_knn_summary.tsv
- hdbscan_primary_assignments.tsv, hdbscan_primary_composition.tsv, hdbscan_all_assignments.tsv, hdbscan_sensitivity.tsv
- package_versions.tsv, provenance_*.json, dependency_install.json, run_*.log

Figures in results/figures are vector PDF (scatter points rasterized for file size) and 600-dpi PNG. Numerical source tables accompany each plot. PDFs embed TrueType text.

Input SHA256 values are recorded before and after execution, with a hard failure if they change. Provenance records script version/hash, arguments, UTC completion time, package versions, platform, Python, seed, and analysis settings. Thread counts are fixed to one.

Methods references:
- https://umap-learn.readthedocs.io/en/latest/reproducibility.html
- https://umap-learn.readthedocs.io/en/latest/parameters.html
- https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html

## Results

All 13 structural checks passed: 12,510 finite 320-D embeddings, exact metadata/CSV alignment, and 2,502 BGCs with five augmentations each. Original input SHA256 values remained unchanged. Independent numerical verification passed for all within-BGC cosines, all cosine neighbors, consensus arithmetic, PCA reconstruction, and saved coordinate alignment.

Among 25,020 pairwise cosines, the median was 0.969034, the mean 0.932363, and the minimum -0.085034. The minimum-pair cutoff of 0.95 flagged 1,655/2,502 BGCs (66.1%). The relative mean-cosine cutoff was 0.810086, flagging 210 BGCs. Padding sensitivity is substantial despite a high overall median.

The five lowest minimum pairwise cosines were:

| BGC | Label | Minimum cosine |
|---|---|---:|
| BGC0001534 | Other | -0.085034 |
| BGC0000871 | Other | -0.061353 |
| BGC0000673 | Terpene | -0.058810 |
| BGC0002249 | Polyketide | 0.009061 |
| BGC0002488 | Other | 0.035763 |

All 2,502 BGCs were retained in the (2502,320) consensus. PC1 explained 39.1% and PC2 27.5% of variance; 5 components reached 90% cumulative variance.

| Saccharide definition | Positive / total | Global prevalence | k | Precision@k | Fold / prevalence |
|---|---:|---:|---:|---:|---:|
| any_saccharide | 195/2502 | 7.79% | 5 | 76.92% | 9.87 |
| any_saccharide | 195/2502 | 7.79% | 10 | 76.62% | 9.83 |
| any_saccharide | 195/2502 | 7.79% | 20 | 74.95% | 9.62 |
| pure_saccharide | 130/2502 | 5.20% | 5 | 92.92% | 17.88 |
| pure_saccharide | 130/2502 | 5.20% | 10 | 93.00% | 17.90 |
| pure_saccharide | 130/2502 | 5.20% | 20 | 92.58% | 17.82 |

Primary HDBSCAN found 6 clusters and classified 843 BGCs (33.7%) as noise. The nine-setting sensitivity table accompanies the exploratory assignments. These results do not establish biological classes.

| Padding QC cohort | BGCs | Flagged | Fraction |
|---|---:|---:|---:|
| all | 2502 | 1655 | 66.1% |
| any_saccharide | 195 | 171 | 87.7% |
| pure_saccharide | 130 | 127 | 97.7% |

This is a descriptive reference analysis. No p-values or held-out performance estimates are claimed. High Saccharide neighborhood enrichment must be considered together with supervised training context and observed augmentation sensitivity.

Run the following after the main pipeline to reproduce independent checks, cohort QC summary, and this results section:

    python -B scripts/verify_and_summarize.py

The full execution arguments and versions are in results/provenance_downstream.json. Six figures were exported as PDF and PNG.
