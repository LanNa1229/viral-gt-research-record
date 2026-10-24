# Original-320D clustering of selected viral GT-associated loci

Read [REPORT.md](REPORT.md) for results, [METHODS_RECORD.md](METHODS_RECORD.md) for executed procedures and a Nature Microbiology Methods draft, and [EXECUTION_EVIDENCE_MATRIX.tsv](EXECUTION_EVIDENCE_MATRIX.tsv) for parameter-to-execution links. [EXECUTION_LOG.md](EXECUTION_LOG.md) indexes successful and failed attempts.

The prespecified baseline identifies five Leiden communities (VC01=68, VC02=55, VC03=74, VC04=35, VC05=55), but the five-way partition is exploratory. HDBSCAN assigns all 287 loci to noise at the fixed 15/5 baseline; cohort and parameter sensitivities materially change some boundaries. Community labels do not imply biological families or functions.

## Key artifacts

- [Figure C1: Leiden on saved UMAP](figures/Figure_C1_saved_UMAP.png) and [SVG](figures/Figure_C1_saved_UMAP.svg).
- [Figure C2: HDBSCAN/noise](figures/Figure_C2_saved_UMAP.png) and [SVG](figures/Figure_C2_saved_UMAP.svg).
- [Figure C3: separate joint-community relationships](figures/Figure_C3_community_relationships.png) and [SVG](figures/Figure_C3_community_relationships.svg).
- [Figure C4: sizes, stability and confounder controls](figures/Figure_C4_size_stability_QC.png) and [SVG](figures/Figure_C4_size_stability_QC.svg).
- [Per-locus primary assignments](tables/per_locus_selected287_clustering.csv), [all-556 records](tables/per_locus_all556_clustering.csv), and [community summary](tables/community_summary.csv).
- [Baseline cross-cohort ARI](tables/cross_cohort_partition_agreement.csv), [matched-grid comparisons](stability/cross_cohort_matched_grid_summary.csv), and [community sensitivity Jaccard](tables/community_sensitivity_Jaccard.csv).
- `assignments/*_leiden_grid.csv`: every seed-specific label for all 81 fits in each of six cohorts. Joint GC/SC IDs are distinct from viral-only VC IDs.
- `stability/`: full grid metrics, within-setting seed ARI, all-grid primary ARI, coassignment frequencies, HDBSCAN grid results and matched-setting cross-cohort ARI/Jaccard.
- `graphs/`: baseline weighted edges, all-k graph diagnostics, node hubness/component statistics and primary directed neighbor ranks.
- `tables/`: reference nearest neighbors/class composition, joint reference/viral ratios, genome diversity, duplicate registry, full15 exclusions, extreme short-window QC and gene-count/hubness strata.
- `runs/`: pre-execution records, exact source/configuration snapshots, distribution/native file fingerprints, stdout/stderr, actual exits, after-execution integrity and output manifests.
- `logs/`: setup attempts, raw relevant Codex command records and their index. Initial failure records are retained.
- [VALIDATION.json](VALIDATION.json): 51 independent arithmetic/alignment/structural checks, without clustering refits.
- `INPUT_MANIFEST.csv`: exact consumed/protected original identities; `DERIVED_INPUT_MANIFEST.csv`: derived inputs used for postprocessing; `OUTPUT_MANIFEST.tsv`: completed source-analysis outputs; `ARCHIVE_SHA256.tsv`: archival bytes (present in the Git archive).

## Scientific inputs and fixed decisions

Primary input is the **previously saved representative embedding** for each of 287 previously selected evaluable loci, drawn from 556 available baseline-locus medoids. All source vectors and representative/window IDs are unchanged. Prior candidate status was the union across 11/15-gene windows of Saccharide classifier positivity (score >=0.5 and >=Other) or Saccharide P@10 >=0.5. Selection was not recalculated. No prediction, P@10, functional annotation, taxonomy or UMAP coordinate was a primary clustering feature or parameter-selection criterion.

Preflight verified the exact hashes against the previous archive at commit `7b35093eefa99a4704e7b64039f51b7564788fd4`. The compatible CORE_ONLY reference has 2,502 rows and 195 any-Saccharide rows, including hybrids. Counts are 287 primary, 556 selection sensitivity, 243 selected full15, 263 unique selected vector/context classes, 2,789 global joint, and 482 Saccharide-focused joint. Inference, model retraining, the original UMAP fit, locus averaging and representative replacement were not performed.

The graph uses exhaustive float64 cosine calculations, union symmetrization, positive cosine weights, stable-ID exact ties and no self-edges. RB-configuration Leiden uses the baseline k=20/resolution=1/seed=42, plus k=10/15/20/30, resolution=0.5/1/1.5/2 and seeds 0-4. Actual optimizer calls/convergence are saved; 100 is the cap, not an assumed executed iteration count. HDBSCAN uses scikit-learn precomputed cosine distance, brute force, EOM, explicit copy=true and the self-inclusive min_samples convention. Configuration was recorded before fitting in `config.json` and `runs/run_001_science/before.json`.

All six cohorts have the same 80 grid fits plus one baseline. Two additional multiplicity-balanced graph controls each have six baseline-setting seeds; they use genome units or existing k21-sketch Jaccard-0.2 components. Relatedness components are coarse proxies, not species or ANI, and do not establish biological independence. Primary biological loci are never silently collapsed.

## Reproduction

The completed execution used the bundled Windows Python 3.12.14 with existing analysis dependencies in `analysis/bgcprophet_reference/.deps` and an isolated graph runtime specified in `config.json`. Exact versions are recorded before execution: NumPy 2.3.5, SciPy 1.18.1, scikit-learn 1.9.0, igraph 1.0.0, leidenalg 0.12.0, texttable 1.7.2, Matplotlib 3.11.1, Pillow 12.3.0. Installed dependency directories/caches are not Git artifacts. Wheel/source fingerprints and acquisition instructions are recorded separately.

A scientific rerun is a **new, separately dated execution**. Work in a new workspace rather than modifying this preserved analysis or the previous UMAP directory. Restore the archived scientific inputs to their workspace-relative locations using `ARCHIVED_INPUT_MANIFEST.tsv` in the Git archive. Restore protected previous-UMAP files from the prior archive's `original/` directory, preserving its hashes; source location is `results/WP07_viral_gt_figure3C_umap_20261024/` at the cited commit. The new archive stores the scientific inputs and immutable identities for protected prior artifacts without duplicating all prior figure/model files.

In the new copy, update only `workspace_root`, `graph_runtime` and dependency locations in its copied `config.json`, and regenerate copied input-manifest absolute paths from their workspace-relative paths without changing hashes. Install the recorded package versions into a dedicated compatible environment. `common.py` prepends the recorded analysis and graph dependency locations. A new installation/platform is not guaranteed to reproduce exact native bytes; its own package and code records must be captured. Do not use today's versions to infer the old runtime.

From the restored workspace root, using the intended Python interpreter:

```powershell
python -B analysis/viral_gt_clustering_20261024/scripts/run_recorded.py science analyze.py
python -B analysis/viral_gt_clustering_20261024/scripts/run_final_recorded.py figures figures_final.py
python -B analysis/viral_gt_clustering_20261024/scripts/run_final_recorded.py validation validate_extended.py
```

The historical scientific invocation used `run_recorded.py science analyze.py`; the newer recording wrapper adds derived-input tracking and reuses the full initial dependency fingerprint while rehashing native files and actual module entrypoints. Both exact wrapper versions and all historical source snapshots are preserved. Additional recorded authoring/saved-partition comparisons use `record_documentation.py` and `record_cross_grid.py`; their source and input hashes are captured in their own before/after records. Read those scripts before rerunning because they write new report artifacts. Scientific fits are never launched merely by importing these records.

The first figure subprocess succeeded but its wrapper caught an edit to an unused validator; its exact assets and code remain under `run_002_figures`. An initial validator exited with a syntax error before checks. The final clean render and corrected successful validation are separate runs. Logs are not rewritten to hide these attempts. There was no repeat of the scientific clustering fit.

Safe archive-only verification is `python reproduce/verify_archive.py` from the Git archive; it checks file sizes and hashes and executes no scientific scripts. New runs may replace files only in their separate reproduction copy. Preserve failed/successful attempts, stdout/stderr and output manifests there as in this record.

## Interpretation and stopping point

All five baseline communities fail at least one operational stability condition. VC03/VC04 are relatively seed-stable at the baseline setting, but are not preserved as a complete set when non-selected loci enter. Full15 and duplicate controls show that short inputs or repeated contexts do not solely explain the baseline pattern; matched-grid comparisons reveal broader parameter dependence. Global and Saccharide-focused joint partitions are separate analyses, not replacements for viral-only labels.

Candidate ascertainment and P@10 share the learned representation and Saccharide labels. The restricted reference comparison is deliberately Saccharide-enriched. Neither reference agreement nor algorithmic agreement is independent biochemical validation. Original UMAP historical gaps, including NOT RECORDED effective epochs, remain unresolved and are not relevant clustering parameters.

No detailed genomic organization analysis, new functional annotation, pathway assignment or representative-neighborhood selection was started. Stop here for joint community-structure review.

Algorithm references: [leidenalg API and objective](https://leidenalg.readthedocs.io/en/stable/reference.html), [scikit-learn HDBSCAN and min_samples convention](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html). The execution's actual versions, explicit values and code/native hashes take precedence over live documentation defaults.
