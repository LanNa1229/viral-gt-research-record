# Baseline 594-member viral GT neighborhood run

The original WP01B cluster A0A3Q8Q3U2 contains exactly 594 unique anchors, including the representative.
Source: `/home/archival-user/projects/viral-gt-research-record/results/WP01B_coarse_clustering/canonical/cluster_membership.tsv`.
Headerless column 1 is representative; column 2 is accession. Baseline SHA256: `cd731e42c1bffd7946427540d4cb61d72f706f494a93df2a8be11831f554fd2f`.
This is the original baseline, not a spike-in rerun. The seven-anchor pilot is a separate cohort.

Latest full run: [`runs/baseline_594_v2`](runs/baseline_594_v2/).
[Commands, dependencies, methods, output schemas and interpretation](checkout/workflows/viral_gt_mining/README.md).
[Frozen inference audit](audit_inference.md). [Exact cohort receipt](inputs/cohort_validation.json).

## Executed and verified

- Input anchors: 594; mapped anchors: 593; unresolved anchors: 1.
- Occurrence-specific windows: 757; complete 15-gene windows: 616; gene rows: 10688.
- Frozen inferences: 737; excluded windows: 20. Overlapping windows remain separate.
- Independent validation: 17029 checks passed; 380 historical artifacts retain their recorded content hashes.
- Extraction: 19.38 s; inference/encoding: 38.31 s; peak inference RSS: 646.6 MiB.
- ESM unique proteins: 7092; new encodings: 871; cached: 6221.
- Development branch: `dev/viral-gt-windows-20261017`. Original repository status is clean. No push/publication performed.

## Classification counts

| Unit | Total | Analyzed at least once | Positive | Fully analyzed negative | Unknown | Fraction among resolved |
|---|---:|---:|---:|---:|---:|---:|
| window | 757 | 737 | 246 | 491 | 20 | 33.4% |
| locus | 565 | 556 | 187 | 367 | 11 | 33.8% |
| genome | 466 | 459 | 166 | 290 | 10 | 36.4% |

A locus/genome is positive if any analyzed constituent window meets the original classifier Saccharide threshold. A negative requires all constituent windows to be analyzed; otherwise its outcome remains unknown.
Identical nucleotide-sequence hash plus source taxid collapses GenBank/RefSeq aliases for genome counts. Loci merge only connected overlapping nucleotide intervals after inference.
Unmapped anchors have no known genome and remain outside genome denominators. Fractions describe this selected structural cohort, not all-virus prevalence. Wilson intervals are descriptive because observations overlap and were selected.

## Mapping and exclusion scope

All remote mappings are bounded by explicit UniProt/UniSave genomic cross-references and exact source-CDS protein sequence matches; this is not an exhaustive public-genome/IPG search.
Every annotated CDS is retained in genomic order. Missing/uncertain translations exclude inference without shifting gene positions. Linear/contig ends remain incomplete; no wrapping or compensation occurs.

- 10 excluded windows: `translation_mismatch`
- 10 excluded windows: `translation_unverified`

Unresolved anchor details are in `anchor_mapping.tsv`, supplemented by `inputs/retrieval_anchor_status.tsv`.

## Pilot and model interpretation

The pilot had 7/7 mapped and inferred windows, six complete and one 11-gene boundary case; 544 independent pilot checks passed. S-H68 QBP05883.1 gives Saccharide 0.999795. The five Catovirus windows span Saccharide scores 0.024441–0.996528, with only two receiving the Saccharide label. The Catovirus ARF08732.1-centered window is not the historical ARF08726–ARF08739 input; see `runs/pilot_v1/manual_region_comparison.tsv`.
Primary representation uses actual window genes, no artificial padding/repetition, the exact frozen checkpoint and the compatible 2,502-BGC CORE_ONLY reference. Seven sigmoid outputs need not sum to one. Original all-actual-position mean pooling equals actual-gene pooling in this input. Reference startup identity/proof gates protect comparability.
Classification and 320D cosine/P@k retrieval remain separate. P@10 enrichment is descriptive over 195/2502 Saccharide-labeled references; no binary localization cutoff or taxonomic enrichment test was chosen. Neither result confirms glycan biosynthesis.
Proteins above 1,022 aa retain complete source sequences but use the original first-1,022-residue encoder convention. Viral domain shift and supervised training-reference geometry remain limitations.

## Files

- [Anchor mappings](runs/baseline_594_v2/anchor_mapping.tsv)
- [Window manifest](runs/baseline_594_v2/viral_gt_windows.tsv)
- [Ordered window genes](runs/baseline_594_v2/viral_gt_window_genes.tsv)
- [Full source protein FASTA](runs/baseline_594_v2/viral_gt_window_proteins.faa)
- [Predictions](runs/baseline_594_v2/viral_bgc_predictions.tsv)
- [Embedding manifest](runs/baseline_594_v2/embedding_manifest.tsv)
- [320D localization](runs/baseline_594_v2/saccharide_localization.tsv)
- [Candidate loci](runs/baseline_594_v2/candidate_loci.tsv)
- [Taxonomy summaries](runs/baseline_594_v2/viral_taxonomy_summary.tsv)
- [Independent validation](runs/baseline_594_v2/validation.json)
- [Output checksums](runs/baseline_594_v2/artifact_checksums.tsv)

## Figures

![Reference/viral overlay](runs/baseline_594_v2/reference_viral_embedding_overlay.png)
![Saccharide-associated fractions](runs/baseline_594_v2/saccharide_associated_fractions.png)
![Taxonomic distribution](runs/baseline_594_v2/viral_taxonomic_distribution.png)

PNG exports are 300 dpi and paired PDFs are vector. Visual review status is recorded alongside the figures. Quantitative retrieval uses the original 320D space, never PCA coordinates.
