# Viral glycan BGC classification and embedding benchmark

Completed pilot for the single user-specified locus on 2026-10-16. No model retraining, SaProt inference, or original-data edits were performed. Results are model recognition and supervised-reference retrieval, not biochemical validation.

## Result

The manually selected S-H68 locus is classified **Saccharide** by both available input strategies. The predicted label and rank are stable between REPEAT_128 and CORE_ONLY. Agreement across all three strategies cannot be assessed because REAL_128 is unavailable. No other viral locus was confirmed by the user for this run; broader cross-locus consistency and biological negative-control performance are therefore not estimated.

| Viral BGC | Strategy | Saccharide score | Rank | Predicted labels | Saccharide P@10 |
|---|---|---:|---:|---|---:|
| S-H68 QBP05878–QBP05897 | REAL_128 | NA | NA | Unavailable: partial contig has 117 genes | NA |
| S-H68 QBP05878–QBP05897 | REPEAT_128 | 0.995225 | 1 | Saccharide | 100% |
| S-H68 QBP05878–QBP05897 | CORE_ONLY | 0.995182 | 1 | Saccharide | 100% |

| Strategy | P@5 | P@10 | P@20 |
|---|---:|---:|---:|
| REPEAT_128 | 100% | 100% | 95% |
| CORE_ONLY | 100% | 100% | 95% |

The within-locus REPEAT_128 / CORE_ONLY embedding cosine is **0.999996485**. This is an input-sensitivity diagnostic; each query's reported neighbors come from its own matched reference space. Saccharide-positive prevalence is 195/2,502 (7.79%) in the biological reference set, counting hybrid labels as positive. P@k is a descriptive positive-neighbor fraction, not a calibrated probability of function.

## Curated input and sequence validation

The user explicitly supplied **QBP05878.1–QBP05897.1 inclusive**: 20 proteins. Expected glycan biosynthesis is the manual annotation from the task, stored separately from model outputs. No prediction was used to choose the core genes.

NCBI protein QBP05878.1 resolves to **Synechococcus phage S-H68, MK016663.1**, an explicitly **partial genome** of 79,639 nt. The available contig has 117 translated CDSs, all with unique protein IDs, valid coordinates, and translations independently verified against their CDS nucleotide sequence. The core spans **22,164–38,671**, 1-based inclusive, and is retained in ascending genomic order. This does not establish the gene count of the complete virus genome.

REAL_128 is unavailable because the source contig has fewer than 128 protein-coding genes. It has an explicit unavailable manifest row and no fabricated embedding or score row. REPEAT_128 is synthetic, repeats the 20-gene ordered core and truncates to 128: the first eight genes appear seven times, and the remaining twelve six times. CORE_ONLY uses exactly the 20 original genes. No core genes are dropped by either available strategy. All-position and core-position means coincide here because every input position is core; the two arrays are retained separately.

Protein embeddings use the existing esm2_t6_8M_UR50D checkpoint, layer 6, evaluation mode, float32 CPU, residue mean excluding BOS/EOS, and the published 1,022-residue truncation. There are 0 truncated pilot proteins. Exact sequence hashes key the cache so repeated positions and strategies reuse the same vector.

Primary sequence record: [NCBI MK016663.1](https://www.ncbi.nlm.nih.gov/nuccore/MK016663.1). Download URLs, hashes, and retrieval date are in download_provenance.tsv. Complete coordinates, translations, accession identities and source hashes are retained in the manifest tables and source GenBank file.

## Frozen model behavior

The official model archive was downloaded into this analysis because Sapelo2's noninteractive SSH connection failed. Its classifier.pt SHA256 is `f45ca298980d4740dd43131a7e69c99cb233368c815470572cfc46f33abf4eb5`, exactly matching the previous validated HPC checkpoint. The local published classifier source also matches the historical HPC source SHA256. Exact source modules are loaded without editing or installing BGC-Prophet.

The original Transformer forward is used for all inputs. Positional encoding is computed but the original unencoded input is passed into the Transformer, preserving published behavior. The published float32 membership mask is retained: values 1 for core and 0 for flanks are **additive key-attention biases**, not boolean exclusion masks. The primary embedding is the all-position mean supplied directly to the final product classifier; an encoder hook also computes a core-position mean. Variable-length inference uses the existing forward directly, without padding CORE_ONLY to 128 and without modifying weights.

Seven outputs are raw sigmoid scores in the published order: Alkaloid, Terpene, NRP, Polyketide, RiPP, Saccharide, Other. They need not sum to one. The threshold is `max(0.5, Other score)`; all scores meeting it are selected, Other is removed when another label also qualifies, and an empty set becomes Unknown. Saccharide rank is 1 plus the number of scores strictly greater than it. The `isBGC=No -> NonBGC` reporting gate is bypassed by direct inference. Thresholds and weights are unchanged.

## Matched MIBiG reference and projections

REAL_128 retains the original 2,502-BGC consensus reference, independently verified as the arithmetic mean of each biological BGC's five original 128-position embeddings. No viral query is plotted in that panel because the strategy is unavailable. This original augmentation-context reference is not claimed to be an authentic-flank distribution match in every respect.

REPEAT_128 and CORE_ONLY references were regenerated from the original published LMDB protein vectors using each BGC's original ordered `sentence` core list. Five augmentations are checked for identical cores and labels, then collapsed to one biological input. Each matched space contains exactly 2,502 biological BGCs. No reference core exceeds 128 genes, so REPEAT_128 does not omit any reference core gene. The existing historical core-context pooling arrays were **not** mistaken for these newly generated input strategies.

Nearest neighbors use cosine similarity in 320 dimensions, with lexicographic accession tie-breaking. The complete top 20 includes top-5/top-10 flags and known labels. PCA is fit on raw, centered reference vectors only, without scaling or whitening, and queries are transformed afterwards. UMAP is fit separately on each reference set: cosine metric, 15 neighbors, min_dist 0.1, 500 epochs, spectral initialization, random_state and transform_seed 1729, one thread. Queries are transformed afterwards; no joint fit includes them. Axes from distinct strategy panels are not directly comparable. Quantitative retrieval never uses 2-D plot distance.

| Strategy | Rank | MIBiG BGC | Known label(s) | Cosine similarity |
|---|---:|---|---|---:|
| REPEAT_128 | 1 | BGC0000778 | Saccharide | 0.956166 |
| REPEAT_128 | 2 | BGC0000753 | Saccharide | 0.954243 |
| REPEAT_128 | 3 | BGC0000723 | Other Saccharide | 0.951991 |
| REPEAT_128 | 4 | BGC0000758 | Saccharide | 0.951951 |
| REPEAT_128 | 5 | BGC0000752 | Saccharide | 0.948575 |
| CORE_ONLY | 1 | BGC0000778 | Saccharide | 0.955848 |
| CORE_ONLY | 2 | BGC0000753 | Saccharide | 0.953843 |
| CORE_ONLY | 3 | BGC0000758 | Saccharide | 0.952108 |
| CORE_ONLY | 4 | BGC0000723 | Other Saccharide | 0.951728 |
| CORE_ONLY | 5 | BGC0000752 | Saccharide | 0.948008 |

![PCA comparison](pca_comparison.png)

![UMAP comparison](umap_comparison.png)

## Quality controls and limits

- Pilot executed before scaling to all 2,502 reference BGCs. A 20-BGC matched-reference pilot and 20 original-row regression passed before full reference generation.
- Hooked and unhooked pilot scores are bitwise identical. Pooled hook output equals the original classifier input. Model weights and checkpoint remained unchanged.
- All saved arrays have aligned sample identities, expected shapes, finite values, and nonzero vectors for cosine comparisons. REAL_128 has no tensor because it is unavailable; REPEAT_128 is (1,128,320), CORE_ONLY (1,20,320).
- Original MIBiG regression maximum absolute differences: embedding **1.2457371e-05**, score **4.2319298e-06**. Cross-runtime tolerance was predeclared as atol=rtol=1e-4; bitwise reproduction of historical runtime outputs is not claimed.
- Uniform-repeat tests using actual viral protein subsets of lengths 1,2,4,8,16 passed at atol=rtol=1e-5. These are numerical controls, not additional curated loci. Across 448 actual MIBiG BGCs whose core lengths divide 128, B/C embeddings pass atol=rtol=1e-4, maximum absolute difference 2.5138259e-05. The 20-gene S-H68 core does not divide 128, so exact equivalence is not expected.
- Original CSV, saved raw MIBiG embeddings/scores, validated extractor, ESM checkpoint, and full 15.1-GB LMDB match their historical SHA256 hashes. LMDB was opened readonly with lock=False/create=False, and its size/mtime stayed unchanged during inference. Original consensus arithmetic also reproduces exactly.
- Historical Catovirus extraction regression outputs were not found in the inspected project resources, so no Catovirus regression is claimed. The saved original MIBiG regression provides a direct checkpoint/extraction compatibility check for this run.
- CORE_ONLY lengths and synthetic repeats differ from the classifier's original 128-gene training construction; viral proteins and the partial-genome source are additional distribution limitations. Matched reference construction controls one input-distribution difference, not all domain shift.
- The MIBiG set was used in training the supervised classifier. Neighborhood enrichment is not held-out performance or independent validation of glycan chemistry.
- Only one viral locus is specified. `viral_pairwise_similarities.tsv` is intentionally header-only: no distinct viral BGC pair exists. No curated non-glycan viral control was supplied, and none was inferred from model predictions.

## Files and reproduction

All nine requested data/figure outputs plus this README are present. The classification and viral embedding files contain the **two available** strategies; input_strategy_manifest.tsv and strategy_comparison.tsv explicitly include the unavailable third strategy. Embedding shape is (2,320), with matching embedding_metadata.tsv. Full seven-score vectors are in classification_scores.tsv. Additional manifests, coordinate tables, reference arrays, numerical controls, download provenance, source checksums and package versions support audit and reruns.

From `/home/rahnn/260908_MiBG` under WSL, use the existing environment (no new packages were installed):

```sh
PY=analysis/bgcprophet_sequence_recovery_20260912/.venv/bin/python
TASK=analysis/viral_glycan_benchmark_20261016
$PY -B $TASK/scripts/benchmark.py prepare
$PY -B $TASK/scripts/benchmark.py esm
$PY -B $TASK/scripts/benchmark.py pilot
$PY -B $TASK/scripts/benchmark.py reference --limit 20
$PY -B $TASK/scripts/benchmark.py reference
```

Then run `scripts/compare.py` from this directory with the existing Windows bundled Python at `C:\Users\rahnn\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe`. It reuses `analysis/bgcprophet_reference/.deps` read-only, sets caches inside this analysis directory, and writes the tables and figures. Exact plotting package versions are in comparison_validation.json.

Finally, in WSL:

```sh
$PY -B $TASK/scripts/validate.py --final
$PY -B $TASK/scripts/report.py
```

The scripts reuse source files in this workspace and the downloaded inputs retained in `sources/`; resource_inventory.tsv lists those dependencies. `benchmark.py` refuses a classifier whose checksum differs from the validated HPC checkpoint and refuses reference generation before the viral pilot passes. All writes target this new analysis directory. Reproduction overwrites generated outputs here. There was no HPC execution or job submission, so an HPC transfer package was not needed.

To extend to additional curated viral loci, supply their explicit ordered core definitions and genome records, then extend the manifest preparation and rerun a validated pilot before scaling. Candidate historical FASTAs were not silently treated as new manual selections.
