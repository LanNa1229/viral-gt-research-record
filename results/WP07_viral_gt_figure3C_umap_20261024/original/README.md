# Exploratory viral locus/reference UMAP — 24 October 2026

This directory contains the initial Figure 3C exploration only: one joint UMAP of existing 320D embeddings, three static views, original-space retrieval diagnostics, and locus/window/anchor audit trails. See `QC_SUMMARY.md` for verified counts and interpretation. No de novo clustering, new biosynthetic-family assignment, model inference, retraining, Git commit, or push is part of this analysis.

## Inputs and compatibility

Paths below are relative to the workspace root (`/home/rahnn/260908_MiBG` in WSL; `\\wsl.localhost\Ubuntu\home\rahnn\260908_MiBG` in Windows).

- Reference vectors: `analysis/viral_glycan_benchmark_20261016/mibig_CORE_ONLY_embeddings.npy`; metadata: `mibig_matched_reference_metadata.tsv` in that directory. These are the same 2,502 biological references used in the existing viral P@k analysis.
- Viral vectors, row manifest, window/gene manifests, anchor mapping, frozen inference caches, and baseline loci: `analysis/viral_gt_mining_20261017/runs/baseline_594_v2/`. In particular, `viral_bgc_embeddings.npy`, `embedding_manifest.tsv`, `viral_gt_windows.tsv`, `viral_gt_window_genes.tsv`, `anchor_mapping.tsv`, and `candidate_loci.tsv`.
- Frozen original anchor cohort: `analysis/viral_gt_mining_20261017/inputs/baseline_594_anchors.txt`.
- Prior 11/15-gene candidate labels: `analysis/viral_gt_figureB_taxonomy_20261018/tables/locus_audit.tsv`, cross-checked against set C in `analysis/viral_gt_window_size_20261017/tables/candidate_recovery_comparison.tsv`. These labels are preserved; the UMAP does not determine candidate status.
- Prior representative-window retrieval values for numerical reproduction: `baseline_594_v2/saccharide_localization.tsv`.
- Original five-augmentation membership/ordered-core data: `BGC_train_dataset_classify.csv`.
- Existing method evidence: `analysis/viral_gt_mining_20261017/audit_inference.md`, frozen `inference_provenance.json`, `analysis/viral_glycan_benchmark_20261016/scripts/benchmark.py`, and its reference validation/regression/pilot/ESM provenance files.

The preferred five-augmentation consensus in `analysis/bgcprophet_reference/results/consensus_embeddings.npy` is **not compatible** with the actual-gene viral inputs: it averages all-position embeddings from 128-position augmented contexts. The contextual-core consensus in `analysis/core_context_consensus_step3_20260911/` also retains contextual effects and is not isolated CORE_ONLY. Neither is used here. Their dimensions alone do not establish compatibility.

The compatible reference vectors were previously generated once per isolated biological core. The five original augmentation rows have identical ordered core genes and labels; the script verifies this against all 12,510 source rows. They are collapsed by core identity, not averaged as contextual vectors. The current task reuses those saved reference vectors without regeneration.

Both datasets use classifier checkpoint SHA256 `f45ca298980d4740dd43131a7e69c99cb233368c815470572cfc46f33abf4eb5` and ESM2 checkpoint SHA256 `46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928`. ESM2 t6/8M layer 6 uses residue means excluding BOS/EOS and the original first-1,022-residue convention. Actual ordered genes enter the original frozen Transformer, with an all-one float32 additive membership mask. The published forward computes but does not use its positional encoding. Its all-position mean equals actual-gene mean for these unpadded inputs. Saved viral vectors are checked against the original cache embeddings and membership-pooled embeddings exactly. Checkpoint/source/proof file identities are verified by SHA256.

Reference protein vectors originated in the published LMDB, whereas viral vectors were generated on CPU. Existing cross-runtime evidence supports near identity, not bitwise identity, at protein-encoding level. No new encoder regression or inference is claimed. The cached window vector identity checks are bitwise comparisons within the saved viral run.

## Biological unit and preprocessing

The original 565 overlap-merged regions remain the locus universe. The frozen `candidate_loci.tsv` establishes the window-to-locus mapping; this task neither merges nor splits regions. Existing downstream summaries provide means/maxima of scalar scores and do not define an observed representative 320D locus embedding. Consequently this task uses the requested cosine medoid.

Within each locus, normalize available 15-gene-workflow vectors in float64, calculate pairwise cosine similarities, and select the observed vector with the highest mean similarity to the other available windows. Including self-similarity would give the same mathematical ranking. A single-window locus uses its existing vector. Ties within an absolute float64 tolerance of 1e-12 use lexicographic stable window IDs. All available occurrence windows participate, including genuine alias/repeated-context observations retained by the frozen mapping; no window weights, class labels, candidate labels, or P@10 values enter this rule. `medoid_window_scores.csv` records all alternatives.

The term “15-gene workflow” includes its original boundary-truncated inputs, which use their actual number of genes. These remain documented in the representative table. All windows are verified to contain the intended contiguous source CDS indexes, without repeated-gene padding, wrapping, invented genes, or compensation for missing translations. No missing embedding is imputed and no 11-gene vector is substituted. Natural duplication of protein sequences is not interpreted as artificial padding.

UMAP receives raw original vectors promoted to a common float64 dtype, with no PCA, feature standardization, whitening, or explicit pre-normalization. UMAP may perform its own internal dtype conversion. The cosine metric handles vector direction. Original-space retrieval uses explicit float64 L2 normalization and a dot product, matching the previous cosine analysis.

## Joint map and plotting

The script fits exactly one joint UMAP, using reference embedding-row order followed by lexicographically sorted viral locus IDs:

```python
umap.UMAP(n_components=2, n_neighbors=30, min_dist=0.1,
          metric='cosine', random_state=42, n_jobs=1)
```

All remaining defaults, including resolved epoch/initialization settings, are captured in `execution_provenance.json`; package versions are in `software_versions.json`. Labels and candidate scores are used only after fitting. All three figures use the identical saved coordinates, plot limits, and scale. Reference BGCs are circles and viral loci are triangles drawn above references.

- Figure A: the seven existing annotations use explicit muted colors. Hybrid references use split-color circles, preserving all memberships without a new class or a forced primary label. Legend counts are memberships and overlap.
- Figure B: any Saccharide label, including hybrid labels, is gold; other references are light gray. Prior candidates are red, other loci are blue. The blue group includes any unresolved prior outcomes; these remain `unresolved` in exported tables and are not asserted to be negatives.
- Figure C: reference BGCs are light gray; viral triangles use viridis for representative-window Saccharide P@10, with a fixed 0–1 color bar. Larger values are drawn last for visibility.

P@10 and nearest-reference statistics are recomputed in the original 320D space against all 2,502 matched references. Reference accession breaks exact cosine ties. All plotted representative P@10 values and nearest-reference identities must reproduce the saved window-level results. No maximum across windows, union score, or 11-gene score is substituted. Reference P@10 fields are blank because Figure C measures viral-query retrieval only.

## Reproduce

The completed analysis uses the existing bundled Windows Python and existing project-local analysis packages. No package installation is required in this workspace. From the workspace root in PowerShell:

```powershell
& 'C:\Users\rahnn\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B analysis/viral_gt_umap_20261024/scripts/exploratory_umap.py
```

Use `--audit-only` to execute input checks, representative selection, and retrieval diagnostics without fitting UMAP. The script detects its workspace from its own path. Windows dependencies are resolved from `analysis/bgcprophet_reference/.deps`; on another system install the versions recorded in `software_versions.json` into the chosen interpreter. Thread settings are fixed to one. Numba and Matplotlib caches are confined to this new directory. Execution may require permission for Windows Python to write to the WSL share. A rerun overwrites only this directory's generated analysis products. All consumed original inputs are hashed before and after the complete run; no source scripts are imported for execution.

## Outputs and verification

- `figures/Figure_A_reference_classes.{png,svg}`
- `figures/Figure_B_candidate_status.{png,svg}`
- `figures/Figure_C_saccharide_P10.{png,svg}`
- `figures/preview_ABC.png`: convenient contact sheet of the three panels.
- `joint_umap_samples.csv`: required sample ID/type, full reference class, viral locus ID, representative window ID, preserved candidate status, P@10, UMAP1, UMAP2.
- `tables/locus_registry.csv`: all baseline loci, inclusion/exclusion, alternative windows and original candidate states.
- `tables/included_locus_representatives.csv`: representative metadata and 320D retrieval diagnostics.
- `tables/window_locus_anchor_crosswalk.csv`: every original window, original cohort anchor(s), source genome/window coordinates, locus, and embedding row.
- `tables/medoid_window_scores.csv`: all observed medoid candidates, scores and deterministic selection.
- `tables/representative_neighbors_top20.csv`: original-space identities, classes and cosine similarities.
- `tables/P10_group_summary.csv`, `tables/P10_group_frequencies.csv`: selected, not-selected, unresolved and combined-other distributions.
- `tables/gene_count_diagnostic.csv`: complete versus boundary-short representative counts and P@10 summaries.
- `tables/representative_15gene_embeddings.npy`, `tables/joint_umap_coordinates.npy`, `joint_umap_model.pkl`: saved representative vectors, coordinates and trusted-local fitted UMAP object. Load pickle files only from trusted sources.
- `input_and_representation_audit.json`, `execution_provenance.json`, `software_versions.json`: dimensions, counts, input gates, numerical reproduction, settings and runtime.
- `tables/input_integrity.csv`: consumed original files unchanged before/after execution.
- `output_checksums.csv`: reproducibility inventory, excluding runtime caches and this self-manifest.
- `QC_SUMMARY.md`, `visual_qa.json`: findings, limitations and rendered figure checks.

UMAP overlap/separation is descriptive and may distort high-dimensional relationships. Candidate selection and P@10 share the same learned representation and reference labels; their agreement is not independent validation. The BGC references were used in supervised model training. Variable gene counts, 128-position training versus isolated-core inference, long-protein truncation, domain shift, genuine repeated contexts across distinct biological loci, and cross-runtime ESM differences remain limitations. This initial run does not test parameter/seed robustness or establish biological families. Stop at these figures and QC for joint inspection.

Additional verification commands (same interpreter): run scripts/validate_outputs.py for independent cosine/export checks; scripts/render_saved_coordinates.py re-renders the saved coordinates without fitting UMAP. QC_DETAILS.md contains the expanded audit narrative behind the concise QC summary.
