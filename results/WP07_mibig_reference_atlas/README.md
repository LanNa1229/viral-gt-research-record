# WP07: BGC-Prophet / MIBiG reference representation record

Integrated 2026-09-11 from `/home/rahnn/260908_MiBG`; this is preservation of completed work, not a new analysis run. Published-style A is retained as the reference representation for future viral-BGC placement. Placement has not been implemented here.

## Representation decision

A is the unmasked mean over all 128 contextual Transformer positions. B pools only contextualized biological core positions, independently derived by protein-ID membership and checked for original protein order and exact agreement with TDlabels (1 = core). Interrupted cores are valid: all 1,900 interrupted rows remain included. N is reconstructed as `(128*A - n*B)/(128-n)`, using float64 arithmetic on saved float32 A and B; it is not an independently extracted exact-precision tensor. Noncore context does not automatically mean native biological flanking genes or zero padding.

Core pooling did not improve augmentation stability. The median per-BGC mean of ten pairwise cosine similarities was A 0.962680, B 0.924871, N 0.961598. Corrected consensus uses all five states (-1, 0, 1, 2, 3), not a single augmentation. Any-Saccharide Precision@20 was A 0.749487, B 0.676154, N 0.740256; pure-Saccharide Precision@20 was A 0.925769, B 0.792692, N 0.924231. B had weaker Saccharide localization, while N retained A-like stability and localization. Widespread B/N perturbation cancellation was not supported; consult the cancellation and correlation tables for the measured scope.

These are supervised training-set representation analyses, not held-out biological validation. Cosine proximity does not prove shared function. UMAP is visualization only; quantitative neighbor analysis used the original 320-dimensional space. The findings support retaining A for downstream reference use without establishing performance on future viral BGCs.

## Preserved evidence

- [Published reference report](published_reference/README.md): original QC, consensus, geometry and Saccharide kNN.
- [Core extraction report](core_extraction/full/README.md) and [pilot report](core_extraction/pilot/README.md): behavior, row alignment and protected-input validation.
- [Corrected consensus report](consensus_ablation/README.md): A/B five-state comparison and reference reproduction.
- [Noncore report](noncore_mechanism/README.md): reconstruction, stability, cancellation and neighborhood comparison.
- [Consensus metadata](canonical/consensus_metadata.tsv): exact rows and augmentation lineage for the external primary matrix.
- [Payload source manifest](../../metadata/manifests/DS-VGT-20260911-001_source_files.tsv): source paths, destination paths, sizes and SHA256 values.
- [Primary matrix identity](../../metadata/manifests/DS-VGT-20260911-001_primary_atlas.json) and [external dependencies](../../metadata/manifests/DS-VGT-20260911-001_external_dependencies.tsv).

Actual historical scripts are under `scripts/archive/WP07_mibig_reference_atlas/`; existing package/version records are under `environment/WP07_mibig_reference_atlas/`. Five selected PDFs are under `figures/WP07_mibig_reference_atlas/`. Historical filenames are retained even where they predate current naming conventions. Original reports and code are copied byte-for-byte, so their local paths, relative links, commands and references to omitted arrays/figures remain historical. Use this README to navigate the curated record. This is not a self-contained runnable software distribution; recreating historical source layouts and providing external inputs is necessary. No portable utilities or replacement scripts were developed.

## Successful extraction and limitations

The preserved full launcher is `run_core_context_full_skylake.sbatch`, associated with successful job 48071887 on c5-14 (Skylake), account nknlab, four CPUs, 16 GB, elapsed 00:04:41, exit 0:0. The alternative validatednode launcher targeted unavailable c4-7 and is not selected. Full validation records 12,510 finite float32 rows of shape (12510,320), all interrupted cores retained, all protected inputs unchanged, and all 12,510 unmasked means exactly equal to the existing A embeddings (maximum absolute difference 0). Failed-run outputs and duplicate transfer bundles are excluded.

The historical published extraction source `/scratch/ln72030/BGC_Prophet/extract_bgc_embeddings.py` is historically validated but not currently recovered locally. Its expected historical SHA256 is `33e6a9eb5d07e7634b3033d96667c450be1b3b7da9dcf8417600dbdd4ff3f616`. No replacement has been created. The extracted LMDB's documented local path is currently absent; its retained archive and historical checksums are recorded separately. Checkpoint identity is a historical HPC checksum, not a new local verification. The artificial MIBiG protein-ID to amino-acid crosswalk remains unresolved; this record makes no sequence reconstruction claim.

No .npy/.npz arrays, checkpoints, LMDB files or archives, complete tensors, dependencies, environments, caches, raw logs, external source trees, or expired signed URLs are included. The primary A consensus array remains at its recorded local source; no public release location or DOI has been assigned to it. Existing source reports and provenance describe prior validations; no model inference or scientific analysis was rerun for integration.
