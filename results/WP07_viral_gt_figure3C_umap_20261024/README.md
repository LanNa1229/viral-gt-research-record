# WP07 Figure 3C: retrospective archive of the initial viral GT joint UMAP

This dedicated result archive preserves the completed `/home/rahnn/260908_MiBG/analysis/viral_gt_umap_20261024/` analysis. It follows the existing `results/WP07_*` convention and keeps the methods, evidence, original code/results/figures and reproducibility records together. No root README, unrelated working-tree files or original analytical files were edited. No new UMAP, clustering, retraining, inference, or scientific metric computation was performed for this archive.

## Start here

- [METHODS_RECORD.md](METHODS_RECORD.md): factual record, Nature Microbiology draft, and explicit gaps in each subsection A–K.
- [EXECUTION_EVIDENCE_MATRIX.tsv](EXECUTION_EVIDENCE_MATRIX.tsv): 121 parameter/decision entries with EXECUTION-VERIFIED, RETROSPECTIVE-ONLY or NOT RECORDED classifications.
- [EXECUTION_LOG.md](EXECUTION_LOG.md): exact original commands, successful and failed attempts, raw log excerpts, line mappings and hashes.
- [CODE_IDENTITY.tsv](CODE_IDENTITY.tsv): recovered executed code versus current preserved source.
- [CLAIM_EVIDENCE_MAP.tsv](CLAIM_EVIDENCE_MAP.tsv): six scientific statements mapped to evidence IDs.
- [Figure A/B/C preview](original/figures/preview_ABC.png), [coordinate CSV](original/joint_umap_samples.csv), [QC summary](original/QC_SUMMARY.md).

The successful fit has script SHA256 `1c5164f368955fc36eec4bede5539be6e8fc7a89fc62a6feab42eed415c88b25`. The recovered fit-time snapshot is [executed_code/exploratory_umap.fit_executed.py](executed_code/exploratory_umap.fit_executed.py). The preserved current main script has SHA256 `c395d57a550494b62e679b629ddfe29382ecca817dd96299ae925aa988aac060` and was used for final figure rendering after a documented color-bar patch. This distinction is explicit; [the diff](evidence/fit_to_render_code.diff) is preserved. The validator, corrected renderer and finalizer are also recovered byte-for-byte from their executed write commands.

## Important retrospective correction

The original README/QC and final edited provenance contain a 500-epoch claim. It came from post-run inspection of installed UMAP source followed by a finalizer edit. The original execution-emitted constructor state records `n_epochs=null`; the actual effective optimization epoch count is **NOT RECORDED**. The numerical fitted kernel a/b values, in contrast, are available in the preserved fitted-model artifact and were recovered by static opcode inspection without unpickling or recomputation. This archive preserves original files unchanged and corrects interpretation in the evidence matrix and new Methods record.

## Contents and integrity

`original/` contains all 34 non-cache source files, byte-for-byte: the four current scripts, original README and QC, coordinate CSV, representative/medoid/crosswalk and nearest-neighbor tables, original manifests/provenance, final PNG/SVG panels, saved representative/coordinate arrays and fitted model. Original Matplotlib/Numba caches are omitted as runtime rebuildable caches; this is explicitly separate from frozen inference caches consumed as audit inputs.

`inputs/workspace/` preserves 629 of the 634 consumed input files, including both reference and viral arrays, metadata, window/gene/anchor/locus registries, prior proof records and 609 frozen inference-cache files. [INPUT_MANIFEST.tsv](INPUT_MANIFEST.tsv) gives every consumed path, size, execution/current SHA256, availability and acquisition route. Five files remain external: the 31.9-MB training CSV, classifier and ESM checkpoints, and two upstream source files for which redistribution permission was not established. They remain locally available with exact hashes and authoritative source pointers in [EXTERNAL_DEPENDENCIES.tsv](EXTERNAL_DEPENDENCIES.tsv). The 15-GB LMDB and raw genome collection were not consumed by this UMAP run and are not necessary for replaying the saved-vector analysis; regeneration of prior embeddings is a separate workflow.

[COPIED_FILES.tsv](COPIED_FILES.tsv) verifies every source-to-archive copy. [ARCHIVE_SHA256.tsv](ARCHIVE_SHA256.tsv) covers all archive files except itself; [AUDIT_VALIDATION.json](AUDIT_VALIDATION.json) records archival checks. The original output manifest is preserved unchanged under `original/`; it is distinct from the new archival manifest. The archive-scoped `.gitattributes` disables text conversion so Git preserves original CRLF/LF bytes. Its `.gitignore` makes narrow exceptions for the small, hash-verified NPY/NPZ/PKL scientific artifacts; global repository ignore rules were not modified.

`evidence/raw/` contains exact original raw scientific command records and relevant tool invocation lines, not a fabricated execution log. Full local session material containing system instructions, encrypted reasoning, unrelated setup or image payloads is intentionally not redistributed. [evidence/record_index.tsv](evidence/record_index.tsv) retains the original session line numbers; full source-prefix identity and exclusion scope are in [evidence/session_identity.json](evidence/session_identity.json). Timestamps identify log records and are not treated as independent proof of the wall clock.

## Reproduction instructions — future run only

The archival checks do not run scientific scripts. To prepare a separately dated future reproduction, use a trusted Python interpreter to restore a **new** workspace:

```bash
python reproduce/restore_layout.py --workspace /absolute/path/to/a/new/figure3c-reproduction --code-stage fit
```

This copies files only and refuses an existing workspace. Retrieve the five external files to their exact workspace-relative locations and verify their SHA256 values from `EXTERNAL_INPUTS_TO_ACQUIRE.tsv`. Use the GitHub blob API's raw-media response (or decode the blob content) for exact upstream source/CSV bytes. For the classifier, extract `classifier.pt` from the recorded official model archive. Do not substitute a latest checkpoint. Do not imply missing inputs can be regenerated from UMAP outputs.

Recreate the recorded Python 3.12.14 environment and distribution versions in `original/software_versions.json`. On Windows the original script prepends `analysis/bgcprophet_reference/.deps`; populate that path consistently or use an isolated interpreter containing the recorded distributions. The complete original transitive environment and wheel hashes are NOT RECORDED, so cross-platform or newly resolved environments are not guaranteed bitwise replicas.

From the restored workspace root, the historical command was:

```powershell
& 'C:\Users\rahnn\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B analysis/viral_gt_umap_20261024/scripts/exploratory_umap.py
```

The recorded command invokes the original fixed parameter set and no script arguments. Only after separately authorizing a future scientific rerun should that command be executed. The `--audit-only` option still recalculates medoids/retrieval and is not an archival integrity command. A fresh run must have a new date, logs and hashes; it cannot fill gaps in the original execution record.

The fit snapshot originally creates the initial figures before the SVG color-bar correction. To reproduce final rendering, restore with `--code-stage render` instead, then run `scripts/render_saved_coordinates.py` against the restored saved coordinate artifacts. It re-renders without fitting, but remains a future rendering operation. The original `validate_outputs.py` recalculates scientific metrics as a validation step; it was not rerun during this archive task. Our archive-only checker uses hashes, structure and syntax, not scientific recomputation. Fitted pickle files are preserved for provenance and were not deserialized; load them only in a trusted environment.

Original relative links, absolute local paths and commands inside `original/` and prior records are historical and can require restoration to the original layout. They were deliberately not rewritten. This new archive README and the manifest paths are the navigation/reproduction contract.

## Unresolved details

- Actual effective optimization epoch count: NOT RECORDED (500 remains RETROSPECTIVE-ONLY).
- Complete inherited environment, transitive package lock and dependency source/wheel/binary hashes: NOT RECORDED.
- Source-analysis Git revision and a standalone analysis log file: NOT RECORDED; exact code hashes and raw Codex command records are available.
- Raw original prior-inference execution transcript: NOT RECORDED in this archive; saved prior-generation proofs and frozen inputs are retained.
- Independent signed/timestamped attestation of the original checksum manifest: NOT RECORDED.

The figure interpretation remains exploratory. Candidate selection and P@10 are not independent validation, and UMAP separation does not define biosynthetic families. Git archival completion does not remove these scientific limitations.

Safe archival verification: run python reproduce/verify_archive.py from any working directory. This uses only file sizes and SHA256 and executes no analytical code. PARAMETERS.execution_verified.json separates recorded constructor values from resolved-state gaps.
