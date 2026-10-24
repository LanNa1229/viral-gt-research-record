# Retrospective methods record: Figure 3C exploratory joint UMAP

Prepared for archival audit on 24 October 2026. This record distinguishes **EXECUTION-VERIFIED**, **RETROSPECTIVE-ONLY** and **NOT RECORDED**. The raw original execution records, not current source defaults or current packages, control the classification. No scientific rerun was performed.

See EXECUTION_EVIDENCE_MATRIX.tsv for one row per material parameter/decision; EXECUTION_LOG.md for raw-record locations and hashes.

## A. Scope and scientific objective

**Factual record.** Completed exploratory Figure 3C joint visualization; hash-matched fit code plus successful command and no clustering call. Evidence: evidence/raw/command_records.jsonl:20 (original session line 213, exec-2869ed2a-e699-4cbd-b31e-8d6ea3f67454); matrix No de novo clustering or retraining.

**Nature Microbiology draft.** To explore the distribution of viral glycosyltransferase-associated genomic loci relative to annotated reference biosynthetic gene clusters, we visualized existing learned embeddings in a joint two-dimensional UMAP. The projection was used descriptively and did not define new biosynthetic families.

**Gaps/qualification.** No additional independent biological validation was performed; no such claim is made.

## B. Biological sampling unit and datasets

**Factual record.** Execution stdout records 594 anchors, 593 mapped anchors, 757 original windows, 737 inferred vectors and 565 baseline loci. One representative from each of 556 usable loci was plotted with 2502 references. Included labels are 287 candidates and 269 not_selected; nine loci unavailable, including one candidate supported only by available 11-gene evidence. Evidence: evidence/raw/command_records.jsonl:20 (original session line 213, exec-2869ed2a-e699-4cbd-b31e-8d6ea3f67454) audit fields; original/tables/locus_registry.csv.

**Nature Microbiology draft.** The biological sampling unit was a previously defined overlap-merged viral locus. We retained 556 of 565 baseline loci with at least one usable embedding from the 15-gene-window workflow and jointly visualized these with 2,502 reference BGCs. Nine loci lacking eligible vectors were excluded without imputation or substitution by 11-gene vectors.

**Gaps/qualification.** No additional sampling/exclusion decisions outside the recorded frozen manifests are inferred.

## C. Embedding generation or reuse

**Factual record.** This run loads existing reference and viral arrays; no model inference is called. Both shapes have 320 features. The successful audit checked stored model hashes, CORE_ONLY masks/pooling/cache vectors and source augmentation-core identity. Prior inference provenance is preserved, but its original terminal execution has not been independently recovered here. Evidence: INPUT_MANIFEST.tsv; evidence/raw/command_records.jsonl:20 (original session line 213, exec-2869ed2a-e699-4cbd-b31e-8d6ea3f67454) checks_passed; executed fit audit().

**Nature Microbiology draft.** Previously generated 320-dimensional BGC-Prophet embeddings were reused. Viral actual-gene CORE_ONLY vectors were compared with the matched isolated-core reference embeddings used in the preceding retrieval analysis. Contextual five-augmentation consensus vectors were excluded because their input context differed from that of the viral vectors. No embeddings or model weights were regenerated for this visualization.

**Gaps/qualification.** Original prior-inference terminal transcript: NOT RECORDED in this archive. Prior model-generation descriptions remain attributable to the preserved prior records, not newly witnessed execution.

## D. Representative-window selection

**Factual record.** Exact recovered represent() code and independent validation record establish cosine medoids over available observed windows, self excluded, lexicographic tie-breaking within 1e-12. Singletons are copied. 130 loci have multiple original windows; 129 have multiple usable windows. Evidence: executed_code/exploratory_umap.fit_executed.py:118; evidence/raw/command_records.jsonl:30 (original session line 286, exec-65cddf17-0eae-4a37-9a33-58649f9e3104) independent_scipy_cosine_medoids.

**Nature Microbiology draft.** For loci with multiple eligible windows, an observed representative was selected as the cosine medoid, maximizing mean similarity to the other windows in the same locus. Ties within 1e-12 were resolved by lexicographic window identifier. Single-window loci used their available vector directly.

**Gaps/qualification.** Missing windows cannot influence representative selection; this limitation is preserved rather than imputed.

## E. Preprocessing and joint UMAP execution

**Factual record.** Hash-matched code concatenates raw vectors, promotes to float64, fits fit_transform(joint) with no y, and records get_params. Metric cosine, neighbors 30, min_dist 0.1, components 2, random_state 42, n_jobs 1 and init spectral are execution-verified. No explicit PCA/scaling/normalization is applied. n_epochs=null is verified; effective 500 was inserted only later. Evidence: evidence/raw/command_records.jsonl:22 (original session line 227, exec-ce9b40c7-4849-40e6-9b26-1b7980bd3661) umap_parameters; executed_code/exploratory_umap.fit_executed.py:231.

**Nature Microbiology draft.** Reference and representative viral embeddings were concatenated and embedded jointly using cosine UMAP with 30 neighbors, min_dist=0.1, two components, spectral initialization and random_state=42, using one job. No class or candidate labels were supplied to fitting, and no PCA or feature standardization preceded UMAP. All panels use the same fitted coordinates.

**Gaps/qualification.** Effective optimization epochs: NOT RECORDED. Numeric fitted kernel a/b were recovered by static inspection of the preserved fitted artifact; see evidence/model_scalar_state.json. The constructor n_epochs/a/b values were null; do not replace those facts with current-library defaults.

## F. 320D nearest-neighbor/P@10 calculations

**Factual record.** Hash-matched neighbors() uses float64 unit vectors and dot products, lexsort by descending cosine then reference accession. Any Saccharide label is positive; top-ten mean is representative-specific. All 556 values and nearest IDs reproduce saved prior results; independent SciPy validator passed. Evidence: evidence/raw/command_records.jsonl:20 (original session line 213, exec-2869ed2a-e699-4cbd-b31e-8d6ea3f67454) P10_reproduction,nearest_reproduction; evidence/raw/command_records.jsonl:30 (original session line 286, exec-65cddf17-0eae-4a37-9a33-58649f9e3104).

**Nature Microbiology draft.** Saccharide P@10 was calculated for each representative viral vector as the fraction of Saccharide-annotated BGCs among its ten nearest references in the original 320-dimensional cosine space. Hybrid labels containing Saccharide were included, and exact similarity ties were resolved by reference accession. These statistics were not calculated from UMAP coordinates.

**Gaps/qualification.** Candidate selection and P@10 are related uses of the same representation and labels; their agreement is not independent validation.

## G. Software, runtime, random seeds, and versions

**Factual record.** The actual successful command specifies the bundled Windows Python. Original pre-finalizer stdout records Python/platform and eight package versions; these were not read from today's installed environment. Seeds and thread assignments are in recorded get_params and recovered code. Evidence: evidence/raw/command_records.jsonl:20 (original session line 213, exec-2869ed2a-e699-4cbd-b31e-8d6ea3f67454) command; evidence/raw/command_records.jsonl:22 (original session line 227, exec-ce9b40c7-4849-40e6-9b26-1b7980bd3661) python/executable/platform/software_versions.

**Nature Microbiology draft.** The recorded execution used Python 3.12.14, umap-learn 0.5.12, NumPy 2.3.5, SciPy 1.18.1 and scikit-learn 1.9.0 on Windows. The UMAP random_state and transform_seed were 42; n_jobs and the explicitly assigned numerical-library thread settings were one.

**Gaps/qualification.** Complete inherited environment, transitive package lock and exact dependency binary/wheel/source checksums: NOT RECORDED. These distribution version records do not prove immutable installed package bytes.

## H. Input/output data provenance and checksums

**Factual record.** The original run checked 634 consumed files before/after and wrote the input integrity table. Finalizer command exited 0 after checking 33 artifact hashes; current archive copies match those entries. All input copies match execution-time hashes. See INPUT_MANIFEST.tsv, COPIED_FILES.tsv, original/output_checksums.csv and ARCHIVE_SHA256.tsv.

**Nature Microbiology draft.** Input files were tracked by SHA256, and the completed run recorded unchanged hashes for all 634 consumed files. Coordinate, representative, figure and provenance artifacts were retained with checksum manifests.

**Gaps/qualification.** No independent signed/timestamped digest of the original output manifest itself was emitted. External dependencies and acquisition limitations are explicitly listed; archive checksums are new retrospective integrity records.

## I. Quality-control procedures

**Factual record.** Recorded checks include dimensions/finite nonzero vectors, metadata alignment, unique source-gene indexes, mask length, cache pooling identity, reference augmentation identity, candidate crosswalks, independent medoids/P@10, CSV round-trip and export XML/resolution. Failure and correction sequence is preserved. Visual QA is an original post-fit review record, not a numerical biology test. Evidence: evidence/raw/command_records.jsonl:20 (original session line 213, exec-2869ed2a-e699-4cbd-b31e-8d6ea3f67454); evidence/raw/command_records.jsonl:24 (original session line 235, exec-a673cccd-1bd4-4401-9c4f-2f7f9829fe41); evidence/raw/command_records.jsonl:28 (original session line 261, exec-8d45d802-2b50-412c-aeb7-f83e452ff21b); evidence/raw/command_records.jsonl:30 (original session line 286, exec-65cddf17-0eae-4a37-9a33-58649f9e3104).

**Nature Microbiology draft.** Quality control verified input dimensions and identifiers, source-gene order and masks, locus crosswalks and representative-specific retrieval values. Independent cosine calculations checked medoid and nearest-neighbor results. Export checks verified coordinate round-trip, PNG resolution and vector SVG output; an SVG color-bar issue was corrected by re-rendering the saved coordinates without refitting.

**Gaps/qualification.** No color-vision-deficiency simulation was recorded. Gene-count sensitivity and additional-seed/parameter runs were not performed.

## J. Limitations and NOT RECORDED items

**Factual record.** Execution evidence supports only one initial map. Boundary-short representatives and exact-context overplotting are recorded QC findings. The effective 500-epoch claim is retrospective. The complete raw session is not redistributed because it includes non-task/system/private content; complete selected raw scientific command records are archived with exact line mapping.

**Nature Microbiology draft.** UMAP was interpreted as an exploratory projection. Unequal gene counts, repeated observed sequence contexts, prior supervised reference training, cross-runtime encoder differences and locus-representative selection may influence its appearance. Candidate-selection/P@10 agreement was not treated as independent biological validation.

**Gaps/qualification.** NOT RECORDED: actual resolved optimization epochs; full inherited environment/transitive lock; exact dependency file hashes; source-analysis Git revision; standalone analysis log; raw original prior-inference execution transcript. Full raw session is intentionally not redistributed, distinct from missing evidence.

## K. Exact mapping of scientific claims to evidence

**Factual record.** CLAIM_EVIDENCE_MAP.tsv maps six material claims to matrix evidence IDs and exact record/code locations. Visual overlap and validation-independence interpretations are RETROSPECTIVE-ONLY; execution-emitted counts and computations are EXECUTION-VERIFIED.

**Nature Microbiology draft.** Viral loci visually overlapped a subset of Saccharide-annotated references, and high-P@10 representatives occurred in that area. This descriptive agreement does not establish new biological clusters or confirm glycan biosynthesis.

**Gaps/qualification.** The visual statement is a draft Results interpretation, not a new quantitative analysis or validated family assignment.
