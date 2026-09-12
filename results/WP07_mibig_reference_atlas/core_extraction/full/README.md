# Full core-context extraction — validated outputs

Use **results_validatednode/** for the completed dataset.
The results/ directory preserves the unsuccessful first attempt's diagnostics; it contains no released embedding array.

## Final validation

Slurm job **48071887**, node c5-14 (Intel Skylake), completed with exit code 0:0 in **00:04:41**.

| Requirement | Result |
|---|---|
| Embedding rows | 12,510 |
| Shape | (12510, 320) |
| dtype | float32 |
| All values finite | PASS |
| Exactly len(sentence) contextual positions pooled | 12,510/12,510 PASS |
| Interrupted rows retained | 1,900 |
| Metadata alignment with existing published-style dataset | Exact |
| Unmasked means equal published-style embeddings | 12,510/12,510 exactly equal |
| Maximum unmasked-reference absolute difference | 0.0 |
| Corresponding core vectors equal successful 20-row pilot | Exactly equal; max difference 0.0 |
| Protected original files unchanged | 33/33 SHA256 checks PASS |

All retrieved file checksums passed. Independent local verification covers every row's identity, metadata, ordered core membership, selected count, retained interruptions, finite values, and saved validation flags. See local_validation.tsv and local_validation_summary.json.

## Deliverables

- results_validatednode/mibig_12510_core_context_embeddings.npy
- results_validatednode/mibig_12510_core_context_metadata.tsv
- results_validatednode/core_context_full_validation.tsv
- results_validatednode/full_summary.json and full_summary.tsv
- results_validatednode/input_integrity_before.tsv and input_integrity_after.tsv
- results_validatednode/reused_audit_provenance.json
- results_validatednode/runtime_provenance.json
- full_retry_48071887.out and full_retry_48071887.err
- job_status.tsv and final_scheduler_record.txt
- OUTPUT_SHA256SUMS

The complete contextual (12510,128,320) tensor was neither accumulated nor saved. Only pooled core vectors and numerical validation/provenance records were retained for the full extraction.

## Method and provenance

The completed membership audit was reused with its CSV and metadata SHA256 bindings. Core=1/noncore=0 TDlabels polarity had already been established independently; the pooling mask comes from protein IDs, not TDlabels. Filtering TDsentence by this membership mask reproduces sentence in exact order, including interrupted cores.

The script uses the validated geneClassifier constructor and original clf.classify() method. It does not reimplement forward, change the model input mask, apply individual-vector normalization, or modify dtype/inference context. CPU, batch_size=128, account nknlab, four CPUs and 16 GB match the validated resource strategy.

A read-only encoder forward hook clones and detaches each batch, computes H[j, derived_membership_mask[j]].mean(dim=0), and releases contextual data after pooling. A classifier pre-hook checks exact equality of the same encoder output's unmasked mean with the classifier input. The model state remains unchanged. All 12,510 unmasked means are compared against the saved published-style array during this same run. No second full baseline inference was performed.

Original CSV, metadata, checkpoint, LMDB data.mdb, reference embeddings, extraction script and all installed package .py source files were hashed before and after. The relevant protected inputs were also checked against the successful pilot. Prior analysis outputs were not modified.

Runtime/package versions and script hashes are in results_validatednode/runtime_provenance.json. No new random seed, precision or backend settings were imposed on the original inference path. The original loader's worker/nested-tensor warnings are retained in stderr.

## CPU reproducibility investigation

The first full attempt, job 48067848 on AMD EPYC Rome node a5-9, completed all rows but failed the reference check: max unmasked-mean difference 0.0005353093147277832, with zero exact matches. It passed shape, finite-value, count and protected-file integrity checks, but released no embedding file. Its results and logs are preserved.

A 20-row diagnostic, job 48071886 on a5-9, confirmed that hooked and unhooked probabilities remain exactly identical (difference 0.0) while cross-run comparison against the original saved outputs differs slightly. This supports CPU-dependent numerical differences; it does not prove the specific low-level implementation cause. No tolerance was relaxed.

The full retry, job 48071887, was scheduled on the successful pilot's CPU family (Intel Skylake) and passed with exact equality to all published-style embeddings. It used the same extraction script and comparison criterion (rtol=1e-5, atol=1e-6; actual differences were zero).

The initially requested exact pilot node c4-7 was unavailable. While pending, the retry's node restriction was cleared and Features=Skylake was set; it ran on c5-14. The small AMD diagnostic's memory request was reduced to 8192 MB while pending. These scheduling changes did not modify model inference.

## Scripts and execution records

Project scripts:
- scripts/run_hpc_core_context_full.py
- scripts/extract_core_context_pilot.py (shared read/audit helpers; audit was not rerun)
- scripts/run_core_context_full.sbatch (first attempt)
- scripts/run_core_context_full_validatednode.sbatch (submitted retry before scheduler override)
- scripts/run_core_context_full_skylake.sbatch (reproducible final CPU-family constraint)
- scripts/run_core_context_node_diagnostic.sbatch
- scripts/verify_retrieved_core_full.py

The initial transfer bundle and SHA256 manifest are preserved in outputs/core_context_full_submission/. The two full attempts used the same Python extraction script hash. The final CPU-family launcher is provided locally for reproducibility; its existing output-directory guard intentionally prevents overwriting these results on rerun.

Remote validated directory:
/scratch/ln72030/BGC_Prophet/core_context_full_20260910/results_validatednode

SSH/SCP reused the authorized ControlMaster socket /home/rahnn/.ssh/cm/sapelo2, with authentication fallback disabled. Slurm stdout/stderr and accounting records for all attempts and the small diagnostic were retrieved.

No PCA, UMAP, HDBSCAN, kNN, training, or downstream analysis was performed.

