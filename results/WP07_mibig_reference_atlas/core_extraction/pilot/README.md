# 20-row core-context behavioral pilot — PASS

Slurm job 48067092 completed with exit code 0:0 in 00:02:04. Only the supplied 20-row pilot was evaluated. No full 12,510-row extraction was launched.

## Results

| Check | Result |
|---|---|
| Encoder output before final mean | (20, 128, 320) |
| Core-context output | (20, 320), float32 |
| Encoder and pooled values finite | 20/20 PASS |
| Selected contextual positions equal sentence length | 20/20 PASS |
| Hooked/unhooked classifier probabilities exactly equal | 20/20 PASS |
| Maximum hooked/unhooked probability difference | 0.0 |
| Captured unmasked mean equals classifier input | 20/20 PASS |
| Interrupted-core rows included | 5/20 |
| Original input/source checksums unchanged | PASS |
| Retrieved-file checksums and local validation | 80/80 PASS |

The completed membership audit was reused, not recalculated. Its training CSV and metadata SHA256 values were checked against HPC inputs. That audit established core=1/noncore=0 for TDlabels, exact agreement on all 12,510 rows, and 10,610 contiguous versus 1,900 interrupted rows across 1,309 BGCs.

Core selection uses only protein-ID membership between sentence and TDsentence. The saved mask retains original protein order and has exactly len(sentence) positions. Pooling is H[j, mask[j]].mean(dim=0), with no individual vector normalization. Start/end bounds are not used for interrupted cores.

The supplied HPC pilot CSV was checked row-for-row against the full source CSV, including sentence, TDsentence and TDlabels. The manifest and metadata record exact source_embedding_row values.

## Validated execution

The existing HPC extract_bgc_embeddings.py, installed classifier source, and prior Slurm launcher were inspected first. The new script constructs the same geneClassifier and invokes its original clf.classify() method for both passes, using CPU and batch_size=128. Both passes therefore have the same 20-row batch. No inference-context, attention-mask, installed code, or precision changes were introduced. Recorder hooks clone detached tensors and return no output replacement; hooks are removed afterward.

The model state was checked for equality after the passes. SHA256 checks cover the checkpoint, LMDB data.mdb, training/pilot CSVs, original metadata, original extraction script, and installed classifier/data source files. All were unchanged. The full contextual tensor was not saved.

Six transferred script/audit files passed checksums before submission; 16 retrieved files passed checksums afterward. The remote job logs are retained. Existing loader warnings about nested tensors and five workers with four allocated CPUs did not prevent successful completion; the installed loader was preserved.

## Additional cross-run comparison

The captured all-position mean embeddings exactly match the corresponding existing full-run embeddings (maximum absolute difference 0.0).
The unhooked pilot probabilities differ from the previous full-run probabilities by at most 5.820766091346741e-11. This separate cross-run comparison is recorded in results/existing_full_run_comparison.tsv. It is not the behavioral hook test: the pilot's hooked versus unhooked difference is exactly 0.0.

## Files

- results/mibig_20_core_context_embeddings.npy
- results/mibig_20_core_context_metadata.tsv
- results/pilot_behavior_validation.tsv
- results/pilot_manifest.tsv
- results/pilot_summary.json
- results/mibig_20_probabilities_unhooked.npy
- results/mibig_20_probabilities_hooked.npy
- results/mibig_20_unmasked_mean_embeddings.npy
- results/existing_full_run_comparison.tsv
- results/input_integrity_before.tsv and input_integrity_after.tsv
- results/reused_audit_provenance.json
- results/runtime_provenance.json
- local_verification.tsv and local_verification_summary.json
- OUTPUT_SHA256SUMS, job_status.tsv, pilot_48067092.out, pilot_48067092.err

Scripts are in the project scripts/ directory:
run_hpc_core_context_pilot.py, extract_core_context_pilot.py, run_core_context_pilot.sbatch, verify_retrieved_core_pilot.py.
The exact transferred versions and SHA256 manifest are preserved in outputs/core_context_pilot_submission/.
Runtime: Python 3.10.20, PyTorch 2.11.0+cu128 (CPU execution), NumPy 2.2.6, pandas 2.3.0, LMDB 2.2.1.
No new random seed or inference settings were imposed on the validated loader.

Remote directory: /scratch/ln72030/BGC_Prophet/core_context_pilot_20260910
Submission used the prepared run_core_context_pilot.sbatch through the existing WSL ControlMaster socket /home/rahnn/.ssh/cm/sapelo2. SSH/SCP used this socket and ProxyCommand=false to prevent authentication fallback. No other job was submitted.

STOP: the 20-row pilot is complete. Full-run inference remains unlaunched.

