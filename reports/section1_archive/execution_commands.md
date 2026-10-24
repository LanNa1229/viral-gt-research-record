# Execution commands and scope

The portable validator and reconstructed plot command are in the root README. All outputs must be directed to a new working directory. Default numerical checks require only Python standard library; bootstrap checks require NumPy; plots require NumPy and Matplotlib>=3.9.

## Safely tested historical postprocessors

1. `build_complete_foldseek_cytoscape_network.py --membership MEMBERSHIP.tsv --m8 aln_result.m8 --outdir NEW_NETWORK_DIR`. Regenerated1727 nodes and966 edges are byte-identical to the archived tables.
2. `summarize_small_cluster_coherence.py --membership selected_cluster_membership.tsv --pairwise pairwise_directed.tsv --fasta full_database_structural_sequences.fasta --outdir NEW_COHERENCE_DIR`. Default TM and coverage thresholds are0.50; the 63-row summary is byte-identical.
3. `summarize_direct_alignment_support.py --all-hits control_vs_baseline_all_hits.tsv --baseline-membership MEMBERSHIP.tsv --control-ids control_ids.txt --representative A0A3Q8Q3U2 --output-dir NEW_CONTROL_DIR`. The12-row coverage figure-input table is byte-identical. Requires pandas.
4. The preserved singleton core's `perform_statistics(master_dataframe,500)` was invoked on the saved QC master table in an isolated process. All30 rows agree numerically;24 pairwise tests are BH-adjusted together; seed 20260712. Requires pandas, numpy, scipy, matplotlib for module import. Structure parsing was not rerun.

The dataset `original_location_manifest.tsv` resolves every original input name to its deduplicated archive path. `restore_original_layout.py` can recreate aliases into a separate working directory; it must not modify the archive. Parameter details are in run_parameter_evidence.tsv. Full command receipts and comparison results remain in the local verification directory.

## Historical expensive workflows, preserved but not executed here

Baseline/spike clustering: `foldseek createdb`, `foldseek cluster --min-seq-id 0.0 -c 0.05 --cov-mode 0 --threads 64`, `foldseek align`, `foldseek convertalis`. The original baseline filename ends `.py` but contains Bash. The missing original `createtsv` invocation must not be invented as recovered history.

Direct control searches: `--exhaustive-search 1 --alignment-type 2 --max-seqs 100000 -e 1000000 -c 0.0 --cov-mode 0`. Six nonredundant control IDs are selected; the seven-entry spike benchmark is a separate population.

Small-cluster alignments: original full-database exhaustive runner and original retained pairwise outputs are archived. No new alignment was run.

Representative annotation: original API client, two original selected batch outputs/ledgers and the separate two-query follow-up are retained. The two original outputs require filtering to the 114non-giant representatives;29359records belong to112matched representatives. The501-row follow-up is separate. The final merge/classification/manual-integration source is missing; the general batch analyzer is not silently substituted for it.

Baseline provenance: the original 447-residue MGAT2 Foldseek search supplied the 1,727 structures used downstream. The separate truncated-core search is outside current release requirements. Its unavailable files are non-blocking; no further search/rerun or independently verified accession-equivalence claim is planned.
