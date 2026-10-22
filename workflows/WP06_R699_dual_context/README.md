# R699 dual-context structural-core mapping

Supplementary Figure X-2-2b: independent C005/C009 reference-core projection onto R699 (9DYT chain A).

## Archive status: incomplete inputs and outputs

The three original analysis/visualization scripts and original Cytoscape annotation table are archived here unchanged. The input structures, reference metrics, alignments, and residue-level analysis outputs have not yet been recovered. This package is NOT yet a complete reproducible release. No original output tables have been reconstructed from reported summary counts.

`reference_members.tsv` is derived from the original Cytoscape annotation CSV: 13 C005 and 15 C009 references. The CSV contains 76 unique nodes: 13 C005 references, 15 C009 references, one R699 query and 47 other neighbors. These categories denote analysis roles, not CAZy assignments for all nodes.

## Essential files

- `map_dual_context_cores_to_r699.py`: freeze reference-core membership from existing metrics and transfer it through medoid residue indices to query-inclusive alignments.
- `RUN_R699_DUAL_CONTEXT_MAPPING.sh`: generate the two query-inclusive FoldMason alignments and run the mapper.
- `make_r699_dual_context_pml_fixed.py`: generate PyMOL selections using the mapped chain, residue number and insertion code.
- `R699_C005_C009_Cytoscape_node_annotation.csv`: original Cytoscape node styling table; join `name` to the network's matching identifier column.

## Required original run directory

Recover the original R699_WF1 directory, including:

- `C005/` (13 coordinate files) and `C009/` (15 coordinate files)
- `query_structure/9DYT_A.pdb`
- `work/00_inputs/pairwise_structural_scores.tsv`
- `work/01_medoids/medoids.tsv`
- `work/02_alignments/alignments/within/C005/C005_aa.fa`
- `work/02_alignments/alignments/within/C009/C009_aa.fa`
- `work/03_metrics/all_residue_metrics.tsv`
- `work/05_query_mapping/` (C005_R699 and C009_R699 alignments)
- `work/06_R699_dual_context_mapping/` (mapping tables, classification, summary and visualization)
- reference-metric generation configuration/commands, original workflow source and logs

The mapper consumes precomputed reference metrics; it does not itself reproduce their calculation. Recover those source commands and quality-filter definitions before declaring end-to-end reproducibility.

## Execution after original inputs are restored

Requires Python 3.10+, pandas, Biopython, FoldMason and PyMOL for visualization. Original Python dependency versions are not known. The analysis notes report FoldMason commit `9a71f1dc2eaa5f6ac5ec79dde0bc614a9dfc2057` and Foldseek commit `755d3fa84b8259d4d30aa92782f0bfc2bfc3363f`; these await confirmation from the original logs.

From this package directory, set RUN_DIR to the recovered run directory and run:

```bash
python -m pip install -r requirements.txt
bash RUN_R699_DUAL_CONTEXT_MAPPING.sh "$RUN_DIR" "$RUN_DIR/query_structure/9DYT_A.pdb" "$PWD/map_dual_context_cores_to_r699.py"
python make_r699_dual_context_pml_fixed.py \
  --input "$RUN_DIR/work/06_R699_dual_context_mapping/R699_dual_context_classification.tsv" \
  --structure "$RUN_DIR/query_structure/9DYT_A.pdb" \
  --output "$RUN_DIR/work/06_R699_dual_context_mapping/view_R699_dual_context.pml"
```

Always supply the structure argument: the unchanged historical shell script defaults to an author-local absolute path. The shell script replaces query-mapping intermediates in the supplied run directory; use a copy of the archived run for reruns. The generated PML references a local structure path and should be regenerated on the machine used for visualization.

## Interpretation and validation

The runner applies `passes_within_quality`, within_lddt >= 0.75 and within_presence >= 0.70. R699 is excluded from reference-core definition but included in the separate alignments used for projection. Supported means alignment correspondence to a reference-core position; it is not a direct measurement of query local geometric conservation, catalytic function, or evolutionary intermediacy.

Historical notes report 145 both-supported, 48 C005-only, 14 C009-only and 235 unsupported positions (442 total). These counts are validation targets only, not independently reproduced results in this package. Confirm that all classified sequence positions map to modeled coordinates before describing the denominator as modeled residues. Domain-specific enrichment and the source of network edges require separate verification.

Completed checks: Python syntax, shell syntax, unique node IDs, reference-set counts. Scientific rerun remains blocked by missing original inputs. SHA256SUMS records the files in this archive.
