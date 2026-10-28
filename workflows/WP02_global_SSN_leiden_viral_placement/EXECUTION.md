# Execution and input/output interfaces

Run commands from the workflow directory after activating the documented environment. Set `DATA` to the extracted companion dataset and choose a new `OUT` directory outside it. The recommended `scripts/reproduce.py` runner resolves paths and passes the recorded configuration explicitly.

## Weighted Leiden

```bash
DATA=/absolute/path/to/section2_ssn_005_data
OUT=../new_leiden_run
python original/analyze_ssn_leiden.py \
  --nodes "$DATA/inputs/metadata/005_GTA_node_metadata_completed_20260727.tsv" \
  --alignments "$DATA/inputs/foldseek/005_GTA_all_vs_all_TM065_cov050_rich.id_corrected.m8" \
  --alignment-columns query,target,fident,alnlen,evalue,bits,alntmscore,qtmscore,ttmscore,qcov,tcov,qlen,tlen,qstart,qend,tstart,tend \
  --outdir "$OUT" --resolutions 1.0 --seeds 1 --primary-resolution 1.0 \
  --n-iterations -1 --skip-plots
```

The supplied corrected rich M8 is headerless with 17 columns; the explicit schema is required. The raw M8 has a header. Node metadata are TSV with unique, case-sensitive `name` IDs and `Group`, `Species`, `GTFam` fields. The script reads metadata as strings, preserves all metadata nodes as vertices, removes self hits and unknown endpoints, and averages available directional row minima for each unordered pair. It does not require reciprocal pairs or impose another coverage/TM cutoff in this configuration. Upstream Foldseek used TM threshold 0.65 and mutual coverage 0.50; these do not require every global-minimum edge weight to exceed 0.65.

The original analyzer also computes unweighted and aligned-core TM comparisons. The preferred exported `community` column is the global-minimum partition when those weights are present. Leiden uses RBConfiguration, resolution 1, seed 1 and convergence iterations (-1). No CAZy annotation is used as an optimization label. Without explicit `--seeds`/`--resolutions`, original defaults run a broader comparison; use the command above for the primary analysis.

Outputs include `edges_clean_undirected.csv`, `cytoscape_nodes_analyzed.csv`, `nodes_all_leiden_assignments.csv`, per-method community/viral tables, `leiden_resolution_summary.csv`, `run_config.json` and `RUN_SUMMARY.md`. `cytoscape_edges_analyzed.csv` duplicates the full edge export; it is not needed as another dataset copy.

## Viral structural placement v3.0.0

```bash
python original/classify_viral_structural_placement_v3.py \
  --nodes "$OUT/cytoscape_nodes_analyzed.csv" \
  --edges "$OUT/edges_clean_undirected.csv" \
  --outdir ../new_placement_run \
  --node-id-column name --community-column community --virus-column is_virus \
  --cazy-family-column CAZy_family --gtfam-column GTFam \
  --weight-column weight_global_min_tm
```

For reference replay, replace the node and edge paths with `$DATA/outputs/primary/leiden/` files. The explicit `--community-column community` is required for these exports; the original classifier's default refers to a different assignment-column spelling. Current default classifier thresholds equal the primary recorded thresholds. The runner nevertheless supplies every threshold from `configs/historical/primary_classification_run.json` explicitly. Historical absolute paths in those archival configs are not used as public runtime paths.

The classifier consumes the full undirected graph, not a display subset. It calculates community composition, direct within-community family support, strength and topology, then emits all-node and viral-only TSV/CSV outputs plus summaries and `classification_run.json`. Protein-level output is `viral_structural_placement.tsv`; the all-node table is `cytoscape_nodes_with_figure2b_classes.tsv`. Neighbour support requires at least 2 annotated neighbours of the same family, weighted annotated-neighbour fraction ≥0.67 and best-minus-second fraction ≥0.20. Concordance compares that family with the modal community family. Composition/topology thresholds and interpretation limits are in Supplementary Methods.

## Figure 2B source tables

```bash
python original/make_figure2b_alluvial.py \
  ../new_placement_run/viral_structural_placement.tsv --output-dir ../new_figure2b
```

This manuscript-specific renderer requires 687 viral rows and produces `figure2b_alluvial_flow_counts.tsv`, `figure2b_node_counts.tsv`, SVG, PDF and PNG. Use the preserved partition to reproduce manuscript counts. A fresh partition's output is a separate result, even when the renderer accepts it. Missing Nimbus Sans causes font fallback only.

## Functional validation

After producing the two runner output directories (with `--render-figure2b` for the preserved run):

```bash
python scripts/validate_execution.py --data-root "$DATA" \
  --preserved-run ../reference_reproduction --fresh-run ../fresh_analysis \
  --report ../execution_validation.json
```

Checks cover CLI/syntax, pinned dependencies, directional aggregation and filters on a small independent fixture, isolates, weighted Leiden, fresh graph identity, all protein-level reference columns, fresh classification consistency and Figure 2B counts. Expected reference counts are assertions, never values assigned to scientific outputs. Validation does not run Foldseek or repeat a parameter sweep.
