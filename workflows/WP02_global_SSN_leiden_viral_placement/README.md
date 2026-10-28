# Cross-domain GT-A structure network and viral placement

This workflow supports Figure 2A/B of the associated viral glycosyltransferase manuscript. It constructs a structure-similarity network across 4,854 curated GT-A-like proteins, detects weighted Leiden communities, and evaluates viral candidates using CAZy family annotations, direct structural neighbours and network topology.

**Workflow:** Foldseek all-versus-all structural comparison → weighted structure-similarity network → Leiden communities → viral structural placement → Figure 2 source tables.

## Installation

Use Linux or WSL with Python 3.11. From this workflow directory:

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r environment/requirements.lock.txt
python -m pip check
```

The tested environment includes NumPy, pandas, SciPy, scikit-learn, python-igraph, leidenalg, Matplotlib and seaborn. See [environment notes](environment/README.md). Foldseek is needed only to recompute structural alignments; supplied alignment results suffice for all commands below.

## Input data

Obtain the companion `section2_ssn_005_data` dataset from **Zenodo: [ZENODO_DOI — insert after approved deposition]** and verify the release SHA256. It contains the 4,854 original coordinate files, raw and corrected Foldseek alignments, node metadata, full graph, primary communities, viral classifications and Figure 2A/B tables. Set `DATA` to the extracted dataset directory:

```bash
DATA=/absolute/path/to/section2_ssn_005_data
python scripts/verify_checksums.py --root "$DATA"
python scripts/verify_release.py --data-root "$DATA" --report ../reference_validation.json
```

To use the original structures, extract and verify them on **case-sensitive Linux storage**; 27 filename pairs differ only by case:

```bash
python scripts/verify_structures.py --data-root "$DATA" \
  --extract-dir "$HOME/GT_A_4854_input" --report ../structure_validation.json
```

## Execution

Reproduce the manuscript classifications and Figure 2 source counts using the preserved primary partition:

```bash
python scripts/reproduce.py --data-root "$DATA" \
  --outdir ../reference_reproduction --mode frozen --render-figure2b
```

Run Leiden and classification from the supplied alignments at **resolution 1.0, seed 1**, using the mean of directional `min(qTM,tTM)` values as the undirected weight:

```bash
python scripts/reproduce.py --data-root "$DATA" \
  --outdir ../fresh_analysis --mode primary
```

Use new output directories. The runner never invokes Foldseek or overwrites the dataset. A fresh run may differ from the manuscript partition in another software environment; exit code 2 reports a completed run with a reference mismatch, while individual step logs identify execution failures. The preserved primary analysis remains the manuscript reference. See [direct script commands and interfaces](EXECUTION.md) and [technical limitations](KNOWN_LIMITATIONS.md).

## Outputs

- `leiden/edges_clean_undirected.csv`: full weighted network; `cytoscape_nodes_analyzed.csv` and `nodes_all_leiden_assignments.csv`: all nodes, including isolates, with community assignments.
- `leiden/leiden_resolution_summary.csv` and `community_morphology_global_min_tm.csv`: network and community statistics.
- `placement/viral_structural_placement.tsv`: protein-level viral evidence and placement; `community_context_summary.tsv`: community composition; additional summary tables describe placement and connectors.
- Reference reproduction: `figure2a/` contains the computational display-edge subset; optional `figure2b/` contains alluvial source counts and rendered files. `reproduction_result.json` records commands, versions and reference comparisons.

There are 130 directly family-supported CAZy-unassigned proteins in the primary analysis; 110 support the modal community family and 20 another family. These are distinct evidence categories, not interchangeable totals. Structural support does not establish enzymatic activity or constitute formal CAZy reassignment.

## Citation and availability

- Manuscript: **[APPROVED_MANUSCRIPT_CITATION / DOI]**.
- Code: [viral-gt-research-record](https://github.com/LanNa1229/viral-gt-research-record), workflow `WP02_global_SSN_leiden_viral_placement`, **[APPROVED_RELEASE_TAG]**.
- Data: **[ZENODO_DOI]**. Replace placeholders when the corresponding releases are approved and identifiers exist.

Original code and verification helpers: MIT. Original data: CC BY 4.0; third-party structures retain their applicable terms. See [data license](DATA_LICENSE.md) and [coordinate rights](RIGHTS_AND_RETRIEVAL.md). This local candidate has not been published.
