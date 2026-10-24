# Restoring a separate reproducibility run

The frozen scientific record is in the current archive; no new scientific run is needed to verify it. Use `python reproduce/verify_archive.py` for archive-only size/SHA256 checking.

To run a separately dated reproduction, first create a dedicated Python environment with the **recorded** package versions in `runs/run_001_science/before.json`, then locate its site-packages directory containing igraph/leidenalg. The graph wheels used originally are external, recorded by exact filenames/SHA256 in `runtime/WHEEL_MANIFEST.tsv` and the successful pip report. Acquire those exact wheel bytes from PyPI when recreating the Windows runtime; Windows wheels cannot be reused as Linux binaries. Other distributions/native files are recorded in the execution fingerprints. A new platform/package build is not guaranteed bitwise identity and must retain its own records.

Locate the previous UMAP archive at `results/WP07_viral_gt_figure3C_umap_20261024/` in the repository, using commit `7b35093eefa99a4704e7b64039f51b7564788fd4` when an immutable source is needed. It supplies protected old analytical files, not new embeddings or inference.

Run the file-restoration helper against a nonexistent new workspace:

```bash
python reproduce/prepare_reproduction.py --workspace-root /absolute/new/workspace --previous-archive /absolute/repository/results/WP07_viral_gt_figure3C_umap_20261024 --graph-runtime /absolute/environment/site-packages
```

It restores the 14 archived scientific/provenance inputs plus protected previous-UMAP artifacts, verifies all 41 original input identities, copies current scripts to a new output directory, and updates only the new copied configuration's locations. It copies no old run logs or result tables into the new clustering directory. If analysis dependencies are not installed under the historical project `.deps` path, the chosen interpreter can provide them from its environment; `--analysis-dependencies` controls the optional prepended relative location. The helper does not install packages, fit models, or execute scientific scripts.

From that new workspace, using the intended environment's interpreter:

```bash
python -B analysis/viral_gt_clustering_20261024/scripts/run_recorded.py science analyze.py
python -B analysis/viral_gt_clustering_20261024/scripts/prepare_derived_inputs.py
python -B analysis/viral_gt_clustering_20261024/scripts/run_final_recorded.py figures figures_final.py
python -B analysis/viral_gt_clustering_20261024/scripts/run_final_recorded.py validation validate_extended.py
python -B analysis/viral_gt_clustering_20261024/scripts/record_documentation.py
python -B analysis/viral_gt_clustering_20261024/scripts/record_cross_grid.py
```

A restoration retains the historical script/directory name but the logger records the new actual UTC execution dates. Treat resulting logs, versions, hashes and results as a new reproduction, never as missing evidence for the original run. The scientific baseline and grid remain fixed. Do not edit the previous UMAP archive or original workspace. The old bootstrap download/install scripts retain historical machine paths and are preserved as records; install into the intended new environment instead of executing those old setup scripts blindly.

The current source `README.md` is preserved byte-for-byte, as are all completed scientific artifacts. `OUTPUT_MANIFEST.tsv` is the original source-analysis manifest. `ARCHIVE_SHA256.tsv` additionally covers archival input copies and these restoration utilities. Logs include unsuccessful attempts; a parent-wrapper integrity failure is distinguished from a successful child renderer. Archived native environments and wheels are external, not silently included or omitted without identities.
