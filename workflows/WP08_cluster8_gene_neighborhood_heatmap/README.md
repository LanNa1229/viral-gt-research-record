# Main Figure 4 — Cluster 8 genomic-neighborhood heatmap

This historical research archive preserves the computational record used to compare the genomic context of 11 viral Cluster 8 glycosyltransferase loci across 33 functional gene classes. The figure summarizes the strongest named structural evidence per locus/class; it is not a count of genes or proof of biochemical activity.

**Lan Na — research conception, methodology, computational analysis, data curation, and visualization.**

All 346 files from the validated October 29 recovery are published individually and unchanged. Additional original query models, saved raw API responses, per-query tables and run records are preserved under `originals/source_runs/`, traced to the exact run-level hit tables used by the heatmap. No Foldseek searches or scientific analyses were rerun for publication. Proposed corrected outputs remain local and are not included in this historical-only version.

## Four-stage computational workflow

1. **Structural search:** `foldseek_batch_api_v2.py` reads query structure manifests, submits structures to the Foldseek API and saves response JSON, run ledgers and hit tables. Original script variants, query models and responses are archived; current online databases cannot recreate historical server state exactly.
2. **Annotation summaries:** `analyze_foldseek_results_v2.py` processes saved `foldseek_all_hits.tsv`, derives readable names and keyword tags, and writes annotation/high-confidence/ORF summary tables including `10_database_top_hit_names.tsv`.
3. **Neighborhood matrix:** `build_cluster8_neighborhood_matrix.py` reads `cluster8_analysis_manifest.fixed.tsv` and the 11 top-hit tables, applies the original 0–3 evidence and class rules, then writes the 11 × 33 matrix and long evidence table. The original focal-ID mismatches and row classifications are retained.
4. **Plotting:** `plot_cluster8_panelB_heatmap.py` and `plot_cluster8_panelB_heatmap_clean.py` render the original and clean heatmaps. No biological curation or figure redesign was applied to this archive.

## Inputs, code and results

- [Original fixed manifest](originals/260630_Temp/cluster8_analysis_manifest.fixed.tsv).
- [Matrix builder](originals/260630_Temp/build_cluster8_neighborhood_matrix.py), [original plotter](originals/260630_Temp/plot_cluster8_panelB_heatmap.py), [clean plotter](originals/260630_Temp/plot_cluster8_panelB_heatmap_clean.py).
- [Batch and annotation script variants / manifests](originals/related_scripts/) and [source query models / complete saved runs](originals/source_runs/).
- `originals/260630_Temp/<candidate_id>/final/`: all 11 original `foldseek_all_hits.tsv` files and `analysis_v2/` tables, including all 11 `10_database_top_hit_names.tsv` files.
- [Figure source matrix](originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_matrix.tsv) and [4,919-record evidence table](originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_evidence_long.tsv).
- [Original and clean PDF/PNG figures](originals/260630_Temp/cluster8_global_neighborhood_panelB/).
- [Main Figure 4 source mapping](docs/FIGURE_SOURCE_MAP.tsv), [Methods and confidence scoring](docs/METHODS_AND_SCORING.md), [all 33 gene-class definitions](docs/GENE_CLASS_ASSIGNMENTS.md).

## Dependencies and reproduction

Offline reconstruction from the preserved top-hit tables requires Python, pandas, NumPy, Matplotlib and Pillow; see [environment.yml](environment.yml). The tested matching environment used Python 3.12.11, pandas 2.3.3, NumPy 2.4.6, Matplotlib 3.10.8, Pillow 12.0.0 and FreeType 2.6.1. The API script additionally requires `requests`, but it is not invoked by offline reproduction.

From this workflow directory, choose a new output directory:

```bash
python -B reproduce.py --out /absolute/path/to/new_validation
```

For the exact rendering check in the matching Linux environment:

```bash
python -B reproduce.py --out /absolute/path/to/new_validation --require-exact-render
```

To plot the immutable historical matrix directly:

```bash
python -B originals/260630_Temp/plot_cluster8_panelB_heatmap_clean.py \
  --matrix originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_matrix.tsv \
  --outdir /absolute/path/to/new_figure_output
```

Expected offline validation outputs include rebuilt matrix/evidence TSVs, validation JSON, and the clean heatmap PDF/PNG in `plot_with_title/`. The wrapper resolves packaged input tables without modifying preserved absolute paths or focal accession IDs. [Detailed reproducibility instructions](docs/REPRODUCIBILITY.md) also describe the original builder CLI. The original API stage is preserved for provenance, not automatically rerun.

## Completed validation

- 11 loci, 33 classes, **363/363 score cells reproduced**.
- **4,919/4,919 evidence records reproduced**; evidence row order is compared as a multiset because the original code iterates sets.
- Original clean **PNG reproduced byte-for-byte**; the reproduced PDF differs only in creation timestamp.
- Original-file hashes and source mtimes were checked during recovery; publication verifies every archived byte against its recorded checksum and its Git blob.
- [Manuscript comparison](validation/manuscript_comparison/AUDIT_REPORT.md): **363/363 cells match** in actual v1.4 and final available v1.5; zero score discrepancies or uncertain cells. Their embedded PNGs are byte-identical. Word crops the bottom 14.262%, hiding the group legend but not the score grid; this is a presentation difference, not a scientific correction.

## Limitations, provenance and rights

Retain the [known historical limitations](docs/KNOWN_LIMITATIONS.md), including A0A6J5M498 row-group fall-through, broad nucleotide-metabolism categories, ambiguous ORF assignments, and five manifest-query identity mismatches. No IDs, scores, annotations, classifications or original figures have been corrected. The earlier recovery reports may describe local proposed corrections; those outputs are outside this historical publication and any future curated version must be separate.

[Recovery report](docs/RECOVERY_REPORT.md), [provenance](docs/PROVENANCE.md), [publication scope and source supplement](docs/PUBLICATION_SCOPE.md), and source checksums document the record. `SHA256SUMS-originals.txt` covers the original 346; `SHA256SUMS-source-supplement.txt` covers additional recovered source files; `SHA256SUMS-package.txt` covers the exact historical publication payload except itself. Historical local paths are execution provenance, not portable defaults.

The repository's existing [MIT notice](licenses/repository_existing/LICENSE-MIT.txt), [CC BY notice](licenses/repository_existing/LICENSE-CC-BY-4.0.md) and [third-party attribution document](licenses/repository_existing/THIRD_PARTY_NOTICES.md) are retained verbatim. Their original Section 1 scope is not expanded to WP08 by this archive. [WP08 rights and attribution](THIRD_PARTY_NOTICES.md) distinguishes existing grants from this authorized historical publication. No Zenodo deposition or release is part of this task.
