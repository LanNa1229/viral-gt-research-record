# Section 1 structural mining and coarse clustering

Version section1-v1.0.0. Sole record creator: Lan Na, University of Georgia, https://orcid.org/0000-0003-2788-6123 . Original code: MIT. Original data/documentation: CC BY 4.0. Third-party terms remain in THIRD_PARTY_NOTICES.md.

This is the selected Section 1 snapshot supporting Figure 1 and Supplementary Figures 1–4, not a license or DOI assignment for the wider repository. The 1,727-structure baseline came from the original 447-residue MGAT2 Foldseek search. Missing truncated-core files are non-blocking historical provenance; no independently verified output equivalence is claimed.

## Find files by figure

Figure numbers follow manuscript v1.5. The [panel evidence map](../../reports/section1_archive/figure_evidence_map.tsv) documents verification and limitations for individual panels.

| Figure | Materials | Location |
| --- | --- | --- |
| Figure 1 | Search provenance and cluster membership | [Search inventory and checksums](../../results/WP01A_initial_search/) · [Canonical membership](../../results/WP01B_coarse_clustering/canonical/cluster_membership.tsv) |
| Supplementary Figure 1 | Network source, editable Cytoscape session and SVG | [Network tables](../../results/WP01B_coarse_clustering/figure_X-1-1/network_source/) · [Cytoscape and SVG](../../figures/WP01B_coarse_clustering/figure_X-1-1/) |
| Supplementary Figure 2 | Positive-control inputs and results | [Control records](../../results/WP01C_positive_controls/) · [Public release exclusions](../../PUBLIC_EXCLUSIONS.tsv) |
| Supplementary Figure 3 | Cluster sizes, annotations and structural coherence | [Clustering records](../../results/WP01B_coarse_clustering/) · [Representative annotation](../../results/WP01D_representative_annotation/) · [Coherence tables](../../results/WP01E_small_cluster_coherence/) |
| Supplementary Figure 4 | Singleton confidence, lengths, thresholds and statistics | [Singleton QC](../../results/WP01F_singleton_qc/) |

## Shared data, code and records

- **Figure source data and assets.** [Frozen plotting inputs and assets](../../results/section1_source_data/) · [Source Data archive on Zenodo](https://doi.org/10.5281/zenodo.23626546).
- **Original outputs.** [Raw clustering assignments and alignments](../../results/WP01B_coarse_clustering/raw/) · [Original results grouped by WP01A–WP01F](../../results/).
- **Code.** [Historical scripts grouped by work package](../../scripts/archive/) · [Section 1 reconstruction and validation scripts](../../scripts/reproducible/WP01_section1/).
- **Reproducibility records.** [Snapshot manifest](snapshot_file_manifest.tsv) · [SHA-256 checksums](snapshot_SHA256SUMS) · [Execution record](../../reports/section1_archive/execution_commands.md) · [Release validation report](../../reports/section1_archive/release_validation_report.md).

## Validation and reconstruction

Run from the repository root (or the root of the selected software ZIP):

```bash
python -B scripts/reproducible/WP01_section1/validate_section1.py --source-data results/section1_source_data/figure_inputs.json --out /path/to/new_validation/result.json --bootstrap
MPLCONFIGDIR=/path/to/new_plot_cache python -B scripts/reproducible/WP01_section1/reconstruct_figures.py --source-data results/section1_source_data/figure_inputs.json --outdir /path/to/new_figures
```

The default numerical validator needs only Python; --bootstrap requires NumPy. Reconstructed figures require NumPy and Matplotlib. The frozen-input run passed 106 checks, including 24 Mann–Whitney/BH comparisons and three displayed 500-resample intervals (seed 20260712). The prior Source Data verification covered 18 sheets and 79,674 data cells without mismatches.

Historical scientific scripts remain byte-identical in scripts/archive/. ORIGINAL execution-supported sources, ARCHIVED surviving sources, RECONSTRUCTED plotting code and AUDIT-ONLY verification code are distinguished in reports/section1_archive/. Historical absolute defaults in original code are retained as necessary provenance; use isolated working directories and explicit inputs/environment variables. Account-specific prefixes are redacted from public metadata, with changed-file identities recorded in REDACTED_METADATA.tsv.

The public dataset withholds institutional AF3 model/source files and structural databases containing them because redistribution authority was not confirmed. The archived model plate is also withheld conservatively; new plotting code shows a disclosure placeholder for panel B while reproducing numerical A/C. Exact exclusions and hashes are in PUBLIC_EXCLUSIONS.tsv. The complete private rc2 archive remains unchanged. Final classifier source, L230 copy lineage and historical-runtime limitations remain disclosed; no new investigation is required for the nonessential truncated-core issue.

Software DOI: https://doi.org/10.5281/zenodo.23626496 . Dataset DOI: https://doi.org/10.5281/zenodo.23626546 . GitHub release: https://github.com/LanNa1229/viral-gt-research-record/releases/tag/section1-v1.0.0 . Record pages establish actual publication status.
