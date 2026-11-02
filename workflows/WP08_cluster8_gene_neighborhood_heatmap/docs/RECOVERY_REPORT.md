# Recovery and validation report — October 29, 2026

The requested completed analysis is recovered, rather than reconstructed.

| Requested item | Recovery location under `originals/` |
|---|---|
| `build_cluster8_neighborhood_matrix.py` | `260630_Temp/` |
| `plot_cluster8_panelB_heatmap_clean.py` | `260630_Temp/` |
| `foldseek_batch_api_v2.py` | `related_scripts/`, multiple source directories |
| `analyze_foldseek_results_v2.py` | `related_scripts/`, multiple source directories |
| `cluster8_analysis_manifest.fixed.tsv` | `260630_Temp/` |
| `cluster8_neighborhood_matrix.tsv` | `260630_Temp/cluster8_global_neighborhood_panelB/` |
| `cluster8_neighborhood_evidence_long.tsv` | same directory |
| original clean PDF and PNG | same directory, `figures_clean/` |
| 11 `10_database_top_hit_names.tsv` inputs | `260630_Temp/<candidate_id>/final/analysis_v2/` |

The non-clean original plotting script and figures, earlier malformed manifest, raw combined hit tables and saved annotation reports were also preserved. Recovery searched the supplied Ubuntu directory, related immediate Foldseek batch directories and the current Windows workspace. It was a targeted local recovery, not an exhaustive disk-forensic search. All essential files were found, so neither code reconstruction nor Foldseek execution was necessary.

## Executed checks

- 11 unique manifest candidates, 33 gene classes, 363 score cells.
- Complete matrix equality including historical `analysis_dir`, row order, scores, grouping, inferred model and group rank.
- All 4,919 evidence records reproduce as a multiset: no missing or extra records. Original set iteration makes evidence row ordering nondeterministic.
- Independent maximum aggregation of the saved evidence exactly matches every gene-class cell.
- Original PNG and reproduced PNG have 4144 × 1958 pixels and identical bytes: SHA-256 `22a09d804340686bab353904471c384eee5a0fec881e504141a6c68c02ee64e7`.
- PDF comparison: identical bytes after removing only the CreationDate value for comparison. Original recorded date: June 30, 2026 at 12:29:46 UTC−04:00. The original PDF is preserved untouched.
- All 346 archived files match their source SHA-256 values; a final source check also confirms original source content and mtimes are unchanged.
- Five exact manifest-query mismatches were found and retained. See `validation/input_audit.json` and limitations.
- Original PNG visually inspected: all 11 labels, 33 columns, legend and score colorbar present. Existing layout crowding near headers and row labels is retained as historical figure content, not silently redesigned.

Historical group counts: 4 dCMP-hydroxymethylase-associated, 5 broad nucleotide-metabolism-associated, 2 additional-GT/glycan-associated. No locus is in the combined C5-MT + 5mYOX class. No locus has a nonzero ThyX score; the short legend label “ThyX/ThyS” therefore cannot be read as proof of ThyX.

Assessment: ready for historical reproducibility review with explicit biological caveats. Release metadata/rights remain a publication prerequisite. Repository source context was read at local HEAD `52dd6147505ce3434a688f7703ed76baec731c3d`; the remote HEAD was not fetched or changed. Existing unrelated working changes are recorded in `audit/repository_discovery.json`.
