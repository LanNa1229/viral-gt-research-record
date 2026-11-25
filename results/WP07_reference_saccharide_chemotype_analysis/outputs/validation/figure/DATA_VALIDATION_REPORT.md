# Data validation report

PASS. The complete machine-readable results are in VALIDATION_RESULTS.json.

| Check | Status | Evidence |
| --- | --- | --- |
| exact_195_unique_original_Saccharide_BGCs | PASS | "All and only original Saccharide members, including hybrids; no duplicate accessions" |
| all_91_original_fields_preserved | PASS | "Every copied field remains exact, including detailed chemistry, confidence, source URLs/hashes and review notes" |
| mutually_exclusive_exhaustive_populations | PASS | "100 GC08 + 24 GC05 + 71 other GCs = 195; GC05 is excluded from Other GCs" |
| original_communities_and_labels_unchanged | PASS | "Exact archived k20_r1_s42 communities and original broad class strings/memberships" |
| hybrids_and_unresolved_retained | PASS | "65 hybrid records, two unresolved and one provisional source assignment retained; unresolved is not a negative" |
| primary_chemotypes_not_reassigned | PASS | "Primary chemistry and confidence match frozen curation decisions; no new annotation" |
| exactly_one_display_group_per_BGC | PASS | "Explicit total mapping of 15 detailed source categories to seven mutually exclusive display groups; no fuzzy classification or missing-value imputation" |
| historical_expectations_independently_reproduced | PASS | {"GC08": [31, 19, 18, 4, 20, 4, 4], "GC05": [4, 2, 3, 0, 7, 5, 3], "Other GCs": [1, 1, 1, 0, 7, 51, 10]} |
| all_21_source_cells_and_100_percent_totals | PASS | "Accession-derived numerators, Saccharide-only denominators and exact rational cumulative widths checked; structural zeros retained" |
| surface_operational_definition_consistent | PASS | "Categories 1-4: 72/100, 9/24, 3/71. GC08 phosphonoglycans remain in Other/unresolved; broader source flags are preserved but not counted" |
| plotted_segment_values_match_source | PASS | "All 21 segment records, including two zeros, match exported source data in stable category order" |
| actual_SVG_segment_geometry | PASS | "Each SVG rectangle width independently matches its percentage; zero segments stay zero; smallest positive segment is 1/71 = 1.40845% (~2.04 mm)" |
| axis_labels_scale_and_order | PASS | "GC08, GC05, Other GCs top-to-bottom; x axis 0-100% at 25-point intervals" |
| typography_and_text_layout | PASS | "Arial 8-8.5 pt; no clipping/overlapping labels and no crowded in-bar text" |
| seven_category_legend_and_fixed_colors | PASS | "Seven solid, distinct HEX colors; four warm surface groups remain consecutive; UMAP Saccharide gold retained for EPS" |
| editable_vector_exports_and_dimensions | PASS | {"dimensions_mm": [182.9999999999947, 95.00000000000804], "fonts": [{"base_font": "/RDHMQP+ArialMT", "type": "/CIDFontType2", "TrueType_embedded": true}], "SVG_text_nodes": 16} |
| 600dpi_PNG | PASS | {"pixels": [4322, 2244], "dpi": [599.9988, 599.9988]} |
| all_frozen_inputs_unchanged | PASS | {"inputs_rehashed": 18, "original_master_SHA256": "5f86b03ff3e071eec4fcafa02eac36c981f1561777a1d1acbe6f69e51315cbc9"} |
| color_vision_and_grayscale_previews | PASS | "Three CVD and grayscale views generated and inspected; exact display simulation matrices and pairwise diagnostics preserved" |
| independent_PDF_render | PASS | "Poppler independently rendered the final vector PDF for final-size visual review" |


## Verified composition

| Chemotype group | GC08 (100) | GC05 (24) | Other GCs (71) |
| --- | --- | --- | --- |
| Capsule | 31 | 4 | 1 |
| EPS/matrix polysaccharide | 19 | 2 | 1 |
| LPS/LOS/O-antigen | 18 | 3 | 1 |
| S-layer/other surface glycan | 4 | 0 | 0 |
| Aminoglycoside/aminocyclitol | 20 | 7 | 7 |
| Glycosylated specialized metabolite | 4 | 5 | 51 |
| Other/unresolved | 4 | 3 | 10 |


No discrepancy from the supplied historical expectations was found. The three populations are disjoint and exhaustive, with Saccharide-only denominators. All 65 hybrid records and all unresolved/provisional source states are retained. All 91 source fields are unchanged at accession level. Every BGC has exactly one display category. All three stacks total exactly 100% under rational accounting.

The final figure uses existing annotations as proposed/source-based categories. This validation establishes provenance and display consistency, not experimental confirmation of those chemical assignments. No new annotation, statistical test or model was performed.
