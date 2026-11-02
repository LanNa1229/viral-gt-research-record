# Targeted manuscript heatmap audit — November 2, 2026

**Conclusion: the final manuscript heatmap faithfully represents the historical TSV. All 363 score cells agree, with no confirmed score discrepancies and no uncertain comparisons. No manual scientific score corrections or row-group corrections are evident.**

The specified v1.4 manuscript and the latest available final v1.5 manuscript contain the same 1646 × 779 PNG, byte-for-byte, with SHA-256 `8a938331cff92bd986ec142b63c18cf316ad6557cf6f2172701ee9178764e5bb`. Both store it as `word/media/image19.png`. The source DOCX paths, file hashes, image inventories and drawing XML are preserved in `document_inventory.json`.

## Scope and sources

- Historical `cluster8_neighborhood_matrix.tsv` and clean PNG from the October 29 byte-preserved recovery of `/home/rahnn/260630_Temp/`; their source provenance and prior preservation check are already established. This audit uses those copies, verifies their identities and does not repeat the completed reproducibility run or execute the original analysis.
- Actual v1.4: `D:/Projects/Viral_GT_Section1_Structural_Clustering_WorkFlow/Viral_GT_draftv1.4_261022_LanNa_references_added.docx`.
- Latest available final v1.5: `C:/Users/rahnn/Downloads/Viral_GT_draftv1.5_261027_LanNa_embedding_supple_final.docx` (local modified time October 29, 2026 12:08:27 Eastern).
- Two additional local drafts were inventoried: v1.3 `...references_added_2.docx` and v1.5 `...SuppleFig_embedding_analysis.docx`. Their media inventories/contact sheets do not contain this target PNG/panel. These are not presumed to be the v1.1/v1.3 versions independently examined by the user. The user's earlier-version identity finding is supplied context, not a newly verified finding here.

## Score comparison

The historical figure is 4144 × 1958 pixels; the manuscript image is 1646 × 779. Each grid was located at its own native resolution. The 33 column labels and 11 row labels/order were checked against the historical plot. The four score colors were independently calibrated from each image's own colorbar and match the original palette exactly.

All cells were sampled in five offset, border-excluding patches. A score was accepted only when every patch assigned at least 99% of pixels to the same palette category and had median RGB distance ≤3. All 363 manuscript cells and all 363 historical-PNG cells passed, and both agree with every corresponding TSV score. Gridline pixels and text were excluded; no ambiguous edge pixel was treated as scientific evidence.

`all_363_score_cell_comparisons.csv` supplies the required locus, functional class, original TSV score and manuscript figure score for every cell, plus sampling evidence. `confirmed_discrepancies.csv` and `uncertain_comparisons.csv` contain headers only because both sets are empty. `image_sampling_details.json` records all patch measurements and calibration. `score_grid_comparison.png` shows the complete TSV/historical-PNG/manuscript decoded grids together.

## Separate grouping, label and cosmetic checks

| Feature | Finding |
|---|---|
| Score cells | 363/363 agree; 0 discrepancies; 0 uncertain |
| Row order and group colors | All 11 rows agree with historical TSV grouping; no regrouping |
| Locus labels | Same 11 labels, same order; no confirmed label edits |
| Class labels, module headers, title and score legend | Same visible text/categories; no confirmed relabeling |
| Embedded raster resolution | Reduced from 4144 × 1958 to 1646 × 779 |
| Word display crop | Both v1.4 and final v1.5 specify `<a:srcRect b="14262"/>`: bottom 14.262% cropped |
| Effect of crop | Hides the bottom neighborhood-class legend; does not crop any score cell, row label or functional-class label |
| Image scaling and antialiasing | Small pixel-level differences around edges/text are compatible with resampling; they do not change decoded values |

The row-group strip is overlapped by locus-label text in the historical layout. Group colors were therefore checked using the modal solid strip fill over a row-interior region, rather than mistaking black label ink for a changed group color. See `row_group_comparisons.csv`.

`word_display_transforms.json` records exact crop and extent settings. Cropping is a document-layout modification; it does not replace the underlying embedded PNG. Native stored images and clearly labeled crop previews are retained separately.

## Was this merely resized or scientifically corrected?

The image is consistent with a downsampled copy of the original plot. All 363 score cells, row-group categories and labels are preserved, including the historically problematic A0A6J5M498 glycan-associated grouping despite TET score 3. A0A3G1L3B7 likewise retains its historical nucleotide-metabolism grouping. The later proposed TET grouping correction from WP08 recovery is not present in either audited manuscript.

Six common resampling filters and several nearest-neighbor coordinate conventions were tested without altering source images. No exact pixel-for-pixel reproduction of the smaller image was obtained; the best whole-image comparison has about 96.2% identical pixels, with small channel error. This establishes neither a unique resize algorithm nor an editing history. Therefore this report does **not** claim forensic proof that no individual cosmetic pixel was ever edited. It does establish that there is no detectable scientific score, group or label correction in the manuscript panel. The confirmed modification is the Word display crop and reduction in resolution.

No changed-score list exists to investigate as manual annotation corrections. Any remembered manual annotation review may have confirmed existing scores, concerned a different figure/version, or left the final panel unchanged; the image alone cannot distinguish those possibilities.

## Preservation and publication status

The original TSV and original clean PNG are copied byte-for-byte under `historical/`. Both manuscript embedded PNGs are preserved under their document-specific paths under `manuscript_images/`. `manuscript_figure_decoded_matrix.tsv` is an independently decoded figure-version table, explicitly separate from the original scientific TSV. DOCX bytes and modification times were checked before/after read-only extraction and again at the end of the audit.

All new records reside in this separate audit directory. Existing WP08 package files and the local Git repository were not edited, staged, committed or published. No Foldseek search, annotation recomputation or original plotting/reproducibility audit was rerun. All scripts used for this targeted comparison and SHA-256 checksums are included.
