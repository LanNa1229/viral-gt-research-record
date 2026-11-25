# Minimalist supplementary UMAP, version 2

A new presentation-only revision with full numerical Cartesian axes, larger Arial text, a dominant scatter panel and a compact vertical legend. The previous figure package is preserved.

## Revised figure

- [Editable SVG](figures/SuppFig_global_leiden_biosynthetic_classes_minimalist_v2.svg)
- [Vector PDF](figures/SuppFig_global_leiden_biosynthetic_classes_minimalist_v2.pdf)
- [600-dpi PNG](figures/SuppFig_global_leiden_biosynthetic_classes_minimalist_v2.png)
- [150-dpi preview](qc/preview_183mm_150dpi.png)
- [96-dpi preview](qc/preview_183mm_96dpi.png)
- [Independent PDF render](qc/PDF_render_300dpi.png)

Canvas 183 x 110 mm; plot width 75%. Arial axis titles 9 pt, ticks 8 pt, legend 8.5 pt. All 2,502 references and 287 selected viral loci remain. The 467 hybrids are neutral hollow circles; viral loci are hollow triangles. Original multilabel annotations and all GC memberships remain in the source data.

## Supporting records

- [Publication-style figure legend](SUPPLEMENTARY_FIGURE_LEGEND.md)
- [Figure source data](FIGURE_SOURCE_DATA.tsv)
- [Original community summaries, outside artwork](COMMUNITY_ANNOTATION_KEY.tsv)
- [All reference-class marginal counts](REFERENCE_CLASS_MARGINAL_COUNTS.tsv)
- [Final HEX palette](COLOR_PALETTE.tsv)
- [Data dictionary](DATA_DICTIONARY.md)
- [Methods record](METHODS_RECORD.md)
- [Visual QC report](VISUAL_QC_REPORT.md)
- [Machine-readable validation](VALIDATION_RESULTS.json)

prepare_figure_inputs.py verifies the prior sealed package and freezes sources. plot_integrated_umap.py renders only the saved data using config.json. validate_figure.py independently checks scientific fields, original coordinates/assignments, actual SVG glyph counts, axes, fonts and export structure. run_stage.py stores contemporaneous code/config snapshots, actual commands, software/platform, stdout/stderr and exit status under executions/.

INPUT_MANIFEST_SHA256.tsv lists 34 input records. PLOTTED_SAMPLE_RECEIPTS.tsv and RENDER_GROUP_RECEIPTS.tsv record every point and display group. LAYOUT_RECEIPTS.json, SOFTWARE_VERSIONS.json, PDF_RENDER_EXECUTION.json and VISUAL_QA.json retain layout and visual evidence. finalize_package.py seals OUTPUT_MANIFEST_SHA256.tsv and CHECKSUM_VERIFICATION.json after all execution files exist.

qc/ contains final-size, CVD and grayscale previews plus the actual publisher image/crop used for direct style comparison. The publisher image is a local reference asset and is not embedded in the scientific exports.

No UMAP/Leiden/model analysis was rerun. No previous figure was overwritten. Nothing was committed or pushed. Any later upload requires explicit scoped approval.
