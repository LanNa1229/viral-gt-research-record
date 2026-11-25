# Visual QC: minimalist supplementary UMAP v2

The source, layout and vector-export checks passed. Final visual-review hashes and any material limitations are recorded in VISUAL_QA.json.

## Direct comparison with Flamholz Figure 3a

The full publisher image was successfully viewed in the browser and saved for local comparison. The panel 3a crop was inspected beside the revised render. Both use a dominant landscape scatter, complete rectangular numerical axes, small functional-category points, a white background and a vertical right-hand legend. The new figure uses this dataset's unchanged saved coordinates, original class vocabulary and existing palette.

The redesigned plot occupies 75% of figure width, with a compact legend instead of a numerical annotation column. Our axes retain equal physical units, preserving the existing point-cloud geometry; its exact aspect and data shape differ from the reference. The prior figure remains intact in its original directory.

## Final-size assessment

- Full X/Y frame and numerical labels are visible. UMAP 1 and UMAP 2 titles are clear; no grid is present.
- Arial 9 pt axis titles, 8 pt tick labels and 8.5 pt legend text are readable in the fixed 183 x 110 mm layout.
- The legend uses seven filled class circles, one hollow hybrid circle and one hollow viral triangle. Long labels wrap within the reserved margin.
- No GC labels, leader lines, community table, percentages, n/N values or explanatory bottom notes remain in the artwork.
- Single-label reference color groups are immediately visible. Hollow hybrid circles reduce the former pie-sector complexity while preserving the location of all 467 multilabel references.
- All 287 selected viral loci are hollow triangles, preserving visibility of reference colors beneath them. Dense exact/near overlaps remain; no samples were moved or jittered.
- Normal, protanopia, deuteranopia, tritanopia and grayscale previews were inspected. Sample types remain distinguishable by shape/fill. Color discrimination among all seven classes is reduced in severe deficiency simulations and grayscale; hybrid memberships must be read from source data.
- The PDF is independently rendered at 300 dpi for final review; 150-dpi and 96-dpi previews retain the same physical layout. Screen previews do not constitute a calibrated print.

## Integrity and interpretation

All 2,502 reference BGCs and 287 selected viral loci are retained, with unchanged saved UMAP coordinates, original reference labels and ten original-space GC assignments. Hybrid is a neutral rendering choice, not a new biological class. All community statistics remain in the separate caption and source tables.

The layout communicates the broad biosynthetic organization of the learned representation more directly than the earlier annotated version. This is descriptive correspondence within a supervised representation; it does not independently validate viral functions. UMAP was not used as a clustering algorithm, and no scientific model was rerun.
