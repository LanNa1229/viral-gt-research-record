# Methods record: minimalist UMAP redesign

## Presentation-only revision

This new sibling directory preserves the previous figure package at analysis/viral_gt_suppfig_global_leiden_20261027/. All 81 files in that package's sealed output manifest were checksum verified before and after this task, including the original PDF/SVG/PNG, plotting code, source data and QC records. No previous file was overwritten. No scientific model, UMAP, Leiden clustering or embeddings were rerun. No Git commit or push was performed.

The input manifest records 34 sources: 21 prior frozen scientific/font inputs plus the prior figure/source/palette/QC/configuration/manifests and the observed publisher reference image. Their SHA-256 values were recorded before rendering. The updated TSV preserves all previous scientific fields exactly. Only rendering metadata change: plot_marker and new rendering_group, plot_edge_or_fill_hex and plot_facecolor fields. All stable IDs, original coordinate strings, class annotations, hybrid memberships and baseline GC IDs remain unchanged.

## Data and lineage

Original UMAP analysis: viral_gt_umap_20261024. The fitted projection contains 3,058 samples: 2,502 references and 556 viral loci. The current visualization uses the same saved coordinates for all 2,502 references and exactly the 287 selected viral locus IDs, giving 2,789 displayed glyphs. The 269 unselected viral loci remain omitted. Filtering does not recompute or transform coordinates.

Global Leiden assignments are read from the existing original-320D baseline column k20_r1_s42 (k=20, resolution=1.0, seed=42). Private archival commit: 284d4f1655cf036c1890bdffcb4aa2dc93f6e41e in LanNa1229/viral-gt-internal-analyses/results/WP07_viral_gt_figure3C_clustering_20261024/. All assignments were independently compared with the archived grid and, for selected viral loci, the integrated viral master. UMAP geometry does not define these communities.

The original seven broad reference labels remain Alkaloid, Terpene, NRP, Polyketide, RiPP, Saccharide and Other. Exactly 2,035 references have a single label and use that label's filled-circle color. All 467 multilabel references use the same neutral hollow-circle rendering. Their original class strings and JSON lists remain unchanged. Hybrid is a rendering category, not a new biological class or a substituted single annotation. The community key and all 70 marginal reference-class composition rows were copied unchanged and independently cross-checked against the original annotation table. They are separate supporting data, absent from the artwork.

## Graphical specifications

Canvas: 183 x 110 mm. One data axes occupies exactly 75% of total figure width (137.25 mm), with an actual height of 80.56 mm after equal-unit aspect adjustment. The plot has a full rectangular black frame with 0.5 pt lines, bottom/left numeric ticks and UMAP 1 / UMAP 2 titles. No background grid. Axes are equal in physical units. Limits are calculated from the displayed saved coordinate extrema with 5.5% padding, expanded to half-unit boundaries: x [-6.5, 16.5], y [-1.5, 12.0]. All displayed sample centers lie strictly inside those limits. No jitter or coordinate movement is used.

Arial typography follows the user's revised readability targets: axis titles 9 pt, tick labels 8 pt, legend/header 8.5 pt. These explicit revised targets take precedence over the earlier figure's 5-7 pt styling. There is no in-artwork title or panel letter.

Filled single-label reference circles have area 6 pt squared, opacity 0.83 and no outlines. Hollow hybrid circles have area 6 pt squared, neutral #666666 outlines 0.35 pt wide and opacity 0.72. Hollow viral triangles have area 8 pt squared, dark-gray #303030 outlines 0.45 pt wide and opacity 0.87. These point sizes encode sample positions only. Hollow markers reduce occlusion without removing overlapping samples.

The vertical right-hand legend has the reference heading and count, seven single-label colors, Hybrid (multiple labels), and one Selected viral GT-associated loci (287) entry. Long labels are wrapped rather than shrunk. The artwork contains no GC01-GC10 labels, leader lines, community annotation table, percentages, n/N values, bottom notes or other headings. Community results are retained in the separate caption and source tables.

| Display category | HEX | Marker |
| --- | --- | --- |
| Alkaloid | #AA79A2 | filled circle |
| Terpene | #009E73 | filled circle |
| NRP | #0072B2 | filled circle |
| Polyketide | #D55E00 | filled circle |
| RiPP | #56B4E9 | filled circle |
| Saccharide | #C8A44C | filled circle |
| Other | #999999 | filled circle |
| Hybrid (multiple labels) | #666666 | hollow circle |
| Selected viral GT-associated loci (287) | #303030 | hollow triangle |


## Reference-style comparison

The current attachment contained the text brief only. The full-resolution published [Flamholz et al. (2024), Figure 3](https://www.nature.com/articles/s41564-023-01584-8/figures/3) was opened in the in-app browser and visually inspected. The observed image was downloaded through its image locator, copied to qc/Flamholz_2024_Figure3_reference.png and hashed. A local QC crop isolates panel 3a for direct comparison. REFERENCE_STYLE_PROVENANCE.json records the actual retrieval and crop.

The revised graphical organization follows panel 3a's full rectangular numeric axes, dominant landscape scatter, small points, white background and simple vertical class legend. The original paper's data, exact colors and other panels are not used in the new scientific figure. Equal-unit axes preserve our saved projection geometry; the data shape and exact frame proportions therefore remain specific to this dataset. The publisher image is a local style-reference asset only and is not embedded in the delivered PDF/SVG/PNG.

An initial direct HTTP image fetch returned a client-challenge page and the public Europe PMC XML endpoint returned an error. These did not supply the comparison image. The normal browser subsequently displayed the full 2,168 x 1,251-pixel publisher image successfully. No CAPTCHA was solved, authentication altered or paywall bypassed.

## Export and accessibility review

Matplotlib exports a fully vector PDF using font type 42, an SVG retaining text nodes (svg.fonttype=none), and a 600-dpi PNG. Export uses the fixed canvas without tight cropping. The PDF page size and embedded Arial TrueType text were inspected independently; no raster images appear in either vector output. Actual SVG point instances were counted per explicit sample group, including all hybrid circles and viral triangles.

Protanopia, deuteranopia and tritanopia previews use the severity-100 Machado et al. matrices already recorded in the prior plotting work, applied in linear sRGB with gamut clipping. A grayscale preview is also retained. The seven class colors remain the previous final HEX map; filled versus hollow circles and triangle shape distinguish single-label references, hybrids and viral loci. Some class colors remain less distinct under severe CVD or grayscale; the display is not an accessibility certification. No class membership is inferred from neutral hybrid color.

The reference image, normal 150/96-dpi previews and all accessibility views were inspected. The final vector PDF is independently rendered at 300 dpi using Poppler, with its actual command and input/output hashes recorded in PDF_RENDER_EXECUTION.json. The final PDF visual review is recorded in VISUAL_QA.json.

## Executed commands and reproducibility

| Frozen execution directory | Start UTC | Actual child argument array |
| --- | --- | --- |
| executions\plot_integrated_umap_20261028T015614Z | 2026-10-28T01:56:14.482860+00:00 | ["C:\\Users\\archival-user\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe", "-B", "\\\\wsl.localhost\\Ubuntu\\home\\archival-user\\260908_MiBG\\analysis\\viral_gt_suppfig_global_leiden_redesign_20261027\\plot_integrated_umap.py"] |
| executions\prepare_figure_inputs_20261028T015459Z | 2026-10-28T01:54:59.128445+00:00 | ["C:\\Users\\archival-user\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe", "-B", "\\\\wsl.localhost\\Ubuntu\\home\\archival-user\\260908_MiBG\\analysis\\viral_gt_suppfig_global_leiden_redesign_20261027\\prepare_figure_inputs.py"] |
| executions\validate_figure_20261028T020238Z | 2026-10-28T02:02:38.560647+00:00 | ["C:\\Users\\archival-user\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe", "-B", "\\\\wsl.localhost\\Ubuntu\\home\\archival-user\\260908_MiBG\\analysis\\viral_gt_suppfig_global_leiden_redesign_20261027\\validate_figure.py"] |


Each run preserves the exact code/configuration, pre-execution command array, Python/platform version, input-manifest hash, stdout, stderr and final exit status under executions/. SOFTWARE_VERSIONS.json records the actual plotting runtime and resolved Arial font. The final config contains all display choices; rendering is deterministic.

VALIDATION_RESULTS.json records the source, glyph, typography, layout and editable-vector checks. INPUT_MANIFEST_SHA256.tsv and the final OUTPUT_MANIFEST_SHA256.tsv preserve source and output hashes. The checksum manifest and its own verification receipt exclude themselves to avoid recursive hashes. Any later GitHub upload requires explicit approval; all current files remain local.
