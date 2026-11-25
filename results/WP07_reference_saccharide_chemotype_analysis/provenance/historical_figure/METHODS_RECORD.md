# Methods record: Saccharide reference chemotype figure

## Scope and controlling sources

This is source-data verification and visualization only. No new chemical annotation, clustering, embeddings, model inference, nearest-neighbor calculation or enrichment test was run. All original analysis/curation files were read without modification. No commit or push was performed.

The controlling source is reference_2502_chemotype_master.tsv, 2,502 biological BGCs x 91 fields, from analysis/viral_gt_reference_chemotype_curation_20261027/. Its frozen SHA-256 is 5f86b03ff3e071eec4fcafa02eac36c981f1561777a1d1acbe6f69e51315cbc9. The dictionary, Methods, reference QC and complete prior validation JSON were inspected. All 10,036 prior checks were programmatically parsed and PASS status confirmed; repeated per-accession check families are summarized in PRIOR_VALIDATION_REVIEW.json rather than treating truncated console output as a complete review.

Every consumed curation artifact was checked against its sealed output manifest, whose own SHA-256 was verified. Existing source decisions were used only to verify unchanged primary chemotypes and confidence. The original global assignment grid and original reference class labels were checked against the pinned archival checksum file, tied to private archival commit 284d4f1655cf036c1890bdffcb4aa2dc93f6e41e. Baseline: k20_r1_s42, original-320D k=20, resolution=1.0, seed=42. Existing MIBiG release identity is MIBiG 3.1; no current live release or new literature source was substituted.

## Populations, categories and denominators

Saccharide selection is membership in the original broad class-token list, including hybrids. Exactly 195 unique accessions qualify: 100 GC08, 24 GC05, 71 in neither GC08 nor GC05. These sets are pairwise disjoint and their union is the selected 195. The earlier 95 non-GC08 references include the 24 GC05 references and are not used as a separate population. Non-Saccharide references, including the five non-Saccharide GC08 members and 29 non-Saccharide GC05 members, are excluded from these denominators.

Existing mutually exclusive primary_chemotype values are mapped through an explicit 15-to-7 lookup, without fuzzy matching or reinterpretation. The exact source spellings, including LPS/LOS abbreviations, are retained in CHEMOTYPE_GROUP_MAPPING.tsv. All 91 original master fields are preserved in ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv; four additional fields record the display population, consolidated chemotype, category order and explicit membership in the operational four-category set. OUTSIDE_DEFINED_SET means outside those named bins, not absence of surface chemistry. Unresolved primary chemotypes have an explicit UNRESOLVED_PRIMARY_CHEMOTYPE membership state.

The operational surface/extracellular definition is categories 1-4 only. The two GC08 phosphonoglycans remain in Other/unresolved; the source's broader 74/100 flag is neither substituted nor overwritten. Glycopeptidolipid and other miscellaneous classes also stay in the requested final bin. Unresolved records remain in denominators and in that bin, without conversion to biological negatives.

The 21 population-category cells were reconstructed from individual accessions. Historical numbers in the request were used solely as validation targets; any discrepancy would stop preparation before plotting. All expected values matched. Bar widths are n/N x 100, using each population's Saccharide BGC count as N. Exact rational accumulation verifies 100% for every row before conversion to display floats. FIGURE_SOURCE_DATA.tsv retains counts, fractions, percentages, cumulative boundaries and exact accession lists for every cell, including zeros.

| Chemotype group | GC08 (100) | GC05 (24) | Other GCs (71) |
| --- | --- | --- | --- |
| Capsule | 31 | 4 | 1 |
| EPS/matrix polysaccharide | 19 | 2 | 1 |
| LPS/LOS/O-antigen | 18 | 3 | 1 |
| S-layer/other surface glycan | 4 | 0 | 0 |
| Aminoglycoside/aminocyclitol | 20 | 7 | 7 |
| Glycosylated specialized metabolite | 4 | 5 | 51 |
| Other/unresolved | 4 | 3 | 10 |


## Graphical design

One 100% stacked horizontal bar chart, 183 x 95 mm, with bars ordered GC08, GC05, Other GCs. Category order is identical across rows. The first four consecutive segments form the operational surface/extracellular set, with dark gold, the previous UMAP Saccharide gold, orange and light amber. Blue, purple and gray complete the stack. No gradients, decorative borders, explanatory in-artwork paragraphs or numerical table are used.

Arial: x-axis title and group labels 8.5 pt, percentage ticks 8 pt, legend 8 pt. The axis is fixed at 0-100%, with labeled 0%, 25%, 50%, 75%, 100% ticks and the title Saccharide reference BGCs (%). White 0.35 pt separators preserve boundaries without inflating segment widths. Bar height is 0.56 category units, leaving clear gaps. A two-column legend below the chart keeps the four warm categories together and uses unwrapped short labels.

Zero-count S-layer/other-surface segments in GC05 and Other GCs remain zero-width, with explicit zero source rows; no minimum-width inflation is applied. The smallest positive segment, 1/71, is 1.40845% or approximately 2.04 mm across the 145 mm plotting width. All exact values remain in the source table rather than labels inside narrow bars.

| Chemotype group | Final HEX |
| --- | --- |
| Capsule | #6F4E0B |
| EPS/matrix polysaccharide | #C8A44C |
| LPS/LOS/O-antigen | #D55E00 |
| S-layer/other surface glycan | #F3D7A5 |
| Aminoglycoside/aminocyclitol | #0072B2 |
| Glycosylated specialized metabolite | #B68BC4 |
| Other/unresolved | #A8A8A8 |


## Palette and export verification

The initial palette was inspected under severity-100 Machado protanopia, deuteranopia and tritanopia matrices in linear sRGB plus grayscale. Blue/purple separation was weak under protanopia, so purple was lightened to #B68BC4 and capsule gold darkened to #6F4E0B. The final palette retains #C8A44C from the previous UMAP. Exact matrices and pairwise CIE76 diagnostics are preserved. The final simulations were visually checked. Diagnostic differences do not certify accessibility for every observer, and grayscale cannot uniquely encode all seven categories; fixed order and segment boundaries remain helpful.

Exports use a fixed canvas, 600-dpi PNG, SVG text nodes and PDF font type 42 with embedded Arial TrueType text. No raster images occur in the vector files. The final PDF is independently rendered at 300 dpi through Poppler, with the actual command and input/output hashes recorded. Both 150-dpi and 96-dpi previews preserve the intended layout. Actual SVG rectangle widths are independently checked against source percentages, including zeros.

## Annotation limitations

The selected cohort contains 168 SOURCE_ANNOTATION_SUPPORTED, 24 LITERATURE_SUPPORTED, one PROVISIONAL_RECORD_LEVEL and two UNRESOLVED annotations. None of these status names establishes uniform independent human manual verification.

- BGC0000262 (Other GCs): aglycone products and hybrid/tailoring context leave product scope unresolved.
- BGC0000727 (GC08): MIBiG indigoidine product versus an auricin/deoxysugar-related sequence description remains an unresolved source conflict.
- BGC0001765 (GC08): the pinned salecan record is strain S33 while cited experimental chemistry concerns ZX09; the EPS category is retained provisionally without transferring experimental claims.

The original detailed notes, statuses, source URLs and record hashes remain intact. ANNOTATION_REVIEW_FLAGS.tsv lists these three records. Broader curation limitations remain in the original curation reports. This task did not attempt to resolve them.

## Actual execution evidence

| Execution folder | Start UTC | Actual child argument array |
| --- | --- | --- |
| executions\plot_chemotype_bars_20261028T022840Z | 2026-10-28T02:28:40.884661+00:00 | ["/redacted-user/.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe", "-B", "/redacted-home/260908_MiBG\\analysis\\viral_gt_suppfig_chemotype_20261027\\plot_chemotype_bars.py"] |
| executions\plot_chemotype_bars_20261028T023154Z | 2026-10-28T02:31:55.111674+00:00 | ["/redacted-user/.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe", "-B", "/redacted-home/260908_MiBG\\analysis\\viral_gt_suppfig_chemotype_20261027\\plot_chemotype_bars.py"] |
| executions\plot_chemotype_bars_20261028T023737Z | 2026-10-28T02:37:38.096463+00:00 | ["/redacted-user/.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe", "-B", "/redacted-home/260908_MiBG\\analysis\\viral_gt_suppfig_chemotype_20261027\\plot_chemotype_bars.py"] |
| executions\prepare_figure_inputs_20261028T022659Z | 2026-10-28T02:26:59.778110+00:00 | ["/redacted-user/.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe", "-B", "/redacted-home/260908_MiBG\\analysis\\viral_gt_suppfig_chemotype_20261027\\prepare_figure_inputs.py"] |
| executions\prepare_figure_inputs_20261028T023429Z | 2026-10-28T02:34:29.369899+00:00 | ["/redacted-user/.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe", "-B", "/redacted-home/260908_MiBG\\analysis\\viral_gt_suppfig_chemotype_20261027\\prepare_figure_inputs.py"] |
| executions\validate_figure_20261028T024156Z | 2026-10-28T02:41:56.362391+00:00 | ["/redacted-user/.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe", "-B", "/redacted-home/260908_MiBG\\analysis\\viral_gt_suppfig_chemotype_20261027\\validate_figure.py"] |


run_stage.py freezes code/configuration and records pre-execution command arrays, Python/platform, input-manifest hash, stdout, stderr and final exit status. The initial plotting attempt produced the artwork but failed to serialize NumPy tick integers into its layout receipt; that implementation error was corrected by converting ticks to native floats, and its unsuccessful log/code are retained. The final palette adjustment has its own preparation and plotting records. No counts or chemistry were changed during these display iterations.

INPUT_MANIFEST_SHA256.tsv freezes all consumed sources before calculations; every source is rehashed after verification. SOFTWARE_VERSIONS.json records actual plotting libraries and resolved Arial. VALIDATION_RESULTS.json contains the independent source, mathematical, geometry, typography and vector-export checks. Output hashes are sealed after all execution files and final visual review records exist. The manifest and its own verification receipt are excluded from recursive hashing.

All files remain local. Any later private archival push requires explicit scoped approval and PRIVATE visibility verification; the public repository is unchanged.
