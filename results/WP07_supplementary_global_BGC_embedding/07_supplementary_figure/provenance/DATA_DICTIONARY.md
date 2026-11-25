# Source-data and display dictionary

FIGURE_SOURCE_DATA.tsv has one row for each of 2,789 unique displayed samples. The twelve non-rendering fields from the preceding figure source are preserved exactly: sample_id, sample_type, global_leiden_community, UMAP1, UMAP2, reference_class_original, reference_class_labels_json, n_reference_class_labels, is_hybrid_reference, viral_locus_id, genome_unit_id and representative_window_id. Stable identifiers are join keys; original UMAP decimal strings and reference label ordering are unmodified.

Rendering fields:
- plot_marker: filled_circle for a single-label reference, hollow_circle for a hybrid reference, hollow_triangle for a selected viral locus.
- rendering_group: original single reference class, Hybrid (multiple labels), or Selected viral GT-associated loci (287). The hybrid label describes presentation only and never replaces original scientific annotations.
- plot_edge_or_fill_hex: fill HEX for single-label references or outline HEX for hollow glyphs.
- plot_facecolor: class HEX or none.

The seven original functional labels and all 467 multilabel memberships remain in reference_class_original and reference_class_labels_json. A hollow hybrid glyph does not indicate missing labels.

COMMUNITY_ANNOTATION_KEY.tsv and REFERENCE_CLASS_MARGINAL_COUNTS.tsv are byte-identical copies of the preceding verified tables. They retain ten original GC labels, per-community reference/viral counts, dominant marginal labels and all 70 class memberships/fractions. Multilabel marginal fractions can exceed 100% in sum. These tables support the separate figure legend and are not rendered as an in-artwork table.

PLOTTED_SAMPLE_RECEIPTS.tsv records each actual plotted ID, display group, marker, unchanged coordinate strings and original class memberships. RENDER_GROUP_RECEIPTS.tsv records the nine plotted groups and their marker parameters. LAYOUT_RECEIPTS.json records actual axes geometry, ticks, font sizes and measured text bounds. COLOR_PALETTE.tsv includes seven annotation colors plus neutral hybrid and viral display colors.
