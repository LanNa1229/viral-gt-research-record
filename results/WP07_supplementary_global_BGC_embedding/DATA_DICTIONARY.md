# Data dictionary

| Stage | Archived material |
| --- | --- |
| Reference embeddings | Matched BGC-Prophet reference vectors and row identifiers, reference class labels and reference input provenance |
| Viral windows | Frozen selected window and CDS records, inference output tables and associated Methods |
| Viral embeddings | Original embedding array, per-window row manifest and retained numeric feature bundles |
| Representative loci | 556 viral representatives, their parent windows and arrays, and the saved joint UMAP coordinates |
| Global Leiden | Selected 287 representative order, baseline k20 graph, global assignment grid, class and membership summaries, graph diagnostics, execution code and optimizer provenance |
| Supplementary figure | Source table, palette, marginal counts, community key, exact historical plotting code, caption and original PNG/PDF/SVG files |

The figure source data contains one row per displayed sample. `sample_id` identifies a reference BGC or a viral locus. `UMAP1` and `UMAP2` are the saved display coordinates. `global_leiden_community` identifies the baseline community from the original embedding coordinates. `reference_class_original` retains the input reference classification. The selected-287 table preserves the row order used to assemble the global clustering input.

The [lineage table](DATA_LINEAGE.tsv) tracks stage-level inputs and outputs. The [public dependency manifest](provenance/PUBLIC_DEPENDENCIES.tsv) maps historical clustering dependencies to public locations. The [source map](provenance/RELEASE_SOURCE_MAP.tsv) records identity, checksums and changes made for release. The [screening ledger](provenance/SCREENED_OMISSIONS.tsv) records omissions and superseded private-archive records by source checksum.
