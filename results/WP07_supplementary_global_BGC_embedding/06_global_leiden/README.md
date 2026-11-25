# Global Leiden reproducibility inputs

The historical global joint baseline contains 2,789 samples with 320 features each, 2,502 reference BGCs and 287 selected viral representatives. It was clustered on the weighted union of the 20 nearest cosine neighbors using resolution 1.0 and seed 42. The saved UMAP coordinates are used for visualization.

The [public dependency manifest](../provenance/PUBLIC_DEPENDENCIES.tsv) identifies the exact public revision or package-relative location for each retained scientific dependency. It covers the original input order, k20 edge and node tables, graph diagnostics, global assignment grid, optimizer record, global stability summaries and executed code snapshots. Two original cross-cohort tables are provided as global-row extracts. Their original and published digests are recorded in the manifest. The Methods and evidence matrix are extracts covering the global workflow. The run receipts and code snapshots preserve the shared execution context.

Follow the [reproduction instructions](../REPRODUCTION.md) to verify the edges and assignments from the archived vectors and optionally replay Leiden with the original engine.
