# Saccharide chemotypes across global Leiden communities

This archive accompanies Supplementary Figure 2 of the Viral GT manuscript. It asks whether Saccharide-associated reference biosynthetic gene clusters (BGCs) in different global Leiden communities represent distinct glycan biosynthetic chemotypes. The archive links every plotted count to a stable MIBiG accession and retains the original figure exports and source table.

The reference population comprises 2,502 biological BGCs matched to MIBiG 3.1 and the BGC-Prophet representation used in the published global embedding workflow. One reference row represents one accession. The 195 Saccharide-associated references include 65 hybrid entries with Saccharide among their original class tokens. Community membership comes from the saved global joint Leiden baseline in the original 320-dimensional representation, using 20 neighbors, resolution 1.0 and seed 42. The metadata, reference embeddings and assignment grid are pinned to immutable public commits.

The historical workflow extracted product, compound, gene and evidence fields from frozen MIBiG JSON and GenBank records. A priority set of 229 accessions combined the 195 Saccharide references with the reference members of GC08 and GC05. Source rules and accession-specific evidence decisions assigned primary categories while retaining original labels, secondary tags and unresolved records. The figure maps the primary categories to seven display groups, assigns each plotted accession to one group and divides the category count by the Saccharide total in its community partition.

The four display categories for capsule, extracellular or matrix polysaccharide, LPS/LOS/O-antigen and S-layer or other surface glycans account for **72/100 in GC08 (72.0%)**, **9/24 in GC05 (37.5%)** and **3/71 in the remaining communities (4.2%)**. This comparison describes the composition of the reference populations. The broader historical surface/extracellular flag includes phosphonoglycans and has a different definition.

- [Inputs](inputs/) holds the reference crosswalk and source acquisition guides. The [public dependency manifest](provenance/PUBLIC_DEPENDENCIES.tsv) records exact versions and digests.
- [Scripts](scripts/) and [annotation evidence](annotation_evidence/) retain historical extraction and curation code, rule versions, accession decisions, evidence cards and literature citations.
- [Outputs](outputs/) holds the reference master, compound and gene tables, composition summaries and validation records.
- [Final PDF, SVG and PNG](figures/final/), [figure source data](figures/source_data/FIGURE_SOURCE_DATA.tsv) and [caption](figures/SUPPLEMENTARY_FIGURE_LEGEND.md) accompany the [methods](METHODS_RECORD.md) and [reproduction commands](REPRODUCTION.md).

The assignments describe reference BGC products. Viral locus products require separate experimental evidence. Automated source rules and AI-assisted interpretation of record and literature evidence were used for computational curation. The archive records the evidence level and uncertainty of individual assignments. [MIBiG 3.1 release](https://zenodo.org/records/13367811).
