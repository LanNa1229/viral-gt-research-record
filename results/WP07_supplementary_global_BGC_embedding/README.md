# Supplementary BGC embedding reproducibility package

This archive accompanies the Viral GT manuscript’s supplementary figure of the joint biosynthetic representation space. The figure overlays 2,502 matched MIBiG/BGC-Prophet reference BGCs and 287 selected viral locus representatives. It uses saved UMAP coordinates and the global Leiden baseline computed in the original 320-dimensional embedding space.

- [Final figures](07_supplementary_figure/figures) · [Caption](07_supplementary_figure/provenance/SUPPLEMENTARY_FIGURE_LEGEND.md) · [Figure source data](07_supplementary_figure/source_data)
- [Reproduction instructions](REPRODUCTION.md) · [Methods record](METHODS_RECORD.md) · [Data dictionary](DATA_DICTIONARY.md)
- [Global Leiden files](06_global_leiden/README.md) · [Public dependency manifest](provenance/PUBLIC_DEPENDENCIES.tsv) · [Data lineage](DATA_LINEAGE.tsv)
- [Release checksum manifest](provenance/ARCHIVE_SHA256.tsv) · [Validation report](provenance/VALIDATION_REPORT.md) · [External inputs](provenance/EXTERNAL_DEPENDENCIES.tsv)

Pin the repository to the commit that contains this README. Relative paths in this package are resolved at that commit. The input manifest gives immutable commits for files reused elsewhere in this repository. The original figure exports and scientific outputs retain their historical bytes, with the two mixed-cohort summary tables published as global-only extracts. The source map records the original hashes, public hashes and publication treatment. The historical global Methods and execution evidence are scoped extracts, and workstation identifiers and personal contact details were removed from copied provenance records.

The archived `BGC_train_dataset_classify.csv` is unchanged from Zhang et al.’s version 1 Figshare release of BGC-Prophet, licensed under CC BY 4.0. The package also includes the MIT notice from the original BGC-Prophet source tree. [Figshare release](https://doi.org/10.6084/m9.figshare.28703642.v1) · [Original CSV](https://ndownloader.figshare.com/files/53346437) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Reference annotations also draw on MIBiG releases under CC BY 4.0. [MIBiG](https://mibig.secondarymetabolites.org/)
