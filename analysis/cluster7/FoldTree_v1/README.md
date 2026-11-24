# Cluster 7 Structural Phylogenetic Analysis — Version 1

Structural phylogenetic relationships of viral glycosyltransferase Cluster 7 (C7) were investigated using FoldTree (Moi et al., 2025).

The dataset comprised 143 non-redundant protein structures, including 23 focal C7 proteins, 98 structural-network neighbors, 15 environmental homologues retrieved from MGnify/ESM Atlas, and eight annotated LPS/O-antigen glycosyltransferase references, with one overlapping accession removed.

Protein structures were compared using Foldseek all-versus-all structural alignments. Three phylogenetic trees were constructed using FoldTree's corrected sequence-identity distance (Fident), alignment-normalized TM-score (AlnTM), and local distance difference test (lDDT).

The analysis was performed on August 3, 2026, using the Snakemake-based FoldTree workflow with `custom_structs=True` and `filter=False`. Full protein structures were used without common-core trimming. Distance matrices were used for tree inference through the workflow's QuickTree-based procedure.

The repository preserves the original phylogenetic trees, distance matrices, protein identifiers, alignment results, and relevant computational records. Version 1 represents the original analysis; future revisions will be archived separately.
