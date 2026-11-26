# Viral Glycosyltransferase Research Record

**Project:** AI-driven viral glycosyltransferase discovery and viral glycan biosynthetic loci  
**Repository initialized:** 2026-07-23  

## 0. Overview 


This repository contains computational workflows used for
the discovery and characterization of viral glycosyltransferases (GTs)
from large-scale viral protein sequence and structure datasets.

The workflow integrates:

- Foldseek structural similarity search
- FoldMason multiple structure alignment
- MMseqs2 clustering
- CAZy annotation
- Cytoscape similarity network analysis
- AlphaFold3 structure prediction
- Genomic neighborhood analysis

The goal is to identify previously uncharacterized viral GT families
and reconstruct potential glycan biosynthesis pathways.



## Manuscript sections and figure data

Figure numbers follow manuscript **v1.5 (2026-10-27)**. The index links to public source data, figure exports, and workflows where available. Some panels were assembled manually. For those panels, source materials may include data tables, structure files, or editable figure files.

### Main figures

| Results section | Figure | Linked materials |
| --- | --- | --- |
| 1. Structure-guided discovery of viral GT candidates | **Figure 1** · Structure-guided discovery and clustering | [Figure folder](workflows/WP01_section1/by_figure/structure_guided_discovery_and_clustering/) · [Section 1 workflow](workflows/WP01_section1/README.md) · [Figure source data on Zenodo](https://doi.org/10.5281/zenodo.23626546) |
| 2. Viral GTs across the GT-A structural landscape | **Figure 2A/B** · Global network and viral placement | [Network and placement workflow](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/workflows/WP02_global_SSN_leiden_viral_placement/README.md) · [Zenodo data record](https://doi.org/10.5281/zenodo.23720756). The full dataset ZIP requires access. |
| 2. Viral GTs across the GT-A structural landscape | **Figure 2C/D** · Local networks and conservation maps | Manual analysis |
| 3. Cluster 8 structural conservation and genomic context | **Figure 3** · Cluster 8 structures and conservation | Manual analysis |
| 3. Cluster 8 structural conservation and genomic context | **Figure 4** · Neighborhood heatmap | [Heatmap package](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/workflows/WP08_cluster8_gene_neighborhood_heatmap/README.md) · [Source matrix, evidence table, plots, and script map](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/workflows/WP08_cluster8_gene_neighborhood_heatmap/docs/FIGURE_SOURCE_MAP.tsv). |
| 4. DNA-associated Cluster 8 loci | **Figures 5–6** · Locus diagrams and structural comparisons | Manual analysis |
| 5. Embedding analysis of viral GT-associated loci | **Figure 7A** · Candidate-screening histogram | [Related screening workflow](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/workflows/viral_gt_candidate_selection/README.md). |
| 5. Embedding analysis of viral GT-associated loci | **Figure 7B** · Taxonomic distribution | [242-genome Sankey workflow](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/workflows/viral_gt_sankey_visualization/README.md) · [Source tables and figure exports on Zenodo](https://doi.org/10.5281/zenodo.23706595). |
| 5. Embedding analysis of viral GT-associated loci | **Figure 7C** · Viral communities and reference glycan chemotypes | [Figure package](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/results/WP07_main_figure7C/README.md) · [Reproduction instructions](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/results/WP07_main_figure7C/REPRODUCTION.md). The package includes the unified source table and PDF, SVG, and PNG exports. |
| 5. Embedding analysis of viral GT-associated loci | **Figure 7D** · Cluster 7 neighborhood comparisons | Manual analysis |

### Supplementary figures

| Results section | Supplementary figure | Linked materials |
| --- | --- | --- |
| 1. Structure-guided discovery of viral GT candidates | **Figure S1** · Cluster membership network | [Figure folder](workflows/WP01_section1/by_figure/cluster_membership_network/) |
| 1. Structure-guided discovery of viral GT candidates | **Figure S2** · Positive-control benchmark | [Figure folder](workflows/WP01_section1/by_figure/positive_control_benchmark/) · [Public exclusions for the S2B model plate](PUBLIC_EXCLUSIONS.tsv) |
| 1. Structure-guided discovery of viral GT candidates | **Figure S3** · Small-cluster distribution and coherence | [Figure folder](workflows/WP01_section1/by_figure/small_cluster_distribution_and_coherence/) |
| 1. Structure-guided discovery of viral GT candidates | **Figure S4** · Singleton length and confidence | [Figure folder](workflows/WP01_section1/by_figure/singleton_length_and_confidence/) |
| 2. Viral GTs across the GT-A structural landscape | **Figure S5** · CAZy-family organization | [Network and placement workflow](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/workflows/WP02_global_SSN_leiden_viral_placement/README.md) (upstream data). |
| 2. Viral GTs across the GT-A structural landscape | **Figure S6A/B** · GT31 and R699 comparisons | [S6A GT31 core workflow](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/workflows/WP06_GT31_conservation/README.md) · [S6B R699 core-mapping archive](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/workflows/WP06_R699_dual_context/README.md). The R699 archive documents missing inputs and outputs for the core-mapping analysis. |
| 2. Viral GTs across the GT-A structural landscape | **Figures S7–S8** · C027 and Cluster 7 phylogenies | [S8 Cluster 7 FoldTree archive](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/analysis/cluster7/FoldTree_v1/README.md) contains trees and distance matrices. |
| 3. Cluster 8 structural conservation and genomic context | **Figures S9–S11** · Cluster 8 clustering and alignments | Manual analysis |
| 4. DNA-associated Cluster 8 loci | **Figures S12–S17** · Precursor chemistry and structural annotations | Manual analysis |
| 5. Embedding analysis of viral GT-associated loci | **Figure S18** · Joint embedding and global communities | [Supplementary embedding package](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/results/WP07_supplementary_global_BGC_embedding/README.md). |
| 5. Embedding analysis of viral GT-associated loci | **Figure S19** · Reference saccharide chemotypes | [Chemotype package](https://github.com/LanNa1229/viral-gt-research-record/blob/c4d46eb5ece3b399225aaaee7d60ef30e506461e/results/WP07_reference_saccharide_chemotype_analysis/README.md). Its README uses the earlier label “Supplementary Figure 2.” |

## 1. Purpose

This repository is the version-controlled research record for computational and integrative work on viral glycosyltransferase discovery, structural bioinformatics, viral glycan biosynthetic loci, and experimental candidate selection.

The repository is designed to:

1. preserve reproducible code, parameters, metadata, notebooks, and analytical decisions;
2. document the chronology and provenance of scientific questions, hypotheses, methods, outputs, and interpretations;
3. connect analyses to manuscript claims, figures, tables, and experimental decisions;
4. identify specific technical and intellectual contributions using the CRediT taxonomy;
5. maintain a professional, orderly, and institutionally accessible research record;


## 2. Scientific scope

Current work packages include:

| Work package | Scope |
|---|---|
| WP01 | Viral glycosyltransferase discovery framework |
| WP02 | Protein-structure curation and metadata normalization |
| WP03 | Foldseek structural comparison and quality control |
| WP04 | Structure-similarity network construction |
| WP05 | Weighted Leiden clustering and community interpretation |
| WP06 | Viral enrichment and cross-domain structural placement |
| WP07 | Genomic-neighborhood and viral glycan-locus analysis |
| WP08 | Cluster 7 and Cluster 8 family-level analysis |
| WP09 | Jumbo-phage DNA modification and sulfation-module analysis |
| WP10 | Phylogenetic or structure-informed evolutionary placement |
| WP11 | AlphaFold structure generation and model quality assessment |
| WP12 | Manuscript claims, figures, tables, and biological interpretation |
| WP13 | Wet-lab validation candidate selection and experimental design |

## 4. Repository architecture

```text
.
├── README.md
├── CONTRIBUTIONS.md
├── PROJECT_TIMELINE.md
├── CREDIT_CONTRIBUTION_LEDGER.tsv
├── CLAIM_FIGURE_EVIDENCE_MATRIX.tsv
├── DATA_TRANSFER_MANIFEST.tsv
├── .gitignore
├── .gitattributes
├── notebooks/
│   ├── template/
│   ├── active/
│   ├── completed/
│   └── retrospective/
├── src/
├── workflows/
├── configs/
├── metadata/
│   ├── schemas/
│   ├── curated/
│   └── manifests/
├── data/
├── results/
│   ├── tables/
│   ├── networks/
│   ├── structures/
│   └── logs/
├── figures/
│   ├── source/
│   ├── panels/
│   ├── final/
│   └── figure_manifests/
├── manuscript_claims/
├── reports/
├── environment/
├── archive_manifests/
├── retrospective_records/
├── transfer_manifests/
├── docs/
└── scripts/
```

## 5. What belongs in Git

Track in Git:

- source code, workflow files, and configuration files;
- human-readable Quarto or Markdown notebooks;
- small curated metadata tables and schemas;
- parameters and command records;
- environment specifications;
- small final summary tables;
- figure-generation scripts;
- vector figures and reasonably sized final figures;
- manuscript claim files;
- CRediT and claim-evidence ledgers;
- transfer manifests;
- archive manifests and checksums;
- README files that explain externally stored datasets.

.

For each externally stored dataset, create a manifest under `metadata/manifests/` containing:

- dataset ID;
- institutional storage location;
- responsible steward;
- creation or acquisition date;
- source;
- file count and total size;
- SHA-256 manifest location;
- access restrictions;
- notebook and commit that used the dataset.





## Section 1 reproducibility release

[Browse all five figure folders](workflows/WP01_section1/by_figure/) · [Panel-level evidence map](reports/section1_archive/figure_evidence_map.tsv)

The [Section 1 workflows](workflows/WP01_section1/README.md) preserve the code,
lightweight canonical results, figure-generation utilities and validation for
Figure 1 and Supplementary Figures 1–4. The downstream baseline comprises
1,727 structures retrieved using the original 447-residue MGAT2 Foldseek query.
The selected Section 1 snapshot uses MIT for original code and CC BY 4.0 for
original data/documentation, subject to its third-party notices and exclusions;
these terms do not relicense unrelated repository content.

- [Versioned Section 1 release](https://github.com/LanNa1229/viral-gt-research-record/releases/tag/section1-v1.0.0)
- [Software archive](https://doi.org/10.5281/zenodo.23626496)
- [Dataset archive and Figure Source Data](https://doi.org/10.5281/zenodo.23626546)

Unconfirmed institutional AF3 model inputs/databases and the qualitative model
plate are excluded from the public subset; numerical control evidence is
retained. Missing truncated-core files are a documented, non-blocking historical
limitation. See the scoped workflow and public exclusion manifest for details.
