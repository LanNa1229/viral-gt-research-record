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


