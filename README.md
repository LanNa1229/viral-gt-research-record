# Viral Glycosyltransferase Research Record

**Repository status:** Private, unpublished research  
**Project:** AI-driven viral glycosyltransferase discovery and viral glycan biosynthetic loci  
**Institution:** University of Georgia  
**Repository initialized:** 2026-07-23  
**Primary purpose:** Reproducibility, research-data stewardship, provenance, manuscript traceability, and transparent contribution documentation

## 1. Purpose

This repository is the version-controlled research record for computational and integrative work on viral glycosyltransferase discovery, structural bioinformatics, viral glycan biosynthetic loci, and experimental candidate selection.

The repository is designed to:

1. preserve reproducible code, parameters, metadata, notebooks, and analytical decisions;
2. document the chronology and provenance of scientific questions, hypotheses, methods, outputs, and interpretations;
3. connect analyses to manuscript claims, figures, tables, and experimental decisions;
4. identify specific technical and intellectual contributions using the CRediT taxonomy;
5. maintain a professional, orderly, and institutionally accessible research record;
6. support appropriate data access and reuse without separating files from their scientific context or provenance.

This repository is not a mechanism for withholding University research data. Research data and supporting records will be shared through authorized channels when requested for a legitimate project purpose, with a transfer manifest, file versions, checksums, and an accompanying README sufficient to preserve reproducibility and provenance.

## 2. Governance statement

- The repository must remain private while it contains unpublished research.
- University research data are treated as institutional research records and are not represented as the personal property of an individual researcher.
- A personal laptop or personal Git-hosting account must not be the sole location of the project record.
- The repository must have an institutionally accessible backup or snapshot on UGA-managed or UGA-approved storage.
- Large raw and intermediate files are stored outside Git. Git tracks their storage location, checksum, version, provenance, and relationship to analyses.
- Git history must not be rewritten to improve chronology. Do not force-push shared branches.
- Historical work must not be backdated. Retrospective records must be explicitly labeled and supported by contemporaneous evidence.
- Contributions are described as specific actions and decisions, not as vague claims of participation.
- CRediT records document contribution roles but do not, by themselves, determine authorship.
- Corrections are made by a new commit that explains the change; prior records are not silently erased.

## 3. Scientific scope

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

Do not ordinarily track in Git:

- raw downloaded proteomes or large FASTA collections;
- full AlphaFold output directories;
- large PDB, CIF, mmCIF, BCIF, or trajectory collections;
- Foldseek databases and large `.m8` alignment tables;
- Cytoscape session files if large or frequently regenerated;
- large all-vs-all edge tables;
- temporary files, caches, checkpoints, or package environments;
- redundant derived outputs that can be regenerated;
- credentials, tokens, private keys, or access-controlled data.

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

## 6. Data tiers

| Tier | Definition | Primary location | Git record |
|---|---|---|---|
| Raw external | Unmodified downloaded source data | UGA-managed storage | source manifest and checksums |
| Raw generated | Original project-generated files | UGA-managed storage | generation record and checksums |
| Intermediate | Regenerable workflow outputs | HPC or UGA storage | workflow, parameters, selected manifests |
| Curated | Reviewed metadata and analysis-ready subsets | Git if small; otherwise UGA storage | table or manifest |
| Final analytical | Tables, statistics, selected networks, figure inputs | Git when reasonably sized | tracked file and notebook |
| Manuscript-ready | Final panels, legends, claim text, release tables | Git plus institutional snapshot | tracked file, tag, snapshot |

## 7. Core record linkage

Every substantive analysis should be linkable through the following chain:

```text
Research question
→ Notebook ID
→ Input dataset ID and checksum
→ Git commit
→ Workflow/configuration
→ Output file and checksum
→ Scientific interpretation
→ Decision
→ Claim ID
→ Figure/table
→ CRediT ledger entry
```

No manuscript claim should be marked `ready` until this chain is complete or its documented limitations are stated.

## 8. Identifier conventions

- Entry: `ENT-VGT-YYYYMMDD-NNN`
- Notebook: `NB-VGT-WP##-YYYYMMDD-NNN`
- Claim: `CLM-VGT-NNN`
- Figure: `FIG##-PANELX`
- Output: `OUT-ENT-VGT-YYYYMMDD-NNN-short-description.ext`
- Retrospective record: `RET-VGT-YYYYMMDD-NNN`
- Transfer: `XFER-VGT-YYYYMMDD-NNN`
- Dataset: `DS-VGT-YYYYMMDD-NNN`
- Milestone tag: `milestone-vMAJOR.MINOR.PATCH-YYYYMMDD`
- Snapshot: `snapshot-YYYYMMDDTHHMMSSZ-<short-commit>`

Use UTC timestamps for machine-generated snapshots and local dates for narrative research entries.

## 9. Branch and commit practice

Protected branch:

- `main`: reviewed, stable research record

Working branches:

- `analysis/WP05-weighted-leiden`
- `analysis/WP09-jumbo-phage-neighborhoods`
- `docs/ENT-VGT-20260723-001-initial-governance`
- `retrospective/RET-VGT-20260723-001-cluster8`
- `figure/FIG03-weighted-leiden`
- `fix/WP03-foldseek-coverage-calculation`

Commit messages should describe the scientific action and object:

- `Identify Cluster 8 viral enrichment using curated structure metadata`
- `Add lineage-aware filtering for jumbo-phage neighborhood analysis`
- `Generate manuscript Figure 3 from weighted Leiden results`
- `Document biological interpretation of sulfation-module candidates`

Avoid messages such as `update`, `fix`, `new results`, `final`, or `final2`.

## 10. Snapshot and backup policy

At each manuscript or analysis milestone:

1. ensure the working tree is clean or explicitly record uncommitted files;
2. create an annotated Git tag;
3. generate a Git bundle containing all refs;
4. generate ZIP and TAR archives from the tagged commit;
5. generate SHA-256 checksum manifests;
6. copy the snapshot to UGA-managed or UGA-approved storage;
7. record the destination, date, responsible person, and verification result in `archive_manifests/`.

Recommended institutional options should be selected according to data classification and project needs. Examples include GACRC storage, UGA OneDrive, Globus-connected UGA storage, Research Institutional File Stores, Secure IFS, or unit-managed storage.

## 11. Retrospective records

Historical work must use the template in:

```text
retrospective_records/RETROSPECTIVE_PROVENANCE_TEMPLATE.md
```

Required opening statement:

> Retrospective provenance entry created on [date]. This entry reconstructs work performed between [date range] using contemporaneous emails, presentation files, manuscript drafts, HPC job records, original file timestamps, and archived outputs.

A retrospective entry must distinguish:

- contemporaneous facts;
- later reconstruction;
- current interpretation;
- unresolved uncertainty.

## 12. Contribution documentation

Contribution records are maintained in:

- `CONTRIBUTIONS.md`
- `CREDIT_CONTRIBUTION_LEDGER.tsv`
- the CRediT section of each notebook
- `CLAIM_FIGURE_EVIDENCE_MATRIX.tsv`

Specific wording is required. For example:

- `Lan conceived the structural-network strategy used to organize sequence-divergent viral GT candidates.`
- `Lan selected the Foldseek coverage and structural-similarity criteria after evaluating control recovery and off-target behavior.`
- `Lan implemented the weighted Leiden workflow and produced the reproducible node, edge, and community audit tables.`
- `Lan interpreted the viral enrichment pattern as evidence for a recurrent virus-associated structural lineage.`
- `Lan connected Cluster 8 neighborhood architecture to the jumbo-phage sulfated-DNA hypothesis.`
- `Lan designed the staged experimental validation strategy for the multi-enzyme DNA-modification module.`

## 13. References and institutional alignment

This repository is intended to support the principles described in:

- UGA Research Data Stewardship policy and guidance: `https://researchdata.uga.edu/policy-guidance`
- UGA Open Scholarship & Data Services: `https://researchdata.uga.edu`
- CRediT Contributor Role Taxonomy: `https://credit.niso.org`

The repository does not replace sponsor requirements, approved data-use agreements, export-control requirements, security rules, laboratory-specific procedures, or formal authorship policies.
