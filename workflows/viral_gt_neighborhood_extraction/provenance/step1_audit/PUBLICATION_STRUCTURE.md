# Proposed Step 1 publication structure — review only

Target repository: **LanNa1229/viral-gt-research-record**, public. No push or publication is authorized by this audit. The inventory recommends destinations; it does not copy, move, delete, commit or upload source files. Preserve unrelated local files and the existing repository history. The internal-analyses repository is not used for this proposed publication-essential package.

## A. GitHub: small code and documentation record

Preserve the exact frozen 15-gene code, its synthetic tests, original workflow documentation, verified dependency-version receipt and the new audit conclusions. Keep original broad launchers under an archive path so readers do not accidentally run inference while trying to reproduce Step 1. Proposed structure:

```text
src/viral_gt_mining/
  __init__.py
  extract.py                         # exact baseline_594_v2 source
  analysis.py                        # original locus/genome grouping definitions
  inference.py                       # conditional legacy import dependency, not models
archive/step1_original/
  workflows/viral_gt_mining/
    prepare_inputs.py
    run_full.py                      # provenance; full pipeline, not Step-1-only
    README.md
  window11/
    experiment.py                    # original prepare_main/shortened implementation
    common.py
  retrieval/
    retrieve.py                      # optional mapping provenance
    supplement.py                    # optional exact-version rescue provenance
tests/step1/test_extract.py           # self-contained synthetic annotated-record fixtures
environment/step1/software_versions.json
workflows/viral_gt_neighborhood_extraction/
  README.md                          # notebook summary + stage-specific commands
  provenance_audit.md
  discrepancies.tsv
  reproducibility_manifest.tsv
  source_and_output_checksums.tsv
  zenodo_data_manifest.tsv           # exact scientific-data identifiers + release link later
results/WP01B_coarse_clustering/       # existing upstream record; preserve/link it
```

The inventory's exact existing source paths and proposed destinations control; this tree is a proposal, not a new extraction implementation. Original source code is 14 A-category files, 157,757 bytes total, including optional historical retrieval/launcher records. The code version cannot be represented by the old base commit alone: it was uncommitted and is fixed by the frozen SHA256. Before publication, choose a clear project-code license and retain correct authorship/third-party attribution; no license permission was inferred during this audit.

## B. Zenodo: immutable scientific inputs, outputs and provenance

Recommend one versioned, publicly accessible data release, linked from GitHub, split into a few archives rather than 1,748 separate uploads:

```text
step1_inputs_and_sources.tar.gz
  inputs/                            # exact 594-list, membership and mapping manifests
  sources/genbank/                   # exact 517 historical files, ~195.68 MB
  sources/alias_evidence/             # 593 raw evidence files, ~13.21 MB, plus URL metadata
step1_outputs_and_provenance.tar.gz
  outputs/baseline_594_v2/            # 757 windows / 10,688 gene rows / mapping / QC
  outputs/window11_594_v1/            # 757 windows / 7,967 gene rows / 757 pairing rows
  provenance/                        # exact code/input identities and execution receipts
  validation/                        # original and new independent audit records
reproducibility_manifest.tsv
SHA256SUMS
README.md
```

These archive names are proposed only; no source-data archives were built. Every recommended scientific source file has an individual size and SHA256 in `reproducibility_manifest.tsv`; no directory is represented as if it had a file checksum. Current B recommendations total 231,234,854 bytes (about 220.52 MiB uncompressed), including optional FASTA exports, receipts and HTTP metadata. Required/conditional/optional status is explicit per row. Source GenBank files and exact alias evidence dominate this footprint; this is far smaller than archiving the workspace, environments or model databases.

Retain all 517 GenBank sources if reproducing the original source manifest and its no-window records exactly. A 515-source subset can reproduce the window-bearing input population but will not reproduce the 517-file manifest; label that distinction. Do not replace .2 historical Tupanvirus records with current .3 annotations. Preserve the raw 594 membership and both UniParc identifiers, including the unresolved anchor. Full gene TSVs already contain available protein sequences, so FASTA exports are convenience duplicates rather than a separate indispensable input.

Existing absolute paths in GenBank and alias manifests require a documented portable relocation map. Preserve original manifests unchanged and record relocated manifests separately after validation. Unmodified alias-manifest execution checks evidence-file hashes, so raw evidence bytes cannot be dropped merely because their current location is a cache. Source accessions, retrieval URLs, versions and provenance should remain attributable. Redistribution/license details for external database records were not independently audited; review them before release rather than treating the archive as automatically CC0.

## Critical distinction: evidence archive versus executable minimal replay

**15-gene extraction:** the existing `extract.py` CLI can run offline from the frozen anchor list, exact GenBank annotations, combined original/supplemental aliases, their raw evidence and a correctly relocated path manifest. It requires only the recorded Python/Biopython environment and stdlib; no model or ESM package is needed. The original historical launcher `run_full.py` is broader and must not be advertised as an extraction-only command. This audit has not generated relocated manifests or executed the CLI.

**11-gene preparation:** the historical `experiment.py prepare` is coupled to three prior preservation inventories (676 baseline, 97 interpretation and 47 retrieval-rescue files), original grouping tables, frozen helper modules and NumPy. Those gates are scientifically useful but pull unrelated later-step artifacts into the exact execution closure. The proposed small release is sufficient to inspect its exact source and validate archived subset identities, but is **not yet a complete unmodified replay package** for that launcher. Do not call it standalone or silently remove checks. Before a fully executable minimal release, either archive the entire historical gate closure or approve a separately identified extraction-only wrapper with regression against all 757 frozen pairs. No new wrapper or extraction was implemented in this forensic audit.

The existing 565 baseline and 575 native-11 locus tables, source code and independent interval-membership checks establish the grouping counts without model reruns. If future readers need all later classifier fields regenerated, that belongs to a later inference release with its own model/reference licensing and dependencies, not this Step-1 audit.

## C. Not needed for a Step-1 upload

Keep these locally; this is not a deletion recommendation:

- The replaceable 119,906,304-byte annotations.sqlite index, which can be rebuilt from exact GenBank inputs.
- ESM/model checkpoints, vector caches, databases, reference arrays, prediction/retrieval outputs and final Sankey artwork that belong to later steps.
- Whole Conda/venv environments, temporary/incomplete downloads, Matplotlib/font caches, Git internals and credentials.
- Identical working copies when frozen snapshots are already the authority.
- Baseline_v1 and other exploratory analyses unless a manuscript claim specifically depends on them; preserve their provenance separately rather than silently treating them as disposable.
- External assets without established redistribution permission. No particular scientific source was determined legally nonredistributable by this audit; unknown permission is recorded as pending review, not guessed.

Six concrete C-category examples have exact sizes/checksums in the inventory. This is a scoped inventory, not exhaustive cache enumeration or permission to delete anything.

## Before any future upload

Review the discrepancy report, exact archive scope, portable path handling, project/third-party license decisions and 11-gene replay dependency gap. Keep the public repository public. Obtain the user's separate explicit approval for the exact GitHub push and any Zenodo release. No upload is performed now.
