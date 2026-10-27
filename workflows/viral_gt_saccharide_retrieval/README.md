# Step 3: reference-based Saccharide nearest-neighbor retrieval

Published v1.0.0 on 2026-10-27 with explicit creator approval. Creator: Lan Na (given name Lan; family name Na).

This package preserves three missing historical scripts byte-for-byte and supplies a separate portable retrieval-only adapter. It references 26 already-public resources by immutable commit/DOI, exact member path and SHA256. No reference matrices, query embeddings, predictions, already-published source files, model weights, source sequences, or earlier 12,510-row augmented collection are included.

The frozen bank contains 2,502 MIBiG-derived CORE_ONLY BGC accessions and 195 source annotations containing the exact whitespace-delimited label Saccharide. Retrieval uses float64 L2-normalized original 320D vectors, descending cosine similarity and ascending accession for ties. Each of the ten neighbors contributes 0.1 to P@10; multilabel references containing Saccharide count as positive. The historical exploratory screen is **P@10 >= 0.5**, combined with the unchanged frozen Saccharide classifier flag using logical OR.

| Frozen result | 15-gene | 11-gene |
|---|---:|---:|
| Evaluated query windows | 737 | 740 |
| P@10 >= 0.5 | 325 | 337 |
| Classifier-only | 6 | 7 |
| Retrieval-only | 85 | 79 |
| Both | 240 | 258 |
| Neither | 406 | 396 |
| Either | 331 | 344 |

Full retrieval reproduction of all 1,477 saved queries passed: 29,540 neighbor ID/rank/label rows and 4,431 P@5/P@10/P@20 values match exactly. Cosine/scalar comparisons use atol=rtol=1e-12, without byte-identity claims. No BGC-Prophet inference or embedding generation ran. This check does not expand Step 2's model-reproduction claim beyond its 12 previously reproduced representative windows.

The separate dataset supplies four missing byte-preserved historical outputs plus new validation and retrospective comparator sensitivity. The historical 245 candidates and 242-genome Sankey remain unchanged; a hypothetical >0.5 rule yields 230/227 with the same three figure exclusions.

Read EXECUTION_GUIDE.md, ENVIRONMENT.md and HISTORICAL_LIMITATIONS.md. Use portable/step3.py; historical mains are broader context and must not be invoked for this retrieval-only workflow.

Step 1: https://doi.org/10.5281/zenodo.23680510
Step 2: https://doi.org/10.5281/zenodo.23682612
Step 3 dataset: https://doi.org/10.5281/zenodo.23683720 ; public record: https://zenodo.org/records/23683720 . Exact approved 24-file software snapshot: https://github.com/LanNa1229/viral-gt-research-record/commit/232437d59a1cf1a7830dc9541300bf95d2132418 .

Approved original-work licenses: MIT applies only to original software; CC BY 4.0 applies only to original derived data/documentation. Third-party rights, attribution and source terms remain separate. See LICENSE_PROPOSAL.md and THIRD_PARTY_NOTICES.md.

## Immutable publication snapshots

The approved software commit above preserves all 24 reviewed file bytes. This scoped follow-up adds the published DOI, exact software commit and license activation, with updated integrity manifests only. Historical scripts, executable retrieval code, validation records and scientific outputs are unchanged.

The Zenodo record contains Step3_frozen_retrieval_outputs_v1.0.0_DATA.zip (22 files; SHA256 d8488e31e1320070a04ecf7c34f6617c77fc8b2ccb5a6c9540afe23f9659f957) and its checksum sidecar. The reviewed software ZIP SHA256 is 8abfae8d0ed2795b946b0e2e861f3570a78ec3ddba8158f81164b9adaaefddf0; its files are represented by the exact GitHub commit, without a redundant Zenodo software archive.

Preparation-time pending wording remains in the immutable archive/snapshot. Public metadata and this current documentation reflect subsequent explicit creator approval. All 26 public dependencies remain pinned in references/public_resources.tsv and resources.json. Historical P@10 >= 0.5 and the 245/242 cohorts remain intact; 230/227 is retrospective sensitivity only. No scientific calculations were rerun during publication.
