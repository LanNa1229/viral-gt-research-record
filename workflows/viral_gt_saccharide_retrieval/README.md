# Step 3: reference-based Saccharide nearest-neighbor retrieval

Local v1.0.0 release candidate; final publication approval pending. Creator: Lan Na (given name Lan; family name Na).

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
Step 3 DOI/code commit: unassigned; add after explicit approval.

Proposed original-work licenses: MIT software, CC BY 4.0 derived data/documentation; third-party rights remain separate. See LICENSE_PROPOSAL.md and THIRD_PARTY_NOTICES.md.
