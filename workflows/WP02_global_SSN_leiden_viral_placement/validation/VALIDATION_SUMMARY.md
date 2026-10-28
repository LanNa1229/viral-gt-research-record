# Functional validation — 2026-10-28

Both unmodified scientific scripts execute successfully in the pinned Python 3.11 environment. Dependency pins match installed versions; `pip check` reports no broken requirements. CLI/import/syntax checks and independent graph fixtures pass.

| Metric | Preserved primary / classifier replay | Fresh execution |
| --- | ---: | ---: |
| Nodes | 4,854 | 4,854 |
| Full undirected edges | 224,825 | 224,825 |
| Isolates | 580 | 580 |
| Multi-protein communities | 80 | 80 |
| Viral proteins | 687 | 687 |
| CAZy-unassigned viral proteins | 383 | 383 |
| Direct family-supported | 130 | 126 |
| Dominant-family-concordant | 110 | 105 |
| Alternative family-supported | 20 | 21 |

All protein-level columns of the preserved classifier output match after numeric serialization tolerance (1e-10 relative, 1e-12 absolute); categorical fields match exactly. The generated Figure 2B count tables match the archived source tables. Fresh node/edge IDs, assignments, degree, strength and direct-family support pass independent graph-based checks. Its full edge table is byte-identical to primary. Fresh-vs-primary partition ARI is 0.9451745569665067; it is a functional new analysis, not exact historical regeneration. Resolution 1, seed 1, original weights and convergence settings were retained.

All 4,854 coordinate files passed a new full extraction on case-sensitive Linux storage with SHA256/size validation; no coordinate was renamed or modified. The original scientific scripts were not changed; no compatibility fix was required. No Foldseek, parameter sweep or interactive visualization audit was run.

The only unresolved data-deposition clearance item is source-specific ownership/redistribution terms for 2,237 coordinate files. AFDB attribution, experimental PDB terms and version-specific local AlphaFold 3 output notices are included. AF3 output sharing through scientific publication is permitted subject to its terms; it is not unrestricted CC BY 4.0. Missing source permission is not assumed to be a grant. Unknown historical dependency versions limit the exact-reexecution claim but do not invalidate the supplied primary dataset or its validated downstream replay.

Machine-readable tests are in `publication_functional_validation_20261028.json`, `publication_reference_validation_20261028.json` and `publication_coordinate_validation_20261028.json`. Test logs and fresh outputs are local validation records outside the reference dataset. No push/upload/publication has occurred.
