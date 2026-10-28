# Scientific and technical limitations

## Community regeneration

The primary weighted graph and community assignments are supplied as the manuscript reference dataset. Both original analytical scripts execute in the pinned Python 3.11 environment. Classification from the supplied partition reproduces all 687 viral rows, including all protein-level fields and Figure 2B counts. Fresh Leiden constructs a byte-identical full graph but yields a different partition (ARI 0.9451745569665067). It retains 4,854 nodes, 224,825 edges and 580 isolates; direct family-support totals are 126/105/21 (all/concordant/alternative), compared with primary 130/110/20. Historical Python/igraph/leidenalg versions are unknown. Exact de novo regeneration of the historical partition has not been demonstrated. Keep fresh results separate; do not replace manuscript reference outputs or claim exact end-to-end replay.

## Structures and redistribution

All 4,854 original coordinate files are included unchanged and pass extraction, size, SHA256 and ID checks. Extraction requires case-sensitive storage because 27 pairs differ only by case; the helper checks this before writing files. The coordinate manifest records actual paths alongside legacy metadata filenames; use actual archived paths for Foldseek. Source-specific redistribution terms remain unestablished for 2,237 files. This is a public-deposition clearance item, not missing scientific data or demonstrated evidence of a prohibited license. See `RIGHTS_AND_RETRIEVAL.md`.

## Interpretation

Structural neighbourhood and community membership do not establish catalytic activity or formal CAZy family membership. The primary 130 direct supports comprise 110 concordant with the modal community family and 20 supporting another family. Figure 2B puts those 20 in its broader limited/ambiguous/alternative category. No curated known-viral-GT evidence table was supplied to the original classifier; its `unclassified_viral` evidence status must not be interpreted as proof that every protein is novel. Topological connectors are screening calls requiring alignment/domain validation for biological interpretations. Unresolved taxonomy remains unresolved.

Interactive figure styling is outside computational validation. Computational source tables are provided; historical styling records are not required to execute these analyses.
