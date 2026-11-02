# Historical publication scope — November 2, 2026

This publication is authorized by Lan Na's explicit request to commit and push the complete historical Main Figure 4 record to `LanNa1229/viral-gt-research-record` on main. It introduces only this workflow directory. No repository-root files or unrelated work are included.

All 346 validated recovered originals are present at the same relative paths and match `SHA256SUMS-originals.txt`. The additional source inventory contains 1130 unchanged original files (4,980,245,796 bytes), located by saved raw-JSON and query-structure references in the historical hit tables. All 11 source-run `foldseek_all_hits.tsv` tables match their corresponding historical copies. No scientific transformation was performed during this curation.

The B2ZYD6 run output was relocated from `260617_Foldseek_batch/260618_foldseek_api_BAG41682_neighbors` to the same run-folder name under `260622_Foldseek_batch`. Its original input structures remain under `260617_Foldseek_batch/structures/`. Both actual source locations are recorded; original path strings are unchanged. The source run for A0A159B7X0 is historically spelled `A0A169B7X0_foldseek`; that spelling is retained, not silently corrected. Additional no-hit query models and saved empty-query responses are retained with their run folders.

`audit/source_supplement_provenance.json` records actual source paths, sizes, mtimes and SHA-256. The largest supplemental file is 88,306,599 bytes, below GitHub's 100 MiB file threshold; ordinary Git blobs preserve all files individually. No LFS pointer or ZIP-only substitute is used.

The completed manuscript audit is included under `validation/manuscript_comparison/`. Contact sheets of unrelated manuscript figures are excluded; the two actual heatmap images, source document hashes/XML, original-figure copies, complete 363-cell comparisons, decoded figure matrix and audit scripts/reports remain. A package-wide checksum covers the selected audit files; its former full-directory checksum remains in the untouched local audit.

Local proposed correction outputs and superseded release/Zenodo drafts are not selected for this historical-only commit. Their exact paths and reasons appear in `audit/publication_scope.json`. These are generated artifacts, not any of the 346 recovered original files. All remain available in the untouched Windows release candidate or as untracked local files. Prior recovery/validation reports describe their own execution-time state and are retained as historical records; the current README and this scope document govern the public payload.

Attribution: Lan Na — research conception, methodology, computational analysis, data curation, and visualization. The figure identification as Main Figure 4 is supplied by the user's publication request. Original code, input data and figures remain byte-preserved. No new biology, annotation correction, figure redesign, Foldseek search, scientific rerun or Zenodo operation is authorized or performed by this publication utility.

The repository has generic ignore rules for CIF models and log files. This authorized historical publication stages only the explicit checksum-verified WP08 path list with `git add -f --pathspec-from-file`; it does not change repository-wide ignore rules, omit source models/logs, or stage any unrelated path.
