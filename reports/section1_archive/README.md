# Section 1 archive and reproducibility record

Scope: Figure 1 and Supplementary Figures 1–4, work packages WP01A–WP01F as defined in work_package_crosswalk.tsv. This record distinguishes ORIGINAL (historical source supported by execution provenance), ARCHIVED (surviving source without exact historical execution proof), RECONSTRUCTED (new implementation/plot assembly), and AUDIT-ONLY (verification utilities).

The initial 96 numerical checks are reused. A new portable frozen-input validation passed106 checks, and isolated historical postprocessors reproduced the network tables,63-cluster coherence summary and12-control-comparison table byte-for-byte. All30 historical statistical rows agree numerically; their serialized missing-value forms differ. No structural search/alignment/model prediction was rerun.

Original scientific code is never overwritten by reconstructed utilities. Some historical scripts contain author-local absolute paths; do not execute them against the archive directories. Use an isolated working copy and explicit inputs. Commands actually tested in this candidate are recorded separately from historical commands.

The baseline script's SHA is historically anchored; its exact May2025 Foldseek executable build is not established. Original search exports and the final PanelC merge/classifier remain unavailable. Truncated-core files are nonessential historical evidence under the current author scope decision. Current-script/output agreement proves present reproducibility, not a historical execution date or original authorship.


## Recovered L230 sources and original spike input availability

The observed POS001 input directory has 1,732 surviving paths: 1,727 baseline PDBs and five control CIFs. Its two L230 CIF copies are absent. Eight historical AF3 model/input/summary records were recovered elsewhere, including both named L230 model sources. These are ARCHIVED source candidates; their byte identity with the absent POS001 copies is not established. The surviving apo model agrees byte-for-byte with the separate direct-control apo copy. The recovered model and apo have 895 matching residue identities but different coordinates (fixed-correspondence rigid-body CA RMSD about 27.293 Angstrom); do not call these recovered candidates structurally identical or silently substitute them as proven original spike copies.

The retained original spike Foldseek database is therefore included as essential analytical evidence rather than an optional cache. Its 14 components total about 3.66 MB; its lookup contains all 1,736 entries and its source map names 1,734 input files, including both L230 identifiers. The historical 7/7 recovery and 483-entry membership claims remain verified independently of the missing CIF-copy lineage. No new clustering was run. The restoration plan distinguishes the 1,732 verified original paths from the two ARCHIVED AF3 source candidates. The original-location manifest preserves 1,727 SHA-identical baseline aliases formerly present in the spike input, without duplicating payload bytes.


## Author scope decision, 24 October 2026

The downstream baseline consists of 1,727 structures retrieved through the original 447-residue MGAT2 Foldseek search. Its file identities and downstream results are independently verified. This provenance is supported by the author's explicit confirmation and the preserved historical records; it is not a newly recovered web-search export. The separate truncated GT-A core search was performed on another computer. Its unavailable files are a documented, non-blocking historical limitation, outside the current release-preparation requirements. No equivalence between the two search outputs is independently verified or claimed. No further search, reconstruction or rerun of the truncated-core workflow is planned.

Current priorities are figure/caption review, scoped GitHub validation, creator/license/redistribution decisions, and Zenodo submission packages. The truncated-core issue does not delay either package. Original manuscript files remain unchanged; precise proposed replacement wording is supplied for author adoption.
