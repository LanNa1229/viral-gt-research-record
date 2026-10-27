# Step1 viral GT genomic-neighborhood extraction — version1.0.0

**Scientific extraction fully reproduced. Final publication approval pending.** The exact historical code is preserved alongside a separate portable adapter. Original input/output files, candidate criteria and final taxonomy figures are unchanged. No BGC-Prophet inference was rerun.

## Scientific results

The original WP01B structural cluster represented by A0A3Q8Q3U2 has594 unique input accessions, including its representative.593 map; UPI001FBC14AE remains unresolved. It is a structure-guided candidate cohort, not the later245/242 selected genomes or proof of GT function for every anchor.

| Metric, fresh = frozen | 15-gene | 11-gene |
|---|---:|---:|
| Matching windows | 757/757 | 757/757 |
| Matching window-CDS rows | 10,688/10,688 | 7,967/7,967 |
| Complete / shortened windows | 616 /141 | 658 /99 |
| Native loci | 565 | 575 |
| Genome units | 466 | 466 |
| Eligible for later inference | 737 | 740 |

All accessions, anchor identities, coordinates, gene orders, completeness, translation statuses and amino-acid sequences match. Both CDS tables and both protein FASTAs are byte-identical. Fixed baseline loci for paired comparisons remain565. Genome units are exact nucleotide hash plus source taxid, not species or independently validated complete genomes. Overlapping window-CDS rows are not unique proteins.

## Evidence and remaining reporting differences

The full standalone package was unpacked into a fresh /tmp directory and completed full-cohort extraction using Python3.12.11/Biopython1.85 with no project PYTHONPATH. The file-access audit observed no undocumented scientific input reads and no reads of frozen expected window/CDS output tables during extraction. The packaged fixed-baseline identity crosswalk is an explicit input; native loci are regenerated. Audit hooks operate at Python level, not as a kernel container or a claim about every platform.

There were zero biological discrepancies. The documented differences are2,789 expected path relocations, one elapsed-time difference,27 omitted auxiliary columns and14 historical auxiliary files not re-emitted. The27 column omissions are16 inference-status fields across structural summaries, two post-extraction retrieval-enrichment fields, and nine redundant paired-map fields; all core window/CDS fields are complete. No model result is fabricated. These reporting/environment differences do not alter the extracted biological content. A complete historical environment lock and literal original shell transcript remain unavailable; the original whole directory is not claimed byte-identical.

Read [current assessment](provenance/REPRODUCIBILITY_ASSESSMENT.md), [full-cohort report](validation/full_cohort/reports/FULL_COHORT_REPRODUCTION_REPORT.md), [comparison summary](validation/full_cohort/tables/table_comparison_summary.tsv), [schema omissions](validation/full_cohort/tables/schema_differences.tsv), and [validation result](validation/validation.json). Full individual-window/CDS/difference and dependency tables accompany the Zenodo data. Older reports under provenance/step1_audit/ and fixture_validation.json are date-specific archival records, superseded as the current assessment by full-cohort evidence.

## Workflow and preserved code

All annotated CDSs count, ordered by genomic start, end, and source feature index. The anchor is relative0; flanks are +/-7 or +/-5 where available. No wrapping, shifting, padding, repetition or overlap merging occurs during extraction. Translation failures remain gene rows and determine later eligibility only. Overlapping anchors and accession aliases remain separate windows; structural grouping occurs afterwards.

- historical/baseline_594_v2/ preserves the exact original implementation snapshot.
- historical/window11/ preserves the original experiment.py/common.py.
- portable/step1.py is a new adapter that imports the hash-pinned original15 extractor and isolates unchanged11 shortening/structural grouping function ASTs. It does not import broad model launchers.
- All820 legacy11 guard dependencies were recovered and hash-verified. The adapter uses explicit input/code hash checks without distributing unrelated model/control artifacts; see [dependency recovery](provenance/DEPENDENCY_RECOVERY.md).
- tests/fixture/ contains eight exact source records and nine requested anchors, with nine frozen window expectations per size. Nine tests passed in two environments and from the standalone code ZIP.
- config.json records fixed historical policy; editing it does not override the hash-pinned algorithm.

The old base commit707f1379770c1c31ebeb5e07baa8551fac7338e1 did not contain the then-uncommitted extraction implementation. Exact frozen SHA256 identities establish its original version. Historical source bytes are never edited or silently relicensed. Current relative manifests are distinct from original absolute-path provenance.

## Tested execution instructions

Use Python>=3.11 and Biopython1.85. Full-cohort validation used Python3.12.11/Biopython1.85/NumPy1.26.4; fixture compatibility was also demonstrated with Python3.11.17/Biopython1.88. requirements.txt provides the minimal dependency target, not a retrospectively invented historical lock. No installation is required when these dependencies are already available.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -B tests/test_portable.py
DATA=/absolute/path/to/extracted_dataset/step1
python -B portable/step1.py verify --data-root "$DATA"
python -B portable/step1.py extract --data-root "$DATA"   --expected-anchors 594 --output /absolute/path/to/NEW_step1_outputs
```

The last command is the now-tested full-cohort extraction operation. It creates window15/,window11/,structural summaries and portable_execution.json; it does not run inference, remote retrieval or taxonomy/structure prediction. Output must be a new directory outside code/data packages. Inputs and historical scripts are hash-checked, and existing outputs cannot be overwritten. The original broad run_full.py/experiment.py execute launchers are archived for history and must not be used as extraction-only entry points. External audit scripts are archived as executed under validation/full_cohort/audit_scripts/ and retain historical local paths; they are evidence, not the portable CLI.

## Creator, licenses and publication references

Creator confirmed by the user:Lan Na, given name Lan, family name Na. Manuscript coauthors are not automatically dataset creators. Third-party database attribution is separate.

Proposed for final approval:MIT for Lan Na's original software and CC-BY-4.0 for original derived data/documentation. NCBI/INSDC records retain their source-policy terms and attribution; UniProt-family records retain upstream CC-BY-4.0 attribution. Embedded source sequences are not silently covered by the project's grant. See [file-level source/license summary](manifests/file_level_sources_licenses.tsv) and [license proposal](LICENSE_PROPOSAL.md).

Existing public project:https://github.com/LanNa1229/viral-gt-research-record . Target directory:workflows/viral_gt_neighborhood_extraction/. No Step1 push or public DOI has occurred at this checkpoint. Publish one scientific dataset linked in Zenodo metadata to the exact approved Step1 GitHub commit, never the old repository base. No duplicate software DOI or automatic whole-repository archiving is proposed. DOI, publication date and public commit are filled only after explicit approval. The immutable data archive can remain hash-stable by placing assigned public identifiers in record metadata and a post-publication reference record.
