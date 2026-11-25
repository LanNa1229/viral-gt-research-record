# Methods record

## Frozen reference input and community membership

The unit is a biological MIBiG accession. The BGC-Prophet reference metadata fixes the 2,502 accessions, their row order and their original single- or multi-label class. The global joint Leiden grid contributes the saved k20_r1_s42 assignment. The [reference community crosswalk](inputs/community_assignments/REFERENCE_COMMUNITY_CROSSWALK.tsv) projects these fields from the master. The independent verification joins the public metadata and grid by accession, checks every matched JSON record digest and derives the Saccharide cohort from the original class tokens.

The [published embedding workflow](https://github.com/LanNa1229/viral-gt-research-record/tree/9bca2d7fa1e97b0c0f98d3c3b5b6c9c76070ead3/results/WP07_supplementary_global_BGC_embedding) documents the upstream embeddings and community inference. The original, frozen reference metadata and vectors are byte-identical at an earlier public revision and the global assignment grid is byte-identical at the release revision. Their direct paths, SHA-256 hashes and sizes appear in [PUBLIC_DEPENDENCIES.tsv](provenance/PUBLIC_DEPENDENCIES.tsv).

## Source extraction and curation

MIBiG 3.1 JSON contains the source cluster, product metadata, class annotations, chemical structures where reported and literature identifiers. The GenBank parser reads the original record and feature annotations and maps source gene identifiers by exact aliases. It retains unresolved mapping and source disagreement. One-to-many compound and gene data remain separate from the one-row-per-accession master. The [historical curation methods](provenance/historical_reference/CURATION_METHODS.md) document source retrieval, exact-identifier joins and evidence scope.

Detailed curation covers the 229-accession union of the Saccharide cohort with GC08 and GC05 reference members. The remaining 2,273 references retain automatically extracted source fields and a NOT_REVIEWED status. The source rules and accession-specific overrides in [the executed classification script](scripts/chemotype_annotation/curate_priority_records.py) encode the decision hierarchy. Direct biochemical and structural evidence applies to the tested enzyme or compound. Genetic evidence supports the reported association. Source records and predictions are recorded as such. Conflicts, provisional record-level decisions and unreviewed annotations remain explicit.

Literature and source-record interpretation was AI-assisted. The relevant task specification, frozen decisions, prior rule versions, evidence cards, bibliographic identifiers and curated feature records preserve the inspectable evidence. The historical record supports neither individual human validation of these decisions nor recovery of the preliminary manual accession ledger. Citation identifiers replace redistributing the 473 retrieved PubMed abstracts. The historical batch URLs and digests are retained because current PubMed responses can change. [NCBI data usage policy](https://www.ncbi.nlm.nih.gov/home/about/policies/).

## Display and historical definitions

The plotted cohort comprises any original class containing a Saccharide token, including hybrids. Each accession contributes once to its primary display category and to exactly one of GC08, GC05 and the other eight communities. The seven-group map pools LPS, LOS and O-antigen, groups S-layer with other cell-surface glycans, and assigns the remaining defined minor and unresolved primary groups to Other/unresolved. The bar denominator is the number of Saccharide references in that partition. The operational surface/extracellular union uses display categories 1–4.

The preliminary GC08 manual totals, subsequent legacy source reconstruction and broader source flag are preserved in [the historical comparison](outputs/community_composition/historical/GC08_HISTORICAL_COUNT_AUDIT.json). The figure's 72/100 differs in category boundaries from the legacy four-group reconstruction, which also totals 72/100. The broader source flag totals 74/100 because it includes two phosphonoglycans. Historical non-GC08 summaries include GC05. The figure's Other GCs excludes both GC08 and GC05.

## Figure and release treatment

Figure 2 is a 183 × 95 mm, 100% stacked horizontal bar chart using the frozen palette and category order. Figure exports are historical PDF, editable SVG and 600 dpi PNG. The [figure Methods](provenance/historical_figure/METHODS_RECORD.md) and [validation records](outputs/validation/figure/) retain the original settings and QC. The archival utilities perform integrity checks, independent composition calculation and scratch-directory replay. They preserve the original scientific results.

Workstation paths, private repository references and contact strings embedded in source comments were redacted in the relevant public copies. The [lineage table](provenance/DATA_LINEAGE.tsv) gives original and public file hashes and release treatment. The official MIBiG archives are external, checksum-pinned inputs. The separate historical viral-neighbor linkage branch is outside the figure workflow. [MIBiG release and attribution](https://zenodo.org/records/13367811).
