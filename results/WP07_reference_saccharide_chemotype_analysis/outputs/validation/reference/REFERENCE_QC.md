# Reference chemotype QC and scientific summary

**PASS**: identity, joins, original labels/rows, priority statuses, normalized tables and saved viral crosswalk verified.

## Dataset accounting

- mibig_matches: 2502
- priority: 229
- priority_assigned: 225
- saccharide_assigned: 193
- priority_unresolved: 4
- background_NOT_REVIEWED: 2273
- priority_provisional: 2
- hybrid_reference_records: 467
- BGCs_with_recorded_structures: 1860
- BGCs_with_explicit_MIBiG_gene_annotations: 396
- raw_GBK_CDS: 41002
- normalized_gene_records: 41002
- genes_with_explicit_activity_assay: 378
- gene_mapping_ambiguous: 0
- compound_relationships: 4145
- viral_links: 5740
- original_vs_mibig_class_disagreements: 0

All 229 priority entries have an explicit curation status; 225 have a proposed category and 4 remain unresolved. 193/195 Saccharide references have a proposed category. NOT_REVIEWED background entries are not biological negatives. Provisional and unresolved entries require further review.

## Chemical composition

GC05 Saccharide references (n=24) are a subset of non-GC08 Saccharide references (n=95). These columns overlap and must not be treated as independent biological samples.

| Primary category | GC08 Saccharide (100) | GC05 Saccharide (24) | Non-GC08 Saccharide (95) |
| --- | ---: | ---: | ---: |
| Capsular polysaccharide | 31 (31.0%) | 4 (16.7%) | 5 (5.3%) |
| Extracellular/matrix polysaccharide | 19 (19.0%) | 2 (8.3%) | 3 (3.2%) |
| Lipopolysaccharide (LPS) | 8 (8.0%) | 1 (4.2%) | 1 (1.1%) |
| Lipooligosaccharide (LOS) | 2 (2.0%) | 0 (0.0%) | 0 (0.0%) |
| O-antigen | 8 (8.0%) | 2 (8.3%) | 3 (3.2%) |
| S-layer glycan | 3 (3.0%) | 0 (0.0%) | 0 (0.0%) |
| Other cell-surface glycan | 1 (1.0%) | 0 (0.0%) | 0 (0.0%) |
| Glycopeptidolipid | 1 (1.0%) | 0 (0.0%) | 0 (0.0%) |
| Teichuronic acid | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) |
| Phosphonoglycan | 2 (2.0%) | 0 (0.0%) | 0 (0.0%) |
| Aminoglycoside / aminocyclitol-related metabolite | 20 (20.0%) | 7 (29.2%) | 14 (14.7%) |
| Glycosylated specialized metabolite | 4 (4.0%) | 5 (20.8%) | 56 (58.9%) |
| Other carbohydrate-related product | 0 (0.0%) | 3 (12.5%) | 11 (11.6%) |
| Non-carbohydrate-associated product | 0 (0.0%) | 0 (0.0%) | 1 (1.1%) |
| Unresolved | 1 (1.0%) | 0 (0.0%) | 1 (1.1%) |

### All reference members, including non-Saccharide exceptions

| Primary category | GC08 all (105) | GC05 all (53) |
| --- | ---: | ---: |
| Capsular polysaccharide | 31 | 4 |
| Extracellular/matrix polysaccharide | 19 | 2 |
| Lipopolysaccharide (LPS) | 8 | 1 |
| Lipooligosaccharide (LOS) | 2 | 0 |
| O-antigen | 8 | 2 |
| S-layer glycan | 3 | 0 |
| Other cell-surface glycan | 1 | 0 |
| Glycopeptidolipid | 1 | 0 |
| Teichuronic acid | 1 | 0 |
| Phosphonoglycan | 2 | 0 |
| Aminoglycoside / aminocyclitol-related metabolite | 20 | 7 |
| Glycosylated specialized metabolite | 6 | 19 |
| Other carbohydrate-related product | 0 | 4 |
| Non-carbohydrate-associated product | 2 | 12 |
| Unresolved | 1 | 2 |

## GC08 historical audit

The four legacy product groups yield 72/100. Source-specific capsule precedence yields 31 CPS and 19 EPS rather than the earlier 29/21 split. BGC0000768 (gellan) and BGC0000797 (S88 sphingan) are explicitly capsular in their cited papers and can also be described as secreted EPS; moving those two to an EPS-first convention would give 29/21 without changing the 72-record legacy union. This numerical reconciliation does not establish the unavailable original accession-level decisions.

The broader source-supported surface/extracellular union is 74/100. The additional BGC0000806 and BGC0000807 are chemically characterized extracellular phosphonoglycans (PMID 24584498). They were outside the earlier four named bins. BGC0000769 glycopeptidolipid has an explicit chemical category but the inspected source did not establish a specific surface location; its surface flag remains NOT_DETERMINED.

## Ambiguities and annotation limitations

- BGC0000727: indigoidine in JSON versus an auricin deoxysugar-tailoring sequence description. Retained as Unresolved; original Saccharide label unchanged.
- BGC0001846: original Saccharide/NRP annotation versus documented valinomycin/montanastatin cyclodepsipeptides; chemistry and model label remain separate.
- BGC0001510: cryptic glycosylation in anisomycin biosynthesis; non-glycosylated final product retained.
- BGC0001747: pseudaminic-acid precursor for a cell-surface protein glycan, not a complete surface polymer.
- BGC0000780: mixed O/K antigen scope; BGC0000793: LPS core oligosaccharide, not automatically a complete O-antigen or LOS pathway.
- BGC0000763 cites an amylovoran paper while the pinned record describes stewartan; BGC0001412 retains historical PS B/C naming ambiguity.
- BGC0001765: cited salecan chemistry concerns ZX09 rather than the S33 sequence. Only provisional record-level chemotype is retained.
- BGC0002450: sarasinoside/glycosidation candidate locus lacks demonstrated BGC-product linkage in the inspected source; provisional.
- Product identity does not demonstrate every enzyme. Gene assay flags are restricted to explicit gene-level source assay evidence. Unknown chemistry/enzymes are not coded as absent.

## Reproducibility and scope

See CHEMOTYPE_DATA_DICTIONARY.md, CURATION_METHODS.md, validation_results.json, SOURCE_MANIFEST_SHA256.tsv, OUTPUT_MANIFEST_SHA256.tsv, execution logs and executed code snapshots. All models/clustering/nearest-neighbor calculations and enrichment tests were excluded. No previous archive or viral master was modified, and no commits/pushes were performed. Proposed private destination is recorded in CURATION_METHODS.md; explicit approval is required before archival push.
