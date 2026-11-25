# Curation methods and evidence limits

## Source snapshot and immutable identity

MIBiG 3.1 JSON and GenBank archives already existed in the workspace and were reused without combining releases. JSON official URL: https://dl.secondarymetabolites.org/mibig/mibig_json_3.1.tar.gz . Original retrieval: 2026-09-08T05:00:37.926106+00:00; SHA-256 87221d7ead87709fbcfe81176d4f1684520017e1d30d5f1f24ef7a39310d1e94. GenBank official URL: https://dl.secondarymetabolites.org/mibig/mibig_gbk_3.1.tar.gz ; SHA-256 db3cb48d74c3143f38f9628f9e2abbfa994bfcff74d60b59a88ffe20ec3fa796. Its original download time was not recovered from the inspected manifest; this limitation is not replaced by filesystem timestamps. The reuse timestamp is recorded in SOURCE_PREFLIGHT.json. Raw tarballs and all raw JSON records are retained.

The exact private archive commit is 284d4f1655cf036c1890bdffcb4aa2dc93f6e41e. Remote Git blob identities were checked against local consumed clustering inputs. The 2,502 original reference accessions, zero-based embedding rows, labels, and core gene counts were preserved; the NPY header/whole-file hash were inspected for identity only. Existing k20_r1_s42 global assignments were selected by exact column name. The 287-locus master and 5,740 saved neighbor identities/ranks/cosines were joined by exact identifiers. No models, nearest-neighbor searches, clustering or UMAP fitting were run.

## Scope and annotation rules

The priority population is the 229-record unique union of all 195 original Saccharide references with all GC08 and GC05 references. Detailed chemotype decisions are independent of community labels and use the same rules throughout. All other 2,273 records receive automatic source extraction and explicit NOT_REVIEWED chemical status. Primary categories are mutually exclusive; secondary tags may overlap. Hybrid original labels are preserved as strings and JSON token lists.

The exact rules, accession-level overrides and explanations are in curate_priority_records.py and CURATION_DECISIONS.json. Obvious product terms and structured source subclasses receive RULE_BASED_SOURCE_SCREEN status; record/structure/literature interpretation receives explicit evidence-guided status. No human manual verification is claimed. No category is chosen to satisfy the 72/100 historical target. A capsule takes primary precedence over a broad extracellular-polysaccharide label when capsule is specifically supported; both descriptions may remain in secondary scope flags. Mixed O/K annotation is retained as Other cell-surface glycan rather than forced into a single O-antigen/LOS/CPS identity.

Gross glycosylated architecture can be curated from a recorded chemical structure. No structure was generated, no molecular mass recalculated, and no sugar stereochemistry, linkage, precursor or tailoring chemistry was invented from broad biosynthetic labels. Normalized names use Unicode NFC and whitespace only, preserving original names independently. Compound relationships are kept separate so congeners/aglycones do not multiply the BGC master.

## Evidence hierarchy

1. Explicit biochemical or chemical-structure experiments in the primary study support only the tested enzyme/compound claim.
2. Genetic or heterologous-expression evidence supports the tested BGC/gene involvement, not every catalytic function.
3. MIBiG product/structure/subclass annotations support record-level chemistry and are retained separately from direct experiment claims.
4. Recorded source genomic predictions, GBK product annotations and MIBiG sequence-based predictions remain predicted/unspecified evidence.
5. Missing, conflicting or unreviewed information is explicit. Product identity never makes all enzymes experimentally characterized.

PubMed metadata/abstracts for all 473 cited PMID records of the priority set were retrieved separately on 2026-10-27. Original source publication IDs are not replaced; supplemental DOI metadata remains a separate field. Exact request URLs, raw XML and hashes are retained. This is record/abstract-guided curation, not exhaustive full-text experimental replication. Supplemental publisher/author source notes are in SUPPLEMENTARY_SOURCE_NOTES.json.

JSON gene annotations and GT specificities are mapped to GBK CDS by unique exact protein/gene/locus aliases. Unmatched or ambiguous aliases remain separately documented, never merged by a fuzzy guess. The standard-library GenBank parser preserves qualifier text/locations and raw archive hashes; the normalized table omits protein translation text because lossless source archives retain it. Existing antiSMASH-related qualifiers in the frozen GBK are source predictions, not newly run analyses.

Chemical claims are accession-specific. Experimental salecan data for ZX09 are not copied to the S33 entry. Phosphonate-containing polysaccharides are not called phosphorylated glycans. A sugar precursor or cryptic glycosylation step is distinguished from a final glycan. Genus/species columns are explicitly lexical components of the recorded source organism; no new taxonomy retrieval was performed.

## Historical GC08 audit

The original accession-level preliminary manual ledger was not supplied or found in the inspected sources. Therefore this is an independent accession-level source reconstruction, not a claimed recovery of every old manual decision. The legacy four-group union is compared with the broader surface/extracellular union separately; overlapping definitions are not independent observations. No enrichment p-values or viral chemistry predictions were calculated.

## Source rights and archive destination

MIBiG identifies CC-BY-4.0 attribution licensing: https://mibig.secondarymetabolites.org/ . Cite the MIBiG contributors and source release. Retrieved publication abstracts may have separate third-party rights; raw literature files are local research-source records and redistribution requires checking applicable terms. No Git commit/push is authorized by completion of this task. Proposed private destination: REDACTED_PRIVATE_REPOSITORY/results/WP07_reference_chemotype_curation_20261027/ .
