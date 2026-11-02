# Foldseek batch result analysis

## Input overview
- Total parsed hit rows: **6152**
- Unique queries: **20**
- Unique clusters: **1**
- Databases: **afdb-swissprot, afdb50, cath50, pdb100**

## High-confidence rule
- probability >= 0.7
- evalue <= 0.001
- q_coverage >= 0.5
- aln_length >= 80
- High-confidence hits: **2507**

## Hit-name search
- Name search terms: glycosyltransferase, glucosyltransferase, glycosyl, glucosyl, galactosyltransferase, mannosyltransferase, xylosyltransferase, rhamnosyltransferase, fucosyltransferase, CAZy, GT-A, GT-B, TET, T-even, T4, hydroxymethyl, 5hmC, 5hmdU, JBP, DNA, methyltransferase, oxygenase, dioxygenase, lipopolysaccharide, LPS, O-antigen, capsular, polysaccharide, hypothetical, uncharacterized, putative
- Name-search hits: **1040**
- Queries with name-search hits: **20**

## Keyword category hit counts
- gt_like: 489 hits, 9 queries
- tet_like: 254 hits, 7 queries
- dna_modification: 250 hits, 17 queries
- lps_capsule_like: 79 hits, 7 queries
- enzyme_characterized_hint: 1959 hits, 19 queries
- hypothetical: 428 hits, 19 queries
- any_keyword: 2637 hits, 20 queries

## Priority database hits
- Priority databases: pdb100
- Priority-database rows: 1397

## Key output files
- `00_all_hits_annotated.tsv`
- `01_top_hits_per_query_database_named.tsv`
- `02_high_confidence_hits.tsv`
- `03_keyword_hits.tsv`
- `04_pdb_priority_hits.tsv`
- `05_query_database_summary.tsv`
- `06_query_overall_summary.tsv`
- `07_cluster_database_summary.tsv`
- `08_pairwise_top_target_overlap.tsv`
- `09_query_cluster_coherence.tsv`
- `10_database_top_hit_names.tsv`
- `11_name_search_hits.tsv`
- `12_name_search_summary_by_query.tsv`

## Interpretation notes
- Start with `10_database_top_hit_names.tsv` for quick readable top-hit inspection.
- Use `11_name_search_hits.tsv` to find hits whose names/descriptions contain terms such as glycosyltransferase, TET, 5hmC, DNA, LPS, or capsular polysaccharide.
- Use `04_pdb_priority_hits.tsv` for experimentally grounded structural homologs.
- Keyword/name hits are triage labels only. Final interpretation should include alignment coverage, active-site conservation, structural superposition, and gene-neighborhood evidence.