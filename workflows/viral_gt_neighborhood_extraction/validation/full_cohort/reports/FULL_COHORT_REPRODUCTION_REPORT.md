# Full-cohort Step 1 reproduction: PASS for scientific extraction

## Conclusion

The full11-gene and15-gene genomic-neighborhood extraction is reproducible from the sealed prepared software and dataset archives in the documented Linux runtime. All757 windows at each size, all18,655 window-CDS rows, native locus memberships and genome units agree with the frozen records. No unexpected biological discrepancy, missing row, extra row, gene-order change or amino-acid sequence change was found.

This is stronger than the earlier representative-fixture result: a fresh full-cohort extraction has now actually run. The previous package/readiness assessment has **not** been modified. Exact byte-for-byte reproduction of every auxiliary historical table or execution receipt is not claimed: the adapter has reporting-schema omissions and relocated paths, detailed below. No BGC-Prophet inference or new candidate/taxonomy analysis occurred.

## Execution and isolation

- Fresh temporary workspace: `/tmp/viral_gt_step1_fullcohort_20261026_2v5b0ebd`.
- Unpacked software: `/tmp/viral_gt_step1_fullcohort_20261026_2v5b0ebd/code/workflows/viral_gt_neighborhood_extraction`.
- Unpacked dataset root: `/tmp/viral_gt_step1_fullcohort_20261026_2v5b0ebd/dataset/step1`.
- Fresh outputs: `/tmp/viral_gt_step1_fullcohort_20261026_2v5b0ebd/regenerated`; no existing output directory was overwritten.
- Input archives: software SHA256 `88f4034343bb47de7e2b81f27676631bd95b933fead0ae94525bef631f142b0a`; dataset SHA256 `d54e2897e22bd04b97925073148a1fed96e7db672b4f9f3ba205ffbd8c59e55f`.
- Interpreter: Python3.12.11/Biopython1.85 from the documented available environment. The package requirements supply Biopython and its transitive runtime dependencies; no environment/package was downloaded or installed.
- Python isolated mode `-I`, no user site and no project PYTHONPATH. The unchanged packaged CLI ran via a separate audit harness; exact process command and environment are recorded in `qc/execution_setup.json`.
- Full extraction plus structural grouping elapsed 19.25 seconds (archive unpacking/comparison excluded).

The harness enforced Python-level file-open/directory/SQLite boundaries. Reads outside unpacked packages were restricted to declared Python/OS runtime roots; writes were restricted to the newly generated output. A deliberately attempted read of the existing project AGENTS.md was blocked as a self-test. There were zero denied external-read attempts by the extraction itself. Model imports, subprocesses and network operations were prohibited and none was observed during extraction. This is an audit-hook test, **not kernel-level sandboxing or a claim of execution on a second physical machine**.

## Full scientific accounting

| Quantity | 15-gene fresh = frozen | 11-gene fresh = frozen |
|---|---:|---:|
| Original requested unique anchors | 594 | Same original set |
| Successfully mapped original anchors | 593 | Same mapped occurrences |
| Unresolved original anchor | UPI001FBC14AE | Preserved |
| Extracted occurrence windows | 757 | 757 |
| Complete windows | 616 | 658 |
| Shortened windows | 141 | 99 |
| Window-CDS rows | 10,688 | 7,967 |
| Eligible for later inference | 737 | 740 |
| Native connected-overlap loci | 565 | 575 |
| Genome units | 466 | 466 |
| Fixed baseline locus identities for paired analysis | 565 | 565 |

Fresh737/740 are **eligibility** counts; model inference was not run. They agree with the historical inferred-row counts. Window-CDS rows include repeated genes in overlapping windows and source aliases;18,655 is not a unique-protein count.

## Exhaustive comparisons

`tables/window_comparison.tsv` contains1,514 individual window results. Every shared field other than source_path is exactly equal, including source accession and nucleotide hash, anchor accession/gene ID, anchor/CDS indexes, strand, coordinates, source CDS count, actual gene count, missing offsets and translation-based eligibility. No size-specific window ID changed.

`tables/CDS_comparison.tsv` contains18,655 individual CDS results. Every field is identical, including accession, gene ID, source feature/index, relative position, order, coordinates, strand, compound location/parts, partial flag, product, amino-acid length/sequence/hash, and translation status/detail. All18,655 stored protein sequences are nonempty and exactly identical, so each has100% amino-acid identity. This matches the provided source translations and does not reclassify invalid nucleotide-to-protein verification cases.

Both regenerated `viral_gt_window_genes.tsv` files and both `viral_gt_window_proteins.faa` files are **byte-identical** to the frozen files. The15-gene extraction_qc.tsv is also byte-identical. FASTA headers, record membership and full protein sequences match exactly.

All565 native15 loci,575 native11 loci and466 genome units at both sizes match on every shared structural/membership field. Genomic segments, locus IDs, genome IDs, source-accession sets, anchor-protein sets, window membership, counts and grouping definitions match. The757 size pairs match all six fields emitted by the adapter. There are519,999 shared-table-cell comparisons:517,211 are exactly equal and2,788 are the verified CSV path relocations described below. Additional run-provenance comparison adds one input-path relocation and one wall-time difference.

## All differences and omissions

1. **Expected path relocation:**2,789 value differences:757 source paths in15 windows,757 in11 windows,757 in15 anchor mappings,517 GenBank manifest paths and one extraction input path. Each source-file relocation is checked against the packaged relocation map and the preserved hash. These are not accession, annotation or sequence changes.
2. **Runtime provenance:**one different `wall_seconds` value. Scientific extraction parameters, input/code hashes, parser policy and recorded Python version agree. Fresh source manifests have the same cache flags as the frozen baseline; no cache-flag discrepancy occurred.
3. **Fresh anchor-map reporting gap:**`retrieval_status` and `retrieval_reason` are not appended to the fresh extraction mapping. The original launcher added these enrichment fields after extraction. All758 original input/occurrence mapping rows (757 mapped plus one unresolved) match on shared fields; original enriched metadata remains in the packaged frozen records and status inputs.
4. **Size-pair map reporting gap:**nine original columns are not emitted in the minimal regenerated pair map: anchor_gene_id, anchor_gene_index, anchor_protein_id, eligible11, eligible15, extraction_verified_from_source_CDS_index, genome_accession, n_removed and source_sha256. The corresponding core identities/eligibility/source hashes are independently checked in the fresh window/CDS records; the historical verification flag is not newly manufactured. The six emitted pairing fields match all757 rows.
5. **No classifier re-emission:**four historical model-status fields (classification_status, fully_analyzed, n_analyzed_windows and n_positive_windows) are intentionally absent from each native-locus/genome summary at both sizes:16 table-column omissions. No inference was performed, so the adapter correctly emits structural summaries without inventing new model statuses.

Together these are27 schema omissions across auxiliary summary/mapping tables, **zero missing scientific fields in either window or CDS table**, and zero unexpected changed shared values. The complete lists are `tables/value_differences.tsv` and `tables/schema_differences.tsv`. `tables/frozen_output_file_inventory.tsv` also lists all30 packaged historical output files;14 auxiliary execution/validation/inference receipts or input-convenience copies are not re-emitted by the extraction-only adapter. They remain unchanged in the archive. New per-run context is recorded in portable_execution.json rather than pretending to recreate old timestamps or validation history.

## Undocumented-local-dependency check

Actual scientific reads were exactly1,114 packaged files:517 GenBank annotations,593 raw alias-evidence files, the anchor list, two portable manifests and the relationship crosswalk. Code reads were confined to packaged software; remaining reads were documented Python/Biopython/transitive runtime files. See `tables/observed_file_dependencies.tsv` and the raw `logs/filesystem_audit.jsonl`.

**No frozen expected window/CDS output table or test expectation was read during extraction.** Comparison against those archived outputs occurred only afterwards in a separate audit script. Thus the matching gene tables/FASTAs were regenerated from source annotations, not copied as expected results.

The packaged relationship crosswalk is an explicit documented dependency for preserving fixed baseline locus links in11-gene window rows. Those fixed labels are reused intentionally; native15/11 loci and genomes are computed fresh from regenerated coordinates/segments and independently match the frozen groups. The package does not claim to rediscover source-genome mappings by fresh online search; original accession crosswalks and exact annotation versions are declared inputs.

This establishes that the complete **scientific extraction** requires no undocumented local research file in the tested environment. The documented interpreter/Biopython/NumPy and OS libraries remain software dependencies external to the data archive. It does not establish exact original environment-lock replay, kernel-isolated execution or portability to every platform/version.

## Preservation and retained results

All1,830 release-package file hashes,1,773 source-copy identities, both sealed archives and their unpacked contents were checked unchanged. Original baseline output/snapshot manifests and finalized242-genome Sankey hashes remain intact. Nothing was published or altered in the release package, historical files or scientific candidate set.

Temporary outputs and the new regenerable annotation cache are retained at `/tmp/viral_gt_step1_fullcohort_20261026_2v5b0ebd/regenerated` for review; no unrelated cleanup was performed. Their sizes/hashes are in `tables/regenerated_output_manifest.tsv`. The separate audit scripts, execution evidence and all comparison results are saved in `/home/rahnn/260908_MiBG/analysis/viral_gt_step1_full_reproduction_20261026` and sealed by SHA256SUMS.tsv.
