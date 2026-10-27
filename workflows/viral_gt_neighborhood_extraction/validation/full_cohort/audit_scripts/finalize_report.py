"""Write final full-cohort reproduction conclusions and seal only the new audit outputs."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import csv,json,hashlib,collections,datetime
D=Path(__file__).resolve().parents[1];ROOT=D.parents[1];REL=ROOT/'analysis/viral_gt_step1_release_20261026';TMP=Path((D/'temporary_workspace.txt').read_text().strip());DATA=TMP/'dataset/step1';OUT=TMP/'regenerated'
def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,rs,fields):
    with Path(p).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rs)
def main():
    compare=json.loads((D/'qc/comparison_result.json').read_text());dep=json.loads((D/'qc/dependency_validation.json').read_text());guard=json.loads((D/'qc/filesystem_guard_summary.json').read_text());run=json.loads((D/'qc/portable_execution.json').read_text());setup=json.loads((D/'qc/execution_setup.json').read_text());assert compare['scientific_extraction_reproduced'] and dep['passed']
    preserved=[]
    for name in ['release_manifest.tsv','source_copy_manifest.tsv']:
        rr=read(REL/name)
        for r in rr:
            p=REL/r['file'];assert sha(p)==r['sha256'],p
            if name=='source_copy_manifest.tsv':assert sha(Path(r['original_path']))==r['sha256'],r['original_path']
        preserved.append(dict(manifest=str(REL/name),files=len(rr),status='PASS'))
    for r in read(REL/'REVIEW_SHA256SUMS.tsv'):assert sha(REL/r['path'])==r['sha256']
    for r in setup['archives']:assert sha(Path(r['path']))==r['sha256']
    for path,digest in json.loads((D/'qc/unpacked_inputs_before.json').read_text()).items():assert sha(Path(path))==digest
    for base in [ROOT/'analysis/viral_gt_mining_20261017/runs/baseline_594_v2',ROOT/'analysis/viral_taxonomy_reconciliation_20261019/figures/final_sankey_242']:
        for name in (['artifact_checksums.tsv','run_snapshot_checksums.tsv'] if base.name=='baseline_594_v2' else ['SHA256SUMS.tsv']):
            rr=read(base/name)
            for r in rr:assert sha(base/r['path'])==r['sha256'],(base,r['path'])
            preserved.append(dict(manifest=str(base/name),files=len(rr),status='PASS'))
    # Index all fresh outputs, including the explicitly regenerable index cache.
    fresh=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),role='regenerable annotation cache, not a biological output' if 'regenerable_cache' in p.parts else 'fresh full-cohort extraction output') for p in sorted(OUT.rglob('*')) if p.is_file()]
    write(D/'tables/regenerated_output_manifest.tsv',fresh,list(fresh[0]))
    runtime=collections.Counter()
    for r in guard['file_accesses']:
        if r['category']!='documented_python_runtime':continue
        s=r['path'];runtime[s.split('/site-packages/')[1].split('/')[0] if '/site-packages/' in s else 'Python standard-library/runtime']+=1
    sources=read(D/'tables/observed_file_dependencies.tsv');assert not [r for r in sources if r['category']=='packaged_data' and r['logical_path'].startswith('outputs/')]
    file_inventory=read(D/'tables/frozen_output_file_inventory.tsv');not_emitted=[r for r in file_inventory if r['regenerated']=='False']
    tablecells=sum(x['compared_cells'] for x in compare['tables']);exactcells=sum(x['equal_cells'] for x in compare['tables'])
    result=dict(passed=True,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),verdict='Full scientific extraction reproduced from the prepared package in the documented runtime; not a byte-for-byte clone of every historical auxiliary report',temporary_workspace=str(TMP),elapsed_seconds=guard['elapsed_seconds'],windows_compared=1514,CDS_rows_compared=18655,all_window_nonpath_fields_match=True,all_CDS_fields_and_order_match=True,AA_sequences_identical=18655,common_table_cells_compared=tablecells,exact_table_cells=exactcells,scientific_discrepancies=0,value_difference_categories=compare['reported_value_differences_by_category'],schema_omissions=27,unregenerated_historical_auxiliary_files=len(not_emitted),no_undocumented_scientific_file_reads_observed=True,frozen_expected_window_or_gene_table_reads_during_extraction=0,documented_runtime_files_by_package=dict(runtime),release_archives_and_package_unchanged=True,historical_inputs_and_outputs_unchanged=True,preservation_checks=preserved,model_inference_calls=0,publication=False,limits=['Python-level audit hooks, not kernel container/syscall isolation','Existing documented Python/Biopython/NumPy and OS runtime remain external software dependencies','Packaged fixed-baseline identity crosswalk is an explicit input; native loci independently regenerated','Auxiliary metadata columns and execution-history receipts are not all re-emitted; no claim of exact original environmental or shell replay'])
    (D/'qc/final_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    report=f'''# Full-cohort Step 1 reproduction: PASS for scientific extraction

## Conclusion

The full11-gene and15-gene genomic-neighborhood extraction is reproducible from the sealed prepared software and dataset archives in the documented Linux runtime. All757 windows at each size, all18,655 window-CDS rows, native locus memberships and genome units agree with the frozen records. No unexpected biological discrepancy, missing row, extra row, gene-order change or amino-acid sequence change was found.

This is stronger than the earlier representative-fixture result: a fresh full-cohort extraction has now actually run. The previous package/readiness assessment has **not** been modified. Exact byte-for-byte reproduction of every auxiliary historical table or execution receipt is not claimed: the adapter has reporting-schema omissions and relocated paths, detailed below. No BGC-Prophet inference or new candidate/taxonomy analysis occurred.

## Execution and isolation

- Fresh temporary workspace: `{TMP}`.
- Unpacked software: `{TMP}/code/workflows/viral_gt_neighborhood_extraction`.
- Unpacked dataset root: `{DATA}`.
- Fresh outputs: `{OUT}`; no existing output directory was overwritten.
- Input archives: software SHA256 `88f4034343bb47de7e2b81f27676631bd95b933fead0ae94525bef631f142b0a`; dataset SHA256 `d54e2897e22bd04b97925073148a1fed96e7db672b4f9f3ba205ffbd8c59e55f`.
- Interpreter: Python3.12.11/Biopython1.85 from the documented available environment. The package requirements supply Biopython and its transitive runtime dependencies; no environment/package was downloaded or installed.
- Python isolated mode `-I`, no user site and no project PYTHONPATH. The unchanged packaged CLI ran via a separate audit harness; exact process command and environment are recorded in `qc/execution_setup.json`.
- Full extraction plus structural grouping elapsed {guard['elapsed_seconds']:.2f} seconds (archive unpacking/comparison excluded).

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

Together these are27 schema omissions across auxiliary summary/mapping tables, **zero missing scientific fields in either window or CDS table**, and zero unexpected changed shared values. The complete lists are `tables/value_differences.tsv` and `tables/schema_differences.tsv`. `tables/frozen_output_file_inventory.tsv` also lists all{len(file_inventory)} packaged historical output files;{len(not_emitted)} auxiliary execution/validation/inference receipts or input-convenience copies are not re-emitted by the extraction-only adapter. They remain unchanged in the archive. New per-run context is recorded in portable_execution.json rather than pretending to recreate old timestamps or validation history.

## Undocumented-local-dependency check

Actual scientific reads were exactly1,114 packaged files:517 GenBank annotations,593 raw alias-evidence files, the anchor list, two portable manifests and the relationship crosswalk. Code reads were confined to packaged software; remaining reads were documented Python/Biopython/transitive runtime files. See `tables/observed_file_dependencies.tsv` and the raw `logs/filesystem_audit.jsonl`.

**No frozen expected window/CDS output table or test expectation was read during extraction.** Comparison against those archived outputs occurred only afterwards in a separate audit script. Thus the matching gene tables/FASTAs were regenerated from source annotations, not copied as expected results.

The packaged relationship crosswalk is an explicit documented dependency for preserving fixed baseline locus links in11-gene window rows. Those fixed labels are reused intentionally; native15/11 loci and genomes are computed fresh from regenerated coordinates/segments and independently match the frozen groups. The package does not claim to rediscover source-genome mappings by fresh online search; original accession crosswalks and exact annotation versions are declared inputs.

This establishes that the complete **scientific extraction** requires no undocumented local research file in the tested environment. The documented interpreter/Biopython/NumPy and OS libraries remain software dependencies external to the data archive. It does not establish exact original environment-lock replay, kernel-isolated execution or portability to every platform/version.

## Preservation and retained results

All1,830 release-package file hashes,1,773 source-copy identities, both sealed archives and their unpacked contents were checked unchanged. Original baseline output/snapshot manifests and finalized242-genome Sankey hashes remain intact. Nothing was published or altered in the release package, historical files or scientific candidate set.

Temporary outputs and the new regenerable annotation cache are retained at `{OUT}` for review; no unrelated cleanup was performed. Their sizes/hashes are in `tables/regenerated_output_manifest.tsv`. The separate audit scripts, execution evidence and all comparison results are saved in `{D}` and sealed by SHA256SUMS.tsv.
'''
    (D/'reports/FULL_COHORT_REPRODUCTION_REPORT.md').write_text(report)
    (D/'README.md').write_text('''# Fresh full-cohort Step1 reproduction

Read [the scientific report](reports/FULL_COHORT_REPRODUCTION_REPORT.md), [per-window comparison](tables/window_comparison.tsv), [per-CDS comparison](tables/CDS_comparison.tsv), [all value differences](tables/value_differences.tsv), [schema omissions](tables/schema_differences.tsv), and [observed dependencies](tables/observed_file_dependencies.tsv). Machine-readable final result:qc/final_validation.json.

Fresh extraction ran from unpacked sealed archives in a separate /tmp directory, with isolated Python and file-access guards. Zero scientific differences were found; all auxiliary differences and omissions are reported explicitly. The release package and original files were not changed. No model inference or publication occurred.

Reproducible audit scripts are under scripts/:run_full_reproduction.py (creates a new temporary directory),compare_outputs.py,audit_dependencies.py,finalize_report.py. guarded_execution.py is a separate audit harness and does not change the packaged CLI. Existing temporary location is recorded in temporary_workspace.txt. Runtime stdout/stderr and file-access trace are under logs/. No new pipeline implementation was created.
''')
    (D/'qc/preservation.json').write_text(json.dumps(dict(passed=True,preservation_checks=preserved,release_review_manifest_unchanged=True,sealed_archives_unchanged=True,unpacked_inputs_unchanged=True),indent=2)+'\n')
    files=sorted(p for p in D.rglob('*') if p.is_file() and p.name not in ['SHA256SUMS.tsv','audit_receipt.json'] and '__pycache__' not in p.parts)
    rr=[dict(path=p.relative_to(D).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in files];write(D/'SHA256SUMS.tsv',rr,list(rr[0]))
    for r in rr:assert sha(D/r['path'])==r['sha256']
    (D/'audit_receipt.json').write_text(json.dumps(dict(passed=True,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=len(rr),bytes=sum(r['bytes'] for r in rr),sha256_manifest=sha(D/'SHA256SUMS.tsv'),temporary_workspace=str(TMP),scientific_discrepancies=0,release_unchanged=True),indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
