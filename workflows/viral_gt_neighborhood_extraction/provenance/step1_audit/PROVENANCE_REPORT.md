# Step 1 provenance and accounting audit

## Scope and authority

This is a forensic read-only audit of the original experiment, not a new extraction pipeline. Only new audit files were written. No extraction function, source parser from the pipeline, BGC-Prophet, ESM encoder, structure prediction, Foldseek, taxonomy inference or original unit test was run. Independent audit code reads existing tables, checks coordinate/index arithmetic and protein hashes, recomputes overlap-component membership from the saved gene locations, and verifies archived source bytes. Candidate selection and the final 242-genome Sankey were not changed.

The original 245/242 candidate sets are downstream selected subsets. They do not define the 594-anchor extraction input. The seven-anchor pilot is separate; it never expands this cohort. The latest corrected primary extraction is baseline_594_v2, not the conservative baseline_594_v1.

## Q1: original anchors

The authoritative original membership is `/home/rahnn/projects/viral-gt-research-record/results/WP01B_coarse_clustering/canonical/cluster_membership.tsv`. It is a headerless two-column table: column 1 is final cluster representative, column 2 member accession. Independent counting verifies 1,727 rows, 761 representatives and exactly 594 unique members for A0A3Q8Q3U2, the largest cluster, including its self-member. The next-stage anchor file is a sorted list of these exact members; there are no duplicate input lines and no candidate-based filtering before extraction.

The WP01B README links this to WP01A's MGAT2/GT-A query against BFVD and the 1,727 retained viral structure entries. It explicitly separates positive-control spike-in analyses. This establishes a structure-guided candidate cohort, not confirmed enzymatic activity for all 594 proteins. The cohort contains 592 UniProt accessions plus UniParc UPI001E761E12 and UPI001FBC14AE; the latter remains unresolved.

| Identity anchor | SHA256 |
|---|---|
| Original cluster membership | cd731e42c1bffd7946427540d4cb61d72f706f494a93df2a8be11831f554fd2f |
| Actual/frozen 594-accession text list | b1d316856a7e60dab0aa6e2414b829e90647f63d316104d2b8cef332e765d928 |
| Frozen anchor TSV | db46efba53ae4374883139862059c0f11c7c168d3afe38f55e3dc58461dc68ab |
| Exact frozen 15-gene extract.py | 9880dab746f4762ba4f7d8aacacb452feec3c0fc2384fcd71f84d42f8d3768e9 |

The actual-run `extraction_provenance.json` binds the input-list hash and extraction-module hash. `run_snapshot_checksums.tsv` binds the frozen code and input copies. The audit verified these hashes against current source bytes, and independently checked exact equality of membership sets. The source list was not inferred from current candidate taxonomy.

Mapping uses explicit UniProt/UniSave/UniParc genomic crosslinks, bounded retrieved annotations and exact full protein hashes. Unversioned source protein accessions may match all versions/occurrences present in supplied annotations; versioned identifiers match exactly. Archived sequence versions can differ from original BFVD entry versions; this audit does not assert universal historical BFVD sequence-version identity. Mapping was not an exhaustive public-genome/IPG search.

All 594 requested rows survive the mapping audit: 593 mapped, one unresolved. UPI001FBC14AE's saved reason is a PDB-only exact sequence with a terminal LEHHHHHH expression tag and no authoritative genomic crosslink. No tag stripping or native-sequence substitution was performed. It contributes no window or known genome denominator. A supplemental, separately documented rescue recovered eight other original inputs, including seven anchors whose proteins occur in historical nucleotide versions .2 but not newer .3 annotations; it did not add anchors to the 594 list.

## Q2: exact implementation and execution evidence

The exact executed 15-gene module is:

`/home/rahnn/260908_MiBG/analysis/viral_gt_mining_20261017/runs/baseline_594_v2/implementation_snapshot/src/viral_gt_mining/extract.py`

Its full SHA256 is `9880dab746f4762ba4f7d8aacacb452feec3c0fc2384fcd71f84d42f8d3768e9`.

The frozen launcher is `.../implementation_snapshot/workflows/viral_gt_mining/run_full.py`, SHA256 `4abf4b65b3c51a62b12c898e7f1129e7fbe51a8d7e5fc8a01363331a4fc18246`. It imports `viral_gt_mining.extract.extract` and calls it directly with the anchor text file, full/supplemental GenBank paths and combined explicit aliases. The launcher first copies implementation/input snapshots and checks their identity. It subsequently runs inference, taxonomy, grouping and validation; it is not a Step-1-only launcher.

The historical workflow README documents this command (not executed during this audit):

```sh
export PYTHONPATH=analysis/viral_gt_mining_20261017/checkout/src
analysis/bgcprophet_sequence_recovery_20260912/.venv/bin/python -u -B \
  analysis/viral_gt_mining_20261017/checkout/workflows/viral_gt_mining/run_full.py \
  --root /home/rahnn/260908_MiBG \
  --task /home/rahnn/260908_MiBG/analysis/viral_gt_mining_20261017 \
  --run-name baseline_594_v2
```

The inputs were `inputs/baseline_594_anchors.txt`, `inputs/full_genbanks.tsv` plus `inputs/supplement_genbanks.tsv`, and `inputs/retrieved_aliases.tsv` plus `inputs/supplement_aliases.tsv`. Every alias's evidence-path SHA is checked. The frozen copies in `input_snapshot/` are the recommended archival authority; the source records themselves live in the paths listed by `source_manifest.tsv`.

Git evidence: the associated repository is LanNa1229/viral-gt-research-record. The analysis clone's origin points to `/home/rahnn/projects/viral-gt-research-record`, branch `dev/viral-gt-windows-20261017`, HEAD/base `707f1379770c1c31ebeb5e07baa8551fac7338e1`. The code/tests/workflows are untracked in that clone, consistent with the historical git_receipt.json. Therefore this is a verified base commit, **not a commit containing the implemented extraction code**. Exact code identity is established by the frozen hash. A literal original terminal argv log was not located among the inspected execution artifacts; the documented command, exact launcher, code/input hashes and completed output records jointly establish execution provenance without proving the literal shell spelling.

The 15-gene code requires Python 3.12.11 and Biopython 1.85 plus standard-library csv/json/hashlib/sqlite3/pathlib/etc. The full historical environment additionally records NumPy 1.26.4, Torch 2.5.1+cpu and fair-esm 2.0.0; Torch/ESM are not required by the 15-gene extraction module. No complete original dependency lockfile was located. The current audit interpreter is separate from the original extraction interpreter and does not redefine its recorded environment.

### Verified algorithm

- Every annotated CDS contributes to the coordinate index, including missing/invalid translations and pseudogenes.
- Order is ascending start, end and original feature index; anchor orientation does not reverse the genomic list.
- The anchor is exactly relative offset zero. For 15 genes, code selects `genes[max(0,anchor_index-7):anchor_index+8]`.
- Missing offsets are explicit. No circular wrapping, opposite-side shifting, padding, repetition, inter-contig join or functional-annotation-based neighbor removal occurs.
- Full source sequences are retained. Translation verification uses source nucleotide CDS, genetic code and codon_start; explicitly fuzzy ends permit a separate exact partial-CDS verification path. Unverified/mismatched translations remain gene rows but prevent inference for that window.
- Window identity binds source accession, nucleotide hash, anchor annotation occurrence and policy. Input duplicates would map to the same window, but this run has none. Repeated lookup targets for one source occurrence are deduplicated. All 757 mapped rows correspond one-to-one to 757 distinct window IDs here.
- Overlapping anchors and GenBank/RefSeq counterpart occurrences remain distinct extracted windows. Only later grouping merges overlap loci and sequence/taxid-equivalent genome aliases.
- A compound CDS crossing a circular origin is preserved as annotated; downstream locus grouping uses its separate nucleotide segments rather than filling the intervening genome. This does not wrap the window around the source-CDS list.

### Later 11-gene extraction

The exact code is `/home/rahnn/260908_MiBG/analysis/viral_gt_window_size_20261017/scripts/experiment.py`, SHA256 `8d6fc8848b720de15b859699bbe4685506d56bdabcce00a3fdf73201a6e182a8`, with `scripts/common.py`. `prepare_main()` verifies original 15-gene membership against source annotations and calls `shortened()` to retain original CDS indexes within ±5. It checks this against `allgenes[max(0,idx-5):idx+6]` and preserves sequence, source accession/hash and anchor identity. It produces new size-specific window IDs plus a bijective 757-row pairing table. It does not merge overlaps.

The METHODS record documents `python -B analysis/viral_gt_window_size_20261017/scripts/experiment.py execute` using the same original recovery environment. `execute` includes model/control work; `prepare` is the extraction stage only. This audit executed neither. The current script hash matches the original extraction_provenance.json. Unlike 15-gene code, it is not in the baseline implementation snapshot; preserve the current hash-matched historical file in a future release.

**Minimal-package limitation:** even the original 11-gene `prepare` calls `verify_prior()`, which requires every artifact named by baseline_files.tsv (676), prior interpretation checksums (97) and retrieval-rescue checksums (47). `common.py` imports the frozen extract, inference and analysis modules plus NumPy, and uses existing baseline genome/locus identity tables. The full original launcher cannot be replayed from only the small Step-1 inputs without those preservation-gate dependencies. This audit identifies that limitation and does not implement a new adapter or remove safeguards. A future release must either supply the full original dependency closure or, with review, add and validate a clearly distinguished portable Step-1-only launcher. Frozen 11-gene outputs, source-validation rows and exact subset checks remain independently auditable now.

## Q3: verified accounting

| Stage / quantity | 15-gene baseline | 11-gene paired experiment |
|---|---:|---:|
| Requested input accessions | 594 | Same original cohort |
| Mapped original accessions | 593 | Same mapped occurrences |
| Distinct extracted occurrence windows | 757 | 757 |
| Complete nominal-size windows | 616 | 658 |
| Shortened windows | 141 | 99 |
| Annotated window-CDS rows | 10,688 | 7,967 |
| Stored inferred windows (later stage) | 737 | 740 |
| Extracted but not inferred | 20 | 17 |
| Native connected-overlap loci | 565 | 575 |
| Baseline loci used for paired comparisons | 565 | 565 fixed identities |
| Source accession versions with windows | 515 | Same 515 |
| Operational genome units | 466 | Same 466 |

474 mapped anchors have one occurrence; 119 have more than one. The complete occurrence multiplicities are 0:1, 1:474, 2:101, 3:6, 4:4, 5:5, 6:2 and 10:1. Thus 593 mapped anchors yield 164 additional occurrences, totaling 757. Different source protein IDs number 757 but their anchor protein sequences number 593; accession/occurrence counts must not be conflated with unique sequence counts.

There are 757 + 757 = **1,514 size-specific extraction rows**, paired into 757 anchor/source occurrences. They are not 1,188 independent genomic neighborhoods. Exactly 23 paired contexts have identical ordered extracted sequence lists despite different nominal sizes; this is an extraction-content fact, not an assertion that every pair was inferred or independently robust.

Of the 565 baseline locus components, 435 contain one window and 130 contain multiple windows (322 windows total). The audit detects 311 overlapping pairs under the original inclusive segmented-interval rule. The difference 757−565=192 is the component reduction, not a count of discarded anchors or independent duplicates. All windows remain in the original extraction/inference tables.

Genomes are grouped by exact whole-source nucleotide SHA256 **and source taxid**; missing identifiers retain accession-specific units. 515 source accession versions reduce to 466 units, with 49 extra counterpart accessions collapsed across 49 multi-accession units. There are only 464 distinct nucleotide hashes: the operational definition keeps sequence-identical records with different source taxids separate. These are not 466 independently verified complete viral genomes, species or evolutionary lineages. No taxonomy was revised in this audit.

The source manifest contains 517 successfully indexed files, but only 515 supplied extracted windows. KY523104.3.gb and MF405918.3.gb are the unused newer versions; supplemental .2 records supplied the exact historical anchor proteins. All 517 source checksums match (195,677,944 bytes total). Aliases depend on 593 distinct frozen evidence files (13,211,629 bytes); these are essential provenance inputs despite living under cache paths.

The 20 uninferred 15-gene windows remain extracted: ten have translation_mismatch and ten translation_unverified. The 11-gene run excludes seven mismatch and ten unverified contexts; three original mismatch cases become eligible when offending outer genes lie outside ±5. It has 737 shared inferred pairs, three newly eligible cases and 17 excluded in both. See `tables/noninferred_and_rescued_windows.tsv` for every identity and reason. Shortening alone does not establish improved biology, and short windows were not automatically excluded.

## Evidence limits and publication decision

`discrepancy_report.tsv` labels claims VERIFIED / PARTIALLY VERIFIED / UNVERIFIED. Verified counts and source hashes are distinguished from unresolved literal shell provenance, uncommitted code, absent exact environment lock, structural-function uncertainty, bounded genomic-mapping scope, portable path relocation and legacy 11-gene closure. See `PUBLICATION_STRUCTURE.md` for the small proposed GitHub code record and Zenodo scientific-data archive. No upload, commit, normalization, relocation, deletion or scientific-file edit has occurred.

