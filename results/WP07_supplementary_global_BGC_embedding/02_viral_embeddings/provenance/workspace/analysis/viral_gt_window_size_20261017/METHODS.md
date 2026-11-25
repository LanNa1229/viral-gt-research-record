# Methods and reproducibility: paired 15- versus 11-gene inputs

## Frozen sources and separate outputs

Data Point 1 is `analysis/viral_gt_mining_20261017/runs/baseline_594_v2`; it is read-only. Data Point 2 is `analysis/viral_gt_mining_20261017/runs/window11_594_v1`. Shortened existing controls are in `analysis/viral_gt_mining_20261017/runs/controls11_existing_v1`. All interpretation, statistics, plots and new protein-cache writes occur in `analysis/viral_gt_window_size_20261017`.

The original 594-accession input is the frozen baseline cluster represented by A0A3Q8Q3U2, from headerless membership column 2 of the original coarse-clustering file, whose SHA256 is `cd731e42c1bffd7946427540d4cb61d72f706f494a93df2a8be11831f554fd2f`. The baseline mapped 593/594 accessions to 757 occurrences; the unresolved accession UPI001FBC14AE is not newly mapped or added. The input cohort table, original anchor crosswalk and existing taxonomy are retained with the new run. No spike-in membership or additional cluster is used.

At the beginning and end, verify all 676 frozen baseline files, 97 prior interpretation artifacts and 47 prior retrieval-rescue artifacts against their existing SHA256 inventories. Cross-platform manifest separators are normalized only for locating the same files. Original file contents, scores, embeddings and taxonomy are never rewritten.

## Occurrence-preserving extraction

Read each of the 515 source annotation files actually used by baseline windows and verify its stored SHA256. Use the frozen original parser, including every annotated CDS, regardless of strand, translation availability or pseudogene status. Sort by 1-based start, end and original feature index as before. Verify the complete archived 15-gene context, anchor gene ID, protein ID, nucleotide hash and CDS indexes against this independently parsed source annotation. Every retained sequence, accession, coordinate, location and product field remains unchanged.

Select indexes `max(0,anchor_index−5)` through `min(n_CDS,anchor_index+6)−1`. Retain boundary windows at their actual length, document absent relative positions, and do not wrap circular records, shift the window, repeat genes or pad positions. Stable new IDs hash source accession, genome nucleotide hash, the identical anchor gene ID, flank size 5 and the original no-wrap policy. A bijective 757-row old/new mapping is exported. Overlapping windows remain separate rows for inference; sequence-identical contexts may reuse the frozen pipeline's content-keyed inference cache without dropping observations.

Eligibility is recomputed on the retained genes using the same verified/full-or-partial translation and viral-lineage rules. An invalid outer gene can therefore cease to exclude a window; the three such newly eligible inputs are separate from the 737 paired observations. Every removed outer gene is exported with original CDS index, relative position, coordinates, sequence hash, sequence and product. Full 15-gene inputs remove four genes at positions −7, −6, +6 and +7; incomplete inputs can remove fewer, including none.

## Exact frozen inference and compatibility gate

Import the immutable baseline implementation snapshot. Before full inference, check model, encoder, wrapper, reference and regression-proof hashes through its original compatibility gate. Recompute four existing baseline inputs spanning short/complete and positive/negative cases, including ARF08732.1. Require all seven scores and every 320D coordinate to match within absolute tolerance 10⁻⁶, label identity to match, and all-actual-position and membership-pooled vectors to be exactly equal. Model tensors must remain unchanged.

The classifier checkpoint SHA256 is `f45ca298980d4740dd43131a7e69c99cb233368c815470572cfc46f33abf4eb5`; the ESM2 checkpoint SHA256 is `46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928`. Reuse exact sequence-hash-keyed per-protein vectors from the existing primary/control caches and audited historical pilots. Verify cache settings, vector shape, finite float32 values and equality of overlapping cache entries, then copy required vectors into a new cache. Only uncached required proteins are encoded; there were 33 for target inputs and zero for controls. Input-length-dependent Transformer representations are newly computed, with no changes to encoder/classifier weights.

ESM2 is the existing t6/8M model, layer 6 residue mean excluding BOS/EOS, first 1,022 residues for long proteins, CPU float32 and the existing batch-size-1 convention. Feed actual n≤11 genes in source-CDS order to the original Transformer, with an all-one float32 additive membership mask and no padding or repetition. Pool the actual contextualized positions by the original all-position mean immediately before the classifier. The source computes a positional-encoding result but passes the original input into the encoder; this behavior is retained. The unchanged classifier labels use threshold=max(Other score,0.5), remove Other when another label qualifies, and permit Unknown. Scores are independent sigmoid outputs and need not sum to one.

For inference outputs, record all seven scores, actual input membership, pooled 320D vectors, model identities, source hashes and content cache keys. Check hooks do not change scores, the pooled vector matches the classifier input, all values are finite, and weights/checkpoint files are unchanged.

## Reference geometry and paired quantities

Use the same 2,502 biological CORE_ONLY BGC vectors and metadata: vector SHA256 `78036e1636eaeeab90f4d157794c868d9b57033195939d83960d8e70bb831e22`, metadata SHA256 `5184f8f3199da2d1fac90e85ff11e2a30b443b70aac60992136718856c18ba4c`. There is one vector per biological BGC, not separate augmentation rows. Exactly 195 references carry Saccharide in their multilabel field.

Normalize saved query/reference vectors in float64. Rank by cosine, breaking ties by reference accession. P@5/P@10/P@20 count Saccharide-labeled BGCs among the k neighbors. Export their identities and cosines for both conditions. Also record nearest-Saccharide cosine/identity, mean cosine to all 195 Saccharide references, cosine to their normalized mean unit-vector centroid, and nearest-Saccharide minus nearest-non-Saccharide cosine. Compare top-k Saccharide sets via gained/lost identities and Jaccard; empty/empty sets are undefined. Record the between-condition cosine and top-20 full-neighbor Jaccard for every paired query.

Only 737 old/new occurrences with valid outputs in both conditions enter delta statistics and transitions. All deltas are 11−15. Score direction uses absolute tolerance 10⁻⁶; retrieval fractions use 10⁻¹² tolerance for arithmetic subtraction. Raw discrete retrieval-delta frequencies, score histograms, quantiles, min/max, sign counts, binary transitions and all per-window values are supplied. Three eligibility rescues and 17 exclusions are separate tables.

## Grouped statistics and candidate units

Primary continuous uncertainty: collapse identical baseline ordered contexts within an exact genome unit, average paired changes within each fixed baseline locus, then average loci within each genome. Resample whole genome means 5,000 times (seed 20261017), reporting the 2.5th and 97.5th percentiles. The complete 15→11 subset isolates removal of four genes; incomplete contexts and unchanged contexts are reported separately. Sensitivities remove generic taxid 2100421 or average genome means within the existing k21 sketch components before resampling components. These components are not species/ANI/phylogenetic lineages.

A=classifier positive; B=classifier negative with P@10≥0.5; C=A∪B. This user-specified retrieval cutoff is exploratory, not biologically calibrated. At entity level, any qualifying member establishes candidate status. Without a candidate, any uninferred member leaves the entity unresolved; otherwise it is observed nonselected. B and A may overlap at locus/genome levels despite disjoint window sets.

Primary entity comparisons keep the same 565 baseline regions and 466 exact-nucleotide/taxid genome units. This avoids counting fragmentation as recovered loci. Separately recompute native 11-gene connected-overlap components using the same baseline algorithm, giving 575 regions, and map each to exactly one frozen baseline parent. Candidate gain/loss registries distinguish observed noncandidate→candidate, prior-unresolved→candidate, candidate→observed noncandidate and candidate→unresolved. Paired-only counts restrict score evidence to the same 737 occurrences while preserving missing members and the original entity universe. Full-run counts additionally include newly eligible observations.

Taxonomy is reused, not refreshed. For recorded and conservative views, report exact analyzed denominators, candidate and unresolved counts for each family/genus. The conservative view withholds inherited ranks for generic source labels matching `uncultured|unclassified|environmental|metagenom|\bsp\.` case-insensitively. Candidates and genomes remain in explicit generic/unresolved bins. All frequencies are conditional on the original GT-anchor-selected mapped cohort.

## Background and removed-gene tests

Reuse the same 1,552 existing control inputs (1,252 random, 239 source-annotation GT/property matched and 61 stricter GT-proxy/property matched). Keep the original center and retain relative positions −5..+5. Previously shorter matched controls stay at their available input length, documented as existing short control inputs; no new genes or random centers are introduced. Full random controls are exactly 15→11 genes. Reverify all source annotation membership and sequences. Run the identical frozen inference/reference pipeline.

Primary background comparison uses the 600 complete target pairs across 356 genomes and their same-genome random controls. Average target changes by unique context→fixed locus→genome, and random controls within the same genome. Difference-in-differences is `(target11−target15)−(control11−control15)`. Resample whole genome contrasts 5,000 times. Nonoverlap sensitivity keeps the original 15-gene nonoverlap flag, a conservative fixed filter; it covers 311 genomes. An all-paired target sensitivity covers 420 random-control genomes. Component-group sensitivity reuses existing nucleotide-sketch groups.

Recheck the original 61 strict pairs at 11 genes before control inference: equal retained gene count, equal GT-proxy count≥1, AA totals within 25%, equal >1,022-residue protein count, reverse-strand count difference≤2 and original nonoverlap. Only 26 pairs across 20 genomes still match. Analyze those separately from all original 61 pairs; never infer specificity from assuming old constraints still hold. Existing source-annotation-only matched controls remain in the raw paired control table, with their known annotation limitations.

GT proxy is source-product GT text OR exact sequence identity to an original cohort anchor. Carbohydrate annotation and unknown/hypothetical annotation are prespecified product-text regexes recorded verbatim in `background_dilution_summary.json`; they are not new or experimentally validated functional assignments. For complete target pairs, examine removal of any GT proxy, any carbohydrate annotation, or ≥3 unknown-annotation genes. Compute within-stratum genome/locus-balanced mean changes and present-minus-absent contrasts. Jointly resample the union of genomes so a genome represented in both strata stays linked. These unadjusted, pointwise exploratory associations do not isolate causality, substrate, glycan production or biological relevance of a removed gene.

## Runtime, commands and verification

The local WSL environment is `analysis/bgcprophet_sequence_recovery_20260912/.venv/bin/python` (Python 3.12.11, NumPy 1.26.4, PyTorch 2.5.1+cpu, fair-esm 2.0.0, Biopython 1.85). No packages are installed. From the repository root, use that Python with `-B` for:

1. `analysis/viral_gt_window_size_20261017/scripts/experiment.py execute`
2. `analysis/viral_gt_window_size_20261017/scripts/paired.py`
3. `analysis/viral_gt_window_size_20261017/scripts/background_dilution.py`
4. `analysis/viral_gt_window_size_20261017/scripts/validate.py`
5. `analysis/viral_gt_window_size_20261017/scripts/report.py`
6. After figure QA: `analysis/viral_gt_window_size_20261017/scripts/validate.py --manifests-only`

The figure script uses the existing bundled Windows Python, `C:\Users\archival-user\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe`, prepending `analysis/bgcprophet_reference/.deps`. Run it with `-B analysis/viral_gt_window_size_20261017/scripts/figures.py`. The existing font cache is copied into this analysis's cache. All PNG/PDF/SVG exports use the same figure objects; inspect rendered PNGs and record fixes. No remote rendering or publication is performed.

Independent audit checks every retained gene against the archived baseline by CDS index, all source hashes, removed/missing membership, eligibility, all old/new pair deltas and entity outcomes, taxonomy partitions, background arithmetic and unambiguous TSV headers. It independently computes all 2,292 new score vectors using only the checkpoint's linear sigmoid head and the saved embeddings (no Transformer rerun), recomputes all new retrieval values, and verifies saved pooled vectors against cache tensors. Maximum head error is 9.43×10⁻⁸. The audit passed 181,614 checks. Complete output manifests include run caches; a separate protein-cache inventory covers every new cache file. Self-manifests and plotting caches are excluded from recursive artifact inventories.
