# Methods: exploratory retrieval rescue

This analysis uses existing frozen predictions and embeddings only. It does not run the ESM encoder, BGC-Prophet classifier, structure predictors, or training. The primary screening result remains `analysis/viral_gt_mining_20261017/runs/baseline_594_v2`. Its 676-file SHA256 inventory and the 97-file inventory of the prior interpretation analysis are verified before and after reanalysis. Outputs are written exclusively to this new directory, apart from the plotting library's temporary font cache if unavailable locally.

## Input population and representation

The original cohort is the baseline 594-member GT-A cluster represented by A0A3Q8Q3U2. We do not expand or redefine its anchors. The existing 757 windows comprise 737 inferred windows and 20 excluded windows. Frozen locus membership merges connected, overlapping inclusive nucleotide segments within an exact-genome/taxid unit; mere adjacency does not merge. There are 565 such loci and 466 genome units; exact genome sequence plus taxid grouping already collapses 49 GenBank/RefSeq aliases. Among these, 556 loci and 459 genomes contain at least one inferred window; 554 loci and 456 genomes have inference for every constituent window. Missing outputs are never imputed as negatives.

The 320D query vectors are the saved means of contextualized actual input genes under CORE_ONLY, using the original all-position pooling and all-one additive mask. All actual n≤15 genes are included without padding or repetition. Their ESM2 provenance is layer 6 residue mean excluding BOS/EOS, truncation at 1,022 residues, CPU float32, as documented in the frozen run. The current analysis reads the pooled arrays and scores; it does not invoke that pipeline again.

References are the same 2,502 biological BGC CORE_ONLY vectors, one row per base BGC, with 195 carrying the Saccharide label. A reference is counted as Saccharide if its whitespace-separated multilabel field contains that exact label. There are no augmentation rows counted separately. The reference vector SHA256 is `78036e1636eaeeab90f4d157794c868d9b57033195939d83960d8e70bb831e22`; metadata SHA256 is `5184f8f3199da2d1fac90e85ff11e2a30b443b70aac60992136718856c18ba4c`.

## Retrieval and exploratory sets

Normalize the saved vectors in float64, calculate cosine similarity, and sort by descending similarity with ascending BGC accession as deterministic tie break. P@k is the fraction of the first k reference BGCs carrying Saccharide, for k=5,10,20. Recomputed values must exactly match the existing fractions to tolerance 1e-12. Export the first 20 neighbors for every query to permit independent verification. Viral query windows have no reference BGC accession identity; references are not duplicated augmentations. This does not establish biological independence or remove all remote/shared-protein homology.

The classifier-positive field is taken directly from the frozen `saccharide_predicted` field, never relabeled by retrieval. “Classifier negative” includes all 491 inferred windows without that field, regardless of Other, Unknown or other class labels. It excludes the 20 uninferred windows.

The original published label rule uses threshold=max(Other score,0.5), selects scores at or above this per-window threshold and removes Other when another class also qualifies. Thus five frozen Saccharide-negative windows have Saccharide score≥0.5 but a larger Other score. The existing labels, their thresholds and Other scores are preserved; a uniform 0.5 score guide would not define the original decision boundary.

- A: frozen Saccharide-positive windows.
- B(t): observed classifier-negative windows with at least 10t Saccharide references in the first ten neighbors, t=0.4,0.5,0.6.
- C(t): A∪B(t).

The cutoffs are user-specified exploratory screening criteria. None is selected using a presumed biological truth panel. Counts at P@10≥0.3 and ≥0.8 are additional descriptive summaries. Other-k fraction cutoffs use the ceiling of k×cutoff: P@5≥0.5 means at least 3/5, not 2.5 neighbors. These are not equivalent resolutions across k. Enrichment=P@10/(195/2502) is descriptive and not a test statistic or probability of glycan biosynthesis.

At locus/genome levels, evidence is the union over constituent windows. An entity is `candidate` if at least one member window is in the relevant set. If none is and any member lacks inference, the status is `unresolved`. Otherwise it is `not_selected`. Consequently B and A are disjoint at window level but may overlap at locus/genome level. Report B counts and B\A separately; the latter is the true gain in distinct entities. A retrieval candidate can resolve a previously unresolved entity without inferring anything about its missing windows. Exact totals, any-inferred, fully-inferred and outcome-resolved denominators are retained separately.

## Taxonomy and redundancy

Reuse the same recorded taxonomy cache and genome metadata as the prior interpretation analysis. There is no taxonomic database refresh. For each cutoff, rank and taxonomy view, every genome is assigned to exactly one recorded taxon or explicit unresolved group. Export the full selected denominator and A/B/C candidate counts for both family and genus, including locus counts and missing-output counts. New taxa require at least one B genome and no A genome in that rank assignment. A family/genus taxid rather than spelling alone defines a lineage for breadth counts.

The conservative view withholds inherited ranks for source labels matching the prior rule `uncultured|unclassified|environmental|metagenom|\bsp\.` (case-insensitive). This sensitivity does not change the original NCBI taxonomy or delete candidates. Such genomes stay in the explicit “Generic source; rank withheld” group and in the selected denominator. A non-generic source label does not itself certify a species identification.

Overlap sensitivity collapses the prior context signature (exact genome unit and ordered protein sequence hashes with genomic positions/strands). Locus and genome identities remain frozen. Complete-window sensitivity restricts to n=15; incomplete windows have n<15. Subsets can share loci/genomes and their entity counts must not be added.

Genome-relatedness sensitivity reuses the existing canonical k21 bottom-1,024 SplitMix64 sketch Jaccard matrix components at cutoffs 0.2,0.5,0.8. For each component count candidate presence and additional components in C absent from A. Also report components containing additional genomes, which may already contain another A genome. These are descriptive operational redundancy groups, not species delimitations, ANI or a phylogeny. No independent-window confidence intervals are assigned to exact counts in this finite selected cohort.

## Cluster 7

Verify original pilot files against its archived checksums, then select only the `CORE_ONLY/primary` row from `cluster7_embeddings.npy`: the actual 14 ARF08726–ARF08739 core genes. Match ESM/classifier checkpoint identities to baseline metadata, all-one mask, layer/truncation, original pooling and passing historical inference QC. Recompute cosine retrieval against the same CORE_ONLY reference vectors and verify against the archived comparison table. Other preprocessing modes remain separate.

Check actual cohort overlap through baseline ordered protein accessions, source genome accession and anchor mapping. Export any overlapping baseline window and its own frozen scores/retrieval. The pilot is a distinct input definition and is not appended to baseline counts, even if its genome or core proteins occur in the baseline. It is an illustrative candidate rather than a validated biological positive or a sensitivity-calibration panel.

## Runtime, reproduction and verification

The scripts use Python's standard library plus the existing repository NumPy/Matplotlib dependencies. This turn used the bundled Windows Python at `C:\Users\archival-user\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe`, prepending the existing `analysis/bgcprophet_reference/.deps` directory. All scripts set `sys.dont_write_bytecode=True` and the plotting script directs Matplotlib caching to this analysis's cache directory. No packages were installed.

From the repository root, run that Python with `-B` and each of:

1. `analysis/viral_gt_retrieval_rescue_20261017/scripts/analyze.py`
2. `analysis/viral_gt_retrieval_rescue_20261017/scripts/figures.py`
3. `analysis/viral_gt_retrieval_rescue_20261017/scripts/validate_and_report.py`

The independent accounting audit rebuilds A/B/C from frozen localization and prediction fields, verifies every entity status and aggregate count, checks taxonomy partitions and explicit missingness, checks all exported neighbor P@k calculations, exercises partially missing status logic, and verifies baseline/prior hashes. It also checks figure file signatures. Visual checks of rendered PNG exports are documented separately; PDF/SVG exports use the same Matplotlib figure objects. `artifact_checksums.tsv` covers deliverables and scripts, excluding caches and itself.

## Interpretation limits

Candidate sets show association with Saccharide references in the model's supervised representation. The score and retrieval are correlated measurements from that representation. There is no verified disjoint held-out panel establishing a biological cutoff or classifier false-negative rate. Prior GT/property-matched controls did not establish added specificity beyond GT content. A multi-enzyme glycan module, pathway direction, product, substrate, autonomous biosynthesis and actual glycan production remain unverified. All taxonomic frequencies are conditional on the original GT-anchor-selected cohort and mapped, eligible windows.
