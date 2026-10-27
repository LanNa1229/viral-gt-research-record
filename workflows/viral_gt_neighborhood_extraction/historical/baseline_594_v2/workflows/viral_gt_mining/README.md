# Viral GT-centered genomic mining

This workflow tests whether neighborhoods of a prespecified structure-derived
viral protein cohort receive Saccharide-associated computational representations.
Neither model classification nor reference proximity establishes glycan synthesis.

## Population and source authority

The full cohort is exactly the original WP01B baseline cluster represented by
`A0A3Q8Q3U2`: 594 distinct accessions including the representative. The headerless
`results/WP01B_coarse_clustering/canonical/cluster_membership.tsv` has representative
in column 1 and member in column 2. Its SHA256 is
`cd731e42c1bffd7946427540d4cb61d72f706f494a93df2a8be11831f554fd2f`.
The original WP01B README explicitly separates this baseline from subsequent
positive-control spike-in analyses. `prepare_inputs.py` checks all these counts and
checksums against the original research checkout and exports an immutable list.
The independently prespecified seven-anchor pilot is separate from the 594 list.

The original scratch workspace is `/home/rahnn/260908_MiBG`. The associated research
repository was clean on `main` at commit `707f1379770c1c31ebeb5e07baa8551fac7338e1`.
Development uses a separate local clone on `dev/viral-gt-windows-20261017` inside
`analysis/viral_gt_mining_20261017/checkout`. No original files are intentionally
modified, and no push, publication or new structure/Foldseek work is performed.

## Environment and execution

Reuse the installed Linux environment (Python 3.12.11, Biopython 1.85, NumPy 1.26.4,
PyTorch 2.5.1+cpu, fair-esm 2.0.0). CUDA is unavailable in this environment. No
training, checkpoint conversion or package installation is required.

From the original scratch workspace in WSL:

```sh
export PYTHONPATH=analysis/viral_gt_mining_20261017/checkout/src
PY=analysis/bgcprophet_sequence_recovery_20260912/.venv/bin/python
TASK=analysis/viral_gt_mining_20261017
ROOT=/home/rahnn/260908_MiBG
$PY -B $TASK/checkout/workflows/viral_gt_mining/prepare_inputs.py --root "$ROOT" --task "$ROOT/$TASK"
$PY -B -m unittest discover -s $TASK/checkout/tests/viral_gt_mining -v
$PY -B -m viral_gt_mining.extract --anchors $TASK/inputs/pilot_anchors.txt \
  --genbank-manifest $TASK/inputs/pilot_genbanks.tsv --run $TASK/runs/pilot_v1 --cache $TASK/cache
$PY -B -m viral_gt_mining.inference --root "$ROOT" --run $TASK/runs/pilot_v1 --localize \
  --sequence-cache-root $TASK/cache/esm
$PY -B $TASK/checkout/workflows/viral_gt_mining/validate_run.py --root "$ROOT" --run $TASK/runs/pilot_v1 --pilot
$PY -u -B -m viral_gt_mining.retrieve --root "$ROOT" --task "$ROOT/$TASK" --workers 4
$PY -u -B -m viral_gt_mining.supplement --root "$ROOT"
$PY -u -B $TASK/checkout/workflows/viral_gt_mining/run_full.py --root "$ROOT" --task "$ROOT/$TASK" --run-name baseline_594_v2
```

Full execution gates on the verified pilot and completed retrieval. Completed
source downloads, annotation indexes, protein representations and individual
window inferences are cached. Interrupted runs resume from these intermediates.
The bounded supplement resolves the seven exact archived Tupanvirus proteins via
their protein-record `coded_by` links to nucleotide versions `.2` (their proteins
are absent from current `.3` annotations). It also resolves the genomic UniParc
anchor via exact protein/genome crosslinks; the PDB-only tagged UniParc construct
remains unresolved. Original requested IDs remain anchors, with no new cohort
members. Supplemental manifests are separate from original retrieval manifests.
To preserve a previous completed result, choose `run_full.py --run-name NAME`.
Repeating the same run name rewrites only derived tables in that run. Do not run
two writers against the same run/cache simultaneously.

The Linux environment lacks matplotlib. Render scientific figures using the
existing bundled Windows Python and `render_windows.py --root UNC_ROOT --run
UNC_RUN`. This reuses the existing Windows packages under
`analysis/bgcprophet_reference/.deps`, without installing dependencies. It writes
300-dpi PNG and vector PDF files and must be followed by visual inspection.

## Mapping and extraction definitions

Prefer saved, verified source GenBank annotations. Remote mapping uses current
UniProt EMBL/RefSeq cross-references, with the latest available archived UniSave
entry for inactive/deleted accessions. Exact full protein-sequence hashes must
match the cross-referenced annotated CDS. All linked genomes and all qualifying
CDS occurrences are retained. This scope is explicit cross-references, not an
exhaustive public-genome sequence search or IPG occurrence census. Archived entry
versions need not be the original BFVD sequence version; their evidence paths,
sequences and retrieval checksums are retained for inspection.

UniProt requests use four workers and a shared domain limiter. NCBI genome
retrieval uses batches of five accessions, at most three concurrent workers and
at most one request start per 0.37 seconds (anonymous rate below 3/s). Transient
errors retry with backoff. Invalid batches split into individual accessions so a
single obsolete accession cannot discard the remaining genomes. Cache files are
atomic and verified by checksum; incomplete/error responses are not valid cache
entries. Failed anchors/genomes remain in explicit status tables.

Every annotated CDS counts, including pseudogenes and genes without usable
translations. CDSs are sorted by genomic start, end, then original feature index.
Windows are exactly index -7 through +7 where available. Coordinates are 1-based
inclusive. Upstream/downstream means increasing genomic coordinates, independent
of anchor strand; amino-acid sequences retain their original coding orientation.
No boundary compensation, inter-contig crossing, artificial genes or circular
wrapping occurs. Compound locations and strands are recorded explicitly; source
translations are checked with Biopython extraction and source genetic code.
Explicitly partial CDSs use a separate verification path: fuzzy biological 5'/3'
ends allow the documented missing start/stop, source codon_start is honored, and
an incomplete trailing nucleotide codon can be ignored only at a fuzzy 3' end.
No amino-acid residue is fabricated or altered; the translated sequence must
exactly match the preserved source protein. These rows are `verified_partial`.
Unverified/mismatched translations exclude the entire window from inference but
do not remove any annotated CDS from the neighborhood. The first conservative
pass is retained as `baseline_594_v1`; the corrected latest analysis uses a new
run name, `baseline_594_v2`, with frozen source-code and input snapshots.

Versioned protein accessions match exactly. Unversioned protein IDs match all
versions/occurrences in supplied annotations. Duplicate input rows remain in
`anchor_mapping.tsv` and link to the same stable window. Different anchor
occurrences and overlapping anchor windows remain separate inference observations.
Window IDs bind genome accession, nucleotide sequence hash, anchor annotation
occurrence and extraction policy. They are stable across input-list reordering
and restarts, and can change when source annotation changes.

## Frozen representation and scientific compatibility

Primary inputs are the actual window genes (up to 15), no padding or repetition.
The exact previously validated `benchmark.py` forward is imported without editing
historical scripts. It verifies classifier `f45ca298980d4740dd43131a7e69c99cb233368c815470572cfc46f33abf4eb5`
and ESM2 `46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928`.
ESM2-8M layer 6 uses residue means excluding BOS/EOS, first 1,022 residues, CPU
float32. Full proteins remain in the extraction FASTA; truncations are explicit
in QC. Unique sequence hashes cache vectors across overlapping windows/runs.

The original Transformer computes positional encoding but does not use it in its
encoder input. Float32 0/1 masks are additive key-attention biases, not exclusion
masks. For these CORE_ONLY inputs all actual positions have membership 1; mean
pooling over actual genes equals the all-input-position mean supplied to the
original classifier. Hooks verify this equality and preserve original scores.
The original seven sigmoid outputs are independent scores and need not sum to 1.
Labels use threshold `max(0.5, Other score)` and the original Other/Unknown rules.
The separate downstream detection reporting gate is bypassed for candidate loci.

The original 128-position construction, genuine context, published all-position
pooling and historical REPEAT_128 sensitivity are documented in
`../../../audit_inference.md` (task-root `audit_inference.md`). No new repetition
strategy is introduced. Direct 128-position viral comparison is unavailable for
short genomes and is deliberately omitted from the primary actual-gene workflow.

Retrieval uses ONLY the existing verified `(2502,320)` CORE_ONLY reference array
`analysis/viral_glycan_benchmark_20261016/mibig_CORE_ONLY_embeddings.npy` and its
row metadata. Startup validates model/source/checkpoint identities, the exact
reference checksums and historical encoder regression proofs before expensive
encoding. The historical encoder comparison shows numerical near-identity, not
bitwise identity, between original LMDB and current CPU encoding. Published
all-position/core-context arrays are incompatible with this isolated-core input
and are not substituted. The 2,502 biological references include 195 with any
Saccharide label; augmentations are not separate biological neighbors.

## Outputs and interpretation

Each versioned `runs/NAME` directory contains:

| File | Grain / definition |
|---|---|
| `anchor_mapping.tsv` | Every input row × matched occurrence; unresolved inputs have a status/reason |
| `viral_gt_windows.tsv` | Stable window ID, source accession/sequence hash, anchor, bounds, missing offsets, taxid, eligibility |
| `viral_gt_window_genes.tsv` | Window × actual annotated CDS, relative position, coordinates, strand, protein sequence/hash and verification |
| `viral_gt_window_proteins.faa` | Available source protein sequences; headers bind window, offset and gene occurrence |
| `extraction_qc.tsv` | Window completeness and inference exclusions |
| `viral_bgc_predictions.tsv` | One inferred window, seven sigmoid scores, original threshold/labels |
| `viral_bgc_embeddings.npy` | Finite float32 `(n_inferred_windows,320)` contextual embeddings |
| `embedding_manifest.tsv` | Exact array row, window ID, ordered genes/hashes, masks and checkpoint identities |
| `inference_qc.tsv` | Every extracted window, inferred or explicitly excluded |
| `saccharide_localization.tsv` | Original-320D cosine neighbors, Saccharide P@5/10/20, descriptive P@10 enrichment over 195/2502 |
| `candidate_loci.tsv` | Connected overlapping genomic intervals, only merged after inference for locus summaries |
| `genome_units.tsv` | Exact nucleotide-hash + source-taxid groups and all GenBank/RefSeq aliases |
| `viral_taxonomy_summary.tsv` | Family/genus/species × level; explicit unresolved groups, counts and denominators |
| `three_level_summary.tsv` | Window/locus/genome positive, fully analyzed negative, unknown and denominator counts |
| `taxonomy_cache.tsv` | NCBI source ranks and stable rank identifiers, raw-XML checksums |
| `validation.json`, `independent_validation.tsv` | Independent source/array/neighbor checks |
| `*_provenance.json`, source manifests | Parameters, runtime, resources, versions and input/output hashes |

Classification positivity is the original model's Saccharide label. Localization
P@k and enrichment are separate descriptive retrieval metrics; no binary
localization cutoff is selected after seeing results. PCA reuses the saved
reference-only CORE_ONLY fit and is a visualization, not a quantitative proximity
test. Scientific exports include reference/viral PCA overlay, fractions and viral
taxonomic distributions.

GenBank/RefSeq counterpart records with identical whole-genome nucleotide hash and
nonempty source taxid count as one genome. Missing taxid/hash keeps accession-
specific units. Loci use connected overlapping inclusive nucleotide intervals
within these genome units; adjacency alone does not merge. All source accessions
remain traceable and all window inference rows remain intact.
For source CDSs that themselves span a circular origin, locus calculations use
separate genomic segments rather than filling the unoccupied middle of the
genome. `candidate_loci.tsv` records these segments explicitly. Extraction never
wraps extra neighboring genes across a contig end.

At locus/genome level, any observed positive makes the unit positive; a negative
requires predictions for every constituent window. Otherwise its outcome is
unknown. Report resolved denominators, missing outcomes and total-population
bounds. Wilson intervals are descriptive: overlapping windows and selected,
dependent viral samples violate independent random-sampling assumptions. No
all-virus prevalence, taxonomic enrichment or experimental confirmation is inferred
from this preselected structural cohort. Unmapped anchors have no known genome
and are reported separately, never treated as negative genomes.

## Performance and integrity

The seven-window local pilot indexed 3 genomes in about 0.5 seconds. First
inference took 4.44 seconds, peak RSS 346 MiB, with 50 existing sequence-cache hits
and 2 new encodings; classifier calls totaled 0.076 seconds. The second run
reused all 7 window results in about 2 seconds. A `(594,320)` float32 output
matrix is only about 0.73 MiB, but multiple genomic mappings can exceed 594 rows.
Protein vectors require 1,280 bytes per unique sequence plus file overhead;
up to 8,910 raw positions before mapping multiplicity/overlap corresponds to
about 10.9 MiB in dense sequence vectors. Genome annotation downloads and cached
whole-genome JSON/SQLite, rather than arrays, dominate disk use.

Do not extrapolate two short pilot encodings to long viral proteins without
qualification. Full inference is CPU-only; four PyTorch threads are used. The
live unique-sequence count, encoded lengths, elapsed times and peak RSS provide
the completed-run benchmark. Retrieval latency depends on genome size and
availability. Keep sources/checkpoints read-only and use `validate_run.py` to
verify historical artifact hashes after execution. Automated tests cover ends,
strand/order, missing CDS proteins, duplicate/versioned accessions, multiple
contigs/occurrences, overlapping anchors, archived mapping, restart identity,
nonfinite caches, missing-outcome denominators and genome counterparts.
