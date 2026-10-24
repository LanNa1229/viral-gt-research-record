# BGC-Prophet inference audit — 2026-10-17

Audit scope: read-only inspection of current local source, latest 2026-10-16 benchmark READMEs/scripts, metadata and numerical array headers/contents. No model inference, retrieval, downloads, training or historical-file mutation was executed in this audit. Paths below are relative to `/home/rahnn/260908_MiBG` (same workspace as the Windows UNC checkout).

## Reusable resources and verified identities

- Exact published model source: `sources/upstream/bgc_prophet/train/classifier.py` and sibling `model.py`. This is loaded directly without installing the upstream package.
- Frozen local checkpoint: `analysis/viral_glycan_benchmark_20261016/sources/classifier.pt`.
- Existing encoder checkpoint: `analysis/bgcprophet_sequence_recovery_20260912/model_cache/hub/checkpoints/esm2_t6_8M_UR50D.pt`.
- Existing WSL runtime: `analysis/bgcprophet_sequence_recovery_20260912/.venv/bin/python`; historical inference metadata records torch 2.5.1+cpu, numpy 1.26.4, fair-esm 2.0.0, biopython 1.85 and lmdb 2.3.0. This audit did not execute that interpreter or probe current GPU availability.
- Reusable verified wrapper: `analysis/viral_glycan_benchmark_20261016/scripts/benchmark.py`, particularly `torch_setup`, `load_model`, `infer`, `predicted`.
- Latest extension pattern: `analysis/viral_glycan_benchmark_cluster7_ARF08726_ARF08739_20261016/scripts/cluster7.py:141` (sequence-hash cache reuse and direct frozen inference).
- Historical all-position export: `recovered_historical_scripts/extract_bgc_embeddings.py` hooks the input to the original classifier head.
- Published reference protein vectors: `analysis/bgcprophet_sequence_recovery_20260912/lmdb_snapshot/data.mdb`, historically verified 15,081,525,248 bytes. Open readonly=True, lock=False, create=False, readahead=False. Do not reconstruct reference vectors from inferred sequences unnecessarily.

The following current file hashes were actually recomputed and equal the historical constants:

| Resource | SHA256 |
|---|---|
| classifier.pt | f45ca298980d4740dd43131a7e69c99cb233368c815470572cfc46f33abf4eb5 |
| ESM2 checkpoint | 46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928 |
| classifier.py | fa4e1d045ac2c94559a9ab41054fd6b28d991d82be63e98585d7bc014564dca8 |

The 15-GB LMDB was not rehashed during this short audit; its historical hash/readonly validation is documented in the recovery README and benchmark source-integrity outputs.

## Original representation and important quirks

The training CSV `BGC_train_dataset_classify.csv` contains 12,510 rows = five augmentations of 2,502 biological BGCs. `sentence` is each ordered biological core; `TDsentence` is the 128-position augmented input. `TDlabels` contains core/flank membership. The loader `sources/upstream/bgc_prophet/train/data.py` fetches each protein's layer-6 `mean_representations` from LMDB in TDsentence order, stacks them, and passes TDlabels as a float32 tensor.

The original forward computes a positional encoding and discards that encoded tensor: the encoder receives `src` directly. Preserve this published behavior. The float32 0/1 tensor passed as `src_key_padding_mask` is an additive key-attention bias: 1=core and 0=flank, with no boolean exclusion of flank positions. It is not a padding mask in the ordinary True/False exclusion sense. The original final head uses the mean of all contextualized positions, followed by a linear map and seven sigmoid outputs.

The general upstream genome split command repeats short genomes until 128 positions and truncates; it even initializes a two-value placeholder TDlabels field before detection. Do not reuse that command for these GT-centered windows. It is a genome-detection input builder, not an appropriate 15-gene-padding implementation.

The validated viral workflow supports variable-length CORE_ONLY inference through the exact original forward. For a GT window use its actual ordered n genes (n<=15), no padding, no repetition, and `[1]*n` membership. This makes actual-window pooling and published all-input-position pooling identical. Boundary-truncated windows remain n-gene inputs. Record variable-length inference as a distribution limitation relative to the original 128-position training construction.

An optional genuine 128-gene context comparison is a different input strategy, available only when 128 authentic genes exist on the same contig. If that input uses the 15 window positions as membership=1 and flanks=0, its contextualized 15-position mean is a supplementary representation; original class scores still come from all 128 positions. Never manufacture class scores for the supplementary pooled vector.

## Concrete callable API

Import the original benchmark module by explicit path or prepend only its script directory, avoiding same-named validator imports. Merely importing it does not execute stages. Keep its module OUT pointing to its original directory so load_model locates the existing classifier. Do not call its prepare/esm/pilot/reference functions, since those write historical outputs.

```python
# Execute with the existing WSL environment from the workspace root.
import importlib.util
from pathlib import Path
p = Path('analysis/viral_glycan_benchmark_20261016/scripts/benchmark.py').resolve()
spec = importlib.util.spec_from_file_location('bgc_frozen_20261016', p)
frozen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(frozen)
torch = frozen.torch_setup()          # once per process; sets CPU thread counts
model = frozen.load_model(torch)      # once; validates source/checkpoint hashes
# vectors is (n,320), ordered exactly like the extracted genomic window.
scores, window_embedding, membership_embedding = frozen.infer(
    torch, model, vectors, [1] * len(vectors), check_hook=True)
labels, threshold = frozen.predicted(scores)
```

The wrapper returns shapes `(7,)`, `(320,)`, `(320,)`; checks finite output and that the recorded all-position mean equals the actual classifier input; check_hook=True verifies hooked/unhooked scores bitwise. Its classifier run deliberately keeps autograd enabled, reproducing the previously validated classify behavior. Do not silently wrap it in inference_mode or change dtype/device until regression confirms equivalence. It handles one window per call; loop with one loaded model for the pilot. Batched optimization would require a separate equivalence check. Mixed variable lengths should be grouped by length if batching is later implemented; unvalidated zero padding changes behavior.

Class order: Alkaloid, Terpene, NRP, Polyketide, RiPP, Saccharide, Other. Saccharide is index 5. These are sigmoid probabilities/scores and need not sum to one. Reporting threshold is max(0.5, Other score); retain all classes meeting it, remove Other if another class qualifies, otherwise Unknown when empty. Direct inference bypasses only the upstream isBGC detection/reporting gate.

## ESM preprocessing and cache reuse

Use the same esm2_t6_8M_UR50D checkpoint, eval mode, float32, layer 6, BOS/EOS excluded, residue mean. Existing converter uses truncation_seq_length=1022; mean slice is `[1:min(len(seq),1022)+1]`. Preserve full sequences and report aa_length, encoded_length and truncation. Do not average BOS/EOS or change the encoder/model layer.

The following current cache pairs were inspected; all arrays are finite float32 and have shape (metadata_rows,320):

| Benchmark directory under analysis | Rows |
|---|---:|
| viral_glycan_benchmark_20261016 | 20 |
| viral_glycan_benchmark_A0A0K0N728_20261016 | 31 |
| viral_glycan_benchmark_cluster7_ARF08726_ARF08739_20261016 | 160 |

Each pair is `protein_embedding_metadata.tsv` + `protein_esm2_embeddings.npy`. Total 211 rows and 211 distinct full-sequence SHA256 hashes. Cache by sequence_sha256, not accession: overlapping windows reuse vectors while preserving separate observations. Validate row alignment, dimensions, finite values, current model hash/preprocessing identity, and duplicate-hash vector consistency before accepting reuse. Persist new vectors under the new analysis, without rewriting these caches. The legacy scripts have only in-process hash dictionaries plus final arrays; a scalable restartable cache should save missing vectors incrementally with checkpoint/preprocessing identity alongside each entry.

Reference LMDB vectors and current CPU ESM extraction use equivalent checkpoint/layer/pooling/truncation, but historical numerical identity is not bitwise: the sequence-recovery pilot documented 50/50 near-identity matches and 0/50 matches at the initial strict tolerance. Its a posteriori relaxed criterion was atol=rtol=1e-3 plus cosine>=0.999999 and expected-ID nearest match. Preserve that caveat; do not claim raw CPU-generated viral vectors are bitwise original-runtime vectors.

## Reference compatibility matrix

All shapes/dtypes below were actually opened and all values checked finite during this audit.

| Reference | Actual shape / dtype | Compatible query construction |
|---|---|---|
| `analysis/viral_glycan_benchmark_20261016/mibig_CORE_ONLY_embeddings.npy` | (2502,320), float32 | Actual ordered window genes only, all membership=1, original forward + mean over those actual genes. Primary new workflow. |
| `analysis/viral_glycan_benchmark_20261016/mibig_REPEAT_128_embeddings.npy` | (2502,320), float32 | Repeated-core input to 128, original all-position mean. Historical comparison only; new task does not authorize artificial repeats as primary. |
| `analysis/bgcprophet_reference/results/consensus_embeddings.npy` | (2502,320), float64 | Original 128-position augmented-context all-position mean, averaged across the five original augmentations. Not compatible with isolated CORE_ONLY vectors. |
| `analysis/core_context_consensus_step3_20260911/results/core_context_consensus_embeddings.npy` | (2502,320), float64 | Mean of contextualized core positions inside original 128-position contexts, averaged across augmentations. Not the isolated CORE_ONLY space. |
| `analysis/bgcprophet_reference/data/mibig_12510_embeddings.npy` | (12510,320), float32 | Individual historical 128-position augmentation rows; do not count as 12,510 independent BGCs. |

CORE_ONLY / REPEAT_128 share `analysis/viral_glycan_benchmark_20261016/mibig_matched_reference_metadata.tsv`. Every biological reference's five ordered core lists and labels were checked identical before collapsing and reinference with original LMDB vectors. `reference_validation.json` records 2,502 rows, original-row regression pass and no repeat truncation. Original/source labels are retained; Saccharide positives include hybrid labels (historically 195/2502).

Original consensus metadata is `analysis/bgcprophet_reference/results/consensus_metadata.tsv`; original row metadata is `analysis/bgcprophet_reference/data/mibig_12510_metadata.tsv`. Check every embedding_row against metadata. Do not pair one array with another strategy's PCA/UMAP fit. The original CORE_ONLY reference is already suitable for actual-window inference: no new reference regeneration is necessary for the primary strategy. A genuine-context 128-position query retains a documented authentic-versus-synthetic context difference even when pooling is matched.

Quantitative neighbor retrieval should use cosine in original 320D with accession tie breaking; PCA/UMAP are displays, not quantitative retrieval coordinates. The references were used in supervised model training; localization is descriptive recognition, not held-out performance or chemical validation.

## Known candidate provenance and broad-list pointers

Existing local genome annotations/sequences: S-H68 `analysis/viral_glycan_benchmark_20261016/sources/MK016663.1.gb` (117 translated CDSs, partial contig); GMA2 `analysis/viral_glycan_benchmark_A0A0K0N728_20261016/sources/KR063281.1.gb` and NC_070761.1 counterpart (126 translated CDSs); Catovirus `analysis/viral_glycan_benchmark_cluster7_ARF08726_ARF08739_20261016/sources/KY684083.1.gb`. Their genome_protein_manifest.tsv files retain complete ordered protein/CDS annotations.

S-H68 structural accession A0A482ICV6 maps by archived exact sequence to QBP05883.1, inside the old manual QBP05878–QBP05897 core. Catovirus structural accession A0A1V0SAJ1 maps by exact archived sequence to ARF08732.1, inside manual ARF08726–ARF08739; current annotation for that protein is hypothetical. Choosing either as the pilot structure-mining anchor must be labeled historical structural association, not experimentally verified GT identity. New +/-7 windows must be compared to old manual cores without assuming equality or reusing old window scores.

The GMA2 extension points outside this workspace to `/home/rahnn/projects/viral-gt-research-record/results/WP01B_coarse_clustering/canonical/cluster_membership.tsv` and `figure_X-1-1/network_source/cytoscape_complete_edges.tsv`. These are broad coarse structural memberships, not automatically the requested user-approved GT accession population. Catovirus SSN provenance points to `/home/rahnn/C7_local_reconcile_merge/C7_core14_ids.txt`, `CLUSTER 7 and neighbor_node.csv`, and `/home/rahnn/260728_Temp/output/preflight_r1_seed1_v3_correct_schema/nodes_all_leiden_assignments.csv`. This audit did not inspect those external files or select a full-scale population from them. The parent extraction audit should reconcile dataset manifests and report absence of a definitive supplied list rather than silently promote all coarse-cluster members.

## Pilot gate and execution limits

A safe next execution is actual-window inference for verified S-H68/Catovirus mappings, importing cached ESM vectors where available and encoding only missing window proteins. Verify n<=15 ordered input, n memberships, no duplication/padding, finite (7,) scores and (320,) vectors, unchanged checkpoint weights/hashes, and one manifest row per output row. All-ones mask means actual-window and original all-position means should agree exactly. Compare manual-region membership by accession, not old predictions.

No new inference success, runtime benchmark, GPU throughput or full-scale biological count is claimed here. The existing environment is CPU-only in historical metadata; current device/throughput and full accession population remain for the parent to verify. For n windows, float32 final embeddings require 1280*n bytes; the raw unique protein cache requires 1280*U bytes plus metadata, sequences and temporary ESM activations. Encoder work, not classifier state, is likely dominant, but runtime must be measured on the pilot.

## Implementation handoff and subsequent cohort clarification

During implementation the user selected the complete 594-member coarse structural cluster with representative A0A3Q8Q3U2 as the baseline cohort. The earlier broad-list caveat describes the initial audit state; the parent now owns cohort materialization and source reconciliation. All 594 members must retain coarse-structural-selection provenance; membership alone does not establish experimentally validated GT function.

Implemented adapters reside in analysis/viral_gt_mining_20261017/checkout/src/viral_gt_mining/inference.py and analysis.py. Inference keeps the frozen original wrapper read-only, caches by full sequence and model settings, emits all four requested B files, and optionally writes separate descriptive matched-reference localization. Analysis preserves raw windows, merges overlapping intervals only for locus summaries, and distinguishes unknown/incomplete inference from negative classifier results. PCA uses the original reference-only CORE_ONLY fit and checks its original file hashes plus reference coordinate reproduction; rendered plots still require visual review. No full inference or analysis run was executed by the preparatory audit agent.
