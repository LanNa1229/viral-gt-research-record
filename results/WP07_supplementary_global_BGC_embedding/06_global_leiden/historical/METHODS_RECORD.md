# Methods record for the historical global 320D clustering

This extract preserves the global analysis methods from the archived run. The release source map records the original and published digests. The contemporaneous baseline settings, software versions and code identities are recorded in the global execution evidence matrix and scientific run receipts.

## Datasets and ascertainment

We reused the archived 556x320 observed 15-gene-workflow medoids, preserving locus/window identifiers and row alignment. Primary clustering used only the 287 previously selected evaluable loci. Candidate status was not recalculated. In the prior screen, selection was the union across 11/15-gene windows of Saccharide classifier score >=0.5 and >=Other score, or Saccharide P@10 >=0.5; those existing labels were reused. This cohort restriction creates ascertainment bias. Existing compatible CORE_ONLY reference embeddings contain 2,502 unique BGCs, including 195 with any Saccharide annotation. Original hashes were matched to archival commit 7b35093eefa99a4704e7b64039f51b7564788fd4. Prior model inference and the original UMAP were not rerun.

## Cosine graph

The executed engine uses float64 L2-normalized vectors and exhaustive pairwise dot products. Cosine similarity is clipped to [-1,1]; distances are 1-similarity, with self and exact-vector duplicate distances set to zero. No feature scaling, PCA or whitening is applied. For each node, k other nodes are selected by ascending distance, breaking exact ties by stable sample ID. The undirected graph is the union of directed neighbor relations, with weight max(0,cosine); zero weights and self-loops are omitted. Every biological locus is retained in the primary graph, including exact duplicates. Stable-ID tie handling can influence duplicate/hub edges; it is explicitly audited. The global baseline graph contains 39,384 edges.

## Leiden and stability

The executed objective is leidenalg.RBConfigurationVertexPartition with singleton initialization, no maximum community-size cap, configured resolution and explicit RNG seed. One optimizer iteration is called repeatedly until improvement <=1e-12, capped at 100 calls; the actual counts/improvements and convergence flag are saved for every fit. The baseline is k=20, resolution=1.0, seed=42. The grid crosses k in 10/15/20/30, resolution in 0.5/1.0/1.5/2.0 and seeds 0-4, giving 80 grid fits plus the baseline for the global cohort. Parameter choice did not use annotations, P@10, desired counts or UMAP appearance. Neutral cluster identifiers are ordered by lexicographically smallest member ID, not by function.

## Software and execution

The successful scientific command, working directory, exact executable source SHA256 and code snapshots are in runs/run_001_science/before.json. Python is 3.12.14 (main, Aug 25 2026, 14:01:42) [MSC v.1944 64 bit (AMD64)]. Distribution versions were recorded at execution: {"igraph": "1.0.0", "leidenalg": "0.12.0", "matplotlib": "3.11.1", "numpy": "2.3.5", "pillow": "12.3.0", "scikit-learn": "1.9.0", "scipy": "1.18.1", "texttable": "1.7.2"}. Full relevant Python/native dependency-file fingerprints are retained for the scientific run. Explicit numerical-library thread variables were set to one.

## NOT RECORDED and limits

The effective optimization epoch count of the **old UMAP** remains NOT RECORDED; it is irrelevant to new 320D clustering and was not reconstructed from installed defaults. A raw original prior-inference execution transcript is NOT RECORDED in this clustering execution; saved source/embedding identities and prior archival proofs are used. No functional identity, pathway identity or new annotation is recorded or inferred. The original archive's historical limitations remain in force. Current clustering iteration counts, parameters and versions are recorded contemporaneously. The inherited full process environment is intentionally not dumped because it can contain credentials; explicit numerical settings, package inventory, code/dependency identities and the actual command are retained.

## Global reference comparison and visualization

The label-free global joint Leiden grid uses 2,789 samples, 287 viral loci and all 2,502 reference BGCs. The GC identifiers are neutral labels. Overlapping reference class annotations are retained. Saved UMAP positions are reused for display. Leiden uses the original embedding space.
