# Initial UMAP QC and interpretation — 24 October 2026

The requested initial exploratory analysis is complete. One label-free joint UMAP contains **2,502 reference BGCs and 556 viral loci (3,058 points)**. Figures A–C use exactly the same coordinates and limits. No new clustering, model inference, family assignments, Git commit, or push was performed.

## Verified inputs and units

| Quantity | Verified count |
|---|---:|
| Original structural-cohort GT anchors | 594 |
| Anchors mapped to at least one original window | 593 |
| Original occurrence-specific 15-gene-workflow windows | 757 |
| Complete 15-gene windows before inference eligibility filtering | 616 |
| Successfully inferred windows, all finite and nonzero | 737 × 320D |
| Excluded original windows | 20 |
| Baseline overlap-merged loci | 565 |
| Loci with usable 15-gene-workflow vectors | 556 |
| Loci excluded from UMAP for missing vectors | 9 |
| Prior candidates in the full baseline universe | 288 |
| Prior not-selected / unresolved loci in full universe | 269 / 8 |
| Candidates / not-selected loci actually plotted | 287 / 269 |
| Unresolved loci plotted | 0 |
| Loci with multiple original windows | 130 |
| Loci with multiple usable windows for medoid selection | 129 |
| Reference BGCs, each represented once | 2,502 × 320D |
| References carrying Saccharide, including hybrids | 195 (7.79%) |
| Reference BGCs with multiple existing class labels | 467 |

The 20 excluded windows comprise 10 translation mismatches and 10 unverified translations. They affect more than the nine completely unavailable loci: a locus can retain another usable window. All nine unavailable loci are listed in `tables/locus_registry.csv`; eight have unresolved prior union outcomes. The ninth is the prior 11-gene-supported Bandra megavirus candidate on MG779342.1, locus `vgl_4304220c157bbfcfc89ea207`, window `vgt_366f3de7fcfac8360f3b2cee`. It has no usable 15-gene embedding and is excluded rather than substituted. The unmapped original anchor remains outside the 565-locus universe.

The five-augmentation 128-position consensus was deliberately excluded because its contextual representation is incompatible with the saved isolated-window viral vectors. The existing matched CORE_ONLY reference set used in the previous viral retrieval analysis was reused instead. Checkpoint, forward-source, reference/proof, embedding, and gene-manifest identities passed the recorded hash gates. Each viral vector was linked to its original input GT anchor(s), source genome/window, gene order, baseline locus, and frozen cache. Cached actual-position and membership-pooled vectors match exactly.

One observed 15-gene-workflow vector per locus was selected by deterministic cosine medoid; scalar locus maxima/means in earlier reports were not treated as representative embeddings. All 556 representative P@10 values and nearest-reference identities reproduced their saved window-level values. Independent SciPy cosine calculations verify medoid choice and retrieval without using UMAP coordinates.

## Quantitative diagnostics in the original 320D space

**222/556 loci (39.93%)** have a nearest reference BGC carrying Saccharide. This is 208/287 selected loci (72.47%) and 14/269 not-selected loci (5.20%). Reference accession resolves exact cosine ties.

| Prior candidate status | Plotted loci | Mean P@10 | Median P@10 | Interquartile range | Range | P@10 ≥0.5 |
|---|---:|---:|---:|---|---|---:|
| Selected by prior 11/15-gene union | 287 | 0.6930 | 0.8 | 0.5–1.0 | 0–1.0 | 229 (79.79%) |
| Not selected by prior union | 269 | 0.0409 | 0 | 0–0.1 | 0–0.4 | 0 |

The 58 selected loci with representative P@10 below 0.5 are compatible with the original rule: selection can arise from classifier evidence, an 11-gene input, or a different 15-gene window at the same locus. Figure C always uses the plotted representative's own 15-gene-workflow P@10. Complete discrete distributions are exported in `tables/P10_group_frequencies.csv`.

Candidate selection and P@10 are related analyses based on the same learned representation and Saccharide reference labels. Their agreement is **not independent validation**. These reference BGCs were also used in supervised model training; the results do not measure held-out accuracy or chemical function.

## Answers from visual inspection

1. **Viral overlap with Saccharide references:** yes, a substantial set of viral loci overlaps the central Saccharide-containing area. Prior selected loci are concentrated there. Viral loci do not cover all Saccharide reference areas; Saccharide references also occur in the reference-dominated left portion with few viral points.
2. **Reference-rich versus separate regions:** many viral loci intermix with reference BGCs in the central and right portions. Some extend along arms/bridges with fewer nearby visible reference points. Most of the large left reference-dominated structure contains few viral loci. These are visual observations, not independently quantified 2D density or distance tests, and not biological clusters.
3. **High-P@10 agreement:** yes. Yellow/green high-P@10 loci concentrate in the central/lower area containing Saccharide references; lower-P@10 loci dominate the upper central and right extensions. This agrees with original-space retrieval, with expected local UMAP distortion. No nearest-neighbor statistic was computed in UMAP space.
4. **Obvious artifacts:** no artificial repeated-gene padding, mask-length mismatch, incompatible consensus mixing, invalid vectors, or row-alignment errors were detected. Important remaining visual limitations are listed below. This initial projection cannot establish that gene-count effects or domain shift are absent.
5. **Next investigation, after reviewing these figures together:** inspect representative genomic neighborhoods from the Saccharide-overlapping area, low-P@10 selected loci, and reference-sparse extensions. Review the short-input and exact-duplicate-context examples alongside their source annotations. If a later analysis is authorized, assess seed/parameter stability and complete-15-gene sensitivity before interpreting projected separation. No such additional fits, sensitivity runs, or clustering were initiated here.

## Artifact and methodological caveats

- **Variable gene counts:** 447 representatives contain all 15 genes; 109 use actual boundary/short-genome inputs of 1–14 genes, including four one-gene representatives. These are original 15-gene-workflow inputs, not substituted 11-gene embeddings. Complete representatives have mean/median P@10 of 0.3904/0.2; short representatives have 0.3248/0.1. This descriptive difference is confounded by cohort composition and does not identify a gene-count cause. Reference cores span 1–115 genes. A gene-count effect cannot be ruled out from this map.
- **Exact overplotting:** 556 distinct loci have 513 unique representative vectors. There are 28 exact-vector duplicate sets involving 71 loci, each sharing the same ordered protein-sequence-hash context. These reflect repeated observed contexts across distinct frozen biological loci, not manufactured padding. Their points can coincide or nearly coincide; visible point count is not locus count. `tables/exact_duplicate_vector_QC.csv` preserves these equivalences. No duplicate loci were silently discarded.
- **Pooling/augmentation:** both sides use isolated actual-gene mean pooling. The original 128-position training construction differs from these inputs; compatible pooling does not remove that domain difference. The five augmented contextual vectors are not averaged into this map.
- **ESM and long proteins:** the historical LMDB versus CPU protein-vector comparison supports near identity rather than bitwise cross-runtime identity. The first-1,022-residue encoder convention remains in force.
- **Locus medoids:** choosing an observed representative avoids synthetic locus averages but cannot capture all within-locus heterogeneity. The medoid is based on available embeddings; missing windows cannot influence it. All alternative window scores and selected IDs remain reviewable.
- **Projection:** UMAP is nonlinear. Inter-island distances, apparent gaps, areas, and point density are not validated biological boundaries. Reference labels were not passed to the fit. Only one seed/parameter setting was requested and executed.

## Reproducibility and export QA

The fit uses cosine distance, n_neighbors=30, min_dist=0.1, n_components=2, random_state=42, n_jobs=1, spectral initialization, and the installed UMAP default of 500 optimization epochs for this sample size. Raw 320D inputs are promoted to float64 without PCA, whitening, feature scaling, or explicit normalization before UMAP; the library performs its internal numerical operations. Retrieval uses float64 unit vectors, consistent with the prior cosine analysis.

Runtime: bundled Windows Python 3.12.14; NumPy 2.3.5; SciPy 1.18.1; scikit-learn 1.9.0; umap-learn 0.5.12; Numba 0.67.0; PyNNDescent 0.6.0; Matplotlib 3.11.1; Pillow 12.3.0. Exact platform and all constructor defaults are recorded in `execution_provenance.json`. All **634 consumed original input/proof/cache files** passed before/after SHA256 integrity checks.

Each PNG is 5,100 × 3,420 pixels at 600 dpi. The SVG scatter marks, text and color bar are vector elements. A color-bar rasterization found during export QA was corrected by re-rendering the same saved coordinates; no UMAP refit occurred. Legends, axis labels, annotations, color scale, point visibility, consistent axes and absence of clipping were inspected in the exported contact sheet and panels. Class color remains essential to distinguishing seven memberships in Figure A; common-color-vision-deficiency simulations were not performed. Reference circles versus viral triangles provide a redundant source-type cue, and Figure C uses a sequential viridis scale.

Stop here for joint inspection and interpretation of the initial figures.
