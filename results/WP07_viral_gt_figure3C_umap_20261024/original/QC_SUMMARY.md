# Initial UMAP QC — 24 October 2026

One joint, label-free cosine UMAP is complete. Figures A–C use the same coordinates. See `README.md` for reproduction and `QC_DETAILS.md` for expanded findings and caveats.

## Counts and representation

| Quantity | Verified count |
|---|---:|
| Original anchors / mapped anchors | 594 / 593 |
| Original 15-gene-workflow windows / successful embeddings | 757 / 737 |
| Complete original windows before eligibility filtering | 616 |
| Baseline merged loci / included loci | 565 / 556 |
| Excluded loci without usable 15-gene embeddings | 9 |
| Full-universe candidates / not selected / unresolved | 288 / 269 / 8 |
| Plotted candidates / not selected / unresolved | 287 / 269 / 0 |
| Loci with multiple original / usable windows | 130 / 129 |
| Reference BGCs / Saccharide-annotated references | 2,502 / 195 |
| Total joint-map points | 3,058 |

Both input arrays have 320 dimensions. The existing CORE_ONLY reference set used in prior viral retrieval is compatible with the actual-gene viral embeddings and was reused. The five-augmentation 128-position consensus is incompatible and was excluded. Original augmentation cores/labels, frozen model/source hashes, gene order, masks, cache vectors and locus/anchor mappings passed audit.

Each included locus uses an observed cosine-medoid window; ties within 1e-12 use stable window IDs. No imputation or 11-gene vector substitution occurred. Twenty original windows lack inference because of translation mismatch (10) or unverified translation (10). Nine loci lack any usable 15-gene-workflow embedding: eight unresolved loci and one prior 11-gene-supported candidate, Bandra megavirus MG779342.1 (`vgl_4304220c157bbfcfc89ea207`). Their original labels remain in the full locus registry.

## Original-space retrieval

The nearest reference is Saccharide-annotated for **222/556 loci (39.93%)**: 208/287 candidates and 14/269 not-selected loci.

| Prior status | n | Mean P@10 | Median | IQR | Range | P@10 ≥0.5 |
|---|---:|---:|---:|---|---|---:|
| Selected | 287 | 0.6930 | 0.8 | 0.5–1.0 | 0–1 | 229 (79.79%) |
| Not selected | 269 | 0.0409 | 0 | 0–0.1 | 0–0.4 | 0 |

All 556 representative P@10 values and nearest-reference IDs reproduce the prior window-level results. Independent SciPy calculations verify retrieval and medoid choice. All neighbor statistics use the original 320D vectors. Selected loci with lower representative P@10 can have classifier, 11-gene, or alternate-window support.

## Findings from the initial figures

1. **Saccharide overlap:** yes. Many viral candidates overlap the central Saccharide-containing area. Viral loci occupy only part of the Saccharide reference distribution; other Saccharide areas have few viral points.
2. **Reference context:** many viral loci intermix with references in the central/right portions, with some extensions containing fewer visible reference points. The large left reference-dominated structure contains few viral loci. These are visual descriptions, not biological clusters or 2D density tests.
3. **P@10 agreement:** high-P@10 loci concentrate in the central/lower Saccharide-overlapping area; lower-P@10 loci dominate several other extensions. This agrees with original-space retrieval, subject to projection distortion.
4. **Artifacts:** no artificial repeated-gene padding, mask mismatch, incompatible consensus mixing, invalid vectors, or row-alignment errors were detected. However, 109 representatives contain only 1–14 actual genes, including four one-gene inputs; 447 contain all 15. Reference cores span 1–115 genes. Gene-count effects remain unresolved. There are 513 unique representative vectors: 71 loci in 28 exact-vector sets share identical ordered sequence-hash contexts, causing overplotting. No biological loci were dropped because their vectors matched.
5. **Next, after joint figure review:** inspect source neighborhoods from the Saccharide-overlapping area, low-P@10 candidates, short inputs and reference-sparse extensions. A later authorized analysis could test seed/parameter stability and complete-15-gene sensitivity. No additional fits or clustering were started.

Candidate selection and P@10 share the same learned representation and Saccharide reference labels; agreement is **not independent validation**. References were used in supervised training. Variable-length isolated inputs versus 128-position training, cross-runtime ESM near identity rather than bitwise identity, first-1,022-residue truncation, missing-window effects on medoids, and viral domain shift remain limitations. UMAP gaps, distances and densities do not establish biosynthetic families.

## Reproducibility and verification

Parameters: n_neighbors=30, min_dist=0.1, cosine, 2 components, seed=42, one thread, spectral initialization, default 500 epochs. Raw 320D vectors enter UMAP without PCA/feature scaling/explicit pre-normalization; retrieval uses float64 L2 normalization. Versions and all defaults are recorded in JSON.

All 634 consumed original files passed before/after hash checks. Independent arithmetic, CSV coordinate round-trip and export checks passed. Each PNG is 5,100 × 3,420 pixels at 600 dpi; all SVG marks, text and color bars are vector. The same saved coordinates were re-rendered to correct SVG color-bar rasterization, without a refit. All three rendered panels were inspected for readable legends, labels, shared axes and clipping. Color-vision-deficiency simulations were not performed.

Stopped after the initial figures and QC. No de novo clustering, retraining, new family assignment, Git commit or push was performed.
