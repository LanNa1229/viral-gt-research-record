# Initial 320D viral-locus clustering report

The prespecified Leiden baseline detects **five exploratory communities** (VC01=68, VC02=55, VC03=74, VC04=35, VC05=55), but the combined evidence does **not support treating this five-way partition as a reproducible set of distinct biological groups**. No community passes every prespecified operational stability flag. The primary graph is connected, HDBSCAN finds mostly noise under the tested settings, and inclusion of unselected loci changes community membership.

| Community | Loci | Genome units | Full 15 genes | Median P@10 | Seed Jaccard | All-556 Jaccard |
|---|---:|---:|---:|---:|---:|---:|
| VC01 | 68 | 63 | 79.4% | 0.9 | 0.588 | 0.840 |
| VC02 | 55 | 52 | 85.5% | 0.6 | 0.703 | 0.662 |
| VC03 | 74 | 63 | 87.8% | 0.5 | 0.923 | 0.526 |
| VC04 | 35 | 35 | 88.6% | 0.8 | 0.971 | 0.306 |
| VC05 | 55 | 54 | 83.6% | 1.0 | 0.782 | 0.553 |


## Stability and algorithm comparison

At k=20 and resolution=1.0, pairwise seed ARI across seeds 0-4 ranges 0.695-1.000, median 0.720. Across the prespecified grid, the primary cohort yields 2-13 communities. VC03 and VC04 have the strongest same-setting seed median best-Jaccard values (0.923 and 0.971), but both change substantially in the all-556 sensitivity. VC01, VC02 and VC05 have lower same-setting seed stability. Higher stability at a different setting was not used to replace the prespecified baseline.

Baseline HDBSCAN (min_cluster_size=15, min_samples=5, sklearn convention including self) finds **zero clusters and 100% noise**. Five of six HDBSCAN settings assign all points to noise. The mcs10/ms5 setting finds three small groups of sizes 16, 14 and 12, with **245/287 (85.4%) noise**. Its ARI to Leiden is 0.0067 when noise is a category. The nonnoise-only ARI is 0.893, but covers only 42 loci (14.6%): it is not broad algorithmic agreement. HDBSCAN membership strengths are not probabilities of biochemical function.

## Sensitivities

- **All 556:** ARI on the 287 shared selected loci is 0.572. VC01 remains intact but acquires other selected members (best Jaccard 0.840); VC02 is partly retained (0.662). VC03, VC04 and VC05 change more (0.526, 0.306 and 0.553). The selection-conditioned five-way partition does not persist as a whole.
- **Full 15 genes:** retaining 243 selected loci, excluding all 44 shorter representatives, gives ARI 0.910; per-community best Jaccard is 0.861-1.000. Short windows do not by themselves explain the principal baseline pattern. Extreme short-input/hubness examples are explicitly listed in tables/extreme_short_window_QC.csv.
- **Duplicate collapse:** 287 loci contain 263 distinct vectors/ordered contexts. The two exact-duplicate equivalence relations coincide. Thirty-eight selected loci belong to 14 repeated-context sets; collapsing and mapping assignments back gives ARI 0.821, with Jaccard 0.763-0.947. Duplicates influence details but do not account for all structure; biological loci were retained separately in the primary analysis.
- **Genome/relatedness:** communities contain 35-63 distinct genome units and 30-58 pre-existing sketch components. No single genome contributes more than 7.3% of a community, and the largest relatedness component contributes at most 14.7%. Genome-balanced and existing k21-sketch-component-balanced graph weights give ARIs 0.876 and 0.782. These controls change some boundaries; they do not establish independence of related contexts. Sketch components are not ANI/species definitions.

## Graph and distance QC

The baseline union kNN graph has 4176 edges, density 0.1018, **one connected component** and no isolated nodes. Weighted modularity is 0.5623; silhouette in original cosine space is 0.200. A positive modularity and modest silhouette do not establish a unique natural partition. Maximum directed kNN in-degree is 56 (mean k=20); in-degree Gini is 0.326. Exact-duplicate edges are 39/4176 (0.93%). Twenty-nine nodes have exact ties at the k=20 boundary, handled by stable IDs. Full similarity, nearest-distance and k-th-distance distributions are saved; absolute cosine magnitude is not interpreted as biological separation.

## Global and Saccharide reference context

The global joint dataset (2,789 samples) has 10 baseline communities; six contain viral loci. Its labels on selected viral loci have ARI 0.380 with viral-only labels. GC08 combines most VC01, all VC05 and part of VC04, with 147 viral loci and 105 references (100 carrying Saccharide). GC05 contains 107 viral loci and 53 references, mixing VC02/VC03 with a VC01 locus. Viral-only and global partitions therefore differ substantially.

The deliberately Saccharide-focused joint dataset (482 samples) has seven baseline communities; six mix viral loci with Saccharide references. ARI on shared viral loci is 0.640. VC03 falls entirely into SC03, VC04 mostly into SC07, and VC05 divides across several Saccharide-focused communities. VC01 and VC05 are comparatively close to Saccharide references (median nearest-Saccharide cosine 0.982 and 0.991); VC02-VC04 have lower medians (0.953-0.959). These are relative descriptive comparisons, not calibrated distance thresholds.

The baseline joint graph's per-community reference contacts and zero-reference-neighbor counts are in tables/community_reference_graph_contacts_summary.csv. Reference-poor local neighborhoods can occur even when the containing joint community includes references; counts are explicitly tied to the chosen kNN graph, not inferred from UMAP islands. Hybrid reference classes are retained as overlapping memberships, so marginal class fractions may sum above one. The Saccharide-only comparison is not an enrichment test.

## Scope and interpretation

Candidate ascertainment already used the 11/15-gene classifier and Saccharide P@10, and the reference BGCs were used in supervised training. Functional labels, P@10, taxonomy, annotations and UMAP positions were not clustering features or parameter-selection criteria. Agreement with Saccharide references and algorithmic agreement are not independent functional validation. The robust result is the **sensitivity pattern**, not a validated five-family classification. Inspect these communities together before selecting any biological follow-up.

No new embeddings, BGC-Prophet inference, UMAP fits, functional annotation, pathway assignment or representative-neighborhood selection was performed. The original UMAP analysis and archive remain unchanged.

## Matched-setting cross-grid comparison

All 80 prespecified grid settings were also compared between selected-only and each corresponding sensitivity/joint cohort on shared selected loci. This uses saved partitions only, with no new fits or parameter selection.

| Comparison | Shared loci | Matched-grid ARI min | Median | Max | Baseline ARI |
|---|---:|---:|---:|---:|---:|
| all556 | 287 | 0.321 | 0.553 | 0.776 | 0.572 |
| selected_full15 | 243 | 0.470 | 0.716 | 0.951 | 0.910 |
| selected_deduplicated | 287 | 0.472 | 0.700 | 0.986 | 0.821 |
| global_joint | 287 | -0.010 | 0.238 | 0.654 | 0.380 |
| saccharide_joint | 287 | 0.469 | 0.678 | 0.950 | 0.640 |
