# Figure B: viral taxonomy of exploratory Saccharide-associated GT neighborhoods

Generated 18 October 2026. All files here are downstream derivatives. No model inference, ESM2 encoding, structure prediction, public upload, or Git push was performed. Frozen inputs were verified by SHA256 before and after analysis.

## Results and recommended manuscript claim

The frozen union contains **288 fixed baseline loci across 245 of 466 genome units**. Genome support is **190 Both / 29 11-only / 26 15-only**. All 190 Both genomes contain at least one same baseline locus supported by both sizes (**217 such loci**); none are supported exclusively at different loci. A baseline locus is a connected overlap region established using the 15-gene windows, not proof of identical module boundaries or an experimentally defined operon.

Of 466 selected units, 461 have at least one inferred window. Union outcomes are 245 candidates, 215 fully analyzed negatives, and 6 unresolved. The resolved denominator is **460**, not 466 or 461. Positive units may have other unresolved windows; a negative requires every relevant window to have inference in both conditions. One 11-only candidate has an unresolved 15-gene outcome, and one 15-only candidate has an unresolved 11-gene outcome. “Only” means observed candidate support, not a confirmed negative in the other condition.

| Interpretation | Selected genome units represented | Candidates | Candidate-bearing taxa | All represented taxa |
|---|---:|---:|---:|---:|
| Conservative family | 156 | 71 | 15 families | 24 families |
| Extended recorded family | 362 | 196 | 19 families | 27 families |
| Non-generic genus | 130 | 65 | 43 genera | 75 genera |

The family and genus coverage sets differ: a verified genus can be retained under its verified higher ancestors despite a missing family rank. Consequently, genus totals are not a partition of the conservative family totals. Recorded candidate genus assignments span 47 genera before generic-source filtering. The genus tree withholds 336 selected units: 260 generic-source units (150 candidates) and 76 non-generic units without a genus rank (30 candidates). Exact exclusion accounting for every view is in `tables/exclusion_summary_all_views.tsv`.

**Supported claim:** Within this GT-A structural-anchor-selected cohort, genomic neighborhoods with Saccharide-associated model predictions or reference retrieval occur across diverse recorded viral taxonomic groups, including 15 conservatively assigned families and 43 non-generic genera. These numbers describe taxonomy coverage of exploratory candidates; they establish neither glycan biosynthesis nor evolutionary homology of the neighborhoods.

**Important uncertainty:** 174/245 candidates are withheld from the conservative family tree: 150 generic-source candidates and 24 non-generic candidates without a family rank. Generic-source exclusions comprise 260 selected units (258 with any inference), while the non-generic missing-family group comprises 50 units (all with inference). Thus 310 selected / 308 partly or fully inferred units are withheld; see `tables/taxonomy_coverage.tsv`. The retained view has 71 positives, 82 negatives, and 3 unresolved. Extended family coverage excludes 104 units, including 49 candidates, because no family rank is recorded. No excluded unit is silently reassigned to a family.

Taxid **2100421** accounts for 148 selected units and 94 candidates. Its generic “uncultured Caudovirales phage” source label inherits Peduoviridae / Maltschvirus in this NCBI snapshot. These placements are excluded from the conservative family and genus trees. The extended tree flags their generic provenance. Sharing this source taxid does not establish close sequence relatedness. Whole-genome exact-hash alias collapse removes exact counterpart duplicates only; remaining closely related genomes are not statistically independent evolutionary replicates. No prevalence confidence intervals or inferred independent lineage counts are presented.

Of the 245 union candidates, 199 have at least one Saccharide classifier-positive window in either condition; **46 have retrieval-only support across both conditions**. Classifier-positive genomes may also contain retrieval-only windows. The optional evidence track uses the mutually exclusive genome partition (199 + 46), not counts of windows. Classifier and embedding retrieval derive from the same trained representation and are not independent biological validations. P@10 ≥ 0.5 remains an exploratory prioritization threshold.

## Source tables and definitions

Paths are relative to the repository root:

- `analysis/viral_gt_mining_20261017/runs/baseline_594_v2/`: frozen 15-gene windows, predictions, retrieval statistics, genome units, fixed loci, anchor mappings, and taxonomy cache.
- `analysis/viral_gt_mining_20261017/runs/window11_594_v1/`: frozen 11-gene predictions, retrieval statistics, and `window15_window11_mapping.tsv`.
- `analysis/viral_gt_window_size_20261017/tables/candidate_recovery_comparison.tsv`: independently checked registry at locus and genome levels.
- `analysis/viral_gt_interpretation_20261017/tables/genome_relatedness_metadata.tsv`: pre-existing generic-source flags.

`tables/master_genome_dataset.tsv` contains exactly one row per frozen genome unit, retaining source accession aliases, organism labels, taxid, recorded ranks, taxonomy hashes, statuses, completeness, evidence, and counts. The original grouping is exact whole-genome nucleotide SHA256 plus source taxid; missing values retain accession-specific units. This preserves 466 units from 515 source accession aliases. Anchor counts are reported separately as unique original cohort accession IDs and genomic centers collapsed across aliases by coordinates plus strand; neither is summed across overlapping windows to infer locus counts. Original cohort IDs can occur in several genomes. The unresolved original cohort anchor UPI001FBC14AE is not an analyzed genome unit and cannot be added to the denominator.

The locus audit retains the original 565 baseline regions, not the 575 native shortened regions. The 11/15 paired mapping links the same genomic anchor occurrences. A candidate at either size has `saccharide_predicted=True` OR `saccharide_fraction_top10 >= 0.5`. Missing inference is unresolved, unless another window supplies positive evidence. Exact physical overlap grouping is inherited unchanged; aliases and multiple anchors do not multiply genome or locus counts.

“Conservative” is an operational taxonomy filter: a non-generic source label and an explicit family/genus rank in the cached NCBI XML. It does not certify an accession's biological identity or current ICTV ratification. Generic flags retain the original case-insensitive rule: `uncultured|unclassified|environmental|metagenom|\bsp\.`. Non-generic missing ranks remain unresolved. Generic-with-missing-family records are counted in the generic exclusion group first, making reasons mutually exclusive.

## Taxonomy and template provenance

All 220 source taxids were checked against five cached **NCBI Taxonomy EFetch XML** batches retrieved **17 October 2026, 05:39:00–05:44:56 UTC**. Every recorded family/genus/species taxid matches its source lineage; shared node names, ranks, and parent paths were checked for consistency. Copied XML, original request metadata, checksums, retrieval timestamps, and URLs are in `sources/taxonomy/`. This is a timestamped live EFetch snapshot; NCBI supplied no numbered taxonomy release. No separate ICTV release was loaded or claimed.

Newick internal labels contain the recorded scientific name and NCBI taxid. Tips are stable `family_<taxid>` or `genus_<taxid>` IDs. Every hierarchy edge follows a verified source lineage. There are no branch lengths, bootstrap values, inferred sequence relationships, or artificial unresolved tips. Unary taxonomy nodes are preserved in Newick; local previews collapse unary nodes for readability while retaining verified ancestry. `sources/taxonomy/nodes.tsv` and `trees.json` retain the full node mapping.

Official [iTOL annotation templates](https://itol.embl.de/help.cgi) were downloaded on 18 October 2026 and archived with retrieval timestamps and hashes in `sources/itol_templates/`. The official version history was checked against iTOL 7.6. Header options, data field widths, leaf identifiers, unique labels, Newick parsing, absence of branch lengths/support, and denominators were validated locally. The previews are local renders, not server-verified iTOL exports.

## Files and iTOL upload sequence

The root folder is the conservative main package. `extended/` contains the recorded-family package; `supplementary_genus/` contains the non-generic genus package. Each package has its own compatible tree, labels, five data tracks, and source table. **Do not mix files from different packages.**

When authorized to use iTOL, upload in this order:

1. Create a tree with `viral_family_taxonomy.nwk` (or `viral_genus_taxonomy.nwk` in the supplement).
2. Apply `itol_family_labels.txt` (or `itol_genus_labels.txt`).
3. Add `itol_taxonomy_confidence.txt`, then `itol_window_support_multibar.txt`.
4. Add `itol_positive_fraction.txt`, followed by `itol_positive_counts.txt`.
5. Optionally add `itol_evidence_multibar.txt`; this partitions classifier-associated and strict retrieval-only candidate genomes.

Use rectangular display, aligned tip labels, unscaled cladogram/ignore branch lengths, no scale bar, no support values, white background, and minimal grey hierarchy lines. Use the multibar's stacked setting (`ALIGN_FIELDS=0`), not separate columns. Counts are Both (blue #0072B2), 11-only (orange #E69F00), 15-only (purple #CC79A7). Keep fraction bars on 0–100%. Exact n/N labels use **candidate / outcome-resolved**, where N = candidates + fully analyzed negatives. The family source table also gives selected total, any-inference count, and unresolved count. N=0 has no fraction bar and an explicit `0/0 (NA)` text label; zero candidate fractions with N>0 remain genuine zeros.

Place confidence immediately after family labels, then stacked support, fraction, and n/N. Confidence colors mean non-generic rank resolved (green), mixed generic/non-generic (amber), generic-source only (vermillion). Keep unresolved and generic exclusions in an inset or accompanying table; do not attach them to an invented branch. Root conservative confidence is uniformly green because generic sources were withheld. Use 8–10 pt leaf text at final size; the 75-genus supplement needs a tall page or several panels. Export vector SVG/PDF from iTOL only after visually checking the final layout.

`figures/` contains local SVG, PDF, and PNG previews for all three variants. The proposed main figure is `Figure_B_conservative`; its inset explicitly accounts for withheld units. Extended family and genus previews retain the same count and fraction semantics. The optional classifier/retrieval track is supplied as iTOL data; quantitative values also appear in every source table.

## Suggested caption

**B. Taxonomic distribution of exploratory Saccharide-associated viral GT neighborhoods.** The family hierarchy follows cached NCBI Taxonomy (17 October 2026); edges do not represent evolutionary distance. Stacked bars count unique genome units with candidate support at both window sizes, 11 genes only, or 15 genes only. Candidates satisfy the frozen Saccharide classifier criterion or Saccharide-reference P@10 ≥ 0.5 in either size. Grey bars and n/N show candidates among outcome-resolved genome units; unresolved outcomes are excluded from N. The conservative tree displays non-generic sources with explicit family assignments, with withheld records summarized separately. All 190 genomes supported by both sizes contain a fixed baseline locus supported by both sizes (217 loci total). Frequencies apply only to the GT-A anchor-selected cohort. Embedding association does not demonstrate glycan biosynthesis or common ancestry of the candidate modules.

## Reproduction and QA

Run `scripts/build_figure_b.py` with standard-library Python in WSL from the repository. It recomputes calls from frozen predictions and reference retrieval without encoding/inference. `scripts/fetch_templates.py` downloads only public templates. `scripts/plot_figure_b.py` uses the existing Windows Matplotlib dependency bundle; it writes only here. `scripts/audit_delivery.py` uses the existing Biopython environment for an independent tree/data/annotation audit. `validation.json` records source and registry checks; `delivery_validation.json` records independent delivery checks. `artifact_checksums.tsv` covers the final downstream files except itself and the ZIP bundle. The ZIP is a convenience copy; the directory remains the authoritative reviewable source.
