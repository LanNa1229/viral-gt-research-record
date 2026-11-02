# Recovered scientific rules

These describe the recovered code, not newly chosen thresholds. Authoritative implementation: `originals/260630_Temp/build_cluster8_neighborhood_matrix.py`.

## Annotation stage versus heatmap stage

The recovered `analyze_foldseek_results_v2.py` names hits from existing target descriptions/IDs, truncates readable names to 300 characters, attaches name/keyword triage tags, and writes the top 20 rows per query/database by rank by default. The saved per-locus reports state high confidence as probability ≥0.70, E-value ≤0.001, query coverage ≥0.50, and alignment length ≥80. That Boolean is distinct from the heatmap's 0–3 score; the matrix builder does not read it. Saved annotation reports provide run-specific priority databases; the script default of `pdb100` must not be assumed to describe every run.

## Heatmap evidence score

The first matching branch applies, followed by database/rank adjustments:

| Initial score | Exact criteria |
|---|---|
| 3 | probability ≥0.90 AND E-value ≤1e−5 AND query coverage ≥0.50 AND alignment length ≥100 |
| 2 | probability ≥0.70 AND E-value ≤1e−3 AND query coverage ≥0.35 AND alignment length ≥80 |
| 1 | E-value ≤1e−2 AND alignment length ≥50; probability and coverage need not pass any cutoff |
| 0 | none of the above |

If the database string contains `cath`, score 3 is capped at 2. If rank >20 and score >1, subtract 1. Target coverage is read but unused. Sequence identity and Foldseek raw score are not scoring criteria. Missing numeric values are converted by `safe_float` to its default 0; an absent E-value field uses 999 but a present NaN becomes 0. This is a historical edge case, not a corrected policy.

The classifier searches lowercase concatenated `hit_name`, `target_id`, `matched_name_terms`, and `matched_keyword_categories`. It allows multiple classes per hit. It ignores score-0 hits and sets each locus/class cell to the **maximum** score among all class-matching hits over all queries and databases. The cell is not a gene count, not a probability of function, and not a vote across independent databases. A zero means no qualifying named match under these rules; it does not establish gene absence.

## Row classification, in order

1. C5_MT ≥2 AND TET_5mYOX ≥2 → `C5-MT + 5mYOX-associated`.
2. dCMP_hydroxymethylase ≥2 → `dCMP hydroxymethylase-associated`.
3. Any of ThyX, ThyS, dCMP_dTMP_enzyme, nucleotidyltransferase ≥2 → `ThyX/ThyS / nucleotide-metabolism-associated`.
4. Any of additional_GT, NDP_sugar_enzyme, epimerase, dehydratase, reductase, aminotransferase, methyltransferase, acetyltransferase, sulfotransferase, APS_PAPS_kinase ≥2 → `additional-GT/glycan-associated`.
5. Otherwise → `ambiguous/no DNA-context`.

Standalone TET support does not trigger a dedicated historical row class. Cluster8_GT_like and TET_associated_GT_like are excluded from the accessory-glycan branch. Row grouping and inferred_model use different criteria.

## Module totals and inferred model

DNA_context_score sums C5_MT, TET_5mYOX, dCMP_hydroxymethylase, hmC_hmdU_related, ThyX, ThyS, dCMP_dTMP_enzyme and the six DNA-metabolism classes. It excludes nucleotidyltransferase even though that class can select the nucleotide-metabolism row group. Glycan score sums the 12 glycan-pathway classes; virion score sums five virion classes. Overlapping class assignment means these are not counts of independent biological activities.

For inferred_model, DNA precursor means any of C5_MT/TET_5mYOX/dCMP_hydroxymethylase/hmC_hmdU_related/ThyX/ThyS ≥2; glycan means any of cluster8_GT_like/TET_associated_GT_like/additional_GT/NDP_sugar_enzyme ≥2; virion means capsid/tail/portal/tape_measure ≥2 (terminase excluded). The ordered outputs are DNA+glycan → hybrid; DNA → DNA-associated; glycan+virion → possible virion/capsid glycosylation; glycan → glycan/no-DNA; otherwise ambiguous. Broad dCMP_dTMP_enzyme or nucleotidyltransferase alone is excluded from this DNA-precursor definition.

## Figure mapping

The recovered clean plot displays 33 classes in six modules. Rows sort by primary-group rank ascending; DNA_context_score, glycan_score and virion_score descending; candidate label ascending. Scores are filled with 0 and clipped to [0,3] for plotting. No matrix class column was missing or required that fallback in validation. The four fixed colors map 0,1,2,3 to `#F2F2F2`, `#C6DBEF`, `#6BAED6`, `#08519C`. The title-enabled invocation reproduces the original; `--no-title` does not.
