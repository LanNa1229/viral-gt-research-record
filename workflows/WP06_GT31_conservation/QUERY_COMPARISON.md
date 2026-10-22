# F2Y1U4 versus GT31: three summary bars

This workflow reproduces the numerical bars in Supplementary Figure X-2-2:
sequence identity **16.42%**, query coverage **99.25%**, and alignment TM-score
(AlnTM) **0.7996**. Code and input data are included in this directory.

## Run

From the repository root, with Python 3.12 (standard library only):

```bash
python workflows/WP06_GT31_conservation/scripts/summarize_query_comparisons.py
```

The script also works from another working directory via its full path.
Use `--outdir path/to/new_folder` to choose a new or empty output directory.
The default is `workflows/WP06_GT31_conservation/query_comparison_output/`.

## Inputs and calculation

1. `inputs/reference_members.tsv` defines SC01 (5 proteins), SC02 (4), SC03 (4),
   and SC04 (4). All 17 references are compared with F2Y1U4, rather than only
   the four representatives used in the local-core analysis.
2. Sequence identity is calculated from
   `inputs/alignments/C012_GT31_reference19_plus_F2Y1U4_aa.fa`. For each reference,
   divide identical amino-acid pairs by alignment columns in which both that
   reference and F2Y1U4 are nongap. The characters `-` and `.` are gaps. The
   other two references remain in the archived MSA but do not enter these bars.
3. Query coverage and AlnTM are read from `qcov` and `alntmscore` in
   `inputs/foldseek/F2Y1U4_vs_reference17.tsv`. These are the 17 query-directed
   rows from the corrected Foldseek run, with F2Y1U4 as query. The reported
   precision is preserved. Reverse-direction scores are not averaged, and
   AlnTM is not replaced with `min(qtmscore, ttmscore)`.
4. For each metric separately, calculate the median within each reference
   subcluster. For the groups with four references, average the two middle
   observations. Then take the unweighted arithmetic mean of the four medians.
   Every subcluster contributes one quarter, despite their unequal sizes.

Identity comes from the FoldMason MSA; coverage and AlnTM come from Foldseek
pairwise comparisons. Their alignment definitions differ. These comparisons
use the full archived input models/alignments, without a core-only residue mask.

## Expected values

| Group | n | Median identity (%) | Median query coverage (%) | Median AlnTM |
|---|---:|---:|---:|---:|
| SC01: B3GALT5-like | 5 | 17.801047 | 98.70 | 0.8067 |
| SC02: Generic B3GALT-like | 4 | 15.306122 | 99.10 | 0.8070 |
| SC03: C1GALT1 | 4 | 18.702423 | 99.60 | 0.7789 |
| SC04: B3GNT3/4 | 4 | 13.861726 | 99.60 | 0.8056 |
| Mean of four medians | | **16.417830** | **99.25** | **0.79955** |

The script aggregates before rounding. Percentage labels are rounded to two
decimal places and AlnTM to four, giving 16.42%, 99.25%, and 0.7996. Decimal
arithmetic from the original observations is used for display rounding so
binary floating-point representations do not change the 0.79955 rounding tie.
AlnTM is a score, not percent sequence identity.

## Outputs and verification

- `per_reference.tsv`: all 17 comparisons, including identical-residue counts
  and the nongap denominators used for sequence identity.
- `subcluster_medians.tsv`: four groups and their three medians, before display rounding.
- `bar_summary.tsv`: three unweighted means plus figure labels and units.
- `validation.json`: success status and validation results.

The script checks the package hashes, reference membership, alignment shape,
and one-to-one query/reference comparisons. Individual identity counts and
group medians must match the archived October audit within 1e-12. Group medians
also must agree with the July table at its reported four-decimal precision
(absolute tolerance 0.00005001). Final display values must match the reported
figure values. Expected files are used only for validation, after calculation.

The saved successful run is `validation/query_comparison_validated_run.json`.
The July values, October audit values, and individual identity counts are in
`expected/query_comparison/`. Sources and extraction hashes are in
`provenance/query_comparison_sources.tsv`. This reruns postprocessing of the
archived inputs; it does not rerun Foldseek, FoldMason, or structure prediction.

Suggested legend wording:

> Bars show the unweighted mean of subcluster medians across four GT31 reference
> groups (SC01–SC04; n = 5, 4, 4 and 4 proteins): sequence identity, 16.42%; query
> coverage, 99.25%; and alignment TM-score (AlnTM), 0.7996. Sequence identity was
> calculated over mutually nongap positions in the FoldMason alignment; query
> coverage and AlnTM were obtained from Foldseek comparisons with F2Y1U4 as query.
