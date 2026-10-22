# Analysis definitions

## Representative selection

The membership table contains 17 GT31 references in SC01–SC04. References must
have at least 100 modeled residues and mean C-alpha B-factor/pLDDT of at least
70 to be eligible as representatives. The distributed inputs are predicted
models whose original B-factor fields were used as pLDDT in the archived
analysis. The driver recomputes these quality measures from the included PDBs.

For each directional Foldseek comparison, the score is
`min(qtmscore, ttmscore)`. The two directional scores for an unordered protein
pair are averaged. For each eligible protein, scores against other eligible
members of its group are summarized by their median. The representative has
the highest median, with ties resolved by mean score, observed-pair fraction,
and protein ID (ascending). Self-comparisons are excluded. The driver requires
the complete 17-by-17 input matrix and both directions for all 136 non-self
pairs; sparse network edges are not used.

This is a median-centrality representative rule, which the historical code
calls a medoid. It is not the conventional minimum-summed-distance medoid.
The archived workflow supplied representatives explicitly. Automatic selection
in this package independently recovers the same four representatives and must
match `expected/medoids.tsv` before proceeding.

## Reference alignment and local geometry

Select the 17 reference rows from the archived reference19 FoldMason AA and
3Di alignments without deleting columns, including columns that become all-gap.
AA alignment columns establish residue correspondence. Mapping aligned
sequences to the PDB residue list requires identity and coverage of at least
0.90. The reference-only alignment has 966 columns.

For a focal residue i in a reference structure R, choose other C-alpha atoms j
with `0 < distance_R(i,j) < 15 angstrom`. No sequence-separation exclusion is
applied. Only neighbors with corresponding residues in the comparison
structure T can be scored. For each such neighbor calculate

`delta_j = abs(distance_R(i,j) - distance_T(i',j'))`.

The local score is the mean of the indicators `delta_j < t`, over all matched
neighbors and the four thresholds `t = 0.5, 1, 2, 4 angstrom`. A comparison
requires at least five matched neighbors. This is a C-alpha-based local
lDDT-like score (0–1), not an all-atom lDDT implementation or pLDDT. Superposition
is not needed because internal distances are compared. Neighborhoods are
defined in the focal reference structure, so reversing a comparison can change
the score.

Within scores average valid local scores between a medoid and the other members
of its own group (4 comparisons for SC01 and 3 for each other group). Cross
scores average valid local scores between that medoid and the other three
medoids. Each comparison contributes equally; neighbor pairs are not pooled
across proteins. Missing scores are omitted from the score mean.

Neighbor coverage is the fraction of reference neighbors with mapped partners,
averaged over comparisons where the focal residue is mapped and has reference
neighbors. Focal presence is the fraction of intended comparison proteins in
which the focal alignment position maps to a residue. These measures are
retained separately from the local score.

## Core definitions

At each alignment column, both tiers require:

- all four medoids to have a mapped residue;
- residue occupancy of at least 0.70 in every reference group;
- pLDDT of at least 70 at the focal residue of every medoid.

A medoid passes the **within** test if its within local score is at least 0.75
and mean neighbor coverage is at least 0.60. It passes the **cross** test if its
cross local score is at least 0.65, mean neighbor coverage is at least 0.40,
and focal presence is at least 0.50. The focal confidence requirement also
applies in each test.

**Strict:** all four medoids pass within and all four pass cross.

**Supported:** at least three medoids pass within and at least three pass cross;
the median within score is at least 0.75 and median cross score at least 0.65.
Within and cross passing sets are counted independently and need not contain
the same three medoids. The shared occupancy and confidence requirements are
not relaxed.

Strict is contained in supported. **Supported-only** means supported minus
strict, giving 147 - 98 = 49 columns. AA and 3Di conservation are descriptive
outputs, not cutoff criteria for strict/supported membership.

## Mapping to F2Y1U4 and visualization

Core membership is defined using only the reference alignment. Each reference
core column is linked through its Q9N294 residue index to the archived
reference19-plus-F2Y1U4 alignment. The query residue in that column is then
assigned the reference class. Positions without a mapped query residue remain
unmapped. This yields 93 mapped strict and 41 mapped supported-only residues.

The query receives the reference class; it is not required to independently
pass the reference within/cross criteria. Query-versus-medoid local geometry is
computed and reported separately. The resulting class PDBs assign 100 to
strict, 60 to supported-only, and 0 to other residues in the B-factor field.
The query writer verifies that no other atom-record fields change.

The sampled GT31 references are predominantly metazoan and small in number.
Core mapping supports structural comparison within this reference set; it does
not establish universal GT31 membership criteria, substrate specificity, or
enzymatic activity. Group membership is annotation-based, not a newly inferred
phylogeny or threshold-based structural clustering.
