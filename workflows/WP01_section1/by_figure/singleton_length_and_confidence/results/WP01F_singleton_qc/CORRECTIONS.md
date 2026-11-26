# v3 corrections after reviewing the first real-data run

## 1. Nested duplicate structure files

The structure directory contained 2,038 PDB files for 1,727 cluster members. The first run marked 311 IDs as ambiguous because the same ID occurred both at the top level and under `selected_hits/`.

v3 resolves this common layout by preferring the unique shallowest path (normally `STRUCTURES_DIR/ID.pdb`). If candidates have equal path priority, byte-identical duplicates are resolved deterministically; genuinely non-identical equal-priority candidates remain ambiguous. Redundant copies associated with a cluster ID are no longer reported as structures outside the membership set.

## 2. QC thresholds

The v2 lower IQR and MAD length fences collapsed to 1 aa because the multi-member length distribution was broad and right-skewed. Similarly, the 5th-percentile pLDDT and 95th-percentile low-confidence-fraction cutoffs were too permissive. This caused 95.4% of singletons to be labeled "unexplained", which is not a defensible interpretation.

v3 uses transparent empirical length strata derived from the complete multi-member population:

- below multi-member P10
- P10 to Q1
- Q1 to Q3
- Q3 to P90
- above P90

Confidence is described using standard median-pLDDT classes:

- <50
- 50 to <70
- 70 to <90
- >=90

The main singleton composition figure now uses a two-axis classification based on multi-member P10 length and median pLDDT 70. These remain QC descriptors, not GT/non-GT assignments.

## 3. Rescue prioritization

- Tier 1: length >= multi-member P10 and median pLDDT >=70
- Tier 2: length < multi-member P10 but median pLDDT >=70
- Tier 3: length >= multi-member P10 but median pLDDT <70
- Tier 4: length < multi-member P10 and median pLDDT <70

Tier 1 is the strongest direct structural-rescue set. Tier 2 should not be discarded because compact viral proteins or isolated domains can be genuine.
