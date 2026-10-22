# Provenance and scope

Packaged on 22 October 2026 from `C012_GT31_research_record_20261021`, which
preserves the July 2026 C012 GT31 workflow and its October recovery/validation.
The packaging date is distinct from the original analysis date.

## Preserved code and data

The five Python files under `src/vgt_constraints/` are byte-identical to the
archived modules. `scripts/make_query_core_classes.py` is the unchanged helper
added in the 21 October research record. `provenance/code_source_manifest.tsv`
records the archived code/configuration paths and SHA-256 hashes. The separate
input bundle contains `provenance/source_manifest.tsv`, with source paths,
hashes, and transformations for the data as well.

The reference membership table retains its original fields and rows, with
structure paths made relative to this package. Its source archive hash is
recorded separately. In the separate input bundle, the 18 original PDBs
(17 references plus F2Y1U4) retain
their original coordinates and confidence values. These are the preserved
analysis inputs; this package does not claim to rerun their model generation
or provide a new structure-database release/version assignment.

The 289-row reference Foldseek table is the complete 17-reference subset of
the corrected 22-by-22 run. The original incorrect/incomplete run is not used.
Raw alignments and expected analysis outputs are copied without modification.

## New packaging code

`scripts/run_analysis.py` performs automatic representative selection, derives
the 17-row reference alignments from the original 19-row alignments, runs the
original pan-core entry point, writes class PDBs, and compares the numerical
outputs with archived tables. It does not execute the separate exploratory
constraint-analysis entry point. Quality filtering, local geometry, core
thresholds, and query mapping use the preserved implementations.

`scripts/make_pymol_pml.py` embeds both class PDBs in one PML to remove external
file-path dependencies. It preserves the dark teal/light blue/grey category
colors. It does not reconstruct the historical structural superposition,
camera angle, or final manuscript layout.

## Upstream alignment record

The retained FoldMason logs identify the software build as
`9a71f1dc2eaa5f6ac5ec79dde0bc614a9dfc2057` (reported under the log's
`MMseqs Version` field). The recorded commands used `easy-msa`,
`--report-mode 1 --refine-iters 0 --pair-threshold 0.0 --threads 16`.
One alignment used 19 references; the other used those references plus F2Y1U4.
The original explicit structure input order is preserved in the input bundle's
`provenance/foldmason_commands.txt`.

Only 17 of those references participate in core definition. The two additional
references remain in the archived alignments. Their structures are not needed
for the distributed postprocessing run and are not included. Rebuilding the
upstream alignments from scratch is outside this package's validated scope;
the exact archived alignments are the reproducibility inputs.

## Validation

The driver verifies package checksums and compares all rows/columns in the
reference core table, medoid local-metric table, and query mapping table with
the archived results. Numeric comparisons use absolute and relative tolerance
of 1e-10; categorical assignments must match. It also checks all summary
values, representative identities, PDB class counts, and preservation of query
atom-record fields outside the B-factor column.

The input bundle's validation record was produced by executing this package
with the listed Python dependency versions. It is not a rerun of Foldseek,
FoldMason, structure prediction, or the PyMOL renderer. The checksum manifests
exclude the validation record itself to avoid a self-referential checksum cycle.
