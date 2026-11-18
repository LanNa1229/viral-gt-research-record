# Sequence- and structure-based assessment of Chimallin homologues in Ralstonia phage phiRSL1
Target genome: **NC_010811.2**. This is an archive of completed work, not a new analysis run.

The complete annotated proteome (343 proteins) was assessed by BLASTP and HMMER using six experimentally supported ChmA references and a 68-accession published homologue set, plus Chimalliviridae core-marker screening. Foldseek local comparisons and global TM-align results are retained separately.

The final structural update adds the representative ORF270 AF3 sample 1 (ranking 0.8052457233370919; source metadata rounds to 0.81). ORF270's maximum resolved-ChmA-normalized global TM-score is **0.131**. The updated dataset comprises **157 proteins and 471 comparisons**; six directed non-self global ChmA controls span **0.8582–0.9263**.

**No convincing canonical ChmA homologue was detected by sequence- or structure-based searches.** This does not establish biological absence of a phage nucleus or demonstrate replacement of physical genome protection by DNA glycosylation.

## Contents
- inputs/: native sequence references, full phiRSL1 proteome and independent NCBI translation source, alignments/profiles, six leave-one-out sets, core-marker references/profiles, 156 AlphaFold DB models, processed experimental ChmA monomers and selected ORF270 AF3 model/confidence records.
- results/frozen/: original frozen tables, raw BLAST/HMM/Foldseek outputs, controls, QC and original interpretation. These historical files correctly describe the original 156-model analysis.
- results/orf270_update/: added-model and harmonized global-control comparisons and 157-model combined table; separate from original results.
- results/execution_evidence/: retained original execution evidence; no analyses rerun for this deposit.
- scripts/: existing analysis scripts and relevant helpers.
- figures/nature_supp_chma_orf270_update/: latest approved PDF/SVG/600 dpi PNG, plotting script, legend, QC, Panel A/B source data and structural-update helper scripts. ORF270 is highlighted in both panels; panel B has no 0.5 heuristic line or label.
- METHODS.md: existing methods, with personal filesystem prefixes redacted.
- FILE_MANIFEST.tsv: archived relative paths, byte sizes and SHA-256 checksums. SOURCE_FILE_MAP.tsv links copied files to original analysis-relative paths and original/archived hashes. PATH_REDACTIONS.tsv discloses path-only redactions.
- THIRD_PARTY_NOTICES.md: applicable data attribution and AF3 output terms. EXCLUSIONS.md documents scope.

This is an archival organization rather than a newly tested executable workflow. Existing scientific code is retained without algorithmic refactoring; personal filesystem prefixes are replaced with semantic placeholders in public copies. The original local files are unchanged. Archived figure exports and numerical data are not regenerated or recalculated. Publication copies of metadata-bearing tables retain all scientific values; filesystem strings alone are redacted.

## Reference sources
Exact accessions are in results/frozen/reference_sequences.tsv and published_chma_accessions.tsv; provenance is retained in inputs/source_metadata/.
- Laughlin et al. (2022), ChmA structures and SI Table 7: https://doi.org/10.1038/s41586-022-05013-4
- Nieweglowska et al. (2023), phiPA3 lattice: https://doi.org/10.1038/s41467-023-36526-9
- Prichard et al. (2023), conserved core and RAY: https://doi.org/10.1016/j.celrep.2023.112432
- Chaikeeratisak et al. (2017): https://doi.org/10.1016/j.celrep.2017.07.064
- Mendoza et al. (2020): https://doi.org/10.1038/s41586-019-1786-y
- Abramson et al. (2024), AlphaFold 3: https://doi.org/10.1038/s41586-024-07487-w
