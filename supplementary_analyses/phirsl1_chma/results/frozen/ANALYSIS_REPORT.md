# Does Ralstonia phage phiRSL1 encode ChmA?
Analysis date: 4 November 2026. Target: NC_010811.2. Local investigation; nothing was pushed or published.

**Conclusion: phiRSL1 is unlikely to use a canonical ChmA-based phage nucleus (moderate confidence). No convincing ChmA homolog was detected with the stated references and thresholds. This does not establish absence of a phage nucleus or an alternative compartmentalization function.**

## Evidence
The fresh complete NCBI annotation contains **343 translated CDSs** on a 231,255 bp genome. All translations independently match the genomic sequence, and the separately retrieved NCBI CDS-protein FASTA agrees with the extracted proteome. No protein was excluded by length or annotation. The RefSeq record is provisional; translation agreement does not validate every gene call experimentally. The record places phiRSL1 in Mieseafarmvirus; this is supporting taxonomy, not experimental evidence about its replication compartment. Genome size was not used to infer nucleus formation.

| Test | Scope | Strongest result | Interpretation |
|---|---|---|---|
| BLASTP, characterized references | Six native ChmA sequences against all 343 proteins | E = 0.40; 23.60% identity over 161 columns; 23.8% reference / 13.6% target span | No convincing hit |
| BLASTP, published homolog set | 68 ChmA accession records against all proteins | E = 0.009; 25.0% identity over 88 columns; 12.7% reference / 37.9% target span | Weak short similarity |
| BLASTP sensitivity check | Six references, SEG and composition statistics off | Best E = 0.065 | No convincing hit |
| HMMER, six characterized sequences | Sensitive search of all proteins | Best E = 0.29 | No convincing hit |
| HMMER, 68 published homologs | Sensitive search of all proteins | ORF270: E = 0.0063, 77/618 profile positions (12.5%) | Weak partial lead, insufficient to assign ChmA |
| Six-frame HMM extension | 4,234 stop-delimited translations, at least 30 aa | Best sequence E = 0.27, same ORF270 region | No significant additional coding region |
| Foldseek, local 3Di+AA | Three experimental monomers x 156 AlphaFold models | Best E = 0.1715; query TM = 0.0556 | No convincing ChmA match |
| Foldseek, global TM-align | All 468 reference/model pairs | Maximum ChmA-normalized TM = 0.277 | No convincing shared full ChmA fold |
| Chimalliviridae-specific core panel | All 21 starred groups in Prichard et al. Table S3 | 0/21 at E <= 1e-5 in BLASTP or HMMER | No family-specific gene-content support |

The operational sequence threshold was E <= 1e-5; coverage and context were also reviewed. HMM E-values are full-sequence values unless explicitly labeled domain E-values. E-values depend on the query/profile and database size and are not probabilities of biological absence.

## Reference set
The six references have experimental evidence of shell localization, nucleus formation, or shell self-assembly in the cited studies. They do not all have equivalent biochemical validation or atomic structures. The 68-record expansion is a **published homology set**, not 68 experimentally verified nucleus-forming phages. It includes Ralstonia phages RSL2, RSF1, RP12, and RP31; these must not be confused with RSL1.

| Phage | Protein | Native query accession | Length | Experimental atomic structure used |
|---|---|---|---:|---|
| 201phi2-1 | gp105 | B3FIW8 | 631 | 7SQR |
| Goslar | gp189 | A0A482GDX1 | 631 | 7SQV |
| phiPA3 | gp53 | F8SJT5 | 602 | 8FNE |
| phiKZ | gp54 | NP_803620.1 | 626 | None used |
| RAY | gp222 | ANH52002.1 | 711 | None used |
| PCH45 | gp202 | QFP93061.1 | 659 | None used |

Native sequence/PDB provenance and exact sequence matches to the published accession set are in [reference_sequences.tsv](results/reference_sequences.tsv). Full source FASTAs, mmCIFs, and original paper supplements are retained.

All six native references were detected when each was withheld from training: leave-one-out E-values range from 1.8e-24 to 1.5e-198. This validates sensitivity to the sampled known family, not arbitrarily divergent proteins. Non-identical experimental ChmA structures matched with E-values below 1e-20 and ChmA-normalized TM-scores of 0.815-0.924 in positive controls.

## Proteins retained for review
No protein qualifies as a supported ChmA candidate. The leading weak observations are retained explicitly:

- **RSL1_ORF270 / YP_001950147.1 / B2ZY82, 265 aa:** strongest expanded-profile lead. HMM positions 18-94/618 align to protein residues 156-248. Full E = 0.0063, domain independent E = 0.012, full score = 9.4 bits. Its best BLASTP match in the 68-record set is to phiKZ ChmA: 24.43% identity over 131 columns, 17.9% reference and 45.3% target span, E = 6.5. No AlphaFold DB model was available. This may be chance similarity or a remotely related fragment/domain; size alone cannot exclude an unusual or split shell. Neighboring annotated ORFs provide no second significant ChmA profile hit. The six-frame extension also found no significant additional fragment.
- **RSL1_ORF229 / YP_001950108.1, 232 aa:** best expanded BLASTP hit, to STsAS ChmA AWN08967.1, 25.0% identity over 88 columns, E = 0.009. Its best global ChmA-normalized TM is only 0.148.
- **RSL1_ORF095 / YP_001949970.1 / B2ZXW4, 111 aa:** lowest local Foldseek E = 0.1715 against 7SQR; 21.1% of resolved ChmA and 82.0% of the short target are covered. ChmA-/target-normalized TM = 0.0556/0.1787; RMSD = 17.88 angstrom; lDDT = 0.244. High target coverage reflects a short protein rather than a shell-sized fold.
- **RSL1_ORF196 / YP_001950075.1 / B2ZY50, 372 aa:** highest global ChmA-normalized TM among uncharacterized models: 0.246 (target-normalized 0.318). No sequence/profile support.
- **RSL1_ORF333 / YP_001950209.1 / B2ZYK3, 882 aa:** highest global score overall: ChmA-/target-normalized TM = 0.277/0.184; annotated DNA topoisomerase II large subunit. This does not support a ChmA assignment.

Detailed identities, coverage, E-values, coordinates, and structural scores are in [candidate_review.tsv](results/candidate_review.tsv), [hmm_chma68_domains.tsv](results/hmm_chma68_domains.tsv), and the complete BLAST/Foldseek outputs. HMMs have no single pairwise sequence identity; the table reports the separate BLAST identity where available.

## Phylogenetic support
The source table contains 68 gene groups; duplications/splits produce the commonly cited 72-gene core for RAY. We recovered all 21 starred family-specific groups and every listed reference protein for those groups from Goslar, phiKZ, and RAY. The broader mapping recovered 211/215 references. Four current phiKZ proteins could not be unambiguously mapped (cg8/gp65, cg15/gp73, cg33/gp118, cg43/gp164); these are outside the starred set and were left unresolved.

All 20 specific groups besides ChmA were negative at E <= 1e-5. The lowest core-profile E was 0.00059 (cg27 to YP_001949904.1), above the threshold; it remains in the marker table. A strong broader-core match was a lysozyme domain: RAY cg55 to phiRSL1 endolysin YP_001950059.1, E = 1.36e-22, 32.65% identity over 196 columns, 8.3% reference and 83.1% target span. A shared lysis domain does not establish Chimalliviridae-specific replication.

These are **gene-content and existing-taxonomy arguments**, not a newly inferred phylogenetic tree. No confident ChmA/family-specific ortholog set was available for a defensible concatenated placement. PhuZ was not used as a necessary marker: experimentally nucleus-forming Asesino lacks it.

## Limits and confidence
- Only **156/343 proteins (45.5%)** had available AlphaFold DB models, including 135 labeled uncharacterized by UniProt. Every model sequence equals its full target protein, but sequence completeness is not structural accuracy.
- **70/156 models have mean pLDDT below 70**; median mean pLDDT is 71.44. The structural negative is weaker than a complete, well-folded structural proteome.
- ORF270 has no model. A highly divergent ChmA, split gene, missed circular-origin ORF, nonhomologous replacement, or another compartment remains possible. The six-frame extension does not join circular-origin translations.
- Profiles and markers sample known Chimalliviridae. The 21 core groups are correlated evolutionary evidence, not 21 independent statistical tests.
- No phiRSL1 infection imaging, candidate localization, biochemical assembly, or genetics was performed or retrieved that directly establishes presence or absence of a nucleus.

**Confidence is high in the bounded result that these searches detected no convincing canonical ChmA homolog; confidence is moderate in the biological inference against the canonical ChmA-based phage nucleus. Biological absence is not established.** Useful next evidence would be reliable models of ORF270 and the unmodeled proteins, then infection microscopy and candidate localization/functional tests.

## Reproduction and evidence files
[METHODS.md](METHODS.md) describes execution order, parameters, normalization, source extraction, and versions. [validation.tsv](results/validation.tsv) records the final checks. Source accessions, retrieval URLs/times, raw outputs, executable scripts, input/output manifests, and SHA-256 checksums are retained. Some earlier exploratory five-reference outputs remain; authoritative final outputs use the six-reference chma_characterized, the chma68 profile, and the final core_cg* tables.

Primary studies:
- [Laughlin et al. 2022: ChmA architecture and self-assembly](https://doi.org/10.1038/s41586-022-05013-4), including SI Table 7.
- [Nieweglowska et al. 2023: phiPA3 shell lattice](https://doi.org/10.1038/s41467-023-36526-9).
- [Chaikeeratisak et al. 2017: conserved Pseudomonas phage nucleus](https://doi.org/10.1016/j.celrep.2017.07.064).
- [Mendoza et al. 2020: PCH45 genome shielding](https://doi.org/10.1038/s41586-019-1786-y).
- [Prichard et al. 2023: core genome and RAY](https://doi.org/10.1016/j.celrep.2023.112432), especially Table S3.
- [Prichard et al. 2025: Asesino nucleus without PhuZ](https://doi.org/10.1038/s41598-024-64095-4).
