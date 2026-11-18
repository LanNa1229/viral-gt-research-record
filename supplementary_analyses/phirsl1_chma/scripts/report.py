#!/usr/bin/env python3
"""Generate the evidence report and source/input/output manifests."""
from pipeline import *
import shutil
def tsv(path):return list(csv.DictReader((ROOT/path).open(),delimiter='\t'))
def write(name,text):(ROOT/name).write_text(text,encoding='utf-8')
def main():
 ref=tsv('results/reference_sequences.tsv')
 write('README.md',"""# Does Ralstonia phage phiRSL1 encode ChmA?
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
"""+''.join(f"| {r['phage']} | {r['protein']} | {r['query_accession']} | {r['length_aa']} | {r['experimental_atomic_pdb'] or 'None used'} |\n" for r in ref)+"""
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
""")
 write('METHODS.md',"""# Reproducible methods
All results are local. No GitHub push, publication, external prediction submission, or destructive cleanup was performed.

## Environment
WSL Ubuntu; Python 3.11.15; Biopython 1.87; BLAST+ 2.12.0+; HMMER 3.4; MUSCLE 5.3; Foldseek hash 755d3fa84b8259d4d30aa92782f0bfc2bfc3363f. Exact paths/system information are in results/environment.json and logs/version_*.stdout. PDF extraction used bundled Windows Python/pypdf; version, page counts, and hashes are in logs/pdf_extraction.json.

An equivalent environment needs Biopython 1.87, pypdf (recorded version), BLAST+, HMMER, MUSCLE, and Foldseek. The two executable constants at the top of scripts/searches.py currently use the installed WSL paths; edit them for another system. BLAST+ and HMMER are found on PATH.

## Input and reference preparation
NCBI EFetch retrieved NC_010811.2 as gbwithparts and independently as fasta_cds_aa. All 343 annotated translations were extracted without length/annotation filtering and verified by genomic translation using annotated table/codon start. Alternative bacterial initiation codons may produce methionine.

PDB entities link native UniProt B3FIW8, A0A482GDX1, and F8SJT5. Native phiKZ NP_803620.1 and RAY ANH52002.1 came from source genomes; PCH45 QFP93061.1 was fetched directly. The 68 accession records from Laughlin et al. 2022 SI Table 7 were all retrieved by protein EFetch and retained as published (including related/identical records). HMMER applies default sequence weighting.

Prichard et al. Table S3 was extracted from the original publisher PDF. PDF pages 5-7 were rendered and visually inspected to verify the columns and all 21 asterisks. Mapping uses exact current/old locus tags, with source genomes NC_048170.1, NC_004629.1, and KU886224.1. The parser checks group cardinalities.

Native phiKZ gp54 is PHIKZ054 / NP_803620.1, 626 aa. The current PHIKZ_p54 is a different 319-aa PHIKZ173 protein and is excluded from final references. Goslar current tags HOV27_gp* map through old Goslar_* tags; RAY labels below 100 have inconsistent padding. Four unmatched general-core phiKZ proteins remain unresolved and are not substituted.

## Sequence analysis
BLASTP searches six characterized references and 68 published homologs against all 343 proteins: reporting E cutoff 100, max_target_seqs 500, four threads, otherwise default protein scoring/SEG/composition statistics. The six-reference sensitivity check uses -seg no -comp_based_stats 0. Raw and headered tables include identity, alignment length, coordinates/lengths, E-values, bit scores, and qcovhsp.

MUSCLE -align supplies alignments for hmmbuild --amino. hmmsearch --max --cpu 4 -E 100 --domE 100 disables heuristic filtering and retains weak results. Separate six-reference and 68-reference profiles search the complete proteome. Each native reference was also withheld, an alignment/profile rebuilt without it, and the held-out sequence searched.

All 68 core groups were screened by BLASTP using 211 mapped references. Profiles were built/searched for all 21 starred groups, with every listed specific-group reference mapped. The source has no phiKZ member for cg14. All 20 specific groups besides cg4/ChmA provide the additional marker evidence.

The operational sequence threshold is full-sequence E <= 1e-5. A threshold-passing match would still need coverage/context review. HMM coverage is aligned HMM span divided by profile length; target coverage is aligned target span divided by target length. BLAST identity uses alignment columns; coverage uses positions on each sequence, not gapped alignment columns. E-values across profiles, queries, and different search database sizes are not directly calibrated biological probabilities.

An annotation sensitivity extension translated both DNA strands in all three frames, split at stops, retaining peptides at least 30 aa without start-codon requirements (4,234 peptides). The chma68 HMM was searched with the same sensitive settings. This does not join circular-origin ORFs or demonstrate translation in vivo.

## Structures
UniProt returned 343 phiRSL1 sequences. Exact sequence matching to the target CDS preceded each AlphaFold DB request. Available models were downloaded; compressed payloads were detected and decompressed while preserving originals. All 156 coordinate-model sequences were independently checked against full targets. Inventory/errors explicitly identify unavailable models. No new structure was predicted.

Experimental references are chain A of 7SQR, 7SQV, and 8FNE. In 8FNE only label_seq_id >=395 is retained to remove the MBP fusion. Monomer files preserve coordinates, use native-relative residue numbers, and assign uniform B-factor 50. Resolved chain lengths are 522/523/525 residues: structural coverage/TM normalization refers to these coordinate chains, not full native 631/631/602-aa sequences.

Foldseek receives explicit model file lists. Both searches exhaust all 3 x 156 pairs with --exhaustive-search 1, four threads, -e 100, --max-seqs 1000, and -a 1. Local search uses default 3Di+AA alignment; global search uses --alignment-type 1 (TM-align). Conversion preserves identity fraction, alignment length, coordinates, query/target lengths and coverage, E-value, bits, alignment-/query-/target-normalized TM, RMSD, and lDDT. The three known ChmA folds were cross-searched as controls.

Local Foldseek E-values and globally optimized normalized TM-scores are distinct evidence. A TM-score around/above 0.5 is a common fold heuristic, not a ChmA functional test (https://doi.org/10.1093/nar/gki524). The strong known-ChmA controls offer an empirical comparison. Poor pLDDT, interdomain orientation, missing models, and divergent proteins limit exclusion. pLDDT is summarized from CA B-factors without excluding low-confidence models.

## Execution order
Run from the repository root using the recorded WSL Python or equivalent environment. Original argv, status, times, stdout/stderr are in logs/commands.jsonl and paired logs. Download URLs, UTC times, sizes, and hashes are in logs/downloads.jsonl; optional failures are retained in logs/download_errors.jsonl.

~~~bash
python -B analysis/phirsl1_chimallin_20261104/scripts/pipeline.py acquire
python -B analysis/phirsl1_chimallin_20261104/scripts/pipeline.py inventory
python -B analysis/phirsl1_chimallin_20261104/scripts/references.py
python -B analysis/phirsl1_chimallin_20261104/scripts/supplements.py
python -B analysis/phirsl1_chimallin_20261104/scripts/extract_pdf_text.py --bootstrap
python -B analysis/phirsl1_chimallin_20261104/scripts/extended_references.py
python -B analysis/phirsl1_chimallin_20261104/scripts/extract_pdf_text.py
python -B analysis/phirsl1_chimallin_20261104/scripts/searches.py sequence
python -B analysis/phirsl1_chimallin_20261104/scripts/searches.py expanded
python -B analysis/phirsl1_chimallin_20261104/scripts/searches.py structure
python -B analysis/phirsl1_chimallin_20261104/scripts/searches.py structure_tm
python -B analysis/phirsl1_chimallin_20261104/scripts/markers.py
python -B analysis/phirsl1_chimallin_20261104/scripts/additional_checks.py
python -B analysis/phirsl1_chimallin_20261104/scripts/summarize.py
python -B analysis/phirsl1_chimallin_20261104/scripts/report.py
python -B analysis/phirsl1_chimallin_20261104/scripts/pipeline.py checksums
~~~

The bootstrap PDF step extracts the 2022 accession table while permitting an unavailable 2023 PDF. The extension script then uses that accession list and downloads the original 2023 publisher PDF/CSV/XLSX. The NIHMS supplementary ZIP was unavailable through Europe PMC, so the publisher route is recorded. Actual extraction used bundled Windows Python/pypdf because WSL lacked pypdf. The extraction script also accepts an explicit analysis root path.

For source visual verification, bundled Poppler pdftoppm rendered pages 5-7 with -f 5 -l 7 -scale-to 1700 -png. Images and the visual validation note are retained in results/.

## Preserved outputs
README.md is the conclusion; summary.json and headered TSVs provide inspectable metrics. candidate_review.tsv joins weak leads and never promotes them to homologs. Full raw BLAST, HMM tbl/domtbl, Foldseek, mappings, model QC, controls, references, and validation remain available.

Earlier five-reference files and failed-attempt logs remain for provenance; final outputs use chma_characterized (six sequences), chma68, and final core_cg* files. work/ contains regenerable search databases, intermediates, and renderings, excluded from stable checksums. Input/output manifests and checksums describe retained stable files. No files outside this analysis directory were modified.
""")
 write('config.json',json.dumps(dict(target='NC_010811.2',sequence_evalue_threshold=1e-5,reporting_evalue=100,threads=4,proteome_length_filter=None,experimental_chma_pdbs=['7SQR','7SQV','8FNE'],specific_core_groups=21,sixframe_min_length=30,global_tm_similarity_heuristic=0.5),indent=2))
 catalog=[]
 for line in (ROOT/'logs/downloads.jsonl').read_text().splitlines():
  r=json.loads(line);p=ROOT/r['path'];preserved=p;transform='none'
  if r['path'].startswith('inputs/afdb/') and (ROOT/'inputs/afdb_compressed'/(p.name+'.gz')).exists():
   preserved=ROOT/'inputs/afdb_compressed'/(p.name+'.gz');transform='gzip decompression; original retained'
  catalog.append(dict(path=r['path'],url=r['url'],utc=r['utc'],download_bytes=r['bytes'],download_sha256=r['sha256'],preserved_download_path=str(preserved.relative_to(ROOT)),preserved_download_matches=(digest(preserved.read_bytes())==r['sha256']),current_path_bytes=p.stat().st_size,current_path_sha256=digest(p.read_bytes()),transformation=transform))
 table('sources/source_catalog.tsv',catalog)
 for n in ['05','06','07']:
  p=ROOT/f'work/core_table-{n}.png'
  if p.exists():shutil.copy2(p,ROOT/f'results/source_core_table_page{n}.png')
 write('results/source_table_visual_validation.txt','Prichard 2023 Table S3, publisher PDF pages 5-7, rendered with Poppler and visually inspected on 2026-11-04. Column assignments and all 21 asterisks agree with the extracted 68-group table. cg27 duplication and cg42 multiple copies retained; phiKZ cg14 n/a is source data. Font warnings did not obscure the rendered table.\n')
 for out,roots in [('input_manifest.tsv',['inputs','sources']),('output_manifest.tsv',['results'])]:
  rows=[]
  for r in roots:
   for p in sorted((ROOT/r).rglob('*')):
    if p.is_file():rows.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=digest(p.read_bytes())))
  table(out,rows)
 print('Report, methods, source catalog and manifests saved',flush=True)
if __name__=='__main__':main()
