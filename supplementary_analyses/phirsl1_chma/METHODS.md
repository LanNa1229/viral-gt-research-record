# Reproducible methods
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
