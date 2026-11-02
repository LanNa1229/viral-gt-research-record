"""Assemble an unpublished WP08 candidate, exclusively outside the existing repository."""
import csv, hashlib, json, shutil, zipfile
from pathlib import Path

BASE=Path(__file__).resolve().parent
STAGE=BASE/'WP08_recovery_staging'
RELEASE=BASE/'WP08_release_candidate'
PKG=RELEASE/'workflows/WP08_cluster8_gene_neighborhood_heatmap'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,text):
    p=PKG/name; p.parent.mkdir(parents=True,exist_ok=True); p.write_text(text.strip()+'\n',encoding='utf-8')
def dump(name,obj): write(name,json.dumps(obj,indent=2,ensure_ascii=False))

def main():
    PKG.mkdir(parents=True,exist_ok=False)
    shutil.copytree(STAGE/'originals',PKG/'originals',copy_function=shutil.copy2)
    shutil.copytree(STAGE/'audit',PKG/'audit',copy_function=shutil.copy2)
    shutil.copy2(STAGE/'SHA256SUMS-originals.txt',PKG/'SHA256SUMS-originals.txt')
    shutil.copy2(STAGE/'validate_recovery.py',PKG/'reproduce.py')
    shutil.copy2(STAGE/'correct_interpretation.py',PKG/'correct_interpretation.py')
    shutil.copy2(BASE/'recover_wp08.py',PKG/'audit/recover_wp08.py')
    shutil.copy2(__file__,PKG/'audit/build_wp08_package.py')
    validation=PKG/'validation'; validation.mkdir()
    for n in ['validation.json','input_audit.json','ambiguous_multiclass_ORFs.json','A0A6J5M498_TET_evidence.tsv','plot_with_title.log','plot_no_title.log']:
        shutil.copy2(STAGE/'validation_linux_final'/n,validation/n)
    shutil.copy2(STAGE/'validation_initial/validation.json',validation/'windows_comparison.json')
    shutil.copytree(STAGE/'validation_linux_final/plot_with_title',validation/'reproduced_figure')
    shutil.copytree(STAGE/'corrected_interpretation_linux',PKG/'corrected_interpretation')
    v=json.loads((validation/'validation.json').read_text())
    provenance=json.loads((PKG/'audit/provenance.json').read_text())
    source_check=json.loads((PKG/'audit/source_preservation_final.json').read_text())
    assert len(provenance['files'])==346 and len(source_check)==346
    assert all(r['sha256_unchanged'] and r['copy_matches'] and r['mtime_unchanged'] for r in source_check)
    for r in provenance['files']: assert sha(PKG/r['archive_path'])==r['sha256']
    write('.gitattributes','* -text\n')
    write('README.md', '''
# WP08: recovered Cluster 8 gene-neighborhood heatmap

The completed June 30, 2026 analysis was recovered from `/home/rahnn/260630_Temp/` on October 29, 2026. All 11 locus annotation tables, the fixed manifest, matrix, evidence table, builder, original and clean plotting scripts, and original PDF/PNG were recovered. Related original Foldseek batch/annotation scripts and query manifests were also preserved. No Foldseek search was rerun and no missing scientific code was reconstructed.

**Validation: historical data and figure reproduced.** All 363 gene-class cells, row metadata, module scores and 4,919 evidence records agree. The clean PNG is byte-for-byte identical. The PDF is identical after ignoring only `/CreationDate` in the comparison; neither archived PDF nor its metadata was edited. Original source and archived-copy SHA-256 values and source mtimes match for 346 files.

**Interpretation requires caveats.** Reproducibility does not establish biochemical function. Read [known limitations](docs/KNOWN_LIMITATIONS.md) before using the figure. A separate review interpretation changes only grouping metadata for two TET-positive loci; the historical results remain authoritative for historical reproduction.

## Files

- `originals/260630_Temp/`: byte-preserved completed analysis, including `cluster8_global_neighborhood_panelB/figures_clean/`.
- `originals/related_scripts/`: original batch/annotation versions and related manifests, retained under their source-directory names.
- `audit/provenance.json`, `SHA256SUMS-originals.txt`: absolute source paths, timestamps, sizes and SHA-256 values.
- `validation/`: executed Ubuntu validation and matching reproduced figure; a Windows comparison documents rendering differences under Matplotlib 3.9.1.
- `docs/`: Methods, confidence rules, gene assignments, limitations, provenance and recovery report.
- `corrected_interpretation/`: explicitly separate interpretation-only matrix, figures and differences.
- `release/`: draft GitHub PR/release text, Zenodo metadata proposal and release checklist.

## Reproduce offline

Run from this workflow directory with Python, pandas, NumPy, Matplotlib and Pillow installed:

```bash
python -B reproduce.py --out /absolute/path/to/new_validation
```

The output directory must be new. This validates the archived inputs, imports the recovered builder, compares the complete matrix and evidence multiset, independently aggregates saved evidence maxima, and runs the recovered clean plotting script in a separate process. It does not use historical absolute input paths and does not invoke the Foldseek API. Originals are hash-checked again after execution.

For the exact-render check on Linux with the matching environment:

```bash
python -B reproduce.py --out /absolute/path/to/new_validation --require-exact-render
```

See [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the tested versions and the original CLI. Numerical data validation can pass under another plotting environment even when raster bytes differ. Strict mode additionally requires the original PNG bytes and PDF equality except for CreationDate.

To reproduce the separate grouping interpretation:

```bash
MPLBACKEND=Agg python -B correct_interpretation.py --out /absolute/path/to/new_corrected_interpretation
```

## Release status

This is a local, unpublished release candidate intended for `LanNa1229/viral-gt-research-record`, under `workflows/WP08_cluster8_gene_neighborhood_heatmap/`. The existing repository was inspected read-only and remains unchanged. No branch, commit, push, PR, GitHub release, Zenodo draft or publication was created. Metadata fields requiring author/rights review are identified in `release/`. This package does not grant a new license to original or third-party material.
''')
    write('docs/RECOVERY_REPORT.md',f'''
# Recovery and validation report — October 29, 2026

The requested completed analysis is recovered, rather than reconstructed.

| Requested item | Recovery location under `originals/` |
|---|---|
| `build_cluster8_neighborhood_matrix.py` | `260630_Temp/` |
| `plot_cluster8_panelB_heatmap_clean.py` | `260630_Temp/` |
| `foldseek_batch_api_v2.py` | `related_scripts/`, multiple source directories |
| `analyze_foldseek_results_v2.py` | `related_scripts/`, multiple source directories |
| `cluster8_analysis_manifest.fixed.tsv` | `260630_Temp/` |
| `cluster8_neighborhood_matrix.tsv` | `260630_Temp/cluster8_global_neighborhood_panelB/` |
| `cluster8_neighborhood_evidence_long.tsv` | same directory |
| original clean PDF and PNG | same directory, `figures_clean/` |
| 11 `10_database_top_hit_names.tsv` inputs | `260630_Temp/<candidate_id>/final/analysis_v2/` |

The non-clean original plotting script and figures, earlier malformed manifest, raw combined hit tables and saved annotation reports were also preserved. Recovery searched the supplied Ubuntu directory, related immediate Foldseek batch directories and the current Windows workspace. It was a targeted local recovery, not an exhaustive disk-forensic search. All essential files were found, so neither code reconstruction nor Foldseek execution was necessary.

## Executed checks

- 11 unique manifest candidates, 33 gene classes, 363 score cells.
- Complete matrix equality including historical `analysis_dir`, row order, scores, grouping, inferred model and group rank.
- All 4,919 evidence records reproduce as a multiset: no missing or extra records. Original set iteration makes evidence row ordering nondeterministic.
- Independent maximum aggregation of the saved evidence exactly matches every gene-class cell.
- Original PNG and reproduced PNG have 4144 × 1958 pixels and identical bytes: SHA-256 `{v['plots'][0]['sha256']}`.
- PDF comparison: identical bytes after removing only the CreationDate value for comparison. Original recorded date: June 30, 2026 at 12:29:46 UTC−04:00. The original PDF is preserved untouched.
- All 346 archived files match their source SHA-256 values; a final source check also confirms original source content and mtimes are unchanged.
- Five exact manifest-query mismatches were found and retained. See `validation/input_audit.json` and limitations.
- Original PNG visually inspected: all 11 labels, 33 columns, legend and score colorbar present. Existing layout crowding near headers and row labels is retained as historical figure content, not silently redesigned.

Historical group counts: 4 dCMP-hydroxymethylase-associated, 5 broad nucleotide-metabolism-associated, 2 additional-GT/glycan-associated. No locus is in the combined C5-MT + 5mYOX class. No locus has a nonzero ThyX score; the short legend label “ThyX/ThyS” therefore cannot be read as proof of ThyX.

Assessment: ready for historical reproducibility review with explicit biological caveats. Release metadata/rights remain a publication prerequisite. Repository source context was read at local HEAD `52dd6147505ce3434a688f7703ed76baec731c3d`; the remote HEAD was not fetched or changed. Existing unrelated working changes are recorded in `audit/repository_discovery.json`.
''')
    write('docs/METHODS_AND_SCORING.md','''
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
''')
    rules=[
    ('C5_MT','C5/cytosine/5-methylcytosine methyltransferase text patterns.'),
    ('TET_5mYOX','5mYOX, 5-methylpyrimidine dioxygenase, tet-like/tet family, 2OG/Fe or 2-oxoglutarate dioxygenase; excludes TET-associated glycosyltransferase text. Does not comprehensively recognize all TET/JBP names.'),
    ('dCMP_hydroxymethylase','dCMP hydroxymethylase, thymidylate synthase/dCMP hydroxymethylase or dUTPase/hydroxymethyl patterns.'),
    ('hmC_hmdU_related','5hmC/5hmdU, hydroxymethylcytosine/uracil/DNA, hmdUTP/hmdCTP.'),
    ('ThyX','thyx word or flavin-dependent thymidylate synthase.'),('ThyS','thys word or thymidylate synthase; can overlap ThyX and dCMP-HMase.'),
    ('dCMP_dTMP_enzyme','dCMP/dTMP/dUMP/dUTP/dCTP/dTTP, pyrimidine, thymidine, deoxycytidine or deoxyuridine.'),
    ('nucleotidyltransferase','Nucleotidyltransferase, pyrophosphorylase, nucleoside-diphosphate, NTP transferase or sugar-1-phosphate nucleotidyltransferase.'),
    ('DNA_polymerase','DNA polymerase, polymerase III/I, family A/B polymerase.'),('helicase','Helicase, DnaB, UvrD or RecA-like helicase.'),
    ('ligase','Any ligase pattern; not limited to DNA ligase.'),('nuclease','Nuclease, endo/exonuclease, RNase H, restriction endonuclease, HNH.'),
    ('primase','Primase.'),('recombinase','Recombinase, integrase, resolvase or broad rec[a-z] pattern.'),
    ('cluster8_GT_like','Exact query_id == manifest cluster8_query_id, plus GT-specific or TET-associated GT text.'),
    ('TET_associated_GT_like','TET-associated glycosyltransferase phrases or broad tet.*gt regex; independent of exact cluster8 query match.'),
    ('additional_GT','GT-specific text on a non-exact cluster8 query, excluding generic transferase/kinase/ligase/polymerase/NTase terms.'),
    ('NDP_sugar_enzyme','NDP/nucleotide/UDP/GDP/CDP/TDP sugar; UDP-glucose/GlcNAc/galactose, GDP-mannose, dTDP or rml/rfb A–D.'),
    ('epimerase','Epimerase/dehydratase/dehydrogenase/isomerase/mutase AND sugar or nucleotide-name pattern; thus broader than epimerases.'),
    ('dehydratase','Any dehydratase.'),('reductase','Any reductase.'),('aminotransferase','Aminotransferase or transaminase.'),
    ('methyltransferase','Any methyltransferase; overlaps C5_MT.'),('acetyltransferase','Acetyltransferase or acyltransferase.'),
    ('sulfotransferase','Sulfotransferase.'),('APS_PAPS_kinase','APS kinase/adenylyl-sulfate kinase/PAPS/sulfate adenylyltransferase/atps/adenosine phosphosulfate; broader than one kinase.'),
    ('capsid','Capsid or HK97.'),('portal','Portal.'),('tail','Tail, baseplate or fiber.'),('tape_measure','Tape measure/tapemeasure.'),
    ('terminase','Terminase.'),('conserved_hypothetical','Conserved hypothetical, hypothetical protein or uncharacterized protein.'),
    ('DUF','DUF followed by digits or domain of unknown function.')]
    assert len(rules)==33
    write('docs/GENE_CLASS_ASSIGNMENTS.md','# Gene-class assignments\n\nHuman-readable audit of all 33 classes. The exact regular expressions and GT exclusion/overlap conditions remain in the recovered `classify_hit` function; this table is explanatory, not replacement code. All rules search the concatenated fields described in Methods.\n\n| Class | Recovered matching logic |\n|---|---|\n'+'\n'.join('| '+c+' | '+r+' |' for c,r in rules))
    write('docs/KNOWN_LIMITATIONS.md','''
# Historical limitations and separate interpretation correction

## TET-positive loci hidden by historical grouping

A0A6J5M498 has TET_5mYOX=3 and C5_MT=0, but is labeled additional-GT/glycan-associated because the first row-class branch requires both C5-MT and TET, and there is no standalone TET branch. Its `inferred_model` already says hybrid DNA-modification + glycan module. This is a mismatch in the grouping scheme, not evidence that the TET-like score should be removed.

The saved evidence specifically identifies UFOVP405_24, afdb-swissprot rank 14, AF-D2W6T1-F1-model_v6, “Tet-like dioxygenase 1”, probability 1.0, E-value 5.287e−9, query coverage 0.9183, alignment length 257. The per-locus table also contains stronger-ranked TET/JBP structural matches with other naming forms, illustrating incomplete name matching. Structural association does not establish a particular DNA substrate or biochemical activity.

A0A3G1L3B7 likewise has TET_5mYOX=2, C5_MT=0; it is historically grouped under broad nucleotide metabolism. The separate correction applies the same rule consistently to both loci. It adds `TET/5mYOX-associated (C5-MT unsupported)` for TET≥2, C5-MT<2 and dCMP-HMase<2. It changes primary_group and display group_rank only (to 2.5); all scores, query IDs, evidence, inferred models and source paths remain unchanged. The corrected plot uses the recovered plotting code with this added group color/order/legend entry. This is an explicitly proposed interpretation, not a new biological validation or replacement of the historical analysis.

## Broad nucleotide and glycan labels

The historical “ThyX/ThyS / nucleotide-metabolism-associated” group accepts a dCMP/dTMP/pyrimidine-related match or a nucleotidyltransferase/pyrophosphorylase match, even without ThyX or ThyS. All 11 ThyX values are zero; all five rows assigned this broad group also have ThyS=0. The shortened figure legend “ThyX/ThyS” is therefore misleading if interpreted literally. Nucleotide/sugar metabolism is not evidence by itself for DNA glycosylation or hydroxymethyl-base production.

Generic reductase, ligase, acyltransferase, epimerase-like and APS/PAPS annotations are broad. A single hit may populate several classes and affect several totals; database redundancy does not constitute independent replication. TET_associated_GT_like uses a broad `tet.*gt` expression as well as explicit phrases, so this column requires hit-level review before mechanistic interpretation.

## Ambiguous ORFs and evidence limitations

Maximum-per-class aggregation can assign mutually alternative functions to the same ORF from different hits, or count broad domains as whole-protein function. Target coverage is unused in the score, and functional specificity, catalytic residues, substrate identity, structural superposition and genomic orientation/distance are not independently tested by the matrix builder. `validation/ambiguous_multiclass_ORFs.json` inventories ORFs with two or more score≥2 classes; multi-class status is a review flag, not proof of an erroneous assignment. Many valid multi-domain/overlapping annotations exist.

Unknown/hypothetical and function-assigned hits can coexist. Zero is a lack of a qualifying named match, not confirmed absence. Neighborhoods contain differing query counts, and this is a selected 11-locus panel rather than an unbiased prevalence estimate. The heatmap summarizes context, not gene order or a demonstrated pathway.

## Accession and query-ID mismatches

| Candidate | Fixed-manifest focal accession | Query namespace in saved table |
|---|---|---|
| B2ZYD6 | BAG41682.1 | BAG_### |
| A0A5P8D8M7 | QFP95377.1 | SEA_YUUY_## |
| A0A3G1L3B7 | ATW62693.1 | SCBWM1_gp# |
| A0A8F2DCD3 | QWS69713.1 | SEA_WILLIAMSTRONG_## |
| A0A159B7X0 | AKJ73160.1 | AKJ_### |

The builder requires exact equality to recognize the focal Cluster 8 ORF. These five focal IDs do not occur in their tables, and all five have historical cluster8_GT_like=0 despite GT-like neighborhood hits. Their focal GT may instead contribute to additional_GT. Related manifests recover local model/query aliases but do not, by themselves, establish all accession-to-query mappings by sequence. No accession was silently remapped. BAG41682/BAG41686 or any other cross-session identity must not be equated without sequence/source evidence.

## Historical provenance limits

The completed raw combined-hit TSVs, annotation outputs and code are recovered; the exact command line and generating script version for every historical Foldseek submission are not conclusively established. Multiple batch script variants exist; all recovered annotation script copies share one hash, while the June 17 batch differs from later copies. Recovery preserves variants and reports the distinction rather than choosing one as the proven generator. Absolute paths in preserved manifests are historical, not portable paths. API database versions and server state cannot be recreated by rerunning an online search today; no new searches were made. Related manifests sometimes name other loci or earlier templates, so directory/name proximity is not treated as confirmed lineage.

## Rendering and publication

Exact PNG reproduction was verified in the existing Linux environment; Windows Matplotlib 3.9.1 produced a one-pixel width difference. The PDF differs only in CreationDate in the matching Linux environment. The recovered layout has some header/row-label crowding and was preserved. This package validates computational reproducibility with caveats, not biochemical function or publication rights. Creators, licensing and final manuscript figure associations still require author confirmation in release materials.
''')
    env=v['environment']
    write('docs/REPRODUCIBILITY.md',f'''
# Reproduction and environment

The successful observed environment was existing Ubuntu/WSL with Python 3.12.11, pandas {env['pandas']}, NumPy {env['numpy']}, Matplotlib {env['matplotlib']}, Pillow {env['pillow']} and Matplotlib FreeType {env['freetype']}. The original PNG records Matplotlib 3.10.8. This is a tested matching rendering environment, not a claim that the entire June environment was recovered. Font/build differences can affect raster bytes even with equal Python package versions.

`environment.yml` records the tested key versions. Additional requests is only for inspecting/running the archived API script; it is deliberately unnecessary for offline reproduction. Do not invoke the API runner for this workflow.

## Offline portable workflow

Run `python -B reproduce.py --out /new/output/path --require-exact-render` from this workflow directory in a matching Linux environment. A data-only portability check can omit `--require-exact-render`; inspect `validation.json` to distinguish exact image comparison from data equality. The code reads each preserved top-hit table by candidate-relative path, calls the recovered functions without edits, retains historical analysis_dir strings for comparison and runs the recovered plot CLI. The original evidence table is preserved byte-for-byte; regenerated evidence is compared as a multiset because set iteration changes row ordering.

The CLI below is also recovered and can be used with a new manifest whose **analysis_dir only** is redirected to the preserved per-locus `final/analysis_v2` folders. Retain candidate IDs and focal query IDs exactly, including known mismatches:

```bash
python -B originals/260630_Temp/build_cluster8_neighborhood_matrix.py --manifest /new/portable_manifest.tsv --outdir /new/matrix_output
python -B originals/260630_Temp/plot_cluster8_panelB_heatmap_clean.py --matrix originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_matrix.tsv --outdir /new/figure_output
```

The new manifest's paths will appear in the generated matrix and differ from preserved historical paths; scientific fields should agree. The supplied `reproduce.py` avoids this metadata difference without editing inputs. No reconstruction of builder/plotter code occurred; the new script is only an offline validation wrapper.

SHA-256 manifests refer to paths relative to this workflow. `SHA256SUMS-originals.txt` covers original bytes; `SHA256SUMS-package.txt` covers all package files except itself. The outer ZIP has its own sidecar checksum. `.gitattributes` disables text conversion within this workflow so Git does not change archived CRLF/LF bytes. The complete recovery includes 346 files and about 315 MB of original bytes, including intermediate/duplicate historical tables; this broad preservation is intentional. No single source file exceeds 7 MB.

Executed validation is in `validation/validation.json`; source preservation in `audit/source_preservation_final.json`; the independent saved-evidence maxima check covers all 363 cells. Original PDF metadata is never modified: CreationDate normalization is comparison-only.
''')
    write('environment.yml',f'''
name: wp08-cluster8-recovery
channels:
  - conda-forge
dependencies:
  - python=3.12.11
  - pandas={env['pandas']}
  - numpy={env['numpy']}
  - matplotlib={env['matplotlib']}
  - pillow={env['pillow']}
''')
    write('docs/PROVENANCE.md','''
# Provenance and preservation

Source authority is the recovered local completed analysis, not a retrospective reconstruction or online annotation update. Original output metadata and source mtimes place completion on June 30, 2026; recovery/validation was performed October 29, 2026. Dates are reported as recorded, without asserting a verified workstation clock history.

`audit/provenance.json` records every absolute Ubuntu source path, resolved path, archive path, size, SHA-256 and source mtime in UTC/nanoseconds. Python `copy2` preserved file bytes and mtimes. Source atimes can change from reading; ownership, source inode/ctime, permissions and Windows filesystem attributes are not asserted as forensic copies. One Zone.Identifier sidecar was excluded; no scientific file in the selected directory was excluded. Other local folders were not exhaustively copied.

The original analysis tree is preserved under `originals/260630_Temp/`; related scripts/manifests are preserved under their original directory names in `originals/related_scripts/`. `audit/script_variants.json` identifies 32 script copies, including different batch versions; `audit/supplemental_manifests.json` identifies 35 related manifests. The 15 recovered annotation-script copies all share SHA-256 `944c08dc08e4ed32722a9fe7bd23479dc663a5ce98c5216e21452a1d69d84d4f`. There are 17 batch-script copies (one older hash and 16 matching later copies); consult the inventory rather than treating duplicate copies as independent evidence.

`audit/recover_wp08.py` and `audit/build_wp08_package.py` are newly authored recovery/packaging utilities, not original scientific scripts. `reproduce.py` and `correct_interpretation.py` are new wrappers with their roles explicitly documented. The historical classifier/plotter is imported from original bytes. No preserved original code was edited.

The source repository inspected was `/home/rahnn/projects/viral-gt-research-record`, remote `https://github.com/LanNa1229/viral-gt-research-record.git`. The existing checkout was not edited; no commit/push or release API was invoked. Its unrelated working changes are recorded, not included or modified. `audit/repository_context/` contains the read-only context used for preparation and is not a claim that the public default branch was refreshed. The private analysis repository was only inspected for identity/context; no uploads were made there.

The local `Zenodo_publication_receipt_23724706.json` in the Windows workspace identifies a separate published structural-analysis record with DOI `10.5281/zenodo.23724706` and creator display name Lan Na. That local receipt supplies draft related-record context only; this task did not change or re-publish that record, and no new DOI was reserved.
''')
    # Verify exact variant counts from recovered inventory instead of assuming a generator.
    variants=json.loads((PKG/'audit/script_variants.json').read_text())
    dump('audit/script_variant_summary.json',{name:dict(copies=len([r for r in variants if Path(r['path']).name==name]),hashes=sorted(set(r['sha256'] for r in variants if Path(r['path']).name==name))) for name in ['foldseek_batch_api_v2.py','analyze_foldseek_results_v2.py']})
    write('release/GITHUB_PR_DRAFT.md','''
# Recover original WP08 Cluster 8 neighborhood analysis and validate historical heatmap

The completed 11-locus Cluster 8 genomic-context heatmap was absent from the packaged research record. This change preserves the recovered original scripts, saved Foldseek inputs/annotations, manifest, matrix, evidence and figures under `workflows/WP08_cluster8_gene_neighborhood_heatmap/`, with SHA-256 provenance and an offline reproduction wrapper.

Validation reproduces all 363 score cells and 4,919 evidence records. The PNG matches byte-for-byte using the existing matching Ubuntu environment; the PDF matches except for CreationDate. All 346 source files retain their contents and source mtimes. No Foldseek jobs were rerun. The documented caveats include TET-positive row-classification fall-through, broad nucleotide-metabolism labels, ambiguous ORFs and five exact accession/query mismatches. A separate interpretation-only output changes group metadata for two TET-positive loci; all historical data and figures are unchanged.

Scope is limited to the new WP08 directory; `.gitattributes` within it preserves original bytes. No root README/ledger or unrelated workflow edits are bundled. This text is a draft; no PR or push has been made. Before publication, verify the current target branch, review rights/creators and obtain scoped push approval under repository AGENTS.md.
''')
    write('release/RELEASE_NOTES.md','''
# WP08 Cluster 8 historical neighborhood analysis — proposed v1.0.0

Recovered completed June 30, 2026 analysis; validated October 29, 2026. Proposed tag: `wp08-cluster8-neighborhood-v1.0.0` (not created).

Includes the original 11-locus ×33-class matrix, 4,919-row evidence table, all saved per-locus annotations and combined hit tables, original PDF/PNG figures, all recovered relevant script variants and related query manifests. SHA-256 inventories, source provenance, Methods, gene-class rules and offline validation are included.

The clean PNG reproduces exactly in the tested Ubuntu environment. The reproduced PDF differs only in creation timestamp. Separate interpretation-only outputs add a standalone TET-positive group for A0A6J5M498 and A0A3G1L3B7; they do not replace historical outputs or alter scores.

Known limits: keyword-derived functional classes, maximum-hit aggregation, broad nucleotide/glycan bins, exact focal-query mismatches and incomplete historical submission lineage. No biochemical claims were newly validated and no online Foldseek searches were rerun. Cite the future WP08 DOI only after it is minted; no WP08 DOI exists in this candidate.
''')
    dump('release/ZENODO_METADATA_PROPOSAL.json',dict(
      status='LOCAL DRAFT ONLY; not an upload-ready API payload; no deposition or DOI created',
      title='Cluster 8 viral glycosyltransferase gene-neighborhood heatmap: recovered historical analysis and reproducibility record',
      version='1.0.0 (proposed)',resource_type='dataset',
      creators=[dict(name='Lan Na',status='Draft display name from prior local structural-record receipt; author list/order/ORCIDs require confirmation')],
      license=None,publication_date=None,
      description='Recovered original June 30, 2026 completed 11-locus genomic-context analysis with unmodified historical code, inputs, outputs and figures; SHA-256 provenance; validated exact PNG reproduction; explicit limitations; separate interpretation-only row grouping correction. No Foldseek searches rerun. Computational reproducibility does not establish biochemical activity.',
      keywords=['viral glycosyltransferase','Cluster 8','genomic neighborhood','Foldseek','DNA modification','reproducibility'],
      related_identifiers=[dict(identifier='https://github.com/LanNa1229/viral-gt-research-record',relation='isSupplementTo',scheme='url'),dict(identifier='10.5281/zenodo.23724706',relation='references',scheme='doi',note='Separate structural record identified by local prior publication receipt; verify relation before publication')],
      required_before_publication=['Confirm complete creators, order, ORCIDs and affiliations','Confirm license and redistribution rights for original scripts, database-derived annotations and all other material','Confirm final title/version and manuscript figure/citation relationship','Decide whether the separate interpretation-only correction should accompany historical-only release','Approve scoped GitHub push and separately approve Zenodo upload/publication'],
      deposit_file='WP08_cluster8_gene_neighborhood_heatmap_v1.0.0_candidate.zip'))
    write('release/PUBLICATION_CHECKLIST.md','''
# Publication handoff (no publication performed)

The scientific reproducibility checks pass. Release metadata and rights are still draft. The candidate tree maps directly to the requested GitHub workflow path, and the candidate ZIP can serve as the Zenodo payload after review.

1. Review recovery report, original checksums, caveats and the separate interpretation-only changes. Confirm authors, rights/license, title/version and manuscript associations; do not use an unreviewed generic license.
2. Re-check candidate checksums and current target branch. Preserve existing unrelated changes in the repository; do not overwrite its checkout or stage unrelated files. Review any deeper repository instructions before applying the tree.
3. After explicit scoped approval, apply only `workflows/WP08_cluster8_gene_neighborhood_heatmap/` on an isolated branch, confirm Git blobs match archived bytes (including CRLF), review the diff and run offline reproduction. Do not change root ledgers/README without including them in an approved scope.
4. Obtain explicit scoped push approval as required by repository `AGENTS.md`, then create the reviewed PR/release if requested. No branch/tag/PR was pre-created here.
5. Review Zenodo metadata proposal and replace unresolved fields. On separate authorization create a new WP08 deposit; do not overwrite the existing structural-analysis record. Upload the checksum-verified candidate archive; attach version-specific GitHub link only after a real tag/commit exists. Reserve/record a DOI only through an authorized API operation.
6. Verify the uploaded archive's content/checksum and metadata; publish only after authorized. Add actual DOI/release URLs to citation and release records in a later scoped update. Until then do not claim the local candidate is deposited or published.

Repository publication policy is preserved in `audit/repository_context/viral-gt-research-record/AGENTS.md`. The candidate was intentionally prepared outside the existing checkout, as requested. No confirmation is needed to inspect or reproduce these local files; publication remains unrequested.
''')
    write('release/RIGHTS_AND_CITATION_DRAFT.md','''
# Rights and citation draft

No new license is assigned by this candidate. Preserve original headers and database/source attribution. No root LICENSE file was found at the exact inspected path; broader rights have not been inferred from GitHub visibility. Creator display name “Lan Na” comes from a prior local structural-record receipt; the complete author list, order and ORCIDs remain for author confirmation.

Proposed citation text, pending author and publication review:

Lan Na [additional creators/order to confirm]. Cluster 8 viral glycosyltransferase gene-neighborhood heatmap: recovered historical analysis and reproducibility record. Proposed version 1.0.0. DOI: pending; no WP08 record has been created.

This record is distinct from the prior structural-analysis DOI `10.5281/zenodo.23724706`. Do not reuse that DOI as though it identifies this recovered heatmap package.
''')
    files=sorted(p for p in PKG.rglob('*') if p.is_file())
    dump('audit/package_validation.json',dict(status='Validated historical reproducibility; metadata/rights review pending before publication',original_files=346,all_original_hashes_match=True,source_bytes_and_mtimes_verified=True,matrix_equal=True,evidence_multiset_equal=True,independent_maxima_equal=True,png_byte_and_pixel_equal=True,pdf_equal_except_creation_date=True,corrected_loci=['A0A3G1L3B7','A0A6J5M498'],existing_repository_modified=False,published=False))
    files=sorted(p for p in PKG.rglob('*') if p.is_file())
    write('SHA256SUMS-package.txt',''.join(f'{sha(p)}  {p.relative_to(PKG).as_posix()}\n' for p in files))
    archive=RELEASE/'WP08_cluster8_gene_neighborhood_heatmap_v1.0.0_candidate.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(p for p in PKG.rglob('*') if p.is_file()): z.write(p,p.relative_to(RELEASE).as_posix())
    archive.with_suffix('.zip.sha256').write_text(f'{sha(archive)}  {archive.name}\n')
    # Read every member back and compare bytes, not only ZIP CRC.
    with zipfile.ZipFile(archive) as z:
        for p in sorted(p for p in PKG.rglob('*') if p.is_file()):
            assert hashlib.sha256(z.read(p.relative_to(RELEASE).as_posix())).hexdigest()==sha(p)
    print(json.dumps(dict(package=str(PKG),archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,files=len(files)+1),indent=2))

if __name__=='__main__': main()
