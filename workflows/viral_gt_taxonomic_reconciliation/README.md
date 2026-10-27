# Viral GT Step 6: taxonomic reconciliation

Step 6 reconciles accession-specific NCBI/ICTV taxonomy with full-genome vConTACT3 and geNomad evidence. It retains source labels, uncertain assignments and the original 466 screened units, 245 candidates (190 Both, 29 11-only, 26 15-only) and separate 242-genome Sankey subset.

## Inputs and setup

Run the unchanged historical scripts in a working copy of the original project layout. Their script parent is the analysis directory; the project root is two levels above it. Restore these inputs before execution:

- `analysis/viral_gt_figureB_taxonomy_20261018/tables/master_genome_dataset.tsv`: original 466-row master, from [Step 4](https://doi.org/10.5281/zenodo.23685561), archive member `outputs/master_genome_dataset.tsv`.
- `analysis/viral_gt_mining_20261017/runs/baseline_594_v2/viral_gt_windows.tsv` and cached source GenBank records: [Step 1](https://doi.org/10.5281/zenodo.23680510). Source paths in the window table must resolve to the restored GenBank files; remap paths only in a working copy if necessary.
- Frozen MSL41.v1 and VMR_MSL41.v1.20260729 workbooks under `sources/ictv/`: [corrected MSL41](https://doi.org/10.5281/zenodo.23376199) and [corrected VMR](https://doi.org/10.5281/zenodo.21694279). `fetch_metadata.py` reuses cached files/responses; without them its `/current` endpoints and live NCBI queries may return newer metadata.
- Existing vConTACT3 v230 database/`230.json` under `db/vcontact3/` ([RefSeq230/GenBank266](https://doi.org/10.5281/zenodo.15598380)); geNomad db1.9 under `db/genomad_db/` ([database](https://doi.org/10.5281/zenodo.14886553)). Reuse these databases; neither is bundled here.

In a separate restored project, place this directory's `scripts/*.py` in `analysis/viral_taxonomy_reconciliation_20261019/scripts/`. Reuse the existing environment where available. For a new installation, this is a version-pinned reconstruction of the historical mamba installation, not a newly validated environment solve:

```sh
PROJECT=/absolute/path/to/restored-project
RUN="$PROJECT/analysis/viral_taxonomy_reconciliation_20261019"
mamba create -y -p "$RUN/env/taxonomy" --override-channels   -c conda-forge -c bioconda python=3.11.17 vcontact3=3.1.6 genomad=1.13.0   numpy=1.26.4 pandas=3.0.6 pyarrow=24.0.0 mmseqs2=18.8cc5c
PY="$RUN/env/taxonomy/bin/python"
"$PY" -m pip install --no-deps --target "$RUN/compat/pandas2"   pandas==2.3.3 pytz==2026.5
```

The historical vConTACT3 launcher uses this pandas/pytz overlay and the local missing-realm-threshold guard. Source code is preserved; no upstream software or models are redistributed.

## Execution

Execute sequentially in the restored working copy. Skip completed stages when reusing their existing outputs. These commands document full analysis; preparing the present release did not run them.

```sh
"$PY" -B "$RUN/scripts/fetch_metadata.py"
"$PY" -B "$RUN/scripts/reconcile_metadata.py"
"$PY" -B "$RUN/scripts/prepare_genomes.py"
# Only when starting from the official cached v230.tar.gz, not an unpacked DB:
# "$PY" -B "$RUN/scripts/prepare_vcontact_reference.py"
"$PY" -B "$RUN/scripts/install_runtime_guard.py"
"$PY" -B "$RUN/scripts/run_classifier.py" genomad
"$PY" -B "$RUN/scripts/summarize_markers.py"
"$PY" -B "$RUN/scripts/run_classifier.py" prokaryotes full
"$PY" -B "$RUN/scripts/run_classifier.py" eukaryotes full
"$PY" -B "$RUN/scripts/run_classifier.py" prokaryotes no_gt
"$PY" -B "$RUN/scripts/run_classifier.py" eukaryotes no_gt
"$PY" -B "$RUN/scripts/prepare_domain_crosschecks.py"
# Historical crosschecks: seven genomes; run only if these inputs were emitted.
"$PY" -B "$RUN/scripts/run_classifier.py" prokaryotes full crosscheck
"$PY" -B "$RUN/scripts/run_classifier.py" prokaryotes no_gt crosscheck
"$PY" -B "$RUN/scripts/audit_reference_identity.py"
"$PY" -B "$RUN/scripts/combine_evidence.py"
```

`prepare_vcontact_reference.py` verifies official archive MD5 `dbf915c9da182e0c3dbec60d2156f8b9` and excludes publisher scratch files during extraction. The genome-preparation script uses full available source records, cached CDS translations, genome lengths and selected GT-anchor coordinates. vConTACT3 runs protein/gene-to-genome input, four threads, `--reduce-memory`, v230, full and selected-anchor-depleted modes. geNomad runs **annotate only**, two threads, 32 splits and `--full-ictv-lineage`; no viral-identification score or virus-score threshold was used.

## Evidence reconciliation and outputs

NCBI ESummary/Taxonomy EFetch supplies versioned accession metadata and virus lineages. VMR accession-token/base matches take priority; specific accession-linked species/genus names are mapped to MSL41, with documented renames. Generic inherited taxonomy is not independent family/genus authority. Alias units follow frozen exact-nucleotide sequence hashes, not shared organism names.

`audit_reference_identity.py` identifies exact reference copies; `combine_evidence.py` excludes self/alias/exact-copy and generic-reference support, reconciles reference names to MSL41, and retains conflicts. Accepted predicted families require full/no-selected-GT agreement, >=2 specific non-self references, purity >=0.8, >=5 family markers, >=3 profiles, bitscore fraction >=0.8, marker span >=0.25 and ambiguous bases <0.05. Fragments require >=10 markers and span >=0.5. Provisional higher ranks require >=3 markers and >=0.5 agreement. “Non-GT” means removal of selected cohort anchors, not all GTs. geNomad db1.9 uses older MSL39 taxonomy; ambiguous family splits are not guessed. No predicted genus/species is promoted to official taxonomy.

Expected outputs include accession evidence/phase1 mapping, full-genome quality and classifier inputs, raw assignment/annotation outputs, marker/reference evidence, and `tables/final_genome_taxonomy_master.tsv` (466 rows) with `final_classification_decisions.tsv`. The separately preserved 245-row, 79-column integrated Master TSV is the final dataset; the historical program assembling that exact export was not recovered, so these scripts do not claim to recreate its 79-column layout. Full classifier execution requires the listed inputs/databases; the one-file Zenodo dataset alone is not an end-to-end rerun bundle.

Reuse [Step 2 predictions/embeddings](https://doi.org/10.5281/zenodo.23682612), [Step 3 retrieval](https://doi.org/10.5281/zenodo.23683720) and the already public [Step 4 selection code](https://github.com/LanNa1229/viral-gt-research-record/blob/a72e62ad7a5efc23b474564a28e8c3c85ed016f7/workflows/viral_gt_candidate_selection/historical/build_figure_b.py); they are not repackaged here.

Original independent scripts are MIT-licensed, using the [existing MIT text](https://github.com/LanNa1229/viral-gt-research-record/blob/7706e5d2e3b1f360ac7dae75454a355931a7d45c/workflows/viral_gt_neighborhood_extraction/LICENSES/MIT.txt). The vConTACT3 guard installer/generated modification is conservatively scoped GPL-3.0; its license is included in `LICENSES/GPL-3.0.txt`. geNomad's installed academic/internal/noncommercial software license remains applicable; its separately CC BY reference database does not relicense the software. No upstream package, database or model is included.
