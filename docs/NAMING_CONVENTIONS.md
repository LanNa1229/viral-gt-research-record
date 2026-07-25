# Naming Conventions

## Identifier formats

| Object | Format | Example |
|---|---|---|
| Research entry | `ENT-VGT-YYYYMMDD-NNN` | `ENT-VGT-20260723-001` |
| Notebook | `NB-VGT-WP##-YYYYMMDD-NNN` | `NB-VGT-WP05-20260723-001` |
| Dataset | `DS-VGT-YYYYMMDD-NNN` | `DS-VGT-20260723-001` |
| Output | `OUT-<EntryID>-<description>` | `OUT-ENT-VGT-20260723-001-community-audit.tsv` |
| Claim | `CLM-VGT-NNN` | `CLM-VGT-003` |
| Retrospective record | `RET-VGT-YYYYMMDD-NNN` | `RET-VGT-20260723-002` |
| Transfer | `XFER-VGT-YYYYMMDD-NNN` | `XFER-VGT-20260723-001` |
| Evidence | `EVD-VGT-YYYYMMDD-NNN` | `EVD-VGT-20260723-001` |

## Notebook filenames

```text
NB-VGT-WP05-20260723-001_weighted-leiden-community-analysis.qmd
NB-VGT-WP09-20260723-001_jumbo-phage-sulfation-neighborhoods.qmd
RET-VGT-20260723-002_cluster8-provenance-reconstruction.qmd
```

## Figure filenames

```text
FIG03-PANELA_weighted-leiden-network_20260723.svg
FIG03-PANELB_viral-placement-classes_20260723.pdf
FIG03_manifest_20260723.tsv
```

Do not use `final`, `final2`, `latest`, or `new`.

## Output filenames

```text
OUT-ENT-VGT-20260723-004_cluster8-neighborhood-classification.tsv
OUT-ENT-VGT-20260723-004_cluster8-neighborhood-heatmap.svg
OUT-ENT-VGT-20260723-004_run-log.txt
```

## Tags

```text
milestone-v0.1.0-20260723
manuscript-v0.3.0-20260815
submission-v1.0.0-YYYYMMDD
revision1-v1.1.0-YYYYMMDD
```

## Branches

```text
analysis/WP05-weighted-leiden
analysis/WP09-jumbo-phage-neighborhoods
figure/FIG03-viral-placement
docs/ENT-VGT-20260723-001-governance
retrospective/RET-VGT-20260723-002-cluster8
fix/WP03-coverage-calculation
```
