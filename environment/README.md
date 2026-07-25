# Software Environment Records

## Purpose

This directory preserves the computational environment needed to reproduce analyses. Environment records should be generated for each major work package and milestone.

## Recommended files

```text
environment/
├── environment_WP05_YYYYMMDD.yml
├── environment_WP05_YYYYMMDD_explicit.txt
├── pip_freeze_WP05_YYYYMMDD.txt
├── system_WP05_YYYYMMDD.txt
├── tool_versions_WP05_YYYYMMDD.tsv
├── container_WP05.def
└── README.md
```

## Capture commands

### Conda

```bash
conda env export --from-history > environment/environment_WP05_$(date +%Y%m%d).yml
conda list --explicit > environment/environment_WP05_$(date +%Y%m%d)_explicit.txt
```

`--from-history` gives a readable specification of intentionally installed packages. The explicit export is more exact but less portable. Preserve both for major milestones.

### Pip

```bash
python -m pip freeze > environment/pip_freeze_WP05_$(date +%Y%m%d).txt
```

### System and tools

```bash
{
  date -u +"captured_utc=%Y-%m-%dT%H:%M:%SZ"
  uname -a
  python --version
  git --version
  quarto --version 2>/dev/null || true
  foldseek version 2>/dev/null || true
} > environment/system_WP05_$(date +%Y%m%d).txt
```

### Git linkage

```bash
git rev-parse HEAD > environment/git_commit_WP05_$(date +%Y%m%d).txt
git status --short > environment/git_status_WP05_$(date +%Y%m%d).txt
```

## Tool-version table

Use tab-separated columns:

```text
tool	version	build_or_container	capture_date_utc	notebook_id	notes
```

## Rules

1. Never commit credentials or license files.
2. Do not commit the full Conda environment directory.
3. Record database versions separately from software versions.
4. For containers, preserve the recipe or digest, not only a mutable tag.
5. Update the notebook with the exact environment filename and Git commit.
6. If an analysis is rerun under a different environment, create a new environment record rather than overwriting the previous file.
