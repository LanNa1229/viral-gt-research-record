# Project Timeline

## Purpose

This timeline records major scientific and administrative milestones for the Viral GT project. It is a high-level index and must link to notebooks, Git commits, archived records, or retrospective provenance entries.

Historical events entered after repository initialization must be labeled `Retrospective` and must not be presented as if they were recorded contemporaneously.

## Status vocabulary

- `Planned`
- `In progress`
- `Completed`
- `Superseded`
- `On hold`
- `Retrospective reconstruction pending`
- `Retrospective entry completed`

## Timeline

| Record date | Work period | Entry ID | Type | Work package | Milestone or decision | Scientific significance | Evidence | Git commit or tag | Status |
|---|---|---|---|---|---|---|---|---|---|
| 2026-07-23 | 2026-07-23 | ENT-VGT-20260723-001 | Contemporaneous | Project-wide | Initialized private Git repository and formal computational notebook system | Establishes prospective provenance, reproducibility, contribution tracking, and institutional backup workflow | README.md; CONTRIBUTIONS.md; initial ledgers | `<initial-commit>` | Completed |
| 2026-07-23 | Prior work; date range TBD | RET-VGT-20260723-001 | Retrospective | WP01–WP06 | Reconstruct origin and development of the structure-first viral GT discovery framework | Links early search design, structural clustering, SSN construction, and weighted Leiden analysis to contemporaneous records | Earliest scripts, HPC logs, presentations, emails, manuscript versions | `<commit-after-entry>` | Retrospective reconstruction pending |
| 2026-07-23 | Prior work; date range TBD | RET-VGT-20260723-002 | Retrospective | WP08–WP09 | Reconstruct Cluster 8 discovery, neighborhood interpretation, jumbo-phage sulfated-DNA hypothesis, and validation design | Preserves the chronology connecting structural placement to pathway-level biological interpretation | Network outputs, neighborhood files, draft text, meeting summaries, email attachments | `<commit-after-entry>` | Retrospective reconstruction pending |
| 2026-07-23 | Prior work; date range TBD | RET-VGT-20260723-003 | Retrospective | WP12 | Reconstruct manuscript figure architecture and claim development | Connects computational outputs to specific manuscript claims and figure panels | Figure plans, presentation slides, manuscript drafts | `<commit-after-entry>` | Retrospective reconstruction pending |

## New-entry template

Copy the following row and complete all fields:

```text
| YYYY-MM-DD | YYYY-MM-DD to YYYY-MM-DD | ENT-or-RET-ID | Contemporaneous or Retrospective | WP## | Milestone or decision | Scientific significance | Evidence location | Commit or tag | Status |
```

## Milestone tag index

| Tag | Date | Scope | Commit | Snapshot manifest | Notes |
|---|---|---|---|---|---|
| `milestone-v0.1.0-20260723` | 2026-07-23 | Repository and research-record system initialized | `<hash>` | `archive_manifests/<snapshot>/SHA256SUMS` | Initial governance milestone |

## Rules

1. Record the date on which the entry is created.
2. For historical work, separately record the reconstructed work period.
3. Do not use file-modification timestamps as the sole evidence of scientific conception.
4. Link milestone statements to the underlying notebook or retrospective entry.
5. When a milestone is revised, add a new row and mark the earlier row `Superseded`; do not silently rewrite history.
