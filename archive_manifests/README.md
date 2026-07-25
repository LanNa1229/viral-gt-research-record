# Archive Manifests and Milestone Snapshots

## Purpose

This directory indexes immutable or read-only milestone snapshots of the research record. Snapshots preserve the repository state and the relationship between Git-tracked files and externally stored research data.

A snapshot is not a substitute for the active repository. It is a point-in-time record used for reproducibility, continuity, institutional access, and manuscript milestones.

## Snapshot contents

Each snapshot directory should contain:

```text
snapshot-YYYYMMDDTHHMMSSZ-<shortcommit>/
├── SNAPSHOT_METADATA.tsv
├── git_status.txt
├── git_log.txt
├── tracked_files.sha256
├── external_datasets.tsv
├── repository.bundle
├── repository_<tag>.tar.gz
├── repository_<tag>.zip
├── SHA256SUMS
└── verification.txt
```

## Snapshot metadata fields

```text
snapshot_id
created_utc
created_by
repository_name
commit
annotated_tag
working_tree_status
institutional_storage_location
copy_verified_date
verified_by
related_manuscript_version
notes
```

## External dataset manifest

The repository archive does not duplicate large raw or intermediate data. `external_datasets.tsv` must list:

```text
dataset_id
data_tier
storage_location
file_manifest
manifest_sha256
size_bytes
access_classification
responsible_steward
related_notebook
related_claim
```

## Creating a snapshot

Use:

```bash
bash scripts/create_snapshot.sh milestone-v0.1.0-20260723 /path/to/uga/backup/location
```

The script creates a Git bundle, TAR archive, ZIP archive, tracked-file checksum manifest, and archive checksum file.

## Verification

After copying to institutional storage:

```bash
sha256sum -c SHA256SUMS
git bundle verify repository.bundle
```

Record the result in `verification.txt` and do not overwrite a prior verification record.

## Institutional backup

The backup location should be accessible through UGA-approved credentials and selected according to data classification. A personal device or personal cloud account must not be the only archived copy.

## Retention

Do not delete milestone snapshots merely because a later snapshot exists. Retention must follow UGA, sponsor, contract, journal, and project requirements. When a hold, audit, dispute, investigation, or other preservation obligation applies, suspend routine deletion and seek appropriate institutional guidance.
