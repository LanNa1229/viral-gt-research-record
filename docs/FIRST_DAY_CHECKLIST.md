# First-Day Checklist

Date: 2026-07-23

- [ ] Create a private remote repository in an authorized UGA-controlled organization or approved private service.
- [ ] Clone or initialize the repository on the Linux/HPC system.
- [ ] Copy this starter scaffold into the repository.
- [ ] Review `README.md` and replace placeholders.
- [ ] Confirm the repository remains private.
- [ ] Configure individual Git identity; do not use a shared account.
- [ ] Commit the governance files.
- [ ] Create `NB-VGT-WP00-20260723-001_repository-initialization.qmd`.
- [ ] Enter the current manuscript scope in `CLAIM_FIGURE_EVIDENCE_MATRIX.tsv`.
- [ ] Review and edit the draft CRediT rows; mark uncorroborated historical rows as retrospective.
- [ ] Create read-only copies of the highest-value pre-existing scripts, manuscript drafts, figures, and outputs.
- [ ] Generate SHA-256 manifests for the archived historical files.
- [ ] Create retrospective records for the structure-first framework and Cluster 8.
- [ ] Create annotated tag `milestone-v0.1.0-20260723`.
- [ ] Run `scripts/create_snapshot.sh`.
- [ ] Copy the snapshot to UGA-managed or UGA-approved storage.
- [ ] Verify the copied snapshot with `sha256sum -c SHA256SUMS`.
- [ ] Record the institutional storage location in the snapshot metadata.
- [ ] Confirm that raw and large intermediate datasets are not accidentally staged in Git.
- [ ] Record who currently has repository access.
