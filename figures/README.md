# Figure Records

- `source/`: small, editable source files and figure-generation metadata.
- `panels/`: generated panels linked to notebooks and commits.
- `final/`: manuscript-ready composite figures.
- `figure_manifests/`: panel-to-input-to-script-to-commit mappings.

Each final figure should have a manifest containing:

```text
figure_id
panel
claim_id
notebook_id
input_files
input_checksums
generation_script
git_commit
generated_date
manual_edits
manual_editor
review_status
```
