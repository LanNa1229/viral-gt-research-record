# Repository Governance and Record Integrity

## Non-destructive history

- Do not use `git push --force` or `git push --force-with-lease` on shared research branches.
- Do not use `git rebase` to alter already shared or tagged research history.
- Do not amend a commit after it has been pushed and used as notebook evidence.
- Correct errors with a new commit that references the earlier commit.
- Do not delete tags associated with milestone snapshots.

## Accounts and access

- Use individual accounts; do not share credentials.
- Use role-appropriate access.
- Keep at least one UGA-accessible institutional copy.
- Record collaborators with repository access and the date access was granted or removed.
- Store access-control records without exposing secrets.

## Signing

Where supported and institutionally appropriate, sign commits and annotated tags with an individual GPG or SSH signing key. Signing can strengthen integrity checking but does not independently prove scientific priority or authorship.

## Sensitive information

Do not commit:

- passwords, API keys, tokens, SSH private keys;
- protected human-subject information;
- export-controlled or CUI data outside an approved environment;
- third-party data whose agreement prohibits the selected repository;
- unpublished collaborator data without authorization.

## Record interpretation

Git records file content, account identity, parent-child commit relationships, and timestamps supplied by the committing system. It is useful provenance but is not infallible and does not by itself prove who first conceived an idea.

Strong provenance relies on triangulation among:

- Git commits and tags;
- institutional snapshots;
- computational notebooks;
- HPC logs;
- emails and presentations;
- manuscript versions;
- file checksums;
- collaborator review;
- claim and contribution ledgers.
