# Portable execution

The original scientific files are unchanged. `scripts/reproduce.py` supplies current paths and recorded arguments, refuses existing output directories and records commands/versions. `scripts/validate_execution.py` checks scientific interfaces and completed runs. Independent reference validation and structure/checksum helpers are separate from analysis code.

Use Linux/WSL, Python 3.11 and the tested dependency lock. Runtime paths in public instructions are user-selected; historical absolute paths remain only in archival configurations/provenance and are replaced by the runner. `--community-column community` and the explicit 17-column rich M8 schema are necessary for direct execution with these exports. See `EXECUTION.md`.

The coordinate helper verifies every file and tests destination case sensitivity before extraction. The Figure 2B renderer is manuscript-specific (687 viral rows); missing Nimbus Sans affects typography only. The later environment produces a different fresh Leiden partition, so preserve the distinction between a successful fresh analysis and exact downstream replay from supplied assignments. No scientific compatibility patch was necessary.
