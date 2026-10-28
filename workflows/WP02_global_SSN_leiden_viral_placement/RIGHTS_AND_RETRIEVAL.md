# Structural inputs: extraction, attribution and redistribution

## Original coordinates

The companion dataset includes `inputs/structures/GT_A_4854_input_structures.zip`: 4,854 files (2,656 PDB, 2,198 CIF), 1,676,507,440 uncompressed bytes. ZIP size: 420,018,827 bytes; SHA256 `cd61755dfebc6a0b48c1795a996155c80de31e76e9d78e1bb18cb99769767c7c`. Original bytes, case and relative paths are preserved. `structures_manifest.tsv` records every ID, path, size, checksum and available source evidence.

Extract with `scripts/verify_structures.py` from the GitHub workflow (or `code_snapshot/scripts/verify_structures.py` from the data package). Use `--data-root` and a new `--extract-dir` on native Linux/WSL storage. The helper rejects case-insensitive destinations: 27 original filename pairs differ only by case. There are also 199 legacy `.pdb` metadata names whose actual originals are `.cif`; protein-ID stems agree. Both names are recorded without renaming original files.

## Source-specific rights and attribution

See [Third-party notices](THIRD_PARTY_NOTICES.md) and the per-file source/license table. Sources are classified as 420 AlphaFold DB, 2,194 locally generated AlphaFold 3 (versions 3.0.0–3.0.2), 3 experimental PDB and 2,237 unresolved. AlphaFold 3 output terms and the required binding notice are included. The complete archive remains local; unresolved-permission files must not be publicly redistributed.

## Optional Foldseek recomputation

Use the original commit `19800e921de4dcadd646feac40d83b4114137515` on PATH and adapt scheduler resources to your system. From the extracted data-package root:

```bash
python code_snapshot/scripts/verify_structures.py --data-root . \
  --extract-dir "$HOME/GT_A_4854_input" --report ../structures_verified.json
sbatch code_snapshot/original/foldseek/viral_all_vs_all_005_rich.sbatch \
  "$HOME/GT_A_4854_input/005_GTA_virus_and_others_curated_20260726" ../new_foldseek_results
```

This is an explicit expensive HPC submission, not a prerequisite for using supplied alignments. The script refuses to overwrite existing expected outputs. Foldseek createdb/all-pairs/structurealign/convertalis settings and the 17-column export schema are in Supplementary Methods. Apply the documented ID-normalization helper as appropriate and verify IDs against metadata. No Foldseek recomputation was run during release preparation.
