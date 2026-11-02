# Reproduction and environment

The successful observed environment was existing Ubuntu/WSL with Python 3.12.11, pandas 2.3.3, NumPy 2.4.6, Matplotlib 3.10.8, Pillow 12.0.0 and Matplotlib FreeType 2.6.1. The original PNG records Matplotlib 3.10.8. This is a tested matching rendering environment, not a claim that the entire June environment was recovered. Font/build differences can affect raster bytes even with equal Python package versions.

`environment.yml` records the tested key versions. Additional requests is only for inspecting/running the archived API script; it is deliberately unnecessary for offline reproduction. Do not invoke the API runner for this workflow.

## Offline portable workflow

Run `python -B reproduce.py --out /new/output/path --require-exact-render` from this workflow directory in a matching Linux environment. A data-only portability check can omit `--require-exact-render`; inspect `validation.json` to distinguish exact image comparison from data equality. The code reads each preserved top-hit table by candidate-relative path, calls the recovered functions without edits, retains historical analysis_dir strings for comparison and runs the recovered plot CLI. The original evidence table is preserved byte-for-byte; regenerated evidence is compared as a multiset because set iteration changes row ordering.

The CLI below is also recovered and can be used with a new manifest whose **analysis_dir only** is redirected to the preserved per-locus `final/analysis_v2` folders. Retain candidate IDs and focal query IDs exactly, including known mismatches:

```bash
python -B originals/260630_Temp/build_cluster8_neighborhood_matrix.py --manifest /new/portable_manifest.tsv --outdir /new/matrix_output
python -B originals/260630_Temp/plot_cluster8_panelB_heatmap_clean.py --matrix originals/260630_Temp/cluster8_global_neighborhood_panelB/cluster8_neighborhood_matrix.tsv --outdir /new/figure_output
```

The new manifest's paths will appear in the generated matrix and differ from preserved historical paths; scientific fields should agree. The supplied `reproduce.py` avoids this metadata difference without editing inputs. No reconstruction of builder/plotter code occurred; the new script is only an offline validation wrapper.

SHA-256 manifests refer to paths relative to this workflow. `SHA256SUMS-originals.txt` covers original bytes; `SHA256SUMS-package.txt` covers all package files except itself. The outer ZIP has its own sidecar checksum. `.gitattributes` disables text conversion within this workflow so Git does not change archived CRLF/LF bytes. The complete recovery includes 346 files and about 315 MB of original bytes, including intermediate/duplicate historical tables; this broad preservation is intentional. No single source file exceeds 7 MB.

Executed validation is in `validation/validation.json`; source preservation in `audit/source_preservation_final.json`; the independent saved-evidence maxima check covers all 363 cells. Original PDF metadata is never modified: CreationDate normalization is comparison-only.
