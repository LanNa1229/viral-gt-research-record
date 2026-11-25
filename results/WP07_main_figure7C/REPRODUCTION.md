# Reproduce the Main Figure 7C baseline and paired plot

From this package directory, use Python 3.12 with the recorded package versions. The historical graph runtime used NumPy 2.3.5, SciPy 1.18.1, scikit-learn 1.9.0, igraph 1.0.0 and leidenalg 0.12.0. The plotting runtime used Matplotlib 3.11.1 and Pillow 12.3.0. See `provenance/HISTORICAL_EXECUTION.json` for the exact font hashes and software versions. Supply licensed local regular and bold Arial files for the closest render. Fonts are separate runtime dependencies.

```sh
python -m pip install numpy==2.3.5 scipy==1.18.1 scikit-learn==1.9.0 igraph==1.0.0 leidenalg==0.12.0 matplotlib==3.11.1 Pillow==12.3.0
python -B reproduce.py ./scratch-clustering --stage clustering
python -B reproduce.py ./scratch-plot --stage plotting --arial /path/to/arial.ttf --arial-bold /path/to/arialbd.ttf
```

Choose fresh scratch directories. The default restorer retrieves immutable GitHub blobs anonymously and checks their SHA-256 digests before use. For offline operation, repeat `--checkout SHA=/path/to/checkout` for both commits in `inputs/PUBLIC_DEPENDENCIES.tsv`. Optional `--deps` flags add existing runtime package directories.

The clustering replay calls the original executed engine's normalization, graph and optimizer functions on the row-aligned historical subset. The recorded historical driver ran a larger project. The portable replay runs only its selected-287 baseline, checks every edge and weight, and compares all 287 canonical assignments. Computation on another numerical stack may differ at floating-point ties or in optimizer behavior. Exact historical graph and assignment files are preserved for comparison.

The plotting replay uses the byte-identical historical script's paired figure branch and its saved configuration, with portable Arial paths. The script originally rendered two alternatives in one execution. The wrapper invokes only the paired branch and writes PDF, SVG and PNG in scratch. The original vector files and PNG remain in `figures/final/`. PDF/SVG serialization and creation metadata may differ between runs.

The unified figure table is the existing historical output with local absolute paths removed from its provenance field. All biological and display fields are unchanged. The replay independently checks joins against public files and the baseline assignment table. The historical integration script is retained as execution evidence, and the portable entry point performs the independently keyed checks. No private GitHub access is required.
