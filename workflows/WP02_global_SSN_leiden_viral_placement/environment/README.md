# Tested scientific environment

Use Python 3.11 on Linux/WSL and `python -m pip install -r environment/requirements.lock.txt` in a virtual environment. The full lock records the installed dependency set; `requirements.reproduction.txt` lists the top-level scientific pins. `requirements.original-minimums.txt` is retained as the original dependency declaration, not a historical lock.

The tested runtime is Python 3.11.15, NumPy 2.2.6, pandas 2.3.3, SciPy 1.15.3, scikit-learn 1.7.2, igraph 0.11.9 and leidenalg 0.10.2; Matplotlib 3.10.7/seaborn 0.13.2 support optional plots. Installed requirements pass `pip check`, all imports and both scientific CLIs. Dependency versions are independently checked by `scripts/validate_execution.py`.

This compatible environment runs both scripts and exactly replays downstream classification from the primary partition. It is not asserted to be the original July environment. With unchanged global-minimum weights, resolution 1 and seed 1, its fresh partition differs from the stored primary partition (ARI 0.9451745569665067). The producing historical package versions were not recorded; changing pins or parameters to seek a desired partition is not part of the release. The compatible igraph 0.11.9 / leidenalg 0.10.2 combination is retained, with measured agreement rather than an unsupported claim of historical equivalence.

The stdlib-only checksum and structure-extraction helpers do not require scientific packages. The analysis itself needs the packages above. Foldseek is optional when using supplied alignments; a full recomputation requires a suitable HPC system and the documented Foldseek commit/settings. No all-versus-all run is performed by the Python reproduction runner.
