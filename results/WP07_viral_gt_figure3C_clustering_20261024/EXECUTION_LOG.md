# Execution log index

Exact commands, cwd, code/configuration/input hashes, runtime versions and fingerprints were saved before fitting. Child stdout/stderr, source/input integrity and actual child exit are in each run. Raw wrapper exits and setup commands are additionally preserved in logs/codex_command_records.jsonl, indexed by logs/codex_command_index.csv.

| Run | Child exit | Sources unchanged | Inputs unchanged | Interpretation |
|---|---:|---|---|---|
| run_001_science | 0 | True | True | Successful recorded execution |
| run_002_figures | 0 | False | True | Child succeeded; parent wrapper failed its source-integrity guard |
| run_003_validation | 1 | True | True | Failed child; source snapshot/stderr retained |
| run_004_figures | 0 | True | True | Successful recorded execution |
| run_005_validation | 0 | True | True | Successful recorded execution |
| run_006_documentation | 0 | True | True | Successful recorded execution |
| run_007_supplementary_stability | 0 | True | True | Successful recorded execution |

- run_001_science is the sole clustering-fitting execution: 486 Leiden grid/baseline fits plus 12 multiplicity-balanced fits, and 24 HDBSCAN fits. Its inputs and executable code remained unchanged.
- run_002_figures produced initial exports, but the wrapper exited 1 after an unused validator was edited. Both code snapshots and exact initial assets remain under figure_exports/. This is not relabeled a clean success.
- run_003_validation exited 1 on syntax before checks. Its failed source and stderr remain.
- run_004_figures is the clean final render, with child exit 0 and unchanged sources/inputs. No scientific clustering or UMAP refit occurred.
- run_005_validation passed 51 checks without refitting clusters.
- run_006_documentation authored the initial Methods/readout from saved records.
- run_007_supplementary_stability compared saved partitions across matching grid settings, without new fits.

Setup attempts are in logs/setup_install*, including the UNC wheel-discovery and sandbox-temp failures followed by successful isolated installation. The Windows pip path normalization wrote initial wheel downloads to a lowercase UNC cache; identical wheel bytes were copied to the intended analysis path and a native Windows runtime location. No original analytical files were removed or changed. Runtime wheels and installed environments are external; exact names/hashes are retained in runtime/WHEEL_MANIFEST.tsv and the pip report.

The original UMAP effective epoch count remains NOT RECORDED. It is not a parameter of this new clustering run. Current clustering iteration counts and actual parameters are recorded, and functional identities are intentionally not assigned. The inherited full environment is not dumped because it can contain credentials; explicit numerical settings, package inventory and relevant code/native identities are retained.

[Leiden objective/API](https://leidenalg.readthedocs.io/en/stable/reference.html) and [HDBSCAN convention](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.HDBSCAN.html) are supporting documentation; actual execution records take precedence over live documentation defaults.
