# Retrieval-only execution

Use Python >=3.11 and NumPy 1.26.4. An existing environment was used; no new environment was installed. Do not execute historical scripts as programs: their mains contain unrelated candidate, background, source-validation or model-dependent work.

Extract the future Step 3 dataset to a separate directory, such as /tmp/step3-data. From workflows/viral_gt_saccharide_retrieval:

~~~bash
python -B portable/step3.py fetch --cache /tmp/step3-public-cache
python -B portable/step3.py verify --cache /tmp/step3-public-cache --dataset /tmp/step3-data
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python -B portable/step3.py reproduce --cache /tmp/step3-public-cache --dataset /tmp/step3-data --output /tmp/step3-reproduced
~~~

The output directory must not exist and must be separate from package/data/cache directories. Fetch is the only network stage. Reproduce checks hashes, AST-loads only the original paired.py reference_metrics function, and independently checks its matrix calculation with a per-query normalized dot product. Original module imports and main functions are never executed.

Reuse existing official archives with:

~~~bash
python -B portable/step3.py fetch --cache /tmp/step3-public-cache --archive step1=/absolute/path/Step1_viral_GT_neighborhoods_v1.0.0_DATA_FINAL.zip --archive step2=/absolute/path/Step2_frozen_11gene_outputs_v1.0.0_DATA.zip
~~~

Entire archive hashes are checked before extracting selected members. Partial/changed files cause failure and remain for inspection. Step 1 is about 74 MB; Step 2 about 1.34 MB. No model weights or large databases are downloaded.

Outputs include reconstructed retrieval/neighbor tables, per-window comparisons, reproduced 18-row historical accounting, and retrospective >= versus > sensitivity. They do not replace frozen candidates or figures.

~~~bash
python -B -m unittest discover -s tests -v
~~~

Exact assertions: reference membership/order and SHA256, row correspondence, ranks/IDs/labels, P@5/P@10/P@20, screening counts. Numerical assertions: cosine/enrichment values at atol=rtol=1e-12. Output formatting/line endings and table byte identity are not reproduction criteria; original archived files remain byte-preserved.

All inputs resolve through resources.json and references/public_resources.tsv. Original absolute paths are provenance, not adapter dependencies. No hidden ESM cache, LMDB, checkpoint, taxonomy database, graphical session or Git checkout is required.
