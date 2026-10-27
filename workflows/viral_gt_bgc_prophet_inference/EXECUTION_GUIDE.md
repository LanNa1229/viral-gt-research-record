# Execution guide

Use Python 3.12.11, torch 2.5.1+cpu, numpy 1.26.4 and fair-esm 2.0.0. Reuse the existing validated environment when available. No packages or environments were installed during release preparation. Commands below are user-run recipes, not evidence that every listed operation was executed.

## Prepare checksum-pinned resources

From this workflow directory:

```sh
PY=/absolute/path/to/the/documented/python
CACHE=/absolute/path/to/NEW_step2_resource_cache
DATA2=/absolute/path/to/extracted_Step2_dataset
DATA1=/absolute/path/to/Step1_viral_GT_neighborhoods_v1.0.0_DATA_FINAL.zip
$PY -B portable/step2.py prepare --cache "$CACHE" --step1-archive "$DATA1" --step2-data "$DATA2"
$PY -B portable/step2.py verify --cache "$CACHE"
```

The prepare command retrieves only pinned public code, frozen 15-gene records, the existing reference/proof files and three small original model resources. Reuse the Step 1 ZIP instead of downloading it again. Only named members are extracted into the cache. Missing or mismatched bytes fail closed; existing mismatched resources are not overwritten. Add --offline when all external files already exist. resources.json lists every URL, exact commit/archive member and SHA256. Upstream author downloads are used for weights; no weights are shipped in the Step 2 dataset.

For a genuinely empty cache, omit --offline; the command can retrieve the Step 1 archive if --step1-archive is omitted. The Step 2 dataset must still be supplied locally until its DOI is assigned. The canonical BGC-Prophet model asset is an 18.4 MB tar file; only its unique classifier.pt member is retained. ESM2 is approximately 30.1 MB plus its 1,511-byte contact-regression sidecar. No large reference LMDB is needed for this validated path.

## Representative reproduction: tested

```sh
OUT=/absolute/path/to/NEW_step2_regression
$PY -B portable/step2.py regress --cache "$CACHE" --output "$OUT"
$PY -B tests/test_release.py --cache "$CACHE" --regression "$OUT" -v
```

This performs six paired occurrences, 12 window results and 66 unique cold-start protein encodings in the recorded preparation. Each BGC input is run once without hooks and once with the original recording hooks, so 12 window results involve 24 classifier forward calls. Full-cohort inference is not triggered. Compare every selected seven-score vector and every 320D coordinate to frozen records, require unchanged labels, exact actual/membership-pooling equality and unchanged model tensors. Tolerance is 1e-6; the observed maximum error was zero for both scores and embeddings. Existing output directories are refused.

## Full 15-gene and 11-gene replay: implemented, not executed

The following commands require separate approval before an agent executes them:

```sh
$PY -B portable/step2.py infer --cache "$CACHE" --size 15 --output /absolute/path/to/NEW_full15 --full-cohort
$PY -B portable/step2.py infer --cache "$CACHE" --size 11 --output /absolute/path/to/NEW_full11 --full-cohort
```

These compose the same pinned historical scientific functions and write separate raw prediction, embedding, manifest and QC files. They do not retrain, run a detector, change thresholds, alter candidate sets, or update historical directories. No full replay was performed in this release preparation. The new full-cohort orchestration/serialization is not covered by a full numerical regression; compare its outputs before claiming full reproduction. The explicit flag protects against accidental expensive execution.

## Historical launchers

The original 15-gene run used run_full.py --root ROOT --task TASK --run-name baseline_594_v2 with PYTHONPATH pointing at the original checkout/src. The original 11-gene experiment uses experiment.py execute, including its gate and infer-main child stages. These files remain at their pinned Step 1 GitHub locations. Those broad historical launchers depend on old workspace layout and cached proofs; do not invoke them from this portable directory. The new adapter replaces that plumbing while preserving the scientific functions. No literal original shell transcript or historical full environment lock has been recovered.
