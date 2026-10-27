# Step 2: viral BGC-Prophet inference and contextual 320D embeddings

Local v1.0.0 release candidate; publication approval pending. Full frozen accounting is verified; representative cold-cache reproduction passed for six paired occurrences (12 windows), with zero score or embedding differences. Full-cohort inference has not been rerun.

This directory contains new portability, verification and documentation only. Already published historical scripts, Step 1 inputs, 15-gene outputs and the matched reference bank are linked by immutable commit/DOI and SHA256 in resources.json and references/public_file_reuse.tsv. Historical source bytes are not duplicated or edited. The separate Step 2 dataset contains only six missing historical 11-gene files plus new release evidence/metadata. No upstream source or model weight is redistributed here.

## Scientific method

Step 1 CDS-table sequences, sorted by relative position, become ESM2 t6/8M layer-6 residue means excluding BOS/EOS; proteins above 1,022 residues use their first 1,022 residues. The frozen BGC-Prophet classifier uses two 320D Transformer layers, five heads, feedforward size 1,280 and a final encoder normalization. The contextual window embedding is the mean of actual gene positions immediately before the linear/sigmoid product-class head. All-one float32 membership is an additive key-attention bias, not Boolean padding exclusion. Actual n genes are supplied; no padding or repetition to 128 occurs. The source computes positional encoding but passes the original src into its encoder; this historical behavior is preserved.

Saccharide uses its original sigmoid score and threshold=max(0.5, Other score). Other is removed if another label qualifies; Unknown remains an analyzed model outcome. These predictions and embedding associations do not verify glycan biosynthesis.

| Frozen accounting | 15-gene | 11-gene |
|---|---:|---:|
| Extracted windows | 757 | 757 |
| Successful predictions | 737 | 740 |
| Saccharide-positive | 246 | 265 |
| Excluded before inference | 20 | 17 |
| Embedding shape | 737 x 320 | 740 x 320 |

## Reproduction and evidence

Read EXECUTION_GUIDE.md and ENVIRONMENT.md. portable/step2.py composes the checksum-pinned original encode_missing, load_model, infer and predicted functions. It replaces absolute-path/cache plumbing only. The validated cold-cache alternative reconstructs protein vectors from preserved sequences, without historical pilot caches or LMDB. Six fixed pairs cover one-gene and five-gene inputs, full positive/negative windows, a long protein and multilabel output.

validation/representative_regression.json, validation/frozen_accounting.json and tests/test_release.py record/check scientific correspondence. The full-cohort infer command is separately gated and has not been executed or validated across all windows. Do not describe this candidate as a complete standalone full-cohort scientific reproduction. See HISTORICAL_LIMITATIONS.md.

Original creator: Lan Na (given Lan; family Na). Proposed MIT covers original software; proposed CC BY 4.0 covers original derived data/documentation. Upstream terms remain separate; see LICENSE_PROPOSAL.md and THIRD_PARTY_NOTICES.md. No Step 2 public commit or dataset DOI is assigned yet. Step 1 dataset: https://doi.org/10.5281/zenodo.23680510 .
