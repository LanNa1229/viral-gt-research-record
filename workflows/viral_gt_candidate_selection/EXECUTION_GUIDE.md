# Reproduce selection without scientific inference

Use Python >=3.11 and only its standard library. No Conda environment, Torch, ESM, NumPy, taxonomy DB or network access during selection is needed.

Unzip the GitHub ZIP and companion data ZIP. From workflows/viral_gt_candidate_selection:

```sh
python3 -I -S -B portable/fetch_resources.py --out /tmp/step4-public-inputs
python3 -I -S -B tests/test_step4.py
python3 -I -S -B validation/run_isolated.py --cache /tmp/step4-public-inputs --dataset /absolute/path/to/unpacked-step4-data --out /tmp/step4-reproduction
```

To reuse existing published ZIPs, add --archive-cache /absolute/path/to/verified-archives to the fetch command. Required ZIP names and SHA256 values are in references/public_archives.tsv. The fetcher accepts only exact checksums and immutable public URLs; it does not use private historical inputs. The cache may contain published source sequences because grouping/exclusion verification needs them, but none is included in either new release archive.

The offline guard allows reads only from the extracted software, verified public cache, expected Step 4 dataset and Python runtime. It denies scientific-output writes outside the new output folder, and denies network/subprocess operations during reproduction. -I -S -B removes user-site/PYTHONPATH influence and bytecode writes.

Outputs: selection_genomes.tsv (466, 33 selection/grouping columns), selection_loci.tsv (565), positive_window_support.tsv (675 positive size-specific evidence rows), candidate_membership245.tsv, sankey_membership242.tsv, reproduction_validation.json and execution_guard.json. These are newly reproduced validation outputs; the canonical original metadata-rich master remains unchanged in the data package.

Expected status: PASS; 15,378 genome-field comparisons; 565 complete locus-row comparisons; 245/242 memberships exact; no discrepancies. See validation/reproduction_validation.json and validation/clean_execution_validation.json for the archive-extracted execution.

For an ordinary unguarded selection call, use portable/step4.py with the same --cache, --dataset and --out arguments. Never run the historical Figure B scripts as part of this workflow.
