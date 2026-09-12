# Published-style BGC-Prophet embedding extraction

The original script was recovered from `/scratch/ln72030/BGC_Prophet/extract_bgc_embeddings.py` and archived here without modification on 2026-09-12. The recovered 4,555-byte file has SHA256 `33e6a9eb5d07e7634b3033d96667c450be1b3b7da9dcf8417600dbdd4ff3f616`, exactly matching the checksum recorded during the original validated extraction workflow. Local recovery source: `/home/rahnn/260908_MiBG/recovered_historical_scripts/extract_bgc_embeddings.py`.

The script instantiates the original pretrained BGC-Prophet `geneClassifier` using the supplied dataset, LMDB, and classifier checkpoint. Before running `clf.classify()`, it registers a temporary forward pre-hook on `clf.model.classifier`. This records the pooled 320-D Transformer representation immediately before the final product classifier; the hook is removed in a `finally` block.

The script does not modify BGC-Prophet inference; it records the classifier input during the normal forward pass. This is published-style BGC-Prophet embedding extraction, not a replacement model. The installed BGC-Prophet pipeline handles construction of the 128-position protein-embedding input from the dataset and LMDB; this wrapper does not manually assemble the 128 protein vectors. Its BGC-level output is one 320-D representation per input row, after the model's published all-128-position pooling.

Outputs are float32 embeddings, seven-class probabilities from the same inference run, and row-aligned metadata containing `embedding_row` and available `ID`, `labels`, and `isBGC` columns. The wrapper checks embedding dimensionality, finite values, and row count agreement. The historical full run produced 12,510 embedding rows, with validated shape (12510,320); probabilities had shape (12510,7). These are historical results, not a new inference run during archival recovery.

## Historical full-run command

The command below is the supplied historical full MIBiG run command. The historical Slurm/SBATCH launcher is not included in this update.

```bash
PY=/scratch/ln72030/conda_envs/bgcprophet/bin/python

$PY /scratch/ln72030/BGC_Prophet/extract_bgc_embeddings.py \
    --datasetPath /scratch/ln72030/BGC_Prophet/reference_data/BGC_train_dataset_classify.csv \
    --classifierPath /scratch/ln72030/BGC_Prophet/models/model/classifier.pt \
    --lmdbPath /scratch/ln72030/BGC_Prophet/reference_data/lmdb_train \
    --outputPath /scratch/ln72030/BGC_Prophet/reference_data/embedding_mibig_full \
    --name mibig_12510 \
    --device cpu \
    --batch_size 128
```

Expected filenames are `mibig_12510_embeddings.npy`, `mibig_12510_probabilities.npy`, and `mibig_12510_metadata.tsv`. Model, LMDB, environment, and arrays remain external; this archival update installs nothing and executes no inference.

See the [recovery manifest](../../../../metadata/manifests/DS-VGT-20260911-001_recovered_extraction.json) for checksum and source identity. Earlier reports and manifests are historical snapshots: their statement that this script was not recovered locally was accurate at first integration and is superseded by this recovery record.
