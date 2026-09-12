#!/usr/bin/env python

import argparse
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd

from bgc_prophet.command.classify import geneClassifier


def main():
    parser = argparse.ArgumentParser(
        description="Export the 320-D BGC-Prophet representation immediately before the product classifier."
    )

    parser.add_argument("--datasetPath", required=True, type=Path)
    parser.add_argument("--classifierPath", required=True, type=Path)
    parser.add_argument("--lmdbPath", required=True, type=Path)
    parser.add_argument("--outputPath", required=True, type=Path)
    parser.add_argument("--name", default="bgc_embedding")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch_size", type=int, default=512)

    args = parser.parse_args()

    args.outputPath.mkdir(parents=True, exist_ok=True)

    classifier_args = Namespace(
        datasetPath=args.datasetPath,
        classifierPath=args.classifierPath,
        lmdbPath=args.lmdbPath,
        outputPath=args.outputPath,
        device=args.device,
        batch_size=args.batch_size,
        name=args.name,
        classify_t=0.5,
    )

    # Load the original pretrained BGC-Prophet classifier.
    clf = geneClassifier(classifier_args)

    captured_embeddings = []

    # The input to clf.model.classifier is the pooled
    # 320-D Transformer representation that we want.
    def capture_classifier_input(module, inputs):
        x = inputs[0]
        captured_embeddings.append(
            x.detach().cpu().numpy().astype(np.float32)
        )

    # Attach a temporary "recorder" immediately before classifier.
    hook = clf.model.classifier.register_forward_pre_hook(
        capture_classifier_input
    )

    try:
        # Run the ORIGINAL BGC-Prophet inference.
        clf.classify()
    finally:
        hook.remove()

    embeddings = np.concatenate(captured_embeddings, axis=0)

    # ----------------------------
    # Basic validation
    # ----------------------------
    if embeddings.ndim != 2:
        raise RuntimeError(
            f"Expected 2-D embeddings, got shape {embeddings.shape}"
        )

    if embeddings.shape[1] != 320:
        raise RuntimeError(
            f"Expected 320 dimensions, got {embeddings.shape[1]}"
        )

    if np.isnan(embeddings).any():
        raise RuntimeError("NaN detected in embeddings.")

    if np.isinf(embeddings).any():
        raise RuntimeError("Inf detected in embeddings.")

    # ----------------------------
    # Save embeddings
    # ----------------------------
    embedding_file = (
        args.outputPath / f"{args.name}_embeddings.npy"
    )

    np.save(
        embedding_file,
        embeddings.astype(np.float32)
    )

    # Also save the 7-class probabilities generated in
    # the same inference run. We will use these for regression checking.
    probability_file = (
        args.outputPath / f"{args.name}_probabilities.npy"
    )

    np.save(
        probability_file,
        clf.results.astype(np.float32)
    )

    # ----------------------------
    # Preserve row IDs / metadata
    # ----------------------------
    df = pd.read_csv(args.datasetPath)

    metadata_cols = [
        col for col in [
            "ID",
            "labels",
            "isBGC",
        ]
        if col in df.columns
    ]

    if metadata_cols:
        metadata = df[metadata_cols].copy()
    else:
        metadata = pd.DataFrame(
            {"original_row": np.arange(len(df))}
        )

    metadata.insert(
        0,
        "embedding_row",
        np.arange(len(metadata))
    )

    if embeddings.shape[0] != len(metadata):
        raise RuntimeError(
            f"Row mismatch: {embeddings.shape[0]} embeddings "
            f"but {len(metadata)} metadata rows."
        )

    metadata_file = (
        args.outputPath / f"{args.name}_metadata.tsv"
    )

    metadata.to_csv(
        metadata_file,
        sep="\t",
        index=False
    )

    print()
    print("=== BGC-Prophet embedding export complete ===")
    print(f"Samples:          {embeddings.shape[0]}")
    print(f"Embedding shape:  {embeddings.shape}")
    print(f"dtype:            {embeddings.dtype}")
    print(f"NaN:              {np.isnan(embeddings).any()}")
    print(f"Inf:              {np.isinf(embeddings).any()}")
    print()
    print(f"Embeddings:       {embedding_file}")
    print(f"Metadata:         {metadata_file}")
    print(f"Probabilities:    {probability_file}")


if __name__ == "__main__":
    main()
