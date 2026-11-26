#!/usr/bin/env python3
"""Select direct representative and best giant-member hits for each control."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd


COLUMNS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen",
    "tlen", "qcov", "tcov", "alntmscore", "qtmscore", "ttmscore",
    "lddt", "prob",
]

NUMERIC = [c for c in COLUMNS if c not in {"query", "target"}]
EXCLUDED_DUPLICATE_IDS = {"Q5UQC3_model"}


def best_hit(frame: pd.DataFrame) -> pd.Series | None:
    if frame.empty:
        return None
    ranked = frame.sort_values(
        ["min_qt_tm", "alntmscore", "lddt", "qcov", "tcov", "bits", "evalue"],
        ascending=[False, False, False, False, False, False, True],
        kind="mergesort",
    )
    return ranked.iloc[0]


def selected_table(
    hits: pd.DataFrame,
    control_ids: list[str],
    selector,
    comparison_type: str,
) -> pd.DataFrame:
    rows = []
    for control_id in control_ids:
        row = best_hit(selector(hits.loc[hits["query"] == control_id]))
        if row is None:
            empty = {column: pd.NA for column in hits.columns}
            empty["query"] = control_id
            empty["comparison_type"] = comparison_type
            empty["status"] = "missing"
            rows.append(empty)
        else:
            out = row.to_dict()
            out["comparison_type"] = comparison_type
            out["status"] = "found"
            rows.append(out)
    columns = ["comparison_type", "status"] + list(hits.columns)
    return pd.DataFrame(rows)[columns]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-hits", type=Path, required=True)
    parser.add_argument("--baseline-membership", type=Path, required=True)
    parser.add_argument("--control-ids", type=Path, required=True)
    parser.add_argument("--representative", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    controls = [line.strip() for line in args.control_ids.read_text().splitlines() if line.strip()]
    if len(controls) != len(set(controls)):
        raise ValueError("control_ids.txt contains duplicate IDs")

    hits = pd.read_csv(args.all_hits, sep="\t", header=None, names=COLUMNS)
    if hits.shape[1] != len(COLUMNS):
        raise ValueError(f"Expected {len(COLUMNS)} columns in {args.all_hits}")
    duplicate_queries = sorted(set(hits["query"]) & EXCLUDED_DUPLICATE_IDS)
    if duplicate_queries:
        raise ValueError(
            "Excluded duplicate control structure(s) were present in the query directory: "
            f"{duplicate_queries}. Remove Q5UQC3_model and retain Q5UQC3_apo_model only."
        )
    for column in NUMERIC:
        hits[column] = pd.to_numeric(hits[column], errors="coerce")
    hits = hits.loc[hits["query"].isin(controls)].copy()
    hits["min_qt_tm"] = hits[["qtmscore", "ttmscore"]].min(axis=1)

    membership = pd.read_csv(
        args.baseline_membership,
        sep="\t",
        header=None,
        names=["representative", "member"],
        dtype=str,
    )
    giant_members = set(
        membership.loc[membership["representative"] == args.representative, "member"]
    )
    if not giant_members:
        raise ValueError(f"No baseline giant members found for {args.representative}")

    representative = selected_table(
        hits,
        controls,
        lambda df: df.loc[df["target"] == args.representative],
        "representative",
    )
    best_giant = selected_table(
        hits,
        controls,
        lambda df: df.loc[df["target"].isin(giant_members)],
        "best_giant_member_including_representative",
    )
    best_nonrep = selected_table(
        hits,
        controls,
        lambda df: df.loc[
            df["target"].isin(giant_members) & (df["target"] != args.representative)
        ],
        "best_nonrepresentative_giant_member",
    )

    representative.to_csv(
        args.output_dir / "control_representative_hits.tsv", sep="\t", index=False
    )
    best_giant.to_csv(
        args.output_dir / "control_best_giant_hits.tsv", sep="\t", index=False
    )
    best_nonrep.to_csv(
        args.output_dir / "control_best_nonrep_giant_hits.tsv", sep="\t", index=False
    )
    figure_input = pd.concat([representative, best_nonrep], ignore_index=True)
    figure_input.to_csv(
        args.output_dir / "coverage_figure_input.tsv", sep="\t", index=False
    )

    qc_rows = []
    for control_id in controls:
        control_hits = hits.loc[hits["query"] == control_id]
        rep_status = representative.loc[representative["query"] == control_id, "status"].iloc[0]
        nonrep_status = best_nonrep.loc[best_nonrep["query"] == control_id, "status"].iloc[0]
        qc_rows.append(
            {
                "control_id": control_id,
                "all_baseline_hits": len(control_hits),
                "baseline_giant_hits": int(control_hits["target"].isin(giant_members).sum()),
                "representative_comparison": rep_status,
                "best_nonrepresentative_comparison": nonrep_status,
            }
        )
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(args.output_dir / "control_search_qc.tsv", sep="\t", index=False)

    missing_rep = qc.loc[qc["representative_comparison"] != "found", "control_id"].tolist()
    missing_nonrep = qc.loc[
        qc["best_nonrepresentative_comparison"] != "found", "control_id"
    ].tolist()
    print(qc.to_string(index=False))
    if missing_rep or missing_nonrep:
        print(f"Missing representative comparisons: {missing_rep}", file=sys.stderr)
        print(f"Missing best non-representative comparisons: {missing_nonrep}", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
