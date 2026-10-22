#!/usr/bin/env python3
"""Reproduce the three F2Y1U4 versus GT31 summary bars from archived inputs.

Python standard library only. Added 22 October 2026; implements the identity
definition preserved in the 21 October audit of the July analysis.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal, ROUND_HALF_UP
import hashlib
import json
import math
from pathlib import Path
import platform
from statistics import mean, median

ROOT = Path(__file__).resolve().parents[1]
QUERY = "F2Y1U4"
GROUP_SIZES = {"SC01": 5, "SC02": 4, "SC03": 4, "SC04": 4}
METRICS = {
    "sequence_identity": "median_AA_identity",
    "query_coverage": "median_query_coverage",
    "alignment_TM_score": "median_Foldseek_AlnTM",
}


def check(condition, message):
    if not condition:
        raise ValueError(message)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def verify_checksums():
    for name in ["SHA256SUMS", "INPUT_SHA256SUMS"]:
        for line in (ROOT / name).read_text(encoding="utf-8").splitlines():
            digest, relative = line.split("  ", 1)
            check(hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest,
                  f"Checksum mismatch: {relative}")


def read_alignment(path):
    sequences = {}
    name = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            check(name not in sequences, f"Duplicate FASTA ID: {name}")
            sequences[name] = ""
        elif line.strip():
            check(name is not None, "Sequence before FASTA header")
            sequences[name] += line.strip()
    check(len(sequences) == 20, "Expected archived 19-reference-plus-query alignment")
    check(len({len(s) for s in sequences.values()}) == 1, "Unequal alignment lengths")
    return sequences


def calculate():
    members = read_tsv(ROOT / "inputs/reference_members.tsv")
    ids = {r["protein_id"] for r in members}
    check(len(members) == len(ids) == 17 and QUERY not in ids,
          "Expected 17 unique reference proteins, excluding the query")
    check(Counter(r["cluster_id"] for r in members) == GROUP_SIZES,
          "Reference group membership changed")
    alignment = read_alignment(
        ROOT / "inputs/alignments/C012_GT31_reference19_plus_F2Y1U4_aa.fa")
    check(ids | {QUERY} <= alignment.keys(), "Missing alignment sequences")
    raw = read_tsv(ROOT / "inputs/foldseek/F2Y1U4_vs_reference17.tsv")
    check(len(raw) == 17 and {r["query"] for r in raw} == {QUERY}
          and {r["target"] for r in raw} == ids,
          "Expected exactly one query-directed Foldseek row per reference")
    foldseek = {r["target"]: r for r in raw}
    individual = []
    for member in sorted(members, key=lambda r: (r["cluster_id"], r["protein_id"])):
        pid = member["protein_id"]
        aligned = [(a, b) for a, b in zip(alignment[QUERY], alignment[pid])
                   if a not in "-." and b not in "-."]
        check(bool(aligned), f"No mutually nongap columns: {pid}")
        matches = sum(a == b for a, b in aligned)
        qcov = float(foldseek[pid]["qcov"])
        alntm = float(foldseek[pid]["alntmscore"])
        check(math.isfinite(qcov) and 0 <= qcov <= 1, f"Invalid qcov: {pid}")
        check(math.isfinite(alntm) and alntm >= 0, f"Invalid AlnTM: {pid}")
        individual.append({"query": QUERY, "reference": pid,
            "unit": member["cluster_id"], "AA_matches": matches,
            "AA_both_nongap_columns": len(aligned), "AA_identity": matches / len(aligned),
            "query_coverage": qcov, "Foldseek_AlnTM": alntm})
    groups = []
    decimal_medians = {metric: [] for metric in METRICS}
    for group in GROUP_SIZES:
        rows = [r for r in individual if r["unit"] == group]
        groups.append({"query": QUERY, "unit": group, "reference_n": len(rows),
            "median_AA_identity": median(r["AA_identity"] for r in rows),
            "median_query_coverage": median(r["query_coverage"] for r in rows),
            "median_Foldseek_AlnTM": median(r["Foldseek_AlnTM"] for r in rows)})
        decimal_medians["sequence_identity"].append(median(
            Decimal(r["AA_matches"]) / Decimal(r["AA_both_nongap_columns"]) for r in rows))
        for metric, source in [("query_coverage", "qcov"),
                               ("alignment_TM_score", "alntmscore")]:
            decimal_medians[metric].append(median(
                Decimal(foldseek[r["reference"]][source]) for r in rows))
    bars = []
    for metric, column in METRICS.items():
        value = mean(r[column] for r in groups)
        is_percent = metric != "alignment_TM_score"
        # Decimal aggregation avoids binary tie-rounding at AlnTM = 0.79955.
        decimal_mean = sum(decimal_medians[metric]) / Decimal(4)
        display = (decimal_mean * (100 if is_percent else 1)).quantize(
            Decimal("0.01" if is_percent else "0.0001"), rounding=ROUND_HALF_UP)
        bars.append({"metric": metric, "mean_subcluster_median": value,
                     "display_value": str(display),
                     "display_unit": "%" if is_percent else "score"})
    return individual, groups, bars


def validate(individual, groups, bars):
    expected = ROOT / "expected/query_comparison"
    archived_individual = read_tsv(expected / "reference_identity.tsv")
    check(len(archived_individual) == 17, "Expected 17 archived identity rows")
    by_reference = {r["reference"]: r for r in archived_individual}
    check(set(by_reference) == {r["reference"] for r in individual}, "Reference IDs changed")
    for row in individual:
        old = by_reference[row["reference"]]
        check(row["query"] == old["query"] and row["unit"] == old["unit"],
              "Archived query/group differs")
        for key in ["AA_matches", "AA_both_nongap_columns", "AA_identity"]:
            check(math.isclose(row[key], float(old[key]), rel_tol=0, abs_tol=1e-12),
                  f"Identity differs from archived audit: {row['reference']} {key}")
    for filename, tolerance in [("subcluster_medians.tsv", 1e-12),
                                ("subcluster_medians_20260731.tsv", 0.00005001)]:
        old_rows = read_tsv(expected / filename)
        archived = {r["unit"]: r for r in old_rows}
        check(len(old_rows) == 4 and set(archived) == set(GROUP_SIZES),
              "Archived subcluster IDs changed")
        for row in groups:
            old = archived[row["unit"]]
            check(old["query"] == QUERY and int(old["reference_n"]) == row["reference_n"],
                  "Archived query or reference count differs")
            for key in METRICS.values():
                check(math.isclose(row[key], float(old[key]), rel_tol=0, abs_tol=tolerance),
                      f"Subcluster median differs from {filename}: {row['unit']} {key}")
    reported = json.loads((expected / "reported_bar_values.json").read_text())
    check({r["metric"]: r["display_value"] for r in bars} == reported,
          "Recalculated bars differ from the reported figure values")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=ROOT / "query_comparison_output")
    args = parser.parse_args()
    verify_checksums()
    out = args.outdir.resolve()
    check(not out.exists() or (out.is_dir() and not any(out.iterdir())),
          "Output directory must be new or empty; choose another --outdir")
    individual, groups, bars = calculate()
    validate(individual, groups, bars)
    out.mkdir(parents=True, exist_ok=True)
    for filename, rows in [("per_reference.tsv", individual),
                           ("subcluster_medians.tsv", groups), ("bar_summary.tsv", bars)]:
        write_tsv(out / filename, rows)
    record = {"status": "PASS", "n_references": 17, "n_subclusters": 4,
              "checksums_verified": True, "individual_identities_match_archive": True,
              "subcluster_medians_match_archive": True, "reported_bars_match": True,
              "aggregation": "unweighted mean of four subcluster medians",
              "bars": bars, "python": platform.python_version(),
              "upstream_foldseek_or_foldmason_rerun": False}
    (out / "validation.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
