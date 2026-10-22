#!/usr/bin/env python3
"""Write strict 100 / supported-only 60 / other 0 without changing PDB coordinates.

Added for the 2026-10-21 research record. This is not the July pan-core caller.
"""
from __future__ import annotations
import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


def truth(value: str) -> bool:
    return str(value).lower() in {"true", "1", "yes"}


def build(pdb_path: Path, mapping_path: Path, output_path: Path) -> dict:
    classes = {}
    with mapping_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            label = row.get("query_residue_label", "")
            if not label or label == "NA":
                continue
            chain, identifier = label.split(":", 1)
            hit = re.fullmatch(r"(-?\d+)([A-Za-z]?)", identifier)
            if not hit:
                raise ValueError(f"Unsupported residue label {label!r}")
            key = (chain, int(hit[1]), hit[2] or " ")
            value = 100 if truth(row["strict_pan_core"]) else 60 if truth(row["supported_pan_core"]) else 0
            if key in classes and classes[key] != value:
                raise ValueError(f"Conflicting class assignments for {key}")
            classes[key] = value

    source = pdb_path.read_text().splitlines(keepends=True)
    result, ca_keys, counts = [], set(), Counter()
    for line in source:
        if line.startswith(("ATOM  ", "HETATM")):
            if len(line.rstrip("\n")) < 66:
                raise ValueError("PDB atom line is shorter than the B-factor field")
            key = (line[21], int(line[22:26]), line[26])
            value = classes.get(key, 0)
            if line[12:16].strip() == "CA" and line.startswith("ATOM  "):
                if key in ca_keys:
                    raise ValueError(f"Duplicate CA residue {key}; resolve alternate models first")
                ca_keys.add(key)
                counts[str(value)] += 1
            modified = line[:60] + f"{value:6.2f}" + line[66:]
            assert line[:60] == modified[:60] and line[66:] == modified[66:]
            result.append(modified)
        else:
            result.append(line)
    missing = {key for key, value in classes.items() if value > 0} - ca_keys
    if missing:
        raise ValueError(f"Mapped core residues absent from PDB: {sorted(missing)}")
    for value in [100, 60]:
        expected = sum(x == value for x in classes.values())
        if counts[str(value)] != expected:
            raise AssertionError((value, counts[str(value)], expected))
    if output_path.resolve() == pdb_path.resolve():
        raise ValueError("Use a new output path so original pLDDT is retained")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(result))
    return {"n_residues": len(ca_keys), "CA_class_counts": dict(counts),
            "only_B_factor_field_changed": True, "output": str(output_path.resolve())}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pdb", type=Path, required=True)
    p.add_argument("--mapping", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(build(args.pdb, args.mapping, args.output), indent=2))


if __name__ == "__main__":
    main()
