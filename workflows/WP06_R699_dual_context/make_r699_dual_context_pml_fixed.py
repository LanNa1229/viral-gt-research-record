#!/usr/bin/env python3
"""Generate a publication-ready PyMOL script from
R699_dual_context_classification.tsv.

Expected columns from map_dual_context_cores_to_r699.py:
  constraint_class, structure_chain, structure_resseq, structure_icode
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

CLASSES = [
    "both_supported",
    "C005_supported",
    "C009_supported",
    "unsupported",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", required=True, type=Path)
    p.add_argument("-s", "--structure", required=True, type=Path)
    p.add_argument("-o", "--output", required=True, type=Path)
    p.add_argument("--object-name", default="R699")
    p.add_argument("--force-chain", default=None)
    return p.parse_args()


def chunks(values: List[str], size: int = 100) -> Iterable[List[str]]:
    for i in range(0, len(values), size):
        yield values[i : i + size]


def load_residue_classes(path: Path, force_chain: str | None) -> Dict[str, List[Tuple[str, str]]]:
    out: Dict[str, List[Tuple[str, str]]] = defaultdict(list)

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"No TSV header found in {path}")

        required = {"constraint_class", "structure_resseq"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Missing required columns: {sorted(missing)}. "
                f"Observed columns: {reader.fieldnames}"
            )

        has_chain = "structure_chain" in reader.fieldnames
        has_icode = "structure_icode" in reader.fieldnames

        for row in reader:
            klass = (row.get("constraint_class") or "").strip()
            if klass not in CLASSES:
                continue

            raw_resseq = (row.get("structure_resseq") or "").strip()
            if not raw_resseq:
                continue

            # Handles values written as either integer strings or 123.0.
            try:
                resseq = str(int(float(raw_resseq)))
            except ValueError:
                resseq = raw_resseq

            icode = (row.get("structure_icode") or "").strip() if has_icode else ""
            residue_token = f"{resseq}{icode}" if icode else resseq

            if force_chain is not None:
                chain = force_chain
            elif has_chain:
                chain = (row.get("structure_chain") or "").strip() or "A"
            else:
                chain = "A"

            out[klass].append((chain, residue_token))

    for klass in CLASSES:
        out[klass] = sorted(set(out[klass]), key=lambda x: (x[0], x[1]))
    return out


def write_selection(handle, object_name: str, selection_name: str, residues: List[Tuple[str, str]]) -> None:
    if not residues:
        handle.write(f"select {selection_name}, none\n")
        return

    by_chain: Dict[str, List[str]] = defaultdict(list)
    for chain, resi in residues:
        by_chain[chain].append(resi)

    handle.write(f"select {selection_name}, none\n")
    for chain, resis in sorted(by_chain.items()):
        unique_resis = sorted(set(resis))
        for part in chunks(unique_resis):
            selector = "+".join(part)
            handle.write(
                f"select {selection_name}, {selection_name} or "
                f"({object_name} and chain {chain} and resi {selector})\n"
            )


def main() -> None:
    args = parse_args()
    if not args.input.exists():
        raise SystemExit(f"ERROR: missing TSV: {args.input}")
    if not args.structure.exists():
        raise SystemExit(f"ERROR: missing structure: {args.structure}")

    classes = load_residue_classes(args.input, args.force_chain)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    obj = args.object_name
    sel = {
        "both_supported": f"{obj}_both_supported",
        "C005_supported": f"{obj}_C005_supported",
        "C009_supported": f"{obj}_C009_supported",
        "unsupported": f"{obj}_unsupported",
    }

    with args.output.open("w", encoding="ascii") as handle:
        handle.write("# R699 dual-context structural-constraint mapping\n")
        handle.write("# Purple: C005 and C009 overlap\n")
        handle.write("# Orange: C005 only\n")
        handle.write("# Blue: C009 only\n")
        handle.write("# Light gray: unsupported\n\n")
        handle.write("reinitialize\n")
        handle.write(f"load {args.structure.resolve()}, {obj}\n")
        handle.write("hide everything, all\n")
        handle.write(f"show cartoon, {obj} and polymer.protein\n")
        handle.write("set_color both_purple, [0.482, 0.196, 0.580]\n")
        handle.write("set_color c005_orange, [0.902, 0.624, 0.000]\n")
        handle.write("set_color c009_blue, [0.000, 0.447, 0.698]\n")
        handle.write("set_color unsupported_gray, [0.850, 0.850, 0.850]\n\n")

        for klass in CLASSES:
            write_selection(handle, obj, sel[klass], classes[klass])
        handle.write("\n")

        # Color unsupported first, then overwrite with supported classes.
        handle.write(f"color unsupported_gray, {sel['unsupported']}\n")
        handle.write(f"color c005_orange, {sel['C005_supported']}\n")
        handle.write(f"color c009_blue, {sel['C009_supported']}\n")
        handle.write(f"color both_purple, {sel['both_supported']}\n\n")

        # Optional ligands/cofactors, if present.
        handle.write(f"select {obj}_ligands, {obj} and not polymer.protein and not solvent\n")
        handle.write(f"show sticks, {obj}_ligands\n")
        handle.write(f"util.cbag {obj}_ligands\n\n")

        handle.write("bg_color white\n")
        handle.write("set orthoscopic, on\n")
        handle.write("set ray_opaque_background, off\n")
        handle.write("set antialias, 2\n")
        handle.write("set cartoon_fancy_helices, on\n")
        handle.write(f"orient {obj}\n")
        handle.write(f"zoom {obj}, complete=1\n\n")

        for klass in CLASSES:
            handle.write(f"count_atoms {sel[klass]} and name CA\n")

        handle.write("\n# Optional publication render:\n")
        handle.write("# ray 2400,1800\n")
        handle.write("# png R699_dual_context_mapping.png, dpi=300\n")

    print(f"Wrote: {args.output}")
    for klass in CLASSES:
        print(f"{klass}: {len(classes[klass])} residues")


if __name__ == "__main__":
    main()
