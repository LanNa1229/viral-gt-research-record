#!/usr/bin/env python3
"""Create a Foldseek batch manifest TSV from a list of PDB/CIF structure files.

Output columns:
query_id, protein_id, cluster, source_genome, model_source, structure_file
"""
import argparse
import csv
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Generate foldseek_manifest TSV from structure files")
    p.add_argument("--structure-dir", required=True, help="Directory containing .pdb/.cif/.mmcif structure files")
    p.add_argument("--out", required=True, help="Output manifest TSV path")
    p.add_argument("--cluster", default="unknown", help="Cluster name, e.g. cluster7 or cluster8")
    p.add_argument("--prefix", default="Q", help="query_id prefix, e.g. C7 or C8")
    p.add_argument("--source-genome", default="unknown", help="Value for source_genome column")
    p.add_argument("--model-source", default="AlphaFold", help="Value for model_source column, e.g. AlphaFoldDB, AF3, ColabFold")
    p.add_argument("--relative-to", default=None, help="Make structure_file paths relative to this directory; default: current working directory")
    return p.parse_args()


def main():
    args = parse_args()
    structure_dir = Path(args.structure_dir).expanduser().resolve()
    if not structure_dir.exists():
        raise SystemExit(f"ERROR: structure directory not found: {structure_dir}")

    files = []
    for pattern in ("*.pdb", "*.cif", "*.mmcif"):
        files.extend(structure_dir.glob(pattern))
    files = sorted(files)
    if not files:
        raise SystemExit(f"ERROR: no .pdb/.cif/.mmcif files found in {structure_dir}")

    base = Path(args.relative_to).expanduser().resolve() if args.relative_to else Path.cwd().resolve()
    out = Path(args.out)

    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["query_id", "protein_id", "cluster", "source_genome", "model_source", "structure_file"],
            delimiter="\t",
        )
        writer.writeheader()
        for i, f in enumerate(files, 1):
            protein_id = f.stem
            try:
                structure_file = f.resolve().relative_to(base).as_posix()
            except ValueError:
                structure_file = f.resolve().as_posix()
            writer.writerow({
                "query_id": f"{args.prefix}_{i:03d}",
                "protein_id": protein_id,
                "cluster": args.cluster,
                "source_genome": args.source_genome,
                "model_source": args.model_source,
                "structure_file": structure_file,
            })

    print(f"Wrote {len(files)} rows to {out}")


if __name__ == "__main__":
    main()
