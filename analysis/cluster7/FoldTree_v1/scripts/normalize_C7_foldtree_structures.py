#!/usr/bin/env python3
"""
Normalize every structure in a FoldTree custom-structure run to a genuine
plain-text PDB file.

Designed for:
    /mnt/d/FoldTree/C7_gated143_foldtree_run

Observed input problem:
    119 files named *.pdb
      - 112 genuine PDBs
      - 7 mmCIF files mislabeled with a .pdb suffix
    24 genuine *.cif files

FoldTree's installed dl_structs.py scans/parses PDB input, so all 143 structures
must be normalized to actual PDB syntax before rerunning the workflow.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, is_aa


SUPPORTED_SUFFIXES = (
    ".pdb", ".cif", ".mmcif",
    ".pdb.gz", ".cif.gz", ".mmcif.gz",
)

CHAIN_IDS = list(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789"
)


def strip_structure_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in sorted(SUPPORTED_SUFFIXES, key=len, reverse=True):
        if lower.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


def is_structure_file(path: Path) -> bool:
    return path.is_file() and path.name.lower().endswith(SUPPORTED_SUFFIXES)


def read_head(path: Path, n: int = 20000) -> str:
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
            return handle.read(n)
    return path.read_text(encoding="utf-8", errors="ignore")[:n]


def materialize_uncompressed(path: Path, temp_dir: Path) -> Path:
    if not path.name.lower().endswith(".gz"):
        return path

    output_name = strip_structure_suffix(path.name)
    suffix = ".cif" if ".cif" in path.name.lower() or ".mmcif" in path.name.lower() else ".pdb"
    output = temp_dir / f"{output_name}{suffix}"

    with gzip.open(path, "rb") as src, output.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return output


def detect_format(path: Path) -> str:
    head = read_head(path).lstrip()
    if (
        head.startswith("data_")
        or "_atom_site." in head
        or "_entry.id" in head
        or "_struct_asym." in head
    ):
        return "mmCIF"
    return "PDB"


def normalize_chain_ids(structure) -> dict[str, str]:
    """
    PDB format supports one-character chain IDs. Rename chains only when needed,
    preserving unique one-character IDs wherever possible.
    """
    chains = []
    for model in structure:
        chains.extend(list(model.get_chains()))

    if len(chains) > len(CHAIN_IDS):
        raise ValueError(
            f"Structure has {len(chains)} chains; plain PDB supports at most "
            f"{len(CHAIN_IDS)} uniquely mapped chain IDs in this converter"
        )

    used = set()
    mapping = {}
    current_ids = [str(chain.id) for chain in chains]

    for old_id in current_ids:
        if len(old_id) == 1 and old_id not in used:
            new_id = old_id
        else:
            new_id = next(candidate for candidate in CHAIN_IDS if candidate not in used)
        mapping[old_id] = new_id
        used.add(new_id)

    for chain, old_id in zip(chains, current_ids):
        chain.id = mapping[old_id]

    return mapping


def parse_structure(path: Path, structure_id: str, fmt: str):
    parser = MMCIFParser(QUIET=True) if fmt == "mmCIF" else PDBParser(QUIET=True)
    return parser.get_structure(structure_id, str(path))


def count_standard_aa(structure) -> int:
    return sum(
        1
        for residue in structure.get_residues()
        if is_aa(residue, standard=True)
    )


def count_chains(structure) -> int:
    return sum(1 for _ in structure.get_chains())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--run-dir",
        default="/mnt/d/FoldTree/C7_gated143_foldtree_run",
    )
    parser.add_argument("--expected-count", type=int, default=143)
    parser.add_argument(
        "--activate",
        action="store_true",
        help="Replace run-dir/structs after all normalized PDBs pass QC",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir)
    source_dir = run_dir / "structs"
    normalized_dir = run_dir / "structs_normalized_pdb"
    temp_dir = run_dir / "00_metadata" / "normalization_tmp"
    report_path = run_dir / "00_metadata" / "structure_normalization_qc.tsv"
    backup_dir = run_dir / "structs_before_pdb_normalization"

    if not source_dir.is_dir():
        raise SystemExit(f"Missing source directory: {source_dir}")

    files = sorted(path for path in source_dir.iterdir() if is_structure_file(path))
    if len(files) != args.expected_count:
        raise SystemExit(
            f"Expected {args.expected_count} structures in {source_dir}; "
            f"observed {len(files)}"
        )

    ids = [strip_structure_suffix(path.name) for path in files]
    duplicates = sorted({identifier for identifier in ids if ids.count(identifier) > 1})
    if duplicates:
        raise SystemExit(
            "Canonical filename collisions detected before conversion: "
            + ", ".join(duplicates)
        )

    if normalized_dir.exists():
        shutil.rmtree(normalized_dir)
    normalized_dir.mkdir(parents=True)

    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True)

    rows = []
    failures = []

    for index, source in enumerate(files, start=1):
        identifier = strip_structure_suffix(source.name)
        detected_format = detect_format(source)
        materialized = materialize_uncompressed(source, temp_dir)
        destination = normalized_dir / f"{identifier}.pdb"

        try:
            structure = parse_structure(materialized, identifier, detected_format)
            input_aa = count_standard_aa(structure)
            input_chains = count_chains(structure)

            if input_aa == 0:
                raise ValueError("No standard amino-acid residues found")

            chain_mapping = normalize_chain_ids(structure)

            writer = PDBIO()
            writer.set_structure(structure)
            writer.save(str(destination))

            output_structure = PDBParser(QUIET=True).get_structure(
                identifier,
                str(destination),
            )
            output_aa = count_standard_aa(output_structure)
            output_chains = count_chains(output_structure)

            if output_aa == 0:
                raise ValueError("Normalized PDB contains no amino-acid residues")
            if output_aa != input_aa:
                raise ValueError(
                    f"Amino-acid residue count changed: input={input_aa}, "
                    f"output={output_aa}"
                )

            status = "converted_from_mmcif" if detected_format == "mmCIF" else "rewritten_pdb"
            row = {
                "identifier": identifier,
                "source_filename": source.name,
                "detected_input_format": detected_format,
                "status": status,
                "input_standard_aa": input_aa,
                "output_standard_aa": output_aa,
                "input_chains": input_chains,
                "output_chains": output_chains,
                "chain_id_mapping": ";".join(
                    f"{old}->{new}" for old, new in chain_mapping.items()
                ),
                "output_filename": destination.name,
                "error": "",
            }
        except Exception as exc:
            failures.append((source.name, type(exc).__name__, str(exc)))
            row = {
                "identifier": identifier,
                "source_filename": source.name,
                "detected_input_format": detected_format,
                "status": "FAILED",
                "input_standard_aa": "",
                "output_standard_aa": "",
                "input_chains": "",
                "output_chains": "",
                "chain_id_mapping": "",
                "output_filename": "",
                "error": f"{type(exc).__name__}: {exc}",
            }

        rows.append(row)
        print(
            f"[{index:3d}/{len(files)}] {source.name}: "
            f"{row['status']}",
            flush=True,
        )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    output_files = sorted(normalized_dir.glob("*.pdb"))
    good_rows = [row for row in rows if row["status"] != "FAILED"]

    print()
    print(f"Input structures:       {len(files)}")
    print(f"Normalized PDB outputs: {len(output_files)}")
    print(f"Successful conversions: {len(good_rows)}")
    print(f"Failures:               {len(failures)}")
    print(f"QC report:              {report_path}")

    if failures:
        print("\nFAILED FILES")
        for filename, error_type, error in failures:
            print(f"{filename}\t{error_type}\t{error}")
        raise SystemExit(
            "Normalization failed. Original structs/ directory was not modified."
        )

    if len(output_files) != args.expected_count:
        raise SystemExit(
            f"Expected {args.expected_count} normalized PDB files; "
            f"observed {len(output_files)}"
        )

    if args.activate:
        if backup_dir.exists():
            raise SystemExit(
                f"Backup already exists: {backup_dir}\n"
                "Refusing to overwrite it."
            )

        source_dir.rename(backup_dir)
        normalized_dir.rename(source_dir)

        print()
        print("ACTIVATED NORMALIZED INPUT")
        print(f"Current FoldTree input: {source_dir}")
        print(f"Original input backup:  {backup_dir}")

        for relative in (
            "sequence_dataset.csv",
            "sequences.fst",
            "finalset.csv",
            "allvall_1.csv",
        ):
            (run_dir / relative).unlink(missing_ok=True)

    else:
        print()
        print("Dry conversion completed; structs/ was not replaced.")
        print("After reviewing the report, activate with:")
        print(
            f"  python {Path(__file__).name} "
            f"--run-dir {run_dir} --expected-count {args.expected_count} --activate"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
