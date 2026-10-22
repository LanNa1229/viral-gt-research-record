#!/usr/bin/env python3
"""Map Workflow-1 reference-only within-community cores from C005 and C009
onto an out-of-sample R699 query through each community medoid.

The core is frozen from all_residue_metrics.tsv (default flag:
passes_within_quality). Query-inclusive FoldMason AA alignments are used only to
transfer medoid sequence positions to R699 sequence positions. The query never
contributes to the reference core definition.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from Bio import Align
from Bio.Data.IUPACData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Map C005 and C009 Workflow-1 within-cores onto R699."
    )
    p.add_argument("--metrics", required=True, type=Path,
                   help="Workflow-1 all_residue_metrics.tsv")
    p.add_argument("--medoids", required=True, type=Path,
                   help="Workflow-1 medoids.tsv")
    p.add_argument("--c005-reference-aa", required=True, type=Path,
                   help="Reference-only C005 FoldMason AA alignment")
    p.add_argument("--c005-query-aa", required=True, type=Path,
                   help="C005 references + R699 FoldMason AA alignment")
    p.add_argument("--c009-reference-aa", required=True, type=Path,
                   help="Reference-only C009 FoldMason AA alignment")
    p.add_argument("--c009-query-aa", required=True, type=Path,
                   help="C009 references + R699 FoldMason AA alignment")
    p.add_argument("--query-id", default="9DYT_A",
                   help="R699 FASTA identifier or normalized identifier")
    p.add_argument("--query-structure", required=True, type=Path,
                   help="R699 PDB/mmCIF structure used in FoldMason")
    p.add_argument("--query-chain", default="A",
                   help="R699 structure chain; default A")
    p.add_argument("--core-flag", default="passes_within_quality",
                   help="Boolean column defining the frozen reference core")
    p.add_argument("--min-within-lddt", type=float, default=None,
                   help="Optional additional within_lddt threshold")
    p.add_argument("--min-within-presence", type=float, default=None,
                   help="Optional additional within_presence threshold")
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def normalize_id(text: str) -> str:
    x = Path(str(text).strip().split()[0]).name
    for suffix in (".pdb.gz", ".cif.gz", ".mmcif.gz", ".pdb", ".cif", ".mmcif"):
        if x.lower().endswith(suffix):
            x = x[:-len(suffix)]
            break
    if "__" in x:
        x = x.split("__", 1)[0]
    return x


def read_fasta(path: Path) -> Dict[str, str]:
    records: Dict[str, List[str]] = {}
    current: Optional[str] = None
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate FASTA ID {current!r} in {path}")
                records[current] = []
            else:
                if current is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                records[current].append(line)
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    return {k: "".join(v) for k, v in records.items()}


def resolve_record(records: Dict[str, str], wanted: str, label: str) -> Tuple[str, str]:
    if wanted in records:
        return wanted, records[wanted]
    nw = normalize_id(wanted).lower()
    hits = [(k, v) for k, v in records.items() if normalize_id(k).lower() == nw]
    if len(hits) == 1:
        return hits[0]
    if not hits:
        available = ", ".join(list(records)[:10])
        raise KeyError(f"Could not find {label} {wanted!r}; first IDs: {available}")
    raise KeyError(f"Ambiguous {label} {wanted!r}: {[k for k, _ in hits]}")


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(
        {"true", "1", "yes", "y", "t"}
    )


def pick_column(df: pd.DataFrame, candidates: Sequence[str], description: str) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"Could not find {description}; tried {list(candidates)}")


def alignment_index_maps(aligned_seq: str) -> Tuple[Dict[int, int], Dict[int, int]]:
    """Return seq_index->column and column->seq_index, all 1-based."""
    seq_to_col: Dict[int, int] = {}
    col_to_seq: Dict[int, int] = {}
    seq_i = 0
    for col_i, aa in enumerate(aligned_seq, start=1):
        if aa not in {"-", "."}:
            seq_i += 1
            seq_to_col[seq_i] = col_i
            col_to_seq[col_i] = seq_i
    return seq_to_col, col_to_seq


def structure_residues(path: Path, chain_id: str) -> Tuple[str, List[Tuple[str, int, str, str]]]:
    lower = path.name.lower()
    parser = MMCIFParser(QUIET=True) if lower.endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
    structure = parser.get_structure("query", str(path))
    model = next(structure.get_models())
    if chain_id not in model:
        chains = [c.id for c in model]
        if len(chains) == 1:
            chain = model[chains[0]]
            print(f"WARNING: chain {chain_id!r} not found; using sole chain {chains[0]!r}", file=sys.stderr)
        else:
            raise KeyError(f"Chain {chain_id!r} not found in {path}; chains={chains}")
    else:
        chain = model[chain_id]

    seq: List[str] = []
    residues: List[Tuple[str, int, str, str]] = []
    table = {k.upper(): v.upper() for k, v in protein_letters_3to1_extended.items()}
    for res in chain:
        if not is_aa(res, standard=False):
            continue
        resname = res.get_resname().strip().upper()
        aa = table.get(resname, "X")
        hetflag, resseq, icode = res.id
        seq.append(aa)
        residues.append((chain.id, int(resseq), str(icode).strip(), resname))
    if not residues:
        raise ValueError(f"No amino-acid residues found in {path}, chain {chain_id}")
    return "".join(seq), residues


def map_query_sequence_to_structure(query_seq: str, structure_seq: str) -> Dict[int, int]:
    """Map 1-based query sequence indices to 1-based structure sequence indices."""
    if query_seq == structure_seq:
        return {i: i for i in range(1, len(query_seq) + 1)}

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    aln = aligner.align(query_seq, structure_seq)[0]
    mapping: Dict[int, int] = {}
    q_blocks, s_blocks = aln.aligned
    for (q0, q1), (s0, s1) in zip(q_blocks, s_blocks):
        length = min(q1 - q0, s1 - s0)
        for offset in range(length):
            mapping[q0 + offset + 1] = s0 + offset + 1
    coverage = len(mapping) / max(1, len(query_seq))
    if coverage < 0.90:
        raise ValueError(
            f"Only {coverage:.1%} of FoldMason query sequence mapped to structure sequence; "
            "check chain and structure identity."
        )
    print(
        f"WARNING: FoldMason query sequence was not identical to structure sequence; "
        f"global alignment mapped {len(mapping)}/{len(query_seq)} residues ({coverage:.1%}).",
        file=sys.stderr,
    )
    return mapping


def medoid_for_context(medoids: pd.DataFrame, context: str) -> str:
    cluster_col = pick_column(medoids, ["cluster_id", "community", "cluster"], "cluster column")
    medoid_col = pick_column(medoids, ["medoid_id", "protein_id", "medoid"], "medoid ID column")
    hit = medoids[medoids[cluster_col].astype(str) == context]
    if len(hit) != 1:
        raise ValueError(f"Expected one medoid row for {context}, found {len(hit)}")
    return str(hit.iloc[0][medoid_col])


def map_context(
    context: str,
    metrics: pd.DataFrame,
    medoid_id: str,
    reference_aa_path: Path,
    query_aa_path: Path,
    query_id: str,
    core_flag: str,
    min_within_lddt: Optional[float],
    min_within_presence: Optional[float],
) -> Tuple[pd.DataFrame, str]:
    cluster_col = pick_column(metrics, ["cluster_id", "community", "cluster"], "metrics cluster column")
    index_col = pick_column(
        metrics,
        ["residue_index", "sequence_index", "seq_index", "medoid_residue_index"],
        "medoid residue sequence-index column",
    )
    if core_flag not in metrics.columns:
        raise KeyError(f"Core flag column {core_flag!r} not present in metrics")

    group = metrics[metrics[cluster_col].astype(str) == context].copy()
    if group.empty:
        raise ValueError(f"No metrics rows for {context}")
    group["_core"] = bool_series(group[core_flag])
    if min_within_lddt is not None:
        if "within_lddt" not in group.columns:
            raise KeyError("--min-within-lddt requested but within_lddt column is absent")
        group["_core"] &= pd.to_numeric(group["within_lddt"], errors="coerce") >= min_within_lddt
    if min_within_presence is not None:
        if "within_presence" not in group.columns:
            raise KeyError("--min-within-presence requested but within_presence column is absent")
        group["_core"] &= pd.to_numeric(group["within_presence"], errors="coerce") >= min_within_presence
    core = group[group["_core"]].copy()

    ref_records = read_fasta(reference_aa_path)
    query_records = read_fasta(query_aa_path)
    _, ref_medoid_aln = resolve_record(ref_records, medoid_id, f"{context} medoid in reference MSA")
    q_medoid_header, q_medoid_aln = resolve_record(query_records, medoid_id, f"{context} medoid in query MSA")
    q_query_header, q_query_aln = resolve_record(query_records, query_id, "R699 query")

    if len(q_medoid_aln) != len(q_query_aln):
        raise ValueError(f"Aligned sequences differ in length in {query_aa_path}")
    ref_medoid_ungapped = ref_medoid_aln.replace("-", "").replace(".", "")
    q_medoid_ungapped = q_medoid_aln.replace("-", "").replace(".", "")
    if ref_medoid_ungapped != q_medoid_ungapped:
        raise ValueError(
            f"{context} medoid sequence differs between reference-only and query-inclusive MSAs; "
            "do not transfer columns directly."
        )

    medoid_len = len(q_medoid_ungapped)
    raw_indices = pd.to_numeric(core[index_col], errors="raise").astype(int)
    if len(raw_indices) and raw_indices.min() == 0 and raw_indices.max() <= medoid_len - 1:
        core["medoid_sequence_index"] = raw_indices + 1
    else:
        core["medoid_sequence_index"] = raw_indices
    if len(core) and (core["medoid_sequence_index"].min() < 1 or core["medoid_sequence_index"].max() > medoid_len):
        raise ValueError(
            f"{context} core residue indices fall outside medoid length {medoid_len}; "
            f"range={core['medoid_sequence_index'].min()}..{core['medoid_sequence_index'].max()}"
        )

    medoid_seq_to_col, _ = alignment_index_maps(q_medoid_aln)
    _, query_col_to_seq = alignment_index_maps(q_query_aln)
    query_ungapped = q_query_aln.replace("-", "").replace(".", "")

    out_rows: List[dict] = []
    metric_cols = [
        c for c in [
            "within_presence", "within_lddt", "within_lddt_sd",
            "within_lddt_coverage", "within_gap_frequency",
            "aa_conservation", "within_3Di_conservation", "local_pLDDT",
        ] if c in core.columns
    ]
    for _, row in core.iterrows():
        m_idx = int(row["medoid_sequence_index"])
        col = medoid_seq_to_col.get(m_idx)
        q_idx = query_col_to_seq.get(col) if col is not None else None
        record = {
            "context": context,
            "medoid_id": normalize_id(medoid_id),
            "medoid_sequence_index": m_idx,
            "query_alignment_column": col if col is not None else "",
            "query_sequence_index": q_idx if q_idx is not None else "",
            "query_aa": query_ungapped[q_idx - 1] if q_idx is not None else "",
            "mapped_to_query": q_idx is not None,
            "query_header": q_query_header,
            "medoid_header": q_medoid_header,
        }
        for c in metric_cols:
            record[c] = row[c]
        out_rows.append(record)

    out = pd.DataFrame(out_rows)
    return out, query_ungapped


def residue_selector(rows: pd.DataFrame) -> str:
    nums = []
    for _, r in rows.iterrows():
        if pd.isna(r.get("structure_resseq")):
            continue
        token = str(int(r["structure_resseq"]))
        icode = str(r.get("structure_icode", "") or "").strip()
        if icode:
            token += icode
        nums.append(token)
    return "+".join(dict.fromkeys(nums)) if nums else "none"


def main() -> None:
    args = parse_args()
    for path in [
        args.metrics, args.medoids, args.c005_reference_aa, args.c005_query_aa,
        args.c009_reference_aa, args.c009_query_aa, args.query_structure,
    ]:
        if not path.exists():
            raise SystemExit(f"ERROR: missing file: {path}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(args.metrics, sep="\t")
    medoids = pd.read_csv(args.medoids, sep="\t")

    c005_medoid = medoid_for_context(medoids, "C005")
    c009_medoid = medoid_for_context(medoids, "C009")

    c005, qseq005 = map_context(
        "C005", metrics, c005_medoid, args.c005_reference_aa,
        args.c005_query_aa, args.query_id, args.core_flag,
        args.min_within_lddt, args.min_within_presence,
    )
    c009, qseq009 = map_context(
        "C009", metrics, c009_medoid, args.c009_reference_aa,
        args.c009_query_aa, args.query_id, args.core_flag,
        args.min_within_lddt, args.min_within_presence,
    )
    if qseq005 != qseq009:
        raise SystemExit("ERROR: R699 ungapped sequence differs between C005 and C009 query-inclusive MSAs")
    query_seq = qseq005

    structure_seq, structure_rows = structure_residues(args.query_structure, args.query_chain)
    q_to_s = map_query_sequence_to_structure(query_seq, structure_seq)

    def add_structure_coordinates(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        chain_vals, resseq_vals, icode_vals, resname_vals = [], [], [], []
        for value in df["query_sequence_index"]:
            if value == "" or pd.isna(value):
                chain_vals.append(""); resseq_vals.append(""); icode_vals.append(""); resname_vals.append("")
                continue
            q_idx = int(value)
            s_idx = q_to_s.get(q_idx)
            if s_idx is None:
                chain_vals.append(""); resseq_vals.append(""); icode_vals.append(""); resname_vals.append("")
                continue
            chain, resseq, icode, resname = structure_rows[s_idx - 1]
            chain_vals.append(chain); resseq_vals.append(resseq); icode_vals.append(icode); resname_vals.append(resname)
        df["structure_chain"] = chain_vals
        df["structure_resseq"] = resseq_vals
        df["structure_icode"] = icode_vals
        df["structure_resname"] = resname_vals
        return df

    c005 = add_structure_coordinates(c005)
    c009 = add_structure_coordinates(c009)
    c005.to_csv(args.outdir / "C005_R699_core_mapping.tsv", sep="\t", index=False)
    c009.to_csv(args.outdir / "C009_R699_core_mapping.tsv", sep="\t", index=False)

    c005_set = set(pd.to_numeric(c005.loc[c005["mapped_to_query"], "query_sequence_index"]).astype(int))
    c009_set = set(pd.to_numeric(c009.loc[c009["mapped_to_query"], "query_sequence_index"]).astype(int))

    rows = []
    for q_idx, aa in enumerate(query_seq, start=1):
        in5 = q_idx in c005_set
        in9 = q_idx in c009_set
        if in5 and in9:
            klass = "both_supported"
        elif in5:
            klass = "C005_supported"
        elif in9:
            klass = "C009_supported"
        else:
            klass = "unsupported"
        s_idx = q_to_s.get(q_idx)
        if s_idx is not None:
            chain, resseq, icode, resname = structure_rows[s_idx - 1]
        else:
            chain, resseq, icode, resname = "", "", "", ""
        rows.append({
            "query_id": normalize_id(args.query_id),
            "query_sequence_index": q_idx,
            "query_aa": aa,
            "structure_chain": chain,
            "structure_resseq": resseq,
            "structure_icode": icode,
            "structure_resname": resname,
            "C005_supported": in5,
            "C009_supported": in9,
            "constraint_class": klass,
        })
    dual = pd.DataFrame(rows)
    dual.to_csv(args.outdir / "R699_dual_context_classification.tsv", sep="\t", index=False)

    summary = (
        dual.groupby("constraint_class", sort=False)
        .size().rename("n_residues").reset_index()
    )
    summary["fraction_of_R699"] = summary["n_residues"] / len(dual)
    summary.to_csv(args.outdir / "R699_dual_context_summary.tsv", sep="\t", index=False)

    object_name = "R699"
    pml = args.outdir / "R699_dual_context_mapping.pml"
    selections = {}
    for klass in ["both_supported", "C005_supported", "C009_supported", "unsupported"]:
        selections[klass] = residue_selector(dual[dual["constraint_class"] == klass])
    with pml.open("w") as handle:
        handle.write("reinitialize\n")
        handle.write(f"load {args.query_structure.resolve()}, {object_name}\n")
        handle.write("hide everything\nshow cartoon, R699\ncolor gray80, R699\n")
        for klass, selector in selections.items():
            handle.write(f"select {klass}, R699 and chain {args.query_chain} and resi {selector}\n")
        handle.write("color purple, both_supported\n")
        handle.write("color tv_orange, C005_supported\n")
        handle.write("color marine, C009_supported\n")
        handle.write("color gray80, unsupported\n")
        handle.write("bg_color white\nset orthoscopic, on\nset ray_opaque_background, off\n")

    print(f"C005 medoid: {c005_medoid}")
    print(f"C009 medoid: {c009_medoid}")
    print(f"C005 frozen-core residues: {len(c005)}; mapped to R699: {int(c005['mapped_to_query'].sum())}")
    print(f"C009 frozen-core residues: {len(c009)}; mapped to R699: {int(c009['mapped_to_query'].sum())}")
    print("\nR699 classification:")
    print(summary.to_string(index=False))
    print(f"\nOutputs: {args.outdir}")


if __name__ == "__main__":
    main()
