#!/usr/bin/env python3
"""Singleton length/pLDDT QC workflow for Foldseek coarse-clustering outputs.

This script reconstructs cluster membership, parses FASTA and predicted structures,
computes sequence/structure-confidence metrics, derives exploratory QC categories,
performs non-parametric statistics, and produces publication-oriented plots.

It never modifies raw input files. All outputs are written beneath --outdir.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import sys
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

try:
    from Bio.PDB import MMCIFParser
except Exception:
    MMCIFParser = None


RNG_SEED = 20260712
GROUP_ORDER = ["giant_cluster", "other_multimember", "singleton"]
PAIRWISE_COMPARISONS = [
    ("singleton", "all_multimember"),
    ("singleton", "giant_cluster"),
    ("singleton", "other_multimember"),
    ("giant_cluster", "other_multimember"),
]
METRICS_FOR_STATS = [
    "fasta_length",
    "median_plddt",
    "mean_plddt",
    "frac_plddt_lt50",
    "frac_plddt_lt70",
    "longest_segment_lt70",
    "structure_coverage_fraction",
]


def log(msg: str, log_handle=None) -> None:
    print(msg, flush=True)
    if log_handle is not None:
        log_handle.write(msg + "\n")
        log_handle.flush()


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return path.open("rt", errors="replace")


def normalize_id(value: str) -> str:
    value = str(value).strip()
    value = value.split()[0] if value else value
    return value


def detect_delimiter(path: Path) -> str:
    lines = []
    with open_text(path) as fh:
        for line in fh:
            line = line.strip("\n\r")
            if not line or line.startswith("#"):
                continue
            lines.append(line)
            if len(lines) >= 10:
                break
    if not lines:
        raise ValueError(f"Cluster file is empty: {path}")
    sample = "\n".join(lines)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,;").delimiter
    except csv.Error:
        return "\t"


def load_cluster_membership(
    path: Path,
    rep_column: Optional[str] = None,
    member_column: Optional[str] = None,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    delimiter = detect_delimiter(path)
    raw_no_header = pd.read_csv(
        path,
        sep=delimiter,
        comment="#",
        header=None,
        dtype=str,
        keep_default_na=False,
    )
    if raw_no_header.shape[1] < 2:
        raise ValueError("Cluster membership file must contain at least two columns.")

    first = [str(x).strip().lower() for x in raw_no_header.iloc[0].tolist()]
    header_tokens = {
        "representative", "representative_id", "rep", "rep_id", "cluster_representative",
        "member", "member_id", "protein", "protein_id", "query", "target"
    }
    header_likely = any(x in header_tokens for x in first)

    if rep_column or member_column or header_likely:
        df = pd.read_csv(
            path,
            sep=delimiter,
            comment="#",
            header=0,
            dtype=str,
            keep_default_na=False,
        )
        df.columns = [str(c).strip() for c in df.columns]
        lower_map = {c.lower(): c for c in df.columns}
        if rep_column is None:
            candidates = [
                "representative_id", "representative", "rep_id", "rep",
                "cluster_representative", "cluster"
            ]
            rep_column = next((lower_map[x] for x in candidates if x in lower_map), None)
        if member_column is None:
            candidates = ["member_id", "member", "protein_id", "protein", "target"]
            member_column = next((lower_map[x] for x in candidates if x in lower_map), None)
        if rep_column is None or member_column is None:
            if df.shape[1] == 2:
                rep_column, member_column = df.columns[:2]
            else:
                raise ValueError(
                    "Could not infer representative/member columns. Pass --rep-column and --member-column."
                )
        cluster = df[[rep_column, member_column]].copy()
    else:
        cluster = raw_no_header.iloc[:, :2].copy()
        cluster.columns = ["representative_id", "member_id"]
        rep_column, member_column = cluster.columns

    cluster.columns = ["representative_id", "member_id"]
    cluster["representative_id"] = cluster["representative_id"].map(normalize_id)
    cluster["member_id"] = cluster["member_id"].map(normalize_id)
    empty_rows = cluster[(cluster["representative_id"] == "") | (cluster["member_id"] == "")]
    cluster = cluster[(cluster["representative_id"] != "") & (cluster["member_id"] != "")].copy()

    exact_duplicate_rows = int(cluster.duplicated(["representative_id", "member_id"]).sum())
    member_rep_counts = cluster.groupby("member_id")["representative_id"].nunique()
    conflicting_members = member_rep_counts[member_rep_counts > 1].index.tolist()
    cluster = cluster.drop_duplicates(["representative_id", "member_id"]).copy()

    audit = {
        "delimiter": repr(delimiter),
        "header_detected": bool(header_likely),
        "representative_column": str(rep_column),
        "member_column": str(member_column),
        "empty_rows_removed": int(len(empty_rows)),
        "exact_duplicate_rows_removed": exact_duplicate_rows,
        "members_assigned_to_multiple_representatives": len(conflicting_members),
        "conflicting_member_ids": conflicting_members,
    }
    return cluster, audit


def parse_fasta(path: Optional[Path]) -> Tuple[Dict[str, str], List[str], Dict[str, str]]:
    if path is None:
        return {}, [], {}
    seqs: Dict[str, str] = {}
    headers: Dict[str, str] = {}
    duplicates: List[str] = []
    current_id: Optional[str] = None
    chunks: List[str] = []

    def commit() -> None:
        nonlocal current_id, chunks
        if current_id is None:
            return
        seq = re.sub(r"[^A-Za-z*]", "", "".join(chunks)).replace("*", "").upper()
        if current_id in seqs:
            duplicates.append(current_id)
        else:
            seqs[current_id] = seq
        chunks = []

    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                commit()
                header = line[1:].strip()
                current_id = normalize_id(header)
                headers[current_id] = header
            else:
                if current_id is None:
                    raise ValueError(f"Sequence encountered before FASTA header in {path}")
                chunks.append(line)
        commit()
    return seqs, duplicates, headers


def strip_structure_suffix(filename: str) -> str:
    name = filename
    for suffix in [".pdb.gz", ".cif.gz", ".mmcif.gz", ".pdb", ".cif", ".mmcif"]:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def discover_structures(structures_dir: Optional[Path]) -> List[Path]:
    if structures_dir is None:
        return []
    valid = (".pdb", ".cif", ".mmcif", ".pdb.gz", ".cif.gz", ".mmcif.gz")
    return sorted(
        p for p in structures_dir.rglob("*")
        if p.is_file() and any(str(p).lower().endswith(ext) for ext in valid)
    )


def load_structure_id_map(path: Optional[Path]) -> Dict[str, Path]:
    if path is None:
        return {}
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    lower = {c.lower(): c for c in df.columns}
    id_col = next((lower[x] for x in ["protein_id", "id", "member_id"] if x in lower), None)
    file_col = next((lower[x] for x in ["structure_file", "path", "filename"] if x in lower), None)
    if id_col is None or file_col is None:
        raise ValueError("Structure ID map must have protein_id/id and structure_file/path columns.")
    return {normalize_id(r[id_col]): Path(r[file_col]).expanduser() for _, r in df.iterrows()}


def match_structure_files(
    protein_ids: Sequence[str],
    structure_files: Sequence[Path],
    explicit_map: Dict[str, Path],
) -> Tuple[Dict[str, Optional[Path]], pd.DataFrame]:
    by_stem: Dict[str, List[Path]] = defaultdict(list)
    for p in structure_files:
        stem = strip_structure_suffix(p.name)
        by_stem[stem].append(p)

    matches: Dict[str, Optional[Path]] = {}
    rows = []
    for pid in protein_ids:
        candidates: List[Path] = []
        method = "none"
        if pid in explicit_map:
            candidates = [explicit_map[pid]]
            method = "explicit_map"
        elif pid in by_stem:
            candidates = by_stem[pid]
            method = "exact_stem"
        else:
            pattern = re.compile(rf"^{re.escape(pid)}(?:[._-]|$)")
            candidates = [p for stem, paths in by_stem.items() if pattern.search(stem) for p in paths]
            if candidates:
                method = "prefix_stem"
        candidates = sorted(set(candidates))
        chosen = candidates[0] if len(candidates) == 1 else None
        matches[pid] = chosen
        rows.append({
            "protein_id": pid,
            "match_method": method,
            "candidate_count": len(candidates),
            "matched_structure": str(chosen) if chosen else "",
            "candidate_structures": ";".join(map(str, candidates)),
            "status": "matched" if chosen else ("ambiguous" if candidates else "missing"),
        })
    return matches, pd.DataFrame(rows)


@dataclass
class StructureMetrics:
    structure_file: str = ""
    structure_format: str = ""
    modeled_residue_count: float = np.nan
    chain_count: float = np.nan
    chain_residue_counts: str = ""
    mean_plddt: float = np.nan
    median_plddt: float = np.nan
    sd_plddt: float = np.nan
    min_plddt: float = np.nan
    max_plddt: float = np.nan
    frac_plddt_lt50: float = np.nan
    frac_plddt_lt70: float = np.nan
    frac_plddt_ge70: float = np.nan
    frac_plddt_ge90: float = np.nan
    longest_segment_lt50: float = np.nan
    longest_segment_lt70: float = np.nan
    n_terminal_low_confidence_length: float = np.nan
    c_terminal_low_confidence_length: float = np.nan
    high_confidence_segment_count_ge20: float = np.nan
    longest_high_confidence_segment: float = np.nan
    plddt_valid_0_100: object = np.nan
    parsing_status: str = "not_parsed"
    notes: str = ""
    residue_profile: Optional[List[Tuple[str, int, float]]] = None


def longest_true_run(mask: Sequence[bool]) -> int:
    best = current = 0
    for value in mask:
        if bool(value):
            current += 1
            best = max(best, current)
        else:
            current = 0
    return best


def terminal_true_run(mask: Sequence[bool], from_end: bool = False) -> int:
    iterable = reversed(mask) if from_end else iter(mask)
    n = 0
    for value in iterable:
        if bool(value):
            n += 1
        else:
            break
    return n


def true_segments(mask: Sequence[bool]) -> List[int]:
    lengths = []
    current = 0
    for value in mask:
        if bool(value):
            current += 1
        elif current:
            lengths.append(current)
            current = 0
    if current:
        lengths.append(current)
    return lengths


def summarize_residue_profile(
    chain_values: Dict[str, List[Tuple[int, float]]],
    structure_file: Path,
    structure_format: str,
) -> StructureMetrics:
    all_values = [v for vals in chain_values.values() for _, v in vals]
    count = len(all_values)
    if count == 0:
        return StructureMetrics(
            structure_file=str(structure_file),
            structure_format=structure_format,
            modeled_residue_count=0,
            chain_count=0,
            parsing_status="no_ca_residues",
            notes="No CA atoms with residue-level confidence values found.",
            residue_profile=[],
        )
    arr = np.asarray(all_values, dtype=float)
    per_chain_lt50 = []
    per_chain_lt70 = []
    per_chain_ge70_segments = []
    n_terminal = []
    c_terminal = []
    profile: List[Tuple[str, int, float]] = []
    for chain, vals in chain_values.items():
        ordered = sorted(vals, key=lambda x: x[0])
        scores = np.asarray([x[1] for x in ordered], dtype=float)
        mask50 = scores < 50
        mask70 = scores < 70
        mask_ge70 = scores >= 70
        per_chain_lt50.append(longest_true_run(mask50))
        per_chain_lt70.append(longest_true_run(mask70))
        per_chain_ge70_segments.extend(true_segments(mask_ge70))
        n_terminal.append(terminal_true_run(mask70, from_end=False))
        c_terminal.append(terminal_true_run(mask70, from_end=True))
        profile.extend((chain, i + 1, float(score)) for i, (_, score) in enumerate(ordered))

    valid = bool(np.nanmin(arr) >= 0 and np.nanmax(arr) <= 100)
    notes = []
    if not valid:
        notes.append("Confidence values fall outside 0-100; B-factor may not represent pLDDT.")
    if np.nanstd(arr, ddof=1) < 1e-6 and len(arr) > 1:
        notes.append("Nearly constant confidence values; verify structure provenance.")
    if len(chain_values) > 1:
        notes.append("Multiple chains detected.")

    return StructureMetrics(
        structure_file=str(structure_file),
        structure_format=structure_format,
        modeled_residue_count=count,
        chain_count=len(chain_values),
        chain_residue_counts=";".join(f"{k}:{len(v)}" for k, v in chain_values.items()),
        mean_plddt=float(np.mean(arr)),
        median_plddt=float(np.median(arr)),
        sd_plddt=float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
        min_plddt=float(np.min(arr)),
        max_plddt=float(np.max(arr)),
        frac_plddt_lt50=float(np.mean(arr < 50)),
        frac_plddt_lt70=float(np.mean(arr < 70)),
        frac_plddt_ge70=float(np.mean(arr >= 70)),
        frac_plddt_ge90=float(np.mean(arr >= 90)),
        longest_segment_lt50=max(per_chain_lt50, default=0),
        longest_segment_lt70=max(per_chain_lt70, default=0),
        n_terminal_low_confidence_length=max(n_terminal, default=0),
        c_terminal_low_confidence_length=max(c_terminal, default=0),
        high_confidence_segment_count_ge20=sum(x >= 20 for x in per_chain_ge70_segments),
        longest_high_confidence_segment=max(per_chain_ge70_segments, default=0),
        plddt_valid_0_100=valid,
        parsing_status="ok" if valid else "suspicious_plddt_range",
        notes=" ".join(notes),
        residue_profile=profile,
    )


def parse_pdb(path: Path) -> StructureMetrics:
    chain_values: Dict[str, List[Tuple[int, float]]] = defaultdict(list)
    seen = set()
    in_first_model = True
    model_seen = False
    with open_text(path) as fh:
        order = 0
        for line in fh:
            rec = line[:6].strip()
            if rec == "MODEL":
                if model_seen:
                    in_first_model = False
                else:
                    model_seen = True
                    in_first_model = True
                continue
            if rec == "ENDMDL" and model_seen:
                break
            if not in_first_model or rec != "ATOM":
                continue
            atom_name = line[12:16].strip()
            altloc = line[16:17].strip()
            if atom_name != "CA" or altloc not in ("", "A", "1"):
                continue
            chain = line[21:22].strip() or "_"
            resseq = line[22:26].strip()
            icode = line[26:27].strip()
            key = (chain, resseq, icode)
            if key in seen:
                continue
            seen.add(key)
            try:
                b = float(line[60:66].strip())
            except ValueError:
                continue
            order += 1
            chain_values[chain].append((order, b))
    return summarize_residue_profile(chain_values, path, "PDB")


def parse_mmcif(path: Path) -> StructureMetrics:
    if MMCIFParser is None:
        raise RuntimeError("Biopython is required for mmCIF parsing.")
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", str(path))
    model = next(structure.get_models())
    chain_values: Dict[str, List[Tuple[int, float]]] = defaultdict(list)
    order = 0
    for chain in model:
        chain_id = str(chain.id).strip() or "_"
        for residue in chain:
            if "CA" not in residue:
                continue
            atom = residue["CA"]
            if atom.is_disordered():
                atom = atom.selected_child
            order += 1
            chain_values[chain_id].append((order, float(atom.get_bfactor())))
    return summarize_residue_profile(chain_values, path, "mmCIF")


def parse_structure(path: Optional[Path]) -> StructureMetrics:
    if path is None:
        return StructureMetrics(parsing_status="missing_structure", notes="No matched structure file.")
    try:
        lower = str(path).lower()
        if lower.endswith((".pdb", ".pdb.gz")):
            return parse_pdb(path)
        if lower.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
            if lower.endswith(".gz"):
                raise RuntimeError("Compressed mmCIF is not supported directly; decompress or map to an uncompressed file.")
            return parse_mmcif(path)
        return StructureMetrics(
            structure_file=str(path), parsing_status="unsupported_format", notes="Unsupported structure extension."
        )
    except Exception as exc:
        return StructureMetrics(
            structure_file=str(path), parsing_status="parse_failure",
            notes=f"{type(exc).__name__}: {exc}"
        )


def mad(values: pd.Series) -> float:
    arr = pd.to_numeric(values, errors="coerce").dropna().to_numpy(float)
    if arr.size == 0:
        return np.nan
    med = np.median(arr)
    return float(np.median(np.abs(arr - med)))


def robust_reference_thresholds(master: pd.DataFrame) -> Dict[str, float]:
    ref = master[master["cluster_group"].isin(["giant_cluster", "other_multimember"])].copy()
    lengths = pd.to_numeric(ref["fasta_length"], errors="coerce").dropna()
    if lengths.empty:
        lengths = pd.to_numeric(ref["modeled_residue_count"], errors="coerce").dropna()
    q1, med, q3 = lengths.quantile([0.25, 0.5, 0.75]).tolist() if len(lengths) else (np.nan,) * 3
    iqr = q3 - q1 if np.isfinite(q1) and np.isfinite(q3) else np.nan
    length_mad = mad(lengths)
    length_lower_iqr = max(1.0, q1 - 1.5 * iqr) if np.isfinite(iqr) else np.nan
    length_upper_iqr = q3 + 1.5 * iqr if np.isfinite(iqr) else np.nan
    if np.isfinite(length_mad) and length_mad > 0:
        length_lower_mad = max(1.0, med - 3.5 * 1.4826 * length_mad)
        length_upper_mad = med + 3.5 * 1.4826 * length_mad
    else:
        length_lower_mad = length_upper_mad = np.nan

    plddt = pd.to_numeric(ref["median_plddt"], errors="coerce").dropna()
    lowfrac = pd.to_numeric(ref["frac_plddt_lt70"], errors="coerce").dropna()
    terminal = pd.concat([
        pd.to_numeric(ref["n_terminal_low_confidence_length"], errors="coerce"),
        pd.to_numeric(ref["c_terminal_low_confidence_length"], errors="coerce"),
    ]).dropna()
    coverage = pd.to_numeric(ref["structure_coverage_fraction"], errors="coerce").dropna()

    return {
        "reference_population_n": int(len(ref)),
        "length_median": float(med) if np.isfinite(med) else np.nan,
        "length_q1": float(q1) if np.isfinite(q1) else np.nan,
        "length_q3": float(q3) if np.isfinite(q3) else np.nan,
        "length_iqr": float(iqr) if np.isfinite(iqr) else np.nan,
        "length_mad": float(length_mad) if np.isfinite(length_mad) else np.nan,
        "length_lower_iqr_fence": float(length_lower_iqr) if np.isfinite(length_lower_iqr) else np.nan,
        "length_upper_iqr_fence": float(length_upper_iqr) if np.isfinite(length_upper_iqr) else np.nan,
        "length_lower_mad_fence": float(length_lower_mad) if np.isfinite(length_lower_mad) else np.nan,
        "length_upper_mad_fence": float(length_upper_mad) if np.isfinite(length_upper_mad) else np.nan,
        "median_plddt_05pct": float(plddt.quantile(0.05)) if len(plddt) else np.nan,
        "median_plddt_reference_median": float(plddt.median()) if len(plddt) else np.nan,
        "frac_lt70_95pct": float(lowfrac.quantile(0.95)) if len(lowfrac) else np.nan,
        "terminal_lt70_95pct": float(terminal.quantile(0.95)) if len(terminal) else np.nan,
        "coverage_05pct": float(coverage.quantile(0.05)) if len(coverage) else np.nan,
    }


def choose_length_fences(th: Dict[str, float]) -> Tuple[float, float]:
    lowers = [th.get("length_lower_iqr_fence"), th.get("length_lower_mad_fence")]
    uppers = [th.get("length_upper_iqr_fence"), th.get("length_upper_mad_fence")]
    lowers = [x for x in lowers if x is not None and np.isfinite(x)]
    uppers = [x for x in uppers if x is not None and np.isfinite(x)]
    return (max(lowers) if lowers else np.nan, min(uppers) if uppers else np.nan)


def classify_qc(master: pd.DataFrame, thresholds: Dict[str, float]) -> pd.DataFrame:
    df = master.copy()
    lower_len, upper_len = choose_length_fences(thresholds)
    low_global_threshold = thresholds.get("median_plddt_05pct", np.nan)
    high_lowfrac_threshold = thresholds.get("frac_lt70_95pct", np.nan)
    terminal_threshold = thresholds.get("terminal_lt70_95pct", np.nan)
    coverage_threshold = thresholds.get("coverage_05pct", np.nan)
    if np.isfinite(terminal_threshold):
        terminal_threshold = max(20.0, terminal_threshold)

    all_flags = []
    primary = []
    for _, row in df.iterrows():
        flags = []
        length = row.get("fasta_length")
        if pd.isna(length):
            length = row.get("modeled_residue_count")
        medp = row.get("median_plddt")
        lowfrac = row.get("frac_plddt_lt70")
        nterm = row.get("n_terminal_low_confidence_length")
        cterm = row.get("c_terminal_low_confidence_length")
        longest_core = row.get("longest_high_confidence_segment")
        coverage = row.get("structure_coverage_fraction")
        parse_status = row.get("parsing_status", "")

        short = bool(np.isfinite(length) and np.isfinite(lower_len) and length < lower_len)
        long = bool(np.isfinite(length) and np.isfinite(upper_len) and length > upper_len)
        low_global = bool(np.isfinite(medp) and np.isfinite(low_global_threshold) and medp < low_global_threshold)
        high_lowfrac = bool(np.isfinite(lowfrac) and np.isfinite(high_lowfrac_threshold) and lowfrac > high_lowfrac_threshold)
        long_terminal = bool(
            np.isfinite(terminal_threshold)
            and max(nterm if np.isfinite(nterm) else 0, cterm if np.isfinite(cterm) else 0) > terminal_threshold
        )
        low_coverage = bool(np.isfinite(coverage) and np.isfinite(coverage_threshold) and coverage < coverage_threshold)
        high_core = bool(
            np.isfinite(longest_core) and np.isfinite(length)
            and longest_core >= max(40, 0.4 * length)
        )
        normal_length = bool(np.isfinite(length) and (not np.isfinite(lower_len) or length >= lower_len)
                             and (not np.isfinite(upper_len) or length <= upper_len))
        high_conf = bool(np.isfinite(medp) and medp >= 70 and np.isfinite(lowfrac) and lowfrac <= 0.30)

        if short:
            flags.append("short_relative_to_clustered")
        if long:
            flags.append("long_relative_to_clustered")
        if low_global:
            flags.append("low_global_confidence")
        if high_lowfrac:
            flags.append("high_low-confidence_fraction")
        if long_terminal:
            flags.append("long_low-confidence_terminal_extension")
        if high_core and long_terminal:
            flags.append("high-confidence_core_with_low-confidence_extension")
        if low_coverage:
            flags.append("low_structure_coverage_relative_to_clustered")
        if normal_length and high_conf:
            flags.append("normal-length_high-confidence")
        if row.get("chain_count", 1) and np.isfinite(row.get("chain_count", np.nan)) and row.get("chain_count", 1) > 1:
            flags.append("unexpected_multichain_structure")
        if parse_status not in ("ok", "suspicious_plddt_range"):
            flags.append("missing_or_unparsed_structure")
        if not flags:
            flags.append("unclassified_by_basic_QC")

        if "missing_or_unparsed_structure" in flags:
            cat = "missing_or_unparsed_structure"
        elif (short and (low_global or high_lowfrac or low_coverage)):
            cat = "probable_fragmented_or_low-confidence"
        elif high_core and long_terminal:
            cat = "high-confidence_core_with_low-confidence_extension"
        elif long and (long_terminal or low_coverage):
            cat = "probable_architecture_or_coverage_sensitive"
        elif low_global or high_lowfrac:
            cat = "low-confidence_model"
        elif short:
            cat = "short_relative_to_clustered"
        elif long:
            cat = "long_relative_to_clustered"
        elif normal_length and high_conf:
            cat = "normal-length_high-confidence"
        else:
            cat = "unclassified_by_basic_QC"
        all_flags.append(";".join(flags))
        primary.append(cat)

    df["QC_flags"] = all_flags
    df["primary_QC_category"] = primary
    return df


def quantile_summary(series: pd.Series) -> Dict[str, float]:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if s.empty:
        return {k: np.nan for k in ["n", "median", "q1", "q3", "iqr", "mean", "sd", "min", "p05", "p10", "p90", "p95", "max"]}
    return {
        "n": int(len(s)),
        "median": float(s.median()),
        "q1": float(s.quantile(0.25)),
        "q3": float(s.quantile(0.75)),
        "iqr": float(s.quantile(0.75) - s.quantile(0.25)),
        "mean": float(s.mean()),
        "sd": float(s.std(ddof=1)) if len(s) > 1 else 0.0,
        "min": float(s.min()),
        "p05": float(s.quantile(0.05)),
        "p10": float(s.quantile(0.10)),
        "p90": float(s.quantile(0.90)),
        "p95": float(s.quantile(0.95)),
        "max": float(s.max()),
    }


def build_group_summary(master: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = [
        "fasta_length", "modeled_residue_count", "structure_coverage_fraction",
        "mean_plddt", "median_plddt", "frac_plddt_lt50", "frac_plddt_lt70",
        "longest_segment_lt50", "longest_segment_lt70",
        "n_terminal_low_confidence_length", "c_terminal_low_confidence_length",
        "longest_high_confidence_segment",
    ]
    for group in GROUP_ORDER:
        sub = master[master["cluster_group"] == group]
        for metric in metrics:
            row = {"cluster_group": group, "metric": metric}
            row.update(quantile_summary(sub[metric]))
            rows.append(row)
    return pd.DataFrame(rows)


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) == 0 or len(y) == 0:
        return np.nan
    u = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic").statistic
    return float(2 * u / (len(x) * len(y)) - 1)


def hodges_lehmann(x: np.ndarray, y: np.ndarray, max_pairs: int = 500_000) -> Tuple[float, str]:
    total = len(x) * len(y)
    if total == 0:
        return np.nan, "unavailable"
    if total <= max_pairs:
        diffs = (x[:, None] - y[None, :]).ravel()
        return float(np.median(diffs)), "exact_all_pairwise_differences"
    rng = np.random.default_rng(RNG_SEED)
    idx_x = rng.integers(0, len(x), size=max_pairs)
    idx_y = rng.integers(0, len(y), size=max_pairs)
    return float(np.median(x[idx_x] - y[idx_y])), f"deterministic_sample_{max_pairs}_pairwise_differences"


def bootstrap_effects(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int,
) -> Tuple[float, float, float, float]:
    if len(x) == 0 or len(y) == 0 or n_boot <= 0:
        return (np.nan,) * 4
    rng = np.random.default_rng(RNG_SEED)
    median_diffs = np.empty(n_boot)
    deltas = np.empty(n_boot)
    for i in range(n_boot):
        xb = x[rng.integers(0, len(x), len(x))]
        yb = y[rng.integers(0, len(y), len(y))]
        median_diffs[i] = np.median(xb) - np.median(yb)
        deltas[i] = cliffs_delta(xb, yb)
    return (
        float(np.quantile(median_diffs, 0.025)),
        float(np.quantile(median_diffs, 0.975)),
        float(np.quantile(deltas, 0.025)),
        float(np.quantile(deltas, 0.975)),
    )


def bh_adjust(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    vals = p[valid]
    if vals.size == 0:
        return out
    order = np.argsort(vals)
    ranked = vals[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    back = np.empty_like(adjusted)
    back[order] = adjusted
    out[valid] = back
    return out


def interpret_delta(delta: float) -> str:
    if not np.isfinite(delta):
        return "unavailable"
    a = abs(delta)
    magnitude = "negligible" if a < 0.147 else "small" if a < 0.33 else "medium" if a < 0.474 else "large"
    direction = "higher in first group" if delta > 0 else "lower in first group" if delta < 0 else "no direction"
    return f"{magnitude}; {direction}"


def perform_statistics(master: pd.DataFrame, n_boot: int) -> pd.DataFrame:
    rows = []
    for metric in METRICS_FOR_STATS:
        arrays = []
        for group in GROUP_ORDER:
            arrays.append(pd.to_numeric(master.loc[master["cluster_group"] == group, metric], errors="coerce").dropna().to_numpy(float))
        valid_arrays = [a for a in arrays if len(a)]
        if len(valid_arrays) == 3:
            pooled = np.concatenate(valid_arrays)
            if len(np.unique(pooled)) < 2:
                h, p = np.nan, np.nan
                eps2 = np.nan
            else:
                h, p = stats.kruskal(*valid_arrays)
                n = sum(map(len, valid_arrays))
                k = len(valid_arrays)
                eps2 = max(0.0, (h - k + 1) / (n - k)) if n > k else np.nan
            rows.append({
                "metric": metric,
                "groups_compared": "giant_cluster vs other_multimember vs singleton",
                "test": "Kruskal-Wallis",
                "statistic": float(h),
                "raw_p_value": float(p),
                "adjusted_p_value": np.nan,
                "effect_size_name": "epsilon_squared",
                "effect_size": eps2,
                "effect_ci_lower": np.nan,
                "effect_ci_upper": np.nan,
                "median_first": np.nan,
                "median_second": np.nan,
                "median_difference_first_minus_second": np.nan,
                "hodges_lehmann": np.nan,
                "hodges_lehmann_method": "",
                "median_difference_ci_lower": np.nan,
                "median_difference_ci_upper": np.nan,
                "interpretation": "Omnibus distributional comparison; inspect pairwise effects.",
            })

        group_data = {
            "singleton": pd.to_numeric(master.loc[master["cluster_group"] == "singleton", metric], errors="coerce").dropna().to_numpy(float),
            "giant_cluster": pd.to_numeric(master.loc[master["cluster_group"] == "giant_cluster", metric], errors="coerce").dropna().to_numpy(float),
            "other_multimember": pd.to_numeric(master.loc[master["cluster_group"] == "other_multimember", metric], errors="coerce").dropna().to_numpy(float),
        }
        group_data["all_multimember"] = np.concatenate([group_data["giant_cluster"], group_data["other_multimember"]])
        for g1, g2 in PAIRWISE_COMPARISONS:
            x, y = group_data[g1], group_data[g2]
            if len(x) == 0 or len(y) == 0:
                continue
            result = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
            delta = cliffs_delta(x, y)
            hl, hl_method = hodges_lehmann(x, y)
            md = float(np.median(x) - np.median(y))
            md_lo, md_hi, d_lo, d_hi = bootstrap_effects(x, y, n_boot)
            rows.append({
                "metric": metric,
                "groups_compared": f"{g1} vs {g2}",
                "test": "Mann-Whitney U",
                "statistic": float(result.statistic),
                "raw_p_value": float(result.pvalue),
                "adjusted_p_value": np.nan,
                "effect_size_name": "Cliffs_delta",
                "effect_size": delta,
                "effect_ci_lower": d_lo,
                "effect_ci_upper": d_hi,
                "median_first": float(np.median(x)),
                "median_second": float(np.median(y)),
                "median_difference_first_minus_second": md,
                "hodges_lehmann": hl,
                "hodges_lehmann_method": hl_method,
                "median_difference_ci_lower": md_lo,
                "median_difference_ci_upper": md_hi,
                "interpretation": interpret_delta(delta),
            })
    stats_df = pd.DataFrame(rows)
    pairwise_mask = stats_df["test"] == "Mann-Whitney U"
    stats_df.loc[pairwise_mask, "adjusted_p_value"] = bh_adjust(stats_df.loc[pairwise_mask, "raw_p_value"])
    return stats_df


def save_figure(fig, outbase: Path, dpi: int = 400) -> None:
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def numeric_group_data(master: pd.DataFrame, metric: str):
    return [
        pd.to_numeric(master.loc[master["cluster_group"] == group, metric], errors="coerce").dropna().to_numpy(float)
        for group in GROUP_ORDER
    ]


def plot_violin(master: pd.DataFrame, metric: str, ylabel: str, title: str, outbase: Path, log_y: bool = False) -> None:
    data = numeric_group_data(master, metric)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    positions = np.arange(1, len(GROUP_ORDER) + 1)
    valid = [(p, d) for p, d in zip(positions, data) if len(d)]
    if valid:
        ax.violinplot([d for _, d in valid], positions=[p for p, _ in valid], showmeans=False, showmedians=False, showextrema=False)
        ax.boxplot([d for _, d in valid], positions=[p for p, _ in valid], widths=0.18, showfliers=False,
                   medianprops={"linewidth": 1.5})
        rng = np.random.default_rng(RNG_SEED)
        for p, d in valid:
            sample = d if len(d) <= 600 else rng.choice(d, 600, replace=False)
            jitter = rng.normal(0, 0.045, size=len(sample))
            ax.scatter(np.full(len(sample), p) + jitter, sample, s=7, alpha=0.25, linewidths=0)
    ax.set_xticks(positions, [g.replace("_", "\n") for g in GROUP_ORDER])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if log_y:
        ax.set_yscale("log")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_figure(fig, outbase)


def plot_scatter(master: pd.DataFrame, x: str, y: str, xlabel: str, ylabel: str, title: str, outbase: Path, log_x: bool = False, singleton_only: bool = False) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    groups = ["singleton"] if singleton_only else GROUP_ORDER
    markers = {"giant_cluster": "o", "other_multimember": "^", "singleton": "."}
    for group in groups:
        sub = master[master["cluster_group"] == group]
        xv = pd.to_numeric(sub[x], errors="coerce")
        yv = pd.to_numeric(sub[y], errors="coerce")
        ok = xv.notna() & yv.notna() & (xv > 0 if log_x else True)
        ax.scatter(xv[ok], yv[ok], s=16 if group != "singleton" else 12, alpha=0.45,
                   marker=markers[group], linewidths=0, label=group.replace("_", " "))
    if log_x:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if not singleton_only:
        ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_figure(fig, outbase)


def plot_singleton_qc(master: pd.DataFrame, outbase: Path) -> None:
    sub = master[master["cluster_group"] == "singleton"]
    counts = sub["primary_QC_category"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(8.0, max(4.0, 0.45 * len(counts) + 1.5)))
    bars = ax.barh(counts.index, counts.values)
    total = len(sub)
    for bar, count in zip(bars, counts.values):
        frac = count / total if total else np.nan
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                f"  {count} ({frac:.1%})", va="center", fontsize=9)
    ax.set_xlabel("Singleton proteins")
    ax.set_title("Singleton composition by exploratory QC category")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_figure(fig, outbase)


def select_profile_representatives(master: pd.DataFrame, max_profiles: int = 8) -> List[str]:
    sub = master[(master["cluster_group"] == "singleton") & (master["parsing_status"] == "ok")].copy()
    selected = []
    priority = [
        "probable_fragmented_or_low-confidence",
        "high-confidence_core_with_low-confidence_extension",
        "probable_architecture_or_coverage_sensitive",
        "low-confidence_model",
        "normal-length_high-confidence",
        "unclassified_by_basic_QC",
    ]
    for category in priority:
        cat = sub[sub["primary_QC_category"] == category].copy()
        if cat.empty:
            continue
        features = ["fasta_length", "median_plddt", "frac_plddt_lt70"]
        for feat in features:
            cat[feat] = pd.to_numeric(cat[feat], errors="coerce")
        med = cat[features].median(numeric_only=True)
        scale = cat[features].mad(numeric_only=True) if hasattr(cat[features], "mad") else cat[features].std()
        scale = scale.replace(0, 1).fillna(1)
        cat["distance"] = ((cat[features] - med).abs() / scale).sum(axis=1)
        selected.append(cat.sort_values(["distance", "protein_id"]).iloc[0]["protein_id"])
        if len(selected) >= max_profiles:
            break
    return selected


def plot_profiles(master: pd.DataFrame, profiles: Dict[str, List[Tuple[str, int, float]]], outbase: Path) -> None:
    selected = [x for x in select_profile_representatives(master) if x in profiles]
    if not selected:
        return
    n = len(selected)
    fig, axes = plt.subplots(n, 1, figsize=(8.5, 2.2 * n), sharex=False)
    if n == 1:
        axes = [axes]
    for ax, pid in zip(axes, selected):
        prof = profiles[pid]
        y = [x[2] for x in prof]
        x = np.arange(1, len(y) + 1)
        ax.plot(x, y, linewidth=1.0)
        ax.axhline(70, linestyle="--", linewidth=0.8)
        ax.axhline(50, linestyle=":", linewidth=0.8)
        row = master.loc[master["protein_id"] == pid].iloc[0]
        ax.set_ylim(0, 100)
        ax.set_ylabel("pLDDT")
        ax.set_title(f"{pid} | {row['primary_QC_category']} | length={row.get('fasta_length', np.nan):.0f}", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlabel("Modeled residue position (chains concatenated only for display)")
    fig.tight_layout()
    save_figure(fig, outbase)


def create_plots(master: pd.DataFrame, profiles: Dict[str, List[Tuple[str, int, float]]], figures_dir: Path) -> None:
    plot_violin(master, "fasta_length", "Sequence length (aa)", "Sequence-length distribution", figures_dir / "plot1_sequence_length_raw")
    positive = master.copy()
    positive["fasta_length"] = pd.to_numeric(positive["fasta_length"], errors="coerce")
    positive = positive[positive["fasta_length"] > 0]
    plot_violin(positive, "fasta_length", "Sequence length (aa; log scale)", "Sequence-length distribution (log scale)", figures_dir / "plot1_sequence_length_log", log_y=True)
    plot_violin(master, "median_plddt", "Median pLDDT", "Predicted-structure confidence", figures_dir / "plot2_median_plddt")
    plot_violin(master, "frac_plddt_lt70", "Fraction of residues with pLDDT < 70", "Low-confidence residue fraction", figures_dir / "plot3_fraction_plddt_lt70")
    plot_scatter(master, "fasta_length", "median_plddt", "Sequence length (aa)", "Median pLDDT", "Sequence length versus median pLDDT", figures_dir / "plot4_length_vs_median_plddt_raw")
    plot_scatter(master, "fasta_length", "median_plddt", "Sequence length (aa; log scale)", "Median pLDDT", "Sequence length versus median pLDDT", figures_dir / "plot4_length_vs_median_plddt_log", log_x=True)
    plot_scatter(master, "fasta_length", "median_plddt", "Sequence length (aa; log scale)", "Median pLDDT", "Singleton-focused: length versus median pLDDT", figures_dir / "plot4_singleton_length_vs_median_plddt", log_x=True, singleton_only=True)
    plot_scatter(master, "fasta_length", "frac_plddt_lt70", "Sequence length (aa; log scale)", "Fraction pLDDT < 70", "Sequence length versus low-confidence fraction", figures_dir / "plot5_length_vs_low_confidence_fraction", log_x=True)
    plot_singleton_qc(master, figures_dir / "plot6_singleton_qc_composition")
    plot_profiles(master, profiles, figures_dir / "plot7_representative_plddt_profiles")


def build_missing_data_audit(
    cluster: pd.DataFrame,
    fasta_ids: set,
    structure_match_df: pd.DataFrame,
    master: pd.DataFrame,
    fasta_duplicates: Sequence[str],
    unused_fasta: Sequence[str],
    unused_structures: Sequence[str],
    cluster_audit: Dict[str, object],
) -> pd.DataFrame:
    rows = []
    cluster_ids = set(cluster["member_id"])
    for pid in sorted(cluster_ids - fasta_ids):
        rows.append({"audit_type": "missing_fasta", "protein_id": pid, "details": "Cluster member absent from FASTA."})
    for pid in sorted(fasta_duplicates):
        rows.append({"audit_type": "duplicated_fasta_id", "protein_id": pid, "details": "Duplicate FASTA identifier; first sequence retained."})
    for pid in sorted(unused_fasta):
        rows.append({"audit_type": "fasta_not_in_cluster_membership", "protein_id": pid, "details": "FASTA entry not present as a cluster member."})
    matched_nonempty = structure_match_df.loc[structure_match_df["matched_structure"] != "", "matched_structure"]
    duplicated_paths = matched_nonempty[matched_nonempty.duplicated(keep=False)]
    for path, count in duplicated_paths.value_counts().items():
        ids = structure_match_df.loc[structure_match_df["matched_structure"] == path, "protein_id"].tolist()
        rows.append({"audit_type": "structure_file_matched_to_multiple_ids", "protein_id": ";".join(ids), "details": f"{path} matched to {count} IDs."})
    for _, r in structure_match_df.iterrows():
        if r["status"] == "missing":
            rows.append({"audit_type": "missing_structure", "protein_id": r["protein_id"], "details": "No structure filename matched."})
        elif r["status"] == "ambiguous":
            rows.append({"audit_type": "ambiguous_structure_match", "protein_id": r["protein_id"], "details": r["candidate_structures"]})
    for _, r in master.iterrows():
        if r["parsing_status"] == "parse_failure":
            rows.append({"audit_type": "structure_parse_failure", "protein_id": r["protein_id"], "details": r["notes"]})
        elif r["parsing_status"] == "suspicious_plddt_range":
            rows.append({"audit_type": "suspicious_plddt_range", "protein_id": r["protein_id"], "details": r["notes"]})
        if pd.notna(r.get("chain_count")) and r.get("chain_count", 1) > 1:
            rows.append({"audit_type": "multiple_chains", "protein_id": r["protein_id"], "details": r.get("chain_residue_counts", "")})
        if pd.notna(r.get("length_difference_fasta_minus_modeled")) and r.get("length_difference_fasta_minus_modeled") < 0:
            rows.append({"audit_type": "modeled_count_exceeds_fasta_length", "protein_id": r["protein_id"], "details": str(r.get("length_difference_fasta_minus_modeled"))})
    for path in sorted(unused_structures):
        rows.append({"audit_type": "structure_not_in_cluster_membership", "protein_id": "", "details": path})
    for pid in cluster_audit.get("conflicting_member_ids", []):
        rows.append({"audit_type": "member_assigned_to_multiple_representatives", "protein_id": pid, "details": "Conflicting cluster assignments retained in audit; master uses first deterministic assignment."})
    return pd.DataFrame(rows, columns=["audit_type", "protein_id", "details"])


def build_prioritized_singletons(master: pd.DataFrame) -> pd.DataFrame:
    sub = master[master["cluster_group"] == "singleton"].copy()
    category_rank = {
        "probable_fragmented_or_low-confidence": 1,
        "low-confidence_model": 1,
        "missing_or_unparsed_structure": 1,
        "probable_architecture_or_coverage_sensitive": 2,
        "high-confidence_core_with_low-confidence_extension": 2,
        "short_relative_to_clustered": 2,
        "long_relative_to_clustered": 2,
        "normal-length_high-confidence": 3,
        "unclassified_by_basic_QC": 3,
    }
    category_label = {
        1: "probable_fragmented_or_low-confidence_candidates",
        2: "probable_architecture_or_coverage-sensitive_candidates",
        3: "normal-length_or_high-confidence_unexplained_singletons",
    }
    sub["priority_group_rank"] = sub["primary_QC_category"].map(category_rank).fillna(3).astype(int)
    sub["priority_group"] = sub["priority_group_rank"].map(category_label)
    sub["direct_rescue_priority"] = np.where(
        sub["priority_group_rank"] == 3,
        "high",
        np.where(sub["priority_group_rank"] == 2, "medium", "low")
    )
    sort_cols = ["priority_group_rank", "median_plddt", "fasta_length", "protein_id"]
    return sub.sort_values(sort_cols, ascending=[True, False, True, True])


def write_summary_and_results(
    master: pd.DataFrame,
    stats_df: pd.DataFrame,
    cluster_summary: Dict[str, object],
    thresholds: Dict[str, float],
    outdir: Path,
) -> None:
    singleton = master[master["cluster_group"] == "singleton"]
    multi = master[master["cluster_group"].isin(["giant_cluster", "other_multimember"])]
    qc_counts = singleton["primary_QC_category"].value_counts()
    n_singleton = len(singleton)
    explained_categories = [
        "probable_fragmented_or_low-confidence",
        "low-confidence_model",
        "probable_architecture_or_coverage_sensitive",
        "high-confidence_core_with_low-confidence_extension",
        "short_relative_to_clustered",
        "long_relative_to_clustered",
    ]
    explained_n = int(qc_counts.reindex(explained_categories).fillna(0).sum())
    normal_hq_n = int(qc_counts.get("normal-length_high-confidence", 0))
    unexplained_n = int(qc_counts.get("unclassified_by_basic_QC", 0) + normal_hq_n)

    def med(df, col):
        return pd.to_numeric(df[col], errors="coerce").median()

    summary_lines = [
        "# Singleton length/pLDDT analysis summary",
        "",
        "## Cluster reconstruction",
        "",
    ]
    for k, v in cluster_summary.items():
        summary_lines.append(f"- **{k}**: {v}")
    summary_lines += [
        "",
        "## Reference-derived exploratory thresholds",
        "",
    ]
    for k, v in thresholds.items():
        summary_lines.append(f"- **{k}**: {v}")
    summary_lines += [
        "",
        "## Descriptive findings",
        "",
        f"- Singleton median sequence length: {med(singleton, 'fasta_length'):.2f} aa.",
        f"- All multi-member median sequence length: {med(multi, 'fasta_length'):.2f} aa.",
        f"- Singleton median pLDDT: {med(singleton, 'median_plddt'):.2f}.",
        f"- All multi-member median pLDDT: {med(multi, 'median_plddt'):.2f}.",
        f"- Singleton categories consistent with atypical length, low confidence, or architecture/coverage sensitivity: {explained_n}/{n_singleton} ({explained_n/n_singleton:.1%})" if n_singleton else "- No singleton proteins were present.",
        f"- Normal-length, high-confidence singletons: {normal_hq_n}/{n_singleton} ({normal_hq_n/n_singleton:.1%})" if n_singleton else "- No singleton proteins were present.",
        f"- Basic length/confidence-unexplained singleton subset: {unexplained_n}/{n_singleton} ({unexplained_n/n_singleton:.1%})" if n_singleton else "- No singleton proteins were present.",
        "",
        "## Interpretation guardrails",
        "",
        "These categories are structure-quality and architecture descriptors, not GT/non-GT assignments. "
        "A low-confidence or atypical-length model may contribute to singleton formation, whereas normal-length high-confidence singletons remain unexplained by these basic metrics and require direct structural rescue and functional annotation.",
    ]
    (outdir / "analysis_summary.md").write_text("\n".join(summary_lines) + "\n")

    # Provisional Results wording, intentionally conditional and conservative.
    result_text = (
        f"Foldseek coarse-clustering partitioned {cluster_summary.get('total_unique_protein_count', 'NA')} unique candidate proteins into "
        f"{cluster_summary.get('total_cluster_count', 'NA')} clusters, including {cluster_summary.get('singleton_cluster_count', 'NA')} singletons. "
        f"Singletons had a median sequence length of {med(singleton, 'fasta_length'):.1f} aa and median model confidence of "
        f"{med(singleton, 'median_plddt'):.1f} pLDDT, compared with {med(multi, 'fasta_length'):.1f} aa and "
        f"{med(multi, 'median_plddt'):.1f} pLDDT among proteins assigned to multi-member clusters. "
        f"Using robust thresholds derived from the multi-member population, {explained_n} of {n_singleton} singletons "
        f"({explained_n/n_singleton:.1%}) showed features consistent with atypical length, globally reduced confidence, incomplete structural coverage, "
        f"or a high-confidence core accompanied by an extended low-confidence region. These observations suggest that sequence/structure quality and architecture "
        f"contribute to singleton formation but do not establish functional status."
        if n_singleton else "No singleton proteins were identified in the supplied cluster membership."
    )
    result_text2 = (
        f"Importantly, {unexplained_n} singletons ({unexplained_n/n_singleton:.1%}) were not explained by the basic length and confidence criteria, including "
        f"{normal_hq_n} normal-length, high-confidence models. This residual population should therefore be retained for direct Foldseek rescue searches, "
        f"domain-aware structural comparisons, and functional annotation rather than being treated as false-positive GT candidates."
        if n_singleton else ""
    )
    (outdir / "provisional_results_text.md").write_text(result_text + "\n\n" + result_text2 + "\n")



# ---------------------------------------------------------------------------
# v3 corrections
# ---------------------------------------------------------------------------
# These definitions intentionally override the v2 implementations above.
# They (i) resolve redundant nested copies of the same structure ID, and
# (ii) replace collapsed outlier fences with empirical length strata plus
# standard pLDDT confidence classes. These are QC descriptors, not functional
# GT/non-GT calls.

def _sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    import hashlib
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            block = fh.read(block_size)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def match_structure_files(
    protein_ids: Sequence[str],
    structure_files: Sequence[Path],
    explicit_map: Dict[str, Path],
) -> Tuple[Dict[str, Optional[Path]], pd.DataFrame]:
    """Match IDs while handling nested duplicate copies conservatively.

    Resolution order for multiple exact/prefix matches:
      1. an explicit mapping;
      2. a unique candidate at the shallowest directory depth;
      3. byte-identical candidates (choose the shallowest deterministically);
      4. otherwise retain an ambiguous match for manual resolution.

    The common case ``STRUCTURES_DIR/ID.pdb`` plus
    ``STRUCTURES_DIR/selected_hits/ID.pdb`` is therefore resolved to the
    top-level canonical file rather than being discarded.
    """
    by_stem: Dict[str, List[Path]] = defaultdict(list)
    for raw in structure_files:
        p = Path(raw)
        by_stem[strip_structure_suffix(p.name)].append(p)

    matches: Dict[str, Optional[Path]] = {}
    rows = []
    for pid in protein_ids:
        candidates: List[Path] = []
        method = "none"
        resolution = ""
        if pid in explicit_map:
            candidates = [Path(explicit_map[pid])]
            method = "explicit_map"
        elif pid in by_stem:
            candidates = by_stem[pid]
            method = "exact_stem"
        else:
            pattern = re.compile(rf"^{re.escape(pid)}(?:[._-]|$)")
            candidates = [p for stem, paths in by_stem.items() if pattern.search(stem) for p in paths]
            if candidates:
                method = "prefix_stem"

        candidates = sorted(set(candidates), key=lambda x: (len(x.parts), str(x)))
        chosen: Optional[Path] = None
        status = "missing"
        if len(candidates) == 1:
            chosen = candidates[0]
            status = "matched"
            resolution = "unique_candidate"
        elif len(candidates) > 1:
            min_depth = min(len(p.parts) for p in candidates)
            shallow = [p for p in candidates if len(p.parts) == min_depth]
            if len(shallow) == 1:
                chosen = shallow[0]
                status = "matched_preferred_shallowest"
                resolution = "unique_shallowest_path"
            else:
                # Only hash the small ambiguous candidate set.
                try:
                    sizes = {p.stat().st_size for p in candidates}
                    identical = len(sizes) == 1 and len({_sha256_file(p) for p in candidates}) == 1
                except OSError:
                    identical = False
                if identical:
                    chosen = candidates[0]
                    status = "matched_identical_duplicates"
                    resolution = "byte_identical_duplicate_files"
                else:
                    status = "ambiguous"
                    resolution = "nonidentical_or_equal_priority_candidates"

        matches[pid] = chosen
        rows.append({
            "protein_id": pid,
            "match_method": method,
            "candidate_count": len(candidates),
            "matched_structure": str(chosen) if chosen else "",
            "candidate_structures": ";".join(map(str, candidates)),
            "status": status,
            "duplicate_resolution": resolution,
        })
    return matches, pd.DataFrame(rows)


def robust_reference_thresholds(master: pd.DataFrame) -> Dict[str, float]:
    """Derive empirical reference strata from all multi-member proteins.

    Unlike Tukey/MAD lower outlier fences, empirical percentiles do not
    collapse to 1 aa when the clustered population is broad and right-skewed.
    """
    ref = master[master["cluster_group"].isin(["giant_cluster", "other_multimember"])].copy()
    lengths = pd.to_numeric(ref["fasta_length"], errors="coerce").dropna()
    if lengths.empty:
        lengths = pd.to_numeric(ref["modeled_residue_count"], errors="coerce").dropna()
    plddt = pd.to_numeric(ref["median_plddt"], errors="coerce").dropna()
    lowfrac = pd.to_numeric(ref["frac_plddt_lt70"], errors="coerce").dropna()
    terminal_fraction = pd.concat([
        pd.to_numeric(ref["n_terminal_low_confidence_length"], errors="coerce") /
        pd.to_numeric(ref["fasta_length"], errors="coerce"),
        pd.to_numeric(ref["c_terminal_low_confidence_length"], errors="coerce") /
        pd.to_numeric(ref["fasta_length"], errors="coerce"),
    ]).replace([np.inf, -np.inf], np.nan).dropna()
    coverage = pd.to_numeric(ref["structure_coverage_fraction"], errors="coerce").dropna()

    def q(series: pd.Series, value: float) -> float:
        return float(series.quantile(value)) if len(series) else np.nan

    return {
        "reference_population_n": int(len(ref)),
        "reference_length_n": int(len(lengths)),
        "reference_confidence_n": int(len(plddt)),
        "length_p05": q(lengths, 0.05),
        "length_p10": q(lengths, 0.10),
        "length_q1": q(lengths, 0.25),
        "length_median": q(lengths, 0.50),
        "length_q3": q(lengths, 0.75),
        "length_p90": q(lengths, 0.90),
        "length_p95": q(lengths, 0.95),
        "median_plddt_p10": q(plddt, 0.10),
        "median_plddt_q1": q(plddt, 0.25),
        "median_plddt_reference_median": q(plddt, 0.50),
        "median_plddt_q3": q(plddt, 0.75),
        "frac_lt70_reference_median": q(lowfrac, 0.50),
        "frac_lt70_q3": q(lowfrac, 0.75),
        "terminal_low_confidence_fraction_p95": q(terminal_fraction, 0.95),
        "coverage_05pct": q(coverage, 0.05),
        "classification_length_anchor": "multi-member empirical P10/Q1/Q3/P90",
        "classification_confidence_anchor": "standard median pLDDT bins: <50, 50-<70, 70-<90, >=90",
    }


def classify_qc(master: pd.DataFrame, thresholds: Dict[str, float]) -> pd.DataFrame:
    df = master.copy()
    p10 = thresholds.get("length_p10", np.nan)
    q1 = thresholds.get("length_q1", np.nan)
    q3 = thresholds.get("length_q3", np.nan)
    p90 = thresholds.get("length_p90", np.nan)
    terminal_p95 = thresholds.get("terminal_low_confidence_fraction_p95", np.nan)

    flags_all: List[str] = []
    primary: List[str] = []
    length_classes: List[str] = []
    confidence_classes: List[str] = []
    combined_classes: List[str] = []

    for _, row in df.iterrows():
        length = pd.to_numeric(pd.Series([row.get("fasta_length")]), errors="coerce").iloc[0]
        medp = pd.to_numeric(pd.Series([row.get("median_plddt")]), errors="coerce").iloc[0]
        lowfrac = pd.to_numeric(pd.Series([row.get("frac_plddt_lt70")]), errors="coerce").iloc[0]
        nterm = pd.to_numeric(pd.Series([row.get("n_terminal_low_confidence_length")]), errors="coerce").iloc[0]
        cterm = pd.to_numeric(pd.Series([row.get("c_terminal_low_confidence_length")]), errors="coerce").iloc[0]
        longest_core = pd.to_numeric(pd.Series([row.get("longest_high_confidence_segment")]), errors="coerce").iloc[0]
        parse_status = str(row.get("parsing_status", ""))
        flags: List[str] = []

        if not np.isfinite(length):
            length_class = "missing_length"
        elif np.isfinite(p10) and length < p10:
            length_class = "below_multimember_p10"
            flags.append("below_multimember_p10_length")
        elif np.isfinite(q1) and length < q1:
            length_class = "multimember_p10_to_q1"
            flags.append("below_multimember_q1_length")
        elif np.isfinite(q3) and length <= q3:
            length_class = "multimember_interquartile_length"
        elif np.isfinite(p90) and length <= p90:
            length_class = "multimember_q3_to_p90"
        else:
            length_class = "above_multimember_p90"
            flags.append("above_multimember_p90_length")

        if parse_status not in ("ok", "suspicious_plddt_range") or not np.isfinite(medp):
            confidence_class = "missing_or_unparsed_confidence"
            flags.append("missing_or_unparsed_structure")
        elif medp < 50:
            confidence_class = "very_low_median_plddt_lt50"
            flags.append("median_plddt_lt50")
        elif medp < 70:
            confidence_class = "low_median_plddt_50_to_lt70"
            flags.append("median_plddt_lt70")
        elif medp < 90:
            confidence_class = "confident_median_plddt_70_to_lt90"
        else:
            confidence_class = "very_high_median_plddt_ge90"
            flags.append("median_plddt_ge90")

        if np.isfinite(lowfrac) and lowfrac > 0.50:
            flags.append("majority_residues_plddt_lt70")

        terminal_fraction = np.nan
        if np.isfinite(length) and length > 0:
            terminal_fraction = max(nterm if np.isfinite(nterm) else 0.0,
                                    cterm if np.isfinite(cterm) else 0.0) / length
        # Use a length-normalized terminal criterion. The empirical p95 can
        # equal 1.0 in very low-confidence reference models and is therefore
        # recorded for provenance but not used as the sole cutoff.
        terminal_extension = bool(
            np.isfinite(terminal_fraction)
            and terminal_fraction >= 0.25
            and max(nterm if np.isfinite(nterm) else 0.0,
                    cterm if np.isfinite(cterm) else 0.0) >= 20
        )
        high_core = bool(
            np.isfinite(longest_core) and np.isfinite(length)
            and longest_core >= max(40.0, 0.40 * length)
        )
        if terminal_extension:
            flags.append("long_low_confidence_terminal_fraction")
        if terminal_extension and high_core:
            flags.append("high_confidence_core_with_low_confidence_extension")

        short_p10 = bool(np.isfinite(length) and np.isfinite(p10) and length < p10)
        confident = bool(np.isfinite(medp) and medp >= 70 and parse_status in ("ok", "suspicious_plddt_range"))
        if "missing_or_unparsed_structure" in flags:
            combined = "missing_or_unparsed_confidence"
        elif short_p10 and not confident:
            combined = "below_multimember_p10_and_low_confidence"
        elif short_p10 and confident:
            combined = "below_multimember_p10_but_confident"
        elif (not short_p10) and not confident:
            combined = "at_or_above_multimember_p10_but_low_confidence"
        else:
            combined = "at_or_above_multimember_p10_and_confident"

        # Primary category retains architecture information when present, but
        # otherwise follows the transparent two-axis classification.
        if "missing_or_unparsed_structure" in flags:
            cat = "missing_or_unparsed_structure"
        elif terminal_extension and high_core:
            cat = "high-confidence_core_with_low-confidence_extension"
        else:
            cat = combined

        flags_all.append(";".join(flags) if flags else "no_basic_qc_flag")
        primary.append(cat)
        length_classes.append(length_class)
        confidence_classes.append(confidence_class)
        combined_classes.append(combined)

    df["relative_length_class"] = length_classes
    df["confidence_class"] = confidence_classes
    df["combined_length_confidence_class"] = combined_classes
    df["QC_flags"] = flags_all
    df["primary_QC_category"] = primary
    return df


def plot_singleton_qc(master: pd.DataFrame, outbase: Path) -> None:
    sub = master[master["cluster_group"] == "singleton"]
    order = [
        "below_multimember_p10_and_low_confidence",
        "below_multimember_p10_but_confident",
        "at_or_above_multimember_p10_but_low_confidence",
        "at_or_above_multimember_p10_and_confident",
        "missing_or_unparsed_confidence",
    ]
    counts = sub["combined_length_confidence_class"].value_counts().reindex(order).dropna()
    counts = counts[counts > 0]
    fig, ax = plt.subplots(figsize=(9.2, max(4.2, 0.55 * len(counts) + 1.8)))
    bars = ax.barh(counts.index, counts.values)
    total = len(sub)
    for bar, count in zip(bars, counts.values):
        frac = count / total if total else np.nan
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
                f"  {int(count)} ({frac:.1%})", va="center", fontsize=9)
    ax.set_xlabel("Singleton proteins")
    ax.set_title("Singleton stratification by multi-member P10 length and median pLDDT 70")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_figure(fig, outbase)


def select_profile_representatives(master: pd.DataFrame, max_profiles: int = 8) -> List[str]:
    sub = master[(master["cluster_group"] == "singleton") & (master["parsing_status"] == "ok")].copy()
    selected: List[str] = []
    priority = [
        "below_multimember_p10_and_low_confidence",
        "below_multimember_p10_but_confident",
        "at_or_above_multimember_p10_but_low_confidence",
        "at_or_above_multimember_p10_and_confident",
        "high-confidence_core_with_low-confidence-extension",
    ]
    for category in priority:
        cat_sub = sub[sub["primary_QC_category"] == category]
        if cat_sub.empty:
            continue
        cat_sub = cat_sub.sort_values(["median_plddt", "fasta_length"], ascending=[False, False])
        selected.append(cat_sub.iloc[0]["protein_id"])
        if len(cat_sub) > 1:
            selected.append(cat_sub.iloc[-1]["protein_id"])
        if len(selected) >= max_profiles:
            break
    return list(dict.fromkeys(selected))[:max_profiles]


def build_prioritized_singletons(master: pd.DataFrame) -> pd.DataFrame:
    sub = master[master["cluster_group"] == "singleton"].copy()
    mapping = {
        "at_or_above_multimember_p10_and_confident": (1, "tier_1_reference_length_or_longer_confident", "high"),
        "below_multimember_p10_but_confident": (2, "tier_2_compact_but_confident", "medium_high"),
        "at_or_above_multimember_p10_but_low_confidence": (3, "tier_3_reference_length_or_longer_low_confidence", "medium"),
        "below_multimember_p10_and_low_confidence": (4, "tier_4_very_short_and_low_confidence", "low"),
        "missing_or_unparsed_confidence": (5, "tier_5_missing_confidence", "low"),
    }
    ranks, labels, priorities = [], [], []
    for value in sub["combined_length_confidence_class"]:
        rank, label, priority = mapping.get(value, (3, "tier_3_other", "medium"))
        ranks.append(rank); labels.append(label); priorities.append(priority)
    sub["priority_group_rank"] = ranks
    sub["priority_group"] = labels
    sub["direct_rescue_priority"] = priorities
    return sub.sort_values(
        ["priority_group_rank", "median_plddt", "fasta_length", "protein_id"],
        ascending=[True, False, False, True],
    )


def write_summary_and_results(
    master: pd.DataFrame,
    stats_df: pd.DataFrame,
    cluster_summary: Dict[str, object],
    thresholds: Dict[str, float],
    outdir: Path,
) -> None:
    singleton = master[master["cluster_group"] == "singleton"].copy()
    multi = master[master["cluster_group"].isin(["giant_cluster", "other_multimember"])].copy()
    other = master[master["cluster_group"] == "other_multimember"].copy()
    n = len(singleton)

    def med(df: pd.DataFrame, col: str) -> float:
        return float(pd.to_numeric(df[col], errors="coerce").median())
    def count(mask) -> int:
        return int(pd.Series(mask, index=singleton.index).fillna(False).sum())
    def pct(value: int) -> str:
        return f"{value/n:.1%}" if n else "NA"

    p10 = thresholds.get("length_p10", np.nan)
    q1 = thresholds.get("length_q1", np.nan)
    s_len = pd.to_numeric(singleton["fasta_length"], errors="coerce")
    s_plddt = pd.to_numeric(singleton["median_plddt"], errors="coerce")
    below_p10 = count(s_len < p10)
    below_q1 = count(s_len < q1)
    low70 = count(s_plddt < 70)
    confident = count(s_plddt >= 70)
    tier1 = count((s_len >= p10) & (s_plddt >= 70))
    compact_confident = count((s_len < p10) & (s_plddt >= 70))
    joint_low = count((s_len < p10) & (s_plddt < 70))

    # Pull the most interpretable singleton-versus-other-multimember effects.
    def stat_row(metric: str):
        hit = stats_df[(stats_df["metric"] == metric) &
                       (stats_df["groups_compared"] == "singleton vs other_multimember")]
        return hit.iloc[0] if len(hit) else None
    length_stat = stat_row("fasta_length")
    plddt_stat = stat_row("median_plddt")

    summary_lines = [
        "# Corrected singleton length/pLDDT analysis summary (v3)", "",
        "## Cluster reconstruction", "",
    ]
    for k, v in cluster_summary.items():
        summary_lines.append(f"- **{k}**: {v}")
    summary_lines += ["", "## Empirical reference thresholds", ""]
    for k, v in thresholds.items():
        summary_lines.append(f"- **{k}**: {v}")
    summary_lines += [
        "", "## Main findings", "",
        f"- Singleton median sequence length: {med(singleton, 'fasta_length'):.1f} aa; all multi-member: {med(multi, 'fasta_length'):.1f} aa; other multi-member: {med(other, 'fasta_length'):.1f} aa.",
        f"- Singleton median pLDDT: {med(singleton, 'median_plddt'):.1f}; all multi-member: {med(multi, 'median_plddt'):.1f}; other multi-member: {med(other, 'median_plddt'):.1f}.",
        f"- {below_p10}/{n} ({pct(below_p10)}) singletons were shorter than the 10th percentile of the multi-member length distribution ({p10:.0f} aa).",
        f"- {below_q1}/{n} ({pct(below_q1)}) singletons were shorter than the first quartile of the multi-member length distribution ({q1:.0f} aa).",
        f"- {low70}/{n} ({pct(low70)}) singletons had median pLDDT <70, whereas {confident}/{n} ({pct(confident)}) had median pLDDT >=70.",
        f"- Joint extreme-short/low-confidence stratum: {joint_low}/{n} ({pct(joint_low)}).",
        f"- Tier-1 direct-rescue stratum (length >= multi-member P10 and median pLDDT >=70): {tier1}/{n} ({pct(tier1)}).",
        f"- Compact but confident stratum (length < multi-member P10 and median pLDDT >=70): {compact_confident}/{n} ({pct(compact_confident)}).",
    ]
    if length_stat is not None:
        summary_lines.append(
            f"- Singleton versus other-multimember length effect: Cliff's delta {float(length_stat['effect_size']):.3f} ({length_stat['interpretation']}), median difference {float(length_stat['median_difference_first_minus_second']):.1f} aa."
        )
    if plddt_stat is not None:
        summary_lines.append(
            f"- Singleton versus other-multimember median-pLDDT effect: Cliff's delta {float(plddt_stat['effect_size']):.3f} ({plddt_stat['interpretation']}), median difference {float(plddt_stat['median_difference_first_minus_second']):.1f}."
        )
    summary_lines += [
        "", "## Interpretation", "",
        "Sequence length is the dominant QC correlate of singleton status. Confidence is lower on average, but the difference relative to other non-giant multi-member clusters is comparatively small. These descriptors do not establish GT/non-GT function. High-confidence compact singletons may represent real small domains or viral proteins, while reference-length/high-confidence singletons are the strongest candidates for direct structural rescue.",
        "", "## Guardrail", "",
        "Do not use the v2 statement that 95.4% of singletons were unexplained. That value arose because the lower length outlier fence collapsed to 1 aa and the confidence cutoffs became excessively permissive in a heterogeneous reference population.",
    ]
    (outdir / "analysis_summary.md").write_text("\n".join(summary_lines) + "\n")

    results = (
        f"Foldseek coarse clustering partitioned {cluster_summary.get('total_unique_protein_count', 'NA')} candidate proteins into "
        f"{cluster_summary.get('total_cluster_count', 'NA')} clusters, including {cluster_summary.get('singleton_cluster_count', 'NA')} singletons. "
        f"Singletons were strongly shifted toward shorter proteins (median {med(singleton, 'fasta_length'):.0f} aa) relative to other non-giant multi-member clusters "
        f"(median {med(other, 'fasta_length'):.0f} aa). {below_p10} of {n} singletons ({pct(below_p10)}) fell below the 10th percentile of the complete multi-member length distribution ({p10:.0f} aa). "
        f"Predicted-structure confidence was also lower among singletons (median pLDDT {med(singleton, 'median_plddt'):.1f}) than among other multi-member proteins "
        f"({med(other, 'median_plddt'):.1f}), although this effect was substantially smaller than the length shift. "
        f"A joint {joint_low}/{n} ({pct(joint_low)}) subset was both shorter than the multi-member 10th percentile and below median pLDDT 70, consistent with short ORFs and/or incompletely resolved folds contributing to singleton formation. "
        f"Conversely, {tier1}/{n} ({pct(tier1)}) singletons were at least as long as the multi-member 10th-percentile boundary and had median pLDDT >=70; these models remain plausible structurally divergent candidates and were prioritized for direct structural rescue. "
        "Length and pLDDT were used only as structure-quality descriptors and not as functional GT/non-GT classifiers."
    )
    (outdir / "provisional_results_text.md").write_text(results + "\n")

def main() -> int:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--cluster-tsv", required=True, type=Path)
    parser.add_argument("--fasta", type=Path, default=None)
    parser.add_argument("--structures-dir", type=Path, default=None)
    parser.add_argument("--structure-id-map", type=Path, default=None,
                        help="Optional TSV/CSV mapping protein_id to structure_file/path.")
    parser.add_argument("--rep-column", default=None)
    parser.add_argument("--member-column", default=None)
    parser.add_argument("--outdir", type=Path, default=Path("singleton_length_plddt_analysis"))
    parser.add_argument("--bootstrap", type=int, default=500,
                        help="Bootstrap replicates for pairwise median-difference and Cliff's-delta CIs.")
    parser.add_argument("--skip-plots", action="store_true")
    args = parser.parse_args()

    outdir = args.outdir.resolve()
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"
    logs_dir = outdir / "logs"
    intermediate_dir = outdir / "intermediate"
    for d in [outdir, tables_dir, figures_dir, logs_dir, intermediate_dir]:
        d.mkdir(parents=True, exist_ok=True)

    with (logs_dir / "run.log").open("w") as runlog:
        log("[1/8] Reading and auditing Foldseek cluster membership", runlog)
        cluster, cluster_audit = load_cluster_membership(args.cluster_tsv, args.rep_column, args.member_column)
        cluster.to_csv(intermediate_dir / "cluster_membership_deduplicated.tsv", sep="\t", index=False)
        (logs_dir / "cluster_input_audit.json").write_text(json.dumps(cluster_audit, indent=2))

        # Resolve conflicting assignments deterministically while retaining the audit.
        cluster_master = cluster.sort_values(["member_id", "representative_id"]).drop_duplicates("member_id", keep="first")
        sizes = cluster_master.groupby("representative_id")["member_id"].nunique().sort_values(ascending=False)
        if sizes.empty:
            raise ValueError("No valid cluster membership rows remain after parsing.")
        max_size = int(sizes.max())
        giant_candidates = sorted(sizes[sizes == max_size].index.tolist())
        giant_rep = giant_candidates[0]
        cluster_master["cluster_size"] = cluster_master["representative_id"].map(sizes).astype(int)
        cluster_master["cluster_group"] = np.where(
            cluster_master["representative_id"] == giant_rep,
            "giant_cluster",
            np.where(cluster_master["cluster_size"] >= 2, "other_multimember", "singleton"),
        )

        log("[2/8] Reading FASTA and matching structure files", runlog)
        seqs, fasta_duplicates, fasta_headers = parse_fasta(args.fasta)
        structure_files = discover_structures(args.structures_dir)
        explicit_map = load_structure_id_map(args.structure_id_map)
        protein_ids = cluster_master["member_id"].tolist()
        matches, structure_match_df = match_structure_files(protein_ids, structure_files, explicit_map)
        structure_match_df.to_csv(intermediate_dir / "structure_id_matching.tsv", sep="\t", index=False)

        log(f"[3/8] Parsing {len(protein_ids)} structures and extracting residue-level pLDDT", runlog)
        rows = []
        profiles: Dict[str, List[Tuple[str, int, float]]] = {}
        errors = []
        for i, row in cluster_master.reset_index(drop=True).iterrows():
            pid = row["member_id"]
            if i % 100 == 0 or i + 1 == len(cluster_master):
                log(f"  processed {i + 1}/{len(cluster_master)}", runlog)
            sm = parse_structure(matches.get(pid))
            if sm.residue_profile is not None:
                profiles[pid] = sm.residue_profile
            fasta_length = len(seqs[pid]) if pid in seqs else np.nan
            modeled = sm.modeled_residue_count
            length_diff = fasta_length - modeled if np.isfinite(fasta_length) and np.isfinite(modeled) else np.nan
            coverage = modeled / fasta_length if np.isfinite(fasta_length) and fasta_length > 0 and np.isfinite(modeled) else np.nan
            note_parts = []
            if pid not in seqs and args.fasta is not None:
                note_parts.append("Missing FASTA entry.")
            if sm.notes:
                note_parts.append(sm.notes)
            record = {
                "protein_id": pid,
                "representative_id": row["representative_id"],
                "cluster_size": int(row["cluster_size"]),
                "cluster_group": row["cluster_group"],
                "fasta_length": fasta_length,
                "modeled_residue_count": modeled,
                "length_difference_fasta_minus_modeled": length_diff,
                "missing_or_unmodeled_residue_count": max(length_diff, 0) if np.isfinite(length_diff) else np.nan,
                "structure_coverage_fraction": coverage,
                "structure_file": sm.structure_file,
                "structure_format": sm.structure_format,
                "chain_count": sm.chain_count,
                "chain_residue_counts": sm.chain_residue_counts,
                "mean_plddt": sm.mean_plddt,
                "median_plddt": sm.median_plddt,
                "sd_plddt": sm.sd_plddt,
                "min_plddt": sm.min_plddt,
                "max_plddt": sm.max_plddt,
                "frac_plddt_lt50": sm.frac_plddt_lt50,
                "frac_plddt_lt70": sm.frac_plddt_lt70,
                "frac_plddt_ge70": sm.frac_plddt_ge70,
                "frac_plddt_ge90": sm.frac_plddt_ge90,
                "longest_segment_lt50": sm.longest_segment_lt50,
                "longest_segment_lt70": sm.longest_segment_lt70,
                "n_terminal_low_confidence_length": sm.n_terminal_low_confidence_length,
                "c_terminal_low_confidence_length": sm.c_terminal_low_confidence_length,
                "high_confidence_segment_count_ge20": sm.high_confidence_segment_count_ge20,
                "longest_high_confidence_segment": sm.longest_high_confidence_segment,
                "plddt_valid_0_100": sm.plddt_valid_0_100,
                "QC_flags": "",
                "primary_QC_category": "",
                "parsing_status": sm.parsing_status,
                "notes": " ".join(note_parts),
            }
            rows.append(record)
            if sm.parsing_status == "parse_failure":
                errors.append(f"{pid}\t{sm.structure_file}\t{sm.notes}")
        (logs_dir / "structure_parse_errors.tsv").write_text("protein_id\tstructure_file\terror\n" + "\n".join(errors) + ("\n" if errors else ""))
        master = pd.DataFrame(rows)

        log("[4/8] Deriving robust QC thresholds and exploratory categories", runlog)
        thresholds = robust_reference_thresholds(master)
        (tables_dir / "qc_thresholds.json").write_text(json.dumps(thresholds, indent=2, allow_nan=True))
        master = classify_qc(master, thresholds)
        master.to_csv(tables_dir / "singleton_length_plddt_master_table.tsv", sep="\t", index=False, na_rep="NA")

        log("[5/8] Building cluster, group-summary, and ID/QC audit tables", runlog)
        cluster_summary = {
            "total_unique_protein_count": int(cluster_master["member_id"].nunique()),
            "total_cluster_count": int(sizes.size),
            "singleton_cluster_count": int((sizes == 1).sum()),
            "singleton_protein_count": int((cluster_master["cluster_group"] == "singleton").sum()),
            "multi_member_cluster_count": int((sizes >= 2).sum()),
            "giant_cluster_representative": giant_rep,
            "giant_cluster_size": max_size,
            "giant_cluster_count": int((cluster_master["cluster_group"] == "giant_cluster").sum()),
            "other_multimember_protein_count": int((cluster_master["cluster_group"] == "other_multimember").sum()),
            "giant_cluster_tie_count": len(giant_candidates),
            "giant_cluster_tied_representatives": ";".join(giant_candidates),
        }
        pd.DataFrame([cluster_summary]).to_csv(tables_dir / "cluster_reconstruction_summary.tsv", sep="\t", index=False)
        group_summary = build_group_summary(master)
        group_summary.to_csv(tables_dir / "group_summary_table.tsv", sep="\t", index=False, na_rep="NA")

        cluster_id_set = set(protein_ids)
        fasta_id_set = set(seqs)
        matched_paths = set(str(p) for p in matches.values() if p is not None)
        unused_structures = [str(p) for p in structure_files if str(p) not in matched_paths]
        audit_df = build_missing_data_audit(
            cluster_master, fasta_id_set, structure_match_df, master,
            fasta_duplicates, sorted(fasta_id_set - cluster_id_set), unused_structures, cluster_audit
        )
        audit_df.to_csv(tables_dir / "missing_data_and_id_audit.tsv", sep="\t", index=False)
        audit_df["audit_type"].value_counts().rename_axis("audit_type").reset_index(name="count").to_csv(
            tables_dir / "missing_data_and_id_audit_summary.tsv", sep="\t", index=False
        )

        prioritized = build_prioritized_singletons(master)
        prioritized.to_csv(tables_dir / "singleton_prioritized_list.tsv", sep="\t", index=False, na_rep="NA")
        prioritized.loc[prioritized["direct_rescue_priority"] == "high", [
            "protein_id", "primary_QC_category", "fasta_length", "median_plddt",
            "frac_plddt_lt70", "structure_coverage_fraction", "QC_flags"
        ]].to_csv(tables_dir / "high_quality_unexplained_singletons_for_direct_rescue.tsv", sep="\t", index=False, na_rep="NA")

        log("[6/8] Running non-parametric statistics and effect-size estimation", runlog)
        stats_df = perform_statistics(master, args.bootstrap)
        stats_df.to_csv(tables_dir / "statistical_comparison_table.tsv", sep="\t", index=False, na_rep="NA")

        if not args.skip_plots:
            log("[7/8] Generating PDF and high-resolution PNG figures", runlog)
            create_plots(master, profiles, figures_dir)
        else:
            log("[7/8] Plot generation skipped", runlog)

        log("[8/8] Writing numerical summary and provisional Results text", runlog)
        write_summary_and_results(master, stats_df, cluster_summary, thresholds, outdir)
        log(f"Completed analysis: {outdir}", runlog)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL: {type(exc).__name__}: {exc}", file=sys.stderr)
        traceback.print_exc()
        raise
