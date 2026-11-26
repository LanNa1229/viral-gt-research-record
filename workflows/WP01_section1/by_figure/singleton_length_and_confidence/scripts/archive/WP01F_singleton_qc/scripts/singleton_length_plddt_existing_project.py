#!/usr/bin/env python3
"""Run singleton length/pLDDT QC directly from existing Foldseek/Cytoscape outputs.

This adapter reuses the analysis/statistics/plotting functions in
singleton_length_plddt_analysis.py, while adding support for:
  * cytoscape_complete_nodes.tsv as a complete cluster-membership source;
  * simultaneous cross-checking against accounting/cluster_membership.tsv;
  * separate original-sequence and Foldseek structure-derived FASTA files;
  * preservation of prior Cytoscape/network metadata in the master table;
  * validation against network_accounting.tsv.

The largest cluster is always determined from reconstructed membership counts.
No cluster ID or expected count is hard-coded.
"""
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

import singleton_length_plddt_analysis as core


NODE_REQUIRED = {"node_id", "cluster_representative"}
NODE_METADATA_COLUMNS = [
    "cluster_size",
    "cluster_class",
    "cluster_visibility_in_original_network",
    "node_role",
    "edge_to_representative",
    "is_cluster_representative",
    "loaded_in_original_direct_edge_network",
    "analysis_group",
]


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", dtype=str, keep_default_na=False)


def load_nodes_table(path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    df = read_table(path)
    df.columns = [str(c).strip() for c in df.columns]
    missing = NODE_REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"Nodes table is missing required columns: {sorted(missing)}")

    df["node_id"] = df["node_id"].map(core.normalize_id)
    df["cluster_representative"] = df["cluster_representative"].map(core.normalize_id)
    empty_mask = (df["node_id"] == "") | (df["cluster_representative"] == "")
    empty_rows = int(empty_mask.sum())
    df = df.loc[~empty_mask].copy()

    exact_duplicate_rows = int(df.duplicated().sum())
    duplicated_node_ids = sorted(
        df.loc[df["node_id"].duplicated(keep=False), "node_id"].unique().tolist()
    )
    conflicting = (
        df.groupby("node_id")["cluster_representative"].nunique().loc[lambda x: x > 1]
    )
    conflicting_ids = sorted(conflicting.index.tolist())

    # Keep one deterministic row per node for metadata; retain all anomalies in audit.
    metadata = (
        df.sort_values(["node_id", "cluster_representative"])
        .drop_duplicates("node_id", keep="first")
        .copy()
    )
    cluster = metadata[["cluster_representative", "node_id"]].copy()
    cluster.columns = ["representative_id", "member_id"]

    actual_sizes = cluster.groupby("representative_id")["member_id"].nunique()
    declared_size_mismatches = 0
    declared_class_mismatches = 0
    if "cluster_size" in metadata.columns:
        declared = pd.to_numeric(metadata["cluster_size"], errors="coerce")
        observed = metadata["cluster_representative"].map(actual_sizes)
        declared_size_mismatches = int((declared.notna() & declared.ne(observed)).sum())
    if "cluster_class" in metadata.columns:
        observed_class = np.where(
            metadata["cluster_representative"].map(actual_sizes).ge(2),
            "multi_member",
            "singleton",
        )
        declared_class_mismatches = int(
            (metadata["cluster_class"].astype(str).str.strip() != observed_class).sum()
        )

    audit = {
        "source": str(path),
        "raw_rows": int(len(df) + empty_rows),
        "valid_rows": int(len(df)),
        "unique_node_ids": int(metadata["node_id"].nunique()),
        "unique_cluster_representatives": int(cluster["representative_id"].nunique()),
        "empty_rows_removed": empty_rows,
        "exact_duplicate_rows": exact_duplicate_rows,
        "duplicated_node_id_count": len(duplicated_node_ids),
        "duplicated_node_ids": duplicated_node_ids,
        "conflicting_node_assignment_count": len(conflicting_ids),
        "conflicting_node_ids": conflicting_ids,
        "declared_cluster_size_mismatch_rows": declared_size_mismatches,
        "declared_cluster_class_mismatch_rows": declared_class_mismatches,
    }
    return cluster, metadata, audit


def membership_map(cluster: pd.DataFrame) -> Dict[str, str]:
    return dict(zip(cluster["member_id"], cluster["representative_id"]))


def compare_memberships(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    ma, mb = membership_map(a), membership_map(b)
    ids = sorted(set(ma) | set(mb))
    rows = []
    for pid in ids:
        ra, rb = ma.get(pid, ""), mb.get(pid, "")
        if ra != rb:
            rows.append({
                "protein_id": pid,
                "membership_representative": ra,
                "nodes_table_representative": rb,
                "status": (
                    "only_in_membership" if ra and not rb else
                    "only_in_nodes_table" if rb and not ra else
                    "representative_disagreement"
                ),
            })
    return pd.DataFrame(
        rows,
        columns=[
            "protein_id", "membership_representative",
            "nodes_table_representative", "status",
        ],
    )


def read_accounting(path: Optional[Path]) -> Dict[str, str]:
    if path is None or not path.exists():
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if not {"metric", "value"}.issubset(df.columns):
        raise ValueError(f"Accounting table must contain metric/value columns: {path}")
    return dict(zip(df["metric"], df["value"]))


def accounting_validation(
    accounting: Dict[str, str],
    cluster_master: pd.DataFrame,
    sizes: pd.Series,
    metadata: Optional[pd.DataFrame],
) -> pd.DataFrame:
    observed: Dict[str, object] = {
        "total_nodes": int(cluster_master["member_id"].nunique()),
        "total_clusters": int(sizes.size),
        "singleton_clusters": int((sizes == 1).sum()),
        "singleton_proteins": int((cluster_master["cluster_size"] == 1).sum()),
        "multi_member_clusters": int((sizes >= 2).sum()),
        "multi_member_proteins": int((cluster_master["cluster_size"] >= 2).sum()),
    }
    if metadata is not None and not metadata.empty:
        if "cluster_visibility_in_original_network" in metadata.columns:
            reps = metadata.drop_duplicates("cluster_representative")
            vis = reps["cluster_visibility_in_original_network"].astype(str)
            observed["visible_multi_member_clusters"] = int(
                ((pd.to_numeric(reps.get("cluster_size"), errors="coerce") >= 2) & vis.eq("visible")).sum()
            )
            observed["absent_multi_member_clusters"] = int(
                ((pd.to_numeric(reps.get("cluster_size"), errors="coerce") >= 2) & vis.eq("absent")).sum()
            )
        if "analysis_group" in metadata.columns:
            vc = metadata["analysis_group"].value_counts()
            observed["group_A_absent_multicluster_proteins"] = int(vc.get("A_absent_multicluster", 0))
            observed["group_B_singleton_proteins"] = int(vc.get("B_singleton", 0))
            observed["groups_A_and_B_combined_proteins"] = int(
                vc.get("A_absent_multicluster", 0) + vc.get("B_singleton", 0)
            )

    rows = []
    for metric, value in observed.items():
        expected = accounting.get(metric, "")
        status = "not_provided"
        if expected != "":
            try:
                status = "match" if float(expected) == float(value) else "mismatch"
            except ValueError:
                status = "match" if str(expected) == str(value) else "mismatch"
        rows.append({
            "metric": metric,
            "reported_value": expected,
            "recomputed_value": value,
            "status": status,
        })
    return pd.DataFrame(rows)


def add_audit(
    rows: List[Dict[str, str]], audit_type: str, protein_id: str, details: str
) -> None:
    rows.append({"audit_type": audit_type, "protein_id": protein_id, "details": details})


def safe_numeric(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan


def build_extended_audit(
    cluster_master: pd.DataFrame,
    metadata: Optional[pd.DataFrame],
    original_seqs: Dict[str, str],
    original_duplicates: Sequence[str],
    foldseek_seqs: Dict[str, str],
    foldseek_duplicates: Sequence[str],
    original_fasta_supplied: bool,
    foldseek_fasta_supplied: bool,
    structure_files: Sequence[Path],
    structure_match_df: pd.DataFrame,
    master: pd.DataFrame,
    membership_comparison: pd.DataFrame,
    validation_df: pd.DataFrame,
) -> pd.DataFrame:
    rows: List[Dict[str, str]] = []
    cluster_ids = set(cluster_master["member_id"])

    if original_fasta_supplied:
        for pid in sorted(cluster_ids - set(original_seqs)):
            add_audit(rows, "missing_original_fasta", pid, "Cluster member absent from original candidate FASTA.")
        for pid in sorted(set(original_seqs) - cluster_ids):
            add_audit(rows, "original_fasta_not_in_membership", pid, "Original FASTA entry is not a cluster member.")
        for pid in sorted(set(original_duplicates)):
            add_audit(rows, "duplicated_original_fasta_id", pid, "Duplicate original FASTA ID; first sequence retained.")

    if foldseek_fasta_supplied:
        for pid in sorted(cluster_ids - set(foldseek_seqs)):
            add_audit(rows, "missing_foldseek_fasta", pid, "Cluster member absent from Foldseek convert2fasta output.")
        for pid in sorted(set(foldseek_seqs) - cluster_ids):
            add_audit(rows, "foldseek_fasta_not_in_membership", pid, "Foldseek FASTA entry is not a cluster member.")
        for pid in sorted(set(foldseek_duplicates)):
            add_audit(rows, "duplicated_foldseek_fasta_id", pid, "Duplicate Foldseek FASTA ID; first sequence retained.")

    if metadata is not None:
        meta_ids = set(metadata["node_id"])
        for pid in sorted(cluster_ids - meta_ids):
            add_audit(rows, "membership_id_missing_from_nodes_table", pid, "No Cytoscape node metadata was found.")
        for pid in sorted(meta_ids - cluster_ids):
            add_audit(rows, "nodes_table_id_missing_from_membership", pid, "Node is absent from canonical membership.")

    for _, r in membership_comparison.iterrows():
        add_audit(rows, f"membership_crosscheck_{r['status']}", r["protein_id"],
                  f"membership={r['membership_representative']}; nodes={r['nodes_table_representative']}")

    for _, r in structure_match_df.iterrows():
        if r["status"] == "missing":
            add_audit(rows, "missing_structure", r["protein_id"], "No structure filename matched.")
        elif r["status"] == "ambiguous":
            add_audit(rows, "ambiguous_structure_match", r["protein_id"], r["candidate_structures"])

    # Treat every candidate path associated with a cluster ID as accounted for,
    # including redundant nested copies such as selected_hits/ID.pdb.
    accounted_paths = set(
        structure_match_df.loc[structure_match_df["matched_structure"] != "", "matched_structure"]
    )
    for value in structure_match_df.get("candidate_structures", pd.Series(dtype=str)).astype(str):
        for candidate in value.split(";"):
            if candidate:
                accounted_paths.add(candidate)
    for path in sorted(str(p) for p in structure_files if str(p) not in accounted_paths):
        add_audit(rows, "structure_not_in_cluster_membership", "", path)

    nonempty = structure_match_df.loc[
        structure_match_df["matched_structure"] != "", "matched_structure"
    ]
    for path, count in nonempty[nonempty.duplicated(keep=False)].value_counts().items():
        ids = structure_match_df.loc[
            structure_match_df["matched_structure"] == path, "protein_id"
        ].tolist()
        add_audit(rows, "structure_file_matched_to_multiple_ids", ";".join(ids), f"{path} matched to {count} IDs.")

    for _, r in master.iterrows():
        pid = r["protein_id"]
        if r["parsing_status"] == "parse_failure":
            add_audit(rows, "structure_parse_failure", pid, str(r.get("notes", "")))
        elif r["parsing_status"] == "suspicious_plddt_range":
            add_audit(rows, "suspicious_plddt_range", pid, str(r.get("notes", "")))
        if pd.notna(r.get("chain_count")) and safe_numeric(r.get("chain_count")) > 1:
            add_audit(rows, "multiple_chains", pid, str(r.get("chain_residue_counts", "")))
        d = safe_numeric(r.get("length_difference_original_minus_modeled"))
        if np.isfinite(d) and d < 0:
            add_audit(rows, "modeled_count_exceeds_original_fasta_length", pid, str(d))
        fd = safe_numeric(r.get("foldseek_sequence_minus_modeled_residue_count"))
        if np.isfinite(fd) and fd != 0:
            add_audit(rows, "foldseek_sequence_modeled_count_difference", pid, str(fd))

    for _, r in validation_df.loc[validation_df["status"] == "mismatch"].iterrows():
        add_audit(rows, "network_accounting_mismatch", "",
                  f"{r['metric']}: reported={r['reported_value']}; recomputed={r['recomputed_value']}")

    return pd.DataFrame(rows, columns=["audit_type", "protein_id", "details"])


def stats_or_empty(master: pd.DataFrame, bootstrap: int) -> pd.DataFrame:
    has_any = any(pd.to_numeric(master.get(metric), errors="coerce").notna().any()
                  for metric in core.METRICS_FOR_STATS if metric in master.columns)
    if not has_any:
        return pd.DataFrame(columns=[
            "metric", "groups_compared", "test", "statistic", "raw_p_value",
            "adjusted_p_value", "effect_size_name", "effect_size",
            "effect_ci_lower", "effect_ci_upper", "median_first", "median_second",
            "median_difference_first_minus_second", "hodges_lehmann",
            "hodges_lehmann_method", "median_difference_ci_lower",
            "median_difference_ci_upper", "interpretation",
        ])
    return core.perform_statistics(master, bootstrap)


def main() -> int:
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--network-dir", type=Path, default=None,
                   help="Directory containing cytoscape_complete_nodes.tsv and optionally network_accounting.tsv.")
    p.add_argument("--nodes-table", type=Path, default=None)
    p.add_argument("--network-accounting", type=Path, default=None)
    p.add_argument("--cluster-tsv", type=Path, default=None,
                   help="Canonical Foldseek representative-member TSV. If supplied with a nodes table, both are cross-checked.")
    p.add_argument("--rep-column", default=None)
    p.add_argument("--member-column", default=None)
    p.add_argument("--original-fasta", type=Path, default=None,
                   help="Original candidate protein FASTA, preferred for true sequence length and structure coverage.")
    p.add_argument("--foldseek-fasta", type=Path, default=None,
                   help="Foldseek convert2fasta output. Treated as structure-derived sequence, not automatically as the original full ORF.")
    p.add_argument("--structures-dir", type=Path, default=None)
    p.add_argument("--structure-id-map", type=Path, default=None)
    p.add_argument("--outdir", type=Path, default=Path("singleton_length_plddt_analysis_v2"))
    p.add_argument("--bootstrap", type=int, default=500)
    p.add_argument("--skip-plots", action="store_true")
    p.add_argument("--allow-membership-mismatch", action="store_true")
    args = p.parse_args()

    if args.network_dir is not None:
        if args.nodes_table is None:
            candidate = args.network_dir / "cytoscape_complete_nodes.tsv"
            if candidate.exists():
                args.nodes_table = candidate
        if args.network_accounting is None:
            candidate = args.network_dir / "network_accounting.tsv"
            if candidate.exists():
                args.network_accounting = candidate

    if args.cluster_tsv is None and args.nodes_table is None:
        p.error("Provide --cluster-tsv, --nodes-table, or --network-dir containing cytoscape_complete_nodes.tsv.")

    outdir = args.outdir.resolve()
    tables_dir = outdir / "tables"
    figures_dir = outdir / "figures"
    logs_dir = outdir / "logs"
    intermediate_dir = outdir / "intermediate"
    for d in [outdir, tables_dir, figures_dir, logs_dir, intermediate_dir]:
        d.mkdir(parents=True, exist_ok=True)

    with (logs_dir / "run.log").open("w") as runlog:
        core.log("[1/9] Reconstructing and cross-checking cluster membership", runlog)
        membership_cluster = None
        membership_audit: Dict[str, object] = {}
        if args.cluster_tsv is not None:
            membership_cluster, membership_audit = core.load_cluster_membership(
                args.cluster_tsv, args.rep_column, args.member_column
            )

        nodes_cluster = None
        node_metadata = None
        nodes_audit: Dict[str, object] = {}
        if args.nodes_table is not None:
            nodes_cluster, node_metadata, nodes_audit = load_nodes_table(args.nodes_table)

        comparison = pd.DataFrame(columns=[
            "protein_id", "membership_representative", "nodes_table_representative", "status"
        ])
        if membership_cluster is not None and nodes_cluster is not None:
            comparison = compare_memberships(membership_cluster, nodes_cluster)
            comparison.to_csv(intermediate_dir / "membership_vs_nodes_table.tsv", sep="\t", index=False)
            if not comparison.empty and not args.allow_membership_mismatch:
                raise ValueError(
                    f"Membership TSV and nodes table disagree for {len(comparison)} proteins. "
                    "Inspect intermediate/membership_vs_nodes_table.tsv or rerun with --allow-membership-mismatch."
                )

        canonical = membership_cluster if membership_cluster is not None else nodes_cluster
        assert canonical is not None
        canonical = canonical.drop_duplicates(["representative_id", "member_id"]).copy()
        canonical.to_csv(intermediate_dir / "cluster_membership_deduplicated.tsv", sep="\t", index=False)
        (logs_dir / "cluster_input_audit.json").write_text(json.dumps({
            "membership_tsv": membership_audit,
            "nodes_table": nodes_audit,
            "membership_nodes_disagreement_count": int(len(comparison)),
            "canonical_source": "cluster_tsv" if membership_cluster is not None else "cytoscape_nodes_table",
        }, indent=2))

        # Deterministic conflict handling is retained only for auditability; a valid Foldseek table should be one-to-one.
        cluster_master = canonical.sort_values(["member_id", "representative_id"]).drop_duplicates("member_id", keep="first")
        sizes = cluster_master.groupby("representative_id")["member_id"].nunique().sort_values(ascending=False)
        if sizes.empty:
            raise ValueError("No valid cluster membership rows remain.")
        max_size = int(sizes.max())
        giant_candidates = sorted(sizes[sizes == max_size].index.tolist())
        giant_rep = giant_candidates[0]
        cluster_master["cluster_size"] = cluster_master["representative_id"].map(sizes).astype(int)
        cluster_master["cluster_group"] = np.where(
            cluster_master["representative_id"] == giant_rep,
            "giant_cluster",
            np.where(cluster_master["cluster_size"] >= 2, "other_multimember", "singleton"),
        )

        if node_metadata is not None:
            meta = node_metadata.copy().rename(columns={
                "node_id": "member_id",
                "cluster_representative": "nodes_table_cluster_representative",
                "cluster_size": "nodes_table_cluster_size",
                "cluster_class": "nodes_table_cluster_class",
            })
            keep = ["member_id"] + [c for c in meta.columns if c != "member_id"]
            cluster_master = cluster_master.merge(meta[keep], on="member_id", how="left", validate="one_to_one")

        core.log("[2/9] Reading original and Foldseek-derived sequence FASTA files", runlog)
        original_seqs, original_duplicates, _ = core.parse_fasta(args.original_fasta)
        foldseek_seqs, foldseek_duplicates, _ = core.parse_fasta(args.foldseek_fasta)

        core.log("[3/9] Matching existing predicted structures", runlog)
        structure_files = core.discover_structures(args.structures_dir)
        explicit_map = core.load_structure_id_map(args.structure_id_map)
        protein_ids = cluster_master["member_id"].tolist()
        matches, structure_match_df = core.match_structure_files(protein_ids, structure_files, explicit_map)
        structure_match_df.to_csv(intermediate_dir / "structure_id_matching.tsv", sep="\t", index=False)

        core.log(f"[4/9] Parsing {len(protein_ids)} structures and residue-level confidence", runlog)
        rows = []
        profiles = {}
        errors = []
        for i, row in cluster_master.reset_index(drop=True).iterrows():
            pid = row["member_id"]
            if i % 100 == 0 or i + 1 == len(cluster_master):
                core.log(f"  processed {i + 1}/{len(cluster_master)}", runlog)
            sm = core.parse_structure(matches.get(pid))
            if sm.residue_profile is not None:
                profiles[pid] = sm.residue_profile

            original_len = len(original_seqs[pid]) if pid in original_seqs else np.nan
            foldseek_len = len(foldseek_seqs[pid]) if pid in foldseek_seqs else np.nan
            modeled = sm.modeled_residue_count
            if np.isfinite(original_len):
                analysis_len = original_len
                length_source = "original_candidate_fasta"
            elif np.isfinite(foldseek_len):
                analysis_len = foldseek_len
                length_source = "foldseek_structure_derived_fasta"
            elif np.isfinite(modeled):
                analysis_len = modeled
                length_source = "modeled_residue_count"
            else:
                analysis_len = np.nan
                length_source = "missing"

            original_diff = (
                original_len - modeled
                if np.isfinite(original_len) and np.isfinite(modeled) else np.nan
            )
            original_coverage = (
                modeled / original_len
                if np.isfinite(original_len) and original_len > 0 and np.isfinite(modeled) else np.nan
            )
            foldseek_diff = (
                foldseek_len - modeled
                if np.isfinite(foldseek_len) and np.isfinite(modeled) else np.nan
            )
            foldseek_consistency = (
                modeled / foldseek_len
                if np.isfinite(foldseek_len) and foldseek_len > 0 and np.isfinite(modeled) else np.nan
            )

            note_parts = []
            if args.original_fasta is not None and pid not in original_seqs:
                note_parts.append("Missing original FASTA entry.")
            if args.foldseek_fasta is not None and pid not in foldseek_seqs:
                note_parts.append("Missing Foldseek-derived FASTA entry.")
            if sm.notes:
                note_parts.append(sm.notes)

            record = {
                "protein_id": pid,
                "representative_id": row["representative_id"],
                "cluster_size": int(row["cluster_size"]),
                "cluster_group": row["cluster_group"],
                # Backward-compatible field used by the core plots/statistics.
                "fasta_length": analysis_len,
                "analysis_sequence_length": analysis_len,
                "analysis_sequence_length_source": length_source,
                "original_fasta_length": original_len,
                "foldseek_structure_sequence_length": foldseek_len,
                "modeled_residue_count": modeled,
                "length_difference_fasta_minus_modeled": original_diff,
                "length_difference_original_minus_modeled": original_diff,
                "missing_or_unmodeled_residue_count": (
                    max(original_diff, 0) if np.isfinite(original_diff) else np.nan
                ),
                "structure_coverage_fraction": original_coverage,
                "coverage_interpretation": (
                    "original_sequence_coverage" if np.isfinite(original_coverage)
                    else "not_computable_without_original_fasta"
                ),
                "foldseek_sequence_minus_modeled_residue_count": foldseek_diff,
                "foldseek_model_consistency_fraction": foldseek_consistency,
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
            # Preserve the old network classifications as provenance, not as replacements for reconstructed groups.
            for col in cluster_master.columns:
                if col not in {"member_id", "representative_id", "cluster_size", "cluster_group"}:
                    record[col] = row.get(col, "")
            rows.append(record)
            if sm.parsing_status == "parse_failure":
                errors.append(f"{pid}\t{sm.structure_file}\t{sm.notes}")

        (logs_dir / "structure_parse_errors.tsv").write_text(
            "protein_id\tstructure_file\terror\n" + "\n".join(errors) + ("\n" if errors else "")
        )
        master = pd.DataFrame(rows)

        core.log("[5/9] Deriving robust reference thresholds and QC categories", runlog)
        thresholds = core.robust_reference_thresholds(master)
        (tables_dir / "qc_thresholds.json").write_text(json.dumps(thresholds, indent=2, allow_nan=True))
        master = core.classify_qc(master, thresholds)
        master.to_csv(tables_dir / "singleton_length_plddt_master_table.tsv", sep="\t", index=False, na_rep="NA")

        core.log("[6/9] Writing reconstruction, accounting, and ID-audit tables", runlog)
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
            "canonical_membership_source": "cluster_tsv" if membership_cluster is not None else "cytoscape_nodes_table",
            "sequence_length_priority": "original FASTA > Foldseek structure-derived FASTA > modeled residue count",
        }
        pd.DataFrame([cluster_summary]).to_csv(
            tables_dir / "cluster_reconstruction_summary.tsv", sep="\t", index=False
        )
        core.build_group_summary(master).to_csv(
            tables_dir / "group_summary_table.tsv", sep="\t", index=False, na_rep="NA"
        )

        accounting = read_accounting(args.network_accounting)
        validation_df = accounting_validation(accounting, cluster_master, sizes, node_metadata)
        validation_df.to_csv(tables_dir / "network_accounting_validation.tsv", sep="\t", index=False)

        audit_df = build_extended_audit(
            cluster_master, node_metadata,
            original_seqs, original_duplicates,
            foldseek_seqs, foldseek_duplicates,
            args.original_fasta is not None, args.foldseek_fasta is not None,
            structure_files, structure_match_df, master,
            comparison, validation_df,
        )
        audit_df.to_csv(tables_dir / "missing_data_and_id_audit.tsv", sep="\t", index=False)
        audit_df["audit_type"].value_counts().rename_axis("audit_type").reset_index(name="count").to_csv(
            tables_dir / "missing_data_and_id_audit_summary.tsv", sep="\t", index=False
        )

        prioritized = core.build_prioritized_singletons(master)
        prioritized.to_csv(tables_dir / "singleton_prioritized_list.tsv", sep="\t", index=False, na_rep="NA")
        rescue_cols = [
            "protein_id", "primary_QC_category", "relative_length_class",
            "confidence_class", "combined_length_confidence_class",
            "analysis_sequence_length", "analysis_sequence_length_source",
            "median_plddt", "frac_plddt_lt70", "longest_high_confidence_segment",
            "structure_coverage_fraction", "foldseek_model_consistency_fraction", "QC_flags",
        ]
        prioritized.loc[prioritized["direct_rescue_priority"] == "high", rescue_cols].to_csv(
            tables_dir / "high_quality_unexplained_singletons_for_direct_rescue.tsv",
            sep="\t", index=False, na_rep="NA"
        )

        core.log("[7/9] Running non-parametric statistics and effect-size estimation", runlog)
        stats_df = stats_or_empty(master, args.bootstrap)
        stats_df.to_csv(tables_dir / "statistical_comparison_table.tsv", sep="\t", index=False, na_rep="NA")

        analyzable = pd.to_numeric(master["fasta_length"], errors="coerce").notna().any()
        confidence_available = pd.to_numeric(master["median_plddt"], errors="coerce").notna().any()
        if not args.skip_plots and analyzable and confidence_available:
            core.log("[8/9] Generating PDF and high-resolution PNG figures", runlog)
            core.create_plots(master, profiles, figures_dir)
        else:
            core.log("[8/9] Plot generation skipped because requested or confidence data are unavailable", runlog)

        core.log("[9/9] Writing numerical summary and provisional Results text", runlog)
        core.write_summary_and_results(master, stats_df, cluster_summary, thresholds, outdir)
        provenance = [
            "# Input provenance and interpretation",
            "",
            f"- Canonical membership: {'cluster TSV' if membership_cluster is not None else 'Cytoscape complete nodes table'}",
            f"- Original candidate FASTA: {args.original_fasta or 'not supplied'}",
            f"- Foldseek structure-derived FASTA: {args.foldseek_fasta or 'not supplied'}",
            f"- Predicted-structure directory: {args.structures_dir or 'not supplied'}",
            "",
            "`database_entries.fasta` generated by `foldseek convert2fasta` is treated as a structure-derived sequence source. "
            "It can support length stratification, but true missing/unmodeled-residue counts and original-sequence coverage are only assigned when an original candidate FASTA is supplied.",
            "",
            "Prior Cytoscape labels such as `analysis_group` and network visibility are retained as provenance columns. "
            "The scientific groups `giant_cluster`, `other_multimember`, and `singleton` are reconstructed independently from actual cluster sizes.",
        ]
        (outdir / "input_provenance.md").write_text("\n".join(provenance) + "\n")
        core.log(f"Completed analysis: {outdir}", runlog)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL: {type(exc).__name__}: {exc}")
        traceback.print_exc()
        raise
