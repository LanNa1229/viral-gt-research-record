#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


M8_COLUMNS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Build a complete Cytoscape network from Foldseek clustering."
    )
    p.add_argument("--membership", required=True, type=Path)
    p.add_argument("--m8", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def read_membership(path):
    cluster_members = defaultdict(list)
    member_to_rep = {}
    order = []

    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_no, row in enumerate(reader, 1):
            if len(row) < 2:
                raise ValueError(f"{path}:{line_no}: expected 2 columns")
            rep, member = row[0].strip(), row[1].strip()
            if member in member_to_rep and member_to_rep[member] != rep:
                raise ValueError(f"{member} assigned to multiple representatives")
            if member not in member_to_rep:
                member_to_rep[member] = rep
                cluster_members[rep].append(member)
                order.append(member)

    return dict(cluster_members), member_to_rep, order


def as_float(value):
    try:
        return float(value)
    except ValueError:
        return float("-inf")


def read_direct_m8(path, member_to_rep):
    direct = {}
    unmatched = []
    nonself_rows = 0

    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_no, row in enumerate(reader, 1):
            if len(row) < 12:
                raise ValueError(f"{path}:{line_no}: expected >=12 columns")
            q, t = row[0].strip(), row[1].strip()
            if q == t:
                continue
            nonself_rows += 1

            if member_to_rep.get(t) == q:
                rep, member = q, t
            elif member_to_rep.get(q) == t:
                rep, member = t, q
            else:
                unmatched.append(row[:12])
                continue

            record = {name: row[i].strip() for i, name in enumerate(M8_COLUMNS)}
            record["raw_query"] = q
            record["raw_target"] = t

            key = (rep, member)
            if key not in direct or as_float(record["bits"]) > as_float(direct[key]["bits"]):
                direct[key] = record

    return direct, unmatched, nonself_rows


def write_tsv(path, fieldnames, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t",
            lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()

    for path in (args.membership, args.m8):
        if not path.is_file():
            raise FileNotFoundError(path)

    args.outdir.mkdir(parents=True, exist_ok=True)

    cluster_members, member_to_rep, order = read_membership(args.membership)
    direct, unmatched, raw_nonself_rows = read_direct_m8(args.m8, member_to_rep)

    cluster_sizes = {rep: len(members) for rep, members in cluster_members.items()}
    multi_reps = {rep for rep, size in cluster_sizes.items() if size >= 2}
    singleton_reps = {rep for rep, size in cluster_sizes.items() if size == 1}
    visible_reps = {rep for rep, member in direct if rep in multi_reps}
    absent_reps = multi_reps - visible_reps

    direct_nodes = set()
    for rep, member in direct:
        direct_nodes.update((rep, member))

    edge_rows = []

    for rep in sorted(multi_reps):
        for member in cluster_members[rep]:
            if member == rep:
                continue

            rec = direct.get((rep, member))
            if rec:
                edge_type = "direct_alignment"
                evidence = "alignment_present_in_m8"
                metrics = {name: rec.get(name, "") for name in M8_COLUMNS[2:]}
                raw_query = rec["raw_query"]
                raw_target = rec["raw_target"]
                is_direct = 1
            else:
                edge_type = "cluster_membership"
                evidence = "same_foldseek_cluster_no_direct_rep_alignment_in_m8"
                metrics = {name: "" for name in M8_COLUMNS[2:]}
                raw_query = ""
                raw_target = ""
                is_direct = 0

            row = {
                "edge_id": f"{rep}__{edge_type}__{member}",
                "source": rep,
                "target": member,
                "interaction": edge_type,
                "edge_type": edge_type,
                "evidence": evidence,
                "is_direct_alignment": is_direct,
                "cluster_representative": rep,
                "cluster_size": cluster_sizes[rep],
                "cluster_visibility_in_original_network": (
                    "visible" if rep in visible_reps else "absent"
                ),
                "raw_m8_query": raw_query,
                "raw_m8_target": raw_target,
            }
            row.update(metrics)
            edge_rows.append(row)

    edge_fields = [
        "edge_id", "source", "target", "interaction", "edge_type",
        "evidence", "is_direct_alignment", "cluster_representative",
        "cluster_size", "cluster_visibility_in_original_network",
        "raw_m8_query", "raw_m8_target", "fident", "alnlen", "mismatch",
        "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits",
    ]
    write_tsv(args.outdir / "cytoscape_complete_edges.tsv", edge_fields, edge_rows)

    node_rows = []

    for node_id in order:
        rep = member_to_rep[node_id]
        size = cluster_sizes[rep]

        if size == 1:
            cluster_class = "singleton"
            visibility = "singleton"
            node_role = "singleton"
            edge_to_rep = "none"
            analysis_group = "B_singleton"
        else:
            cluster_class = "multi_member"
            visibility = "visible" if rep in visible_reps else "absent"
            analysis_group = (
                "A_absent_multicluster"
                if rep in absent_reps
                else "C_visible_multicluster"
            )

            if node_id == rep:
                node_role = "representative"
                edge_to_rep = "self"
            elif (rep, node_id) in direct:
                node_role = "direct_member"
                edge_to_rep = "direct_alignment"
            else:
                node_role = "indirect_member"
                edge_to_rep = "cluster_membership"

        node_rows.append({
            "node_id": node_id,
            "cluster_representative": rep,
            "cluster_size": size,
            "cluster_class": cluster_class,
            "cluster_visibility_in_original_network": visibility,
            "node_role": node_role,
            "edge_to_representative": edge_to_rep,
            "is_cluster_representative": int(node_id == rep),
            "loaded_in_original_direct_edge_network": int(node_id in direct_nodes),
            "analysis_group": analysis_group,
        })

    node_fields = [
        "node_id", "cluster_representative", "cluster_size", "cluster_class",
        "cluster_visibility_in_original_network", "node_role",
        "edge_to_representative", "is_cluster_representative",
        "loaded_in_original_direct_edge_network", "analysis_group",
    ]
    write_tsv(args.outdir / "cytoscape_complete_nodes.tsv", node_fields, node_rows)

    style_rows = [
        {
            "edge_type": "direct_alignment",
            "biological_meaning": "Direct representative-member alignment present in m8",
            "recommended_line_style": "SOLID",
            "recommended_width": "2.0",
            "recommended_color": "#2166AC",
        },
        {
            "edge_type": "cluster_membership",
            "biological_meaning": "Same Foldseek cluster; no direct final-representative alignment in m8",
            "recommended_line_style": "LONG_DASH",
            "recommended_width": "1.0",
            "recommended_color": "#969696",
        },
    ]
    write_tsv(
        args.outdir / "cytoscape_edge_style_legend.tsv",
        ["edge_type", "biological_meaning", "recommended_line_style",
         "recommended_width", "recommended_color"],
        style_rows,
    )

    with (args.outdir / "complete_multimember_network.sif").open("w") as handle:
        for row in edge_rows:
            handle.write(f"{row['source']}\t{row['edge_type']}\t{row['target']}\n")

    # Optional SIF containing all 1,727 nodes. SIF permits one-column lines
    # for isolated nodes, so singleton proteins can be imported as isolated nodes.
    with (args.outdir / "complete_all_nodes_network.sif").open("w") as handle:
        for row in edge_rows:
            handle.write(f"{row['source']}\t{row['edge_type']}\t{row['target']}\n")
        for singleton in sorted(singleton_reps):
            handle.write(f"{singleton}\n")

    group_a = [r for r in node_rows if r["analysis_group"] == "A_absent_multicluster"]
    group_a_reps = [r for r in group_a if r["is_cluster_representative"] == 1]
    group_b = [r for r in node_rows if r["analysis_group"] == "B_singleton"]

    write_tsv(args.outdir / "group_A_absent_multicluster_proteins.tsv",
              node_fields, group_a)
    write_tsv(args.outdir / "group_A_absent_multicluster_representatives.tsv",
              node_fields, group_a_reps)
    write_tsv(args.outdir / "group_B_singletons.tsv",
              node_fields, group_b)
    write_tsv(args.outdir / "groups_A_and_B_combined.tsv",
              node_fields, group_a + group_b)

    with (args.outdir / "unmatched_m8_nonself.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(M8_COLUMNS)
        writer.writerows(unmatched)

    direct_edges = sum(r["edge_type"] == "direct_alignment" for r in edge_rows)
    membership_edges = sum(r["edge_type"] == "cluster_membership" for r in edge_rows)
    edge_nodes = {r["source"] for r in edge_rows} | {r["target"] for r in edge_rows}

    accounting = [
        ("total_nodes", len(member_to_rep)),
        ("total_clusters", len(cluster_members)),
        ("singleton_clusters", len(singleton_reps)),
        ("singleton_proteins", len(singleton_reps)),
        ("multi_member_clusters", len(multi_reps)),
        ("multi_member_proteins", sum(cluster_sizes[r] for r in multi_reps)),
        ("visible_multi_member_clusters", len(visible_reps)),
        ("absent_multi_member_clusters", len(absent_reps)),
        ("original_direct_edge_network_nodes", len(direct_nodes)),
        ("complete_multimember_network_nodes", len(edge_nodes)),
        ("direct_alignment_edges", direct_edges),
        ("cluster_membership_edges", membership_edges),
        ("total_complete_edges", len(edge_rows)),
        ("group_A_absent_multicluster_proteins", len(group_a)),
        ("group_A_absent_multicluster_representatives", len(group_a_reps)),
        ("group_B_singleton_proteins", len(group_b)),
        ("groups_A_and_B_combined_proteins", len(group_a) + len(group_b)),
        ("raw_nonself_m8_rows", raw_nonself_rows),
        ("unmatched_nonself_m8_rows", len(unmatched)),
    ]

    with (args.outdir / "network_accounting.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(accounting)

    expected_edges = sum(cluster_sizes[r] for r in multi_reps) - len(multi_reps)
    expected_nodes = sum(cluster_sizes[r] for r in multi_reps)

    problems = []
    if len(edge_rows) != expected_edges:
        problems.append(f"edge count {len(edge_rows)} != expected {expected_edges}")
    if len(edge_nodes) != expected_nodes:
        problems.append(f"node count {len(edge_nodes)} != expected {expected_nodes}")

    print("=" * 68)
    print("Foldseek complete-network accounting")
    print("=" * 68)
    for key, value in accounting:
        print(f"{key:46s} {value}")

    if problems:
        print("\nWARNING:", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        return 2

    print("\nAll integrity checks passed.")
    print(f"Output directory: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
