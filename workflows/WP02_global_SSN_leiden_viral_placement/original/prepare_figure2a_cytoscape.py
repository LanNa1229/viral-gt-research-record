#!/usr/bin/env python3
"""Prepare compact Cytoscape inputs for Viral GT manuscript Figure 2A.

The analytical graph is not altered. In addition to the full edge table, this
script creates a display-only backbone consisting of:

1. a maximum-weight spanning forest, which preserves connectivity within every
   non-isolated connected component; and
2. the strongest K incident edges for every node.

The backbone is intended only to reduce edge overplotting in the publication
view. All community assignments and statistics remain those calculated from
the full network.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


PALETTE = [
    "#4477AA",  # blue
    "#EE6677",  # rose
    "#228833",  # green
    "#CCBB44",  # gold
    "#66CCEE",  # cyan
    "#AA3377",  # purple
    "#EE8866",  # orange
    "#44AA99",  # teal
    "#999933",  # olive
    "#882255",  # wine
    "#6699CC",  # light blue
    "#AA4466",  # mauve
    "#117733",  # dark green
    "#332288",  # indigo
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare Figure 2A Cytoscape import and styling tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--edges", required=True)
    parser.add_argument("--communities", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--node-id-column", default="name")
    parser.add_argument("--source-column", default="source")
    parser.add_argument("--target-column", default="target")
    parser.add_argument("--community-column", default="community")
    parser.add_argument("--weight-column", default="weight_global_min_tm")
    parser.add_argument("--top-edges-per-node", type=int, default=3)
    parser.add_argument("--top-cazy-families", type=int, default=14)
    parser.add_argument("--major-community-min-n", type=int, default=40)
    return parser.parse_args()


def read_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", keep_default_na=False)


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_bool(value: object) -> bool:
    return str(value).strip().casefold() in {
        "1",
        "true",
        "t",
        "yes",
        "y",
        "virus",
        "viral",
    }


def normalize_cazy(value: object) -> str:
    text = str(value).strip().upper()
    match = re.match(r"^(GT\d+)(?:$|[-_].*)", text)
    return match.group(1) if match else "Unassigned"


class DisjointSet:
    def __init__(self, items: list[str]) -> None:
        self.parent = {item: item for item in items}
        self.rank = {item: 0 for item in items}

    def find(self, item: str) -> str:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: str, right: str) -> bool:
        root_l = self.find(left)
        root_r = self.find(right)
        if root_l == root_r:
            return False
        if self.rank[root_l] < self.rank[root_r]:
            root_l, root_r = root_r, root_l
        self.parent[root_r] = root_l
        if self.rank[root_l] == self.rank[root_r]:
            self.rank[root_l] += 1
        return True


def build_backbone(
    edges: pd.DataFrame,
    node_ids: list[str],
    source: str,
    target: str,
    weight: str,
    top_k: int,
) -> tuple[pd.DataFrame, dict[str, int]]:
    work = edges.copy()
    work[source] = work[source].astype(str)
    work[target] = work[target].astype(str)
    work[weight] = pd.to_numeric(work[weight], errors="raise")
    work["edge_key_internal"] = work.apply(
        lambda row: "||".join(sorted((row[source], row[target]))), axis=1
    )
    work = work.sort_values(
        [weight, source, target], ascending=[False, True, True]
    ).drop_duplicates("edge_key_internal", keep="first")

    # Maximum-weight spanning forest.
    dsu = DisjointSet(node_ids)
    forest_keys: set[str] = set()
    for row in work.itertuples(index=False):
        left = str(getattr(row, source))
        right = str(getattr(row, target))
        if dsu.union(left, right):
            forest_keys.add(getattr(row, "edge_key_internal"))

    # Strongest K incident edges for every node.
    strongest_keys: set[str] = set()
    incident_count: dict[str, int] = defaultdict(int)
    for row in work.itertuples(index=False):
        left = str(getattr(row, source))
        right = str(getattr(row, target))
        key = getattr(row, "edge_key_internal")
        take_left = incident_count[left] < top_k
        take_right = incident_count[right] < top_k
        if take_left or take_right:
            strongest_keys.add(key)
        if take_left:
            incident_count[left] += 1
        if take_right:
            incident_count[right] += 1

    chosen = forest_keys | strongest_keys
    backbone = work.loc[work["edge_key_internal"].isin(chosen)].copy()
    backbone["backbone_msf"] = backbone["edge_key_internal"].isin(forest_keys)
    backbone["backbone_topk"] = backbone["edge_key_internal"].isin(strongest_keys)
    backbone["interaction"] = "structural_similarity"
    backbone["edge_display_width"] = (
        0.25
        + 1.25
        * (backbone[weight] - backbone[weight].min())
        / max(backbone[weight].max() - backbone[weight].min(), 1e-12)
    )
    backbone["edge_display_transparency"] = (
        35
        + 90
        * (backbone[weight] - backbone[weight].min())
        / max(backbone[weight].max() - backbone[weight].min(), 1e-12)
    ).round().astype(int)
    backbone = backbone.drop(columns=["edge_key_internal"])

    stats = {
        "full_unique_edges": int(len(work)),
        "maximum_spanning_forest_edges": int(len(forest_keys)),
        "strongest_incident_edge_union": int(len(strongest_keys)),
        "display_backbone_edges": int(len(backbone)),
        "top_edges_per_node": int(top_k),
    }
    return backbone, stats


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    nodes = read_table(args.nodes)
    edges = read_table(args.edges)
    communities = read_table(args.communities)

    required_nodes = {
        args.node_id_column,
        args.community_column,
        "CAZy_family",
        "is_virus",
        "is_giant_component",
        "weighted_degree",
    }
    missing_nodes = sorted(required_nodes - set(nodes.columns))
    if missing_nodes:
        raise ValueError(f"Node table missing columns: {missing_nodes}")
    required_edges = {args.source_column, args.target_column, args.weight_column}
    missing_edges = sorted(required_edges - set(edges.columns))
    if missing_edges:
        raise ValueError(f"Edge table missing columns: {missing_edges}")

    nodes[args.node_id_column] = nodes[args.node_id_column].astype(str).str.strip()
    edges[args.source_column] = edges[args.source_column].astype(str).str.strip()
    edges[args.target_column] = edges[args.target_column].astype(str).str.strip()
    node_ids = nodes[args.node_id_column].tolist()
    node_set = set(node_ids)
    if len(node_set) != len(node_ids):
        raise ValueError("Node IDs are not unique")
    unknown_edge_nodes = (
        set(edges[args.source_column]) | set(edges[args.target_column])
    ) - node_set
    if unknown_edge_nodes:
        raise ValueError(
            f"Edges reference {len(unknown_edge_nodes)} unknown nodes; "
            f"examples: {sorted(unknown_edge_nodes)[:10]}"
        )

    nodes["cazy_base"] = nodes["CAZy_family"].map(normalize_cazy)
    cazy_counts = (
        nodes.loc[nodes["cazy_base"].ne("Unassigned"), "cazy_base"]
        .value_counts()
        .sort_values(ascending=False)
    )
    featured_families = cazy_counts.head(args.top_cazy_families).index.tolist()
    family_colors = dict(zip(featured_families, PALETTE))
    family_colors["Other assigned"] = "#C7C7C7"
    family_colors["Unassigned"] = "#FFFFFF"

    nodes["figure2a_fill_group"] = nodes["cazy_base"].where(
        nodes["cazy_base"].isin(featured_families), "Other assigned"
    )
    nodes.loc[
        nodes["cazy_base"].eq("Unassigned"), "figure2a_fill_group"
    ] = "Unassigned"
    nodes["figure2a_fill_color"] = nodes["figure2a_fill_group"].map(family_colors)
    nodes["figure2a_is_virus"] = nodes["is_virus"].map(parse_bool)
    nodes["figure2a_border_color"] = nodes["figure2a_is_virus"].map(
        {True: "#111111", False: "#FFFFFF"}
    )
    nodes["figure2a_border_width"] = nodes["figure2a_is_virus"].map(
        {True: 2.4, False: 0.3}
    )
    nodes["figure2a_node_size"] = nodes["figure2a_is_virus"].map(
        {True: 13.0, False: 7.0}
    )
    nodes["figure2a_node_shape"] = nodes["figure2a_is_virus"].map(
        {True: "ELLIPSE", False: "ELLIPSE"}
    )

    # Community labels are attached to one high-strength representative each.
    community_lookup = communities.set_index("community").to_dict("index")
    nodes["figure2a_label"] = ""
    nodes["figure2a_major_community"] = False
    nodes["figure2a_virus_enriched_community"] = False
    for community, group in nodes.groupby(args.community_column, sort=False):
        context = community_lookup.get(str(community), {})
        community_n = int(context.get("community_n", len(group)))
        enrichment = str(context.get("viral_enrichment_class", ""))
        major = community_n >= args.major_community_min_n
        special = enrichment == "virus_enriched"
        if not (major or special):
            continue
        representative = group.loc[
            pd.to_numeric(group["within_community_strength"], errors="coerce")
            .fillna(-math.inf)
            .idxmax()
        ]
        identity = str(context.get("dominant_cazy_family", "")).strip()
        identity_class = str(context.get("cazy_identity_class", "")).strip()
        purity = float(context.get("dominant_cazy_purity", 0) or 0)
        if identity_class == "single_family_dominant":
            descriptor = identity
        elif identity:
            descriptor = f"{identity}-leading mixed"
        else:
            descriptor = "CAZy-unassigned"
        label = f"{community}  {descriptor}  n={community_n}"
        nodes.loc[
            nodes[args.node_id_column].eq(representative[args.node_id_column]),
            "figure2a_label",
        ] = label
        nodes.loc[
            nodes[args.community_column].eq(community),
            "figure2a_major_community",
        ] = major
        nodes.loc[
            nodes[args.community_column].eq(community),
            "figure2a_virus_enriched_community",
        ] = special

    compact_columns = [
        args.node_id_column,
        "figure2a_label",
        args.community_column,
        "connected_component",
        "is_giant_component",
        "is_graph_isolate",
        "figure2a_is_virus",
        "Group_corrected",
        "Species",
        "cazy_base",
        "figure2a_fill_group",
        "figure2a_fill_color",
        "figure2a_border_color",
        "figure2a_border_width",
        "figure2a_node_size",
        "figure2a_node_shape",
        "weighted_degree",
        "within_community_strength",
        "participation_coefficient",
        "placement_major",
        "figure2b_group",
        "figure2a_major_community",
        "figure2a_virus_enriched_community",
    ]
    compact_columns = [column for column in compact_columns if column in nodes.columns]
    compact = nodes[compact_columns].copy()
    compact = compact.rename(columns={args.node_id_column: "id"})

    backbone, backbone_stats = build_backbone(
        edges,
        node_ids,
        args.source_column,
        args.target_column,
        args.weight_column,
        args.top_edges_per_node,
    )
    backbone = backbone.rename(
        columns={
            args.source_column: "source",
            args.target_column: "target",
        }
    )
    full_edges = edges.rename(
        columns={
            args.source_column: "source",
            args.target_column: "target",
        }
    ).copy()
    full_edges["interaction"] = "structural_similarity"

    connected = compact.loc[~compact["is_graph_isolate"].map(parse_bool)].copy()
    isolates = compact.loc[compact["is_graph_isolate"].map(parse_bool)].copy()

    compact.to_csv(outdir / "figure2a_nodes_all.csv", index=False)
    connected.to_csv(outdir / "figure2a_nodes_connected.csv", index=False)
    isolates.to_csv(outdir / "figure2a_nodes_isolates.csv", index=False)
    full_edges.to_csv(outdir / "figure2a_edges_full.csv", index=False)
    backbone.to_csv(outdir / "figure2a_edges_display_backbone.csv", index=False)

    legend_rows = []
    for family in featured_families:
        legend_rows.append(
            {
                "fill_group": family,
                "color": family_colors[family],
                "node_n": int((nodes["cazy_base"] == family).sum()),
            }
        )
    for family in ["Other assigned", "Unassigned"]:
        legend_rows.append(
            {
                "fill_group": family,
                "color": family_colors[family],
                "node_n": int((nodes["figure2a_fill_group"] == family).sum()),
            }
        )
    pd.DataFrame(legend_rows).to_csv(
        outdir / "figure2a_cazy_color_legend.csv", index=False
    )

    run = {
        "input_files": {
            "nodes": str(Path(args.nodes).resolve()),
            "edges": str(Path(args.edges).resolve()),
            "communities": str(Path(args.communities).resolve()),
        },
        "input_sha256": {
            "nodes": sha256(args.nodes),
            "edges": sha256(args.edges),
            "communities": sha256(args.communities),
        },
        "counts": {
            "nodes": int(len(nodes)),
            "full_edges": int(len(edges)),
            "connected_nodes": int(len(connected)),
            "isolates": int(len(isolates)),
            "viral_nodes": int(nodes["figure2a_is_virus"].sum()),
            "viral_isolates": int(
                (
                    nodes["figure2a_is_virus"]
                    & nodes["is_graph_isolate"].map(parse_bool)
                ).sum()
            ),
            **backbone_stats,
        },
        "style": {
            "featured_cazy_families": featured_families,
            "family_colors": family_colors,
            "viral_encoding": "black 2.4-width border and 13-unit node",
            "nonviral_encoding": "white 0.3-width border and 7-unit node",
        },
        "important_note": (
            "The display backbone is a rendering-only subset. Community "
            "detection and all quantitative results use the full edge set."
        ),
    }
    with open(outdir / "figure2a_preparation_manifest.json", "w") as handle:
        json.dump(run, handle, indent=2)

    print(json.dumps(run["counts"], indent=2))


if __name__ == "__main__":
    main()
