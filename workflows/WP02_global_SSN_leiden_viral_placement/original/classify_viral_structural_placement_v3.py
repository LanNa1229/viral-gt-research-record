#!/usr/bin/env python3
"""Classify viral proteins in a Leiden-partitioned structural network.

Version 3 separates five questions that must not be collapsed into one label:

1. Which Leiden community contains the protein?
2. How well annotated is that community and what is its CAZy-family composition?
3. Is an assigned node concordant with the community's dominant CAZy family?
4. Does an unassigned node have direct, family-specific structural-neighbour support?
5. Is the node internal, peripheral, or an intercommunity connector?

The script also supports curated known-viral-GT metadata, ID lookup, refined
bridge subtypes, and both CSV and TSV outputs.

Required Python packages: pandas, numpy, scipy
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


SCRIPT_VERSION = "3.0.0"
TRUE_VALUES = {"1", "true", "t", "yes", "y", "virus", "viral"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "cellular", "nonvirus", "non-virus", ""}
RELIABLE_FAMILY_COMPOSITIONS = {
    "single_family",
    "dominant_family_with_minor_admixture",
}
UNANNOTATED_COMPOSITIONS = {"annotation_depleted"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Classify viral-protein placement from a Leiden node table and a "
            "cleaned undirected structural-similarity edge table."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--nodes", required=True, help="Node CSV/TSV with community assignments.")
    parser.add_argument("--edges", required=True, help="Undirected edge CSV/TSV.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--node-id-column", default="name")
    parser.add_argument("--source-column", default="source")
    parser.add_argument("--target-column", default="target")
    parser.add_argument("--community-column", default="community_global_min_tm_r1")
    parser.add_argument("--virus-column", default="is_virus")
    parser.add_argument("--cazy-family-column", default="CAZy_family")
    parser.add_argument(
        "--gtfam-column",
        default="GTFam",
        help="Row-wise fallback used when the CAZy-family value is absent or unparsable.",
    )
    parser.add_argument(
        "--weight-column",
        default="weight_global_min_tm",
        help="Edge-weight column; use 'none' for an unweighted analysis.",
    )

    composition = parser.add_argument_group("community-composition thresholds")
    composition.add_argument(
        "--known-family-min-annotated",
        type=int,
        default=5,
        help="Minimum CAZy-annotated members for a sufficiently supported composition call.",
    )
    composition.add_argument("--known-family-min-purity", type=float, default=0.80)
    composition.add_argument("--mixed-min-second-family-n", type=int, default=2)
    composition.add_argument("--mixed-min-second-family-fraction", type=float, default=0.20)
    composition.add_argument(
        "--small-community-max",
        type=int,
        default=4,
        help="Maximum size of a non-singleton community labeled small.",
    )
    composition.add_argument(
        "--large-community-min",
        type=int,
        default=20,
        help="Minimum size of a community labeled large.",
    )
    composition.add_argument("--virus-enrichment-q", type=float, default=0.05)

    topology = parser.add_argument_group("node-topology thresholds")
    topology.add_argument(
        "--peripheral-quantile",
        type=float,
        default=0.25,
        help="Within-community-strength quantile below which a family-context node is peripheral.",
    )
    topology.add_argument("--bridge-min-degree", type=int, default=3)
    topology.add_argument("--bridge-min-external-neighbors", type=int, default=1)
    topology.add_argument("--bridge-min-external-strength-fraction", type=float, default=0.20)
    topology.add_argument("--bridge-min-participation", type=float, default=0.30)
    topology.add_argument(
        "--bridge-hub-min-external-communities",
        type=int,
        default=3,
        help="External-community count used to call a multi-community hub.",
    )

    support = parser.add_argument_group("direct family-support thresholds")
    support.add_argument("--family-support-min-neighbors", type=int, default=2)
    support.add_argument("--family-support-min-strength-fraction", type=float, default=0.67)
    support.add_argument("--family-support-min-margin", type=float, default=0.20)

    evidence = parser.add_argument_group("viral evidence and lookup")
    evidence.add_argument(
        "--viral-reference-metadata",
        help=(
            "Optional CSV/TSV. Rows matching the node-ID column are curated viral references. "
            "Optional columns: viral_evidence_class, known_gt_name, evidence_type, aliases."
        ),
    )
    evidence.add_argument("--reference-id-column", default="name")
    evidence.add_argument("--reference-class-column", default="viral_evidence_class")
    evidence.add_argument("--reference-name-column", default="known_gt_name")
    evidence.add_argument("--reference-evidence-column", default="evidence_type")
    evidence.add_argument("--reference-alias-column", default="aliases")
    evidence.add_argument(
        "--unlisted-viral-class",
        default="new_viral_candidate",
        help=(
            "Class assigned to viral nodes absent from supplied reference metadata. "
            "Without metadata, previously unclassified viral nodes remain unclassified_viral."
        ),
    )
    evidence.add_argument(
        "--query-id",
        action="append",
        default=[],
        help="Exact node ID or curated alias to retrieve; repeat for multiple IDs.",
    )
    evidence.add_argument(
        "--query-ids-file",
        help="Optional text/CSV/TSV file; the first field of each nonempty line is queried.",
    )
    return parser.parse_args()


def read_delimited(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", keep_default_na=False)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_bool(value: object) -> bool:
    text = str(value).strip().casefold()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise ValueError(f"Cannot parse boolean value {value!r}")


def parse_cazy_family(value: object) -> str:
    """Return the base CAZy GT family, accepting labels such as GT2-DPM_like."""
    text = str(value).strip().upper()
    match = re.match(r"^(GT\d+)(?=$|[\s,;:_/\-])", text)
    return match.group(1) if match else ""


def fdr_bh(pvalues: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvalues), dtype=float)
    q = np.full_like(p, np.nan)
    finite = np.isfinite(p)
    if not finite.any():
        return q
    idx = np.where(finite)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(ranked, 1.0)
    return q


def normalized_entropy(counts: pd.Series) -> float:
    counts = counts[counts > 0]
    if len(counts) <= 1:
        return 0.0
    proportions = counts / counts.sum()
    return float(-(proportions * np.log(proportions)).sum() / np.log(len(counts)))


def effective_family_number(counts: pd.Series) -> float:
    counts = counts[counts > 0]
    if counts.empty:
        return 0.0
    proportions = counts / counts.sum()
    return float(np.exp(-(proportions * np.log(proportions)).sum()))


def validate_probability(name: str, value: float) -> None:
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must lie between 0 and 1; received {value}")


def community_size_class(n: int, args: argparse.Namespace) -> str:
    if n == 1:
        return "singleton"
    if n <= args.small_community_max:
        return "small"
    if n >= args.large_community_min:
        return "large"
    return "medium"


def community_composition_class(
    annotated_n: int,
    family_count: int,
    purity: float,
    second_n: int,
    second_fraction: float,
    args: argparse.Namespace,
) -> str:
    if annotated_n == 0:
        return "annotation_depleted"
    if annotated_n < args.known_family_min_annotated:
        return "sparsely_annotated"
    if family_count == 1:
        return "single_family"
    if purity >= args.known_family_min_purity:
        return "dominant_family_with_minor_admixture"
    if (
        second_n >= args.mixed_min_second_family_n
        and second_fraction >= args.mixed_min_second_family_fraction
    ):
        return "genuine_mixed_family"
    return "composition_ambiguous"


def prepare_nodes(args: argparse.Namespace) -> pd.DataFrame:
    nodes = read_delimited(args.nodes)
    required = [args.node_id_column, args.community_column, args.virus_column]
    missing = [column for column in required if column not in nodes.columns]
    if missing:
        raise ValueError(f"Node table lacks required columns: {missing}")

    nodes[args.node_id_column] = nodes[args.node_id_column].astype(str).str.strip()
    if nodes[args.node_id_column].eq("").any():
        raise ValueError("Node table contains empty node IDs")
    if nodes[args.node_id_column].duplicated().any():
        examples = nodes.loc[
            nodes[args.node_id_column].duplicated(), args.node_id_column
        ].head().tolist()
        raise ValueError(f"Node IDs must be unique; duplicate examples: {examples}")

    nodes[args.community_column] = nodes[args.community_column].astype(str).str.strip()
    if nodes[args.community_column].eq("").any():
        raise ValueError("Community assignments contain empty values")
    nodes["_is_virus"] = nodes[args.virus_column].map(parse_bool)

    if args.cazy_family_column in nodes.columns:
        primary_raw = nodes[args.cazy_family_column].astype(str).str.strip()
        primary = primary_raw.map(parse_cazy_family)
    else:
        primary_raw = pd.Series("", index=nodes.index, dtype=object)
        primary = pd.Series("", index=nodes.index, dtype=object)

    if args.gtfam_column in nodes.columns:
        fallback_raw = nodes[args.gtfam_column].astype(str).str.strip()
        fallback = fallback_raw.map(parse_cazy_family)
    else:
        fallback_raw = pd.Series("", index=nodes.index, dtype=object)
        fallback = pd.Series("", index=nodes.index, dtype=object)

    if args.cazy_family_column not in nodes.columns and args.gtfam_column not in nodes.columns:
        raise ValueError(
            f"Neither {args.cazy_family_column!r} nor fallback {args.gtfam_column!r} "
            "exists in the node table"
        )

    use_fallback = primary.eq("") & fallback.ne("")
    family = primary.where(~use_fallback, fallback)
    nodes["_cazy_family"] = family
    nodes["cazy_family_normalized"] = family
    nodes["cazy_family_parse_source"] = np.select(
        [primary.ne(""), use_fallback],
        [args.cazy_family_column, args.gtfam_column],
        default="unassigned",
    )
    nodes["cazy_status"] = np.where(family.ne(""), "CAZy_assigned", "CAZy_unassigned")
    nodes["cazy_family_raw_primary"] = primary_raw
    nodes["cazy_family_raw_fallback"] = fallback_raw
    return nodes


def merge_viral_reference_metadata(
    nodes: pd.DataFrame, args: argparse.Namespace
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    result = nodes.copy()
    defaults = {
        "viral_evidence_class": "",
        "known_gt_name": "",
        "evidence_type": "",
        "known_gt_aliases": "",
    }
    for column, default in defaults.items():
        if column not in result.columns:
            result[column] = default
        else:
            result[column] = result[column].astype(str).str.strip()

    result.loc[~result["_is_virus"], "viral_evidence_class"] = "nonviral_reference"
    metadata_unmatched = pd.DataFrame(columns=[args.reference_id_column])
    stats = {
        "reference_metadata_rows": 0,
        "reference_metadata_matched_nodes": 0,
        "reference_metadata_unmatched_rows": 0,
    }

    if args.viral_reference_metadata:
        metadata = read_delimited(args.viral_reference_metadata)
        if args.reference_id_column not in metadata.columns:
            raise ValueError(
                f"Reference metadata lacks ID column {args.reference_id_column!r}"
            )
        metadata[args.reference_id_column] = (
            metadata[args.reference_id_column].astype(str).str.strip()
        )
        if metadata[args.reference_id_column].eq("").any():
            raise ValueError("Reference metadata contains empty IDs")
        if metadata[args.reference_id_column].duplicated().any():
            examples = metadata.loc[
                metadata[args.reference_id_column].duplicated(), args.reference_id_column
            ].head().tolist()
            raise ValueError(f"Reference metadata IDs must be unique; examples: {examples}")

        node_ids = set(result[args.node_id_column])
        matched = metadata[args.reference_id_column].isin(node_ids)
        metadata_unmatched = metadata.loc[~matched].copy()
        matched_metadata = metadata.loc[matched].copy()
        stats = {
            "reference_metadata_rows": len(metadata),
            "reference_metadata_matched_nodes": len(matched_metadata),
            "reference_metadata_unmatched_rows": len(metadata_unmatched),
        }

        metadata_index = matched_metadata.set_index(args.reference_id_column)
        node_index = result[args.node_id_column]
        is_curated = node_index.isin(metadata_index.index)
        curated_nonviral = is_curated & ~result["_is_virus"]
        if curated_nonviral.any():
            examples = result.loc[
                curated_nonviral, args.node_id_column
            ].head().tolist()
            raise ValueError(
                "Curated viral-reference metadata matched nodes marked nonviral; "
                f"examples: {examples}"
            )

        def metadata_map(column: str, default: str = "") -> pd.Series:
            if column not in metadata_index.columns:
                return pd.Series(default, index=result.index, dtype=object)
            return node_index.map(metadata_index[column]).fillna(default).astype(str).str.strip()

        curated_class = metadata_map(args.reference_class_column, "known_viral_GT")
        curated_class = curated_class.mask(
            is_curated & curated_class.eq(""), "known_viral_GT"
        )
        result.loc[is_curated, "viral_evidence_class"] = curated_class[is_curated]
        result.loc[is_curated, "known_gt_name"] = metadata_map(
            args.reference_name_column
        )[is_curated]
        result.loc[is_curated, "evidence_type"] = metadata_map(
            args.reference_evidence_column
        )[is_curated]
        result.loc[is_curated, "known_gt_aliases"] = metadata_map(
            args.reference_alias_column
        )[is_curated]

        unlisted_viral = result["_is_virus"] & ~is_curated
        result.loc[
            unlisted_viral & result["viral_evidence_class"].eq(""),
            "viral_evidence_class",
        ] = args.unlisted_viral_class
    else:
        result.loc[
            result["_is_virus"] & result["viral_evidence_class"].eq(""),
            "viral_evidence_class",
        ] = "unclassified_viral"

    result["is_known_viral_gt"] = (
        result["viral_evidence_class"].str.casefold().eq("known_viral_gt")
    )
    return result, metadata_unmatched, stats


def prepare_edges(
    args: argparse.Namespace, valid_nodes: set[str]
) -> tuple[pd.DataFrame, dict[str, int]]:
    edges = read_delimited(args.edges)
    required = [args.source_column, args.target_column]
    missing = [column for column in required if column not in edges.columns]
    if missing:
        raise ValueError(f"Edge table lacks required columns: {missing}")

    source = edges[args.source_column].astype(str).str.strip()
    target = edges[args.target_column].astype(str).str.strip()
    if args.weight_column.casefold() == "none":
        weight = pd.Series(np.ones(len(edges), dtype=float), index=edges.index)
    else:
        if args.weight_column not in edges.columns:
            raise ValueError(f"Edge table lacks weight column {args.weight_column!r}")
        weight = pd.to_numeric(edges[args.weight_column], errors="coerce")
        if weight.isna().any():
            raise ValueError(
                f"Weight column {args.weight_column!r} contains missing/non-numeric values"
            )
        if (weight < 0).any():
            raise ValueError("Negative edge weights are not supported")

    work = pd.DataFrame({"source": source, "target": target, "weight": weight})
    n_input = len(work)
    missing_node = ~work["source"].isin(valid_nodes) | ~work["target"].isin(valid_nodes)
    n_missing = int(missing_node.sum())
    work = work.loc[~missing_node].copy()
    self_hit = work["source"].eq(work["target"])
    n_self = int(self_hit.sum())
    work = work.loc[~self_hit].copy()

    source_array = work["source"].to_numpy(dtype=str)
    target_array = work["target"].to_numpy(dtype=str)
    work["u"] = np.where(source_array <= target_array, source_array, target_array)
    work["v"] = np.where(source_array <= target_array, target_array, source_array)
    n_before_collapse = len(work)
    work = work.groupby(["u", "v"], as_index=False, sort=False)["weight"].mean()
    work = work.rename(columns={"u": "source", "v": "target"})
    stats = {
        "edge_rows_input": n_input,
        "edge_rows_missing_node_removed": n_missing,
        "self_edges_removed": n_self,
        "duplicate_directional_or_undirected_rows_collapsed": n_before_collapse
        - len(work),
        "unique_undirected_edges_analyzed": len(work),
    }
    return work, stats


def add_within_family_support(
    result: pd.DataFrame,
    source_idx: np.ndarray,
    target_idx: np.ndarray,
    weight: np.ndarray,
    same: np.ndarray,
) -> pd.DataFrame:
    n = len(result)
    family = result["_cazy_family"].to_numpy(dtype=str)
    directed = pd.concat(
        [
            pd.DataFrame(
                {
                    "node_index": source_idx[same],
                    "neighbor_family": family[target_idx[same]],
                    "strength": weight[same],
                }
            ),
            pd.DataFrame(
                {
                    "node_index": target_idx[same],
                    "neighbor_family": family[source_idx[same]],
                    "strength": weight[same],
                }
            ),
        ],
        ignore_index=True,
    )
    directed = directed.loc[directed["neighbor_family"].ne("")].copy()
    grouped = (
        directed.groupby(["node_index", "neighbor_family"], as_index=False)
        .agg(neighbor_n=("strength", "size"), strength=("strength", "sum"))
        .sort_values(
            ["node_index", "strength", "neighbor_n", "neighbor_family"],
            ascending=[True, False, False, True],
        )
    )
    grouped["rank"] = grouped.groupby("node_index").cumcount() + 1

    annotated_neighbor_n = np.zeros(n, dtype=int)
    annotated_neighbor_strength = np.zeros(n, dtype=float)
    if not grouped.empty:
        total = grouped.groupby("node_index").agg(
            neighbor_n=("neighbor_n", "sum"), strength=("strength", "sum")
        )
        idx = total.index.to_numpy(dtype=int)
        annotated_neighbor_n[idx] = total["neighbor_n"].to_numpy(dtype=int)
        annotated_neighbor_strength[idx] = total["strength"].to_numpy(dtype=float)

    result["within_cazy_annotated_neighbor_n"] = annotated_neighbor_n
    result["within_cazy_annotated_neighbor_strength"] = annotated_neighbor_strength
    for rank, label in [(1, "best"), (2, "second")]:
        family_values = np.full(n, "", dtype=object)
        neighbor_values = np.zeros(n, dtype=int)
        strength_values = np.zeros(n, dtype=float)
        rows = grouped.loc[grouped["rank"].eq(rank)]
        if not rows.empty:
            idx = rows["node_index"].to_numpy(dtype=int)
            family_values[idx] = rows["neighbor_family"].to_numpy(dtype=str)
            neighbor_values[idx] = rows["neighbor_n"].to_numpy(dtype=int)
            strength_values[idx] = rows["strength"].to_numpy(dtype=float)
        result[f"within_{label}_cazy_family"] = family_values
        result[f"within_{label}_cazy_family_neighbor_n"] = neighbor_values
        result[f"within_{label}_cazy_family_strength"] = strength_values
        result[f"within_{label}_cazy_family_strength_fraction"] = np.divide(
            strength_values,
            annotated_neighbor_strength,
            out=np.zeros(n, dtype=float),
            where=annotated_neighbor_strength > 0,
        )

    result["within_best_vs_second_strength_margin"] = np.divide(
        result["within_best_cazy_family_strength"].to_numpy(dtype=float)
        - result["within_second_cazy_family_strength"].to_numpy(dtype=float),
        annotated_neighbor_strength,
        out=np.zeros(n, dtype=float),
        where=annotated_neighbor_strength > 0,
    )
    return result


def compute_metrics(
    nodes: pd.DataFrame,
    edges: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    n = len(nodes)
    node_id = args.node_id_column
    community_column = args.community_column
    id_to_index = pd.Series(np.arange(n), index=nodes[node_id]).to_dict()
    source_idx = edges["source"].map(id_to_index).to_numpy(dtype=int)
    target_idx = edges["target"].map(id_to_index).to_numpy(dtype=int)
    weight = edges["weight"].to_numpy(dtype=float)
    community = nodes[community_column].to_numpy(dtype=str)

    degree = np.bincount(source_idx, minlength=n) + np.bincount(target_idx, minlength=n)
    total_strength = (
        np.bincount(source_idx, weights=weight, minlength=n)
        + np.bincount(target_idx, weights=weight, minlength=n)
    )
    same = community[source_idx] == community[target_idx]
    within_degree = (
        np.bincount(source_idx[same], minlength=n)
        + np.bincount(target_idx[same], minlength=n)
    )
    within_strength = (
        np.bincount(source_idx[same], weights=weight[same], minlength=n)
        + np.bincount(target_idx[same], weights=weight[same], minlength=n)
    )
    external_degree = degree - within_degree
    external_strength = total_strength - within_strength
    within_fraction = np.divide(
        within_strength, total_strength, out=np.zeros(n, dtype=float), where=total_strength > 0
    )
    external_fraction = np.divide(
        external_strength, total_strength, out=np.zeros(n, dtype=float), where=total_strength > 0
    )

    left = pd.DataFrame(
        {
            "node_index": source_idx,
            "neighbor_community": community[target_idx],
            "strength": weight,
        }
    )
    right = pd.DataFrame(
        {
            "node_index": target_idx,
            "neighbor_community": community[source_idx],
            "strength": weight,
        }
    )
    neighbor_by_community = (
        pd.concat([left, right], ignore_index=True)
        .groupby(["node_index", "neighbor_community"], as_index=False, sort=False)
        .agg(strength=("strength", "sum"), edge_n=("strength", "size"))
    )
    denominator = total_strength[
        neighbor_by_community["node_index"].to_numpy(dtype=int)
    ]
    neighbor_by_community["fraction_squared"] = np.divide(
        neighbor_by_community["strength"].to_numpy(dtype=float),
        denominator,
        out=np.zeros(len(neighbor_by_community), dtype=float),
        where=denominator > 0,
    ) ** 2
    participation = np.zeros(n, dtype=float)
    summed_squares = neighbor_by_community.groupby("node_index")[
        "fraction_squared"
    ].sum()
    participation[summed_squares.index.to_numpy(dtype=int)] = (
        1.0 - summed_squares.to_numpy()
    )
    participation = np.clip(participation, 0.0, 1.0)

    own_community = community[
        neighbor_by_community["node_index"].to_numpy(dtype=int)
    ]
    external_neighbor_profile = neighbor_by_community.loc[
        neighbor_by_community["neighbor_community"].to_numpy(dtype=str) != own_community
    ].copy()
    external_community_count = np.zeros(n, dtype=int)
    if not external_neighbor_profile.empty:
        counts = external_neighbor_profile.groupby("node_index")[
            "neighbor_community"
        ].nunique()
        external_community_count[counts.index.to_numpy(dtype=int)] = counts.to_numpy(
            dtype=int
        )

    cross_source = source_idx[~same]
    cross_target = target_idx[~same]
    cross_source_comm = community[cross_source]
    cross_target_comm = community[cross_target]
    pair_a = np.where(
        cross_source_comm <= cross_target_comm, cross_source_comm, cross_target_comm
    )
    pair_b = np.where(
        cross_source_comm <= cross_target_comm, cross_target_comm, cross_source_comm
    )
    cross = pd.DataFrame(
        {
            "source_index": cross_source,
            "target_index": cross_target,
            "pair_a": pair_a,
            "pair_b": pair_b,
        }
    )
    min_pair_edge_count = np.zeros(n, dtype=int)
    if not cross.empty:
        pair_counts = cross.groupby(["pair_a", "pair_b"]).size().rename("pair_edge_count")
        cross = cross.join(pair_counts, on=["pair_a", "pair_b"])
        for endpoint in ("source_index", "target_index"):
            endpoint_min = cross.groupby(endpoint)["pair_edge_count"].min()
            idx = endpoint_min.index.to_numpy(dtype=int)
            values = endpoint_min.to_numpy(dtype=int)
            existing = min_pair_edge_count[idx]
            min_pair_edge_count[idx] = np.where(
                existing == 0, values, np.minimum(existing, values)
            )

    result = nodes.copy()
    result["degree"] = degree
    result["weighted_degree"] = total_strength
    result["within_community_degree"] = within_degree
    result["within_community_strength"] = within_strength
    result["external_degree"] = external_degree
    result["external_strength"] = external_strength
    result["within_strength_fraction"] = within_fraction
    result["external_strength_fraction"] = external_fraction
    result["participation_coefficient"] = participation
    result["external_community_count"] = external_community_count
    result["minimum_incident_community_pair_edge_count"] = min_pair_edge_count
    result["sole_direct_community_pair_connector"] = min_pair_edge_count == 1
    result = add_within_family_support(
        result, source_idx, target_idx, weight, same
    )

    community_codes, unique_communities = pd.factorize(
        result[community_column], sort=False
    )
    n_communities = len(unique_communities)
    internal_edge_count = np.bincount(
        community_codes[source_idx[same]], minlength=n_communities
    )
    internal_strength = np.bincount(
        community_codes[source_idx[same]],
        weights=weight[same],
        minlength=n_communities,
    )
    external_edge_count = (
        np.bincount(community_codes[source_idx[~same]], minlength=n_communities)
        + np.bincount(community_codes[target_idx[~same]], minlength=n_communities)
    )
    external_community_strength = (
        np.bincount(
            community_codes[source_idx[~same]],
            weights=weight[~same],
            minlength=n_communities,
        )
        + np.bincount(
            community_codes[target_idx[~same]],
            weights=weight[~same],
            minlength=n_communities,
        )
    )
    edge_summary = pd.DataFrame(
        {
            "community": unique_communities,
            "internal_edges": internal_edge_count,
            "external_edges": external_edge_count,
            "internal_strength": internal_strength,
            "external_strength": external_community_strength,
        }
    ).set_index("community")

    total_nodes = len(result)
    total_virus = int(result["_is_virus"].sum())
    rows: list[dict[str, object]] = []
    for comm, subset in result.groupby(community_column, sort=False):
        family_counts = subset.loc[
            subset["_cazy_family"].ne(""), "_cazy_family"
        ].value_counts()
        annotated_n = int(family_counts.sum())
        dominant = str(family_counts.index[0]) if annotated_n else ""
        dominant_n = int(family_counts.iloc[0]) if annotated_n else 0
        second = str(family_counts.index[1]) if len(family_counts) > 1 else ""
        second_n = int(family_counts.iloc[1]) if len(family_counts) > 1 else 0
        purity = dominant_n / annotated_n if annotated_n else np.nan
        second_fraction = second_n / annotated_n if annotated_n else 0.0
        virus_n = int(subset["_is_virus"].sum())
        n_nodes = len(subset)
        pvalue = (
            float(hypergeom.sf(virus_n - 1, total_nodes, total_virus, n_nodes))
            if total_virus
            else np.nan
        )
        possible_edges = n_nodes * (n_nodes - 1) / 2
        internal_edges = int(edge_summary.loc[comm, "internal_edges"])
        internal_s = float(edge_summary.loc[comm, "internal_strength"])
        external_s = float(edge_summary.loc[comm, "external_strength"])
        composition = community_composition_class(
            annotated_n,
            len(family_counts),
            purity,
            second_n,
            second_fraction,
            args,
        )
        rows.append(
            {
                "community": comm,
                "community_n": n_nodes,
                "community_size_class": community_size_class(n_nodes, args),
                "internal_edges": internal_edges,
                "external_edges": int(edge_summary.loc[comm, "external_edges"]),
                "internal_strength": internal_s,
                "external_strength": external_s,
                "density": internal_edges / possible_edges if possible_edges else 0.0,
                "conductance": external_s / (2 * internal_s + external_s)
                if (2 * internal_s + external_s) > 0
                else 0.0,
                "cazy_annotated_n": annotated_n,
                "cazy_annotation_coverage": annotated_n / n_nodes,
                "cazy_family_count": int(len(family_counts)),
                "dominant_cazy_family": dominant,
                "dominant_cazy_n": dominant_n,
                "second_cazy_family": second,
                "second_cazy_n": second_n,
                "second_cazy_fraction": second_fraction,
                "cazy_purity": purity,
                "cazy_normalized_entropy": normalized_entropy(family_counts),
                "cazy_effective_family_n": effective_family_number(family_counts),
                "community_composition_class": composition,
                "has_reliable_family_context": composition
                in RELIABLE_FAMILY_COMPOSITIONS,
                "virus_n": virus_n,
                "virus_fraction": virus_n / n_nodes,
                "virus_enrichment_p": pvalue,
            }
        )
    community_summary = pd.DataFrame(rows)
    community_summary["virus_enrichment_q"] = fdr_bh(
        community_summary["virus_enrichment_p"]
    )
    community_summary["community_context"] = community_summary[
        "community_composition_class"
    ] + "_neighborhood"
    dark = community_summary["community_composition_class"].eq(
        "annotation_depleted"
    ) & (community_summary["virus_enrichment_q"] < args.virus_enrichment_q)
    community_summary.loc[
        dark, "community_context"
    ] = "virus_enriched_dark_neighborhood"
    return result, community_summary, external_neighbor_profile


def add_external_profiles(
    result: pd.DataFrame,
    community_summary: pd.DataFrame,
    external_profile: pd.DataFrame,
) -> pd.DataFrame:
    n = len(result)
    lookup = community_summary.set_index("community")
    profile = external_profile.copy()
    if profile.empty:
        result["top_external_community"] = ""
        result["top_external_community_strength"] = 0.0
        result["top_external_community_edge_n"] = 0
        result["top_external_community_composition_class"] = ""
        result["top_external_dominant_cazy_family"] = ""
        result["external_community_profile"] = ""
        return result

    profile["composition"] = profile["neighbor_community"].map(
        lookup["community_composition_class"]
    )
    profile["dominant_family"] = profile["neighbor_community"].map(
        lookup["dominant_cazy_family"]
    )
    profile["annotated_n"] = profile["neighbor_community"].map(
        lookup["cazy_annotated_n"]
    )
    profile["purity"] = profile["neighbor_community"].map(lookup["cazy_purity"])
    profile = profile.sort_values(
        ["node_index", "strength", "edge_n", "neighbor_community"],
        ascending=[True, False, False, True],
    )
    profile["rank"] = profile.groupby("node_index").cumcount() + 1
    top = profile.loc[profile["rank"].eq(1)]

    top_comm = np.full(n, "", dtype=object)
    top_strength = np.zeros(n, dtype=float)
    top_edge_n = np.zeros(n, dtype=int)
    top_composition = np.full(n, "", dtype=object)
    top_family = np.full(n, "", dtype=object)
    idx = top["node_index"].to_numpy(dtype=int)
    top_comm[idx] = top["neighbor_community"].to_numpy(dtype=str)
    top_strength[idx] = top["strength"].to_numpy(dtype=float)
    top_edge_n[idx] = top["edge_n"].to_numpy(dtype=int)
    top_composition[idx] = top["composition"].fillna("").to_numpy(dtype=str)
    top_family[idx] = top["dominant_family"].fillna("").to_numpy(dtype=str)
    result["top_external_community"] = top_comm
    result["top_external_community_strength"] = top_strength
    result["top_external_community_edge_n"] = top_edge_n
    result["top_external_community_composition_class"] = top_composition
    result["top_external_dominant_cazy_family"] = top_family

    def render_profile(group: pd.DataFrame) -> str:
        parts = []
        for row in group.itertuples(index=False):
            family = row.dominant_family if row.dominant_family else "NA"
            composition = row.composition if row.composition else "NA"
            parts.append(
                f"{row.neighbor_community}|{composition}|{family}|"
                f"edges={int(row.edge_n)}|strength={float(row.strength):.6g}"
            )
        return ";".join(parts)

    rendered = {
        int(node_index): render_profile(group.drop(columns=["node_index"]))
        for node_index, group in profile.groupby("node_index", sort=False)
    }
    result["external_community_profile"] = (
        result.index.to_series().map(rendered).fillna("")
    )
    return result


def classify_nodes(
    result: pd.DataFrame,
    community_summary: pd.DataFrame,
    external_profile: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    community_column = args.community_column
    lookup = community_summary.set_index("community")
    mapped_columns = [
        "community_n",
        "community_size_class",
        "cazy_annotated_n",
        "cazy_annotation_coverage",
        "cazy_family_count",
        "dominant_cazy_family",
        "dominant_cazy_n",
        "second_cazy_family",
        "second_cazy_n",
        "second_cazy_fraction",
        "cazy_purity",
        "cazy_normalized_entropy",
        "cazy_effective_family_n",
        "community_composition_class",
        "has_reliable_family_context",
        "virus_n",
        "virus_fraction",
        "virus_enrichment_q",
        "community_context",
    ]
    for column in mapped_columns:
        result[column] = result[community_column].map(lookup[column])
    result = add_external_profiles(result, community_summary, external_profile)

    cutoff = result.groupby(community_column)["within_community_strength"].quantile(
        args.peripheral_quantile
    )
    result["community_peripheral_strength_cutoff"] = result[community_column].map(
        cutoff
    )
    result["is_graph_isolate"] = result["degree"].eq(0)
    result["is_leiden_singleton"] = result["community_n"].eq(1)
    result["is_singleton_or_isolate"] = (
        result["is_graph_isolate"] | result["is_leiden_singleton"]
    )

    result["bridge_candidate"] = (
        ~result["is_singleton_or_isolate"]
        & (result["degree"] >= args.bridge_min_degree)
        & (result["within_community_degree"] >= 1)
        & (result["external_degree"] >= args.bridge_min_external_neighbors)
        & (
            result["external_strength_fraction"]
            >= args.bridge_min_external_strength_fraction
        )
        & (result["participation_coefficient"] >= args.bridge_min_participation)
    )
    result["bridge_support"] = np.select(
        [
            result["bridge_candidate"] & (result["external_degree"] >= 2),
            result["bridge_candidate"] & result["external_degree"].eq(1),
        ],
        ["multiple_external_edges", "single_external_edge_candidate"],
        default="not_bridge_candidate",
    )

    result["cazy_relation_to_community"] = np.select(
        [
            result["_cazy_family"].eq(""),
            result["dominant_cazy_family"].eq(""),
            result["_cazy_family"].eq(result["dominant_cazy_family"]),
        ],
        ["not_applicable_CAZy_unassigned", "community_without_dominant_CAZy", "concordant"],
        default="discordant",
    )
    result["assigned_family_relation_to_direct_support"] = np.select(
        [
            result["_cazy_family"].eq(""),
            result["within_best_cazy_family"].eq(""),
            result["_cazy_family"].eq(result["within_best_cazy_family"]),
        ],
        ["not_applicable_CAZy_unassigned", "no_direct_annotated_neighbor_support", "concordant"],
        default="discordant",
    )

    support_pass = (
        (result["within_best_cazy_family_neighbor_n"] >= args.family_support_min_neighbors)
        & (
            result["within_best_cazy_family_strength_fraction"]
            >= args.family_support_min_strength_fraction
        )
        & (
            result["within_best_vs_second_strength_margin"]
            >= args.family_support_min_margin
        )
    )
    unassigned = result["_cazy_family"].eq("")
    result["unassigned_family_support_call"] = np.select(
        [
            ~unassigned,
            result["within_cazy_annotated_neighbor_n"].eq(0),
            result["within_best_cazy_family_neighbor_n"]
            < args.family_support_min_neighbors,
            unassigned & ~support_pass,
            unassigned
            & support_pass
            & result["within_best_cazy_family"].eq(
                result["dominant_cazy_family"]
            ),
            unassigned & support_pass,
        ],
        [
            "not_applicable_CAZy_assigned",
            "no_direct_annotated_neighbor_support",
            "insufficient_independent_annotated_neighbors",
            "ambiguous_multi_family_direct_support",
            "direct_support_concordant_with_community",
            "direct_support_discordant_with_community",
        ],
        default="not_assessed",
    )
    result["directly_supported_cazy_family"] = np.where(
        unassigned & support_pass, result["within_best_cazy_family"], ""
    )
    result["direct_supports_community_dominant_family"] = (
        unassigned
        & support_pass
        & result["within_best_cazy_family"].eq(result["dominant_cazy_family"])
        & result["dominant_cazy_family"].ne("")
    )

    peripheral = (
        result["has_reliable_family_context"]
        & ~result["is_singleton_or_isolate"]
        & (
            result["within_community_strength"]
            < result["community_peripheral_strength_cutoff"]
        )
    )
    result["is_peripheral_family_context_node"] = peripheral
    result["topology_role"] = np.select(
        [
            result["is_singleton_or_isolate"],
            result["bridge_candidate"],
            peripheral,
        ],
        ["singleton_or_isolate", "intercommunity_connector", "peripheral"],
        default="embedded_or_internal",
    )

    own_composition = result["community_composition_class"]
    external_composition = result["top_external_community_composition_class"]
    own_reliable = own_composition.isin(RELIABLE_FAMILY_COMPOSITIONS)
    external_reliable = external_composition.isin(RELIABLE_FAMILY_COMPOSITIONS)
    own_mixed = own_composition.eq("genuine_mixed_family")
    external_mixed = external_composition.eq("genuine_mixed_family")
    own_unannotated = own_composition.isin(UNANNOTATED_COMPOSITIONS)
    external_unannotated = external_composition.isin(UNANNOTATED_COMPOSITIONS)
    bridge = result["bridge_candidate"]
    result["bridge_subtype"] = np.select(
        [
            ~bridge,
            bridge & result["external_degree"].eq(1),
            bridge
            & (
                result["external_community_count"]
                >= args.bridge_hub_min_external_communities
            ),
            bridge
            & own_reliable
            & external_reliable
            & result["dominant_cazy_family"].eq(
                result["top_external_dominant_cazy_family"]
            ),
            bridge
            & own_reliable
            & external_reliable
            & result["dominant_cazy_family"].ne(
                result["top_external_dominant_cazy_family"]
            ),
            bridge
            & ((own_reliable & external_mixed) | (own_mixed & external_reliable)),
            bridge
            & (
                (own_reliable & external_unannotated)
                | (own_unannotated & external_reliable)
            ),
            bridge & own_unannotated & external_unannotated,
            bridge
            & (
                own_composition.eq("sparsely_annotated")
                | external_composition.eq("sparsely_annotated")
            ),
        ],
        [
            "not_bridge_candidate",
            "weak_single_edge_connector",
            "multi_community_hub",
            "same_family_fragmentation_bridge",
            "inter_family_bridge_candidate",
            "family_to_mixed_connector",
            "family_to_unannotated_connector",
            "unannotated_community_connector",
            "sparse_annotation_connector",
        ],
        default="other_intercommunity_connector",
    )
    result["bridge_family_pair"] = np.where(
        bridge,
        result["dominant_cazy_family"].replace("", "NA")
        + "->"
        + result["top_external_dominant_cazy_family"].replace("", "NA"),
        "",
    )
    result["bridge_alignment_validation_status"] = np.where(
        bridge, "not_assessed_requires_alignment_coordinates", "not_applicable"
    )

    result["family_context_interpretation"] = np.select(
        [
            ~result["has_reliable_family_context"],
            result["_cazy_family"].ne("")
            & result["cazy_relation_to_community"].eq("concordant"),
            result["_cazy_family"].ne("")
            & result["cazy_relation_to_community"].eq("discordant"),
            result["direct_supports_community_dominant_family"],
            unassigned & support_pass,
            unassigned,
        ],
        [
            "community_without_reliable_single_or_dominant_family_context",
            "assigned_family_concordant_with_community",
            "assigned_family_discordant_with_community",
            "unassigned_directly_supported_for_community_dominant_family",
            "unassigned_directly_supported_for_different_family",
            "unassigned_community_membership_without_decisive_direct_family_support",
        ],
        default="not_assessed",
    )

    major = np.full(len(result), "nonviral_not_classified", dtype=object)
    virus = result["_is_virus"].to_numpy(dtype=bool)
    remaining = virus.copy()

    def assign(mask: np.ndarray, label: str) -> None:
        nonlocal remaining
        selected = remaining & mask
        major[selected] = label
        remaining &= ~selected

    assign(result["is_singleton_or_isolate"].to_numpy(dtype=bool), "singleton_or_isolate")
    assign(result["bridge_candidate"].to_numpy(dtype=bool), "intercommunity_bridge")
    assign(
        peripheral.to_numpy(dtype=bool),
        "peripheral_established_family_context",
    )
    assign(
        result["has_reliable_family_context"].to_numpy(dtype=bool),
        "established_family_context",
    )
    assign(
        result["community_composition_class"]
        .eq("genuine_mixed_family")
        .to_numpy(dtype=bool),
        "genuine_mixed_family_community",
    )
    assign(
        result["community_composition_class"]
        .eq("sparsely_annotated")
        .to_numpy(dtype=bool),
        "sparsely_annotated_community",
    )
    assign(
        result["community_context"]
        .eq("virus_enriched_dark_neighborhood")
        .to_numpy(dtype=bool),
        "virus_enriched_dark_community",
    )
    assign(
        result["community_composition_class"]
        .eq("annotation_depleted")
        .to_numpy(dtype=bool),
        "annotation_depleted_community",
    )
    assign(
        result["community_composition_class"]
        .eq("composition_ambiguous")
        .to_numpy(dtype=bool),
        "composition_ambiguous_community",
    )
    major[remaining] = "unresolved_structural_context"
    result["placement_major"] = major
    result["figure2b_class"] = np.where(
        result["_is_virus"],
        result["placement_major"]
        + "__"
        + result["cazy_status"].str.lower()
        + "__"
        + result["viral_evidence_class"].str.lower(),
        "",
    )
    result["figure2b_group"] = np.where(
        result["_is_virus"],
        result["community_composition_class"]
        + "__"
        + result["topology_role"]
        + "__"
        + result["cazy_status"].str.lower(),
        "",
    )
    return result


def load_query_ids(args: argparse.Namespace) -> list[str]:
    queries = [str(value).strip() for value in args.query_id if str(value).strip()]
    if args.query_ids_file:
        path = Path(args.query_ids_file)
        for line in path.read_text(encoding="utf-8").splitlines():
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            value = re.split(r"[\t,]", text, maxsplit=1)[0].strip()
            if value and value.casefold() not in {"name", "id", "query_id"}:
                queries.append(value)
    return list(dict.fromkeys(queries))


def query_placements(
    viral: pd.DataFrame, args: argparse.Namespace
) -> pd.DataFrame:
    queries = load_query_ids(args)
    if not queries:
        return pd.DataFrame()

    key_to_indices: dict[str, set[int]] = {}
    for idx, row in viral.iterrows():
        keys = [str(row[args.node_id_column]).strip()]
        aliases = str(row.get("known_gt_aliases", "")).strip()
        if aliases:
            keys.extend(
                token.strip()
                for token in re.split(r"[;,|]", aliases)
                if token.strip()
            )
        for key in keys:
            key_to_indices.setdefault(key.casefold(), set()).add(idx)

    output_rows: list[dict[str, object]] = []
    for query in queries:
        matches = sorted(key_to_indices.get(query.casefold(), set()))
        if not matches:
            output_rows.append(
                {"query_id": query, "query_match_status": "not_found"}
            )
            continue
        status = "unique" if len(matches) == 1 else "ambiguous_multiple_matches"
        for idx in matches:
            record = viral.loc[idx].to_dict()
            output_rows.append(
                {
                    "query_id": query,
                    "query_match_status": status,
                    **record,
                }
            )
    return pd.DataFrame(output_rows)


def write_table_both(df: pd.DataFrame, stem: Path) -> None:
    df.to_csv(stem.with_suffix(".csv"), index=False)
    df.to_csv(stem.with_suffix(".tsv"), sep="\t", index=False)


def write_outputs(
    result: pd.DataFrame,
    community_summary: pd.DataFrame,
    metadata_unmatched: pd.DataFrame,
    metadata_stats: dict[str, int],
    edge_stats: dict[str, int],
    args: argparse.Namespace,
) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    public_result = result.drop(columns=["_is_virus", "_cazy_family"])
    viral = public_result.loc[
        public_result[args.virus_column].map(parse_bool)
    ].copy()

    write_table_both(
        public_result, outdir / "cytoscape_nodes_with_figure2b_classes"
    )
    write_table_both(viral, outdir / "viral_structural_placement")
    write_table_both(community_summary, outdir / "community_context_summary")

    detailed_summary = (
        viral.groupby(
            [
                "viral_evidence_class",
                "cazy_status",
                "community_composition_class",
                "topology_role",
                "bridge_subtype",
                "placement_major",
            ],
            dropna=False,
        )
        .size()
        .rename("n_viral_proteins")
        .reset_index()
        .sort_values("n_viral_proteins", ascending=False)
    )
    detailed_summary["fraction_of_viral_proteins"] = (
        detailed_summary["n_viral_proteins"] / len(viral)
        if len(viral)
        else np.nan
    )
    write_table_both(detailed_summary, outdir / "viral_placement_detailed_summary")

    cross_tab = (
        viral.groupby(
            [
                "viral_evidence_class",
                "cazy_status",
                "community_composition_class",
                "topology_role",
            ],
            dropna=False,
        )
        .size()
        .rename("n_viral_proteins")
        .reset_index()
    )
    write_table_both(cross_tab, outdir / "known_vs_candidate_placement_crosstab")

    known = viral.loc[viral["is_known_viral_gt"]].copy()
    known = known.sort_values(
        ["known_gt_name", args.node_id_column], kind="stable"
    )
    write_table_both(known, outdir / "known_viral_gt_anchor_table")

    bridges = viral.loc[viral["bridge_candidate"]].copy()
    bridges = bridges.sort_values(
        [
            "external_community_count",
            "participation_coefficient",
            "external_strength_fraction",
            "weighted_degree",
        ],
        ascending=[False, False, False, False],
    )
    write_table_both(bridges, outdir / "bridge_candidates")

    query_results = query_placements(viral, args)
    if not query_results.empty:
        write_table_both(query_results, outdir / "viral_gt_query_results")
    if not metadata_unmatched.empty:
        write_table_both(metadata_unmatched, outdir / "unmatched_reference_metadata")

    config = {
        "script_version": SCRIPT_VERSION,
        "inputs": {
            "nodes": str(Path(args.nodes).resolve()),
            "nodes_sha256": sha256_file(args.nodes),
            "edges": str(Path(args.edges).resolve()),
            "edges_sha256": sha256_file(args.edges),
            "viral_reference_metadata": (
                str(Path(args.viral_reference_metadata).resolve())
                if args.viral_reference_metadata
                else None
            ),
            "viral_reference_metadata_sha256": (
                sha256_file(args.viral_reference_metadata)
                if args.viral_reference_metadata
                else None
            ),
        },
        "columns": {
            "node_id": args.node_id_column,
            "community": args.community_column,
            "virus": args.virus_column,
            "cazy_family": args.cazy_family_column,
            "gtfam_fallback": args.gtfam_column,
            "weight": args.weight_column,
        },
        "thresholds": {
            "known_family_min_annotated": args.known_family_min_annotated,
            "known_family_min_purity": args.known_family_min_purity,
            "mixed_min_second_family_n": args.mixed_min_second_family_n,
            "mixed_min_second_family_fraction": args.mixed_min_second_family_fraction,
            "small_community_max": args.small_community_max,
            "large_community_min": args.large_community_min,
            "peripheral_quantile": args.peripheral_quantile,
            "virus_enrichment_q": args.virus_enrichment_q,
            "bridge_min_degree": args.bridge_min_degree,
            "bridge_min_external_neighbors": args.bridge_min_external_neighbors,
            "bridge_min_external_strength_fraction": args.bridge_min_external_strength_fraction,
            "bridge_min_participation": args.bridge_min_participation,
            "bridge_hub_min_external_communities": args.bridge_hub_min_external_communities,
            "family_support_min_neighbors": args.family_support_min_neighbors,
            "family_support_min_strength_fraction": args.family_support_min_strength_fraction,
            "family_support_min_margin": args.family_support_min_margin,
        },
        "counts": {
            "nodes": len(public_result),
            "viral_nodes": len(viral),
            "communities": len(community_summary),
            "graph_isolates": int(public_result["is_graph_isolate"].sum()),
            "viral_graph_isolates": int(viral["is_graph_isolate"].sum()),
            "known_viral_gt_nodes": len(known),
            "bridge_candidates": len(bridges),
            **metadata_stats,
            **edge_stats,
        },
        "interpretation_guardrails": [
            (
                "Leiden community membership is a network-partition result; it does not by "
                "itself establish family membership, function, or homology."
            ),
            (
                "A CAZy-unassigned protein is assigned a directly supported family only when "
                "family-specific neighbour count, strength fraction, and support-margin "
                "thresholds all pass."
            ),
            (
                "Bridge subtypes are topological screening calls, not graph-theory articulation "
                "points or evidence of evolutionary intermediacy."
            ),
            (
                "Full-length versus domain-mediated bridge hypotheses require Foldseek alignment "
                "coordinates or an equivalent domain-resolved analysis."
            ),
            (
                "Virus-enriched annotation-depleted communities require structural, domain, and "
                "genomic-context validation before any novel-family claim."
            ),
        ],
    }
    (outdir / "classification_run.json").write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def validate_args(args: argparse.Namespace) -> None:
    for name in [
        "known_family_min_purity",
        "mixed_min_second_family_fraction",
        "peripheral_quantile",
        "virus_enrichment_q",
        "bridge_min_external_strength_fraction",
        "bridge_min_participation",
        "family_support_min_strength_fraction",
        "family_support_min_margin",
    ]:
        validate_probability(name.replace("_", "-"), getattr(args, name))
    if args.known_family_min_annotated < 1:
        raise ValueError("--known-family-min-annotated must be at least 1")
    if args.mixed_min_second_family_n < 1:
        raise ValueError("--mixed-min-second-family-n must be at least 1")
    if args.small_community_max < 2:
        raise ValueError("--small-community-max must be at least 2")
    if args.large_community_min <= args.small_community_max:
        raise ValueError("--large-community-min must exceed --small-community-max")
    if args.bridge_min_degree < 1 or args.bridge_min_external_neighbors < 1:
        raise ValueError("Bridge degree thresholds must be at least 1")
    if args.bridge_hub_min_external_communities < 2:
        raise ValueError("--bridge-hub-min-external-communities must be at least 2")
    if args.family_support_min_neighbors < 1:
        raise ValueError("--family-support-min-neighbors must be at least 1")


def main() -> None:
    args = parse_args()
    validate_args(args)
    nodes = prepare_nodes(args)
    nodes, metadata_unmatched, metadata_stats = merge_viral_reference_metadata(
        nodes, args
    )
    edges, edge_stats = prepare_edges(args, set(nodes[args.node_id_column]))
    result, community_summary, external_profile = compute_metrics(
        nodes, edges, args
    )
    result = classify_nodes(
        result, community_summary, external_profile, args
    )
    write_outputs(
        result,
        community_summary,
        metadata_unmatched,
        metadata_stats,
        edge_stats,
        args,
    )

    viral_n = int(result["_is_virus"].sum())
    known_n = int((result["_is_virus"] & result["is_known_viral_gt"]).sum())
    bridge_n = int((result["_is_virus"] & result["bridge_candidate"]).sum())
    print(
        f"Classified {viral_n:,} viral nodes; known viral GTs: {known_n:,}; "
        f"bridge candidates: {bridge_n:,}"
    )
    print(f"Outputs written to: {Path(args.outdir).resolve()}")


if __name__ == "__main__":
    main()
