#!/usr/bin/env python3
"""Reproducible weighted-Leiden analysis for a Foldseek structure network.

The pipeline:
  1. audits/corrects viral origin metadata without overwriting source columns;
  2. removes self hits and collapses directional Foldseek rows to undirected edges;
  3. compares unweighted, alnTM-weighted, and min(qTM,tTM)-weighted Leiden;
  4. quantifies CAZy recapitulation, graph morphology, heterogeneous-community
     transitions, and viral-protein distributions;
  5. exports Cytoscape-ready node and edge tables plus audit tables.

Weighted modes are enabled only when the alignment table contains the required
Foldseek columns. A standard 12-column m8 still runs the unweighted baseline.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ssn_leiden_matplotlib")

try:
    import numpy as np
    import pandas as pd
    import igraph as ig
    import leidenalg as la
    from scipy.stats import hypergeom
    from sklearn.metrics import (
        adjusted_rand_score,
        completeness_score,
        homogeneity_score,
        normalized_mutual_info_score,
        v_measure_score,
    )
except ImportError as exc:  # pragma: no cover - dependency failure path
    missing = getattr(exc, "name", str(exc))
    raise SystemExit(
        f"Missing Python dependency: {missing}\n"
        "Install dependencies with:\n"
        "  python -m pip install pandas numpy scipy scikit-learn "
        "python-igraph leidenalg matplotlib seaborn"
    ) from exc


STANDARD_M8_COLUMNS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
]

RICH_M8_COLUMNS = [
    "query", "target", "fident", "alnlen", "evalue", "bits",
    "alntmscore", "qtmscore", "ttmscore", "qcov", "tcov", "qlen", "tlen",
]

NUMERIC_ALIGNMENT_COLUMNS = {
    "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend",
    "tstart", "tend", "evalue", "bits", "alntmscore", "qtmscore",
    "ttmscore", "qcov", "tcov", "qlen", "tlen", "lddt",
}

VIRUS_PATTERN = re.compile(
    r"virus|viral|phage|bacteriophage|virophage|viridae|virales|viricetes|viricota|provirus|prophage",
    flags=re.IGNORECASE,
)

TRUE_VALUES = {"1", "true", "t", "yes", "y", "virus", "viral"}
FALSE_VALUES = {"0", "false", "f", "no", "n", "cellular", "nonvirus", "non-virus"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit metadata and compare weighted Leiden partitions of a Foldseek SSN.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--nodes", required=True, help="Node CSV/TSV.")
    parser.add_argument("--alignments", required=True, help="Foldseek m8/TSV, preferably rich m8.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--node-id-column", default="name")
    parser.add_argument("--group-column", default="Group")
    parser.add_argument("--species-column", default="Species")
    parser.add_argument("--gtfam-column", default="GTFam")
    parser.add_argument(
        "--alignment-columns",
        default=None,
        help="Comma-separated column names for a headerless nonstandard alignment table.",
    )
    parser.add_argument(
        "--virus-overrides",
        default=None,
        help="Optional CSV/TSV with name,is_virus,Group_corrected,note. Overrides take precedence.",
    )
    parser.add_argument(
        "--focus-clusters",
        default=None,
        help="Optional CSV/TSV with name,focus_cluster for manually selected Cytoscape regions.",
    )
    parser.add_argument("--resolutions", default="0.5,0.75,1.0,1.25,1.5,2.0")
    parser.add_argument("--seeds", default="1,2,3,4,5")
    parser.add_argument("--primary-resolution", type=float, default=1.0)
    parser.add_argument("--n-iterations", type=int, default=-1)
    parser.add_argument("--reciprocal-only", action="store_true")
    parser.add_argument(
        "--analysis-min-coverage",
        type=float,
        default=None,
        help="Optional additional mutual-coverage filter; input Foldseek filtering is otherwise retained.",
    )
    parser.add_argument(
        "--analysis-min-alntmscore",
        type=float,
        default=None,
        help="Optional additional alnTM filter for rich m8.",
    )
    parser.add_argument("--heterogeneous-purity", type=float, default=0.70)
    parser.add_argument("--high-purity", type=float, default=0.80)
    parser.add_argument("--min-annotated", type=int, default=5)
    parser.add_argument("--skip-plots", action="store_true")
    return parser.parse_args()


def log(message: str) -> None:
    print(f"[SSN] {message}", flush=True)


def parse_number_list(text: str, cast) -> list:
    values = [cast(x.strip()) for x in text.split(",") if x.strip()]
    if not values:
        raise ValueError("Expected at least one comma-separated value")
    return values


def resolution_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def read_delimited(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep=None,
        engine="python",
        dtype=str,
        keep_default_na=False,
    )


def read_alignment_table(path: str | Path, explicit_columns: str | None) -> pd.DataFrame:
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        first = next((line.rstrip("\n\r") for line in handle if line.strip()), "")
    if not first:
        raise ValueError(f"Alignment file is empty: {path}")

    fields = first.split("\t")
    normalized = {x.strip().lower() for x in fields}
    has_header = "query" in normalized and "target" in normalized

    if has_header:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        df.columns = [str(c).strip().lower() for c in df.columns]
    else:
        if explicit_columns:
            columns = [x.strip().lower() for x in explicit_columns.split(",") if x.strip()]
        elif len(fields) == len(STANDARD_M8_COLUMNS):
            columns = STANDARD_M8_COLUMNS
        elif len(fields) == len(RICH_M8_COLUMNS):
            columns = RICH_M8_COLUMNS
        else:
            raise ValueError(
                f"Headerless alignment has {len(fields)} columns. Supply --alignment-columns."
            )
        if len(columns) != len(fields):
            raise ValueError(
                f"--alignment-columns supplies {len(columns)} names but file has {len(fields)} columns"
            )
        df = pd.read_csv(path, sep="\t", header=None, names=columns, dtype=str, keep_default_na=False)

    required = {"query", "target"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Alignment table lacks required columns: {sorted(missing)}")

    df["query"] = df["query"].astype(str).str.strip()
    df["target"] = df["target"].astype(str).str.strip()
    for column in NUMERIC_ALIGNMENT_COLUMNS.intersection(df.columns):
        df[column] = pd.to_numeric(df[column], errors="coerce")

    for column in ("qcov", "tcov"):
        if column in df.columns and df[column].dropna().gt(1.0).any():
            plausible_percent = df[column].dropna().between(0, 100).all()
            if plausible_percent:
                warnings.warn(f"{column} appears to be percent; dividing by 100")
                df[column] = df[column] / 100.0
    return df


def parse_optional_bool(value: object) -> bool | None:
    text = str(value).strip().lower()
    if not text:
        return None
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise ValueError(f"Cannot parse boolean value: {value!r}")


def column_lookup(df: pd.DataFrame) -> dict[str, str]:
    return {str(c).strip().lower(): c for c in df.columns}


def audit_and_correct_nodes(
    node_path: str | Path,
    node_id_column: str,
    group_column: str,
    species_column: str,
    gtfam_column: str,
    override_path: str | Path | None,
    focus_path: str | Path | None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    nodes = read_delimited(node_path)
    required = [node_id_column]
    for column in required:
        if column not in nodes.columns:
            raise ValueError(f"Node table lacks required column: {column}")

    for column in (group_column, species_column, gtfam_column):
        if column not in nodes.columns:
            nodes[column] = ""

    nodes[node_id_column] = nodes[node_id_column].astype(str).str.strip()
    if nodes[node_id_column].eq("").any():
        raise ValueError("Node IDs contain empty values")
    if nodes[node_id_column].duplicated().any():
        dup = nodes.loc[nodes[node_id_column].duplicated(), node_id_column].head().tolist()
        raise ValueError(f"Node IDs are not unique; examples: {dup}")

    original_group = nodes[group_column].fillna("").astype(str).str.strip()
    species = nodes[species_column].fillna("").astype(str).str.strip()
    from_group = original_group.str.casefold().eq("virus")
    from_species = species.str.contains(VIRUS_PATTERN, na=False)

    nodes["Group_original"] = original_group
    nodes["Species_original"] = species
    nodes["virus_evidence_group"] = from_group
    nodes["virus_evidence_species_keyword"] = from_species
    nodes["is_virus"] = from_group | from_species
    nodes["virus_annotation_source"] = np.select(
        [from_group & from_species, from_group, from_species],
        ["group+species", "original_group", "species_keyword"],
        default="no_viral_evidence",
    )
    nodes["Group_corrected"] = np.where(
        nodes["is_virus"],
        "Virus",
        original_group.replace("", "Unknown"),
    )
    nodes["virus_override_applied"] = False
    nodes["virus_override_note"] = ""

    warnings_list: list[str] = []
    if override_path:
        overrides = read_delimited(override_path)
        lookup = column_lookup(overrides)
        id_key = node_id_column.lower()
        if id_key not in lookup:
            if "name" in lookup:
                id_col = lookup["name"]
            else:
                raise ValueError("Virus override table requires a name/node ID column")
        else:
            id_col = lookup[id_key]
        override_map = nodes.set_index(node_id_column).index
        for _, row in overrides.iterrows():
            node_id = str(row[id_col]).strip()
            if not node_id:
                continue
            if node_id not in override_map:
                warnings_list.append(f"Virus override ID not found in node table: {node_id}")
                continue
            mask = nodes[node_id_column].eq(node_id)
            if "is_virus" in lookup:
                parsed = parse_optional_bool(row[lookup["is_virus"]])
                if parsed is not None:
                    nodes.loc[mask, "is_virus"] = parsed
            if "group_corrected" in lookup:
                value = str(row[lookup["group_corrected"]]).strip()
                if value:
                    nodes.loc[mask, "Group_corrected"] = value
            if "note" in lookup:
                nodes.loc[mask, "virus_override_note"] = str(row[lookup["note"]]).strip()
            nodes.loc[mask, "virus_override_applied"] = True
            nodes.loc[mask, "virus_annotation_source"] = "manual_override"
            if bool(nodes.loc[mask, "is_virus"].iloc[0]):
                nodes.loc[mask, "Group_corrected"] = "Virus"

    nodes["Group_annotation_changed"] = (
        nodes["Group_original"].replace("", "Unknown") != nodes["Group_corrected"]
    )
    nodes["virus_metadata_needs_review"] = (
        (nodes["virus_evidence_group"] & ~nodes["virus_evidence_species_keyword"])
        | nodes["Group_annotation_changed"]
    ) & ~nodes["virus_override_applied"]

    gtfam = nodes[gtfam_column].fillna("").astype(str).str.strip()
    nodes["CAZy_subgroup"] = gtfam
    nodes["CAZy_family"] = gtfam.str.extract(r"^(GT\d+)", flags=re.IGNORECASE, expand=False)
    nodes["CAZy_family"] = nodes["CAZy_family"].fillna("").str.upper()

    if focus_path:
        focus = read_delimited(focus_path)
        lookup = column_lookup(focus)
        if node_id_column.lower() not in lookup and "name" not in lookup:
            raise ValueError("Focus-cluster table requires a name/node ID column")
        if "focus_cluster" not in lookup:
            raise ValueError("Focus-cluster table requires focus_cluster column")
        fid = lookup.get(node_id_column.lower(), lookup.get("name"))
        focus_small = focus[[fid, lookup["focus_cluster"]]].copy()
        focus_small.columns = [node_id_column, "focus_cluster"]
        focus_small[node_id_column] = focus_small[node_id_column].astype(str).str.strip()
        nodes = nodes.merge(focus_small, on=node_id_column, how="left")
        nodes["focus_cluster"] = nodes["focus_cluster"].fillna("")
    else:
        nodes["focus_cluster"] = ""

    audit_columns = [
        node_id_column, "Species_original", "Group_original", "Group_corrected",
        "is_virus", "virus_annotation_source", "virus_evidence_group",
        "virus_evidence_species_keyword", "virus_override_applied",
        "virus_override_note", "Group_annotation_changed", "virus_metadata_needs_review",
    ]
    audit = nodes.loc[
        nodes["Group_annotation_changed"] | nodes["virus_metadata_needs_review"] | nodes["is_virus"],
        audit_columns,
    ].copy()
    return nodes, audit, warnings_list


def clean_undirected_edges(
    alignments: pd.DataFrame,
    node_ids: set[str],
    reciprocal_only: bool,
    min_coverage: float | None,
    min_alntmscore: float | None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    stats: dict[str, int] = {"alignment_rows_input": len(alignments)}
    df = alignments.copy()
    known = df["query"].isin(node_ids) & df["target"].isin(node_ids)
    stats["alignment_rows_missing_node"] = int((~known).sum())
    df = df.loc[known].copy()

    stats["self_hits_removed"] = int(df["query"].eq(df["target"]).sum())
    df = df.loc[~df["query"].eq(df["target"])].copy()

    if {"qcov", "tcov"}.issubset(df.columns):
        df["row_mutual_coverage"] = df[["qcov", "tcov"]].min(axis=1)
        if min_coverage is not None:
            before = len(df)
            df = df.loc[df["row_mutual_coverage"] >= min_coverage].copy()
            stats["rows_removed_analysis_coverage"] = before - len(df)
    elif min_coverage is not None:
        raise ValueError("--analysis-min-coverage requires qcov and tcov")

    if "alntmscore" in df.columns and min_alntmscore is not None:
        before = len(df)
        df = df.loc[df["alntmscore"] >= min_alntmscore].copy()
        stats["rows_removed_analysis_alntmscore"] = before - len(df)
    elif min_alntmscore is not None:
        raise ValueError("--analysis-min-alntmscore requires alntmscore")

    if {"qtmscore", "ttmscore"}.issubset(df.columns):
        df["row_global_min_tm"] = df[["qtmscore", "ttmscore"]].min(axis=1)

    q = df["query"].to_numpy(dtype=str)
    t = df["target"].to_numpy(dtype=str)
    df["source"] = np.where(q <= t, q, t)
    df["target_undirected"] = np.where(q <= t, t, q)
    df["direction_token"] = df["query"] + "\t" + df["target"]

    grouped = df.groupby(["source", "target_undirected"], sort=False, observed=True)
    edges = grouped.size().rename("directional_rows").reset_index()
    orientations = grouped["direction_token"].nunique().rename("orientation_count").reset_index(drop=True)
    edges["orientation_count"] = orientations
    edges["reciprocal"] = edges["orientation_count"] >= 2

    aggregations = {
        "evalue": ("evalue", "min"),
        "bits_mean": ("bits", "mean"),
        "bits_max": ("bits", "max"),
        "weight_alntmscore": ("alntmscore", "mean"),
        "weight_global_min_tm": ("row_global_min_tm", "mean"),
        "mutual_coverage_mean": ("row_mutual_coverage", "mean"),
        "mutual_coverage_min": ("row_mutual_coverage", "min"),
        "fident_mean": ("fident", "mean"),
        "alnlen_mean": ("alnlen", "mean"),
    }
    for output_col, (input_col, operation) in aggregations.items():
        if input_col not in df.columns:
            continue
        values = getattr(grouped[input_col], operation)().reset_index(drop=True)
        edges[output_col] = values

    edges = edges.rename(columns={"target_undirected": "target"})
    stats["unique_undirected_edges_before_reciprocal_filter"] = len(edges)
    stats["reciprocal_edges"] = int(edges["reciprocal"].sum())
    stats["one_way_edges"] = int((~edges["reciprocal"]).sum())
    if reciprocal_only:
        edges = edges.loc[edges["reciprocal"]].copy()
    stats["unique_undirected_edges_analysis"] = len(edges)
    return edges.reset_index(drop=True), stats


def build_graph(nodes: pd.DataFrame, edges: pd.DataFrame, node_id_column: str) -> tuple[ig.Graph, np.ndarray, np.ndarray]:
    node_names = nodes[node_id_column].astype(str).tolist()
    index = {name: i for i, name in enumerate(node_names)}
    source_idx = edges["source"].map(index).to_numpy(dtype=int)
    target_idx = edges["target"].map(index).to_numpy(dtype=int)
    graph = ig.Graph(n=len(node_names), edges=list(zip(source_idx, target_idx)), directed=False)
    graph.vs["name"] = node_names
    return graph, source_idx, target_idx


def relabel_membership_by_size(membership: Sequence[int]) -> tuple[np.ndarray, dict[int, str]]:
    membership = np.asarray(membership, dtype=int)
    counts = pd.Series(membership).value_counts()
    ordered = sorted(counts.index.tolist(), key=lambda x: (-counts[x], int(x)))
    mapping = {old: f"C{rank:03d}" for rank, old in enumerate(ordered, start=1)}
    labels = np.asarray([mapping[x] for x in membership], dtype=object)
    return labels, mapping


def pairwise_ari(partitions: Sequence[Sequence[int]]) -> tuple[float, float]:
    if len(partitions) < 2:
        return 1.0, 1.0
    values = [adjusted_rand_score(a, b) for a, b in itertools.combinations(partitions, 2)]
    return float(np.mean(values)), float(np.min(values))


def label_agreement_metrics(membership: Sequence[object], labels: pd.Series) -> dict[str, float]:
    comm = np.asarray(membership, dtype=object)
    label_array = labels.fillna("").astype(str).to_numpy()
    keep = label_array != ""
    if keep.sum() == 0:
        return {
            "annotated_n": 0, "purity": np.nan, "recovery": np.nan,
            "homogeneity": np.nan, "completeness": np.nan,
            "v_measure": np.nan, "ari": np.nan, "nmi": np.nan,
            "median_family_community_count": np.nan,
        }
    c = comm[keep]
    y = label_array[keep]
    table = pd.crosstab(pd.Series(c, name="community"), pd.Series(y, name="label"))
    purity = table.max(axis=1).sum() / table.to_numpy().sum()
    recovery = table.max(axis=0).sum() / table.to_numpy().sum()
    family_fragmentation = (table > 0).sum(axis=0)
    return {
        "annotated_n": int(keep.sum()),
        "purity": float(purity),
        "recovery": float(recovery),
        "homogeneity": float(homogeneity_score(y, c)),
        "completeness": float(completeness_score(y, c)),
        "v_measure": float(v_measure_score(y, c)),
        "ari": float(adjusted_rand_score(y, c)),
        "nmi": float(normalized_mutual_info_score(y, c)),
        "median_family_community_count": float(family_fragmentation.median()),
    }


def partition_summary(
    membership: Sequence[object],
    nodes: pd.DataFrame,
) -> dict[str, float]:
    membership = np.asarray(membership, dtype=object)
    sizes = pd.Series(membership).value_counts()
    family = label_agreement_metrics(membership, nodes["CAZy_family"])
    subgroup = label_agreement_metrics(membership, nodes["CAZy_subgroup"])
    result: dict[str, float] = {
        "n_communities": int(len(sizes)),
        "singleton_communities": int((sizes == 1).sum()),
        "largest_community_n": int(sizes.max()),
    }
    result.update({f"family_{k}": v for k, v in family.items()})
    result.update({f"subgroup_{k}": v for k, v in subgroup.items()})
    if "is_giant_component" in nodes.columns:
        giant = nodes["is_giant_component"].astype(bool).to_numpy()
        giant_membership = membership[giant]
        giant_family = label_agreement_metrics(
            giant_membership, nodes.loc[giant, "CAZy_family"].reset_index(drop=True)
        )
        giant_subgroup = label_agreement_metrics(
            giant_membership, nodes.loc[giant, "CAZy_subgroup"].reset_index(drop=True)
        )
        result["giant_component_n"] = int(giant.sum())
        result["giant_n_communities"] = int(pd.Series(giant_membership).nunique())
        result.update({f"giant_family_{k}": v for k, v in giant_family.items()})
        result.update({f"giant_subgroup_{k}": v for k, v in giant_subgroup.items()})
    return result


def run_leiden_grid(
    graph: ig.Graph,
    nodes: pd.DataFrame,
    methods: dict[str, np.ndarray | None],
    resolutions: Sequence[float],
    seeds: Sequence[int],
    n_iterations: int,
) -> tuple[pd.DataFrame, dict[tuple[str, float], np.ndarray]]:
    summaries: list[dict] = []
    best_memberships: dict[tuple[str, float], np.ndarray] = {}
    for method, weights in methods.items():
        log(f"Leiden method={method}")
        weights_list = None if weights is None else np.asarray(weights, dtype=float).tolist()
        for resolution in resolutions:
            memberships: list[np.ndarray] = []
            qualities: list[float] = []
            for seed in seeds:
                partition = la.find_partition(
                    graph,
                    la.RBConfigurationVertexPartition,
                    weights=weights_list,
                    resolution_parameter=float(resolution),
                    n_iterations=n_iterations,
                    seed=int(seed),
                )
                memberships.append(np.asarray(partition.membership, dtype=int))
                qualities.append(float(partition.quality()))
            best_index = int(np.argmax(qualities))
            best_raw = memberships[best_index]
            best_labels, _ = relabel_membership_by_size(best_raw)
            best_memberships[(method, float(resolution))] = best_labels
            stability_mean, stability_min = pairwise_ari(memberships)
            row = {
                "method": method,
                "resolution": float(resolution),
                "best_seed": int(seeds[best_index]),
                "leiden_quality": float(qualities[best_index]),
                "seed_stability_ari_mean": stability_mean,
                "seed_stability_ari_min": stability_min,
            }
            row.update(partition_summary(best_labels, nodes))
            summaries.append(row)
    return pd.DataFrame(summaries), best_memberships


def fdr_bh(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
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
    p = counts / counts.sum()
    return float(-(p * np.log(p)).sum() / np.log(len(counts)))


def community_morphology(
    method: str,
    membership: np.ndarray,
    nodes: pd.DataFrame,
    source_idx: np.ndarray,
    target_idx: np.ndarray,
    weights: np.ndarray | None,
) -> pd.DataFrame:
    edge_weights = np.ones(len(source_idx), dtype=float) if weights is None else np.asarray(weights, dtype=float)
    communities = pd.Series(membership).drop_duplicates().tolist()
    node_groups = pd.Series(np.arange(len(nodes))).groupby(membership).apply(list).to_dict()
    internal_count = defaultdict(int)
    external_count = defaultdict(int)
    internal_strength = defaultdict(float)
    external_strength = defaultdict(float)
    for s, t, weight in zip(source_idx, target_idx, edge_weights):
        cs, ct = membership[s], membership[t]
        if cs == ct:
            internal_count[cs] += 1
            internal_strength[cs] += float(weight)
        else:
            external_count[cs] += 1
            external_count[ct] += 1
            external_strength[cs] += float(weight)
            external_strength[ct] += float(weight)

    total_nodes = len(nodes)
    total_viral = int(nodes["is_virus"].astype(bool).sum())
    rows: list[dict] = []
    for community in communities:
        idx = node_groups[community]
        subset = nodes.iloc[idx]
        n = len(idx)
        families = subset.loc[subset["CAZy_family"] != "", "CAZy_family"].value_counts()
        annotated_n = int(families.sum())
        dominant = families.index[0] if annotated_n else ""
        dominant_n = int(families.iloc[0]) if annotated_n else 0
        purity = dominant_n / annotated_n if annotated_n else np.nan
        virus_n = int(subset["is_virus"].astype(bool).sum())
        possible_edges = n * (n - 1) / 2
        density = internal_count[community] / possible_edges if possible_edges else 0.0
        denominator = 2 * internal_strength[community] + external_strength[community]
        conductance = external_strength[community] / denominator if denominator else 0.0
        pvalue = float(hypergeom.sf(virus_n - 1, total_nodes, total_viral, n)) if total_viral else np.nan
        rows.append({
            "method": method,
            "community": community,
            "n_nodes": n,
            "internal_edges": internal_count[community],
            "external_edges": external_count[community],
            "internal_strength": internal_strength[community],
            "external_strength": external_strength[community],
            "internal_external_strength_ratio": internal_strength[community] / (external_strength[community] + 1e-12),
            "density": density,
            "conductance": conductance,
            "CAZy_annotated_n": annotated_n,
            "CAZy_family_count": int(len(families)),
            "dominant_CAZy_family": dominant,
            "dominant_CAZy_n": dominant_n,
            "CAZy_purity": purity,
            "CAZy_normalized_entropy": normalized_entropy(families),
            "virus_n": virus_n,
            "virus_fraction": virus_n / n,
            "unassigned_virus_n": int((subset["is_virus"].astype(bool) & subset["CAZy_family"].eq("")).sum()),
            "virus_enrichment_p": pvalue,
        })
    result = pd.DataFrame(rows).sort_values(["n_nodes", "community"], ascending=[False, True])
    result["virus_enrichment_q"] = fdr_bh(result["virus_enrichment_p"].to_numpy())
    return result


def node_position_metrics(
    method: str,
    membership: np.ndarray,
    nodes: pd.DataFrame,
    morphology: pd.DataFrame,
    source_idx: np.ndarray,
    target_idx: np.ndarray,
    weights: np.ndarray | None,
    high_purity: float,
    min_annotated: int,
) -> pd.DataFrame:
    n = len(nodes)
    edge_weights = np.ones(len(source_idx), dtype=float) if weights is None else np.asarray(weights, dtype=float)
    total = np.zeros(n, dtype=float)
    within = np.zeros(n, dtype=float)
    neighbor_comm_strength: list[dict[str, float]] = [defaultdict(float) for _ in range(n)]
    for s, t, weight in zip(source_idx, target_idx, edge_weights):
        weight = float(weight)
        total[s] += weight
        total[t] += weight
        neighbor_comm_strength[s][membership[t]] += weight
        neighbor_comm_strength[t][membership[s]] += weight
        if membership[s] == membership[t]:
            within[s] += weight
            within[t] += weight
    participation = np.zeros(n, dtype=float)
    for i in range(n):
        if total[i] > 0:
            participation[i] = 1.0 - sum((x / total[i]) ** 2 for x in neighbor_comm_strength[i].values())

    result = nodes.copy()
    result["analysis_method"] = method
    result["community"] = membership
    result["weighted_degree"] = total
    result["within_community_strength"] = within
    result["within_strength_fraction"] = np.divide(within, total, out=np.zeros_like(within), where=total > 0)
    result["participation_coefficient"] = participation

    m = morphology.set_index("community")
    result["community_dominant_CAZy"] = result["community"].map(m["dominant_CAZy_family"])
    result["community_CAZy_purity"] = result["community"].map(m["CAZy_purity"])
    result["community_CAZy_annotated_n"] = result["community"].map(m["CAZy_annotated_n"])
    result["community_virus_enrichment_q"] = result["community"].map(m["virus_enrichment_q"])
    result["community_n"] = result["community"].map(m["n_nodes"])

    q25 = result.groupby("community")["within_community_strength"].quantile(0.25)
    provisional: list[str] = []
    for _, row in result.iterrows():
        if not bool(row["is_virus"]):
            provisional.append("")
        elif row["weighted_degree"] == 0 or row["community_n"] == 1:
            provisional.append("singleton_or_isolate")
        elif (
            row["community_CAZy_annotated_n"] >= min_annotated
            and pd.notna(row["community_CAZy_purity"])
            and row["community_CAZy_purity"] >= high_purity
        ):
            threshold = q25.get(row["community"], 0.0)
            if row["within_community_strength"] >= threshold:
                provisional.append("known_family_embedded")
            else:
                provisional.append("peripheral_family_affiliate")
        elif row["community_CAZy_annotated_n"] > 0:
            provisional.append("mixed_structural_neighborhood")
        elif (
            pd.notna(row["community_virus_enrichment_q"])
            and row["community_virus_enrichment_q"] < 0.05
        ):
            provisional.append("virus_enriched_unassigned_community")
        else:
            provisional.append("unassigned_structural_neighborhood")
    result["provisional_viral_class"] = provisional
    return result


def heterogeneous_transitions(
    baseline_membership: np.ndarray,
    target_membership: np.ndarray,
    baseline_morphology: pd.DataFrame,
    target_morphology: pd.DataFrame,
    nodes: pd.DataFrame,
    hetero_purity: float,
    min_annotated: int,
) -> pd.DataFrame:
    baseline_info = baseline_morphology.set_index("community")
    heterogeneous = baseline_morphology.loc[
        (baseline_morphology["CAZy_annotated_n"] >= min_annotated)
        & (baseline_morphology["CAZy_family_count"] >= 2)
        & (baseline_morphology["CAZy_purity"] < hetero_purity),
        "community",
    ].tolist()
    frame = pd.DataFrame({
        "baseline_community": baseline_membership,
        "target_community": target_membership,
        "CAZy_family": nodes["CAZy_family"].to_numpy(),
        "is_virus": nodes["is_virus"].to_numpy(),
    })
    frame = frame.loc[frame["baseline_community"].isin(heterogeneous)].copy()
    if frame.empty:
        return pd.DataFrame()
    result = (
        frame.groupby(["baseline_community", "target_community"], observed=True)
        .agg(
            n=("target_community", "size"),
            virus_n=("is_virus", "sum"),
            CAZy_annotated_n=("CAZy_family", lambda x: int((x != "").sum())),
            CAZy_families=("CAZy_family", lambda x: ";".join(sorted(set(v for v in x if v)))),
        )
        .reset_index()
    )
    totals = result.groupby("baseline_community")["n"].transform("sum")
    result["fraction_of_baseline"] = result["n"] / totals
    target_info = target_morphology.set_index("community")
    result["target_dominant_CAZy"] = result["target_community"].map(target_info["dominant_CAZy_family"])
    result["target_CAZy_purity"] = result["target_community"].map(target_info["CAZy_purity"])
    result["baseline_CAZy_purity"] = result["baseline_community"].map(baseline_info["CAZy_purity"])
    return result.sort_values(["baseline_community", "n"], ascending=[True, False])


def focus_cluster_transitions(
    nodes: pd.DataFrame,
    memberships: dict[str, np.ndarray],
) -> pd.DataFrame:
    if "focus_cluster" not in nodes or nodes["focus_cluster"].eq("").all():
        return pd.DataFrame()
    frame = pd.DataFrame({"focus_cluster": nodes["focus_cluster"]})
    frame = frame.loc[frame["focus_cluster"] != ""].copy()
    for method, membership in memberships.items():
        frame[method] = membership[frame.index]
    grouping = ["focus_cluster"] + list(memberships)
    return frame.groupby(grouping, observed=True).size().rename("n").reset_index()


def plot_summary(summary: pd.DataFrame, outdir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        warnings.warn("matplotlib/seaborn unavailable; skipping plots")
        return
    sns.set_theme(style="whitegrid")
    family_prefix = "giant_family" if "giant_family_purity" in summary.columns else "family"
    metrics = [
        (f"{family_prefix}_purity", "CAZy family purity (giant component)"),
        (f"{family_prefix}_recovery", "CAZy family recovery (giant component)"),
        (f"{family_prefix}_v_measure", "CAZy family V-measure (giant component)"),
        ("seed_stability_ari_mean", "Seed stability (ARI)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for ax, (metric, title) in zip(axes.flat, metrics):
        sns.lineplot(data=summary, x="resolution", y=metric, hue="method", marker="o", ax=ax)
        ax.set_title(title)
        ax.set_ylim(0, 1.02)
    fig.savefig(outdir / "leiden_method_resolution_summary.png", dpi=220)
    plt.close(fig)


def plot_cazy_heatmap(
    method: str,
    membership: np.ndarray,
    nodes: pd.DataFrame,
    outdir: Path,
    max_communities: int = 30,
    max_families: int = 30,
) -> None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        return
    frame = pd.DataFrame({"community": membership, "family": nodes["CAZy_family"]})
    frame = frame.loc[frame["family"] != ""]
    if frame.empty:
        return
    table = pd.crosstab(frame["community"], frame["family"])
    top_rows = table.sum(axis=1).nlargest(max_communities).index
    top_cols = table.sum(axis=0).nlargest(max_families).index
    table = table.loc[top_rows, top_cols]
    normalized = table.div(table.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig_w = max(9, 0.35 * len(top_cols) + 4)
    fig_h = max(6, 0.28 * len(top_rows) + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
    sns.heatmap(normalized, cmap="mako", vmin=0, vmax=1, ax=ax)
    ax.set_title(f"{method}: community × CAZy family (row-normalized)")
    fig.savefig(outdir / f"cazy_heatmap_{method}.png", dpi=220)
    plt.close(fig)


def write_override_template(path: Path, node_id_column: str) -> None:
    template = pd.DataFrame(columns=[node_id_column, "is_virus", "Group_corrected", "note"])
    template.to_csv(path, index=False)


def write_run_report(
    path: Path,
    args: argparse.Namespace,
    stats: dict[str, int],
    methods: Sequence[str],
    nodes: pd.DataFrame,
    warnings_list: Sequence[str],
) -> None:
    lines = [
        "# Structure-similarity network analysis run\n",
        f"- Nodes: {len(nodes):,}",
        f"- Corrected viral-origin nodes: {int(nodes['is_virus'].sum()):,}",
        f"- Group annotations changed automatically/manually: {int(nodes['Group_annotation_changed'].sum()):,}",
        f"- Viral metadata entries flagged for review: {int(nodes['virus_metadata_needs_review'].sum()):,}",
        f"- Analysis methods: {', '.join(methods)}",
        f"- Primary resolution: {args.primary_resolution:g}",
        "",
        "## Edge cleaning",
    ]
    lines.extend(f"- {key}: {value:,}" for key, value in stats.items())
    lines.extend([
        "",
        "## Interpretation guardrails",
        "- CAZy labels were used only for post-hoc evaluation, not Leiden optimization.",
        "- Group_original is preserved; Group_corrected and is_virus are auditable derived fields.",
        "- provisional_viral_class is a screening label and requires biological review.",
        "- Cytoscape layout appearance is not used as a clustering metric.",
    ])
    if warnings_list:
        lines.extend(["", "## Warnings"])
        lines.extend(f"- {message}" for message in warnings_list)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    resolutions = parse_number_list(args.resolutions, float)
    seeds = parse_number_list(args.seeds, int)
    if args.primary_resolution not in resolutions:
        resolutions.append(args.primary_resolution)
        resolutions = sorted(set(resolutions))

    log("Auditing and correcting node metadata")
    nodes, audit, warnings_list = audit_and_correct_nodes(
        args.nodes,
        args.node_id_column,
        args.group_column,
        args.species_column,
        args.gtfam_column,
        args.virus_overrides,
        args.focus_clusters,
    )
    nodes.to_csv(outdir / "nodes_metadata_corrected.csv", index=False)
    audit.to_csv(outdir / "virus_metadata_audit.csv", index=False)
    write_override_template(outdir / "virus_overrides_template.csv", args.node_id_column)

    log("Reading Foldseek alignments and cleaning undirected edges")
    alignments = read_alignment_table(args.alignments, args.alignment_columns)
    edges, edge_stats = clean_undirected_edges(
        alignments,
        set(nodes[args.node_id_column]),
        args.reciprocal_only,
        args.analysis_min_coverage,
        args.analysis_min_alntmscore,
    )
    edges.to_csv(outdir / "edges_clean_undirected.csv", index=False)

    graph, source_idx, target_idx = build_graph(nodes, edges, args.node_id_column)
    raw_components = np.asarray(graph.connected_components(mode="weak").membership, dtype=int)
    component_labels, _ = relabel_membership_by_size(raw_components)
    nodes["connected_component"] = component_labels
    component_sizes = pd.Series(component_labels).value_counts()
    nodes["connected_component_size"] = pd.Series(component_labels).map(component_sizes).to_numpy()
    nodes["is_giant_component"] = component_labels == "C001"
    # Re-save after graph-derived component metadata have been added.
    nodes.to_csv(outdir / "nodes_metadata_corrected.csv", index=False)
    methods: dict[str, np.ndarray | None] = {"unweighted": None}
    if "weight_alntmscore" in edges and edges["weight_alntmscore"].notna().all():
        methods["alntmscore"] = edges["weight_alntmscore"].to_numpy(dtype=float)
    else:
        warnings_list.append("alntmscore-weighted Leiden skipped: rich m8 column unavailable/incomplete")
    if "weight_global_min_tm" in edges and edges["weight_global_min_tm"].notna().all():
        methods["global_min_tm"] = edges["weight_global_min_tm"].to_numpy(dtype=float)
    else:
        warnings_list.append("min(qTM,tTM)-weighted Leiden skipped: q/ttmscore unavailable/incomplete")

    log(f"Running methods: {', '.join(methods)}")
    summary, best_memberships = run_leiden_grid(
        graph, nodes, methods, resolutions, seeds, args.n_iterations
    )
    summary.to_csv(outdir / "leiden_resolution_summary.csv", index=False)

    assignments = nodes.copy()
    for (method, resolution), membership in best_memberships.items():
        assignments[f"community_{method}_r{resolution_label(resolution)}"] = membership
    assignments.to_csv(outdir / "nodes_all_leiden_assignments.csv", index=False)

    primary_memberships = {
        method: best_memberships[(method, float(args.primary_resolution))]
        for method in methods
    }
    primary = nodes.copy()
    for method, membership in primary_memberships.items():
        primary[f"community_{method}"] = membership

    agreement_rows = []
    for (method_a, part_a), (method_b, part_b) in itertools.combinations(primary_memberships.items(), 2):
        agreement_rows.append({
            "method_a": method_a,
            "method_b": method_b,
            "ARI": adjusted_rand_score(part_a, part_b),
            "NMI": normalized_mutual_info_score(part_a, part_b),
        })
    pd.DataFrame(agreement_rows).to_csv(outdir / "primary_partition_agreement.csv", index=False)

    morphology_by_method: dict[str, pd.DataFrame] = {}
    viral_tables: list[pd.DataFrame] = []
    for method, membership in primary_memberships.items():
        morphology = community_morphology(
            method, membership, nodes, source_idx, target_idx, methods[method]
        )
        morphology_by_method[method] = morphology
        morphology.to_csv(outdir / f"community_morphology_{method}.csv", index=False)
        viral = node_position_metrics(
            method, membership, nodes, morphology, source_idx, target_idx,
            methods[method], args.high_purity, args.min_annotated,
        )
        viral = viral.loc[viral["is_virus"].astype(bool)].copy()
        viral.to_csv(outdir / f"viral_distribution_{method}.csv", index=False)
        viral_tables.append(viral)
        if not args.skip_plots:
            plot_cazy_heatmap(method, membership, nodes, outdir)

    if "unweighted" in primary_memberships:
        for method in primary_memberships:
            if method == "unweighted":
                continue
            transitions = heterogeneous_transitions(
                primary_memberships["unweighted"], primary_memberships[method],
                morphology_by_method["unweighted"], morphology_by_method[method],
                nodes, args.heterogeneous_purity, args.min_annotated,
            )
            transitions.to_csv(
                outdir / f"heterogeneous_transitions_unweighted_to_{method}.csv", index=False
            )

    focus = focus_cluster_transitions(nodes, primary_memberships)
    if not focus.empty:
        focus.to_csv(outdir / "manual_focus_cluster_transitions.csv", index=False)

    preferred = "global_min_tm" if "global_min_tm" in methods else (
        "alntmscore" if "alntmscore" in methods else "unweighted"
    )
    preferred_viral = node_position_metrics(
        preferred,
        primary_memberships[preferred],
        nodes,
        morphology_by_method[preferred],
        source_idx,
        target_idx,
        methods[preferred],
        args.high_purity,
        args.min_annotated,
    )
    cytoscape_nodes = preferred_viral.copy()
    for method, membership in primary_memberships.items():
        cytoscape_nodes[f"community_{method}"] = membership
    cytoscape_nodes.to_csv(outdir / "cytoscape_nodes_analyzed.csv", index=False)
    edges.to_csv(outdir / "cytoscape_edges_analyzed.csv", index=False)

    if not args.skip_plots:
        plot_summary(summary, outdir)

    config = {
        "arguments": vars(args),
        "resolutions": resolutions,
        "seeds": seeds,
        "methods_run": list(methods),
        "preferred_method_for_cytoscape_node_metrics": preferred,
        "edge_stats": edge_stats,
    }
    (outdir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    write_run_report(
        outdir / "RUN_SUMMARY.md", args, edge_stats, list(methods), nodes, warnings_list
    )
    log(f"Complete. Outputs written to {outdir.resolve()}")


if __name__ == "__main__":
    main()
