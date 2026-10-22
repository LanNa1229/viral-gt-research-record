from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .common import (
    load_structure_residues,
    map_alignment_sequence_to_structure,
    normalized_entropy_conservation,
    read_fasta,
    read_tsv,
    require_columns,
    write_metric_pdb,
)
from .constraints import aggregate_local_lddt, contiguous_groups


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA")


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def _resolve_alignment_to_protein_ids(
    path: str | Path,
    members: pd.DataFrame,
    required_protein_ids: Iterable[str],
) -> dict[str, str]:
    records = read_fasta(path)
    required = [str(x) for x in required_protein_ids]
    by_protein: dict[str, list[str]] = {}
    for row in members.itertuples(index=False):
        protein_id = str(row.protein_id)
        candidates = [protein_id]
        if hasattr(row, "foldmason_id") and str(row.foldmason_id):
            candidates.append(str(row.foldmason_id))
        by_protein[protein_id] = candidates

    resolved: dict[str, str] = {}
    used_headers: set[str] = set()
    for protein_id in required:
        candidates = by_protein.get(protein_id, [protein_id])
        hits = [candidate for candidate in candidates if candidate in records]
        if not hits:
            normalized_hits = [
                header
                for header in records
                if Path(header).stem in candidates
                or Path(header).name in candidates
                or any(header.startswith(f"{candidate}__") for candidate in candidates)
            ]
            hits = normalized_hits
        if len(hits) != 1:
            raise ValueError(
                f"Could not uniquely resolve {protein_id} in {path}; candidates={candidates}, hits={hits}"
            )
        header = hits[0]
        if header in used_headers:
            raise ValueError(f"Alignment header {header} resolved to more than one protein")
        resolved[protein_id] = records[header]
        used_headers.add(header)
    return resolved


def _prepare_members(members_path: str | Path, medoids_path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    members = read_tsv(members_path)
    medoids = read_tsv(medoids_path)
    require_columns(members, ["protein_id", "cluster_id"], str(members_path))
    require_columns(medoids, ["cluster_id", "medoid_id"], str(medoids_path))

    members = members.copy()
    if "exclude" in members:
        members = members[~_as_bool(members["exclude"])].copy()
    structure_col = next(
        (x for x in ("staged_structure_path", "structure_path", "pdb_path", "cif_path") if x in members),
        None,
    )
    if structure_col is None:
        raise ValueError("members table needs staged_structure_path, structure_path, pdb_path, or cif_path")
    members["_structure_path"] = members[structure_col].astype(str)
    if "chain_id" not in members:
        members["chain_id"] = ""

    medoids = medoids[["cluster_id", "medoid_id"]].copy()
    medoids["cluster_id"] = medoids["cluster_id"].astype(str)
    medoids["medoid_id"] = medoids["medoid_id"].astype(str)
    member_ids = set(members["protein_id"].astype(str))
    missing = sorted(set(medoids["medoid_id"]) - member_ids)
    if missing:
        raise ValueError(f"Medoids absent from members table: {missing}")
    represented = set(medoids["cluster_id"])
    members = members[members["cluster_id"].astype(str).isin(represented)].copy()
    counts = members.groupby("cluster_id").size()
    too_small = counts[counts < 2]
    if not too_small.empty:
        raise ValueError(f"Every pan-core subcluster needs at least two members: {too_small.to_dict()}")
    return members.reset_index(drop=True), medoids.sort_values("cluster_id").reset_index(drop=True)


def _build_maps(
    aa: dict[str, str],
    manifest: pd.DataFrame,
    config: dict,
) -> tuple[dict[str, dict[int, object]], dict[str, list[object]], pd.DataFrame]:
    maps: dict[str, dict[int, object]] = {}
    residues_by_id: dict[str, list[object]] = {}
    rows = []
    q = config["quality"]
    manifest_by_id = manifest.set_index("protein_id")
    for protein_id, aligned_sequence in aa.items():
        row = manifest_by_id.loc[protein_id]
        _, residues = load_structure_residues(
            row["_structure_path"], str(row.get("chain_id", "")) or None
        )
        column_map, qc = map_alignment_sequence_to_structure(
            aligned_sequence,
            residues,
            float(q["min_sequence_mapping_identity"]),
            float(q["min_sequence_mapping_coverage"]),
        )
        maps[protein_id] = column_map
        residues_by_id[protein_id] = residues
        rows.append({"protein_id": protein_id, **qc})
    return maps, residues_by_id, pd.DataFrame(rows)


def _cluster_column_metrics(
    column: int,
    cluster_ids: list[str],
    members_by_cluster: dict[str, list[str]],
    aa: dict[str, str],
    di: dict[str, str],
) -> dict:
    occupancies = []
    within_di_conservation = []
    within_di_consensus_fraction = []
    cluster_consensus_states = []
    for cluster_id in cluster_ids:
        ids = members_by_cluster[cluster_id]
        aa_states = [aa[x][column] for x in ids if aa[x][column] not in "-."]
        di_states = [di[x][column] for x in ids if di[x][column] not in "-."]
        occupancies.append(len(aa_states) / len(ids))
        stats = normalized_entropy_conservation(di_states, 20)
        within_di_conservation.append(stats["conservation"])
        within_di_consensus_fraction.append(stats["dominant_fraction"])
        if stats["dominant"]:
            cluster_consensus_states.append(stats["dominant"])

    all_ids = [x for ids in members_by_cluster.values() for x in ids]
    pan_aa = normalized_entropy_conservation(
        [aa[x][column] for x in all_ids if aa[x][column] not in "-."], 20
    )
    pan_di = normalized_entropy_conservation(
        [di[x][column] for x in all_ids if di[x][column] not in "-."], 20
    )
    cluster_di = normalized_entropy_conservation(cluster_consensus_states, 20)

    def finite_min(values: list[float]) -> float:
        finite = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
        return float(finite.min()) if len(finite) else np.nan

    def finite_mean(values: list[float]) -> float:
        finite = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
        return float(finite.mean()) if len(finite) else np.nan

    return {
        "min_subcluster_occupancy": float(np.min(occupancies)),
        "mean_subcluster_occupancy": float(np.mean(occupancies)),
        "min_subcluster_3di_conservation": finite_min(
            within_di_conservation
        ),
        "mean_subcluster_3di_conservation": finite_mean(
            within_di_conservation
        ),
        "min_subcluster_3di_consensus_fraction": finite_min(
            within_di_consensus_fraction
        ),
        "pan_aa_conservation": pan_aa["conservation"],
        "pan_aa_consensus": pan_aa["dominant"],
        "pan_aa_consensus_fraction": pan_aa["dominant_fraction"],
        "pan_aa_n": pan_aa["n"],
        "pan_3di_conservation": pan_di["conservation"],
        "pan_3di_consensus": pan_di["dominant"],
        "pan_3di_consensus_fraction": pan_di["dominant_fraction"],
        "pan_3di_n": pan_di["n"],
        "cluster_consensus_3di_conservation": cluster_di["conservation"],
        "cluster_consensus_3di": cluster_di["dominant"],
        "cluster_consensus_3di_fraction": cluster_di["dominant_fraction"],
        "n_cluster_consensus_states": cluster_di["n"],
    }


def _anchor_metrics(
    medoids: pd.DataFrame,
    members_by_cluster: dict[str, list[str]],
    maps: dict[str, dict[int, object]],
    config: dict,
) -> pd.DataFrame:
    rows = []
    medoid_ids = medoids["medoid_id"].astype(str).tolist()
    for medoid in medoids.itertuples(index=False):
        cluster_id = str(medoid.cluster_id)
        medoid_id = str(medoid.medoid_id)
        within_ids = [x for x in members_by_cluster[cluster_id] if x != medoid_id]
        cross_ids = [x for x in medoid_ids if x != medoid_id]
        within = aggregate_local_lddt(medoid_id, within_ids, maps, config).rename(
            columns={
                "lddt": "within_lddt",
                "lddt_coverage": "within_lddt_coverage",
                "presence": "within_presence",
                "n_lddt_comparisons": "within_n_lddt_comparisons",
            }
        )
        cross = aggregate_local_lddt(medoid_id, cross_ids, maps, config).rename(
            columns={
                "lddt": "cross_lddt",
                "lddt_coverage": "cross_lddt_coverage",
                "presence": "cross_presence",
                "n_lddt_comparisons": "cross_n_lddt_comparisons",
            }
        )
        keep_cross = [
            "alignment_column",
            "residue_index",
            "cross_lddt",
            "cross_lddt_coverage",
            "cross_presence",
            "cross_n_lddt_comparisons",
        ]
        table = within[
            [
                "alignment_column",
                "residue_index",
                "chain_id",
                "resseq",
                "icode",
                "residue_label",
                "residue_name",
                "aa",
                "local_pLDDT",
                "within_lddt",
                "within_lddt_coverage",
                "within_presence",
                "within_n_lddt_comparisons",
            ]
        ].merge(cross[keep_cross], on=["alignment_column", "residue_index"], how="left")
        table.insert(0, "medoid_id", medoid_id)
        table.insert(0, "cluster_id", cluster_id)
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def _safe_min(values: pd.Series) -> float:
    finite = pd.to_numeric(values, errors="coerce").dropna()
    return float(finite.min()) if len(finite) else np.nan


def _safe_median(values: pd.Series) -> float:
    finite = pd.to_numeric(values, errors="coerce").dropna()
    return float(finite.median()) if len(finite) else np.nan


def call_pan_core_columns(
    aa: dict[str, str],
    di: dict[str, str],
    members: pd.DataFrame,
    medoids: pd.DataFrame,
    maps: dict[str, dict[int, object]],
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cluster_ids = medoids["cluster_id"].astype(str).tolist()
    members_by_cluster = {
        cluster_id: members.loc[
            members["cluster_id"].astype(str) == cluster_id, "protein_id"
        ].astype(str).tolist()
        for cluster_id in cluster_ids
    }
    anchors = _anchor_metrics(medoids, members_by_cluster, maps, config)
    pc = config["pan_core"]
    n_clusters = len(cluster_ids)
    supported_min = int(
        math.ceil(float(pc["supported_min_cluster_pass_fraction"]) * n_clusters)
    )
    length = len(next(iter(aa.values())))
    rows = []
    for column in range(length):
        focal = anchors[anchors["alignment_column"] == column + 1].copy()
        medoids_present = int(len(focal))
        if medoids_present:
            within_pass = (
                (focal["local_pLDDT"] >= float(pc["min_local_plddt"]))
                & (focal["within_lddt"] >= float(pc["within_high"]))
                & (
                    focal["within_lddt_coverage"]
                    >= float(pc["min_within_lddt_coverage"])
                )
            )
            cross_pass = (
                (focal["local_pLDDT"] >= float(pc["min_local_plddt"]))
                & (focal["cross_lddt"] >= float(pc["cross_high"]))
                & (
                    focal["cross_lddt_coverage"]
                    >= float(pc["min_cross_lddt_coverage"])
                )
                & (focal["cross_presence"] >= float(pc["min_cross_presence"]))
            )
            n_within_pass = int(within_pass.sum())
            n_cross_pass = int(cross_pass.sum())
        else:
            n_within_pass = 0
            n_cross_pass = 0
        column_metrics = _cluster_column_metrics(
            column, cluster_ids, members_by_cluster, aa, di
        )
        all_medoids_present = medoids_present == n_clusters
        occupancy_pass = (
            column_metrics["min_subcluster_occupancy"]
            >= float(pc["min_subcluster_occupancy"])
        )
        confidence_pass = (
            medoids_present == n_clusters
            and _safe_min(focal["local_pLDDT"]) >= float(pc["min_local_plddt"])
        )
        strict = bool(
            all_medoids_present
            and occupancy_pass
            and confidence_pass
            and n_within_pass == n_clusters
            and n_cross_pass == n_clusters
        )
        supported = bool(
            all_medoids_present
            and occupancy_pass
            and confidence_pass
            and n_within_pass >= supported_min
            and n_cross_pass >= supported_min
            and _safe_median(focal["within_lddt"]) >= float(pc["within_high"])
            and _safe_median(focal["cross_lddt"]) >= float(pc["cross_high"])
        )
        pass_fraction = (
            min(n_within_pass, n_cross_pass) / n_clusters if n_clusters else np.nan
        )
        plddt_weight = np.clip(
            (_safe_median(focal["local_pLDDT"]) - 50.0) / 40.0, 0.0, 1.0
        )
        support_score = (
            _safe_median(focal["within_lddt"])
            * _safe_median(focal["cross_lddt"])
            * column_metrics["min_subcluster_occupancy"]
            * pass_fraction
            * plddt_weight
            if medoids_present
            else np.nan
        )
        rows.append(
            {
                "reference_alignment_column": column + 1,
                "n_subclusters_total": n_clusters,
                "n_medoids_present": medoids_present,
                "n_subclusters_passing_within": n_within_pass,
                "n_subclusters_passing_cross": n_cross_pass,
                "min_medoid_pLDDT": _safe_min(focal["local_pLDDT"]),
                "median_medoid_pLDDT": _safe_median(focal["local_pLDDT"]),
                "min_within_lddt": _safe_min(focal["within_lddt"]),
                "median_within_lddt": _safe_median(focal["within_lddt"]),
                "min_cross_lddt": _safe_min(focal["cross_lddt"]),
                "median_cross_lddt": _safe_median(focal["cross_lddt"]),
                "min_within_lddt_coverage": _safe_min(
                    focal["within_lddt_coverage"]
                ),
                "min_cross_lddt_coverage": _safe_min(
                    focal["cross_lddt_coverage"]
                ),
                "min_cross_presence": _safe_min(focal["cross_presence"]),
                **column_metrics,
                "strict_pan_core": strict,
                "supported_pan_core": supported,
                "pan_core_class": (
                    "strict pan-shared core"
                    if strict
                    else "supported pan-shared core"
                    if supported
                    else "not pan-core"
                ),
                "pan_core_support_score": support_score,
            }
        )
    return pd.DataFrame(rows), anchors


def _anchor_mapping(
    columns: pd.DataFrame,
    anchor_id: str,
    anchor_map: dict[int, object],
) -> pd.DataFrame:
    rows = []
    for row in columns.itertuples(index=False):
        column0 = int(row.reference_alignment_column) - 1
        residue = anchor_map.get(column0)
        rows.append(
            {
                "reference_alignment_column": int(row.reference_alignment_column),
                "anchor_medoid_id": anchor_id,
                "anchor_residue_index": residue.seq_index if residue else np.nan,
                "anchor_residue_label": residue.residue_label if residue else "",
                "anchor_residue_name": residue.resname if residue else "",
                "anchor_local_pLDDT": residue.plddt if residue else np.nan,
            }
        )
    return pd.DataFrame(rows)


def summarize_core_regions(
    core_columns: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    selected = core_columns[core_columns["strict_pan_core"]].copy()
    selected = selected[selected["anchor_residue_index"].notna()].sort_values(
        "anchor_residue_index"
    )
    indices = selected["anchor_residue_index"].astype(int).tolist()
    groups = contiguous_groups(
        indices, int(config["pan_core"]["region_max_bridge_residues"])
    )
    rows = []
    by_index = selected.set_index("anchor_residue_index")
    region_id = 0
    for group in groups:
        window = selected[
            (selected["anchor_residue_index"] >= min(group))
            & (selected["anchor_residue_index"] <= max(group))
        ]
        if len(window) < int(config["pan_core"]["region_min_length"]):
            continue
        region_id += 1
        first = by_index.loc[min(group)]
        last = by_index.loc[max(group)]
        rows.append(
            {
                "region_id": f"PCR{region_id:02d}",
                "anchor_medoid_id": first["anchor_medoid_id"],
                "anchor_residue_start": min(group),
                "anchor_residue_end": max(group),
                "anchor_residue_range": (
                    f"{first['anchor_residue_label']}-{last['anchor_residue_label']}"
                ),
                "n_strict_core_residues": len(group),
                "sequence_span_residues": max(group) - min(group) + 1,
                "mean_min_subcluster_occupancy": window[
                    "min_subcluster_occupancy"
                ].mean(),
                "mean_median_within_lddt": window["median_within_lddt"].mean(),
                "mean_median_cross_lddt": window["median_cross_lddt"].mean(),
                "mean_pan_3di_conservation": window["pan_3di_conservation"].mean(),
                "mean_pan_aa_conservation": window["pan_aa_conservation"].mean(),
                "mean_pan_core_support_score": window[
                    "pan_core_support_score"
                ].mean(),
            }
        )
    return pd.DataFrame(rows)


def spatial_core_components(
    core_columns: pd.DataFrame,
    anchor_map: dict[int, object],
    config: dict,
) -> tuple[pd.DataFrame, dict[int, int]]:
    strict = core_columns[core_columns["strict_pan_core"]].copy()
    residues = []
    for row in strict.itertuples(index=False):
        residue = anchor_map.get(int(row.reference_alignment_column) - 1)
        if residue is not None:
            residues.append(residue)
    if not residues:
        return pd.DataFrame(), {}
    coords = np.asarray([x.ca for x in residues], dtype=float)
    delta = coords[:, None, :] - coords[None, :, :]
    distances = np.sqrt(np.sum(delta * delta, axis=-1))
    cutoff = float(config["pan_core"]["contact_cutoff_A"])
    unseen = set(range(len(residues)))
    components = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        stack = [start]
        component = [start]
        while stack:
            i = stack.pop()
            neighbors = [
                j
                for j in list(unseen)
                if distances[i, j] <= cutoff
                or abs(residues[i].seq_index - residues[j].seq_index) == 1
            ]
            for j in neighbors:
                unseen.remove(j)
                stack.append(j)
                component.append(j)
        components.append(sorted(component))
    components.sort(key=lambda x: (-len(x), min(residues[i].seq_index for i in x)))
    rows = []
    residue_to_component: dict[int, int] = {}
    total = len(residues)
    for component_id, component in enumerate(components, 1):
        indices = sorted(residues[i].seq_index for i in component)
        for index in indices:
            residue_to_component[index] = component_id
        rows.append(
            {
                "component_id": f"PCC{component_id:02d}",
                "n_core_residues": len(indices),
                "fraction_of_strict_core": len(indices) / total,
                "anchor_residue_index_min": min(indices),
                "anchor_residue_index_max": max(indices),
                "anchor_residue_indices": ",".join(str(x) for x in indices),
                "contact_cutoff_A": cutoff,
            }
        )
    return pd.DataFrame(rows), residue_to_component


def threshold_sensitivity(
    core_columns: pd.DataFrame,
    anchor_metrics: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    pc = config["pan_core"]
    n_clusters = int(core_columns["n_subclusters_total"].max())
    rows = []
    for occupancy in (0.50, 0.70, 0.80):
        for within in (0.70, 0.75, 0.80):
            for cross in (0.60, 0.65, 0.70):
                for required_fraction in (0.75, 1.0):
                    required = int(math.ceil(required_fraction * n_clusters))
                    anchor = anchor_metrics.copy()
                    anchor["passes_candidate_within"] = (
                        (
                            anchor["local_pLDDT"]
                            >= float(pc["min_local_plddt"])
                        )
                        & (anchor["within_lddt"] >= within)
                        & (
                            anchor["within_lddt_coverage"]
                            >= float(pc["min_within_lddt_coverage"])
                        )
                    )
                    anchor["passes_candidate_cross"] = (
                        (
                            anchor["local_pLDDT"]
                            >= float(pc["min_local_plddt"])
                        )
                        & (anchor["cross_lddt"] >= cross)
                        & (
                            anchor["cross_lddt_coverage"]
                            >= float(pc["min_cross_lddt_coverage"])
                        )
                        & (
                            anchor["cross_presence"]
                            >= float(pc["min_cross_presence"])
                        )
                    )
                    passing = (
                        anchor.groupby("alignment_column", as_index=False)
                        .agg(
                            candidate_n_within_pass=(
                                "passes_candidate_within",
                                "sum",
                            ),
                            candidate_n_cross_pass=(
                                "passes_candidate_cross",
                                "sum",
                            ),
                        )
                        .rename(
                            columns={
                                "alignment_column": "reference_alignment_column"
                            }
                        )
                    )
                    evaluated = core_columns.merge(
                        passing, on="reference_alignment_column", how="left"
                    ).fillna(
                        {
                            "candidate_n_within_pass": 0,
                            "candidate_n_cross_pass": 0,
                        }
                    )
                    mask = (
                        (evaluated["n_medoids_present"] == n_clusters)
                        & (
                            evaluated["min_subcluster_occupancy"]
                            >= occupancy
                        )
                        & (
                            evaluated["min_medoid_pLDDT"]
                            >= float(pc["min_local_plddt"])
                        )
                        & (
                            evaluated["min_within_lddt_coverage"]
                            >= float(pc["min_within_lddt_coverage"])
                        )
                        & (
                            evaluated["min_cross_lddt_coverage"]
                            >= float(pc["min_cross_lddt_coverage"])
                        )
                        & (
                            evaluated["min_cross_presence"]
                            >= float(pc["min_cross_presence"])
                        )
                        & (
                            evaluated["candidate_n_within_pass"]
                            >= required
                        )
                        & (
                            evaluated["candidate_n_cross_pass"]
                            >= required
                        )
                    )
                    rows.append(
                        {
                            "min_subcluster_occupancy": occupancy,
                            "within_lddt_threshold": within,
                            "cross_lddt_threshold": cross,
                            "required_subcluster_fraction": required_fraction,
                            "n_pan_core_columns": int(mask.sum()),
                        }
                    )
    return pd.DataFrame(rows)


def _consensus_excluding(
    records: dict[str, str],
    ids: list[str],
    subject_id: str,
    column: int,
) -> str:
    states = [
        records[x][column]
        for x in ids
        if x != subject_id and records[x][column] not in "-."
    ]
    if not states:
        return ""
    counts = Counter(states)
    return sorted(counts.items(), key=lambda x: (-x[1], x[0]))[0][0]


def calibrate_reference_members(
    core_columns: pd.DataFrame,
    aa: dict[str, str],
    di: dict[str, str],
    members: pd.DataFrame,
    medoids: pd.DataFrame,
    maps: dict[str, dict[int, object]],
    config: dict,
) -> pd.DataFrame:
    core = core_columns[core_columns["strict_pan_core"]].copy()
    core_columns0 = core["reference_alignment_column"].astype(int).sub(1).tolist()
    reference_ids = members["protein_id"].astype(str).tolist()
    medoid_ids = medoids["medoid_id"].astype(str).tolist()
    cluster_by_id = members.set_index("protein_id")["cluster_id"].astype(str).to_dict()
    rows = []
    for protein_id in reference_ids:
        comparisons = [x for x in medoid_ids if x != protein_id]
        local = aggregate_local_lddt(protein_id, comparisons, maps, config)
        local_by_column = local.set_index("alignment_column")["lddt"].to_dict()
        present = [c for c in core_columns0 if aa[protein_id][c] not in "-."]
        di_evaluable = [c for c in present if di[protein_id][c] not in "-."]
        di_matches = [
            di[protein_id][c]
            == _consensus_excluding(di, reference_ids, protein_id, c)
            for c in di_evaluable
        ]
        aa_matches = [
            aa[protein_id][c]
            == _consensus_excluding(aa, reference_ids, protein_id, c)
            for c in present
        ]
        lddt = [
            local_by_column.get(c + 1, np.nan)
            for c in present
            if np.isfinite(local_by_column.get(c + 1, np.nan))
        ]
        rows.append(
            {
                "protein_id": protein_id,
                "cluster_id": cluster_by_id[protein_id],
                "is_query": False,
                "n_strict_pan_core_columns": len(core_columns0),
                "n_core_columns_present": len(present),
                "pan_core_coverage": (
                    len(present) / len(core_columns0) if core_columns0 else np.nan
                ),
                "n_3di_evaluable": len(di_evaluable),
                "pan_core_3di_consensus_match_fraction": (
                    float(np.mean(di_matches)) if di_matches else np.nan
                ),
                "n_aa_evaluable": len(present),
                "pan_core_aa_consensus_match_fraction": (
                    float(np.mean(aa_matches)) if aa_matches else np.nan
                ),
                "n_local_lddt_evaluable": len(lddt),
                "median_pan_core_local_lddt": (
                    float(np.median(lddt)) if lddt else np.nan
                ),
                "mean_pan_core_local_lddt": (
                    float(np.mean(lddt)) if lddt else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def map_query_to_pan_core(
    core_columns: pd.DataFrame,
    anchor_id: str,
    reference_maps: dict[str, dict[int, object]],
    query_aa_path: str | Path,
    query_di_path: str | Path,
    query_id: str,
    query_structure: str | Path,
    query_chain: str,
    members: pd.DataFrame,
    medoids: pd.DataFrame,
    config: dict,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reference_ids = members["protein_id"].astype(str).tolist()
    query_manifest = members.copy()
    query_row = {column: "" for column in query_manifest.columns}
    query_row.update(
        {
            "protein_id": query_id,
            "cluster_id": "QUERY",
            "_structure_path": str(Path(query_structure).resolve()),
            "chain_id": query_chain,
        }
    )
    query_manifest = pd.concat(
        [query_manifest, pd.DataFrame([query_row])], ignore_index=True
    )
    required = reference_ids + [query_id]
    query_aa = _resolve_alignment_to_protein_ids(
        query_aa_path, query_manifest, required
    )
    query_di = _resolve_alignment_to_protein_ids(
        query_di_path, query_manifest, required
    )
    query_maps, query_residues, mapping_qc = _build_maps(
        query_aa, query_manifest, config
    )

    ref_anchor_by_residue = {
        residue.seq_index: column
        for column, residue in reference_maps[anchor_id].items()
    }
    query_anchor_by_residue = {
        residue.seq_index: column for column, residue in query_maps[anchor_id].items()
    }
    medoid_ids = medoids["medoid_id"].astype(str).tolist()
    query_lddt = aggregate_local_lddt(
        query_id, medoid_ids, query_maps, config
    ).set_index("alignment_column")

    rows = []
    for row in core_columns.itertuples(index=False):
        reference_column0 = int(row.reference_alignment_column) - 1
        anchor_residue = reference_maps[anchor_id].get(reference_column0)
        query_column0 = (
            query_anchor_by_residue.get(anchor_residue.seq_index)
            if anchor_residue is not None
            else None
        )
        query_residue = (
            query_maps[query_id].get(query_column0)
            if query_column0 is not None
            else None
        )
        local_lddt = (
            query_lddt.loc[query_column0 + 1, "lddt"]
            if query_column0 is not None
            and query_column0 + 1 in query_lddt.index
            else np.nan
        )
        query_3di = (
            query_di[query_id][query_column0]
            if query_column0 is not None
            and query_di[query_id][query_column0] not in "-."
            else ""
        )
        rows.append(
            {
                "query_id": query_id,
                "reference_alignment_column": int(row.reference_alignment_column),
                "pan_core_class": row.pan_core_class,
                "strict_pan_core": bool(row.strict_pan_core),
                "supported_pan_core": bool(row.supported_pan_core),
                "anchor_medoid_id": anchor_id,
                "anchor_residue_index": (
                    anchor_residue.seq_index if anchor_residue else np.nan
                ),
                "query_alignment_column": (
                    query_column0 + 1 if query_column0 is not None else np.nan
                ),
                "query_residue_index": (
                    query_residue.seq_index if query_residue else np.nan
                ),
                "query_residue_label": (
                    query_residue.residue_label if query_residue else ""
                ),
                "query_aa": (
                    query_aa[query_id][query_column0]
                    if query_column0 is not None
                    and query_aa[query_id][query_column0] not in "-."
                    else ""
                ),
                "query_3di": query_3di,
                "pan_aa_consensus": row.pan_aa_consensus,
                "pan_3di_consensus": row.pan_3di_consensus,
                "query_matches_pan_aa_consensus": bool(
                    query_column0 is not None
                    and query_aa[query_id][query_column0] not in "-."
                    and query_aa[query_id][query_column0]
                    == row.pan_aa_consensus
                ),
                "query_matches_pan_3di_consensus": bool(
                    query_3di and query_3di == row.pan_3di_consensus
                ),
                "query_vs_four_medoids_local_lddt": local_lddt,
                "pan_core_support_score": row.pan_core_support_score,
            }
        )
    mapping = pd.DataFrame(rows)
    summaries = []
    for tier, mask in [
        ("strict pan-shared core", mapping["strict_pan_core"]),
        ("supported pan-shared core", mapping["supported_pan_core"]),
    ]:
        selected = mapping[mask].copy()
        aligned = selected[selected["query_residue_index"].notna()]
        lddt = pd.to_numeric(
            aligned["query_vs_four_medoids_local_lddt"], errors="coerce"
        ).dropna()
        summaries.append(
            {
                "query_id": query_id,
                "core_tier": tier,
                "n_reference_core_columns": len(selected),
                "n_query_aligned": len(aligned),
                "query_core_coverage": (
                    len(aligned) / len(selected) if len(selected) else np.nan
                ),
                "n_query_3di_evaluable": int((aligned["query_3di"] != "").sum()),
                "n_query_aa_evaluable": int((aligned["query_aa"] != "").sum()),
                "query_pan_aa_consensus_match_fraction": (
                    aligned.loc[
                        aligned["query_aa"] != "",
                        "query_matches_pan_aa_consensus",
                    ].mean()
                    if (aligned["query_aa"] != "").any()
                    else np.nan
                ),
                "query_pan_3di_consensus_match_fraction": (
                    aligned.loc[
                        aligned["query_3di"] != "",
                        "query_matches_pan_3di_consensus",
                    ].mean()
                    if (aligned["query_3di"] != "").any()
                    else np.nan
                ),
                "n_query_local_lddt_evaluable": len(lddt),
                "median_query_vs_medoids_local_lddt": (
                    float(lddt.median()) if len(lddt) else np.nan
                ),
                "mean_query_vs_medoids_local_lddt": (
                    float(lddt.mean()) if len(lddt) else np.nan
                ),
                "fraction_query_local_lddt_ge_cross_threshold": (
                    float((lddt >= float(config["pan_core"]["cross_high"])).mean())
                    if len(lddt)
                    else np.nan
                ),
                "fraction_query_local_lddt_ge_within_threshold": (
                    float((lddt >= float(config["pan_core"]["within_high"])).mean())
                    if len(lddt)
                    else np.nan
                ),
            }
        )
    summary = pd.DataFrame(summaries)

    strict = mapping[mapping["strict_pan_core"]].copy()
    _, query_structure_residues = load_structure_residues(
        query_structure, query_chain or None
    )
    coverage_values = {}
    lddt_values = {}
    di_values = {}
    for row in strict.itertuples(index=False):
        if pd.isna(row.query_residue_index):
            continue
        residue = query_structure_residues[int(row.query_residue_index) - 1]
        key = (residue.chain_id, residue.resseq, residue.icode)
        coverage_values[key] = 1.0
        if np.isfinite(row.query_vs_four_medoids_local_lddt):
            lddt_values[key] = float(row.query_vs_four_medoids_local_lddt)
        di_values[key] = 1.0 if row.query_matches_pan_3di_consensus else 0.25

    query_component_rows = []
    aligned_strict = strict[strict["query_residue_index"].notna()].copy()
    if len(aligned_strict):
        query_indices = sorted(
            aligned_strict["query_residue_index"].astype(int).unique().tolist()
        )
        query_core_residues = [
            query_structure_residues[index - 1] for index in query_indices
        ]
        coords = np.asarray([x.ca for x in query_core_residues], dtype=float)
        delta = coords[:, None, :] - coords[None, :, :]
        distances = np.sqrt(np.sum(delta * delta, axis=-1))
        cutoff = float(config["pan_core"]["contact_cutoff_A"])
        unseen = set(range(len(query_core_residues)))
        components = []
        while unseen:
            start = min(unseen)
            unseen.remove(start)
            stack = [start]
            component = [start]
            while stack:
                i = stack.pop()
                neighbors = [
                    j
                    for j in list(unseen)
                    if distances[i, j] <= cutoff
                    or abs(
                        query_core_residues[i].seq_index
                        - query_core_residues[j].seq_index
                    )
                    == 1
                ]
                for j in neighbors:
                    unseen.remove(j)
                    stack.append(j)
                    component.append(j)
            components.append(sorted(component))
        components.sort(
            key=lambda x: (
                -len(x),
                min(query_core_residues[i].seq_index for i in x),
            )
        )
        for component_id, component in enumerate(components, 1):
            indices = sorted(
                query_core_residues[i].seq_index for i in component
            )
            query_component_rows.append(
                {
                    "component_id": f"QPCC{component_id:02d}",
                    "n_query_core_residues": len(indices),
                    "fraction_of_aligned_strict_core": (
                        len(indices) / len(query_core_residues)
                    ),
                    "query_residue_index_min": min(indices),
                    "query_residue_index_max": max(indices),
                    "query_residue_indices": ",".join(str(x) for x in indices),
                    "contact_cutoff_A": cutoff,
                }
            )
    query_components = pd.DataFrame(query_component_rows)
    write_tsv(
        query_components,
        output_dir / f"{query_id}_strict_pan_core_spatial_components.tsv",
    )
    if len(query_components):
        strict_index = summary["core_tier"] == "strict pan-shared core"
        summary.loc[strict_index, "n_query_spatial_components"] = len(
            query_components
        )
        summary.loc[
            strict_index, "largest_query_component_fraction"
        ] = query_components.iloc[0]["fraction_of_aligned_strict_core"]
    output_dir.mkdir(parents=True, exist_ok=True)
    write_metric_pdb(
        query_structure,
        output_dir / f"{query_id}_strict_pan_core_coverage.pdb",
        coverage_values,
    )
    write_metric_pdb(
        query_structure,
        output_dir / f"{query_id}_strict_pan_core_local_lddt.pdb",
        lddt_values,
    )
    write_metric_pdb(
        query_structure,
        output_dir / f"{query_id}_strict_pan_core_3di_match.pdb",
        di_values,
    )
    pml = [
        "bg_color white",
        "set ray_opaque_background, off",
        f"load {query_id}_strict_pan_core_coverage.pdb, pan_core_coverage",
        f"load {query_id}_strict_pan_core_local_lddt.pdb, pan_core_lddt",
        f"load {query_id}_strict_pan_core_3di_match.pdb, pan_core_3di",
        "hide everything",
        "show cartoon, pan_core_coverage",
        "color grey80, pan_core_coverage",
        "select strict_core, pan_core_coverage and b > 99",
        "color marine, strict_core",
        "show sticks, strict_core",
        "spectrum b, red_white_blue, pan_core_lddt, minimum=0, maximum=100",
        "spectrum b, grey70_white_blue, pan_core_3di, minimum=0, maximum=100",
        "disable pan_core_lddt",
        "disable pan_core_3di",
        "orient pan_core_coverage",
    ]
    (output_dir / f"view_{query_id}_pan_core.pml").write_text("\n".join(pml) + "\n")
    return mapping, summary, mapping_qc


def add_query_to_calibration(
    calibration: pd.DataFrame,
    query_summary: pd.DataFrame,
    query_id: str,
) -> pd.DataFrame:
    strict = query_summary[
        query_summary["core_tier"] == "strict pan-shared core"
    ].iloc[0]
    row = {
        "protein_id": query_id,
        "cluster_id": "QUERY",
        "is_query": True,
        "n_strict_pan_core_columns": strict["n_reference_core_columns"],
        "n_core_columns_present": strict["n_query_aligned"],
        "pan_core_coverage": strict["query_core_coverage"],
        "n_3di_evaluable": strict["n_query_3di_evaluable"],
        "pan_core_3di_consensus_match_fraction": strict[
            "query_pan_3di_consensus_match_fraction"
        ],
        "n_aa_evaluable": strict["n_query_aa_evaluable"],
        "pan_core_aa_consensus_match_fraction": strict[
            "query_pan_aa_consensus_match_fraction"
        ],
        "n_local_lddt_evaluable": strict["n_query_local_lddt_evaluable"],
        "median_pan_core_local_lddt": strict[
            "median_query_vs_medoids_local_lddt"
        ],
        "mean_pan_core_local_lddt": strict[
            "mean_query_vs_medoids_local_lddt"
        ],
        "fraction_local_lddt_ge_cross_threshold": strict[
            "fraction_query_local_lddt_ge_cross_threshold"
        ],
        "fraction_local_lddt_ge_within_threshold": strict[
            "fraction_query_local_lddt_ge_within_threshold"
        ],
    }
    return pd.concat([calibration, pd.DataFrame([row])], ignore_index=True)


def plot_pan_core_summary(
    core_columns: pd.DataFrame,
    calibration: pd.DataFrame,
    sensitivity: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="ticks", context="paper")
    colors = {
        "SC01": "#0072B2",
        "SC02": "#56B4E9",
        "SC03": "#E69F00",
        "SC04": "#009E73",
        "QUERY": "#CC3377",
    }
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.35))

    class_counts = (
        core_columns["pan_core_class"]
        .value_counts()
        .reindex(
            [
                "strict pan-shared core",
                "supported pan-shared core",
                "not pan-core",
            ],
            fill_value=0,
        )
    )
    axes[0].bar(
        ["Strict", "Supported", "Other"],
        class_counts.values,
        color=["#225EA8", "#9ECAE1", "#D9D9D9"],
        edgecolor="none",
    )
    axes[0].set_ylabel("Common-MSA columns")
    axes[0].set_title("Pan-core calls")

    metrics = [
        ("pan_core_coverage", "Core coverage"),
        (
            "pan_core_3di_consensus_match_fraction",
            "3Di consensus match",
        ),
        ("median_pan_core_local_lddt", "Median local lDDT"),
    ]
    long = calibration.melt(
        id_vars=["protein_id", "cluster_id", "is_query"],
        value_vars=[x[0] for x in metrics],
        var_name="metric",
        value_name="value",
    )
    label_map = dict(metrics)
    long["metric_label"] = long["metric"].map(label_map)
    for metric_label, group in long.groupby("metric_label", sort=False):
        x = list(label_map.values()).index(metric_label)
        refs = group[~group["is_query"]]
        jitter = np.linspace(-0.14, 0.14, len(refs))
        for offset, row in zip(jitter, refs.itertuples(index=False)):
            axes[1].scatter(
                x + offset,
                row.value,
                s=23,
                color=colors.get(row.cluster_id, "#777777"),
                alpha=0.85,
                linewidth=0,
            )
        query = group[group["is_query"]]
        if len(query):
            axes[1].scatter(
                x,
                query.iloc[0]["value"],
                s=75,
                marker="D",
                color=colors["QUERY"],
                edgecolor="black",
                linewidth=0.6,
                zorder=5,
            )
    axes[1].set_xticks(range(len(metrics)))
    axes[1].set_xticklabels([x[1] for x in metrics], rotation=20, ha="right")
    axes[1].set_ylim(0, 1.03)
    axes[1].set_ylabel("Fraction or score")
    axes[1].set_title("Reference calibration and query")

    subset = sensitivity[
        (sensitivity["min_subcluster_occupancy"] == 0.70)
        & (sensitivity["required_subcluster_fraction"] == 1.0)
    ].copy()
    pivot = subset.pivot(
        index="within_lddt_threshold",
        columns="cross_lddt_threshold",
        values="n_pan_core_columns",
    )
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".0f",
        cmap="Blues",
        cbar_kws={"label": "Core columns"},
        ax=axes[2],
    )
    axes[2].set_xlabel("Cross-subcluster lDDT threshold")
    axes[2].set_ylabel("Within-subcluster lDDT threshold")
    axes[2].set_title("Strict-call sensitivity")
    sns.despine(ax=axes[0])
    sns.despine(ax=axes[1])
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(
            output_dir / f"pan_shared_core_overview.{suffix}",
            dpi=dpi if suffix == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)


def compute_pan_core_workflow(
    members_path: str | Path,
    medoids_path: str | Path,
    reference_aa_path: str | Path,
    reference_di_path: str | Path,
    output_dir: str | Path,
    config: dict,
    anchor_medoid_id: str | None = None,
    query_id: str | None = None,
    query_aa_path: str | Path | None = None,
    query_di_path: str | Path | None = None,
    query_structure: str | Path | None = None,
    query_chain: str = "",
) -> dict:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    members, medoids = _prepare_members(members_path, medoids_path)
    reference_ids = members["protein_id"].astype(str).tolist()
    aa = _resolve_alignment_to_protein_ids(
        reference_aa_path, members, reference_ids
    )
    di = _resolve_alignment_to_protein_ids(
        reference_di_path, members, reference_ids
    )
    if set(aa) != set(di):
        raise ValueError("Reference AA and 3Di alignments resolve to different proteins")
    for protein_id in reference_ids:
        if len(aa[protein_id]) != len(di[protein_id]):
            raise ValueError(f"AA and 3Di lengths differ for {protein_id}")
        if not np.array_equal(
            np.asarray(list(aa[protein_id])) == "-",
            np.asarray(list(di[protein_id])) == "-",
        ):
            raise ValueError(f"AA and 3Di gaps differ for {protein_id}")

    maps, _, mapping_qc = _build_maps(aa, members, config)
    core_columns, anchor_metrics = call_pan_core_columns(
        aa, di, members, medoids, maps, config
    )
    anchor_id = anchor_medoid_id or str(medoids.iloc[0]["medoid_id"])
    if anchor_id not in set(medoids["medoid_id"]):
        raise ValueError(f"Anchor {anchor_id} is not one of the supplied medoids")
    anchor_table = _anchor_mapping(core_columns, anchor_id, maps[anchor_id])
    core_columns = core_columns.merge(
        anchor_table, on="reference_alignment_column", how="left"
    )
    all_component_tables = []
    component_maps: dict[str, dict[int, int]] = {}
    for medoid_id in medoids["medoid_id"].astype(str):
        component_table, component_map = spatial_core_components(
            core_columns, maps[medoid_id], config
        )
        if not component_table.empty:
            component_table.insert(0, "medoid_id", medoid_id)
            all_component_tables.append(component_table)
        component_maps[medoid_id] = component_map
    components = (
        pd.concat(all_component_tables, ignore_index=True)
        if all_component_tables
        else pd.DataFrame()
    )
    residue_to_component = component_maps[anchor_id]
    core_columns["spatial_component_id"] = core_columns[
        "anchor_residue_index"
    ].map(
        lambda x: (
            f"PCC{residue_to_component.get(int(x), 0):02d}"
            if pd.notna(x) and int(x) in residue_to_component
            else ""
        )
    )
    regions = summarize_core_regions(core_columns, config)
    sensitivity = threshold_sensitivity(core_columns, anchor_metrics, config)
    calibration = calibrate_reference_members(
        core_columns, aa, di, members, medoids, maps, config
    )

    query_mapping = pd.DataFrame()
    query_summary = pd.DataFrame()
    query_qc = pd.DataFrame()
    if query_id:
        if not all(
            x is not None
            for x in (query_aa_path, query_di_path, query_structure)
        ):
            raise ValueError(
                "query-id requires query-aa, query-3di, and query-structure"
            )
        query_output = output_dir / "query_mapping"
        query_mapping, query_summary, query_qc = map_query_to_pan_core(
            core_columns,
            anchor_id,
            maps,
            query_aa_path,
            query_di_path,
            query_id,
            query_structure,
            query_chain,
            members,
            medoids,
            config,
            query_output,
        )
        calibration = add_query_to_calibration(
            calibration, query_summary, query_id
        )

    write_tsv(mapping_qc, output_dir / "reference_alignment_mapping_qc.tsv")
    write_tsv(core_columns, output_dir / "pan_shared_core_columns.tsv")
    write_tsv(anchor_metrics, output_dir / "pan_core_medoid_anchor_metrics.tsv")
    write_tsv(regions, output_dir / "pan_shared_core_regions.tsv")
    write_tsv(components, output_dir / "pan_shared_core_spatial_components.tsv")
    write_tsv(sensitivity, output_dir / "pan_core_threshold_sensitivity.tsv")
    write_tsv(calibration, output_dir / "pan_core_reference_calibration.tsv")
    if not query_mapping.empty:
        write_tsv(
            query_mapping,
            output_dir / "query_mapping" / f"{query_id}_pan_core_mapping.tsv",
        )
        write_tsv(
            query_summary,
            output_dir / "query_mapping" / f"{query_id}_pan_core_summary.tsv",
        )
        write_tsv(
            query_qc,
            output_dir / "query_mapping" / f"{query_id}_alignment_mapping_qc.tsv",
        )

    members_index = members.set_index("protein_id")
    for medoid_id in medoids["medoid_id"].astype(str):
        medoid_values = {}
        for row in core_columns.itertuples(index=False):
            residue = maps[medoid_id].get(
                int(row.reference_alignment_column) - 1
            )
            if residue is None:
                continue
            value = (
                1.0
                if row.strict_pan_core
                else 0.60
                if row.supported_pan_core
                else 0.0
            )
            medoid_values[
                (residue.chain_id, residue.resseq, residue.icode)
            ] = value
        write_metric_pdb(
            members_index.loc[medoid_id, "_structure_path"],
            output_dir / f"{medoid_id}_pan_core_classes.pdb",
            medoid_values,
        )

    plot_pan_core_summary(
        core_columns,
        calibration,
        sensitivity,
        output_dir / "figures",
        int(config["visualization"]["dpi"]),
    )
    summary = {
        "n_reference_proteins": len(reference_ids),
        "n_subclusters": len(medoids),
        "reference_alignment_length": len(core_columns),
        "n_strict_pan_core_columns": int(core_columns["strict_pan_core"].sum()),
        "n_supported_pan_core_columns": int(
            core_columns["supported_pan_core"].sum()
        ),
        "n_supported_only_pan_core_columns": int(
            (
                core_columns["supported_pan_core"]
                & ~core_columns["strict_pan_core"]
            ).sum()
        ),
        "n_strict_core_regions": len(regions),
        "spatial_components_per_medoid": (
            components.groupby("medoid_id").size().astype(int).to_dict()
            if len(components)
            else {}
        ),
        "largest_spatial_component_fraction_per_medoid": (
            components.groupby("medoid_id")[
                "fraction_of_strict_core"
            ].max().astype(float).to_dict()
            if len(components)
            else {}
        ),
        "anchor_medoid_id": anchor_id,
        "query_id": query_id,
    }
    if not query_summary.empty:
        strict = query_summary[
            query_summary["core_tier"] == "strict pan-shared core"
        ].iloc[0]
        summary.update(
            {
                "query_strict_core_coverage": strict["query_core_coverage"],
                "query_strict_core_3di_match": strict[
                    "query_pan_3di_consensus_match_fraction"
                ],
                "query_strict_core_aa_match": strict[
                    "query_pan_aa_consensus_match_fraction"
                ],
                "query_strict_core_median_local_lddt": strict[
                    "median_query_vs_medoids_local_lddt"
                ],
                "query_fraction_local_lddt_ge_cross_threshold": strict[
                    "fraction_query_local_lddt_ge_cross_threshold"
                ],
                "query_fraction_local_lddt_ge_within_threshold": strict[
                    "fraction_query_local_lddt_ge_within_threshold"
                ],
            }
        )
    (output_dir / "pan_core_run_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n"
    )
    return summary
