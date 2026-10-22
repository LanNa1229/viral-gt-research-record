from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .common import (
    load_structure_residues,
    map_alignment_sequence_to_structure,
    normalized_entropy_conservation,
    read_fasta,
    read_tsv,
    require_columns,
    resolve_alignment_records,
    sanitize_group,
    write_metric_pdb,
)
from .medoids import lookup_pair, pair_lookup, standardize_pairwise, translate_pairwise_ids


def load_alignment_pair(prefix: str | Path, expected_ids: Iterable[str]) -> tuple[dict[str, str], dict[str, str]]:
    prefix = Path(prefix)
    aa_path = Path(f"{prefix}_aa.fa")
    di_path = Path(f"{prefix}_3di.fa")
    if not aa_path.exists() or not di_path.exists():
        raise FileNotFoundError(f"Expected FoldMason outputs not found: {aa_path} and {di_path}")
    aa = resolve_alignment_records(read_fasta(aa_path), expected_ids)
    di = resolve_alignment_records(read_fasta(di_path), expected_ids)
    if set(aa) != set(di):
        raise ValueError(f"AA and 3Di record IDs differ for {prefix}")
    for key in aa:
        if len(aa[key]) != len(di[key]):
            raise ValueError(f"AA and 3Di alignment lengths differ for {key} in {prefix}")
        aa_gaps = np.asarray(list(aa[key])) == "-"
        di_gaps = np.asarray(list(di[key])) == "-"
        if not np.array_equal(aa_gaps, di_gaps):
            raise ValueError(f"AA and 3Di gap patterns differ for {key} in {prefix}")
    return aa, di


def build_coordinate_maps(
    aa_records: dict[str, str],
    manifest: pd.DataFrame,
    config: dict,
) -> tuple[dict[str, dict[int, object]], dict[str, list[object]], pd.DataFrame]:
    q = config["quality"]
    maps = {}
    residues_by_id = {}
    qc_rows = []
    manifest_by_id = manifest.set_index("foldmason_id")
    for foldmason_id, aligned_sequence in aa_records.items():
        row = manifest_by_id.loc[foldmason_id]
        _, residues = load_structure_residues(
            row["staged_structure_path"], str(row.get("chain_id", "")) or None
        )
        column_map, qc = map_alignment_sequence_to_structure(
            aligned_sequence,
            residues,
            float(q["min_sequence_mapping_identity"]),
            float(q["min_sequence_mapping_coverage"]),
        )
        maps[foldmason_id] = column_map
        residues_by_id[foldmason_id] = residues
        qc_rows.append({"foldmason_id": foldmason_id, "protein_id": row["protein_id"], **qc})
    return maps, residues_by_id, pd.DataFrame(qc_rows)


def pairwise_distance_matrix(coords: np.ndarray) -> np.ndarray:
    delta = coords[:, None, :] - coords[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=-1))


def aggregate_local_lddt(
    reference_id: str,
    comparison_ids: list[str],
    coordinate_maps: dict[str, dict[int, object]],
    config: dict,
) -> pd.DataFrame:
    ref_map = coordinate_maps[reference_id]
    columns = sorted(ref_map)
    col_to_i = {c: i for i, c in enumerate(columns)}
    ref_coords = np.asarray([ref_map[c].ca for c in columns], dtype=float)
    ref_dist = pairwise_distance_matrix(ref_coords)
    radius = float(config["lddt"]["inclusion_radius_A"])
    thresholds = np.asarray(config["lddt"]["thresholds_A"], dtype=float)
    min_neighbors = int(config["lddt"]["min_neighbors"])
    neighbor_mask = (ref_dist > 0.0) & (ref_dist < radius)
    neighbor_counts = neighbor_mask.sum(axis=1)

    lddt_values: list[list[float]] = [[] for _ in columns]
    coverage_values: list[list[float]] = [[] for _ in columns]
    focal_present = np.zeros(len(columns), dtype=int)

    for alt_id in comparison_ids:
        alt_map = coordinate_maps[alt_id]
        common_cols = [c for c in columns if c in alt_map]
        if not common_cols:
            continue
        common_indices = np.asarray([col_to_i[c] for c in common_cols], dtype=int)
        alt_coords = np.asarray([alt_map[c].ca for c in common_cols], dtype=float)
        alt_dist = pairwise_distance_matrix(alt_coords)
        alt_pos_by_ref_i = {ref_i: pos for pos, ref_i in enumerate(common_indices)}
        common_set = set(common_indices.tolist())
        for ref_i in common_indices:
            focal_present[ref_i] += 1
            total_neighbors = int(neighbor_counts[ref_i])
            if total_neighbors == 0:
                continue
            valid_neighbor_indices = [
                j for j in np.flatnonzero(neighbor_mask[ref_i]) if int(j) in common_set
            ]
            coverage = len(valid_neighbor_indices) / total_neighbors
            coverage_values[ref_i].append(coverage)
            if len(valid_neighbor_indices) < min_neighbors:
                continue
            alt_i = alt_pos_by_ref_i[int(ref_i)]
            ref_d = ref_dist[ref_i, valid_neighbor_indices]
            alt_js = [alt_pos_by_ref_i[int(j)] for j in valid_neighbor_indices]
            alt_d = alt_dist[alt_i, alt_js]
            delta = np.abs(ref_d - alt_d)
            score = float(np.mean(delta[:, None] < thresholds[None, :]))
            lddt_values[ref_i].append(score)

    rows = []
    n_total = len(comparison_ids)
    for i, column in enumerate(columns):
        residue = ref_map[column]
        lvals = np.asarray(lddt_values[i], dtype=float)
        cvals = np.asarray(coverage_values[i], dtype=float)
        rows.append(
            {
                "alignment_column": column + 1,
                "residue_index": residue.seq_index,
                "chain_id": residue.chain_id,
                "resseq": residue.resseq,
                "icode": residue.icode.strip(),
                "residue_label": residue.residue_label,
                "residue_name": residue.resname,
                "aa": residue.aa,
                "local_pLDDT": residue.plddt,
                "lddt": float(np.mean(lvals)) if len(lvals) else np.nan,
                "lddt_sd": float(np.std(lvals, ddof=1)) if len(lvals) > 1 else (0.0 if len(lvals) == 1 else np.nan),
                "lddt_coverage": float(np.mean(cvals)) if len(cvals) else 0.0,
                "presence": (focal_present[i] / n_total) if n_total else np.nan,
                "gap_frequency": (1.0 - focal_present[i] / n_total) if n_total else np.nan,
                "n_lddt_comparisons": len(lvals),
                "n_focal_present": int(focal_present[i]),
                "n_comparators_total": n_total,
                "n_reference_neighbors": int(neighbor_counts[i]),
            }
        )
    return pd.DataFrame(rows)


def kabsch_superpose(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    p = mobile - mobile_center
    q = reference - reference_center
    u, _, vt = np.linalg.svd(p.T @ q)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    translation = reference_center - mobile_center @ rotation
    return rotation, translation, mobile @ rotation + translation


def within_ca_dispersion(
    reference_id: str,
    comparison_ids: list[str],
    coordinate_maps: dict[str, dict[int, object]],
    trim_cutoff_A: float = 4.0,
    trim_iterations: int = 2,
) -> pd.DataFrame:
    ref_map = coordinate_maps[reference_id]
    columns = sorted(ref_map)
    transformed_by_col: dict[int, list[np.ndarray]] = {c: [ref_map[c].ca] for c in columns}
    for alt_id in comparison_ids:
        alt_map = coordinate_maps[alt_id]
        common = [c for c in columns if c in alt_map]
        if len(common) < 10:
            continue
        mobile = np.asarray([alt_map[c].ca for c in common], dtype=float)
        reference = np.asarray([ref_map[c].ca for c in common], dtype=float)
        keep = np.ones(len(common), dtype=bool)
        rotation = np.eye(3)
        translation = np.zeros(3)
        for _ in range(trim_iterations + 1):
            rotation, translation, _ = kabsch_superpose(mobile[keep], reference[keep])
            transformed = mobile @ rotation + translation
            residual = np.linalg.norm(transformed - reference, axis=1)
            new_keep = residual <= trim_cutoff_A
            if new_keep.sum() < 10 or np.array_equal(new_keep, keep):
                break
            keep = new_keep
        transformed = mobile @ rotation + translation
        for c, xyz in zip(common, transformed):
            transformed_by_col[c].append(xyz)

    rows = []
    for c in columns:
        coords = np.asarray(transformed_by_col[c], dtype=float)
        centroid = coords.mean(axis=0)
        sq = np.sum((coords - centroid) ** 2, axis=1)
        rows.append(
            {
                "alignment_column": c + 1,
                "ca_variance_A2": float(np.mean(sq)) if len(coords) > 1 else np.nan,
                "ca_rms_dispersion_A": float(np.sqrt(np.mean(sq))) if len(coords) > 1 else np.nan,
                "n_ca_coordinates": len(coords),
            }
        )
    return pd.DataFrame(rows)


def conservation_by_reference_column(
    reference_id: str,
    selected_ids: list[str],
    aa_records: dict[str, str],
    di_records: dict[str, str],
    reference_map: dict[int, object],
) -> pd.DataFrame:
    rows = []
    for column in sorted(reference_map):
        aa = normalized_entropy_conservation([aa_records[x][column] for x in selected_ids], 20)
        di = normalized_entropy_conservation([di_records[x][column] for x in selected_ids], 20)
        rows.append(
            {
                "alignment_column": column + 1,
                "aa_conservation": aa["conservation"],
                "aa_entropy": aa["entropy"],
                "aa_consensus": aa["dominant"],
                "aa_consensus_fraction": aa["dominant_fraction"],
                "aa_n": aa["n"],
                "within_3Di_conservation": di["conservation"],
                "within_3Di_entropy": di["entropy"],
                "within_3Di_consensus": di["dominant"],
                "within_3Di_consensus_fraction": di["dominant_fraction"],
                "within_3Di_n": di["n"],
            }
        )
    return pd.DataFrame(rows)


def select_comparators(
    focal_medoid_id: str,
    medoids_in_group: pd.DataFrame,
    pairwise_standard: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    tm = pair_lookup(pairwise_standard, "symmetric_tm")
    cov = pair_lookup(pairwise_standard, "reciprocal_coverage")
    scope = config["between_scope"]
    rows = []
    for row in medoids_in_group.itertuples(index=False):
        if str(row.medoid_id) == str(focal_medoid_id):
            continue
        score = lookup_pair(tm, focal_medoid_id, str(row.medoid_id))
        coverage = lookup_pair(cov, focal_medoid_id, str(row.medoid_id))
        passes_tm = np.isfinite(score) and score >= float(scope["min_symmetric_tm"])
        passes_cov = (not np.isfinite(coverage)) or coverage >= float(scope["min_reciprocal_coverage"])
        rows.append(
            {
                "medoid_id": str(row.medoid_id),
                "foldmason_id": str(row.foldmason_id),
                "cluster_id": str(row.cluster_id),
                "symmetric_tm": score,
                "reciprocal_coverage": coverage,
                "eligible": bool(passes_tm and passes_cov),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out[out["eligible"]].sort_values(["symmetric_tm", "medoid_id"], ascending=[False, True])
    return out.head(int(scope["nearest_k"])).reset_index(drop=True)


def between_3di_metrics(
    focal_id: str,
    comparator_ids: list[str],
    focal_residue_index_to_column: dict[int, int],
    between_di: dict[str, str],
    within_metrics: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    within_by_index = within_metrics.set_index("residue_index")
    for residue_index, column in focal_residue_index_to_column.items():
        state = str(within_by_index.loc[residue_index, "within_3Di_consensus"])
        states = [between_di[x][column] for x in comparator_ids if between_di[x][column] not in "-."]
        if states and state not in ("", "nan"):
            background_frequency = sum(x == state for x in states) / len(states)
            divergence = 1.0 - background_frequency
        else:
            background_frequency = np.nan
            divergence = np.nan
        entropy = normalized_entropy_conservation(states, 20)
        rows.append(
            {
                "residue_index": residue_index,
                "between_3Di_alignment_column": column + 1,
                "between_3Di_background_frequency": background_frequency,
                "between_3Di_divergence": divergence,
                "between_3Di_entropy": entropy["entropy"],
                "between_3Di_conservation": entropy["conservation"],
                "between_3Di_n": len(states),
            }
        )
    return pd.DataFrame(rows)


def prefix_columns(df: pd.DataFrame, prefix: str, keep: set[str]) -> pd.DataFrame:
    return df.rename(columns={c: f"{prefix}{c}" for c in df.columns if c not in keep})


def add_scores_and_categories(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    df = df.copy()
    scores = config["scores"]
    floor = float(scores["plddt_floor"])
    full = float(scores["plddt_full_weight"])
    df["pLDDT_weight"] = np.clip((df["local_pLDDT"] - floor) / (full - floor), 0.0, 1.0)
    min_comp = int(config["between_scope"]["min_comparators_for_call"])
    df["between_evidence_weight"] = (
        np.sqrt(df["between_lddt_coverage"].clip(0.0, 1.0))
        * np.minimum(1.0, df["between_n_lddt_comparisons"] / min_comp)
    )
    raw = (
        df["within_lddt"]
        * (1.0 - df["between_lddt"])
        * df["within_lddt_coverage"]
        * df["pLDDT_weight"]
        * (1.0 - df["within_gap_frequency"])
    )
    df["cluster_specific_constraint_score_raw"] = raw
    if bool(scores["include_between_evidence_weight"]):
        df["cluster_specific_constraint_score"] = raw * df["between_evidence_weight"]
    else:
        df["cluster_specific_constraint_score"] = raw
    df["cluster_specific_3Di_score"] = (
        df["within_3Di_conservation"]
        * df["between_3Di_divergence"]
        * df["within_presence"]
        * np.minimum(1.0, df["between_3Di_n"] / min_comp)
    )
    df["structured_insertion_score"] = np.where(
        df["between_n_comparators_total"] >= min_comp,
        df["within_lddt"]
        * df["within_lddt_coverage"]
        * df["pLDDT_weight"]
        * (1.0 - df["within_gap_frequency"])
        * (1.0 - df["between_presence"].fillna(0.0)),
        np.nan,
    )

    c = config["classification"]
    quality_within = (
        (df["local_pLDDT"] >= float(config["quality"]["min_local_plddt"]))
        & (df["within_presence"] >= float(c["min_within_presence"]))
        & (df["within_lddt_coverage"] >= float(c["min_within_lddt_coverage"]))
        & (df["within_gap_frequency"] <= float(c["max_within_gap_frequency"]))
    )
    between_evaluable = (
        (df["between_n_comparators_total"] >= min_comp)
        & (df["between_presence"] >= float(c["min_between_presence"]))
        & (df["between_lddt_coverage"] >= float(c["min_between_lddt_coverage"]))
        & df["between_lddt"].notna()
    )
    df["passes_within_quality"] = quality_within
    df["passes_between_evidence"] = between_evaluable
    df["passes_constraint_call_qc"] = quality_within & between_evaluable

    within_high = df["within_lddt"] >= float(c["within_high"])
    between_high = df["between_lddt"] >= float(c["between_high"])
    category = np.full(len(df), "insufficient evidence", dtype=object)
    usable = df["passes_constraint_call_qc"].to_numpy()
    wh = within_high.to_numpy()
    bh = between_high.to_numpy()
    category[usable & wh & bh] = "shared GT core"
    category[usable & wh & ~bh] = "cluster-specific structural constraint"
    category[usable & ~wh & ~bh] = "variable or poorly constrained"
    category[usable & ~wh & bh] = "ambiguous or alignment artifact"
    df["structural_category"] = category
    df["structured_insertion_candidate"] = (
        quality_within
        & within_high
        & (df["between_n_comparators_total"] >= min_comp)
        & (df["between_presence"] < float(c["min_between_presence"]))
    )
    df["sequence_divergent_but_structurally_constrained"] = (
        (df["aa_conservation"] < 0.50)
        & (df["within_3Di_conservation"] >= 0.80)
        & within_high
        & quality_within
    )
    return df


def attach_active_site_annotations(df: pd.DataFrame, annotations: pd.DataFrame | None) -> pd.DataFrame:
    df = df.copy()
    df["residue_annotation"] = ""
    df["structural_role"] = ""
    if annotations is None or annotations.empty:
        return df
    require_columns(annotations, ["protein_id"], "active_site_residues.tsv")
    ann = annotations.copy()
    if "residue_index" in ann.columns:
        ann["residue_index"] = pd.to_numeric(ann["residue_index"], errors="coerce")
        keys = ["protein_id", "residue_index"]
    elif "resseq" in ann.columns:
        ann["resseq"] = pd.to_numeric(ann["resseq"], errors="coerce")
        keys = ["protein_id", "resseq"]
        if "chain_id" in ann.columns:
            keys.append("chain_id")
    else:
        raise ValueError("active_site_residues.tsv needs residue_index or resseq")
    label_col = next((c for c in ("annotation", "role", "residue_role") if c in ann.columns), None)
    if label_col is None:
        ann["annotation"] = "annotated residue"
        label_col = "annotation"
    role_col = "structural_role" if "structural_role" in ann.columns else None
    keep = keys + [label_col] + ([role_col] if role_col else [])
    ann = ann[keep].rename(columns={label_col: "_annotation", **({role_col: "_role"} if role_col else {})})
    aggregations = {
        "_annotation": lambda s: ";".join(sorted(set(str(x) for x in s if str(x)))),
    }
    if role_col:
        aggregations["_role"] = lambda s: ";".join(sorted(set(str(x) for x in s if str(x))))
    ann = ann.groupby(keys, as_index=False, dropna=False).agg(aggregations)
    merged = df.merge(ann, on=keys, how="left")
    merged["residue_annotation"] = merged["_annotation"].fillna("")
    merged["structural_role"] = merged.get("_role", "").fillna("") if "_role" in merged else ""
    return merged.drop(columns=[c for c in ("_annotation", "_role") if c in merged])


def contiguous_groups(indices: list[int], max_bridge: int) -> list[list[int]]:
    if not indices:
        return []
    groups = [[indices[0]]]
    for value in indices[1:]:
        if value - groups[-1][-1] <= max_bridge + 1:
            groups[-1].append(value)
        else:
            groups.append([value])
    return groups


def summarize_regions(residues: pd.DataFrame, config: dict) -> pd.DataFrame:
    region_cfg = config["regions"]
    rows = []
    for (cluster_id, medoid_id), cdf in residues.groupby(["cluster_id", "medoid_id"], sort=True):
        cdf = cdf.sort_values("residue_index").copy()
        constraint_candidate = (
            (cdf["structural_category"] == "cluster-specific structural constraint")
            & (cdf["cluster_specific_constraint_score"] >= float(region_cfg["min_constraint_score"]))
        )
        candidate_sets = [
            ("divergent_aligned_region", cdf.loc[constraint_candidate, "residue_index"].astype(int).tolist()),
            ("structured_cluster_specific_insertion", cdf.loc[cdf["structured_insertion_candidate"], "residue_index"].astype(int).tolist()),
        ]
        for region_type, indices in candidate_sets:
            for group in contiguous_groups(indices, int(region_cfg["max_bridge_residues"])):
                start, end = min(group), max(group)
                window = cdf[(cdf["residue_index"] >= start) & (cdf["residue_index"] <= end)]
                candidate_fraction = len(group) / len(window)
                if len(window) < int(region_cfg["min_region_length"]):
                    continue
                if candidate_fraction < float(region_cfg["min_candidate_fraction"]):
                    continue
                annotated = window[window["residue_annotation"].astype(str) != ""]
                nearby = cdf[
                    (cdf["residue_index"] >= start - int(region_cfg["active_site_window"]))
                    & (cdf["residue_index"] <= end + int(region_cfg["active_site_window"]))
                    & (cdf["residue_annotation"].astype(str) != "")
                ]
                roles = sorted(set(x for x in window["structural_role"].astype(str) if x))
                first = window.iloc[0]
                last = window.iloc[-1]
                rows.append(
                    {
                        "cluster_id": cluster_id,
                        "medoid_id": medoid_id,
                        "region_type": region_type,
                        "residue_index_start": start,
                        "residue_index_end": end,
                        "residue_range": f"{first.chain_id}:{int(first.resseq)}{first.icode}-{last.chain_id}:{int(last.resseq)}{last.icode}",
                        "region_length": len(window),
                        "candidate_fraction": candidate_fraction,
                        "mean_within_lddt": window["within_lddt"].mean(),
                        "mean_between_lddt": window["between_lddt"].mean(),
                        "mean_cluster_specific_constraint_score": window["cluster_specific_constraint_score"].mean(),
                        "mean_structured_insertion_score": window["structured_insertion_score"].mean(),
                        "mean_AA_conservation": window["aa_conservation"].mean(),
                        "mean_3Di_conservation": window["within_3Di_conservation"].mean(),
                        "mean_cluster_specific_3Di_score": window["cluster_specific_3Di_score"].mean(),
                        "mean_pLDDT": window["local_pLDDT"].mean(),
                        "mean_CA_rms_dispersion_A": window["ca_rms_dispersion_A"].mean(),
                        "structural_role": ";".join(roles) if roles else "unassigned",
                        "annotated_residues_in_region": ";".join(
                            f"{r.residue_label}:{r.residue_annotation}" for r in annotated.itertuples()
                        ),
                        "nearby_active_site_or_binding_residues": ";".join(
                            f"{r.residue_label}:{r.residue_annotation}" for r in nearby.itertuples()
                        ),
                    }
                )
    return pd.DataFrame(rows)


def score_digit(value: float) -> str:
    if not np.isfinite(value):
        return " "
    return str(int(np.clip(np.floor(value * 10), 0, 9)))


def write_alignment_blocks(
    regions: pd.DataFrame,
    residue_metrics: pd.DataFrame,
    alignment_context: dict,
    outdir: Path,
    config: dict,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    top_n = int(config["regions"]["top_alignment_blocks_per_cluster"])
    if regions.empty:
        return
    rank_column = "mean_cluster_specific_constraint_score"
    for cluster_id, rdf in regions.groupby("cluster_id", sort=True):
        context = alignment_context[str(cluster_id)]
        selected = rdf.sort_values(rank_column, ascending=False).head(top_n)
        lines = []
        metrics = residue_metrics[residue_metrics["cluster_id"] == cluster_id].set_index("residue_index")
        within_ref_map = context["within_ref_map"]
        residue_to_within_col = {r.seq_index: c for c, r in within_ref_map.items()}
        between_ref_map = context["between_ref_map"]
        residue_to_between_col = {r.seq_index: c for c, r in between_ref_map.items()}
        for region in selected.itertuples(index=False):
            lines.append(
                f"## {cluster_id} {region.region_type} residues {region.residue_index_start}-{region.residue_index_end} "
                f"score={region.mean_cluster_specific_constraint_score:.3f}"
            )
            within_cols = [
                residue_to_within_col[i]
                for i in range(int(region.residue_index_start), int(region.residue_index_end) + 1)
                if i in residue_to_within_col
            ]
            if within_cols:
                lo, hi = min(within_cols), max(within_cols) + 1
                lines.append("# within-cluster AA")
                for ident, seq in context["within_aa"].items():
                    lines.append(f"AA  {ident:<28} {seq[lo:hi]}")
                lines.append("# within-cluster 3Di")
                for ident, seq in context["within_di"].items():
                    lines.append(f"3DI {ident:<28} {seq[lo:hi]}")
                annotation = []
                lddt_line = []
                score_line = []
                for col in range(lo, hi):
                    residue = within_ref_map.get(col)
                    if residue is None or residue.seq_index not in metrics.index:
                        annotation.append(" ")
                        lddt_line.append(" ")
                        score_line.append(" ")
                    else:
                        row = metrics.loc[residue.seq_index]
                        annotation.append("C" if row["structural_category"] == "cluster-specific structural constraint" else ".")
                        lddt_line.append(score_digit(float(row["within_lddt"])))
                        score_line.append(score_digit(float(row["cluster_specific_constraint_score"])))
                lines.append(f"ANN {'':<28} {''.join(annotation)}  C=constraint call")
                lines.append(f"WLD {'':<28} {''.join(lddt_line)}  0-9 within-cluster lDDT")
                lines.append(f"SCR {'':<28} {''.join(score_line)}  0-9 final score")

            between_cols = [
                residue_to_between_col[i]
                for i in range(int(region.residue_index_start), int(region.residue_index_end) + 1)
                if i in residue_to_between_col
            ]
            if between_cols:
                lo, hi = min(between_cols), max(between_cols) + 1
                display_ids = [context["focal_foldmason_id"], *context["comparator_foldmason_ids"]]
                lines.append("# between-cluster medoid AA")
                for ident in display_ids:
                    lines.append(f"AA  {ident:<28} {context['between_aa'][ident][lo:hi]}")
                lines.append("# between-cluster medoid 3Di")
                for ident in display_ids:
                    lines.append(f"3DI {ident:<28} {context['between_di'][ident][lo:hi]}")
            lines.append("")
        (outdir / f"{sanitize_group(cluster_id)}_top_regions.txt").write_text("\n".join(lines) + "\n")


def write_bfactor_maps(cluster_df: pd.DataFrame, structure_path: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics = {
        "within_lddt": "within_lddt",
        "between_lddt": "between_lddt",
        "constraint_score": "cluster_specific_constraint_score",
        "aa_conservation": "aa_conservation",
        "3di_conservation": "within_3Di_conservation",
        "3di_signature": "cluster_specific_3Di_score",
    }
    created = []
    for label, column in metrics.items():
        values = {
            (str(r.chain_id), int(r.resseq), str(r.icode or "")): float(getattr(r, column))
            for r in cluster_df.itertuples(index=False)
            if pd.notna(getattr(r, column))
        }
        path = output_dir / f"{label}.pdb"
        write_metric_pdb(structure_path, path, values, scale=100.0)
        created.append((label, path))
    pml = ["bg_color white", "set ray_opaque_background, off"]
    colors = {
        "within_lddt": "blue_white_red",
        "between_lddt": "red_white_blue",
        "constraint_score": "white_yellow_red",
        "aa_conservation": "white_cyan_blue",
        "3di_conservation": "white_green",
        "3di_signature": "white_magenta",
    }
    for label, path in created:
        object_name = label.replace("-", "_")
        pml.extend(
            [
                f"load {path.name}, {object_name}",
                f"hide everything, {object_name}",
                f"show cartoon, {object_name}",
                f"spectrum b, {colors[label]}, {object_name}, minimum=0, maximum=100",
            ]
        )
    pml.append("disable all")
    pml.append("enable constraint_score")
    pml.append("orient constraint_score")
    (output_dir / "view_metrics.pml").write_text("\n".join(pml) + "\n")


def compute_constraint_workflow(
    prepared_members_path: str | Path,
    medoids_path: str | Path,
    pairwise_path: str | Path,
    jobs_path: str | Path,
    output_dir: str | Path,
    config: dict,
    active_site_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    members = read_tsv(prepared_members_path)
    medoids = read_tsv(medoids_path)
    jobs = read_tsv(jobs_path)
    pairwise_raw = translate_pairwise_ids(read_tsv(pairwise_path), members)
    pairwise_standard = standardize_pairwise(pairwise_raw)
    annotations = read_tsv(active_site_path) if active_site_path else None
    require_columns(
        members,
        ["protein_id", "cluster_id", "analysis_group", "foldmason_id", "staged_structure_path", "selected_for_within_alignment", "is_medoid"],
        "prepared_members.tsv",
    )
    for col in ("selected_for_within_alignment", "is_medoid"):
        members[col] = members[col].astype(str).str.lower().isin({"true", "1", "yes"})
    jobs_index = {}
    for row in jobs.itertuples(index=False):
        key = (str(row.scope), str(row.cluster_id) if str(row.scope) == "within" else str(row.analysis_group))
        jobs_index[key] = str(row.output_prefix)

    residue_tables = []
    mapping_qc_tables = []
    comparator_tables = []
    alignment_context = {}

    for medoid in medoids.itertuples(index=False):
        cluster_id = str(medoid.cluster_id)
        group = str(medoid.analysis_group)
        focal_protein_id = str(medoid.medoid_id)
        focal_foldmason_id = str(medoid.foldmason_id)
        within_manifest = members[(members["cluster_id"] == cluster_id) & members["selected_for_within_alignment"]]
        if len(within_manifest) < 2:
            print(f"SKIP cluster {cluster_id}: fewer than two aligned within-cluster structures")
            continue
        within_ids = within_manifest["foldmason_id"].astype(str).tolist()
        within_aa, within_di = load_alignment_pair(jobs_index[("within", cluster_id)], within_ids)
        within_maps, _, mapping_qc = build_coordinate_maps(within_aa, within_manifest, config)
        mapping_qc["scope"] = "within"
        mapping_qc["cluster_id"] = cluster_id
        mapping_qc_tables.append(mapping_qc)
        within_comparison_ids = [x for x in within_ids if x != focal_foldmason_id]
        within = aggregate_local_lddt(focal_foldmason_id, within_comparison_ids, within_maps, config)
        within = prefix_columns(
            within,
            "within_",
            {"alignment_column", "residue_index", "chain_id", "resseq", "icode", "residue_label", "residue_name", "aa", "local_pLDDT"},
        )
        ca = within_ca_dispersion(focal_foldmason_id, within_comparison_ids, within_maps)
        conservation = conservation_by_reference_column(
            focal_foldmason_id, within_ids, within_aa, within_di, within_maps[focal_foldmason_id]
        )
        table = within.merge(ca, on="alignment_column", how="left").merge(
            conservation, on="alignment_column", how="left"
        )

        between_manifest = members[(members["analysis_group"] == group) & members["is_medoid"]]
        between_ids = between_manifest["foldmason_id"].astype(str).tolist()
        if len(between_manifest) < 2:
            raise ValueError(f"Analysis group {group} has fewer than two medoids")
        between_aa, between_di = load_alignment_pair(jobs_index[("between", group)], between_ids)
        between_maps, _, mapping_qc = build_coordinate_maps(between_aa, between_manifest, config)
        mapping_qc["scope"] = "between"
        mapping_qc["cluster_id"] = cluster_id
        mapping_qc_tables.append(mapping_qc[mapping_qc["foldmason_id"] == focal_foldmason_id])
        medoids_in_group = medoids[medoids["analysis_group"] == group]
        comparator_df = select_comparators(focal_protein_id, medoids_in_group, pairwise_standard, config)
        if not comparator_df.empty:
            comparator_df["focal_cluster_id"] = cluster_id
            comparator_df["focal_medoid_id"] = focal_protein_id
            comparator_tables.append(comparator_df)
        comparator_ids = comparator_df["foldmason_id"].tolist() if not comparator_df.empty else []
        between = aggregate_local_lddt(focal_foldmason_id, comparator_ids, between_maps, config)
        between = prefix_columns(
            between,
            "between_",
            {"residue_index", "chain_id", "resseq", "icode", "residue_label", "residue_name", "aa", "local_pLDDT"},
        )
        between = between.drop(
            columns=[c for c in ("chain_id", "resseq", "icode", "residue_label", "residue_name", "aa", "local_pLDDT") if c in between],
            errors="ignore",
        )
        table = table.merge(between, on="residue_index", how="left")

        between_ref_map = between_maps[focal_foldmason_id]
        residue_to_between_col = {r.seq_index: c for c, r in between_ref_map.items()}
        di_metrics = between_3di_metrics(
            focal_foldmason_id,
            comparator_ids,
            residue_to_between_col,
            between_di,
            table,
        )
        table = table.merge(di_metrics, on="residue_index", how="left")
        table.insert(0, "medoid_id", focal_protein_id)
        table.insert(0, "cluster_id", cluster_id)
        table.insert(2, "analysis_group", group)
        table["n_cluster_members"] = int(medoid.cluster_size)
        table["n_within_alignment_members"] = len(within_ids)
        table["n_eligible_between_comparators"] = len(comparator_ids)
        table = add_scores_and_categories(table, config)
        table = attach_active_site_annotations(table, annotations)
        residue_tables.append(table)

        source_structure = str(medoid.structure_path)
        write_bfactor_maps(
            table,
            source_structure,
            output_dir / "mapped_structures" / sanitize_group(cluster_id),
        )
        alignment_context[cluster_id] = {
            "within_aa": within_aa,
            "within_di": within_di,
            "within_ref_map": within_maps[focal_foldmason_id],
            "between_aa": between_aa,
            "between_di": between_di,
            "between_ref_map": between_ref_map,
            "focal_foldmason_id": focal_foldmason_id,
            "comparator_foldmason_ids": comparator_ids,
        }

    if not residue_tables:
        raise ValueError("No clusters produced residue metrics")
    residues = pd.concat(residue_tables, ignore_index=True)
    mapping_qc = pd.concat(mapping_qc_tables, ignore_index=True) if mapping_qc_tables else pd.DataFrame()
    comparators = pd.concat(comparator_tables, ignore_index=True) if comparator_tables else pd.DataFrame()
    regions = summarize_regions(residues, config)

    summary_rows = []
    medoid_index = medoids.set_index("cluster_id")
    for cluster_id, cdf in residues.groupby("cluster_id", sort=True):
        m = medoid_index.loc[cluster_id]
        member_cluster = members[members["cluster_id"] == cluster_id]
        organism_col = next((c for c in ("organism_type", "superkingdom", "source_domain") if c in member_cluster.columns), None)
        viral_fraction = (
            member_cluster[organism_col].astype(str).str.contains("virus", case=False, regex=False).mean()
            if organism_col else 1.0
        )
        summary_rows.append(
            {
                "cluster_id": cluster_id,
                "analysis_group": m["analysis_group"],
                "medoid_id": m["medoid_id"],
                "cluster_size": int(m["cluster_size"]),
                "median_within_TM": pd.to_numeric(m["median_within_TM"], errors="coerce"),
                "median_within_lddt": cdf["within_lddt"].median(),
                "median_between_lddt": cdf["between_lddt"].median(),
                "mean_cluster_specific_score": cdf["cluster_specific_constraint_score"].mean(),
                "mean_AA_conservation": cdf["aa_conservation"].mean(),
                "mean_3Di_conservation": cdf["within_3Di_conservation"].mean(),
                "annotation_fraction": (cdf["residue_annotation"].astype(str) != "").mean(),
                "viral_fraction": viral_fraction,
                "mean_pLDDT": cdf["local_pLDDT"].mean(),
                "constraint_residue_fraction": (cdf["structural_category"] == "cluster-specific structural constraint").mean(),
                "shared_core_residue_fraction": (cdf["structural_category"] == "shared GT core").mean(),
                "n_ranked_regions": int((regions["cluster_id"] == cluster_id).sum()) if not regions.empty else 0,
            }
        )
    summary = pd.DataFrame(summary_rows)

    residues.to_csv(output_dir / "all_residue_metrics.tsv", sep="\t", index=False, na_rep="NA")
    regions.to_csv(output_dir / "ranked_regions.tsv", sep="\t", index=False, na_rep="NA")
    summary.to_csv(output_dir / "cluster_summary.tsv", sep="\t", index=False, na_rep="NA")
    mapping_qc.to_csv(output_dir / "alignment_mapping_qc.tsv", sep="\t", index=False, na_rep="NA")
    comparators.to_csv(output_dir / "between_comparators.tsv", sep="\t", index=False, na_rep="NA")
    counts = (
        residues.groupby(["cluster_id", "structural_category"], dropna=False)
        .size()
        .rename("n_residues")
        .reset_index()
    )
    counts.to_csv(output_dir / "structural_category_counts.tsv", sep="\t", index=False)
    for cluster_id, cdf in residues.groupby("cluster_id"):
        cluster_dir = output_dir / "per_cluster"
        cluster_dir.mkdir(parents=True, exist_ok=True)
        cdf.to_csv(cluster_dir / f"{sanitize_group(cluster_id)}_residue_metrics.tsv", sep="\t", index=False, na_rep="NA")
    write_alignment_blocks(regions, residues, alignment_context, output_dir / "alignment_blocks", config)
    return residues, regions, summary, mapping_qc
