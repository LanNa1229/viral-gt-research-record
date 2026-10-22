from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .common import (
    json_dumps,
    load_config,
    read_tsv,
    require_columns,
    safe_id,
    sanitize_group,
    structure_qc,
)


QUERY_TM_COLUMNS = ("qtmscore", "qtm", "q_tm", "query_tm")
TARGET_TM_COLUMNS = ("ttmscore", "ttm", "t_tm", "target_tm")
SINGLE_TM_COLUMNS = ("tm_score", "alntmscore", "tm", "tmscore")
QUERY_COV_COLUMNS = ("qcov", "query_coverage", "q_cov")
TARGET_COV_COLUMNS = ("tcov", "target_coverage", "t_cov")


def first_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace("", np.nan), errors="coerce")


def standardize_pairwise(df: pd.DataFrame) -> pd.DataFrame:
    require_columns(df, ["query", "target"], "pairwise_structural_scores.tsv")
    qtm = first_column(df, QUERY_TM_COLUMNS)
    ttm = first_column(df, TARGET_TM_COLUMNS)
    stm = first_column(df, SINGLE_TM_COLUMNS)
    if qtm and ttm:
        q = numeric(df[qtm])
        t = numeric(df[ttm])
        symmetric_tm = np.minimum(q, t)
        mean_directional_tm = (q + t) / 2.0
    elif stm:
        symmetric_tm = numeric(df[stm])
        mean_directional_tm = symmetric_tm
    else:
        raise ValueError(
            "Pairwise table needs qtmscore+ttmscore or one of: " + ", ".join(SINGLE_TM_COLUMNS)
        )

    qcov = first_column(df, QUERY_COV_COLUMNS)
    tcov = first_column(df, TARGET_COV_COLUMNS)
    if qcov and tcov:
        reciprocal_coverage = np.minimum(numeric(df[qcov]), numeric(df[tcov]))
    elif qcov:
        reciprocal_coverage = numeric(df[qcov])
    elif tcov:
        reciprocal_coverage = numeric(df[tcov])
    else:
        reciprocal_coverage = np.nan

    out = pd.DataFrame(
        {
            "query": df["query"].astype(str),
            "target": df["target"].astype(str),
            "symmetric_tm": symmetric_tm,
            "mean_directional_tm": mean_directional_tm,
            "reciprocal_coverage": reciprocal_coverage,
        }
    )
    out = out[out["query"] != out["target"]].dropna(subset=["symmetric_tm"])
    out["pair_a"] = out[["query", "target"]].min(axis=1)
    out["pair_b"] = out[["query", "target"]].max(axis=1)
    return (
        out.groupby(["pair_a", "pair_b"], as_index=False)
        .agg(
            symmetric_tm=("symmetric_tm", "mean"),
            mean_directional_tm=("mean_directional_tm", "mean"),
            reciprocal_coverage=("reciprocal_coverage", "mean"),
            n_directional_rows=("query", "size"),
        )
    )


def translate_pairwise_ids(df: pd.DataFrame, members: pd.DataFrame) -> pd.DataFrame:
    """Translate staged FoldMason/Foldseek IDs back to biological protein IDs."""
    require_columns(df, ["query", "target"], "pairwise_structural_scores.tsv")
    mapping = dict(zip(members["foldmason_id"].astype(str), members["protein_id"].astype(str)))
    biological = set(members["protein_id"].astype(str))
    out = df.copy()
    for column in ("query", "target"):
        normalized = out[column].astype(str).map(lambda x: Path(x).name)
        for suffix in (".pdb.gz", ".cif.gz", ".mmcif.gz", ".pdb", ".cif", ".mmcif", ".gz"):
            normalized = normalized.map(lambda x, s=suffix: x[: -len(s)] if x.lower().endswith(s) else x)
        out[column] = normalized.map(lambda x: mapping.get(x, x))
    observed = set(out["query"]) | set(out["target"])
    unmatched = sorted(observed - biological)
    if unmatched:
        raise ValueError(
            "Pairwise IDs do not match protein_id or staged foldmason_id; first unmatched IDs: "
            + ", ".join(unmatched[:5])
        )
    return out


def pair_lookup(pairwise: pd.DataFrame, value_column: str = "symmetric_tm") -> dict[tuple[str, str], float]:
    return {
        (str(r.pair_a), str(r.pair_b)): float(getattr(r, value_column))
        for r in pairwise.itertuples()
        if pd.notna(getattr(r, value_column))
    }


def lookup_pair(lookup: dict[tuple[str, str], float], a: str, b: str, missing: float = np.nan) -> float:
    key = tuple(sorted((str(a), str(b))))
    return lookup.get(key, missing)


def load_and_qc_members(members_path: str | Path, config: dict) -> pd.DataFrame:
    members = read_tsv(members_path)
    require_columns(members, ["protein_id", "cluster_id"], "cluster_members.tsv")
    path_column = "structure_path"
    if path_column not in members.columns:
        for candidate in ("pdb_path", "cif_path"):
            if candidate in members.columns:
                members[path_column] = members[candidate]
                break
    require_columns(members, [path_column], "cluster_members.tsv")
    if members["protein_id"].duplicated().any():
        dup = members.loc[members["protein_id"].duplicated(), "protein_id"].tolist()[:5]
        raise ValueError(f"protein_id must be unique; duplicates include: {dup}")
    members["analysis_group"] = members.get("analysis_group", "ALL").replace("", "ALL")
    members["chain_id"] = members.get("chain_id", "").fillna("")

    qc_rows = []
    for row in members.itertuples(index=False):
        path = Path(str(row.structure_path)).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Structure not found for {row.protein_id}: {path}")
        qc = structure_qc(path, str(row.chain_id) or None)
        qc_rows.append({"protein_id": str(row.protein_id), "structure_path": str(path), **qc})
    qc_df = pd.DataFrame(qc_rows)
    members = members.drop(columns=["structure_path"]).merge(qc_df, on="protein_id", how="left")

    supplied_plddt_col = next((c for c in ("mean_pLDDT", "pLDDT", "mean_plddt") if c in members.columns), None)
    if supplied_plddt_col:
        members["mean_pLDDT"] = numeric(members[supplied_plddt_col])
        members["confidence_source"] = f"table:{supplied_plddt_col}"
    else:
        members["mean_pLDDT"] = numeric(members["mean_structure_bfactor"])
        members["confidence_source"] = "structure_CA_B_factor"

    q = config["quality"]
    members["passes_model_plddt"] = members["mean_pLDDT"] >= float(q["min_model_plddt"])
    members["passes_min_length"] = members["n_residues"] >= int(q["min_model_residues"])
    if "exclude" in members.columns:
        excluded = members["exclude"].astype(str).str.lower().isin({"1", "true", "yes", "y"})
    else:
        excluded = pd.Series(False, index=members.index)
    members["explicitly_excluded"] = excluded
    members["eligible_for_medoid"] = (
        members["passes_model_plddt"] & members["passes_min_length"] & ~members["explicitly_excluded"]
    )
    members["foldmason_id"] = members["protein_id"].map(safe_id)
    return members


def choose_medoids(
    members: pd.DataFrame,
    pairwise: pd.DataFrame,
    config: dict,
    manual_representatives: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tm_lookup = pair_lookup(pairwise, "symmetric_tm")
    manual_map: dict[str, str] = {}
    if manual_representatives is not None:
        require_columns(manual_representatives, ["cluster_id", "medoid_id"], "manual_representatives.tsv")
        manual_map = dict(zip(manual_representatives["cluster_id"], manual_representatives["medoid_id"]))

    centrality_rows = []
    medoid_rows = []
    for cluster_id, cluster_df in members.groupby("cluster_id", sort=True):
        eligible = cluster_df[cluster_df["eligible_for_medoid"]].copy()
        if eligible.empty:
            raise ValueError(f"Cluster {cluster_id} has no models passing medoid QC")
        ids = eligible["protein_id"].astype(str).tolist()
        for protein_id in ids:
            others = [x for x in ids if x != protein_id]
            observed = [lookup_pair(tm_lookup, protein_id, x) for x in others]
            observed_non_na = [x for x in observed if np.isfinite(x)]
            complete = [0.0 if not np.isfinite(x) else x for x in observed]
            centrality_rows.append(
                {
                    "cluster_id": cluster_id,
                    "protein_id": protein_id,
                    "eligible_cluster_size": len(ids),
                    "n_expected_pairs": len(others),
                    "n_observed_pairs": len(observed_non_na),
                    "pairwise_observation_fraction": (len(observed_non_na) / len(others)) if others else 1.0,
                    "median_within_TM": float(np.median(complete)) if complete else 1.0,
                    "mean_within_TM": float(np.mean(complete)) if complete else 1.0,
                    "median_observed_within_TM": float(np.median(observed_non_na)) if observed_non_na else (1.0 if not others else np.nan),
                }
            )
        cdf = pd.DataFrame([x for x in centrality_rows if x["cluster_id"] == cluster_id])
        if str(cluster_id) in manual_map:
            selected = str(manual_map[str(cluster_id)])
            if selected not in set(cluster_df["protein_id"]):
                raise ValueError(f"Manual medoid {selected} is not a member of cluster {cluster_id}")
            selection_mode = "manual"
        else:
            selected = (
                cdf.sort_values(
                    ["median_within_TM", "mean_within_TM", "pairwise_observation_fraction", "protein_id"],
                    ascending=[False, False, False, True],
                )
                .iloc[0]["protein_id"]
            )
            selection_mode = "highest_median_symmetric_TM"
        selected_member = cluster_df.loc[cluster_df["protein_id"] == selected].iloc[0]
        selected_centrality = cdf.loc[cdf["protein_id"] == selected].iloc[0]
        medoid_rows.append(
            {
                "cluster_id": cluster_id,
                "analysis_group": selected_member["analysis_group"],
                "medoid_id": selected,
                "foldmason_id": selected_member["foldmason_id"],
                "cluster_size": len(cluster_df),
                "eligible_cluster_size": len(eligible),
                "median_within_TM": selected_centrality["median_within_TM"],
                "mean_within_TM": selected_centrality["mean_within_TM"],
                "pairwise_observation_fraction": selected_centrality["pairwise_observation_fraction"],
                "mean_pLDDT": selected_member["mean_pLDDT"],
                "n_residues": selected_member["n_residues"],
                "structure_path": selected_member["structure_path"],
                "chain_id": selected_member["chain_id"],
                "selection_mode": selection_mode,
            }
        )
    return pd.DataFrame(medoid_rows), pd.DataFrame(centrality_rows)


def maximin_sample(ids: list[str], medoid_id: str, n: int, tm_lookup: dict[tuple[str, str], float]) -> list[str]:
    if len(ids) <= n:
        return ids
    selected = [medoid_id]
    remaining = sorted(set(ids) - {medoid_id})
    while remaining and len(selected) < n:
        def min_distance(candidate: str) -> float:
            distances = [1.0 - lookup_pair(tm_lookup, candidate, x, missing=0.0) for x in selected]
            return min(distances)

        best = max(remaining, key=lambda x: (min_distance(x), x))
        selected.append(best)
        remaining.remove(best)
    return selected


def suffix_for_structure(path: Path) -> str:
    lower = path.name.lower()
    for suffix in (".pdb.gz", ".cif.gz", ".mmcif.gz", ".pdb", ".cif", ".mmcif"):
        if lower.endswith(suffix):
            return suffix
    raise ValueError(f"Unsupported structure suffix: {path}")


def stage_member_structures(members: pd.DataFrame, outdir: str | Path) -> pd.DataFrame:
    outdir = Path(outdir)
    structure_dir = outdir / "structures"
    structure_dir.mkdir(parents=True, exist_ok=True)
    members = members.copy()
    staged_paths = []
    for row in members.itertuples(index=False):
        source = Path(str(row.structure_path))
        destination = structure_dir / f"{row.foldmason_id}{suffix_for_structure(source)}"
        if not destination.exists() or destination.stat().st_size != source.stat().st_size:
            shutil.copy2(source, destination)
        staged_paths.append(str(destination.resolve()))
    members["initial_staged_structure_path"] = staged_paths
    return members


def prepare_alignment_jobs(
    members: pd.DataFrame,
    medoids: pd.DataFrame,
    pairwise: pd.DataFrame,
    outdir: str | Path,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    outdir = Path(outdir)
    staged_dir = outdir / "staged_structures"
    staged_dir.mkdir(parents=True, exist_ok=True)
    members = members.copy()
    medoid_by_cluster = dict(zip(medoids["cluster_id"], medoids["medoid_id"]))
    tm_lookup = pair_lookup(pairwise, "symmetric_tm")
    max_members = int(config["medoid"]["max_within_members"])

    selected_ids: set[str] = set()
    for cluster_id, cdf in members.groupby("cluster_id", sort=True):
        eligible_ids = cdf.loc[cdf["eligible_for_medoid"], "protein_id"].astype(str).tolist()
        selected_ids.update(maximin_sample(eligible_ids, medoid_by_cluster[cluster_id], max_members, tm_lookup))
    members["selected_for_within_alignment"] = members["protein_id"].isin(selected_ids)
    members["is_medoid"] = members.apply(
        lambda r: str(r["protein_id"]) == str(medoid_by_cluster[r["cluster_id"]]), axis=1
    )

    staged_paths = []
    for row in members.itertuples(index=False):
        source = Path(str(row.structure_path))
        destination = staged_dir / f"{row.foldmason_id}{suffix_for_structure(source)}"
        if not destination.exists() or destination.stat().st_size != source.stat().st_size:
            shutil.copy2(source, destination)
        staged_paths.append(str(destination.resolve()))
    members["staged_structure_path"] = staged_paths

    jobs = []
    for cluster_id, cdf in members[members["selected_for_within_alignment"]].groupby("cluster_id", sort=True):
        base = outdir / "alignments" / "within" / sanitize_group(cluster_id)
        prefix = base / sanitize_group(cluster_id)
        jobs.append(
            {
                "scope": "within",
                "analysis_group": str(cdf["analysis_group"].iloc[0]),
                "cluster_id": str(cluster_id),
                "n_structures": len(cdf),
                "structure_paths_json": json_dumps(cdf["staged_structure_path"]),
                "foldmason_ids_json": json_dumps(cdf["foldmason_id"]),
                "output_prefix": str(prefix.resolve()),
                "tmp_dir": str((outdir / "tmp" / "within" / sanitize_group(cluster_id)).resolve()),
            }
        )
    medoid_members = members[members["is_medoid"]].copy()
    for group, gdf in medoid_members.groupby("analysis_group", sort=True):
        base = outdir / "alignments" / "between" / sanitize_group(group)
        prefix = base / sanitize_group(group)
        jobs.append(
            {
                "scope": "between",
                "analysis_group": str(group),
                "cluster_id": "",
                "n_structures": len(gdf),
                "structure_paths_json": json_dumps(gdf["staged_structure_path"]),
                "foldmason_ids_json": json_dumps(gdf["foldmason_id"]),
                "output_prefix": str(prefix.resolve()),
                "tmp_dir": str((outdir / "tmp" / "between" / sanitize_group(group)).resolve()),
            }
        )
    return members, pd.DataFrame(jobs)


def run_foldmason_jobs(jobs_path: str | Path, config_path: str | Path | None, dry_run: bool = False) -> None:
    config = load_config(config_path)
    fm = config["foldmason"]
    jobs = read_tsv(jobs_path)
    require_columns(jobs, ["scope", "structure_paths_json", "output_prefix", "tmp_dir"], "alignment_jobs.tsv")
    executable = str(fm["executable"])
    for row in jobs.itertuples(index=False):
        structures = json.loads(row.structure_paths_json)
        if len(structures) < 2:
            print(f"SKIP {row.scope} {getattr(row, 'cluster_id', '')}: fewer than two structures")
            continue
        Path(row.output_prefix).parent.mkdir(parents=True, exist_ok=True)
        Path(row.tmp_dir).mkdir(parents=True, exist_ok=True)
        command = [
            executable,
            "easy-msa",
            *structures,
            row.output_prefix,
            row.tmp_dir,
            "--report-mode",
            str(fm["report_mode"]),
            "--refine-iters",
            str(fm["refine_iters"]),
            "--pair-threshold",
            str(fm["pair_threshold"]),
            "--threads",
            str(fm["threads"]),
        ]
        print(" ".join(command))
        if not dry_run:
            subprocess.run(command, check=True)
