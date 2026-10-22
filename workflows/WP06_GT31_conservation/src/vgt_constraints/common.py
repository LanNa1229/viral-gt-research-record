from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
import yaml
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import is_aa


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}


DEFAULT_CONFIG = {
    "quality": {
        "min_model_plddt": 70.0,
        "min_model_residues": 100,
        "min_sequence_mapping_identity": 0.90,
        "min_sequence_mapping_coverage": 0.90,
        "min_local_plddt": 70.0,
    },
    "medoid": {
        "score_mode": "min_directional_tm",
        "min_pairwise_observations": 1,
        "max_within_members": 50,
    },
    "foldmason": {
        "executable": "foldmason",
        "report_mode": 1,
        "refine_iters": 0,
        "pair_threshold": 0.0,
        "threads": 8,
    },
    "lddt": {
        "inclusion_radius_A": 15.0,
        "thresholds_A": [0.5, 1.0, 2.0, 4.0],
        "min_neighbors": 5,
    },
    "between_scope": {
        "min_symmetric_tm": 0.40,
        "min_reciprocal_coverage": 0.40,
        "nearest_k": 20,
        "min_comparators_for_call": 3,
    },
    "classification": {
        "within_high": 0.75,
        "between_high": 0.65,
        "min_within_presence": 0.70,
        "min_within_lddt_coverage": 0.60,
        "min_between_presence": 0.50,
        "min_between_lddt_coverage": 0.40,
        "max_within_gap_frequency": 0.30,
    },
    "pan_core": {
        "min_subcluster_occupancy": 0.70,
        "min_local_plddt": 70.0,
        "within_high": 0.75,
        "cross_high": 0.65,
        "min_within_lddt_coverage": 0.60,
        "min_cross_lddt_coverage": 0.40,
        "min_cross_presence": 0.50,
        "supported_min_cluster_pass_fraction": 0.75,
        "region_min_length": 4,
        "region_max_bridge_residues": 1,
        "contact_cutoff_A": 10.0,
    },
    "scores": {
        "plddt_floor": 50.0,
        "plddt_full_weight": 90.0,
        "include_between_evidence_weight": True,
    },
    "regions": {
        "min_constraint_score": 0.12,
        "min_region_length": 4,
        "max_bridge_residues": 1,
        "min_candidate_fraction": 0.60,
        "active_site_window": 5,
        "top_alignment_blocks_per_cluster": 5,
    },
    "visualization": {
        "dpi": 300,
        "max_scatter_points_per_category": 20000,
    },
}


def deep_merge(base: dict, override: Mapping) -> dict:
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(out.get(key), Mapping):
            out[key] = deep_merge(dict(out[key]), value)
        else:
            out[key] = value
    return out


def load_config(path: str | Path | None) -> dict:
    if path is None:
        return DEFAULT_CONFIG
    with open(path) as handle:
        user = yaml.safe_load(handle) or {}
    return deep_merge(DEFAULT_CONFIG, user)


def require_columns(df: pd.DataFrame, required: Iterable[str], source: str) -> None:
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(f"{source} is missing required columns: {', '.join(missing)}")


def read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def safe_id(value: str) -> str:
    core = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_.") or "protein"
    digest = hashlib.sha1(str(value).encode()).hexdigest()[:8]
    return f"{core[:80]}__{digest}"


def sanitize_group(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_.") or "ALL"


def read_fasta(path: str | Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current = None
    chunks: list[str] = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                if current is None:
                    raise ValueError(f"Sequence encountered before FASTA header in {path}")
                chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks)
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Alignment sequences have unequal lengths in {path}: {sorted(lengths)}")
    return records


def normalize_alignment_id(header_id: str) -> str:
    name = Path(header_id).name
    for suffix in (".pdb.gz", ".cif.gz", ".mmcif.gz", ".pdb", ".cif", ".mmcif", ".gz"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def resolve_alignment_records(records: Mapping[str, str], expected_ids: Iterable[str]) -> dict[str, str]:
    expected = set(expected_ids)
    out: dict[str, str] = {}
    for raw_id, seq in records.items():
        candidates = [raw_id, normalize_alignment_id(raw_id)]
        hit = next((x for x in candidates if x in expected), None)
        if hit is not None:
            if hit in out:
                raise ValueError(f"Duplicate alignment record resolved to {hit}")
            out[hit] = seq
    missing = sorted(expected - set(out))
    if missing:
        raise ValueError(
            f"Alignment is missing {len(missing)} expected staged IDs; first missing: {missing[:5]}"
        )
    return out


@dataclass(frozen=True)
class ResidueRecord:
    seq_index: int
    chain_id: str
    resseq: int
    icode: str
    resname: str
    aa: str
    ca: np.ndarray
    plddt: float

    @property
    def residue_label(self) -> str:
        ins = self.icode.strip()
        return f"{self.chain_id}:{self.resseq}{ins}"


def structure_parser(path: Path):
    lower = path.name.lower()
    if lower.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def load_structure_residues(path: str | Path, chain_id: str | None = None) -> tuple[object, list[ResidueRecord]]:
    path = Path(path)
    structure = structure_parser(path).get_structure(path.stem, str(path))
    model = next(structure.get_models())
    residues: list[ResidueRecord] = []
    for chain in model:
        if chain_id and chain.id != chain_id:
            continue
        for residue in chain:
            if "CA" not in residue or not is_aa(residue, standard=False):
                continue
            resname = residue.get_resname().upper()
            aa = AA3_TO_1.get(resname, "X")
            ca_atom = residue["CA"]
            hetflag, resseq, icode = residue.id
            residues.append(
                ResidueRecord(
                    seq_index=len(residues) + 1,
                    chain_id=str(chain.id),
                    resseq=int(resseq),
                    icode=str(icode or ""),
                    resname=resname,
                    aa=aa,
                    ca=np.asarray(ca_atom.coord, dtype=float),
                    plddt=float(ca_atom.bfactor),
                )
            )
    if not residues:
        raise ValueError(f"No amino-acid CA atoms found in {path}")
    return structure, residues


def structure_qc(path: str | Path, chain_id: str | None = None) -> dict:
    _, residues = load_structure_residues(path, chain_id)
    chains = sorted({r.chain_id for r in residues})
    plddt = np.asarray([r.plddt for r in residues], dtype=float)
    return {
        "n_residues": len(residues),
        "mean_structure_bfactor": float(np.nanmean(plddt)),
        "median_structure_bfactor": float(np.nanmedian(plddt)),
        "n_chains": len(chains),
        "chains": ",".join(chains),
    }


def map_alignment_sequence_to_structure(
    aligned_sequence: str,
    residues: list[ResidueRecord],
    min_identity: float,
    min_coverage: float,
) -> tuple[dict[int, ResidueRecord], dict]:
    ungapped = "".join(x for x in aligned_sequence if x not in "-.").upper()
    structure_sequence = "".join(r.aa for r in residues)
    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(structure_sequence, ungapped)[0]
    ungapped_to_structure: dict[int, int] = {}
    matches = 0
    mapped = 0
    for (s0, s1), (u0, u1) in zip(alignment.aligned[0], alignment.aligned[1]):
        block = min(s1 - s0, u1 - u0)
        for offset in range(block):
            si = int(s0 + offset)
            ui = int(u0 + offset)
            ungapped_to_structure[ui] = si
            mapped += 1
            if structure_sequence[si] == ungapped[ui] or "X" in (structure_sequence[si], ungapped[ui]):
                matches += 1
    identity = matches / mapped if mapped else 0.0
    coverage = mapped / len(ungapped) if ungapped else 0.0
    if identity < min_identity or coverage < min_coverage:
        raise ValueError(
            f"MSA-to-structure mapping failed QC: identity={identity:.3f}, coverage={coverage:.3f}; "
            f"required {min_identity:.3f}/{min_coverage:.3f}"
        )
    column_map: dict[int, ResidueRecord] = {}
    u_index = 0
    for column, char in enumerate(aligned_sequence):
        if char in "-.":
            continue
        if u_index in ungapped_to_structure:
            column_map[column] = residues[ungapped_to_structure[u_index]]
        u_index += 1
    return column_map, {
        "mapping_identity": identity,
        "mapping_coverage": coverage,
        "msa_ungapped_length": len(ungapped),
        "structure_length": len(residues),
    }


def normalized_entropy_conservation(chars: Iterable[str], alphabet_size: int = 20) -> dict:
    values = [str(x) for x in chars if str(x) not in ("-", ".", "", "nan")]
    if not values:
        return {"conservation": np.nan, "entropy": np.nan, "dominant": "", "dominant_fraction": np.nan, "n": 0}
    counts = Counter(values)
    probs = np.asarray(list(counts.values()), dtype=float) / len(values)
    entropy = float(-(probs * np.log(probs)).sum())
    max_entropy = math.log(alphabet_size)
    conservation = 1.0 - min(1.0, entropy / max_entropy)
    dominant, count = counts.most_common(1)[0]
    return {
        "conservation": conservation,
        "entropy": entropy,
        "dominant": dominant,
        "dominant_fraction": count / len(values),
        "n": len(values),
    }


def json_dumps(values: Iterable[str]) -> str:
    return json.dumps(list(values), separators=(",", ":"))


class ResidueMetricSelect(Select):
    def __init__(self, values: Mapping[tuple[str, int, str], float], scale: float = 100.0):
        self.values = values
        self.scale = scale

    def accept_atom(self, atom):
        residue = atom.get_parent()
        chain = residue.get_parent()
        _, resseq, icode = residue.id
        value = self.values.get((str(chain.id), int(resseq), str(icode or "")), np.nan)
        atom.set_bfactor(0.0 if not np.isfinite(value) else float(np.clip(value * self.scale, 0.0, 100.0)))
        return 1


def write_metric_pdb(
    structure_path: str | Path,
    output_path: str | Path,
    metric_values: Mapping[tuple[str, int, str], float],
    scale: float = 100.0,
) -> None:
    structure = structure_parser(Path(structure_path)).get_structure("mapped", str(structure_path))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path), ResidueMetricSelect(metric_values, scale=scale))
