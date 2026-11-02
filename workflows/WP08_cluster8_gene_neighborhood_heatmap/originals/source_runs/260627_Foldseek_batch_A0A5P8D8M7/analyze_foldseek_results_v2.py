#!/usr/bin/env python3
"""
Analyze unified TSV output from foldseek_batch_api_v2.py.

Main input:
    final/foldseek_all_hits.tsv

Purpose:
    This script is optimized for viral GT / CAZy-family discovery workflows.
    It does not only show top Foldseek hits; it also extracts readable hit names,
    searches hit names/descriptions for biologically relevant terms, and writes
    database-level, query-level, and cluster-level summary tables.

Key outputs:
    00_all_hits_annotated.tsv
    01_top_hits_per_query_database_named.tsv
    02_high_confidence_hits.tsv
    03_keyword_hits.tsv
    04_pdb_priority_hits.tsv
    05_query_database_summary.tsv
    06_query_overall_summary.tsv
    07_cluster_database_summary.tsv
    08_pairwise_top_target_overlap.tsv
    09_query_cluster_coherence.tsv
    10_database_top_hit_names.tsv
    11_name_search_hits.tsv
    12_name_search_summary_by_query.tsv
    keyword_<category>_hits.tsv

Notes:
    - The script searches locally within Foldseek output fields such as target,
      target_id, and target_description. It does not call UniProt/PDB online.
    - If your Foldseek output lacks detailed descriptions, use PDB/UniProt target
      IDs from the output to perform a separate external annotation lookup.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

try:
    import pandas as pd
except ImportError as exc:
    raise SystemExit("Missing dependency: pandas. Install with: pip install pandas") from exc


DEFAULT_KEYWORDS: Dict[str, List[str]] = {
    "gt_like": [
        "glycosyltransferase", "glycosyl transferase", "glucosyltransferase",
        "glucosyl transferase", "galactosyltransferase", "mannosyltransferase",
        "xylosyltransferase", "rhamnosyltransferase", "fucosyltransferase",
        "sialyltransferase", "galactofuranosyltransferase", "glycosylation",
        "transferase family", "gt-", "gt ", "cazy", "udp-glucose",
        "udp glucose", "udp-glc", "udp-gal", "udp-glcna", "gdp-mannose",
        "glycosyl", "glucosyl", "gtr", "gt-a", "gt-b",
        "alpha-glucosyltransferase", "beta-glucosyltransferase",
        "dna alpha-glucosyltransferase", "dna beta-glucosyltransferase",
    ],
    "tet_like": [
        "tet", "t-even", "t4", "t2", "t6", "hydroxymethyl", "5-hydroxymethyl",
        "5hm", "5hmc", "5hmd", "5hmdU", "hmdu", "hydroxymethylcytosine",
        "hydroxymethyldeoxyuridine", "jbp", "base j", "ten-eleven",
    ],
    "dna_modification": [
        "dna", "methyltransferase", "restriction", "modification", "modified base",
        "base modification", "cytosine", "thymidine", "deoxyuridine", "nucleotide",
        "phage dna", "genome modification", "epigenetic", "oxygenase", "dioxygenase",
        "2-oxoglutarate", "2og", "fe(ii)", "iron-dependent",
    ],
    "lps_capsule_like": [
        "lipopolysaccharide", "lps", "o-antigen", "o antigen", "capsular",
        "capsule", "exopolysaccharide", "polysaccharide", "cell wall", "glycan",
        "peptidoglycan", "teichoic", "capsular polysaccharide",
    ],
    "enzyme_characterized_hint": [
        "crystal structure", "complex with", "bound to", "udp", "gdp", "acceptor",
        "donor", "catalytic", "mutant", "enzyme", "transferase", "substrate",
        "product", "inhibitor", "ligand",
    ],
    "hypothetical": [
        "hypothetical", "uncharacterized", "unknown function", "putative",
        "predicted protein", "domain of unknown function", "duf",
    ],
}

# These are shown as exact matched words/phrases in a separate name-search table.
DEFAULT_NAME_SEARCH_TERMS = [
    "glycosyltransferase", "glucosyltransferase", "glycosyl", "glucosyl",
    "galactosyltransferase", "mannosyltransferase", "xylosyltransferase",
    "rhamnosyltransferase", "fucosyltransferase", "CAZy", "GT-A", "GT-B",
    "TET", "T-even", "T4", "hydroxymethyl", "5hmC", "5hmdU", "JBP",
    "DNA", "methyltransferase", "oxygenase", "dioxygenase",
    "lipopolysaccharide", "LPS", "O-antigen", "capsular", "polysaccharide",
    "hypothetical", "uncharacterized", "putative",
]

NUMERIC_COLUMNS = [
    "rank", "seq_identity", "aln_length", "mismatches", "gaps_opened",
    "q_start", "q_end", "t_start", "t_end", "probability", "evalue",
    "score", "q_len", "t_len", "q_coverage", "t_coverage",
]

PREFERRED_COLUMNS = [
    "query_id", "protein_id", "cluster", "database", "rank", "target_id",
    "hit_name", "hit_annotation", "target_description", "matched_name_terms",
    "matched_keyword_categories", "probability", "evalue", "score", "seq_identity",
    "aln_length", "q_coverage", "t_coverage", "q_len", "t_len",
    "q_start", "q_end", "t_start", "t_end", "high_confidence",
]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze foldseek_all_hits.tsv and search hit names/descriptions for biological terms."
    )
    parser.add_argument("--input", required=True, help="Path to foldseek_all_hits.tsv")
    parser.add_argument("--outdir", default=None, help="Output directory. Default: <input_dir>/analysis")
    parser.add_argument("--top-n", type=int, default=20, help="Top N hits per query/database for summaries")
    parser.add_argument("--highconf-prob", type=float, default=0.70, help="High-confidence probability cutoff")
    parser.add_argument("--highconf-evalue", type=float, default=1e-3, help="High-confidence E-value cutoff")
    parser.add_argument("--highconf-qcov", type=float, default=0.50, help="High-confidence query coverage cutoff")
    parser.add_argument("--highconf-alnlen", type=int, default=80, help="High-confidence alignment length cutoff")
    parser.add_argument("--target-databases", nargs="*", default=["pdb100"],
                        help="Databases to treat as experimentally/curated high-priority. Default: pdb100")
    parser.add_argument("--keyword-file", default=None,
                        help="Optional TSV: category<TAB>keyword. Adds/overrides keyword categories if --replace-keywords is used.")
    parser.add_argument("--replace-keywords", action="store_true",
                        help="Use only --keyword-file keywords instead of defaults.")
    parser.add_argument("--name-terms", nargs="*", default=None,
                        help="Extra terms to search in hit names/descriptions, e.g. --name-terms glycosyltransferase TET 5hmC")
    parser.add_argument("--name-terms-file", default=None,
                        help="Optional one-term-per-line file for hit-name search terms.")
    parser.add_argument("--replace-name-terms", action="store_true",
                        help="Use only terms from --name-terms/--name-terms-file instead of defaults.")
    parser.add_argument("--outlier-top-k", type=int, default=20,
                        help="Top K target IDs per query/database used for target-overlap coherence analysis")
    return parser.parse_args(argv)


def read_keyword_file(path: Optional[str], replace_defaults: bool) -> Dict[str, List[str]]:
    keywords = {} if replace_defaults else {k: list(v) for k, v in DEFAULT_KEYWORDS.items()}
    if not path:
        return keywords
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Keyword file not found: {p}")
    with p.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"Keyword file line {line_no} must be category<TAB>keyword: {line}")
            category, keyword = parts[0].strip(), parts[1].strip()
            if category and keyword:
                keywords.setdefault(category, []).append(keyword)
    return keywords


def read_name_terms(args: argparse.Namespace) -> List[str]:
    terms: List[str] = [] if args.replace_name_terms else list(DEFAULT_NAME_SEARCH_TERMS)
    if args.name_terms:
        terms.extend(args.name_terms)
    if args.name_terms_file:
        p = Path(args.name_terms_file)
        if not p.exists():
            raise FileNotFoundError(f"Name terms file not found: {p}")
        with p.open() as handle:
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#"):
                    terms.append(line)
    # preserve order while removing duplicates case-insensitively
    seen = set()
    unique = []
    for t in terms:
        key = t.lower()
        if key not in seen:
            seen.add(key)
            unique.append(t)
    return unique


def ensure_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Input TSV is missing required columns: {missing}. Found: {list(df.columns)}")


def first_existing_text(row: pd.Series, cols: Sequence[str]) -> str:
    for c in cols:
        if c in row.index:
            val = row.get(c, "")
            if pd.notna(val) and str(val).strip():
                return str(val).strip()
    return ""


def clean_text(x: str) -> str:
    x = str(x or "")
    x = re.sub(r"\s+", " ", x).strip()
    return x


def derive_hit_name(row: pd.Series) -> str:
    """Create a concise display name from Foldseek target fields."""
    desc = clean_text(first_existing_text(row, ["target_description", "description", "target"] ))
    tid = clean_text(first_existing_text(row, ["target_id", "target"]))

    # Foldseek/PDB descriptions are often already useful. Trim very long strings.
    if desc:
        # Remove duplicated target id at the beginning if present.
        if tid and desc.lower().startswith(tid.lower()):
            desc = desc[len(tid):].strip(" |;:-_\t")
        # Cut common pipe-formatted records to readable part if needed.
        desc = re.sub(r"^\S+\|\S+\|", "", desc)
        return desc[:300]

    return tid[:300]


def derive_hit_annotation(row: pd.Series) -> str:
    """A longer annotation text searched for terms."""
    pieces = []
    for c in ["target", "target_id", "target_description", "hit_name"]:
        if c in row.index:
            v = clean_text(row.get(c, ""))
            if v and v not in pieces:
                pieces.append(v)
    return " | ".join(pieces)


def load_hits(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input TSV not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    ensure_columns(df, ["query_id", "database", "rank", "target", "target_id", "target_description"])

    for col in NUMERIC_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col].replace("", pd.NA), errors="coerce")

    for col in ["protein_id", "cluster", "source_genome", "model_source"]:
        if col not in df.columns:
            df[col] = ""

    # readable names / annotation text
    df["hit_name"] = df.apply(derive_hit_name, axis=1)
    df["hit_annotation"] = df.apply(derive_hit_annotation, axis=1)
    df["target_text"] = df["hit_annotation"].fillna("").astype(str).str.replace(r"\s+", " ", regex=True)

    df = df.sort_values(["query_id", "database", "rank"], kind="mergesort").reset_index(drop=True)
    return df


def compile_keyword_regex(words: List[str]) -> re.Pattern:
    escaped = []
    for w in words:
        w = w.strip()
        if not w:
            continue
        # Let user regex through only if wrapped in /.../
        if len(w) >= 2 and w.startswith("/") and w.endswith("/"):
            escaped.append(f"(?:{w[1:-1]})")
        else:
            # Word-ish boundary for short terms reduces false positives: TET should not match "tetramer" unless wanted.
            if re.fullmatch(r"[A-Za-z0-9_+-]{2,8}", w):
                escaped.append(rf"(?<![A-Za-z0-9_]){re.escape(w)}(?![A-Za-z0-9_])")
            else:
                escaped.append(re.escape(w))
    if not escaped:
        return re.compile(r"a^", re.IGNORECASE)
    return re.compile("|".join(escaped), re.IGNORECASE)


def add_keyword_flags(df: pd.DataFrame, keywords: Dict[str, List[str]]) -> pd.DataFrame:
    out = df.copy()
    for category, words in keywords.items():
        pattern = compile_keyword_regex(words)
        out[f"flag_{category}"] = out["target_text"].fillna("").map(lambda x: bool(pattern.search(str(x))))
    flag_cols = [c for c in out.columns if c.startswith("flag_")]
    out["flag_any_keyword"] = out[flag_cols].any(axis=1) if flag_cols else False
    out["matched_keyword_categories"] = out.apply(
        lambda row: ",".join([c.replace("flag_", "") for c in flag_cols if bool(row.get(c, False))]),
        axis=1,
    )
    return out


def add_name_search_terms(df: pd.DataFrame, terms: List[str]) -> pd.DataFrame:
    out = df.copy()
    term_patterns = [(t, compile_keyword_regex([t])) for t in terms]

    def match_terms(text: str) -> str:
        found = []
        for term, pattern in term_patterns:
            if pattern.search(str(text)):
                found.append(term)
        return ",".join(found)

    out["matched_name_terms"] = out["target_text"].fillna("").map(match_terms)
    out["flag_name_search_hit"] = out["matched_name_terms"].astype(str).ne("")
    return out


def add_high_conf_flag(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    out = df.copy()
    prob_ok = out["probability"].ge(args.highconf_prob) if "probability" in out else False
    e_ok = out["evalue"].le(args.highconf_evalue) if "evalue" in out else False
    qcov_ok = out["q_coverage"].ge(args.highconf_qcov) if "q_coverage" in out else False
    aln_ok = out["aln_length"].ge(args.highconf_alnlen) if "aln_length" in out else False
    out["high_confidence"] = prob_ok & e_ok & qcov_ok & aln_ok
    return out


def select_columns(df: pd.DataFrame, extra: Optional[List[str]] = None) -> pd.DataFrame:
    extra = extra or []
    cols = [c for c in PREFERRED_COLUMNS + extra if c in df.columns]
    rest = [c for c in df.columns if c not in cols and c != "target_text"]
    return df[cols + rest]


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def first_or_blank(series: pd.Series) -> str:
    for x in series:
        if pd.notna(x) and str(x) != "":
            return str(x)
    return ""


def build_top_hits(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    return (
        df.sort_values(["query_id", "database", "rank"], kind="mergesort")
          .groupby(["query_id", "database"], group_keys=False, dropna=False)
          .head(top_n)
          .reset_index(drop=True)
    )


def build_query_database_summary(df: pd.DataFrame) -> pd.DataFrame:
    flag_cols = [c for c in df.columns if c.startswith("flag_") and c not in {"flag_any_keyword", "flag_name_search_hit"}]
    rows = []
    for (query_id, database), g in df.groupby(["query_id", "database"], dropna=False, sort=True):
        g_sorted = g.sort_values("rank", kind="mergesort")
        top = g_sorted.head(1)
        name_hits = g[g["flag_name_search_hit"]].sort_values("rank", kind="mergesort")
        keyword_hits = g[g["flag_any_keyword"]].sort_values("rank", kind="mergesort")
        row = {
            "query_id": query_id,
            "protein_id": first_or_blank(g.get("protein_id", pd.Series(dtype=str))),
            "cluster": first_or_blank(g.get("cluster", pd.Series(dtype=str))),
            "database": database,
            "n_hits": len(g),
            "n_high_confidence": int(g.get("high_confidence", pd.Series(False, index=g.index)).sum()),
            "n_keyword_hits": int(g.get("flag_any_keyword", pd.Series(False, index=g.index)).sum()),
            "n_name_search_hits": int(g.get("flag_name_search_hit", pd.Series(False, index=g.index)).sum()),
            "best_target_id": first_or_blank(top.get("target_id", pd.Series(dtype=str))),
            "best_hit_name": first_or_blank(top.get("hit_name", pd.Series(dtype=str))),
            "best_target_description": first_or_blank(top.get("target_description", pd.Series(dtype=str))),
            "best_probability": top["probability"].iloc[0] if "probability" in top else pd.NA,
            "best_evalue": top["evalue"].iloc[0] if "evalue" in top else pd.NA,
            "best_score": top["score"].iloc[0] if "score" in top else pd.NA,
            "best_q_coverage": top["q_coverage"].iloc[0] if "q_coverage" in top else pd.NA,
            "best_aln_length": top["aln_length"].iloc[0] if "aln_length" in top else pd.NA,
            "best_name_search_rank": name_hits["rank"].iloc[0] if len(name_hits) else pd.NA,
            "best_name_search_target_id": first_or_blank(name_hits.get("target_id", pd.Series(dtype=str))),
            "best_name_search_hit_name": first_or_blank(name_hits.get("hit_name", pd.Series(dtype=str))),
            "best_name_search_terms": first_or_blank(name_hits.get("matched_name_terms", pd.Series(dtype=str))),
            "best_keyword_rank": keyword_hits["rank"].iloc[0] if len(keyword_hits) else pd.NA,
            "best_keyword_target_id": first_or_blank(keyword_hits.get("target_id", pd.Series(dtype=str))),
            "best_keyword_hit_name": first_or_blank(keyword_hits.get("hit_name", pd.Series(dtype=str))),
            "best_keyword_categories": first_or_blank(keyword_hits.get("matched_keyword_categories", pd.Series(dtype=str))),
        }
        for col in flag_cols:
            row[f"n_{col.replace('flag_', '')}"] = int(g[col].sum())
        rows.append(row)
    return pd.DataFrame(rows)


def build_query_overall_summary(df: pd.DataFrame, target_databases: Sequence[str]) -> pd.DataFrame:
    flag_cols = [c for c in df.columns if c.startswith("flag_") and c not in {"flag_any_keyword", "flag_name_search_hit"}]
    rows = []
    for query_id, g in df.groupby("query_id", dropna=False, sort=True):
        g_ranked = g.sort_values(["rank"], kind="mergesort")
        priority = g[g["database"].isin(target_databases)].sort_values(["database", "rank"], kind="mergesort")
        best_priority = priority.head(1)
        best_any = g_ranked.head(1)
        highconf = g[g.get("high_confidence", pd.Series(False, index=g.index))]
        keyword = g[g.get("flag_any_keyword", pd.Series(False, index=g.index))]
        name_hits = g[g.get("flag_name_search_hit", pd.Series(False, index=g.index))]

        row = {
            "query_id": query_id,
            "protein_id": first_or_blank(g.get("protein_id", pd.Series(dtype=str))),
            "cluster": first_or_blank(g.get("cluster", pd.Series(dtype=str))),
            "source_genome": first_or_blank(g.get("source_genome", pd.Series(dtype=str))),
            "model_source": first_or_blank(g.get("model_source", pd.Series(dtype=str))),
            "n_databases": g["database"].nunique(),
            "n_total_hits": len(g),
            "n_high_confidence_hits": len(highconf),
            "n_keyword_hits": len(keyword),
            "n_name_search_hits": len(name_hits),
            "best_any_database": first_or_blank(best_any.get("database", pd.Series(dtype=str))),
            "best_any_target_id": first_or_blank(best_any.get("target_id", pd.Series(dtype=str))),
            "best_any_hit_name": first_or_blank(best_any.get("hit_name", pd.Series(dtype=str))),
            "best_any_target_description": first_or_blank(best_any.get("target_description", pd.Series(dtype=str))),
            "best_priority_database": first_or_blank(best_priority.get("database", pd.Series(dtype=str))),
            "best_priority_target_id": first_or_blank(best_priority.get("target_id", pd.Series(dtype=str))),
            "best_priority_hit_name": first_or_blank(best_priority.get("hit_name", pd.Series(dtype=str))),
            "best_priority_target_description": first_or_blank(best_priority.get("target_description", pd.Series(dtype=str))),
        }
        for col in flag_cols:
            category = col.replace("flag_", "")
            hits = g[g[col]].sort_values(["database", "rank"], kind="mergesort")
            row[f"has_{category}"] = bool(len(hits))
            if len(hits):
                htop = hits.head(1)
                row[f"best_{category}_database"] = first_or_blank(htop.get("database", pd.Series(dtype=str)))
                row[f"best_{category}_rank"] = htop["rank"].iloc[0] if "rank" in htop else pd.NA
                row[f"best_{category}_target_id"] = first_or_blank(htop.get("target_id", pd.Series(dtype=str)))
                row[f"best_{category}_hit_name"] = first_or_blank(htop.get("hit_name", pd.Series(dtype=str)))
                row[f"best_{category}_description"] = first_or_blank(htop.get("target_description", pd.Series(dtype=str)))
            else:
                row[f"best_{category}_database"] = ""
                row[f"best_{category}_rank"] = pd.NA
                row[f"best_{category}_target_id"] = ""
                row[f"best_{category}_hit_name"] = ""
                row[f"best_{category}_description"] = ""
        rows.append(row)
    return pd.DataFrame(rows)


def build_cluster_summary(df: pd.DataFrame) -> pd.DataFrame:
    flag_cols = [c for c in df.columns if c.startswith("flag_") and c not in {"flag_any_keyword", "flag_name_search_hit"}]
    rows = []
    for (cluster, database), g in df.groupby(["cluster", "database"], dropna=False, sort=True):
        row = {
            "cluster": cluster if cluster else "NA",
            "database": database,
            "n_queries": g["query_id"].nunique(),
            "n_hits": len(g),
            "n_high_confidence_hits": int(g.get("high_confidence", pd.Series(False, index=g.index)).sum()),
            "n_keyword_hits": int(g.get("flag_any_keyword", pd.Series(False, index=g.index)).sum()),
            "n_name_search_hits": int(g.get("flag_name_search_hit", pd.Series(False, index=g.index)).sum()),
            "queries_with_high_confidence": int(g[g.get("high_confidence", pd.Series(False, index=g.index))]["query_id"].nunique()),
            "queries_with_any_keyword": int(g[g.get("flag_any_keyword", pd.Series(False, index=g.index))]["query_id"].nunique()),
            "queries_with_name_search_hit": int(g[g.get("flag_name_search_hit", pd.Series(False, index=g.index))]["query_id"].nunique()),
        }
        for col in flag_cols:
            category = col.replace("flag_", "")
            hits = g[g[col]]
            row[f"hits_{category}"] = int(len(hits))
            row[f"queries_with_{category}"] = int(hits["query_id"].nunique())
        rows.append(row)
    return pd.DataFrame(rows)


def build_name_search_summary_by_query(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    hits = df[df["flag_name_search_hit"]].copy()
    for (query_id, database), g in df.groupby(["query_id", "database"], dropna=False, sort=True):
        gh = hits[(hits["query_id"].astype(str) == str(query_id)) & (hits["database"].astype(str) == str(database))]
        gh = gh.sort_values("rank", kind="mergesort")
        top = gh.head(1)
        terms = sorted(set(",".join(gh["matched_name_terms"].dropna().astype(str)).split(",")) - {""})
        rows.append({
            "query_id": query_id,
            "protein_id": first_or_blank(g.get("protein_id", pd.Series(dtype=str))),
            "cluster": first_or_blank(g.get("cluster", pd.Series(dtype=str))),
            "database": database,
            "n_name_search_hits": len(gh),
            "all_matched_name_terms": ",".join(terms),
            "best_name_search_rank": top["rank"].iloc[0] if len(top) else pd.NA,
            "best_name_search_target_id": first_or_blank(top.get("target_id", pd.Series(dtype=str))),
            "best_name_search_hit_name": first_or_blank(top.get("hit_name", pd.Series(dtype=str))),
            "best_name_search_matched_terms": first_or_blank(top.get("matched_name_terms", pd.Series(dtype=str))),
            "best_name_search_probability": top["probability"].iloc[0] if len(top) and "probability" in top else pd.NA,
            "best_name_search_evalue": top["evalue"].iloc[0] if len(top) and "evalue" in top else pd.NA,
            "best_name_search_q_coverage": top["q_coverage"].iloc[0] if len(top) and "q_coverage" in top else pd.NA,
            "best_name_search_aln_length": top["aln_length"].iloc[0] if len(top) and "aln_length" in top else pd.NA,
        })
    return pd.DataFrame(rows)


def build_database_top_hit_names(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    rows = []
    for (query_id, database), g in df.groupby(["query_id", "database"], dropna=False, sort=True):
        qtop = g.sort_values("rank", kind="mergesort").head(top_n)
        for _, row in qtop.iterrows():
            rows.append({
                "query_id": row.get("query_id", ""),
                "protein_id": row.get("protein_id", ""),
                "cluster": row.get("cluster", ""),
                "database": row.get("database", ""),
                "rank": row.get("rank", pd.NA),
                "target_id": row.get("target_id", ""),
                "hit_name": row.get("hit_name", ""),
                "matched_name_terms": row.get("matched_name_terms", ""),
                "matched_keyword_categories": row.get("matched_keyword_categories", ""),
                "probability": row.get("probability", pd.NA),
                "evalue": row.get("evalue", pd.NA),
                "score": row.get("score", pd.NA),
                "seq_identity": row.get("seq_identity", pd.NA),
                "aln_length": row.get("aln_length", pd.NA),
                "q_coverage": row.get("q_coverage", pd.NA),
                "t_coverage": row.get("t_coverage", pd.NA),
                "target_description": row.get("target_description", ""),
            })
    return pd.DataFrame(rows)


def build_top_overlap(df: pd.DataFrame, top_k: int) -> pd.DataFrame:
    rows = []
    for database, dg in df.groupby("database", dropna=False, sort=True):
        target_sets: Dict[str, set] = {}
        meta: Dict[str, Dict[str, str]] = {}
        for query_id, qg in dg.groupby("query_id", dropna=False, sort=True):
            qtop = qg.sort_values("rank", kind="mergesort").head(top_k)
            target_sets[str(query_id)] = set(qtop["target_id"].dropna().astype(str)) - {""}
            meta[str(query_id)] = {
                "cluster": first_or_blank(qg.get("cluster", pd.Series(dtype=str))),
                "protein_id": first_or_blank(qg.get("protein_id", pd.Series(dtype=str))),
            }
        qids = sorted(target_sets)
        for i, q1 in enumerate(qids):
            for q2 in qids[i + 1:]:
                s1, s2 = target_sets[q1], target_sets[q2]
                union = s1 | s2
                inter = s1 & s2
                rows.append({
                    "database": database,
                    "query_id_1": q1,
                    "protein_id_1": meta[q1]["protein_id"],
                    "cluster_1": meta[q1]["cluster"],
                    "query_id_2": q2,
                    "protein_id_2": meta[q2]["protein_id"],
                    "cluster_2": meta[q2]["cluster"],
                    "top_k": top_k,
                    "n_overlap": len(inter),
                    "n_union": len(union),
                    "jaccard_top_targets": len(inter) / len(union) if union else math.nan,
                    "overlapping_target_ids": ",".join(sorted(inter)),
                })
    return pd.DataFrame(rows)


def build_query_coherence(df: pd.DataFrame, top_k: int) -> pd.DataFrame:
    rows = []
    for (cluster, database), g in df.groupby(["cluster", "database"], dropna=False, sort=True):
        target_sets = {}
        for query_id, qg in g.groupby("query_id", dropna=False, sort=True):
            qtop = qg.sort_values("rank", kind="mergesort").head(top_k)
            target_sets[str(query_id)] = set(qtop["target_id"].dropna().astype(str)) - {""}
        qids = sorted(target_sets)
        if len(qids) < 2:
            continue
        for qid in qids:
            scores = []
            for other in qids:
                if other == qid:
                    continue
                s1, s2 = target_sets[qid], target_sets[other]
                union = s1 | s2
                scores.append(len(s1 & s2) / len(union) if union else math.nan)
            valid = [s for s in scores if not pd.isna(s)]
            rows.append({
                "cluster": cluster if cluster else "NA",
                "database": database,
                "query_id": qid,
                "top_k": top_k,
                "n_cluster_queries": len(qids),
                "mean_jaccard_to_cluster_mates": sum(valid) / len(valid) if valid else math.nan,
                "min_jaccard_to_cluster_mates": min(valid) if valid else math.nan,
                "max_jaccard_to_cluster_mates": max(valid) if valid else math.nan,
            })
    return pd.DataFrame(rows)


def infer_pdb_code(target_id: str) -> str:
    if not isinstance(target_id, str):
        return ""
    m = re.search(r"\b([0-9][A-Za-z0-9]{3})(?:[_\-.:][A-Za-z0-9])?\b", target_id)
    return m.group(1).upper() if m else ""


def add_pdb_code(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["pdb_code_inferred"] = out["target_id"].fillna("").astype(str).map(infer_pdb_code)
    return out


def write_markdown_report(
    path: Path,
    df: pd.DataFrame,
    highconf: pd.DataFrame,
    keyword_hits: pd.DataFrame,
    name_hits: pd.DataFrame,
    priority_hits: pd.DataFrame,
    args: argparse.Namespace,
    name_terms: Sequence[str],
) -> None:
    n_queries = df["query_id"].nunique()
    n_clusters = df["cluster"].replace("", pd.NA).dropna().nunique() if "cluster" in df else 0
    dbs = ", ".join(map(str, sorted(df["database"].dropna().unique())))
    flag_cols = [c for c in df.columns if c.startswith("flag_") and c not in {"flag_any_keyword", "flag_name_search_hit"}]

    lines = []
    lines.append("# Foldseek batch result analysis")
    lines.append("")
    lines.append("## Input overview")
    lines.append(f"- Total parsed hit rows: **{len(df)}**")
    lines.append(f"- Unique queries: **{n_queries}**")
    lines.append(f"- Unique clusters: **{n_clusters}**")
    lines.append(f"- Databases: **{dbs}**")
    lines.append("")
    lines.append("## High-confidence rule")
    lines.append(f"- probability >= {args.highconf_prob}")
    lines.append(f"- evalue <= {args.highconf_evalue:g}")
    lines.append(f"- q_coverage >= {args.highconf_qcov}")
    lines.append(f"- aln_length >= {args.highconf_alnlen}")
    lines.append(f"- High-confidence hits: **{len(highconf)}**")
    lines.append("")
    lines.append("## Hit-name search")
    lines.append(f"- Name search terms: {', '.join(name_terms[:40])}{' ...' if len(name_terms) > 40 else ''}")
    lines.append(f"- Name-search hits: **{len(name_hits)}**")
    lines.append(f"- Queries with name-search hits: **{name_hits['query_id'].nunique() if len(name_hits) else 0}**")
    lines.append("")
    lines.append("## Keyword category hit counts")
    for col in flag_cols:
        lines.append(f"- {col.replace('flag_', '')}: {int(df[col].sum())} hits, {df[df[col]]['query_id'].nunique()} queries")
    lines.append(f"- any_keyword: {len(keyword_hits)} hits, {keyword_hits['query_id'].nunique() if len(keyword_hits) else 0} queries")
    lines.append("")
    lines.append("## Priority database hits")
    lines.append(f"- Priority databases: {', '.join(args.target_databases)}")
    lines.append(f"- Priority-database rows: {len(priority_hits)}")
    lines.append("")
    lines.append("## Key output files")
    for name in [
        "00_all_hits_annotated.tsv",
        "01_top_hits_per_query_database_named.tsv",
        "02_high_confidence_hits.tsv",
        "03_keyword_hits.tsv",
        "04_pdb_priority_hits.tsv",
        "05_query_database_summary.tsv",
        "06_query_overall_summary.tsv",
        "07_cluster_database_summary.tsv",
        "08_pairwise_top_target_overlap.tsv",
        "09_query_cluster_coherence.tsv",
        "10_database_top_hit_names.tsv",
        "11_name_search_hits.tsv",
        "12_name_search_summary_by_query.tsv",
    ]:
        lines.append(f"- `{name}`")
    lines.append("")
    lines.append("## Interpretation notes")
    lines.append("- Start with `10_database_top_hit_names.tsv` for quick readable top-hit inspection.")
    lines.append("- Use `11_name_search_hits.tsv` to find hits whose names/descriptions contain terms such as glycosyltransferase, TET, 5hmC, DNA, LPS, or capsular polysaccharide.")
    lines.append("- Use `04_pdb_priority_hits.tsv` for experimentally grounded structural homologs.")
    lines.append("- Keyword/name hits are triage labels only. Final interpretation should include alignment coverage, active-site conservation, structural superposition, and gene-neighborhood evidence.")
    path.write_text("\n".join(lines))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    input_path = Path(args.input).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else input_path.parent / "analysis"
    outdir.mkdir(parents=True, exist_ok=True)

    keywords = read_keyword_file(args.keyword_file, args.replace_keywords)
    name_terms = read_name_terms(args)

    df = load_hits(input_path)
    df = add_keyword_flags(df, keywords)
    df = add_name_search_terms(df, name_terms)
    df = add_high_conf_flag(df, args)
    df = add_pdb_code(df)

    annotated_cols = [
        "high_confidence", "matched_keyword_categories", "matched_name_terms",
        "flag_name_search_hit", "pdb_code_inferred",
    ]

    # 00: annotated all hits
    write_tsv(select_columns(df, annotated_cols), outdir / "00_all_hits_annotated.tsv")

    # 01: top hits per query/database with names up front
    top_hits = build_top_hits(df, args.top_n)
    write_tsv(select_columns(top_hits, annotated_cols), outdir / "01_top_hits_per_query_database_named.tsv")

    # 02: high confidence hits
    highconf = df[df["high_confidence"]].copy()
    write_tsv(select_columns(highconf, annotated_cols), outdir / "02_high_confidence_hits.tsv")

    # 03: keyword hits
    keyword_hits = df[df["flag_any_keyword"]].copy()
    write_tsv(select_columns(keyword_hits, annotated_cols), outdir / "03_keyword_hits.tsv")

    # 04: priority/PDB hits
    priority_hits = df[df["database"].isin(args.target_databases)].copy()
    write_tsv(select_columns(priority_hits, annotated_cols), outdir / "04_pdb_priority_hits.tsv")

    # Per-category keyword tables.
    for col in [c for c in df.columns if c.startswith("flag_") and c not in {"flag_any_keyword", "flag_name_search_hit"}]:
        category = col.replace("flag_", "")
        sub = df[df[col]].copy()
        write_tsv(select_columns(sub, annotated_cols), outdir / f"keyword_{category}_hits.tsv")

    # 05-07: summaries
    qdb_summary = build_query_database_summary(df)
    write_tsv(qdb_summary, outdir / "05_query_database_summary.tsv")

    q_summary = build_query_overall_summary(df, args.target_databases)
    write_tsv(q_summary, outdir / "06_query_overall_summary.tsv")

    cluster_summary = build_cluster_summary(df)
    write_tsv(cluster_summary, outdir / "07_cluster_database_summary.tsv")

    # 08-09: coherence/outlier analyses
    overlap = build_top_overlap(df, args.outlier_top_k)
    write_tsv(overlap, outdir / "08_pairwise_top_target_overlap.tsv")

    coherence = build_query_coherence(df, args.outlier_top_k)
    write_tsv(coherence, outdir / "09_query_cluster_coherence.tsv")

    # 10: compact readable top-hit names
    db_top_names = build_database_top_hit_names(df, args.top_n)
    write_tsv(db_top_names, outdir / "10_database_top_hit_names.tsv")

    # 11: explicit name search hits
    name_hits = df[df["flag_name_search_hit"]].copy()
    write_tsv(select_columns(name_hits, annotated_cols), outdir / "11_name_search_hits.tsv")

    # 12: query/database-level name search summary
    name_summary = build_name_search_summary_by_query(df)
    write_tsv(name_summary, outdir / "12_name_search_summary_by_query.tsv")

    write_markdown_report(
        outdir / "foldseek_analysis_report.md",
        df=df,
        highconf=highconf,
        keyword_hits=keyword_hits,
        name_hits=name_hits,
        priority_hits=priority_hits,
        args=args,
        name_terms=name_terms,
    )

    print(f"Input rows: {len(df)}")
    print(f"Unique queries: {df['query_id'].nunique()}")
    print(f"Name-search hits: {len(name_hits)}")
    print(f"Output directory: {outdir}")
    print(f"Readable top-hit table: {outdir / '10_database_top_hit_names.tsv'}")
    print(f"Name-search hit table: {outdir / '11_name_search_hits.tsv'}")
    print(f"Main report: {outdir / 'foldseek_analysis_report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
