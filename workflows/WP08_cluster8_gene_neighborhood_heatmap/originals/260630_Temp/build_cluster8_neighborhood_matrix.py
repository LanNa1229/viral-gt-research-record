#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd


GENE_CLASSES = [
    # DNA base precursor / hydroxymethylation
    "C5_MT",
    "TET_5mYOX",
    "dCMP_hydroxymethylase",
    "hmC_hmdU_related",

    # pyrimidine / nucleotide metabolism
    "ThyX",
    "ThyS",
    "dCMP_dTMP_enzyme",
    "nucleotidyltransferase",

    # DNA metabolism
    "DNA_polymerase",
    "helicase",
    "ligase",
    "nuclease",
    "primase",
    "recombinase",

    # glycan pathway
    "cluster8_GT_like",
    "TET_associated_GT_like",
    "additional_GT",
    "NDP_sugar_enzyme",
    "epimerase",
    "dehydratase",
    "reductase",
    "aminotransferase",
    "methyltransferase",
    "acetyltransferase",
    "sulfotransferase",
    "APS_PAPS_kinase",

    # virion structural
    "capsid",
    "portal",
    "tail",
    "tape_measure",
    "terminase",

    # unknown / dark ORFs
    "conserved_hypothetical",
    "DUF",
]


def safe_float(x, default=0.0):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default


def normalize_text(*fields):
    text = " ".join(str(x) for x in fields if pd.notna(x))
    return text.lower()


def evidence_score(row):
    """
    Convert Foldseek evidence to 0/1/2/3.
    0 = absent
    1 = weak / keyword-level / low coverage
    2 = medium structural/domain evidence
    3 = strong structural evidence
    """
    db = str(row.get("database", "")).lower()
    prob = safe_float(row.get("probability", 0))
    evalue = safe_float(row.get("evalue", 999))
    qcov = safe_float(row.get("q_coverage", 0))
    tcov = safe_float(row.get("t_coverage", 0))
    aln = safe_float(row.get("aln_length", 0))
    rank = safe_float(row.get("rank", 999))

    # Strong: good Foldseek probability, coverage, and alignment
    if prob >= 0.90 and evalue <= 1e-5 and qcov >= 0.50 and aln >= 100:
        score = 3
    elif prob >= 0.70 and evalue <= 1e-3 and qcov >= 0.35 and aln >= 80:
        score = 2
    elif evalue <= 1e-2 and aln >= 50:
        score = 1
    else:
        score = 0

    # PDB/SwissProt are more interpretable; CATH fold-only hits should not overdrive function
    if "cath" in db and score > 2:
        score = 2

    # Very low-rank hits are less reliable for functional assignment
    if rank > 20 and score > 1:
        score -= 1

    return max(score, 0)


def classify_hit(row, cluster8_query_id=None):
    """
    Classify one Foldseek hit into one or more gene classes.
    Uses hit_name plus target_id and keyword categories.
    """
    query_id = str(row.get("query_id", ""))
    text = normalize_text(
        row.get("hit_name", ""),
        row.get("target_id", ""),
        row.get("matched_name_terms", ""),
        row.get("matched_keyword_categories", ""),
    )

    classes = set()

    is_cluster8_orf = cluster8_query_id is not None and query_id == str(cluster8_query_id)

    # -------------------------
    # DNA base precursor / DNA modification
    # -------------------------
    if re.search(r"\bc5\b.*methyltransferase|cytosine.*methyltransferase|dna.*cytosine.*methyltransferase|5-methylcytosine.*methyltransferase", text):
        classes.add("C5_MT")

    if re.search(r"5m yox|5myox|5-methylpyrimidine.*dioxygenase|tet-like|tet family|2og.*fe.*dioxygenase|2-oxoglutarate.*dioxygenase", text):
        # Avoid counting "TET-associated glycosyltransferase" as oxygenase.
        if not re.search(r"tet-associated.*glycosyltransferase|tet associated.*glycosyltransferase", text):
            classes.add("TET_5mYOX")

    if re.search(r"dcmp.*hydroxymethylase|hydroxymethyl.*dcmp|thymidylate synthase/dcmp hydroxymethylase|dutpase.*hydroxymethyl", text):
        classes.add("dCMP_hydroxymethylase")

    if re.search(r"5hmc|5-hmc|5hmd.?u|5-hmd.?u|hydroxymethylcytosine|hydroxymethyluracil|hydroxymethyl.*dna|hmdutp|hmdctp", text):
        classes.add("hmC_hmdU_related")

    # -------------------------
    # Pyrimidine / nucleotide metabolism
    # -------------------------
    if re.search(r"\bthyx\b|flavin-dependent thymidylate synthase", text):
        classes.add("ThyX")

    if re.search(r"\bthys\b|thymidylate synthase", text):
        # dCMP hydroxymethylase can be thymidylate-synthase-like; keep both if text supports it.
        classes.add("ThyS")

    if re.search(r"dcmp|dtmp|dump|dutp|dctp|dttp|pyrimidine|thymidine|deoxycytidine|deoxyuridine", text):
        classes.add("dCMP_dTMP_enzyme")

    if re.search(r"nucleotidyltransferase|pyrophosphorylase|nucleoside-diphosphate|nucleoside diphosphate|ntp transferase|sugar-1-phosphate nucleotidyltransferase", text):
        classes.add("nucleotidyltransferase")

    # -------------------------
    # DNA metabolism
    # -------------------------
    if re.search(r"dna polymerase|polymerase iii|polymerase i|family a polymerase|family b polymerase", text):
        classes.add("DNA_polymerase")

    if re.search(r"helicase|dnab|uvrd|reca-like helicase", text):
        classes.add("helicase")

    if re.search(r"ligase|dna ligase", text):
        classes.add("ligase")

    if re.search(r"nuclease|endonuclease|exonuclease|ribonuclease h|restriction endonuclease|hnh", text):
        classes.add("nuclease")

    if re.search(r"primase|dna primase", text):
        classes.add("primase")

    if re.search(r"recombinase|integrase|resolvase|rec[a-z]", text):
        classes.add("recombinase")

    # -------------------------
    # Glycan pathway
    # -------------------------
    tet_assoc_gt = re.search(
        r"tet-associated.*glycosyltransferase|tet associated.*glycosyltransferase|tet.*gt",
        text
    )

    gt_specific = re.search(
        r"glycosyltransferase|glycosyl transferase|glucosyltransferase|galactosyltransferase|"
        r"mannosyltransferase|rhamnosyltransferase|xylosyltransferase|fucosyltransferase|"
        r"hexosyltransferase|gt2\b|gt4\b|gt-a|gt-b|cazy.*gt",
        text
    )

    # Avoid generic transferases being called GT
    generic_not_gt = re.search(
        r"sulfotransferase|methyltransferase|phosphotransferase|aminotransferase|"
        r"acetyltransferase|kinase|ligase|polymerase|nucleotidyltransferase",
        text
    )

    if is_cluster8_orf and (tet_assoc_gt or gt_specific):
        classes.add("cluster8_GT_like")

    if tet_assoc_gt:
        classes.add("TET_associated_GT_like")

    if gt_specific and not generic_not_gt:
        if is_cluster8_orf:
            classes.add("cluster8_GT_like")
        else:
            classes.add("additional_GT")

    if re.search(r"ndp-sugar|nucleotide sugar|udp-sugar|gdp-sugar|cdp-sugar|tdp-sugar|"
                 r"udp-glucose|udp-glcnac|udp-galactose|gdp-mannose|dtdp|rml[a-d]|rfb[a-d]", text):
        classes.add("NDP_sugar_enzyme")

    if re.search(r"epimerase|dehydratase|dehydrogenase|isomerase|mutase", text):
        if re.search(r"sugar|glucose|galactose|mannose|rhamnose|fucose|nucleotide|ndp|udp|gdp|cdp|tdp|dtdp", text):
            classes.add("epimerase")

    if re.search(r"dehydratase", text):
        classes.add("dehydratase")

    if re.search(r"reductase", text):
        classes.add("reductase")

    if re.search(r"aminotransferase|transaminase", text):
        classes.add("aminotransferase")

    if re.search(r"methyltransferase", text):
        classes.add("methyltransferase")

    if re.search(r"acetyltransferase|acyltransferase", text):
        classes.add("acetyltransferase")

    if re.search(r"sulfotransferase", text):
        classes.add("sulfotransferase")

    if re.search(r"adenylyl-sulfate kinase|aps kinase|paps|sulfate adenylyltransferase|atps|adenosine.*phosphosulfate", text):
        classes.add("APS_PAPS_kinase")

    # -------------------------
    # Virion structural
    # -------------------------
    if re.search(r"capsid|major capsid|minor capsid|hk97", text):
        classes.add("capsid")

    if re.search(r"portal", text):
        classes.add("portal")

    if re.search(r"tail|baseplate|fiber", text):
        classes.add("tail")

    if re.search(r"tape measure|tapemeasure", text):
        classes.add("tape_measure")

    if re.search(r"terminase|large terminase|small terminase", text):
        classes.add("terminase")

    # -------------------------
    # Unknown / dark proteins
    # -------------------------
    if re.search(r"conserved hypothetical|hypothetical protein|uncharacterized protein", text):
        classes.add("conserved_hypothetical")

    if re.search(r"\bduf\d+|domain of unknown function", text):
        classes.add("DUF")

    return classes


def infer_scores_for_candidate(top_hits, candidate_id, cluster8_query_id=None):
    row = {c: 0 for c in GENE_CLASSES}
    evidence_records = []

    for _, hit in top_hits.iterrows():
        score = evidence_score(hit)
        if score == 0:
            continue

        classes = classify_hit(hit, cluster8_query_id=cluster8_query_id)
        if not classes:
            continue

        for cls in classes:
            if cls not in row:
                continue

            if score > row[cls]:
                row[cls] = score

            evidence_records.append({
                "candidate_id": candidate_id,
                "query_id": hit.get("query_id", ""),
                "protein_id": hit.get("protein_id", ""),
                "gene_class": cls,
                "score": score,
                "database": hit.get("database", ""),
                "rank": hit.get("rank", ""),
                "target_id": hit.get("target_id", ""),
                "hit_name": hit.get("hit_name", ""),
                "probability": hit.get("probability", ""),
                "evalue": hit.get("evalue", ""),
                "q_coverage": hit.get("q_coverage", ""),
                "t_coverage": hit.get("t_coverage", ""),
                "aln_length": hit.get("aln_length", ""),
            })

    return row, evidence_records


def compute_module_scores(row):
    dna_cols = [
        "C5_MT",
        "TET_5mYOX",
        "dCMP_hydroxymethylase",
        "hmC_hmdU_related",
        "ThyX",
        "ThyS",
        "dCMP_dTMP_enzyme",
        "DNA_polymerase",
        "helicase",
        "ligase",
        "nuclease",
        "primase",
        "recombinase",
    ]

    glycan_cols = [
        "cluster8_GT_like",
        "TET_associated_GT_like",
        "additional_GT",
        "NDP_sugar_enzyme",
        "epimerase",
        "dehydratase",
        "reductase",
        "aminotransferase",
        "methyltransferase",
        "acetyltransferase",
        "sulfotransferase",
        "APS_PAPS_kinase",
    ]

    virion_cols = [
        "capsid",
        "portal",
        "tail",
        "tape_measure",
        "terminase",
    ]

    dna_score = sum(row.get(c, 0) for c in dna_cols)
    glycan_score = sum(row.get(c, 0) for c in glycan_cols)
    virion_score = sum(row.get(c, 0) for c in virion_cols)

    return dna_score, glycan_score, virion_score


def infer_primary_group(row):
    """
    Row grouping for Panel B.
    Important: cluster8_GT_like itself should not be enough to call DNA context.
    """
    if row.get("C5_MT", 0) >= 2 and row.get("TET_5mYOX", 0) >= 2:
        return "C5-MT + 5mYOX-associated"

    if row.get("dCMP_hydroxymethylase", 0) >= 2:
        return "dCMP hydroxymethylase-associated"

    if (
        row.get("ThyX", 0) >= 2
        or row.get("ThyS", 0) >= 2
        or row.get("dCMP_dTMP_enzyme", 0) >= 2
        or row.get("nucleotidyltransferase", 0) >= 2
    ):
        return "ThyX/ThyS / nucleotide-metabolism-associated"

    # Here, use accessory glycan context, not cluster8 itself.
    accessory_glycan = max(
        row.get("additional_GT", 0),
        row.get("NDP_sugar_enzyme", 0),
        row.get("epimerase", 0),
        row.get("dehydratase", 0),
        row.get("reductase", 0),
        row.get("aminotransferase", 0),
        row.get("methyltransferase", 0),
        row.get("acetyltransferase", 0),
        row.get("sulfotransferase", 0),
        row.get("APS_PAPS_kinase", 0),
    )

    if accessory_glycan >= 2:
        return "additional-GT/glycan-associated"

    return "ambiguous/no DNA-context"


def infer_model(row):
    has_dna_precursor = (
        row.get("C5_MT", 0) >= 2
        or row.get("TET_5mYOX", 0) >= 2
        or row.get("dCMP_hydroxymethylase", 0) >= 2
        or row.get("hmC_hmdU_related", 0) >= 2
        or row.get("ThyX", 0) >= 2
        or row.get("ThyS", 0) >= 2
    )

    has_glycan = (
        row.get("cluster8_GT_like", 0) >= 2
        or row.get("TET_associated_GT_like", 0) >= 2
        or row.get("additional_GT", 0) >= 2
        or row.get("NDP_sugar_enzyme", 0) >= 2
    )

    has_virion = (
        row.get("capsid", 0) >= 2
        or row.get("tail", 0) >= 2
        or row.get("portal", 0) >= 2
        or row.get("tape_measure", 0) >= 2
    )

    if has_dna_precursor and has_glycan:
        return "hybrid DNA-modification + glycan module"
    if has_dna_precursor:
        return "DNA-modification-associated"
    if has_glycan and has_virion:
        return "possible virion/capsid glycosylation context"
    if has_glycan:
        return "glycan-associated but no DNA context"
    return "ambiguous"


def main():
    parser = argparse.ArgumentParser(
        description="Build cluster8_neighborhood_matrix.tsv from multiple analyze_foldseek_results_v2.py output folders."
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="TSV with columns: candidate_id, cluster8_query_id, analysis_dir"
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--top-hit-file",
        default="10_database_top_hit_names.tsv",
        help="Foldseek analysis table to read from each analysis_dir"
    )
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    matrix_rows = []
    all_evidence = []

    required_cols = {"candidate_id", "analysis_dir"}
    missing = required_cols - set(manifest.columns)
    if missing:
        raise ValueError(f"Manifest missing required columns: {missing}")

    for _, rec in manifest.iterrows():
        candidate_id = str(rec["candidate_id"])
        analysis_dir = Path(str(rec["analysis_dir"]))
        cluster8_query_id = rec.get("cluster8_query_id", None)

        top_file = analysis_dir / args.top_hit_file
        if not top_file.exists():
            print(f"[WARN] Missing file for {candidate_id}: {top_file}")
            continue

        df = pd.read_csv(top_file, sep="\t")

        # Ensure numeric columns
        for col in ["rank", "probability", "evalue", "score", "seq_identity", "aln_length", "q_coverage", "t_coverage"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        class_scores, evidence = infer_scores_for_candidate(
            df,
            candidate_id=candidate_id,
            cluster8_query_id=cluster8_query_id
        )

        dna_score, glycan_score, virion_score = compute_module_scores(class_scores)

        row = {
            "candidate_id": candidate_id,
            "cluster8_query_id": cluster8_query_id,
            "analysis_dir": str(analysis_dir),
            **class_scores,
            "DNA_context_score": dna_score,
            "glycan_score": glycan_score,
            "virion_score": virion_score,
        }

        row["primary_group"] = infer_primary_group(row)
        row["inferred_model"] = infer_model(row)

        if "notes" in manifest.columns:
            row["notes"] = rec.get("notes", "")

        matrix_rows.append(row)
        all_evidence.extend(evidence)

    matrix = pd.DataFrame(matrix_rows)

    group_order = {
        "C5-MT + 5mYOX-associated": 1,
        "dCMP hydroxymethylase-associated": 2,
        "ThyX/ThyS / nucleotide-metabolism-associated": 3,
        "additional-GT/glycan-associated": 4,
        "ambiguous/no DNA-context": 5,
    }

    matrix["group_rank"] = matrix["primary_group"].map(group_order).fillna(99)

    matrix = matrix.sort_values(
        by=["group_rank", "DNA_context_score", "glycan_score", "virion_score", "candidate_id"],
        ascending=[True, False, False, False, True]
    )

    matrix_out = outdir / "cluster8_neighborhood_matrix.tsv"
    matrix.to_csv(matrix_out, sep="\t", index=False)

    evidence_out = outdir / "cluster8_neighborhood_evidence_long.tsv"
    pd.DataFrame(all_evidence).to_csv(evidence_out, sep="\t", index=False)

    print(f"[OK] Wrote: {matrix_out}")
    print(f"[OK] Wrote: {evidence_out}")
    print(f"[OK] Rows: {len(matrix)}")


if __name__ == "__main__":
    main()
