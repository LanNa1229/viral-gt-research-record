"""Summarize frozen viral-window inference without treating missing scores as negatives.

Fractions describe a selected structural/GT candidate cohort, never all-virus
prevalence. Retrieval P@k stays a separate descriptive metric. Original window
observations remain unchanged; overlap merging is used only for locus summaries.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
from pathlib import Path
import time

from .inference import digest, read_tsv, save_json, write_tsv


CI_NOTE = "95% Wilson intervals are descriptive; selected and overlapping observations violate independent random-sampling assumptions"
SCOPE = "selected structure/GT-associated viral candidate cohort; not prevalence among all viruses or taxonomic enrichment"
SUMMARY_FIELDS = ["level", "n_units_total", "n_any_analyzed", "n_fully_analyzed", "n_positive", "n_negative_fully_analyzed", "n_unknown", "n_resolved", "positive_fraction_resolved", "wilson95_low", "wilson95_high", "positive_fraction_total_lower_bound", "positive_fraction_total_upper_bound", "n_windows_total", "n_windows_analyzed", "denominator_definition", "uncertainty_note", "population_scope"]
GENOME_DEFINITION = "exact whole-genome nucleotide SHA256 + source taxid; missing taxid or nucleotide hash keeps accession-specific units"
LOCUS_DEFINITION = "connected overlapping inclusive nucleotide window segments within a genome unit; source origin-spanning CDS gaps excluded; adjacency alone does not merge"
LOCUS_FIELDS = ["locus_id", "genome_unit_id", "genome_accession", "source_genome_accessions", "counterpart_alias_count", "genome_sequence_sha256", "taxid", "organism", "source_organism_labels", "start_gene_index", "end_gene_index", "start_1based", "end_1based", "n_windows", "n_analyzed_windows", "n_positive_windows", "classification_status", "fully_analyzed", "window_ids", "anchor_protein_ids", "genome_grouping_definition", "merge_definition"]
GENOME_FIELDS = ["genome_unit_id", "genome_sequence_sha256", "taxid", "source_genome_accessions", "counterpart_alias_count", "source_organism_labels", "n_windows", "n_analyzed_windows", "n_positive_windows", "classification_status", "fully_analyzed", "window_ids", "genome_grouping_definition"]


LOCUS_FIELDS += ['genomic_segments_1based']


def as_bool(value) -> bool:
    return str(value).lower() in {"true", "1", "yes"}


def wilson(k: int, n: int) -> tuple[float | str, float | str]:
    if not n:
        return "", ""
    z = 1.959963984540054
    p = k / n
    denom = 1 + z*z/n
    center = (p + z*z/(2*n)) / denom
    width = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
    return max(0.0, center - width), min(1.0, center + width)


def unit_status(window_ids, predictions: dict) -> dict:
    ids = list(window_ids)
    analyzed = [wid for wid in ids if wid in predictions]
    positive = [wid for wid in analyzed if as_bool(predictions[wid]["saccharide_predicted"])]
    complete = len(analyzed) == len(ids)
    status = "positive" if positive else "negative" if complete and ids else "unknown"
    return dict(window_ids=ids, n_windows=len(ids), n_analyzed=len(analyzed), n_positive=len(positive), complete=complete, status=status)


def genome_key(window: dict) -> tuple[str, str, str]:
    sequence_hash = window.get("genome_sequence_sha256", "").strip()
    taxid = window.get("taxid", "").strip()
    if sequence_hash and taxid:
        return "sequence_and_taxid", sequence_hash, taxid
    return "accession_specific", window["genome_accession"], sequence_hash + ":" + taxid


def genome_unit_id(key: tuple) -> str:
    return "vgg_" + hashlib.sha256(json.dumps(key, separators=(",", ":")).encode()).hexdigest()[:24]


def group_genomes(windows: list[dict], predictions: dict) -> tuple[list[dict], list[dict]]:
    grouped = collections.defaultdict(list)
    for window in windows:
        grouped[genome_key(window)].append(window)
    units, rows = [], []
    for key, members in sorted(grouped.items()):
        state = unit_status(sorted(w["window_id"] for w in members), predictions)
        accessions = sorted({w["genome_accession"] for w in members})
        labels = sorted({w.get("organism", "") for w in members} - {""})
        state.update(genome_key=key, unit_id=genome_unit_id(key), source_genome_accessions=accessions)
        units.append(state)
        rows.append(dict(genome_unit_id=state["unit_id"], genome_sequence_sha256=members[0].get("genome_sequence_sha256", ""), taxid=members[0].get("taxid", ""), source_genome_accessions=";".join(accessions), counterpart_alias_count=len(accessions)-1, source_organism_labels=";".join(labels), n_windows=state["n_windows"], n_analyzed_windows=state["n_analyzed"], n_positive_windows=state["n_positive"], classification_status=state["status"], fully_analyzed=state["complete"], window_ids=";".join(state["window_ids"]), genome_grouping_definition=GENOME_DEFINITION))
    return units, rows


def summarize(level: str, units: list[dict]) -> dict:
    total = len(units)
    positive = sum(u["status"] == "positive" for u in units)
    negative = sum(u["status"] == "negative" for u in units)
    resolved = positive + negative
    unknown = total - resolved
    lo, hi = wilson(positive, resolved)
    return dict(level=level, n_units_total=total, n_any_analyzed=sum(u["n_analyzed"] > 0 for u in units), n_fully_analyzed=sum(u["complete"] for u in units), n_positive=positive, n_negative_fully_analyzed=negative, n_unknown=unknown, n_resolved=resolved, positive_fraction_resolved=positive/resolved if resolved else "", wilson95_low=lo, wilson95_high=hi, positive_fraction_total_lower_bound=positive/total if total else "", positive_fraction_total_upper_bound=(positive+unknown)/total if total else "", n_windows_total=sum(u["n_windows"] for u in units), n_windows_analyzed=sum(u["n_analyzed"] for u in units), denominator_definition="observed positive units plus fully analyzed units with no positive; partially analyzed units without a positive remain unknown", uncertainty_note=CI_NOTE, population_scope=SCOPE)


def combine_segments(segments):
    result=[]
    for start,end in sorted(segments):
        if not result or start>result[-1][1]:result.append([start,end])
        else:result[-1][1]=max(result[-1][1],end)
    return result


def window_segments(window,members):
    """Avoid filling the chromosome-sized gap in a source origin-spanning CDS."""
    regular=[];wrapped=[]
    for gene in members:
        parts=json.loads(gene.get('location_parts','[]'))
        crossing=(len(parts)>1 and ((parts[0][2]==1 and parts[0][0]>parts[-1][0])
            or (parts[0][2]==-1 and parts[0][0]<parts[-1][0])))
        if crossing:wrapped.extend((int(p[0]),int(p[1])) for p in parts)
        elif gene.get('start_1based'):regular.append((int(gene['start_1based']),int(gene['end_1based'])))
    if not wrapped:return [[int(window['start_1based']),int(window['end_1based'])]]
    if regular:wrapped.append((min(p[0] for p in regular),max(p[1] for p in regular)))
    return combine_segments(wrapped)


def merge_loci(windows: list[dict], genes: list[dict], predictions: dict) -> tuple[list[dict], list[dict]]:
    grouped_genes = collections.defaultdict(list)
    for gene in genes:
        grouped_genes[gene["window_id"]].append(gene)
    by_genome = collections.defaultdict(list)
    segmented={}
    for window in windows:
        members = grouped_genes[window["window_id"]]
        if not members:
            raise ValueError(f"Window lacks gene rows: {window['window_id']}")
        indices = [int(g["gene_index"]) for g in members]
        if len(indices) != int(window["n_genes"]) or len(set(indices)) != len(indices):
            raise ValueError(f"Window gene rows disagree with manifest: {window['window_id']}")
        start, end = int(window["start_1based"]), int(window["end_1based"])
        if start < 1 or end < start:
            raise ValueError(f"Invalid nucleotide interval: {window['window_id']}")
        by_genome[genome_key(window)].append((start, end, window))
        segmented[window['window_id']]=window_segments(window,members)
    loci, units = [], []
    for key, intervals in sorted(by_genome.items()):
        parent=list(range(len(intervals)))
        def find(i):
            while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
            return i
        for i,(_,_,left) in enumerate(intervals):
            for j in range(i):
                right=intervals[j][2]
                if any(a<=d and c<=b for a,b in segmented[left['window_id']] for c,d in segmented[right['window_id']]):
                    parent[find(i)]=find(j)
        groups=collections.defaultdict(list)
        for i,item in enumerate(intervals):groups[find(i)].append(item)
        components=[[min(x[0] for x in values),max(x[1] for x in values),[x[2] for x in values]] for values in groups.values()]
        for start, end, members in components:
            ids = sorted(w["window_id"] for w in members)
            regions=combine_segments([part for wid in ids for part in segmented[wid]])
            locus_id = "vgl_" + hashlib.sha256(json.dumps([key, regions], separators=(",", ":")).encode()).hexdigest()[:24]
            state = unit_status(ids, predictions)
            state["unit_id"] = locus_id
            state["genome_key"] = key
            units.append(state)
            organisms = sorted({w.get("organism", "") for w in members} - {""})
            taxids = sorted({w.get("taxid", "") for w in members} - {""})
            if len(taxids) > 1:
                raise ValueError(f"Conflicting source taxonomy inside {locus_id}")
            accessions = sorted({w["genome_accession"] for w in members})
            indices = [int(g["gene_index"]) for wid in ids for g in grouped_genes[wid]]
            loci.append(dict(locus_id=locus_id, genome_unit_id=genome_unit_id(key), genome_accession=accessions[0], source_genome_accessions=";".join(accessions), counterpart_alias_count=len(accessions)-1, genome_sequence_sha256=members[0].get("genome_sequence_sha256", ""), taxid=";".join(taxids), organism=organisms[0] if organisms else "", source_organism_labels=";".join(organisms), start_gene_index=min(indices) if len(accessions)==1 else "", end_gene_index=max(indices) if len(accessions)==1 else "", start_1based=start, end_1based=end, n_windows=len(ids), n_analyzed_windows=state["n_analyzed"], n_positive_windows=state["n_positive"], classification_status=state["status"], fully_analyzed=state["complete"], window_ids=";".join(ids), anchor_protein_ids=";".join(sorted({w["anchor_protein_id"] for w in members})), genome_grouping_definition=GENOME_DEFINITION, merge_definition=LOCUS_DEFINITION))
            loci[-1]['genomic_segments_1based']=json.dumps(regions)
    return loci, units


def taxonomy_cache(path: Path | None) -> dict:
    if path is None:
        return {}
    rows = read_tsv(path)
    cache = {}
    for row in rows:
        taxid = row.get("taxid", row.get("taxonomy_id", ""))
        if not taxid:
            raise ValueError(f"Taxonomy cache row lacks stable taxid: {path}")
        if taxid in cache and row != cache[taxid]:
            raise ValueError(f"Conflicting taxonomy cache rows for {taxid}")
        cache[taxid] = row
    return cache


def make_taxonomy(windows: list[dict], loci_units: list[dict], genome_units: list[dict], predictions: dict, cache: dict) -> list[dict]:
    rows = []
    for rank in ["family", "genus", "species"]:
        groups = collections.defaultdict(list)
        for window in windows:
            info = cache.get(window.get("taxid", ""), {})
            name = info.get(rank, "")
            identifier = info.get(rank + "_taxid", "")
            groups[(name or "Unresolved", identifier, "rank_cache" if name else "rank_unresolved")].append(window)
        for (name, identifier, source), members in sorted(groups.items()):
            member_ids = {w["window_id"] for w in members}
            # Whole locus/genome membership belongs to one taxon. Never count
            # only a scored subfragment as a fully observed negative unit.
            levels = dict(window=[unit_status([w["window_id"]], predictions) for w in members], locus=[u for u in loci_units if set(u["window_ids"]) <= member_ids], genome=[u for u in genome_units if set(u["window_ids"]) <= member_ids])
            for level, units in levels.items():
                rows.append(dict(rank=rank, taxon=name, rank_taxid=identifier, rank_source=source, source_taxids=";".join(sorted({w.get("taxid", "") for w in members} - {""})), **summarize(level, units)))
    return rows


def plots(root: Path, run: Path, windows: list[dict], predictions: dict, genome_units: list[dict], summaries: list[dict]) -> dict:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.ticker import PercentFormatter, MaxNLocator
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "pdf.fonttype": 42})
    base = root / "analysis/viral_glycan_benchmark_20261016"
    required = ["mibig_CORE_ONLY_embeddings.npy", "mibig_matched_reference_metadata.tsv", "pca_CORE_ONLY_fit.npz", "projection_coordinates.tsv"]
    historical_checksums = {r["path"]: r["sha256"] for r in read_tsv(base / "artifact_checksums.tsv")}
    for name in required:
        if name not in historical_checksums or digest(base / name) != historical_checksums[name]:
            raise ValueError(f"Reference/PCA artifact no longer matches verified manifest: {name}")
    reference = np.load(base / required[0], allow_pickle=False)
    metadata = read_tsv(base / required[1])
    if reference.shape != (2502, 320) or len(metadata) != 2502 or not np.isfinite(reference).all():
        raise ValueError("Invalid reference population for PCA overlay")
    if any(int(r["embedding_row"]) != i for i, r in enumerate(metadata)):
        raise ValueError("Reference/PCA metadata row mismatch")
    with np.load(base / required[2], allow_pickle=False) as fit:
        components, mean, variance = (fit[name].copy() for name in ["components", "mean", "explained_variance_ratio"])
    if components.shape != (2, 320) or mean.shape != (320,) or variance.shape != (2,) or not all(np.isfinite(v).all() for v in [components, mean, variance]):
        raise ValueError("Invalid saved PCA fit")
    if not np.allclose(mean, reference.mean(0), atol=1e-5, rtol=1e-5):
        raise ValueError("Saved PCA mean mismatches current CORE_ONLY references")
    reference_xy = (reference - mean) @ components.T
    oldcoords = {r["sample_id"]: (float(r["x"]), float(r["y"])) for r in read_tsv(base / required[3]) if r["method"] == "PCA" and r["strategy"] == "CORE_ONLY" and r["kind"] == "reference"}
    expected_xy = np.asarray([oldcoords[r["BGC"]] for r in metadata])
    if not np.allclose(reference_xy, expected_xy, atol=1e-4, rtol=1e-4):
        raise ValueError("Saved CORE_ONLY PCA reference coordinates fail reproduction")
    query = np.load(run / "viral_bgc_embeddings.npy", allow_pickle=False)
    manifest = read_tsv(run / "embedding_manifest.tsv")
    if query.shape != (len(manifest), 320) or not np.isfinite(query).all() or any(int(r["embedding_row"]) != i for i, r in enumerate(manifest)):
        raise ValueError("Viral embedding row mismatch")
    if {r["window_id"] for r in manifest} != set(predictions):
        raise ValueError("Prediction/embedding observation populations differ")
    query_xy = (query - mean) @ components.T
    sac = np.asarray(["Saccharide" in r["label"].split() for r in metadata])
    query_positive = np.asarray([as_bool(predictions[r["window_id"]]["saccharide_predicted"]) for r in manifest], dtype=bool)
    coordinates = [dict(kind="reference", sample_id=r["BGC"], x=float(xy[0]), y=float(xy[1]), saccharide_classifier_positive="", saccharide_source_label=bool(flag)) for r, xy, flag in zip(metadata, reference_xy, sac)]
    coordinates += [dict(kind="viral_window", sample_id=r["window_id"], x=float(xy[0]), y=float(xy[1]), saccharide_classifier_positive=bool(flag), saccharide_source_label="") for r, xy, flag in zip(manifest, query_xy, query_positive)]
    write_tsv(run / "pca_core_only_coordinates.tsv", coordinates, ["kind", "sample_id", "x", "y", "saccharide_classifier_positive", "saccharide_source_label"])
    generated = []

    def save(fig, name):
        for suffix in ["png", "pdf"]:
            path = run / f"{name}.{suffix}"
            fig.savefig(path, bbox_inches="tight", facecolor="white")
            generated.append(path.name)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 5.8), layout="constrained")
    ax.scatter(reference_xy[~sac, 0], reference_xy[~sac, 1], color="#B9C0C6", s=9, alpha=.45, linewidths=0, label=f"Other reference labels (n={(~sac).sum():,})")
    ax.scatter(reference_xy[sac, 0], reference_xy[sac, 1], color="#087F8C", s=15, alpha=.8, linewidths=0, label=f"Saccharide reference labels (n={sac.sum():,})")
    for flag, marker, color, label in [(False, "o", "#DDAF35", "Window: other classifier result"), (True, "*", "#D1495B", "Window: Saccharide classifier result")]:
        selected = query_positive == flag
        if selected.any():
            ax.scatter(query_xy[selected, 0], query_xy[selected, 1], marker=marker, color=color, s=80 if not flag else 135, edgecolor="#202629", linewidth=.65, zorder=5, label=f"{label} (n={selected.sum():,})")
    ax.set(xlabel=f"PC1 ({100*variance[0]:.1f}% reference variance)", ylabel=f"PC2 ({100*variance[1]:.1f}% reference variance)", title="Viral windows in the matched CORE_ONLY reference space")
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("Reference-only PCA fit; quantitative retrieval uses 320D cosine", fontsize=9, y=1.025)
    save(fig, "reference_viral_embedding_overlay")

    fig, ax = plt.subplots(figsize=(7.4, 5), layout="constrained")
    for i, row in enumerate(summaries):
        if row["n_resolved"]:
            p, lo, hi = (float(row[k]) for k in ["positive_fraction_resolved", "wilson95_low", "wilson95_high"])
            ax.errorbar(i, p, yerr=[[p-lo], [hi-p]], fmt="o", color="#087F8C", capsize=5, markersize=7)
            ax.annotate(f"{row['n_positive']}/{row['n_resolved']} resolved\n{row['n_unknown']} unknown", (i, p), xytext=(0, 12), textcoords="offset points", ha="center", fontsize=9)
        else:
            ax.text(i, .5, "No resolved units", ha="center", fontsize=9)
    ax.set_xticks(range(len(summaries)), [r["level"].capitalize() for r in summaries])
    ax.set(xlim=(-.5, len(summaries)-.5), ylim=(-.04, 1.22), ylabel="Fraction with a Saccharide classifier result", title="Saccharide-associated fractions in the selected candidate cohort")
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_yticks([0, .25, .5, .75, 1])
    fig.supxlabel("95% Wilson intervals are descriptive; observations overlap and are selected.\nUnknown units are excluded from the resolved denominator.", fontsize=9)
    save(fig, "saccharide_associated_fractions")

    by_id = {w["window_id"]: w for w in windows}
    taxon_counts = collections.defaultdict(collections.Counter)
    for unit in genome_units:
        source = by_id[unit["window_ids"][0]]
        label = f"{source.get('organism') or 'Unknown organism'} [taxid {source.get('taxid') or 'unresolved'}]"
        taxon_counts[label][unit["status"]] += 1
    ranked = sorted(taxon_counts, key=lambda name: (-sum(taxon_counts[name].values()), name))
    shown = ranked[:12]
    if len(ranked) > 12:
        remaining = collections.Counter()
        for name in ranked[12:]:
            remaining.update(taxon_counts[name])
        taxon_counts["Remaining source taxa"] = remaining
        shown.append("Remaining source taxa")
    fig, ax = plt.subplots(figsize=(9, max(3.8, .48*len(shown)+1.6)), layout="constrained")
    left = np.zeros(len(shown))
    for status, color, label, hatch in [("positive", "#087F8C", "Observed Saccharide result", ""), ("negative", "#B9C0C6", "Fully analyzed; no Saccharide result", ""), ("unknown", "#F1E2B9", "Unknown / incomplete inference", "//")]:
        values = np.asarray([taxon_counts[name][status] for name in shown])
        ax.barh(range(len(shown)), values, left=left, color=color, label=label, edgecolor="#485158", linewidth=.4, hatch=hatch)
        left += values
    ax.set_yticks(range(len(shown)), shown)
    ax.invert_yaxis()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set(xlabel="Genome units (whole nucleotide sequence + source taxid)", title="Source-taxonomy distribution of the selected candidate genomes", xlim=(0, max(float(left.max()) if len(left) else 0, 1)*1.12))
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(.5, -.18), fontsize=8)
    fig.supxlabel("Organism labels and taxids come from source records; family/genus/species ranks require a taxonomy cache.", fontsize=9)
    save(fig, "viral_taxonomic_distribution")
    return dict(generated_files=generated, pca_fit=str(base / required[2]), pca_fit_sha256=digest(base / required[2]), pca_reference_coordinate_regression="PASS atol=rtol=1e-4", reference_n=len(metadata), query_n=len(manifest), matplotlib=matplotlib.__version__, visual_review="not yet performed; inspect rendered exports before claiming publication quality")


def analyze(root: Path, run: Path, taxonomy_path: Path | None = None, make_plots: bool = True) -> dict:
    started = time.monotonic()
    windows = read_tsv(run / "viral_gt_windows.tsv")
    genes = read_tsv(run / "viral_gt_window_genes.tsv")
    scores = read_tsv(run / "viral_bgc_predictions.tsv")
    if len({w["window_id"] for w in windows}) != len(windows):
        raise ValueError("Duplicate extraction window identities")
    predictions = {r["window_id"]: r for r in scores}
    if len(predictions) != len(scores) or not set(predictions) <= {w["window_id"] for w in windows}:
        raise ValueError("Duplicate or unmatched inference windows")
    # Independently verify original threshold and label decision for all scores.
    for row in scores:
        probability = float(row["saccharide_probability"])
        threshold = max(.5, float(row["probability_Other"]))
        if not 0 <= probability <= 1 or abs(threshold - float(row["threshold"])) > 1e-7 or (probability >= threshold) != as_bool(row["saccharide_predicted"]):
            raise ValueError(f"Classifier decision mismatch: {row['window_id']}")
    loci, locus_units = merge_loci(windows, genes, predictions)
    genome_units, genome_rows = group_genomes(windows, predictions)
    summaries = [summarize("window", [unit_status([w["window_id"]], predictions) for w in windows]), summarize("locus", locus_units), summarize("genome", genome_units)]
    if taxonomy_path is None and (run / "taxonomy_ranks.tsv").exists():
        taxonomy_path = run / "taxonomy_ranks.tsv"
    taxonomic = make_taxonomy(windows, locus_units, genome_units, predictions, taxonomy_cache(taxonomy_path))
    write_tsv(run / "candidate_loci.tsv", loci, LOCUS_FIELDS)
    write_tsv(run / "genome_units.tsv", genome_rows, GENOME_FIELDS)
    write_tsv(run / "three_level_summary.tsv", summaries, SUMMARY_FIELDS)
    write_tsv(run / "viral_taxonomy_summary.tsv", taxonomic, ["rank", "taxon", "rank_taxid", "rank_source", "source_taxids"] + SUMMARY_FIELDS)
    plots_provenance = plots(root, run, windows, predictions, genome_units, summaries) if make_plots else None
    inputs = ["viral_gt_windows.tsv", "viral_gt_window_genes.tsv", "viral_bgc_predictions.tsv", "embedding_manifest.tsv", "viral_bgc_embeddings.npy"]
    result = dict(status="PASS", window_n=len(windows), analyzed_window_n=len(predictions), locus_n=len(loci), genome_n=len(genome_units), source_genome_accession_n=len({w["genome_accession"] for w in windows}), counterpart_aliases_collapsed=sum(r["counterpart_alias_count"] for r in genome_rows), unknown_window_n=len(windows)-len(predictions), taxonomy_rank_source=str(taxonomy_path) if taxonomy_path else "unresolved; no rank cache supplied", taxonomy_rank_cache_sha256=digest(taxonomy_path) if taxonomy_path else None, population_scope=SCOPE, uncertainty_note=CI_NOTE, input_checksums={name: digest(run / name) for name in inputs}, plots=plots_provenance, elapsed_seconds=time.monotonic()-started, definitions=dict(genome=GENOME_DEFINITION, locus=LOCUS_DEFINITION, counterpart_alias_count="number of additional distinct source accessions beyond the first representative; requires whole-genome sequence identity and same nonempty source taxid", gene_indices="retained as within-record annotation coordinates only; blank for multi-accession loci because annotation indices may differ", negative="every constituent window has a prediction and none meets original Saccharide threshold", positive="at least one constituent window meets original frozen classifier threshold", unknown="no observed positive and at least one window has no prediction", retrieval="P@k is separate and descriptive; no retrieval-positive cutoff declared"))
    output_names = ["candidate_loci.tsv", "genome_units.tsv", "three_level_summary.tsv", "viral_taxonomy_summary.tsv"]
    if plots_provenance:
        output_names += plots_provenance["generated_files"] + ["pca_core_only_coordinates.tsv"]
    result["output_checksums"] = {name: digest(run / name) for name in output_names}
    save_json(run / "analysis_provenance.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--taxonomy-cache", type=Path)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    result = analyze(args.root.resolve(), args.run.resolve(), args.taxonomy_cache, not args.no_plots)
    print(json.dumps({key: result[key] for key in ["status", "window_n", "analyzed_window_n", "locus_n", "genome_n", "unknown_window_n", "elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
