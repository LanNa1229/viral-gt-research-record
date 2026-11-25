"""Verify frozen inputs, replay only the historical baseline and render panel B in scratch."""
import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys
from urllib.parse import quote
from urllib.request import urlopen
import xml.etree.ElementTree as ET

PACKAGE = Path(__file__).resolve().parent
STEM = "Figure3C_B_paired_umaps"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t" if path.suffix == ".tsv" else ","))

def keyed(rows, column):
    result = {row[column]: row for row in rows}
    assert len(result) == len(rows) and all(result)
    return result

def restore(scratch, checkouts):
    pins = read(PACKAGE / "inputs/PUBLIC_DEPENDENCIES.tsv")
    files = {}
    for pin in pins:
        dest = scratch / "dependencies" / (pin["role"] + Path(pin["repository_path"]).suffix)
        dest.parent.mkdir(parents=True, exist_ok=True)
        if checkouts:
            shutil.copyfile(checkouts[pin["commit_sha"]] / pin["repository_path"], dest)
        else:
            url = f'https://raw.githubusercontent.com/{pin["repository"]}/{pin["commit_sha"]}/{quote(pin["repository_path"])}'
            with urlopen(url, timeout=90) as source, dest.open("wb") as target:
                shutil.copyfileobj(source, target)
        assert dest.stat().st_size == int(pin["size_bytes"]) and sha(dest) == pin["sha256"], pin["role"]
        files[pin["role"]] = dest
    return files

def verify_data(files):
    import numpy as np
    original_rows = read(files["joint_umap"])
    original = keyed(original_rows, "sample_id")
    points = np.array([[float(r["UMAP1"]), float(r["UMAP2"])] for r in original_rows])
    assert len(original) == 3058 and np.array_equal(points, np.load(files["coordinates"], allow_pickle=False))
    refs = {k:r for k,r in original.items() if r["sample_type"] == "reference_BGC"}
    assert len(refs) == 2502
    order = read(PACKAGE / "inputs/viral_clustering_inputs/selected_287_input_order.csv")
    assert order == read(files["selected_order"])
    representatives = read(files["representative_manifest"])
    assert len(representatives) == 556
    assert [int(r["original_556_row"]) for r in order] == [i for i,r in enumerate(representatives) if r["candidate_status"] == "candidate"]
    gc = keyed(read(files["global_assignments"]), "sample_id")
    vc = keyed(read(PACKAGE / "outputs/viral_community_assignments/selected287_baseline.csv"), "sample_id")
    sources = keyed(read(PACKAGE / "figures/source_data/FIGURE_SOURCE_DATA.tsv"), "sample_id")
    ref_meta = keyed(read(files["reference_metadata"]), "BGC")
    sac = {k for k,r in refs.items() if "Saccharide" in r["reference_class"].split()}
    master = keyed(read(files["chemotype_master"]), "mibig_accession")
    assignments = keyed(read(files["chemotype_accessions"]), "mibig_accession")
    mapping = keyed(read(files["chemotype_mapping"]), "source_primary_chemotype")
    assert set(refs) == set(master) == set(ref_meta)
    assert len(sac) == 195 and sac == set(assignments)
    assert set(sources) == set(refs) | set(vc) == set(gc)
    assert len(sources) == 2789 and len(vc) == 287
    for i, item in enumerate(order):
        reference = representatives[int(item["original_556_row"])]
        assert int(item["selected_row"]) == i and item["viral_locus_id"] in vc
        for key in ("viral_locus_id", "representative_window_id"):
            assert item[key] == reference[key] == original[item["viral_locus_id"]][key]
    for key, row in sources.items():
        for field in ("UMAP1", "UMAP2", "sample_type", "viral_locus_id"):
            assert row[field] == original[key][field]
        assert row["global_leiden_community"] == gc[key]["k20_r1_s42"]
        if key in refs:
            assert row["reference_class_original"] == original[key]["reference_class"] == ref_meta[key]["label"]
            assert row["reference_primary_chemotype"] == master[key]["primary_chemotype"]
            assert (row["is_saccharide_reference"] == "True") == (key in sac)
            if key in sac:
                primary = assignments[key]["primary_chemotype"]
                assert row["reference_primary_chemotype"] == primary
                assert row["reference_chemotype_visualization_group"] == assignments[key]["display_chemotype_group"] == mapping[primary]["display_chemotype_group"]
        else:
            assert row["viral_only_leiden_community"] == vc[key]["k20_r1_s42"]
    counts = dict(sorted(Counter(r["k20_r1_s42"] for r in vc.values()).items()))
    assert counts == dict(zip(("VC01", "VC02", "VC03", "VC04", "VC05"), (68,55,74,35,55)))
    assert np.load(files["reference_vectors"], allow_pickle=False).shape == (2502, 320)
    return order, vc, {"references": 2502, "saccharide": 195, "other_references": 2307, "viral": 287, "counts": counts}

def replay_cluster(scratch, files, order, vc):
    code = scratch / "clustering/scripts"
    shutil.copytree(PACKAGE / "scripts/clustering", code)
    shutil.copyfile(PACKAGE / "config.json", code.parent / "config.json")
    sys.path.insert(0, str(code))
    import engine
    import numpy as np
    all_vectors = np.load(files["representative_vectors"], allow_pickle=False)
    assert all_vectors.shape == (556,320)
    ids = [r["viral_locus_id"] for r in order]
    x = all_vectors[[int(r["original_556_row"]) for r in order]]
    sim, distance, inverse = engine.distances(x)
    graph, diagnostics, _, _ = engine.graph(sim, distance, ids, 20, inverse)
    saved = read(PACKAGE / "outputs/viral_community_assignments/selected287_k20_edges.csv")
    edges = [(ids[e.source], ids[e.target]) for e in graph.es]
    assert edges == [(r["source"], r["target"]) for r in saved]
    error = max(abs(e["weight"] - float(r["weight"])) for e,r in zip(graph.es, saved))
    assert error <= 1e-12
    labels, info = engine.leiden(graph, ids, 1.0, 42, "VC")
    assert all(label == vc[key]["k20_r1_s42"] for key,label in zip(ids, labels))
    return {"edges":len(saved), "max_weight_error":error, "components":diagnostics["n_components"],
            "exact_memberships":len(ids), "optimizer_calls":info["actual_optimizer_calls"]}

def verify_svg(path):
    ns = {"s":"http://www.w3.org/2000/svg"}
    root = ET.parse(path).getroot()
    counts = Counter()
    for element in root.findall(".//s:g", ns):
        gid = element.get("id", "")
        for panel in ("B_left", "B_right"):
            if gid.startswith(panel + "_"):
                counts[panel] += len(element.findall(".//s:use", ns)) or len(element.findall("s:path", ns))
    assert counts == {"B_left":2789, "B_right":2789}
    assert not root.findall(".//s:image", ns)
    return dict(counts)

def replay_plot(scratch, arial, bold):
    from PIL import Image, ImageChops
    target = scratch / "plot"
    target.mkdir()
    for name in ("plot_designs.py", "config.json"):
        shutil.copyfile(PACKAGE / "scripts/plotting" / name, target / name)
    shutil.copyfile(PACKAGE / "figures/source_data/FIGURE_SOURCE_DATA.tsv", target / "FIGURE_SOURCE_DATA.tsv")
    config_path = target / "config.json"
    config = json.loads(config_path.read_text())
    historical = json.loads((PACKAGE / "provenance/HISTORICAL_EXECUTION.json").read_text())
    for key,path in {"font_path":arial, "bold_font_path":bold}.items():
        config[key] = str(path.resolve())
    exact_fonts = all(sha(Path(config[k])) == v for k,v in historical["font_sha256"].items())
    config_path.write_text(json.dumps(config), encoding="utf-8")
    spec = importlib.util.spec_from_file_location("historical_plot", target / "plot_designs.py")
    plot = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plot)
    fig, receipts, groups, layout = plot.build(plot.read("FIGURE_SOURCE_DATA.tsv"), "B")
    panels = {name:{r["sample_id"]:(r["UMAP1"],r["UMAP2"]) for r in receipts if r["panel"] == name}
              for name in ("B_left","B_right")}
    assert len(receipts) == 5578 and panels["B_left"] == panels["B_right"] and len(panels["B_left"]) == 2789
    for ext in ("png", "svg", "pdf"):
        fig.savefig(target / f"{STEM}.{ext}", dpi=600 if ext == "png" else 100)
    counts = verify_svg(target / f"{STEM}.svg")
    original = PACKAGE / "figures/final" / f"{STEM}.png"
    with Image.open(original) as saved, Image.open(target / f"{STEM}.png") as generated:
        assert saved.size == generated.size
        same = ImageChops.difference(saved.convert("RGB"), generated.convert("RGB")).getbbox() is None
    return {"panels":counts, "matched_font_bytes":exact_fonts, "PNG_pixels_identical":same, "coordinates_equal_between_panels":True}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scratch", type=Path)
    parser.add_argument("--checkout", action="append", default=[], metavar="SHA=PATH")
    parser.add_argument("--arial", type=Path)
    parser.add_argument("--arial-bold", type=Path)
    parser.add_argument("--deps", action="append", default=[], help="Existing Python package directory")
    parser.add_argument("--stage", choices=["all","data","clustering","plotting"], default="all")
    args = parser.parse_args()
    args.scratch = args.scratch.resolve()
    if args.stage in ("all","plotting") and not (args.arial and args.arial_bold):
        parser.error("plotting requires --arial and --arial-bold")
    sys.path.extend(args.deps)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
        os.environ[key] = "1"
    args.scratch.mkdir(parents=True, exist_ok=False)
    checkouts = {k:Path(v) for k,v in (a.split("=",1) for a in args.checkout)}
    files = restore(args.scratch, checkouts)
    order, vc, totals = verify_data(files)
    result = {"data":totals,"dependency_hashes_verified":len(files),
              "dependency_access":"local pinned checkout" if checkouts else "anonymous immutable HTTPS",
              "original_SVG":verify_svg(PACKAGE / "figures/final" / f"{STEM}.svg")}
    if args.stage in ("all","clustering"):
        result["clustering"] = replay_cluster(args.scratch, files, order, vc)
    if args.stage in ("all","plotting"):
        result["plotting"] = replay_plot(args.scratch, args.arial, args.arial_bold)
    (args.scratch / "verification.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
