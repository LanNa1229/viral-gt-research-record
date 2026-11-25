"""Verify or replay the frozen reference-only supplementary workflow."""
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import hashlib
import json
import math
import re
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request


ROOT = Path(__file__).resolve().parents[1]
PROV = ROOT / "provenance"
MASTER = ROOT / "outputs/reference_chemotype_master"
EVIDENCE = ROOT / "annotation_evidence/evidence_tables"
DECISIONS = ROOT / "annotation_evidence/curation_records"
SOURCE = ROOT / "figures/source_data"
FIGURE = ROOT / "figures/final"
STEM = "SuppFig2_saccharide_reference_chemotype_distribution"


def read(path, delimiter="\t"):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def index(rows, key):
    indexed = {r[key]: r for r in rows}
    assert len(indexed) == len(rows), key
    return indexed


def check(path, row):
    assert path.stat().st_size == int(row["size_bytes"]), path
    assert sha(path) == row["sha256"], path


def check_archive():
    manifest = PROV / "ARCHIVE_SHA256.tsv"
    rows = index(read(manifest), "path")
    files = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*") if p.is_file() and p != manifest}
    assert files == set(rows), (files - set(rows), set(rows) - files)
    for name, row in rows.items():
        check(ROOT / name, row)
    return {"files": len(files) + 1, "bytes": sum(int(r["size_bytes"]) for r in rows.values()) + manifest.stat().st_size,
            "manifest_sha256": sha(manifest)}


def pull(url, dst, row):
    if not dst.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(url, headers={"User-Agent": "Saccharide-chemotype-reproduction/1"})
        with tempfile.NamedTemporaryFile(dir=dst.parent, delete=False) as handle:
            tmp = Path(handle.name)
            try:
                with urllib.request.urlopen(request, timeout=90) as reply:
                    shutil.copyfileobj(reply, handle)
                handle.flush()
                check(tmp, row)
            except BaseException:
                tmp.unlink()
                raise
        tmp.rename(dst)
    check(dst, row)
    return dst


def dependencies(work):
    paths = {}
    for row in read(PROV / "PUBLIC_DEPENDENCIES.tsv"):
        url = f"https://raw.githubusercontent.com/{row['repository']}/{row['commit']}/{row['path']}"
        paths[row["role"]] = pull(url, work / Path(row["path"]).name, row)
    return paths


def mibig_snapshot(work):
    row = next(r for r in read(PROV / "EXTERNAL_DEPENDENCIES.tsv") if r["role"] == "mibig_json")
    return pull(row["url"], work / Path(row["url"]).name, row)


def composition(deps, json_path):
    master = index(read(MASTER / "reference_2502_chemotype_master.tsv"), "mibig_accession")
    metadata = index(read(deps["matched_reference_metadata"]), "BGC")
    assignments = index(read(deps["global_leiden_assignments"], ","), "sample_id")
    crosswalk = index(read(ROOT / "inputs/community_assignments/REFERENCE_COMMUNITY_CROSSWALK.tsv"), "mibig_accession")
    assert len(master) == 2502 and set(master) == set(metadata) == set(crosswalk) <= set(assignments)
    assert len(assignments) == 2789
    decisions = index(load(DECISIONS / "CURATION_DECISIONS.json")["decisions"], "mibig_accession")
    mapping = index(read(SOURCE / "CHEMOTYPE_GROUP_MAPPING.tsv"), "source_primary_chemotype")
    assert {int(r["category_order"]) for r in mapping.values()} == set(range(1, 8))
    assert all((r["in_four_category_surface_definition"] == "True") == (int(r["category_order"]) <= 4)
               for r in mapping.values())
    trace = index(read(SOURCE / "ACCESSION_LEVEL_CHEMOTYPE_MAPPING.tsv"), "mibig_accession")
    selected = {}
    buckets = defaultdict(list)
    hybrids = 0
    frozen_json = {}
    priority = set()
    with tarfile.open(json_path, "r:gz") as archive:
        for member in archive:
            acc = Path(member.name).stem
            if member.isfile() and member.name.endswith(".json") and acc in master:
                assert acc not in frozen_json
                raw = archive.extractfile(member).read()
                assert json.loads(raw)["cluster"]["mibig_accession"] == acc
                frozen_json[acc] = hashlib.sha256(raw).hexdigest()
    assert set(frozen_json) == set(master)
    for acc, ref in master.items():
        origin = metadata[acc]
        label = origin["label"]
        tokens = label.split()
        community = assignments[acc]["k20_r1_s42"]
        assert ref["source_record_SHA256"] == frozen_json[acc]
        assert origin["embedding_row"] == ref["original_embedding_row"]
        assert label == ref["original_bgcprophet_class"] and tokens == json.loads(ref["original_class_multilabel"])
        assert community == ref["global_leiden_community"]
        assert all(crosswalk[acc][key] == ref[key] for key in crosswalk[acc])
        has_saccharide = "Saccharide" in tokens
        assert (ref["is_saccharide_reference"] == "True") == has_saccharide
        if has_saccharide or community in ("GC08", "GC05"):
            priority.add(acc)
            assert decisions[acc]["primary_chemotype"] == ref["primary_chemotype"]
        if not has_saccharide:
            continue
        group = community if community in ("GC08", "GC05") else "Other GCs"
        entry = mapping[ref["primary_chemotype"]]
        rank = int(entry["category_order"])
        assert all(trace[acc][k] == ref[k] for k in ref)
        assert trace[acc]["figure_population"] == group
        assert trace[acc]["display_chemotype_group"] == entry["display_chemotype_group"]
        assert trace[acc]["display_category_order"] == entry["category_order"]
        membership = trace[acc]["four_category_surface_membership"]
        assert (membership == "IN_DEFINED_SET") == (entry["in_four_category_surface_definition"] == "True"), membership
        selected[acc] = group
        buckets[group, rank].append(acc)
        hybrids += len(tokens) > 1
    assert set(decisions) == priority and len(priority) == 229
    assert set(selected) == set(trace) and len(selected) == 195 and hybrids == 65

    source = {(r["population"], int(r["category_order"])): r for r in read(SOURCE / "FIGURE_SOURCE_DATA.tsv")}
    assert len(source) == 21
    expected = {"GC08": [31, 19, 18, 4, 20, 4, 4],
                "GC05": [4, 2, 3, 0, 7, 5, 3], "Other GCs": [1, 1, 1, 0, 7, 51, 10]}
    summary = {}
    recalculated_rows = []
    by_rank = {}
    for entry in mapping.values():
        rank = int(entry["category_order"])
        pair = (entry["display_chemotype_group"], entry["legend_label"])
        assert by_rank.setdefault(rank, pair) == pair
    for position, (group, target) in enumerate(expected.items(), 1):
        n = sum(g == group for g in selected.values())
        counts = [len(buckets[group, rank]) for rank in range(1, 8)]
        assert counts == target and sum(counts) == n
        left = 0
        for rank, count in enumerate(counts, 1):
            row = source[group, rank]
            assert int(row["population_order"]) == position
            assert (row["chemotype_group"], row["legend_label"]) == by_rank[rank]
            assert int(row["n_BGCs"]) == count and int(row["denominator_saccharide_BGCs"]) == n
            checks = {"fraction": Fraction(count, n), "percentage": Fraction(100 * count, n),
                      "left_percentage": Fraction(100 * left, n),
                      "right_percentage": Fraction(100 * (left + count), n)}
            for field, value in checks.items():
                assert math.isclose(float(row[field]), float(value), rel_tol=0, abs_tol=1e-9), row
            accessions = json.loads(row["accessions_json"])
            assert len(accessions) == count and sorted(accessions) == sorted(buckets[group, rank])
            recalculated_rows.append(dict(population=group, category=by_rank[rank][0], count=count, denominator=n,
                                          fraction=float(checks["fraction"]), percentage=float(checks["percentage"])))
            left += count
        top = sum(counts[:4])
        summary[group] = dict(counts=counts, denominator=n, surface_count=top, surface_fraction=top / n)
    assert [summary[group]["surface_count"] for group in expected] == [72, 9, 3]
    return dict(references=len(master), unique_saccharide=len(selected), hybrids=hybrids, priority=len(priority),
                rows_validated=len(source), populations=summary, recalculated_rows=recalculated_rows)


def dated_run(path, cwd):
    with (cwd / (path.stem + ".replay.stdout.txt")).open("wb") as log:
        subprocess.run([sys.executable, "-B", str(path)], cwd=cwd, check=True, stdout=log)


def replay_figure(work, font):
    font_row = next(r for r in read(PROV / "EXTERNAL_DEPENDENCIES.tsv") if r["role"] == "original_arial")
    assert font.is_file(), font
    font_matches = sha(font) == font_row["sha256"]
    work.mkdir(exist_ok=False)
    copied = ROOT / "scripts/figure_generation"
    shutil.copyfile(copied / "plot_chemotype_bars.py", work / "plot_chemotype_bars.py")
    shutil.copyfile(SOURCE / "FIGURE_SOURCE_DATA.tsv", work / "FIGURE_SOURCE_DATA.tsv")
    shutil.copyfile(ROOT / "outputs/validation/figure/INPUT_VALIDATION.json", work / "INPUT_VALIDATION.json")
    config = load(copied / "config.json")
    config["font_path"] = str(font.resolve())
    write_json(work / "config.json", config)
    dated_run(work / "plot_chemotype_bars.py", work)

    import xml.etree.ElementTree as ET
    import re
    ns = {"s": "http://www.w3.org/2000/svg"}
    def geometry(path):
        root = ET.parse(path).getroot()
        segments = {}
        for group in root.findall(".//s:g", ns):
            tag = group.get("id", "")
            if not tag.startswith("segment_"):
                continue
            geom = group.find("s:path", ns)
            values = re.findall(r"[-+]?\d+(?:\.\d+)?", geom.attrib["d"])
            segments[tag] = (values, geom.attrib["style"])
        text = [("".join(item.itertext()), item.get("style")) for item in root.findall(".//s:text", ns)]
        return root.attrib["viewBox"], segments, text
    reproduced = work / "figures" / STEM
    original = FIGURE / STEM
    layout, segments, labels = geometry(reproduced.with_suffix(".svg"))
    old_layout, old_segments, old_labels = geometry(original.with_suffix(".svg"))
    assert layout == old_layout and len(segments) == 21 and segments == old_segments
    assert labels == old_labels
    assert all(reproduced.with_suffix(f".{ext}").stat().st_size > 0 for ext in ("pdf", "svg", "png"))
    return dict(original_font_sha256_matches=font_matches, matching_svg_segments=len(segments),
                matching_svg_labels=len(labels), png_byte_match=sha(reproduced.with_suffix(".png")) == sha(original.with_suffix(".png")),
                reproduced_files=["figures/" + STEM + "." + ext for ext in ("pdf", "svg", "png")])


def publication_bytes(path):
    """Hash a reconstructed table after removing source-record contact strings."""
    data = re.sub(rb"(?i)\b[a-z0-9_.%+-]+@[a-z0-9.-]+\.[a-z]{2,}\b", b"[contact redacted]", path.read_bytes())
    return data


def replay_reference(work, deps, json_path, gbk_path):
    work.mkdir(exist_ok=False)
    raw = work / "raw_sources"
    raw.mkdir()
    json_dest = raw / "mibig_json_3.1"
    json_dest.mkdir()
    shutil.copyfile(json_path, raw / "mibig_json_3.1.tar.gz")
    shutil.copyfile(gbk_path, raw / gbk_path.name)
    meta = index(read(deps["matched_reference_metadata"]), "BGC")
    with tarfile.open(json_path, "r:gz") as archive:
        for member in archive:
            name = Path(member.name).name
            if member.isfile() and name.endswith(".json") and Path(name).stem in meta:
                (json_dest / name).write_bytes(archive.extractfile(member).read())
    prep = ROOT / "scripts/preprocessing"
    annot = ROOT / "scripts/chemotype_annotation"
    for filename, source in (("extract_gbk_annotations.py", prep), ("curate_priority_records.py", annot),
                             ("curation_common.py", annot)):
        shutil.copyfile(source / filename, work / filename)
    shutil.copyfile(EVIDENCE / "PRIORITY_SOURCE_RECORDS.json", work / "PRIORITY_SOURCE_RECORDS.json")
    shutil.copyfile(EVIDENCE / "PUBMED_CITATIONS.json", work / "PUBMED_PRIMARY_RECORDS.json")
    shutil.copyfile(DECISIONS / "CURATED_CHEMICAL_FEATURES.json", work / "CURATED_CHEMICAL_FEATURES.json")
    dated_run(work / "extract_gbk_annotations.py", work)
    for name, frozen in (("PRIORITY_EVIDENCE_CARDS.json", EVIDENCE / "PRIORITY_EVIDENCE_CARDS.json"),
                         ("GBK_RECORD_MANIFEST.tsv", PROV / "historical_reference/GBK_RECORD_MANIFEST.tsv")):
        assert publication_bytes(work / name).replace(b"\r\n", b"\n") == frozen.read_bytes().replace(b"\r\n", b"\n"), name
    dated_run(work / "curate_priority_records.py", work)
    ran = work / "CURATION_DECISIONS.json"
    original = DECISIONS / ran.name
    assert ran.read_bytes().replace(b"\r\n", b"\n") == original.read_bytes().replace(b"\r\n", b"\n")

    hist = (annot / "build_reference_database.py").read_text(encoding="utf-8")
    # Use the original reference construction. Replace the two workstation bindings and
    # omit only input provenance and the separate, unpublished viral-link calculation.
    start = hist.index(" archived=read_json(OUT/'SOURCE_PREFLIGHT.json')")
    end = hist.index("\n", hist.index("source_manifest.append(dict(role='existing_viral_master_input'", start))
    code = hist[:start] + hist[end:]
    start = code.index(" # Join original Top-20 relationships")
    end = code.index(" write_table('reference_2502_chemotype_master.tsv'", start)
    code = code[:start] + code[end:]
    code = code.replace(";write_table('viral_287_reference_chemotype_links.tsv',links)", "")
    start = code.index(" print('Built'")
    end = code.index("\nif __name__", start)
    code = code[:start] + code[end:]
    code = code.replace(" meta=read_table(ARCHIVE/'inputs/workspace/analysis/viral_glycan_benchmark_20261016/mibig_matched_reference_metadata.tsv')",
                        " meta=read_table(OUT/'mibig_matched_reference_metadata.tsv')")
    code = code.replace("read_table(ARCHIVE/'assignments/global_joint_leiden_grid.csv')",
                        "read_table(OUT/'global_joint_leiden_grid.csv')")
    for role in ("matched_reference_metadata", "global_leiden_assignments"):
        shutil.copyfile(deps[role], work / deps[role].name)
    script = work / "replay_reference_build.py"
    script.write_text(code, encoding="utf-8")
    dated_run(script, work)
    names = sorted(path.name for path in MASTER.glob("*.tsv"))
    for name in names:
        assert hashlib.sha256(publication_bytes(work / name)).hexdigest() == sha(MASTER / name), name
    assert sha(work / "chemotype_curation_audit.tsv") == sha(EVIDENCE / "chemotype_curation_audit.tsv")
    return dict(evidence_cards_match_after_contact_redaction_and_lf_normalization=True,
                gbk_record_manifest_match=True, frozen_curated_decisions_match_after_lf_normalization=True,
                output_tables_match_after_contact_redaction=names + ["chemotype_curation_audit.tsv"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("archive", "composition", "figure", "reference"))
    parser.add_argument("work", type=Path, nargs="?")
    parser.add_argument("--font", type=Path)
    parser.add_argument("--gbk", type=Path)
    options = parser.parse_args()
    if options.action == "archive":
        result = check_archive()
    else:
        if options.work is None:
            parser.error("provide a scratch directory outside the archive")
        work = options.work.resolve()
        if work.is_relative_to(ROOT):
            parser.error("scratch directory must be outside the archive")
        work.mkdir(parents=True, exist_ok=True)
        if options.action == "figure":
            if options.font is None:
                parser.error("provide --font with the path to licensed Arial regular")
            result = replay_figure(work / "figure", options.font)
        else:
            deps = dependencies(work / "dependencies")
            json_path = mibig_snapshot(work / "dependencies")
            if options.action == "composition":
                result = composition(deps, json_path)
            else:
                row = next(r for r in read(PROV / "EXTERNAL_DEPENDENCIES.tsv") if r["role"] == "mibig_gbk")
                gbk = options.gbk or pull(row["url"], work / "dependencies/mibig_gbk_3.1.tar.gz", row)
                check(gbk, row)
                assert gbk.name == "mibig_gbk_3.1.tar.gz"
                result = replay_reference(work / "reference", deps, json_path, gbk)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
