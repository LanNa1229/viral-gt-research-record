"""Restartable frozen BGC-Prophet CORE_ONLY inference; historical inputs stay read-only.

Run with the existing WSL recovery environment and PYTHONPATH=src. No model is
downloaded, trained or modified. The original validated forward is imported by
file path; actual window genes are neither padded nor repeated.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import time
from typing import Any

try:
    import resource
except ImportError:  # Analysis can reuse these file helpers on Windows.
    resource = None


CLASSIFIER_SHA = "f45ca298980d4740dd43131a7e69c99cb233368c815470572cfc46f33abf4eb5"
ESM_SHA = "46f002a9870c9bdecd0ea887acb1f9a38a6b561e8f8bf8a6990b679b9d31b928"
SOURCE_SHA = "fa4e1d045ac2c94559a9ab41054fd6b28d991d82be63e98585d7bc014564dca8"
ENCODER_SOURCE_SHA = "fa4ad0c224707ee316b3ab3f37c102796ac176b7203fad357b99c331178a65f0"
REFERENCE_FILES = {
    "mibig_CORE_ONLY_embeddings.npy": "78036e1636eaeeab90f4d157794c868d9b57033195939d83960d8e70bb831e22",
    "mibig_matched_reference_metadata.tsv": "5184f8f3199da2d1fac90e85ff11e2a30b443b70aac60992136718856c18ba4c",
    "reference_validation.json": "a75043f212d4df3297e2dfa37145127f5e890f866d4c36cc6284eeed1780306e",
    "mibig_original_regression.tsv": "4622855130a83bb97c51f1680aa7d423b0cdba2482236cc5ac341c735b7d6cef",
    "pilot_validation.json": "8ca48c1c1cf2e587aaed694dd66f62811f7d811f3f0432387348a3e2110c01b7",
    "esm_provenance.json": "f1d149751ff26261b97251d99f4a8b0a9b98351a2713104923a41a968c715fb7",
    "scripts/benchmark.py": "d0b7608341ef31ad29bae58b147804c7a6a0ecddcbc6b8b43801e255ed4b9233",
}
LABELS = ["Alkaloid", "Terpene", "NRP", "Polyketide", "RiPP", "Saccharide", "Other"]
HISTORICAL = [
    "viral_glycan_benchmark_20261016",
    "viral_glycan_benchmark_A0A0K0N728_20261016",
    "viral_glycan_benchmark_cluster7_ARF08726_ARF08739_20261016",
]
PREDICTION_FIELDS = ["window_id", "genome_accession", "anchor_protein_id", "n_genes", "strategy"] + [
    "probability_" + name for name in LABELS
] + ["saccharide_probability", "saccharide_rank", "predicted_labels", "threshold", "saccharide_predicted", "probability_semantics"]
MANIFEST_FIELDS = ["embedding_row", "window_id", "genome_accession", "anchor_protein_id", "n_genes", "strategy", "pooling", "membership_mask", "ordered_protein_ids", "ordered_sequence_sha256", "input_cache_key", "cache_reused", "classifier_sha256", "esm_sha256"]
QC_FIELDS = ["window_id", "status", "reason", "n_genes", "input_shape", "finite", "pooled_matches_classifier_input", "actual_and_published_pooling_equal", "cache_reused", "n_truncated_proteins", "n_sequence_cache_hits"]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def save_json(path: Path, value: Any) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def atomic_numpy(path: Path, value: Any) -> None:
    import numpy as np
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as stream:
        np.save(stream, value, allow_pickle=False)
    os.replace(temporary, path)


def load_frozen(root: Path):
    source = root / "analysis/viral_glycan_benchmark_20261016/scripts/benchmark.py"
    spec = importlib.util.spec_from_file_location("viral_gt_frozen_20261016", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def valid_vector(vector: Any, shape: tuple[int, ...]) -> bool:
    import numpy as np
    return vector.shape == shape and vector.dtype == np.float32 and bool(np.isfinite(vector).all())


def seed_cache(root: Path) -> tuple[dict, list[dict]]:
    """Read only the three audited caches, retaining model provenance checks."""
    import numpy as np
    cache: dict = {}
    provenance = []
    for name in HISTORICAL:
        directory = root / "analysis" / name
        meta = directory / "protein_embedding_metadata.tsv"
        array = directory / "protein_esm2_embeddings.npy"
        proof_path = directory / ("esm_provenance.json" if name == HISTORICAL[0] else "inference_validation.json")
        proof = json.loads(proof_path.read_text())
        proof_hash = proof.get("esm_checkpoint_sha256", proof.get("sha256"))
        if proof_hash != ESM_SHA:
            raise ValueError(f"Unverified historical encoder identity: {proof_path}")
        rows, vectors = read_tsv(meta), np.load(array, allow_pickle=False)
        if not valid_vector(vectors, (len(rows), 320)):
            raise ValueError(f"Invalid historical cache: {array}")
        for index, row in enumerate(rows):
            if int(row["embedding_row"]) != index:
                raise ValueError(f"Historical cache row misalignment: {meta}")
            key = row["sequence_sha256"]
            if key in cache and not np.array_equal(cache[key], vectors[index]):
                raise ValueError(f"Conflicting existing vectors for sequence {key}")
            cache[key] = vectors[index]
        provenance.append(dict(metadata=str(meta), metadata_sha256=digest(meta), array=str(array), array_sha256=digest(array), rows=len(rows), proof_sha256=digest(proof_path)))
    return cache, provenance


def prepare_windows(run: Path) -> tuple[list[tuple[dict, list[dict]]], list[dict], dict]:
    windows = read_tsv(run / "viral_gt_windows.tsv")
    genes = read_tsv(run / "viral_gt_window_genes.tsv")
    grouped: dict[str, list] = collections.defaultdict(list)
    for gene in genes:
        grouped[gene["window_id"]].append(gene)
    if len({w["window_id"] for w in windows}) != len(windows):
        raise ValueError("Duplicate window_id in viral_gt_windows.tsv")
    unknown = set(grouped) - {w["window_id"] for w in windows}
    if unknown:
        raise ValueError(f"Gene rows reference unknown windows: {sorted(unknown)[:5]}")
    eligible, qc, sequences = [], [], {}
    for window in windows:
        reason = ""
        ordered = sorted(grouped[window["window_id"]], key=lambda g: int(g["relative_position"]))
        if window.get("inference_eligible", "false").lower() not in {"true", "1", "yes"}:
            reason = "extraction_marked_ineligible"
        elif not 1 <= len(ordered) <= 15 or int(window["n_genes"]) != len(ordered):
            reason = "gene_count_mismatch_or_not_1_to_15"
        elif len({int(g["relative_position"]) for g in ordered}) != len(ordered):
            reason = "duplicate_relative_position"
        elif len({g["gene_index"] for g in ordered}) != len(ordered):
            reason = "duplicate_genomic_occurrence_in_window"
        else:
            for gene in ordered:
                sequence = gene["sequence"]
                actual = hashlib.sha256(sequence.encode()).hexdigest()
                if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO") or actual != gene["sequence_sha256"]:
                    reason = "missing_invalid_or_hash_mismatched_sequence"
                    break
        if reason:
            qc.append(dict(window_id=window["window_id"], status="excluded", reason=reason, n_genes=len(ordered)))
            continue
        eligible.append((window, ordered))
        for gene in ordered:
            sequences[gene["sequence_sha256"]] = gene["sequence"]
    return eligible, qc, sequences


def encode_missing(torch, esm_checkpoint: Path, sequences: dict, cache: dict,
                   cache_dir: Path, batch_size: int) -> dict:
    """Persist each new vector atomically; restart reuses completed sequences."""
    import numpy as np
    missing = []
    local_hits = 0
    initial_hits = sum(key in cache for key in sequences)
    for key, sequence in sequences.items():
        if key in cache:
            continue
        path = cache_dir / (key + ".npy")
        if path.exists():
            vector = np.load(path, allow_pickle=False)
            if valid_vector(vector, (320,)):
                cache[key] = vector
                local_hits += 1
                continue
            raise ValueError(f"Invalid sequence cache; review before replacement: {path}")
        missing.append((key, sequence))
    started = time.monotonic()
    if missing:
        import esm
        model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(esm_checkpoint))
        model.eval().float()
        converter = alphabet.get_batch_converter(truncation_seq_length=1022)
        with torch.inference_mode():
            for start in range(0, len(missing), batch_size):
                batch = missing[start:start + batch_size]
                _, _, tokens = converter(batch)
                hidden = model(tokens, repr_layers=[6], return_contacts=False)["representations"][6]
                for index, (key, sequence) in enumerate(batch):
                    vector = hidden[index, 1:min(len(sequence), 1022) + 1].mean(0).cpu().numpy()
                    if not valid_vector(vector, (320,)):
                        raise ValueError(f"Invalid generated ESM vector: {key}")
                    atomic_numpy(cache_dir / (key + ".npy"), vector)
                    cache[key] = vector
                print(f"ESM completed {min(start + batch_size, len(missing))}/{len(missing)} new sequences", flush=True)
        del model
    return dict(unique_sequences=len(sequences), historical_cache_hits=initial_hits, persistent_cache_hits=local_hits, newly_encoded=len(missing), esm_seconds=time.monotonic() - started, esm_batch_size=batch_size)


def load_window_cache(path: Path, key: str):
    import numpy as np
    if not path.exists():
        return None
    with np.load(path, allow_pickle=False) as result:
        if str(result["input_key"].item()) != key:
            raise ValueError(f"Window cache identity mismatch: {path}")
        values = tuple(result[name].copy() for name in ("scores", "embedding", "membership_embedding"))
    if not all(valid_vector(v, shape) for v, shape in zip(values, [(7,), (320,), (320,)])):
        raise ValueError(f"Invalid cached inference result: {path}")
    if not ((values[0] >= 0).all() and (values[0] <= 1).all() and np.array_equal(values[1], values[2])):
        raise ValueError(f"Cached inference probability or pooling mismatch: {path}")
    return values


def save_window_cache(path: Path, key: str, values: tuple) -> None:
    import numpy as np
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as stream:
        np.savez(stream, input_key=np.asarray(key), scores=values[0], embedding=values[1], membership_embedding=values[2])
    os.replace(temporary, path)


def verify_reference(root: Path) -> dict:
    """Fail closed before retrieval if the audited representation space changed."""
    base = root / "analysis/viral_glycan_benchmark_20261016"
    historical = {r["path"]: r["sha256"] for r in read_tsv(base / "artifact_checksums.tsv")}
    for name, expected in REFERENCE_FILES.items():
        if historical.get(name) != expected or digest(base / name) != expected:
            raise ValueError(f"Audited CORE_ONLY reference/proof checksum mismatch: {name}")
    fixed_resources = [
        (base / "sources/classifier.pt", CLASSIFIER_SHA),
        (root / "sources/upstream/bgc_prophet/train/classifier.py", SOURCE_SHA),
        (root / "sources/upstream/bgc_prophet/train/model.py", ENCODER_SOURCE_SHA),
        (root / "analysis/bgcprophet_sequence_recovery_20260912/model_cache/hub/checkpoints/esm2_t6_8M_UR50D.pt", ESM_SHA),
    ]
    for path, expected in fixed_resources:
        if digest(path) != expected:
            raise ValueError(f"Reference/query checkpoint or forward source mismatch: {path}")
    proof = json.loads((base / "reference_validation.json").read_text())
    pilot = json.loads((base / "pilot_validation.json").read_text())
    encoder = json.loads((base / "esm_provenance.json").read_text())
    if proof.get("status") != "PASS" or proof.get("n_biological_bgcs") != 2502 or not proof.get("original_regression_passed"):
        raise ValueError("CORE_ONLY biological population or regression proof failed")
    if pilot.get("status") != "PASS" or pilot.get("checkpoint_sha256") != CLASSIFIER_SHA or pilot.get("classifier_source_sha256") != SOURCE_SHA or not pilot.get("weights_unchanged"):
        raise ValueError("CORE_ONLY frozen checkpoint/forward proof failed")
    if encoder.get("sha256") != ESM_SHA or encoder.get("layer") != 6 or encoder.get("truncation") != 1022 or encoder.get("pooling") != "residue mean excluding BOS/EOS":
        raise ValueError("CORE_ONLY encoder/preprocessing proof failed")
    regression = read_tsv(base / "mibig_original_regression.tsv")
    if len(regression) != 20 or not all(r["passed"].lower() == "true" and float(r["atol"]) == 1e-4 and float(r["rtol"]) == 1e-4 for r in regression):
        raise ValueError("Original-checkpoint reference regression proof failed")
    return dict(status="PASS", fixed_artifact_sha256=REFERENCE_FILES, classifier_sha256=CLASSIFIER_SHA, esm_sha256=ESM_SHA, classifier_source_sha256=SOURCE_SHA, encoder_source_sha256=ENCODER_SOURCE_SHA, pooling="actual isolated core genes; original all-position mean; all-one additive membership", reference_rows=2502)


def localize(root: Path, run: Path, embeddings, manifest: list[dict], predictions: list[dict]) -> dict:
    """Descriptive 320D retrieval with the matched isolated-core references."""
    import numpy as np
    compatibility = verify_reference(root)
    base = root / "analysis/viral_glycan_benchmark_20261016"
    array = base / "mibig_CORE_ONLY_embeddings.npy"
    metadata_path = base / "mibig_matched_reference_metadata.tsv"
    reference, metadata = np.load(array, allow_pickle=False), read_tsv(metadata_path)
    if not valid_vector(reference, (2502, 320)) or len(metadata) != 2502:
        raise ValueError("Invalid CORE_ONLY reference population")
    if any(int(row["embedding_row"]) != i for i, row in enumerate(metadata)):
        raise ValueError("Reference metadata row misalignment")
    norms = np.linalg.norm(reference.astype(np.float64), axis=1)
    if (norms == 0).any():
        raise ValueError("Zero reference vector cannot be used for cosine retrieval")
    normalized = reference.astype(np.float64) / norms[:, None]
    ids = np.asarray([row["BGC"] for row in metadata])
    positive = np.asarray(["Saccharide" in row["label"].split() for row in metadata])
    rows = []
    for query, row, prediction in zip(embeddings, manifest, predictions):
        norm = np.linalg.norm(query.astype(np.float64))
        if norm == 0:
            raise ValueError(f"Zero query vector: {row['window_id']}")
        similarity = normalized @ (query.astype(np.float64) / norm)
        order = np.lexsort((ids, -similarity))
        nearest_positive = next(int(i) for i in order if positive[i])
        rows.append(dict(window_id=row["window_id"], embedding_row=row["embedding_row"], strategy="CORE_ONLY", predicted_saccharide=prediction["saccharide_predicted"], saccharide_probability=prediction["saccharide_probability"], nearest_reference=ids[order[0]], nearest_reference_label=metadata[order[0]]["label"], nearest_reference_cosine=float(similarity[order[0]]), saccharide_fraction_top5=float(positive[order[:5]].mean()), saccharide_fraction_top10=float(positive[order[:10]].mean()), saccharide_fraction_top20=float(positive[order[:20]].mean()), descriptive_saccharide_enrichment_top10=float(positive[order[:10]].mean()/positive.mean()), nearest_saccharide_reference=ids[nearest_positive], nearest_saccharide_cosine=float(similarity[nearest_positive]), reference_saccharide_fraction=float(positive.mean()), interpretation="descriptive supervised-reference localization and top10/global-reference ratio; no significance inference, positivity cutoff or chemical validation"))
    fields = ["window_id", "embedding_row", "strategy", "predicted_saccharide", "saccharide_probability", "nearest_reference", "nearest_reference_label", "nearest_reference_cosine", "saccharide_fraction_top5", "saccharide_fraction_top10", "saccharide_fraction_top20", "descriptive_saccharide_enrichment_top10", "nearest_saccharide_reference", "nearest_saccharide_cosine", "reference_saccharide_fraction", "interpretation"]
    write_tsv(run / "saccharide_localization.tsv", rows, fields)
    return dict(reference_array=str(array), reference_array_sha256=digest(array), reference_metadata_sha256=digest(metadata_path), compatibility_gate=compatibility, n_reference_bgcs=len(metadata), n_saccharide_positive=int(positive.sum()), pooling="isolated-core all-actual-position mean", retrieval="cosine in 320D; accession tie break; descriptive P@k", n_queries=len(rows))


def run_inference(root: Path, run: Path, esm_batch_size: int = 1, with_localization: bool = False, sequence_cache_root: Path | None = None) -> dict:
    import numpy as np
    started = time.monotonic()
    root, run = root.resolve(), run.resolve()
    run.mkdir(parents=True, exist_ok=True)
    if esm_batch_size < 1:
        raise ValueError("esm_batch_size must be positive")
    compatibility = verify_reference(root)
    frozen = load_frozen(root)
    classifier = root / "analysis/viral_glycan_benchmark_20261016/sources/classifier.pt"
    esm_checkpoint = root / "analysis/bgcprophet_sequence_recovery_20260912/model_cache/hub/checkpoints/esm2_t6_8M_UR50D.pt"
    source = root / "sources/upstream/bgc_prophet/train/classifier.py"
    for path, expected in [(classifier, CLASSIFIER_SHA), (esm_checkpoint, ESM_SHA), (source, SOURCE_SHA)]:
        if digest(path) != expected:
            raise ValueError(f"Checkpoint/source checksum mismatch: {path}")
    versions = {name: importlib.metadata.version(name) for name in ["torch", "numpy", "fair-esm"]}
    settings = dict(strategy="CORE_ONLY", gene_order="relative_position ascending, genomic coordinate order", membership="all float32 ones; original additive bias", pooling="mean of contextualized actual window genes; original classifier input", padding="none", repetition="none", classifier_sha256=CLASSIFIER_SHA, classifier_source_sha256=SOURCE_SHA, encoder_source_sha256=digest(root / "sources/upstream/bgc_prophet/train/model.py"), frozen_wrapper_sha256=digest(Path(frozen.__file__)), esm_checkpoint_sha256=ESM_SHA, esm_layer=6, esm_max_residues=1022, esm_pooling="residue mean excluding BOS/EOS", device="cpu", dtype="float32", versions=versions)
    cache_root = run / "inference_cache"
    protein_settings = {key: settings[key] for key in ["esm_checkpoint_sha256", "esm_layer", "esm_max_residues", "esm_pooling", "device", "dtype", "versions"]}
    protein_settings["esm_batch_size"] = esm_batch_size
    protein_base = sequence_cache_root.resolve() if sequence_cache_root else cache_root / "proteins"
    protein_dir = protein_base / json_digest(protein_settings)
    window_dir = cache_root / "windows"
    protein_dir.mkdir(parents=True, exist_ok=True)
    window_dir.mkdir(parents=True, exist_ok=True)
    if (protein_dir / "settings.json").exists() and json.loads((protein_dir / "settings.json").read_text()) != protein_settings:
        raise ValueError(f"Existing sequence-cache settings mismatch: {protein_dir}")
    save_json(protein_dir / "settings.json", protein_settings)
    eligible, qc, sequences = prepare_windows(run)
    cache, cache_provenance = seed_cache(root)
    historical_keys = set(cache)
    torch = frozen.torch_setup()
    encoding = encode_missing(torch, esm_checkpoint, sequences, cache, protein_dir, esm_batch_size)
    model, initial_weights = None, None
    scores_rows, manifest, embeddings = [], [], []
    inference_seconds, completed, reused = 0.0, 0, 0
    for window, genes in eligible:
        key = json_digest(dict(settings=settings, protein_settings=protein_settings, ordered_sequence_sha256=[g["sequence_sha256"] for g in genes]))
        values = load_window_cache(window_dir / (key + ".npz"), key)
        cache_reused = values is not None
        if not cache_reused:
            if model is None:
                model = frozen.load_model(torch)
                initial_weights = {name: value.detach().clone() for name, value in model.state_dict().items()}
            tick = time.monotonic()
            values = tuple(frozen.infer(torch, model, [cache[g["sequence_sha256"]] for g in genes], [1] * len(genes), True))
            inference_seconds += time.monotonic() - tick
            if not np.array_equal(values[1], values[2]):
                raise ValueError("Actual-gene and all-position pooling differ with all-one membership")
            save_window_cache(window_dir / (key + ".npz"), key, values)
            completed += 1
        else:
            reused += 1
        probabilities, embedding, _ = values
        labels, threshold = frozen.predicted(probabilities)
        common = dict(window_id=window["window_id"], genome_accession=window.get("genome_accession", ""), anchor_protein_id=window.get("anchor_protein_id", ""), n_genes=len(genes), strategy="CORE_ONLY")
        scores_rows.append(dict(**common, **{"probability_" + label: float(value) for label, value in zip(LABELS, probabilities)}, saccharide_probability=float(probabilities[5]), saccharide_rank=1 + int((probabilities > probabilities[5]).sum()), predicted_labels=labels, threshold=threshold, saccharide_predicted="Saccharide" in labels.split(), probability_semantics="raw frozen-model sigmoid outputs; seven scores need not sum to one"))
        manifest.append(dict(**common, embedding_row=len(embeddings), pooling=settings["pooling"], membership_mask=" ".join(["1"] * len(genes)), ordered_protein_ids=" ".join(g["protein_id"] for g in genes), ordered_sequence_sha256=" ".join(g["sequence_sha256"] for g in genes), input_cache_key=key, cache_reused=cache_reused, classifier_sha256=CLASSIFIER_SHA, esm_sha256=ESM_SHA))
        embeddings.append(embedding)
        qc.append(dict(window_id=window["window_id"], status="PASS", reason="", n_genes=len(genes), input_shape=f"(1,{len(genes)},320)", finite=True, pooled_matches_classifier_input=True, actual_and_published_pooling_equal=True, cache_reused=cache_reused, n_truncated_proteins=sum(len(g["sequence"]) > 1022 for g in genes), n_sequence_cache_hits=sum(g["sequence_sha256"] in historical_keys for g in genes)))
        print(f"Inference {len(embeddings)}/{len(eligible)} {window['window_id']} {'cache' if cache_reused else 'computed'}", flush=True)
    if model is not None and not all(torch.equal(initial_weights[k], value) for k, value in model.state_dict().items()):
        raise ValueError("Model state changed during inference")
    for path, expected in [(classifier, CLASSIFIER_SHA), (esm_checkpoint, ESM_SHA), (source, SOURCE_SHA)]:
        if digest(path) != expected:
            raise ValueError(f"Checkpoint/source changed during inference: {path}")
    array = np.stack(embeddings).astype(np.float32) if embeddings else np.empty((0, 320), dtype=np.float32)
    atomic_numpy(run / "viral_bgc_embeddings.npy", array)
    write_tsv(run / "viral_bgc_predictions.tsv", scores_rows, PREDICTION_FIELDS)
    write_tsv(run / "embedding_manifest.tsv", manifest, MANIFEST_FIELDS)
    write_tsv(run / "inference_qc.tsv", qc, QC_FIELDS)
    localization = localize(root, run, array, manifest, scores_rows) if with_localization else None
    provenance = dict(status="PASS" if eligible else "NO_ELIGIBLE_WINDOWS", settings=settings, root=str(root), run=str(run), historical_sequence_cache=cache_provenance, input_checksums={name: digest(run / name) for name in ["viral_gt_windows.tsv", "viral_gt_window_genes.tsv"]}, n_success=len(embeddings), n_excluded=len(qc) - len(embeddings), n_new_window_inferences=completed, n_resumed_window_inferences=reused, encoding=encoding, classifier_seconds=inference_seconds, elapsed_seconds=time.monotonic() - started, peak_process_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if resource else None, localization=localization, limitations=["Variable-length CORE_ONLY inputs differ from original 128-position training construction", "Training-reference localization is descriptive, not held-out accuracy or chemical validation", "Published LMDB reference vectors and current CPU ESM vectors have historical near-identity rather than bitwise-identity evidence", "Long proteins use the original first-1022-residue convention", "Window observations preserve overlapping anchors and are not independent loci/genomes"])
    provenance["sequence_cache_directory"] = str(protein_dir)
    provenance["sequence_cache_shared_explicitly"] = sequence_cache_root is not None
    provenance["reference_compatibility_gate"] = compatibility
    provenance["output_checksums"] = {name: digest(run / name) for name in ["viral_bgc_predictions.tsv", "viral_bgc_embeddings.npy", "embedding_manifest.tsv", "inference_qc.tsv"] + (["saccharide_localization.tsv"] if with_localization else [])}
    save_json(run / "inference_provenance.json", provenance)
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Original MiBG workspace containing verified resources")
    parser.add_argument("--run", type=Path, required=True, help="Versioned extraction output directory")
    parser.add_argument("--esm-batch-size", type=int, default=1, help="Default 1 reproduces the validated CPU encoder calls")
    parser.add_argument("--localize", action="store_true", help="Also emit descriptive cosine retrieval against matched CORE_ONLY references")
    parser.add_argument("--sequence-cache-root", type=Path, help="Optional shared protein-vector cache root across pilot/full runs; keyed by frozen encoder settings")
    args = parser.parse_args()
    result = run_inference(args.root, args.run, args.esm_batch_size, args.localize, args.sequence_cache_root)
    print(json.dumps({k: result[k] for k in ["status", "n_success", "n_excluded", "n_new_window_inferences", "n_resumed_window_inferences", "elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
