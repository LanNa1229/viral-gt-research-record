"""Bounded, isolated repair of archived Tupanvirus and UniParc crosslinks.

This never changes the baseline HTTP/genome caches, input manifests or runs.
Only exact original-anchor sequence hashes may enter the supplemental aliases.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request

from .inference import digest, read_tsv, save_json, write_tsv

ARCHIVED = ["A0A6N1NGR4", "A0A6N1NK79", "A0A6N1NM99", "A0A6N1NPT3", "A0A6N1NQZ0", "A0A6N1NUE9", "A0A6N1NUS0"]
UNIPARC = ["UPI001E761E12", "UPI001FBC14AE"]


def supplement(root: Path) -> dict:
    from Bio import SeqIO
    from .extract import parse_record
    task = root / "analysis/viral_gt_mining_20261017"
    cache = task / "cache/uniparc_supplement"
    inputs = task / "inputs"
    cache.mkdir(parents=True, exist_ok=True)
    names = ["supplement_aliases.tsv", "supplement_genbanks.tsv", "supplement_status.tsv", "supplement_anchor_sequences.faa", "supplement_provenance.json"]
    if any((inputs/name).exists() for name in names):
        raise ValueError("Supplement outputs already exist; preserve them and use a new version rather than overwrite")
    retrievals = []

    def fetch(url: str, name: str) -> Path:
        path = cache / name
        cached = path.exists()
        if not cached:
            for attempt in range(4):
                try:
                    req = urllib.request.Request(url, headers={"User-Agent": "viral-gt-mining-supplement/1.0"})
                    with urllib.request.urlopen(req, timeout=45) as response:
                        payload = response.read()
                    if not payload:
                        raise ValueError("Empty API response")
                    path.write_bytes(payload)
                    time.sleep(.4)
                    break
                except Exception:
                    if attempt == 3:
                        raise
                    time.sleep(2**attempt)
        retrievals.append(dict(url=url, path=str(path), sha256=digest(path), bytes=path.stat().st_size, reused=cached))
        print(f"Supplement {'reused' if cached else 'fetched'} {name}", flush=True)
        return path

    def efetch(database: str, accession: str, name: str) -> Path:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(dict(db=database, id=accession, rettype="gbwithparts" if database == "nuccore" else "gb", retmode="text"))
        return fetch(url, name)

    old_aliases = {r["input_accession"]: r for r in read_tsv(inputs / "retrieved_aliases.tsv") if r["input_accession"] in ARCHIVED}
    requests = []
    sequences = {}
    for anchor in ARCHIVED:
        row = old_aliases[anchor]
        requests.append(dict(input_accession=anchor, protein_id=row["protein_id"], expected_sha256=row["sequence_sha256"], taxid=row["taxid"], anchor_evidence_path=row["evidence_path"], anchor_evidence_sha256=row["evidence_sha256"], method="exact_archived_unisave_and_ncbi_coded_by"))
    uniparc_data = {}
    for anchor in UNIPARC:
        path = fetch(f"https://rest.uniprot.org/uniparc/{anchor}.json", anchor + ".json")
        obj = json.loads(path.read_text())
        if obj["uniParcId"] != anchor:
            raise ValueError("UniParc response identity mismatch")
        sequence = obj["sequence"]["value"]
        expected = hashlib.sha256(sequence.encode()).hexdigest()
        sequences[anchor] = sequence
        uniparc_data[anchor] = obj
        if anchor == "UPI001FBC14AE":
            # PDB-only sequence has a terminal LEHHHHHH expression tag.
            # Do not strip it or replace the requested sequence with a native one.
            continue
        for reference in obj["uniParcCrossReferences"]:
            if reference["database"] not in {"RefSeq", "EMBL"} or not reference.get("active"):
                continue
            version = reference.get("version")
            if not isinstance(version, int):
                raise ValueError("UniParc protein cross-reference lacks explicit sequence version")
            protein_id = reference["id"] + "." + str(version)
            requests.append(dict(input_accession=anchor, protein_id=protein_id, expected_sha256=expected, taxid=str(reference["organism"]["taxonId"]), anchor_evidence_path=str(path), anchor_evidence_sha256=digest(path), method="exact_uniparc_sequence_and_ncbi_coded_by"))
    proteins_file = efetch("protein", ",".join(sorted({r["protein_id"] for r in requests})), "supplement_proteins.gb")
    proteins = {r.id: r for r in SeqIO.parse(proteins_file, "genbank")}
    candidates, failures = [], []
    for request in requests:
        protein = proteins.get(request["protein_id"])
        if protein is None:
            failures.append(dict(input_accession=request["input_accession"], reason="exact_versioned_protein_not_returned", protein_id=request["protein_id"]))
            continue
        sequence = str(protein.seq).upper()
        if hashlib.sha256(sequence.encode()).hexdigest() != request["expected_sha256"]:
            failures.append(dict(input_accession=request["input_accession"], reason="original_anchor_vs_protein_sequence_mismatch", protein_id=request["protein_id"]))
            continue
        sequences[request["input_accession"]] = sequence
        coded = [f.qualifiers["coded_by"][0] for f in protein.features if f.type == "CDS" and f.qualifiers.get("coded_by")]
        sources = {acc for text in coded for acc in re.findall(r"([A-Z][A-Z_]*\d+\.\d+):", text)}
        if len(coded) != 1 or len(sources) != 1:
            failures.append(dict(input_accession=request["input_accession"], reason="ambiguous_or_missing_coded_by_genome", protein_id=request["protein_id"]))
            continue
        candidates.append(dict(request, genome_accession=next(iter(sources)), coded_by=coded[0]))
    genome_data = {}
    for accession in sorted({c["genome_accession"] for c in candidates}):
        path = efetch("nuccore", accession, accession + ".gb")
        record = SeqIO.read(path, "genbank")
        if record.id != accession:
            raise ValueError(f"Requested archived genome version not returned: {accession} -> {record.id}")
        info, genes = parse_record(record)
        genome_data[accession] = (path, record, info, genes)
    aliases, accepted_genomes = [], set()
    for candidate in candidates:
        path, record, info, genes = genome_data[candidate["genome_accession"]]
        matches = [g for g in genes if g["protein_id"] == candidate["protein_id"] and g["sequence_sha256"] == candidate["expected_sha256"]]
        if len(matches) != 1 or info["taxid"] != candidate["taxid"]:
            failures.append(dict(input_accession=candidate["input_accession"], reason="exact_CDS_identity_or_taxid_failed", protein_id=candidate["protein_id"]))
            continue
        gene = matches[0]
        if gene["translation_status"] not in {"verified", "verified_partial"}:
            failures.append(dict(input_accession=candidate["input_accession"], reason="source_CDS_translation_not_verified:"+gene["translation_status"], protein_id=candidate["protein_id"]))
            continue
        # Independently check the protein-coded coordinates and strand too.
        feature = record.features[int(gene["source_feature_index"])]
        source_numbers = list(map(int, re.findall(r"\d+", candidate["coded_by"].split(":", 1)[1])))
        if len(feature.location.parts) != 1 or source_numbers != [gene["start_1based"], gene["end_1based"]] or ("complement" in candidate["coded_by"]) != (gene["strand"] == -1):
            failures.append(dict(input_accession=candidate["input_accession"], reason="protein_coded_by_vs_CDS_coordinate_mismatch", protein_id=candidate["protein_id"]))
            continue
        aliases.append(dict(input_accession=candidate["input_accession"], genome_accession=record.id, protein_id=gene["protein_id"], sequence_sha256=gene["sequence_sha256"], evidence_path=candidate["anchor_evidence_path"], evidence_sha256=candidate["anchor_evidence_sha256"], taxid=info["taxid"], mapping_method=candidate["method"], genome_path=str(path), genome_sha256=digest(path), protein_evidence_path=str(proteins_file), protein_evidence_sha256=digest(proteins_file), coded_by=candidate["coded_by"], translation_status=gene["translation_status"]))
        accepted_genomes.add(record.id)
    statuses = []
    for anchor in ARCHIVED + UNIPARC:
        accepted = [a for a in aliases if a["input_accession"] == anchor]
        reason = ";".join(sorted({f["reason"] for f in failures if f["input_accession"] == anchor}))
        if anchor == "UPI001FBC14AE":
            reason = "PDB-only exact sequence; terminal LEHHHHHH expression tag; no authoritative genomic crosslink; no tag stripping or native-sequence substitution"
        statuses.append(dict(input_accession=anchor, status="exact_mapping_recovered" if accepted else "unresolved", n_aliases=len(accepted), exact_sequence_sha256=hashlib.sha256(sequences[anchor].encode()).hexdigest() if anchor in sequences else "", aa_length=len(sequences[anchor]) if anchor in sequences else "", genome_accessions=";".join(sorted({a["genome_accession"] for a in accepted})), reason=reason))
    alias_fields = ["input_accession", "genome_accession", "protein_id", "sequence_sha256", "evidence_path", "evidence_sha256", "taxid", "mapping_method", "genome_path", "genome_sha256", "protein_evidence_path", "protein_evidence_sha256", "coded_by", "translation_status"]
    write_tsv(inputs / names[0], aliases, alias_fields)
    write_tsv(inputs / names[1], [dict(path=str(genome_data[acc][0]), genome_accession=acc, sha256=digest(genome_data[acc][0])) for acc in sorted(accepted_genomes)], ["path", "genome_accession", "sha256"])
    write_tsv(inputs / names[2], statuses, ["input_accession", "status", "n_aliases", "exact_sequence_sha256", "aa_length", "genome_accessions", "reason"])
    with (inputs / names[3]).open("w") as stream:
        for anchor, sequence in sorted(sequences.items()):
            stream.write(f">{anchor} sequence_sha256={hashlib.sha256(sequence.encode()).hexdigest()}\n{sequence}\n")
    provenance = dict(status="completed", original_cohort_preserved=594, requested_inputs=ARCHIVED+UNIPARC, recovered_inputs=sum(r["status"] == "exact_mapping_recovered" for r in statuses), recovered_aliases=len(aliases), genome_records=len(accepted_genomes), retrievals=retrievals, failures=failures, output_sha256={name:digest(inputs/name) for name in names[:-1]}, policy="Only exact original requested sequence and source CDS translation with matching protein coded_by, archived genome version and source taxid accepted. Baseline files/runs/caches unchanged.")
    save_json(inputs / names[4], provenance)
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    result = supplement(args.root.resolve())
    print(json.dumps({k:result[k] for k in ["status", "recovered_inputs", "recovered_aliases", "genome_records"]}, indent=2))


if __name__ == "__main__":
    main()
