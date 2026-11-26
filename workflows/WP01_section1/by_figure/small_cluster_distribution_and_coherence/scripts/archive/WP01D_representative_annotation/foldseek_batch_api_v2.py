#!/usr/bin/env python3
"""
Batch Foldseek Search Server API runner.

Submit many PDB/mmCIF files to the Foldseek web/API server, poll jobs,
save raw JSON results, and write a unified TSV table.

Manifest TSV required columns:
    query_id    structure_file

Useful optional columns:
    protein_id  cluster  source_genome  model_source
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import requests
except ImportError as exc:
    raise SystemExit("Missing dependency: requests. Install with `pip install requests`.") from exc


FOLDSEEK_BASE = "https://search.foldseek.com/api"
REQUIRED_MANIFEST_COLUMNS = {"query_id", "structure_file"}

LEDGER_FIELDS = [
    "run_id", "query_id", "structure_file", "mode", "databases", "ticket_id",
    "status", "submitted_at", "completed_at", "raw_json_file",
    "download_tar_file", "error",
]

BASE_OUTPUT_FIELDS = [
    "run_id", "query_id", "protein_id", "cluster", "source_genome", "model_source",
    "structure_file", "mode", "databases_submitted", "database", "rank",
    "target", "target_id", "target_description", "seq_identity", "aln_length",
    "mismatches", "gaps_opened", "q_start", "q_end", "t_start", "t_end",
    "probability", "evalue", "score", "q_len", "t_len", "q_coverage",
    "t_coverage", "ticket_id", "raw_json_file", "status", "submitted_at",
    "completed_at", "error",
]
ALIGNMENT_OUTPUT_FIELDS = ["q_aln", "t_aln"]


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def setup_logging(outdir: Path, verbose: bool = False) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    handlers = [logging.StreamHandler(sys.stderr), logging.FileHandler(outdir / "run.log")]
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
    )


def read_manifest(path: Path) -> List[Dict[str, str]]:
    """Read a tab-separated manifest.

    Robust behavior:
    - skips completely blank lines;
    - skips comment rows where query_id starts with '#';
    - raises a clear error for partially filled rows.
    """
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Manifest has no header: {path}")
        fieldnames = [f.strip() if isinstance(f, str) else f for f in reader.fieldnames]
        reader.fieldnames = fieldnames
        missing = REQUIRED_MANIFEST_COLUMNS - set(fieldnames)
        if missing:
            raise ValueError(
                f"Manifest is missing required columns: {sorted(missing)}. "
                f"Found columns: {fieldnames}"
            )
        rows: List[Dict[str, str]] = []
        seen = set()
        for line_no, row in enumerate(reader, start=2):
            clean = {k: (v.strip() if isinstance(v, str) else "") for k, v in row.items()}

            # Skip fully blank rows, including rows that contain only tabs/spaces.
            if not any(clean.values()):
                logging.warning("Skipping blank manifest line %s", line_no)
                continue

            query_id = clean.get("query_id", "")
            structure_file = clean.get("structure_file", "")

            # Allow comment rows in the manifest.
            if query_id.startswith("#"):
                logging.warning("Skipping comment manifest line %s", line_no)
                continue

            if not query_id or not structure_file:
                raise ValueError(
                    f"Manifest line {line_no} is partially filled. "
                    f"Both query_id and structure_file are required. Row={clean}"
                )
            if query_id in seen:
                raise ValueError(f"Duplicate query_id in manifest: {query_id}")
            seen.add(query_id)
            rows.append(clean)
        return rows


def resolve_structure_path(structure_file: str, manifest_path: Path) -> Path:
    p = Path(structure_file)
    if p.is_absolute():
        return p
    cwd_path = Path.cwd() / p
    if cwd_path.exists():
        return cwd_path
    return manifest_path.parent / p


def validate_manifest(rows: List[Dict[str, str]], manifest_path: Path) -> None:
    missing = []
    for row in rows:
        p = resolve_structure_path(row["structure_file"], manifest_path)
        if not p.exists():
            missing.append((row["query_id"], str(p)))
        elif p.suffix.lower() not in {".pdb", ".cif", ".mmcif", ".ent"}:
            logging.warning("Unusual structure file extension for %s: %s", row["query_id"], p)
    if missing:
        preview = "\n".join(f"  {qid}: {path}" for qid, path in missing[:20])
        extra = "" if len(missing) <= 20 else f"\n  ... and {len(missing) - 20} more"
        raise FileNotFoundError(f"Missing structure files:\n{preview}{extra}")


def read_ledger(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.exists():
        return {}
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["query_id"]: row for row in reader if row.get("query_id")}


def write_ledger(path: Path, ledger: Dict[str, Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=LEDGER_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for query_id in sorted(ledger):
            writer.writerow({field: ledger[query_id].get(field, "") for field in LEDGER_FIELDS})
    tmp.replace(path)


def safe_filename(text: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in {"-", "_", "."} else "_" for ch in text)
    return safe.strip("_") or "query"


def submit_foldseek_job(
    session: requests.Session,
    structure_path: Path,
    databases: Sequence[str],
    mode: str,
    base_url: str,
    timeout: int,
) -> Dict[str, Any]:
    data: List[Tuple[str, str]] = [("mode", mode)] + [("database[]", db) for db in databases]

    max_retries = 12
    base_sleep = 60

    last_response = None

    for attempt in range(1, max_retries + 1):
        with structure_path.open("rb") as handle:
            files = {"q": (structure_path.name, handle, "application/octet-stream")}
            response = session.post(
                f"{base_url}/ticket",
                files=files,
                data=data,
                timeout=timeout,
            )

        last_response = response

        if response.status_code != 429:
            response.raise_for_status()
            return response.json()

        retry_after = response.headers.get("Retry-After")
        if retry_after and retry_after.isdigit():
            sleep_s = int(retry_after)
        else:
            sleep_s = min(base_sleep * attempt, 900)

        logging.warning(
            "Foldseek API returned 429 Too Many Requests for %s. "
            "Sleeping %s seconds before retry %d/%d.",
            structure_path.name,
            sleep_s,
            attempt,
            max_retries,
        )

        response.close()
        time.sleep(sleep_s)

    if last_response is not None:
        last_response.raise_for_status()

    raise RuntimeError(f"Foldseek submission failed after {max_retries} retries: {structure_path}")

def get_job_status(session: requests.Session, ticket_id: str, base_url: str, timeout: int) -> Dict[str, Any]:
    response = session.get(f"{base_url}/ticket/{ticket_id}", timeout=timeout)
    response.raise_for_status()
    return response.json()


def wait_for_completion(
    session: requests.Session,
    ticket_id: str,
    base_url: str,
    poll_interval: int,
    max_wait: int,
    timeout: int,
) -> Dict[str, Any]:
    start = time.time()
    while True:
        status_data = get_job_status(session, ticket_id, base_url, timeout)
        status = status_data.get("status", "")
        logging.info("ticket=%s status=%s", ticket_id, status)
        if status in {"COMPLETE", "ERROR", "UNKNOWN"}:
            return status_data
        if time.time() - start > max_wait:
            raise TimeoutError(
                f"Foldseek job did not finish within {max_wait} seconds: "
                f"ticket={ticket_id}, last_status={status}"
            )
        time.sleep(poll_interval)


def fetch_result_json(
    session: requests.Session,
    ticket_id: str,
    base_url: str,
    result_index: int,
    timeout: int,
) -> Dict[str, Any]:
    response = session.get(f"{base_url}/result/{ticket_id}/{result_index}", timeout=timeout)
    response.raise_for_status()
    return response.json()


def download_result_tar(
    session: requests.Session,
    ticket_id: str,
    base_url: str,
    out_path: Path,
    timeout: int,
) -> None:
    with session.get(f"{base_url}/result/download/{ticket_id}", stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with out_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 64):
                if chunk:
                    handle.write(chunk)


def flatten_alignments(obj: Any) -> Iterable[Dict[str, Any]]:
    """Yield alignment dictionaries from nested Foldseek alignment structures."""
    if isinstance(obj, dict):
        yield obj
    elif isinstance(obj, list):
        for item in obj:
            yield from flatten_alignments(item)


def split_target(target: str) -> Tuple[str, str]:
    parts = (target or "").split(maxsplit=1)
    if not parts:
        return "", ""
    return parts[0], parts[1] if len(parts) > 1 else ""


def to_int(value: Any) -> Optional[int]:
    if value in {None, ""}:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def coverage(start: Any, end: Any, length: Any) -> str:
    s, e, l = to_int(start), to_int(end), to_int(length)
    if s is None or e is None or l in {None, 0}:
        return ""
    return f"{max(0, e - s + 1) / l:.4f}"


def parse_foldseek_result(
    result_data: Dict[str, Any],
    manifest_row: Dict[str, str],
    ledger_row: Dict[str, str],
    top_n: int,
    include_alignments: bool,
) -> List[Dict[str, Any]]:
    output_rows: List[Dict[str, Any]] = []
    results = result_data.get("results", [])
    if isinstance(results, dict):
        results = [results]
    if not isinstance(results, list):
        logging.warning("Unexpected result structure for %s", manifest_row.get("query_id"))
        return output_rows

    for db_index, db_result in enumerate(results):
        if not isinstance(db_result, dict):
            continue
        db_name = db_result.get("db") or db_result.get("database") or f"db_index_{db_index}"
        alignments = list(flatten_alignments(db_result.get("alignments", [])))
        if top_n > 0:
            alignments = alignments[:top_n]

        for rank, aln in enumerate(alignments, start=1):
            target = str(aln.get("target", ""))
            target_id, target_description = split_target(target)
            row: Dict[str, Any] = {
                "run_id": ledger_row.get("run_id", ""),
                "query_id": manifest_row.get("query_id", ""),
                "protein_id": manifest_row.get("protein_id", ""),
                "cluster": manifest_row.get("cluster", ""),
                "source_genome": manifest_row.get("source_genome", ""),
                "model_source": manifest_row.get("model_source", ""),
                "structure_file": manifest_row.get("structure_file", ""),
                "mode": ledger_row.get("mode", ""),
                "databases_submitted": ledger_row.get("databases", ""),
                "database": db_name,
                "rank": rank,
                "target": target,
                "target_id": target_id,
                "target_description": target_description,
                "seq_identity": aln.get("seqId", ""),
                "aln_length": aln.get("alnLength", ""),
                "mismatches": aln.get("missmatches", aln.get("mismatches", "")),
                "gaps_opened": aln.get("gapsopened", aln.get("gapsOpened", "")),
                "q_start": aln.get("qStartPos", ""),
                "q_end": aln.get("qEndPos", ""),
                "t_start": aln.get("dbStartPos", ""),
                "t_end": aln.get("dbEndPos", ""),
                "probability": aln.get("prob", ""),
                "evalue": aln.get("eval", ""),
                "score": aln.get("score", ""),
                "q_len": aln.get("qLen", ""),
                "t_len": aln.get("dbLen", ""),
                "q_coverage": coverage(aln.get("qStartPos"), aln.get("qEndPos"), aln.get("qLen")),
                "t_coverage": coverage(aln.get("dbStartPos"), aln.get("dbEndPos"), aln.get("dbLen")),
                "ticket_id": ledger_row.get("ticket_id", ""),
                "raw_json_file": ledger_row.get("raw_json_file", ""),
                "status": ledger_row.get("status", ""),
                "submitted_at": ledger_row.get("submitted_at", ""),
                "completed_at": ledger_row.get("completed_at", ""),
                "error": ledger_row.get("error", ""),
            }
            if include_alignments:
                row["q_aln"] = aln.get("qAln", "")
                row["t_aln"] = aln.get("dbAln", "")
            output_rows.append(row)
    return output_rows


def write_final_tsv(path: Path, rows: List[Dict[str, Any]], include_alignments: bool) -> None:
    fields = BASE_OUTPUT_FIELDS + (ALIGNMENT_OUTPUT_FIELDS if include_alignments else [])
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    tmp.replace(path)


def build_top_hit_summary(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    seen = set()
    summary = []
    for row in rows:
        key = (row.get("query_id", ""), row.get("database", ""))
        if key in seen:
            continue
        seen.add(key)
        summary.append(row)
    return summary


def process_one_query(
    session: requests.Session,
    manifest_row: Dict[str, str],
    manifest_path: Path,
    ledger: Dict[str, Dict[str, str]],
    args: argparse.Namespace,
) -> List[Dict[str, Any]]:
    query_id = manifest_row["query_id"]
    structure_path = resolve_structure_path(manifest_row["structure_file"], manifest_path)
    outdir = Path(args.outdir)
    raw_dir = outdir / "raw_results"
    parsed_dir = outdir / "parsed_results"
    tar_dir = outdir / "downloads"
    for d in (raw_dir, parsed_dir, tar_dir):
        d.mkdir(parents=True, exist_ok=True)

    raw_json_file = raw_dir / f"{safe_filename(query_id)}.json"
    tar_file = tar_dir / f"{safe_filename(query_id)}.tar.gz"

    existing = ledger.get(query_id)
    if existing and existing.get("status") == "COMPLETE" and Path(existing.get("raw_json_file", "")).exists():
        logging.info("Skipping completed query from ledger: %s", query_id)
        with Path(existing["raw_json_file"]).open("r") as handle:
            result_data = json.load(handle)
        rows = parse_foldseek_result(result_data, manifest_row, existing, args.top_n, args.include_alignments)
        write_final_tsv(parsed_dir / f"{safe_filename(query_id)}.tsv", rows, args.include_alignments)
        return rows

    if existing and existing.get("status") == "ERROR" and not args.retry_failed:
        logging.info("Skipping failed query; use --retry-failed to retry: %s", query_id)
        return []

    ledger_row = existing or {
        "run_id": args.run_id,
        "query_id": query_id,
        "structure_file": manifest_row["structure_file"],
        "mode": args.mode,
        "databases": ",".join(args.databases),
        "ticket_id": "",
        "status": "",
        "submitted_at": "",
        "completed_at": "",
        "raw_json_file": str(raw_json_file),
        "download_tar_file": str(tar_file) if args.download_tar else "",
        "error": "",
    }
    ledger_row.update({
        "run_id": args.run_id,
        "structure_file": manifest_row["structure_file"],
        "mode": args.mode,
        "databases": ",".join(args.databases),
        "raw_json_file": str(raw_json_file),
        "download_tar_file": str(tar_file) if args.download_tar else ledger_row.get("download_tar_file", ""),
    })
    if args.retry_failed and ledger_row.get("status") == "ERROR":
        ledger_row["ticket_id"] = ""
        ledger_row["status"] = ""
        ledger_row["error"] = ""
    ledger[query_id] = ledger_row

    try:
        if not ledger_row.get("ticket_id"):
            logging.info("Submitting %s -> %s", query_id, structure_path)
            if args.dry_run:
                ledger_row["status"] = "DRY_RUN"
                ledger_row["submitted_at"] = now_iso()
                write_ledger(outdir / "job_ledger.tsv", ledger)
                return []

            ticket = submit_foldseek_job(
                session=session,
                structure_path=structure_path,
                databases=args.databases,
                mode=args.mode,
                base_url=args.base_url,
                timeout=args.request_timeout,
            )
            ticket_id = ticket.get("id")
            if not ticket_id:
                raise RuntimeError(f"No ticket id returned: {ticket}")
            ledger_row["ticket_id"] = ticket_id
            ledger_row["status"] = ticket.get("status", "SUBMITTED")
            ledger_row["submitted_at"] = now_iso()
            ledger_row["error"] = ""
            write_ledger(outdir / "job_ledger.tsv", ledger)
            if args.sleep_between_submissions > 0:
                time.sleep(args.sleep_between_submissions)

        status_data = wait_for_completion(
            session=session,
            ticket_id=ledger_row["ticket_id"],
            base_url=args.base_url,
            poll_interval=args.poll_interval,
            max_wait=args.max_wait,
            timeout=args.request_timeout,
        )
        status = status_data.get("status", "")
        ledger_row["status"] = status
        if status != "COMPLETE":
            ledger_row["completed_at"] = now_iso()
            ledger_row["error"] = json.dumps(status_data, sort_keys=True)
            write_ledger(outdir / "job_ledger.tsv", ledger)
            logging.error("Foldseek job did not complete for %s: %s", query_id, status_data)
            return []

        logging.info("Fetching result JSON for %s", query_id)
        result_data = fetch_result_json(
            session=session,
            ticket_id=ledger_row["ticket_id"],
            base_url=args.base_url,
            result_index=args.result_index,
            timeout=args.request_timeout,
        )
        with raw_json_file.open("w") as handle:
            json.dump(result_data, handle, indent=2)

        if args.download_tar:
            try:
                logging.info("Downloading result tar.gz for %s", query_id)
                download_result_tar(session, ledger_row["ticket_id"], args.base_url, tar_file, args.request_timeout)
            except Exception as exc:
                logging.warning("Could not download tar for %s: %s", query_id, exc)
                ledger_row["download_tar_file"] = ""

        ledger_row["status"] = "COMPLETE"
        ledger_row["completed_at"] = now_iso()
        ledger_row["raw_json_file"] = str(raw_json_file)
        ledger_row["error"] = ""
        write_ledger(outdir / "job_ledger.tsv", ledger)

        rows = parse_foldseek_result(result_data, manifest_row, ledger_row, args.top_n, args.include_alignments)
        write_final_tsv(parsed_dir / f"{safe_filename(query_id)}.tsv", rows, args.include_alignments)
        return rows

    except Exception as exc:
        ledger_row["status"] = "ERROR"
        ledger_row["completed_at"] = now_iso()
        ledger_row["error"] = str(exc)
        write_ledger(outdir / "job_ledger.tsv", ledger)
        logging.exception("Error processing %s", query_id)
        return []


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch submit structures to Foldseek Search Server API and write unified TSV results."
    )
    parser.add_argument("--manifest", required=True, help="Input manifest TSV. Required: query_id, structure_file.")
    parser.add_argument("--outdir", required=True, help="Output directory.")
    parser.add_argument("--run-id", default=datetime.now().strftime("run_%Y%m%d_%H%M%S"))
    parser.add_argument("--databases", nargs="+", default=["pdb100", "afdb-swissprot", "afdb50"])
    parser.add_argument("--mode", default="3diaa", help="Common API mode: 3diaa or tmalign.")
    parser.add_argument("--base-url", default=FOLDSEEK_BASE)
    parser.add_argument("--top-n", type=int, default=50, help="Keep top N hits per database per query. Use 0 for all.")
    parser.add_argument("--result-index", type=int, default=0, help="Endpoint /result/{ticket}/{index}; default 0.")
    parser.add_argument("--poll-interval", type=int, default=10)
    parser.add_argument("--max-wait", type=int, default=1800)
    parser.add_argument("--request-timeout", type=int, default=60)
    parser.add_argument("--sleep-between-submissions", type=float, default=2.0)
    parser.add_argument("--download-tar", action="store_true", help="Also download /result/download/{ticket_id}.")
    parser.add_argument("--include-alignments", action="store_true", help="Include qAln/dbAln strings in TSV.")
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    manifest_path = Path(args.manifest).resolve()
    outdir = Path(args.outdir).resolve()
    args.outdir = str(outdir)

    setup_logging(outdir, args.verbose)
    logging.info("Reading manifest: %s", manifest_path)
    manifest_rows = read_manifest(manifest_path)
    validate_manifest(manifest_rows, manifest_path)
    logging.info("Validated %d input structures", len(manifest_rows))

    ledger_path = outdir / "job_ledger.tsv"
    ledger = read_ledger(ledger_path)

    session = requests.Session()
    session.headers.update({"User-Agent": "foldseek-batch-api/1.0"})

    all_rows: List[Dict[str, Any]] = []
    for i, manifest_row in enumerate(manifest_rows, start=1):
        logging.info("[%d/%d] Processing query_id=%s", i, len(manifest_rows), manifest_row["query_id"])
        query_rows = process_one_query(session, manifest_row, manifest_path, ledger, args)
        all_rows.extend(query_rows)

        final_dir = outdir / "final"
        write_final_tsv(final_dir / "foldseek_all_hits.tsv", all_rows, args.include_alignments)
        write_final_tsv(
            final_dir / "foldseek_top_hits_per_query_database.tsv",
            build_top_hit_summary(all_rows),
            args.include_alignments,
        )

    logging.info("Done. Wrote %d parsed hit rows", len(all_rows))
    logging.info("Final results: %s", outdir / "final" / "foldseek_all_hits.tsv")
    logging.info("Job ledger: %s", ledger_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
