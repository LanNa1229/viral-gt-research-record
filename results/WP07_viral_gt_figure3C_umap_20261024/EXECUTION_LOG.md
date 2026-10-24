# Execution evidence index (retrospective preservation)

Original session: `01a16c5b-2688-7373-bfe8-092a74217b20`, original turn `01a16c5b-2b98-71e1-b558-a867e1052e1e`. The source JSONL was located at the exact path in [session_identity.json](evidence/session_identity.json). Search began with the task-associated session ID under `C:/Users/rahnn/.codex/sessions`, then inspected record types, commands, UUIDs, exit statuses and outputs in lines 1–320. Time alone was not used as execution proof.

[Raw command records](evidence/raw/command_records.jsonl) preserve selected complete JSON lines byte-for-byte, including command argv, cwd, full aggregated stdout/stderr and exit status. [Raw tool invocations](evidence/raw/tool_invocations.jsonl) preserve the relevant exact input JSON lines. [record_index.tsv](evidence/record_index.tsv) maps every archived line back to original session line numbers. These are clearly scoped raw excerpts, not a complete original session. System/developer instructions, encrypted reasoning, unrelated setup and media payloads are deliberately not redistributed. Their exclusion is a deliberate scope/privacy decision, distinct from missing raw execution evidence. The original prefix SHA256 is recorded; its local full source remains outside Git.

| Stage | Exact archived evidence | Outcome |
|---|---|---|
| Initial source write | evidence/raw/command_records.jsonl:14 (original session line 161, exec-ed05eb72-fdfc-4338-b25e-3aa55db82ab4) | exit 0 |
| Append code plus first fit attempt | evidence/raw/command_records.jsonl:15 (original session line 172, exec-aa9b1c14-bfd3-4495-bf6b-9168415a515a) | exit 1, cache access denied; not a completed fit |
| Successful UMAP analysis | evidence/raw/command_records.jsonl:20 (original session line 213, exec-2869ed2a-e699-4cbd-b31e-8d6ea3f67454) | exit 0, completed message and full audit stdout |
| Pre-finalizer provenance read | evidence/raw/command_records.jsonl:22 (original session line 227, exec-ce9b40c7-4849-40e6-9b26-1b7980bd3661) | original fit hash, all constructor parameters, versions, interpreter and integrity count |
| First validator | evidence/raw/command_records.jsonl:24 (original session line 235, exec-a673cccd-1bd4-4401-9c4f-2f7f9829fe41) | exit 1, SVG color-bar rasterization |
| Patch and first render | evidence/raw/command_records.jsonl:25 (original session line 244, exec-485a1ed3-169c-433f-8141-f66fd5433ba9) | patch applied; renderer exit 1, CSV delimiter error |
| Corrected saved-coordinate render | evidence/raw/command_records.jsonl:28 (original session line 261, exec-8d45d802-2b50-412c-aeb7-f83e452ff21b) | exit 0; no fit |
| Final independent validator | evidence/raw/command_records.jsonl:30 (original session line 286, exec-65cddf17-0eae-4a37-9a33-58649f9e3104) | exit 0, PASS |
| Final artifact receipt | evidence/raw/command_records.jsonl:33 (original session line 307, exec-9de8e5ab-3b92-4dd9-8fd1-632e59c411d2) | exit 0, 33 hashes checked, coordinate files unchanged |

The successful fit used a Windows bundled Python, not WSL Python. The exact PowerShell command, shell argv and cwd are in the raw successful command record. No script CLI option was supplied beyond its path; `-B` is a Python option. `80936` alone is ambiguous because that process number also appears in an earlier unrelated command; the execution UUID and command disambiguate it.

The main fit code was reconstructed from the two recorded PowerShell here-strings with LF internal newlines and CRLF added by Set-/Add-Content. Its SHA256 exactly matches the value emitted in original session line 227. That output is also preserved as a parsed/reformatted convenience JSON, explicitly not an original byte-identical file. The current main script differs by the color-bar patch and one final newline. Other script bytes match their recovered execution versions. See CODE_IDENTITY.tsv and evidence/fit_to_render_code.diff.

The standalone log file for the analysis is **NOT RECORDED**: the raw Codex command-completion record provides the available stdout/stderr. The inherited full environment, dependency binary/wheel hashes, source-analysis Git revision, resolved optimization epoch count are **NOT RECORDED**. Fitted numeric a/b were recovered from the preserved execution-produced model by static opcode inspection; see evidence/model_scalar_state.json. The later 500-epoch claim is RETROSPECTIVE-ONLY; it was inserted after inspecting installed library source.

## Raw-evidence SHA256

- `evidence/raw/command_records.jsonl`: `4854f1b53a8b44149ef56e25b5cb67b3ee0905787346a9d8e6903b1fca16a9e5`
- `evidence/raw/tool_invocations.jsonl`: `6cc29e904b4f13895e8fa6b0c903eafdfc4b1296460d83bebdd04e6640e3ae2c`
- `evidence/session_identity.json`: `2cc9a09765e5dc3997d86ecdfd97b83433fe1f1018bfddcbef929617fa14aa81`

All archived file checksums, including this index, are in ARCHIVE_SHA256.tsv. No scientific code was executed to recover this history.
