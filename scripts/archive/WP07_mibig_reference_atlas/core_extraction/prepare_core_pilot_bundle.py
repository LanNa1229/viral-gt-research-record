from pathlib import Path
import shutil,hashlib,tarfile
root=Path.cwd()
stage=root/"outputs/core_context_pilot_submission"
stage.mkdir(exist_ok=False)
for name in ("run_hpc_core_context_pilot.py","extract_core_context_pilot.py","run_core_context_pilot.sbatch"):
    shutil.copyfile(root/"scripts"/name,stage/name)
audit=stage/"saved_membership_audit";audit.mkdir()
for name in ("provenance.json","core_position_summary.json","core_position_validation.tsv"):
    shutil.copyfile(root/"outputs/core_context_membership_audit"/name,audit/name)
files=sorted(p for p in stage.rglob("*") if p.is_file())
(stage/"SHA256SUMS").write_text("".join(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+p.relative_to(stage).as_posix()+"\n" for p in files))
with tarfile.open(stage/"pilot_bundle.tar.gz","w:gz") as t:
    for p in files+[stage/"SHA256SUMS"]:t.add(p,arcname=p.relative_to(stage).as_posix())
print((stage/"SHA256SUMS").read_text())

