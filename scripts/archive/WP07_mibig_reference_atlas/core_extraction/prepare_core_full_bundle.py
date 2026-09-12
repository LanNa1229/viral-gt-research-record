from pathlib import Path
import ast,shutil,tarfile,hashlib
root=Path.cwd();stage=root/"outputs/core_context_full_submission";stage.mkdir()
script=root/"scripts/run_hpc_core_context_full.py"
script.write_text(script.read_text(encoding="utf-8-sig"))
ast.parse(script.read_text())
sb=(root/"scripts/run_core_context_pilot.sbatch").read_text().replace("mibig_core_p20","mibig_core_full").replace("00:20:00","02:00:00").replace("core_context_pilot_20260910","core_context_full_20260910").replace("pilot_%j","full_%j").replace("run_hpc_core_context_pilot.py","run_hpc_core_context_full.py")
(root/"scripts/run_core_context_full.sbatch").write_text(sb)
for name in ("run_hpc_core_context_full.py","extract_core_context_pilot.py","run_core_context_full.sbatch"):
    shutil.copyfile(root/"scripts"/name,stage/name)
a=stage/"saved_membership_audit";a.mkdir()
for name in ("provenance.json","core_position_summary.json","core_position_validation.tsv"):
    shutil.copyfile(root/"outputs/core_context_membership_audit"/name,a/name)
files=sorted(p for p in stage.rglob("*") if p.is_file())
(stage/"SHA256SUMS").write_text("".join(hashlib.sha256(p.read_bytes()).hexdigest()+"  "+p.relative_to(stage).as_posix()+"\n" for p in files))
with tarfile.open(stage/"full_bundle.tar.gz","w:gz") as t:
    for p in files+[stage/"SHA256SUMS"]:t.add(p,arcname=p.relative_to(stage).as_posix())
print((stage/"SHA256SUMS").read_text())

