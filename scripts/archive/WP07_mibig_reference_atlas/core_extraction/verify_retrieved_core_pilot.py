"""Independent validation of retrieved pilot outputs; no model inference."""
import csv,json,hashlib,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"outputs/core_context_pilot_20260910";R=BASE/"results"
def read(path,delimiter="\t"):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter=delimiter))
checks=[]
def check(name,ok):
    checks.append(dict(check=name,passed=bool(ok)))
    if not ok:raise ValueError(name)
for line in (BASE/"OUTPUT_SHA256SUMS").read_text().splitlines():
    h,name=line.split("  ",1)
    check("sha256:"+name,hashlib.sha256((BASE/name).read_bytes()).hexdigest()==h)
x=np.load(R/"mibig_20_core_context_embeddings.npy",allow_pickle=False)
p=np.load(R/"mibig_20_probabilities_hooked.npy",allow_pickle=False)
q=np.load(R/"mibig_20_probabilities_unhooked.npy",allow_pickle=False)
check("embedding_shape_dtype_finite",x.shape==(20,320) and x.dtype==np.float32 and np.isfinite(x).all())
check("probabilities_equal_finite",p.shape==q.shape==(20,7) and np.array_equal(p,q) and np.isfinite(p).all())
meta=read(R/"mibig_20_core_context_metadata.tsv");manifest=read(R/"pilot_manifest.tsv")
validation=read(R/"pilot_behavior_validation.tsv")
source=read(ROOT/"BGC_train_dataset_classify.csv",",")
check("row_count_and_order",len(meta)==len(manifest)==len(validation)==20 and [int(r["embedding_row"]) for r in meta]==list(range(20)))
for j,r in enumerate(meta):
    i=int(r["source_embedding_row"]);original=source[i]
    check("identity:"+r["ID"],all(r[k]==original[k] for k in ("ID","labels","isBGC")) and manifest[j]["ID"]==validation[j]["ID"]==r["ID"])
    core=original["sentence"].split();td=original["TDsentence"].split()
    indices=[k for k,v in enumerate(td) if v in set(core)]
    check("pooled_positions:"+r["ID"],[td[k] for k in indices]==core and len(indices)==int(validation[j]["selected_positions"])==len(core)
        and indices==[int(k) for k in manifest[j]["core_positions_0based"].split()])
    check("behavior:"+r["ID"],all(validation[j][k]=="True" for k in ("selected_count_matches_core","tdlabels_agree","encoder_finite","pooled_finite","unmasked_mean_matches_classifier_input","probabilities_exactly_identical"))
          and validation[j]["encoder_shape"]=="(20, 128, 320)" and validation[j]["pooled_shape"]=="(20, 320)")
check("original_file_integrity",all(r["unchanged"]=="True" and r["sha256_before"]==r["sha256_after"] for r in read(R/"input_integrity_after.tsv")))
comparison=read(R/"existing_full_run_comparison.tsv")
maximum=max(float(r["mean_max_abs_difference"]) for r in comparison)
probmax=max(float(r["probability_max_abs_difference"]) for r in comparison)
with (BASE/"local_verification.tsv").open("x",newline="") as f:
    w=csv.DictWriter(f,["check","passed"],delimiter="\t");w.writeheader();w.writerows(checks)
summary=dict(local_checks=len(checks),passed=True,existing_full_mean_max_abs_difference=maximum,
    existing_full_probability_max_abs_difference=probmax,interrupted_pilot_rows=sum(r["core_interrupted"]=="True" for r in validation))
(BASE/"local_verification_summary.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))

