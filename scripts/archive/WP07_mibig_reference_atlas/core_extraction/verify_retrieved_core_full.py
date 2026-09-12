"""Read-only numerical verification of retrieved full core-context outputs."""
import csv,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"outputs/core_context_full_20260910";R=BASE/"results_validatednode"
def read(p,d="\t"):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter=d))
for line in (BASE/"OUTPUT_SHA256SUMS").read_text().splitlines():
    digest,name=line.split("  ",1)
    assert hashlib.sha256((BASE/name).read_bytes()).hexdigest()==digest,name
x=np.load(R/"mibig_12510_core_context_embeddings.npy",allow_pickle=False)
assert x.shape==(12510,320) and x.dtype==np.float32 and np.isfinite(x).all()
meta=read(R/"mibig_12510_core_context_metadata.tsv")
original_meta=read(ROOT/"analysis/bgcprophet_reference/data/mibig_12510_metadata.tsv")
validation=read(R/"core_context_full_validation.tsv")
source=read(ROOT/"BGC_train_dataset_classify.csv",",")
assert meta==original_meta and len(validation)==len(source)==12510
checks=[]
for i,(m,v,s) in enumerate(zip(meta,validation,source)):
    assert int(m["embedding_row"])==int(v["embedding_row"])==i and m["ID"]==v["ID"]==s["ID"]
    assert all(m[k]==s[k] for k in ("labels","isBGC"))
    core=s["sentence"].split();td=s["TDsentence"].split()
    positions=[j for j,p in enumerate(td) if p in set(core)]
    assert [td[j] for j in positions]==core
    assert len(positions)==len(core)==int(v["core_length"])==int(v["selected_positions"])
    interrupted=positions!=list(range(positions[0],positions[-1]+1))
    assert interrupted==(v["core_interrupted"]=="True")
    for k in ("selected_count_matches_core","tdlabels_agree","encoder_finite","pooled_finite",
              "metadata_aligned","unmasked_reference_close","unmasked_mean_matches_classifier_input"):
        assert v[k]=="True",(i,k)
    checks.append(dict(embedding_row=i,ID=m["ID"],passed=True,core_length=len(core),
                       selected_positions=len(positions),core_interrupted=interrupted))
assert sum(r["core_interrupted"] for r in checks)==1900
integrity=read(R/"input_integrity_after.tsv")
assert all(r["unchanged"]=="True" and r["sha256_before"]==r["sha256_after"] for r in integrity)
for path in (ROOT/"BGC_train_dataset_classify.csv",ROOT/"analysis/bgcprophet_reference/data/mibig_12510_metadata.tsv"):
    name=path.name
    matching=[r for r in integrity if Path(r["path"]).name==name]
    assert len(matching)==1 and hashlib.sha256(path.read_bytes()).hexdigest()==matching[0]["sha256_before"]
pilotdir=ROOT/"outputs/core_context_pilot_20260910/results"
pm=read(pilotdir/"mibig_20_core_context_metadata.tsv")
px=np.load(pilotdir/"mibig_20_core_context_embeddings.npy",allow_pickle=False)
pilot_delta=float(np.max(np.abs(x[[int(r["source_embedding_row"]) for r in pm]]-px)))
summary=dict(passed=True,rows=12510,shape=list(x.shape),dtype=str(x.dtype),all_finite=True,
    metadata_exact=True,interrupted_rows=1900,protected_files_unchanged=len(integrity),
    unmasked_reference_exact_rows=sum(r["unmasked_reference_exact"]=="True" for r in validation),
    unmasked_reference_max_abs_difference=max(float(r["unmasked_reference_max_abs_difference"]) for r in validation),
    successful_pilot_core_max_abs_difference=pilot_delta)
with (BASE/"local_validation.tsv").open("x",newline="") as f:
    w=csv.DictWriter(f,list(checks[0]),delimiter="\t");w.writeheader();w.writerows(checks)
(BASE/"local_validation_summary.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))

