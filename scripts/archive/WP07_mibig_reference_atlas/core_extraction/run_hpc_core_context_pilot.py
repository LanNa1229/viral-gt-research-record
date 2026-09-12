"""HPC-only 20-row core-context pilot, using original geneClassifier.classify."""
import argparse,sys,json,inspect,datetime,os
from pathlib import Path
from argparse import Namespace
import numpy as np
import torch
from bgc_prophet.command.classify import geneClassifier
from extract_core_context_pilot import audit,read,tsv,sha

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--work",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();w=a.work;out=a.out
    out.mkdir(parents=True,exist_ok=False)
    train=w/"reference_data/BGC_train_dataset_classify.csv"
    pilot=w/"reference_data/BGC_train_dataset_classify_pilot20.csv"
    checkpoint=w/"models/model/classifier.pt";lmdb=w/"reference_data/lmdb_train"
    fullmeta=w/"reference_data/embedding_mibig_full/mibig_12510_metadata.tsv"
    source=Path(inspect.getfile(geneClassifier))
    paths=[train,pilot,checkpoint,lmdb/"data.mdb",fullmeta,w/"extract_bgc_embeddings.py",
           source,source.parents[1]/"train/classifier.py",source.parents[1]/"train/data.py"]
    print("Hashing original inputs before pilot",flush=True)
    before={str(q):sha(q) for q in paths}
    tsv(out/"input_integrity_before.tsv",[dict(path=q,sha256=h) for q,h in before.items()])
    # Reuse the completed audit, bound to the exact original input checksums.
    saved=Path(__file__).parent/"saved_membership_audit"
    ap=json.loads((saved/"provenance.json").read_text())
    summary=json.loads((saved/"core_position_summary.json").read_text())
    assert ap["inputs_unchanged"] and summary["structural_failures"]==0 and summary["mask_agreements"]==12510
    assert sha(train)==ap["input_sha256"]["BGC_train_dataset_classify.csv"]
    assert sha(fullmeta)==ap["input_sha256"]["analysis/bgcprophet_reference/data/mibig_12510_metadata.tsv"]
    rows=read(train,",");report=read(saved/"core_position_validation.tsv","\t")
    assert len(rows)==len(report)==12510
    for i,r in enumerate(report):
        assert int(r["embedding_row"])==i and r["ID"]==rows[i]["ID"] and r["structural_pass"]=="True"
        for key in ("tdlabels_agree","core_interrupted","core_contiguous"):
            r[key]=r[key]=="True"
    (out/"reused_audit_provenance.json").write_text(json.dumps(dict(
        source=str(saved),audit_provenance=ap,summary=summary,
        validation_sha256=sha(saved/"core_position_validation.tsv"),structural_audit_recomputed=False),indent=2))
    small=read(pilot,",");lookup={r["ID"]:i for i,r in enumerate(rows)}
    if len(small)!=20 or len({r["ID"] for r in small})!=20:raise ValueError("Pilot CSV must contain exactly 20 distinct rows")
    indices=[]
    for r in small:
        i=lookup[r["ID"]]
        if r!=rows[i]:raise ValueError("Pilot CSV row differs from full original CSV: "+r["ID"])
        indices.append(i)
    masks=np.array([[v=="1" for v in report[i]["derived_core_mask"].split()] for i in indices],dtype=bool)
    lengths=np.array([len(r["sentence"].split()) for r in small])
    assert masks.shape==(20,128) and np.array_equal(masks.sum(1),lengths)
    tsv(out/"pilot_manifest.tsv",[dict(embedding_row=j,source_embedding_row=i,ID=small[j]["ID"],
        core_length=int(lengths[j]),core_positions_0based=report[i]["core_positions_0based"],
        core_interrupted=report[i]["core_interrupted"]) for j,i in enumerate(indices)])
    args=Namespace(datasetPath=pilot,classifierPath=checkpoint,lmdbPath=lmdb,outputPath=out,
                   device="cpu",batch_size=128,name="mibig_20_core_context",classify_t=.5)
    clf=geneClassifier(args)
    assert not any(m.training for m in clf.model.modules())
    state={k:v.detach().clone() for k,v in clf.model.state_dict().items()}
    print("Unhooked original classify()",flush=True)
    clf.classify();baseline=clf.results.copy()
    captures=[];classifier_inputs=[]
    def capture(module,inputs,output):
        if not isinstance(output,torch.Tensor):raise ValueError("Non-tensor encoder output")
        captures.append(output.detach().clone())
    def capture_mean(module,inputs):
        classifier_inputs.append(inputs[0].detach().clone())
    h=clf.model.encoder.register_forward_hook(capture)
    h2=clf.model.classifier.register_forward_pre_hook(capture_mean)
    print("Hooked original classify()",flush=True)
    try:clf.classify()
    finally:h.remove();h2.remove()
    hooked=clf.results.copy()
    if len(captures)!=len(classifier_inputs):raise ValueError("Hook call mismatch")
    vectors=[];unmasked=[];checks=[];offset=0
    for H,cin in zip(captures,classifier_inputs):
        n=H.shape[0]
        if tuple(H.shape)!=(n,128,320):raise ValueError("Encoder shape mismatch")
        mm=torch.as_tensor(masks[offset:offset+n],device=H.device,dtype=torch.bool)
        pooled=torch.stack([H[j,mm[j]].mean(dim=0) for j in range(n)])
        whole=H.mean(dim=1).squeeze(1)
        for j in range(n):
            k=offset+j;i=indices[k]
            checks.append(dict(embedding_row=k,source_embedding_row=i,ID=small[k]["ID"],
                encoder_shape=str(tuple(H.shape)),pooled_shape=str(tuple(pooled.shape)),
                core_length=int(lengths[k]),selected_positions=int(mm[j].sum()),
                selected_count_matches_core=int(mm[j].sum())==int(lengths[k]),
                tdlabels_agree=report[i]["tdlabels_agree"],core_interrupted=report[i]["core_interrupted"],
                encoder_finite=bool(torch.isfinite(H[j]).all()),pooled_finite=bool(torch.isfinite(pooled[j]).all()),
                unmasked_mean_matches_classifier_input=bool(torch.equal(whole[j],cin[j])),
                probabilities_exactly_identical=bool(np.array_equal(baseline[k],hooked[k])),
                max_probability_difference=float(np.max(np.abs(baseline[k]-hooked[k])))))
        vectors.append(pooled.cpu().numpy().astype(np.float32))
        unmasked.append(whole.cpu().numpy().astype(np.float32));offset+=n
    tsv(out/"pilot_behavior_validation.tsv",checks)
    valid_fields=("selected_count_matches_core","tdlabels_agree","encoder_finite","pooled_finite",
                  "unmasked_mean_matches_classifier_input","probabilities_exactly_identical")
    passed=(offset==20 and baseline.shape==hooked.shape==(20,7) and np.isfinite(baseline).all()
        and np.isfinite(hooked).all() and all(all(r[f] for f in valid_fields) for r in checks)
        and all(torch.equal(v,clf.model.state_dict()[k]) for k,v in state.items()))
    print("Hashing original inputs after pilot",flush=True)
    integrity=[dict(path=q,sha256_before=h,sha256_after=sha(Path(q))) for q,h in before.items()]
    for r in integrity:r["unchanged"]=r["sha256_before"]==r["sha256_after"]
    tsv(out/"input_integrity_after.tsv",integrity)
    passed=passed and all(r["unchanged"] for r in integrity)
    summary=dict(pilot_pass=bool(passed),samples=20,encoder_shapes=[list(h.shape) for h in captures],
        core_embedding_shape=[20,320],max_probability_difference=float(np.max(np.abs(baseline-hooked))),
        probabilities_exactly_identical=bool(np.array_equal(baseline,hooked)),
        all_original_inputs_unchanged=all(r["unchanged"] for r in integrity),
        interrupted_pilot_rows=sum(r["core_interrupted"] for r in checks))
    (out/"pilot_summary.json").write_text(json.dumps(summary,indent=2))
    if not passed:raise RuntimeError("Pilot failed: no core embeddings released")
    core=np.concatenate(vectors);assert core.shape==(20,320) and core.dtype==np.float32
    np.save(out/"mibig_20_core_context_embeddings.npy",core,allow_pickle=False)
    np.save(out/"mibig_20_unmasked_mean_embeddings.npy",np.concatenate(unmasked),allow_pickle=False)
    np.save(out/"mibig_20_probabilities_unhooked.npy",baseline,allow_pickle=False)
    np.save(out/"mibig_20_probabilities_hooked.npy",hooked,allow_pickle=False)
    tsv(out/"mibig_20_core_context_metadata.tsv",[dict(embedding_row=j,source_embedding_row=i,
        ID=small[j]["ID"],labels=small[j]["labels"],isBGC=small[j]["isBGC"]) for j,i in enumerate(indices)])
    # Additional comparison to existing full-run outputs; not the hook equivalence test.
    original=np.load(w/"reference_data/embedding_mibig_full/mibig_12510_embeddings.npy",allow_pickle=False)
    oldp=np.load(w/"reference_data/embedding_mibig_full/mibig_12510_probabilities.npy",allow_pickle=False)
    tsv(out/"existing_full_run_comparison.tsv",[dict(embedding_row=j,source_embedding_row=i,
        mean_max_abs_difference=float(np.max(np.abs(np.concatenate(unmasked)[j]-original[i]))),
        probability_max_abs_difference=float(np.max(np.abs(baseline[j]-oldp[i])))) for j,i in enumerate(indices)])
    import importlib.metadata
    (out/"runtime_provenance.json").write_text(json.dumps(dict(script_sha256=sha(Path(__file__)),
        audit_script_sha256=sha(Path(__file__).with_name("extract_core_context_pilot.py")),
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,
        packages={k:importlib.metadata.version(k) for k in ("torch","numpy","pandas","lmdb")},
        command=sys.argv,job_id=os.environ.get("SLURM_JOB_ID"),device="cpu",batch_size=128,
        original_classify_used=True,full_tensor_saved=False,full_run_launched=False),indent=2))
    print(json.dumps(summary,indent=2),flush=True)
if __name__=="__main__":main()

