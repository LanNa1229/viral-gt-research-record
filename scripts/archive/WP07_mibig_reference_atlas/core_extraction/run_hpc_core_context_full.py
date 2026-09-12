"""Full core-context extraction: original classify(), streaming read-only hooks."""
import argparse,sys,json,inspect,datetime,os
from pathlib import Path
from argparse import Namespace
import numpy as np
import torch
from bgc_prophet.command.classify import geneClassifier
from extract_core_context_pilot import read,tsv,sha

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--work",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args();w=args.work;out=args.out
    out.mkdir(parents=True,exist_ok=False)
    train=w/"reference_data/BGC_train_dataset_classify.csv"
    checkpoint=w/"models/model/classifier.pt";lmdb=w/"reference_data/lmdb_train"
    reference_dir=w/"reference_data/embedding_mibig_full"
    fullmeta=reference_dir/"mibig_12510_metadata.tsv"
    reference_path=reference_dir/"mibig_12510_embeddings.npy"
    source=Path(inspect.getfile(geneClassifier));package=source.parents[1]
    saved=Path(__file__).parent/"saved_membership_audit"
    paths=[train,checkpoint,lmdb/"data.mdb",fullmeta,reference_path,w/"extract_bgc_embeddings.py",
           *sorted(package.rglob("*.py"))]
    print("Hashing protected inputs and installed source",flush=True)
    before={str(p):sha(p) for p in paths}
    tsv(out/"input_integrity_before.tsv",[dict(path=p,sha256=h) for p,h in before.items()])
    previous=read(w/"core_context_pilot_20260910/results/input_integrity_before.tsv","\t")
    if any(before[r["path"]]!=r["sha256"] for r in previous if r["path"] in before):
        raise ValueError("Protected input differs from successful pilot")
    ap=json.loads((saved/"provenance.json").read_text())
    aq=json.loads((saved/"core_position_summary.json").read_text())
    assert ap["inputs_unchanged"] and aq["structural_failures"]==0 and aq["mask_agreements"]==12510
    assert aq["tdlabels_core_value"]==1 and aq["interrupted_rows"]==1900
    assert before[str(train)]==ap["input_sha256"]["BGC_train_dataset_classify.csv"]
    assert before[str(fullmeta)]==ap["input_sha256"]["analysis/bgcprophet_reference/data/mibig_12510_metadata.tsv"]
    rows=read(train,",");meta=read(fullmeta,"\t");report=read(saved/"core_position_validation.tsv","\t")
    assert len(rows)==len(meta)==len(report)==12510
    for i,(r,m,a) in enumerate(zip(rows,meta,report)):
        assert int(m["embedding_row"])==int(a["embedding_row"])==i
        assert all(m[k]==r[k] for k in ("ID","labels","isBGC")) and a["ID"]==r["ID"]
        assert a["structural_pass"]==a["tdlabels_agree"]=="True"
    masks=np.array([[v=="1" for v in a["derived_core_mask"].split()] for a in report],dtype=bool)
    lengths=np.array([len(r["sentence"].split()) for r in rows])
    assert masks.shape==(12510,128) and np.array_equal(masks.sum(1),lengths)
    reference=np.load(reference_path,allow_pickle=False)
    assert reference.shape==(12510,320) and np.isfinite(reference).all()
    (out/"reused_audit_provenance.json").write_text(json.dumps(dict(audit=ap,summary=aq,
        validation_sha256=sha(saved/"core_position_validation.tsv"),structural_audit_recomputed=False),indent=2))
    clf=geneClassifier(Namespace(datasetPath=train,classifierPath=checkpoint,lmdbPath=lmdb,
        outputPath=out,device="cpu",batch_size=128,name="mibig_12510_core_context",classify_t=.5))
    assert not any(m.training for m in clf.model.modules())
    state={k:v.detach().clone() for k,v in clf.model.state_dict().items()}
    pooled=np.empty((12510,320),dtype=np.float32);checks=[];pending=[];offset=0;calls=0
    # No forward changes. Pool each detached batch and retain no contextual tensor.
    def capture(module,inputs,output):
        nonlocal offset,calls
        if pending:raise ValueError("Unconsumed encoder output")
        if not isinstance(output,torch.Tensor):raise ValueError("Non-tensor encoder output")
        H=output.detach().clone();n=H.shape[0];start=offset;end=start+n
        if tuple(H.shape)!=(n,128,320) or end>12510:raise ValueError("Encoder shape/row mismatch")
        mask=torch.as_tensor(masks[start:end],device=H.device,dtype=torch.bool)
        vectors=torch.stack([H[j,mask[j]].mean(dim=0) for j in range(n)])
        whole=H.mean(dim=1).squeeze(1)
        assert vectors.shape==(n,320) and torch.isfinite(H).all() and torch.isfinite(vectors).all()
        whole_np=whole.cpu().numpy().astype(np.float32)
        pooled[start:end]=vectors.cpu().numpy().astype(np.float32)
        for j,i in enumerate(range(start,end)):
            count=int(mask[j].sum());assert count==int(lengths[i])
            delta=float(np.max(np.abs(whole_np[j]-reference[i])))
            agrees=bool(np.allclose(whole_np[j],reference[i],rtol=1e-5,atol=1e-6))
            checks.append(dict(embedding_row=i,ID=rows[i]["ID"],core_length=int(lengths[i]),
                selected_positions=count,selected_count_matches_core=True,
                core_interrupted=report[i]["core_interrupted"]=="True",tdlabels_agree=True,
                encoder_shape=str(tuple(H.shape)),pooled_shape=str(tuple(vectors.shape)),
                encoder_finite=True,pooled_finite=True,metadata_aligned=True,
                unmasked_reference_exact=np.array_equal(whole_np[j],reference[i]),
                unmasked_reference_close=agrees,unmasked_reference_max_abs_difference=delta,
                unmasked_mean_matches_classifier_input=False))
        pending.append((start,end,whole));offset=end;calls+=1
        if calls%10==0:print(f"Captured {offset}/12510 rows",flush=True)
    def capture_classifier(module,inputs):
        if len(pending)!=1:raise ValueError("Classifier/encoder call mismatch")
        start,end,whole=pending.pop()
        if not torch.equal(whole,inputs[0].detach()):raise ValueError("Encoder mean differs from classifier input")
        for i in range(start,end):checks[i]["unmasked_mean_matches_classifier_input"]=True
    h=clf.model.encoder.register_forward_hook(capture)
    h2=clf.model.classifier.register_forward_pre_hook(capture_classifier)
    print("Single original classify() over 12510 rows",flush=True)
    try:clf.classify()
    finally:h.remove();h2.remove()
    tsv(out/"core_context_full_validation.tsv",checks)
    valid=(offset==12510 and len(checks)==12510 and not pending and pooled.shape==(12510,320)
        and pooled.dtype==np.float32 and np.isfinite(pooled).all()
        and sum(r["core_interrupted"] for r in checks)==1900
        and all(r["unmasked_reference_close"] and r["unmasked_mean_matches_classifier_input"] for r in checks)
        and all(torch.equal(v,clf.model.state_dict()[k]) for k,v in state.items())
        and clf.results.shape==(12510,7) and np.isfinite(clf.results).all())
    print("Checking protected inputs after extraction",flush=True)
    integrity=[dict(path=p,sha256_before=h,sha256_after=sha(Path(p))) for p,h in before.items()]
    for r in integrity:r["unchanged"]=r["sha256_before"]==r["sha256_after"]
    tsv(out/"input_integrity_after.tsv",integrity)
    valid=valid and all(r["unchanged"] for r in integrity)
    summary=dict(passed=bool(valid),rows=offset,shape=list(pooled.shape),dtype=str(pooled.dtype),
        all_finite=bool(np.isfinite(pooled).all()),interrupted_rows=sum(r["core_interrupted"] for r in checks),
        unmasked_reference_exact_rows=sum(bool(r["unmasked_reference_exact"]) for r in checks),
        unmasked_reference_max_abs_difference=max(r["unmasked_reference_max_abs_difference"] for r in checks),
        reference_rtol=1e-5,reference_atol=1e-6,all_protected_inputs_unchanged=all(r["unchanged"] for r in integrity),
        encoder_batches=calls,baseline_inference_runs=0,original_classify_runs=1,full_tensor_saved=False)
    (out/"full_summary.json").write_text(json.dumps(summary,indent=2))
    tsv(out/"full_summary.tsv",[dict(metric=k,value=v) for k,v in summary.items()])
    if not valid:raise RuntimeError("Full extraction failed validation; embeddings not released")
    np.save(out/"mibig_12510_core_context_embeddings.npy",pooled,allow_pickle=False)
    tsv(out/"mibig_12510_core_context_metadata.tsv",meta)
    import importlib.metadata
    (out/"runtime_provenance.json").write_text(json.dumps(dict(script_sha256=sha(Path(__file__)),
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),python=sys.version,
        packages={k:importlib.metadata.version(k) for k in ("torch","numpy","pandas","lmdb")},
        job_id=os.environ.get("SLURM_JOB_ID"),command=sys.argv,device="cpu",batch_size=128,
        model_forward_reimplemented=False,individual_vectors_normalized=False,
        pooling="H[j, derived_membership_mask[j]].mean(dim=0)",random_settings="Unchanged from validated original loader",
        validated_pilot_job="48067092"),indent=2))
    print(json.dumps(summary,indent=2),flush=True)
if __name__=="__main__":main()

