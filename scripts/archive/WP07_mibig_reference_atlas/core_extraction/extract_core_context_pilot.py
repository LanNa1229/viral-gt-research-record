"""Core-context ablation: validate all rows, then run a 20-row HPC pilot only.
No full-run mode is provided. No UMAP/PCA, training, or package modification.
Adapter contract is documented in README_core_context_pilot.md.
"""
import argparse,csv,hashlib,json,re,sys,importlib.util,datetime
from pathlib import Path
VERSION="2.0.0"
def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def tsv(path,rows,fields=None):
    with open(path,"x",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fields or list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
def read(path,delim):
    with open(path,encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter=delim))
def audit(csv_path,metadata_path,out):
    rows=read(csv_path,",");meta=read(metadata_path,"\t")
    if len(rows)!=12510 or len(meta)!=12510:raise ValueError("Expected 12510 CSV and metadata rows")
    indices=[int(r["embedding_row"]) for r in meta]
    if sorted(indices)!=list(range(12510)):raise ValueError("Invalid embedding_row bijection")
    meta=sorted(meta,key=lambda r:int(r["embedding_row"]))
    if len({r["ID"] for r in rows})!=12510:raise ValueError("Duplicate sample IDs")
    report=[];groups={}
    for i,(r,m) in enumerate(zip(rows,meta)):
        errors=[]
        if any(r[k]!=m[k] for k in ("ID","labels","isBGC")):errors.append("metadata_alignment")
        match=re.fullmatch(r"(BGC[0-9]{7})_(-1|0|1|2|3)",r["ID"])
        if not match:errors.append("invalid_sample_ID")
        else:groups.setdefault(match[1],[]).append((int(match[2]),i))
        core=r["sentence"].split();td=r["TDsentence"].split();raw=r["TDlabels"].split()
        if not (1<=len(core)<=128 and len(td)==len(raw)==128):errors.append("inconsistent_lengths")
        if len(set(core))!=len(core):errors.append("duplicate_core_ID")
        if match and any(not re.fullmatch(re.escape(match[1])+r"_protein_[0-9]+",p) for p in core):errors.append("invalid_core_ID")
        if any(v not in ("0","1") for v in raw):errors.append("nonbinary_TDlabels")
        mask=[p in set(core) for p in td]
        positions=[j for j,v in enumerate(mask) if v]
        occurrences=[td.count(p) for p in core]
        if any(n==0 for n in occurrences):errors.append("missing_core_ID")
        if any(n>1 for n in occurrences):errors.append("duplicate_core_occurrence")
        ordered=[p for p,v in zip(td,mask) if v]==core
        if not ordered:errors.append("core_order_disagreement")
        if sum(mask)!=len(core):errors.append("mask_count_disagreement")
        if sum(mask)>len(core):errors.append("extra_membership_matches")
        contiguous=bool(positions) and positions==list(range(positions[0],positions[-1]+1))
        # Test both binary polarities independently; establish one global polarity.
        one_match=len(raw)==len(mask) and all(v in ("0","1") for v in raw) and mask==[v=="1" for v in raw]
        zero_match=len(raw)==len(mask) and all(v in ("0","1") for v in raw) and mask==[v=="0" for v in raw]
        report.append(dict(embedding_row=i,ID=r["ID"],BGC=match[1] if match else "",
            labels=r["labels"],core_length=len(core),td_length=len(td),tdlabels_length=len(raw),
            each_core_ID_once=all(n==1 for n in occurrences),ordered_core_agrees=ordered,
            selected_position_count=sum(mask),mask_count_agrees=sum(mask)==len(core),
            core_contiguous=contiguous,core_interrupted=bool(positions) and not contiguous,
            core_start_0based=positions[0] if positions else None,
            core_end_exclusive_0based=positions[-1]+1 if positions else None,
            core_positions_0based=" ".join(map(str,positions)),
            derived_core_mask=" ".join(str(int(v)) for v in mask),
            tdlabels_one_is_core_matches=one_match,tdlabels_zero_is_core_matches=zero_match,
            structural_pass=not errors,errors=";".join(errors)))
    candidates=[value for value,field in ((1,"tdlabels_one_is_core_matches"),(0,"tdlabels_zero_is_core_matches"))
                if all(r[field] for r in report)]
    polarity=candidates[0] if len(candidates)==1 else None
    for r in report:
        r["tdlabels_core_value"]=polarity
        r["tdlabels_agree"]=polarity is not None and r["tdlabels_one_is_core_matches" if polarity==1 else "tdlabels_zero_is_core_matches"]
        if not r["tdlabels_agree"]:
            r["errors"]+=";global_TDlabels_polarity_or_mask_disagreement"
            r["structural_pass"]=False
    group_errors=[]
    if len(groups)!=2502:group_errors.append("Expected 2502 BGCs")
    for b,g in groups.items():
        if sorted(s for s,i in g)!=[-1,0,1,2,3]:group_errors.append(b+": augmentation suffixes")
        if len({(rows[i]["sentence"],rows[i]["labels"]) for s,i in g})!=1:group_errors.append(b+": inconsistent core/labels")
    tsv(out/"core_position_validation.tsv",report)
    summary=dict(rows=len(rows),bgcs=len(groups),structural_failures=sum(not r["structural_pass"] for r in report),
        tdlabels_core_value=polarity,tdlabels_noncore_value=1-polarity if polarity is not None else None,
        one_is_core_matching_rows=sum(r["tdlabels_one_is_core_matches"] for r in report),
        zero_is_core_matching_rows=sum(r["tdlabels_zero_is_core_matches"] for r in report),
        mask_agreements=sum(r["tdlabels_agree"] for r in report),
        mask_disagreements=sum(not r["tdlabels_agree"] for r in report),
        contiguous_rows=sum(r["core_contiguous"] for r in report),
        interrupted_rows=sum(r["core_interrupted"] for r in report),
        interrupted_unique_bgcs=len({r["BGC"] for r in report if r["core_interrupted"]}),
        core_length_min=min(r["core_length"] for r in report),core_length_max=max(r["core_length"] for r in report),
        group_errors=group_errors)
    (out/"core_position_summary.json").write_text(json.dumps(summary,indent=2))
    tsv(out/"core_position_summary.tsv",[dict(metric=k,value=v) for k,v in summary.items()])
    if summary["structural_failures"] or group_errors:raise ValueError("Core structure ambiguous; see audit tables")
    # Four distinct biological BGCs, all five augmentations each:
    # BGC0000001, shortest, longest, pure Saccharide (deterministic ties).
    chosen=["BGC0000001"]
    length={b:report[g[0][1]]["core_length"] for b,g in groups.items()}
    for ordering in (sorted(groups,key=lambda b:(length[b],b)),
                     sorted(groups,key=lambda b:(-length[b],b)),
                     sorted(b for b,g in groups.items() if rows[g[0][1]]["labels"].strip()=="Saccharide")):
        chosen.append(next(b for b in ordering if b not in chosen))
    pilot=sorted(i for b in chosen for s,i in groups[b])
    assert len(pilot)==20
    tsv(out/"pilot_rows.tsv",[dict(pilot_embedding_row=j,source_embedding_row=i,ID=rows[i]["ID"],
         labels=rows[i]["labels"],core_length=report[i]["core_length"],
         core_positions_0based=report[i]["core_positions_0based"],
         core_contiguous=report[i]["core_contiguous"]) for j,i in enumerate(pilot)])
    print(json.dumps(summary,indent=2),flush=True)
    return rows,meta,report,pilot
def run_pilot(args,rows,meta,report,pilot):
    # Import only in the validated HPC environment, never in audit-only mode.
    import numpy as np
    import torch
    spec=importlib.util.spec_from_file_location("validated_extraction_adapter",args.adapter)
    adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
    runtime=adapter.build_runtime(args.adapter_config)
    model=runtime["model"];forward=runtime["forward_rows"]
    if model.training or any(m.training for m in model.modules()):raise ValueError("Adapter must return eval-mode model")
    encoder=model.get_submodule(args.encoder_module)
    classifier=model.get_submodule(args.classifier_module)
    baseline=[];behavior=[];pooled=[];allmeans=[]
    # Pilot runs with the exact same row grouping with and without hooks.
    # No seeding, backend, precision, mask, autocast, or inference-context changes here:
    # forward_rows must preserve the validated extraction's execution context.
    batches=[pilot[j:j+args.batch_size] for j in range(0,20,args.batch_size)]
    for indices in batches:
        p=forward(indices)
        if not isinstance(p,torch.Tensor) or p.shape!=(len(indices),7) or not torch.isfinite(p).all():
            raise ValueError("Invalid baseline probabilities")
        baseline.append(p.detach().clone())
    for indices,p0 in zip(batches,baseline):
        captured=[];classifier_inputs=[]
        def encoder_hook(module,inputs,output):
            if not isinstance(output,torch.Tensor):raise ValueError("Encoder output is not a tensor")
            captured.append(output.detach().clone())
            # No return: do not substitute or mutate the forward output.
        def classifier_pre_hook(module,inputs):
            classifier_inputs.append(inputs[0].detach().clone())
        h=encoder.register_forward_hook(encoder_hook)
        h2=classifier.register_forward_pre_hook(classifier_pre_hook)
        try:p1=forward(indices)
        finally:h.remove();h2.remove()
        if len(captured)!=1 or len(classifier_inputs)!=1:raise ValueError("Expected exactly one encoder and classifier call")
        H=captured[0]
        shape_ok=tuple(H.shape)==(len(indices),128,320)
        if not shape_ok:raise ValueError(f"Unexpected encoder shape: {tuple(H.shape)}")
        finite=bool(torch.isfinite(H).all())
        same=bool(torch.equal(p0,p1))
        # Also verify that the captured tensor really feeds the published mean.
        full_mean=H.mean(dim=1).squeeze(1)
        mean_matches=bool(torch.equal(full_mean,classifier_inputs[0]))
        masks=torch.tensor([[v=="1" for v in report[i]["derived_core_mask"].split()] for i in indices],dtype=torch.bool,device=H.device)
        counts=masks.sum(dim=1)
        expected_counts=torch.tensor([len(rows[i]["sentence"].split()) for i in indices],device=H.device)
        counts_ok=bool(torch.equal(counts,expected_counts))
        if not counts_ok:raise ValueError("Selected contextual position count differs from sentence length")
        vectors=torch.stack([H[j,masks[j]].mean(dim=0) for j in range(len(indices))])
        valid=tuple(vectors.shape)==(len(indices),320) and bool(torch.isfinite(vectors).all())
        for j,i in enumerate(indices):
            behavior.append(dict(source_embedding_row=i,ID=rows[i]["ID"],encoder_shape=str(tuple(H.shape)),
                core_pooled_shape=str(tuple(vectors.shape)),encoder_finite=finite,core_pooled_finite=valid,
                core_length=int(expected_counts[j].item()),selected_position_count=int(counts[j].item()),
                selected_count_matches_core=counts_ok,core_interrupted=report[i]["core_interrupted"],
                probabilities_exactly_equal=bool(torch.equal(p0[j],p1[j])),
                max_probability_difference=float((p0[j]-p1[j]).abs().max().item()),
                captured_mean_equals_classifier_input=mean_matches))
        tpath=args.out/"pilot_behavior_validation.tsv"
        if tpath.exists():tpath.unlink() # Only this new run's partial pilot report.
        tsv(tpath,behavior)
        if not(same and finite and valid and mean_matches):raise ValueError("Behavioral pilot FAILED; no embeddings saved")
        pooled.append(vectors.detach().cpu().to(torch.float32).numpy())
        allmeans.append(full_mean.detach().cpu().to(torch.float32).numpy())
    result=np.concatenate(pooled)
    assert result.shape==(20,320) and result.dtype==np.float32 and np.isfinite(result).all()
    # Pilot filenames explicitly distinguish these 20 rows from the full dataset.
    np.save(args.out/"mibig_20_core_context_embeddings.npy",result,allow_pickle=False)
    np.save(args.out/"mibig_20_unmasked_mean_embeddings.npy",np.concatenate(allmeans),allow_pickle=False)
    tsv(args.out/"mibig_20_core_context_metadata.tsv",[dict(embedding_row=j,source_embedding_row=i,
        ID=meta[i]["ID"],labels=meta[i]["labels"],isBGC=meta[i]["isBGC"]) for j,i in enumerate(pilot)])
    # Read-only reference comparison: diagnoses differences in HPC call context.
    if args.reference_embeddings:
        ref=np.load(args.reference_embeddings,allow_pickle=False)
        assert ref.shape==(12510,320)
        observed=np.concatenate(allmeans)
        tsv(args.out/"pilot_existing_representation_comparison.tsv",[dict(source_embedding_row=i,
            ID=rows[i]["ID"],exactly_equal=bool(np.array_equal(observed[j],ref[i])),
            max_abs_difference=float(np.max(np.abs(observed[j]-ref[i])))) for j,i in enumerate(pilot)])
    (args.out/"pilot_runtime.json").write_text(json.dumps(dict(torch_version=torch.__version__,
        numpy_version=np.__version__,adapter_sha256=sha(args.adapter),runtime_provenance=runtime["provenance"],
        encoder_module=args.encoder_module,classifier_module=args.classifier_module,full_context_saved=False,
        pilot_pass=True,rows=20),indent=2,default=str))
    print("20-row behavioral pilot PASS. STOP: no full-run command or full dataset output created.",flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv",type=Path,required=True);p.add_argument("--metadata",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True);p.add_argument("--mode",choices=["audit","pilot"],default="audit")
    p.add_argument("--adapter",type=Path);p.add_argument("--adapter-config",type=Path)
    p.add_argument("--encoder-module",default="encoder");p.add_argument("--classifier-module",default="classifier")
    p.add_argument("--batch-size",type=int,default=5)
    p.add_argument("--reference-embeddings",type=Path)
    args=p.parse_args()
    if args.out.exists():raise FileExistsError("Output directory must be new; existing data are never overwritten")
    if args.mode=="pilot" and (not args.adapter or not args.adapter_config):p.error("Pilot requires validated HPC adapter and configuration")
    if args.batch_size<1:p.error("Batch size must be positive")
    args.out.mkdir(parents=True)
    hashes={str(path):sha(path) for path in (args.csv,args.metadata)}
    try:
        rows,meta,report,pilot=audit(args.csv,args.metadata,args.out)
        if args.mode=="pilot":run_pilot(args,rows,meta,report,pilot)
    finally:
        unchanged=all(sha(Path(path))==h for path,h in hashes.items())
        (args.out/"provenance.json").write_text(json.dumps(dict(script_version=VERSION,script_sha256=sha(Path(__file__)),
           utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),command=sys.argv,input_sha256=hashes,
           inputs_unchanged=unchanged,mode=args.mode),indent=2))
        if not unchanged:raise RuntimeError("Source inputs changed")

