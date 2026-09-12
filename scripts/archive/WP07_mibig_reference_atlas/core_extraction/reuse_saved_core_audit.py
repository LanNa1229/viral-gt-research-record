from pathlib import Path
p=Path("scripts/run_hpc_core_context_pilot.py")
s=p.read_text(encoding="utf-8-sig")
old='''    auditout=out/"all_rows_audit";auditout.mkdir()
    rows,meta,report,_=audit(train,fullmeta,auditout)
    (auditout/"pilot_rows.tsv").rename(auditout/"coverage_candidate_rows_not_executed.tsv")'''
new='''    # Reuse the completed audit, bound to the exact original input checksums.
    saved=Path(__file__).parent/"saved_membership_audit"
    ap=json.loads((saved/"provenance.json").read_text())
    summary=json.loads((saved/"core_position_summary.json").read_text())
    assert ap["inputs_unchanged"] and summary["structural_failures"]==0 and summary["mask_agreements"]==12510
    assert sha(train)==ap["input_sha256"]["BGC_train_dataset_classify.csv"]
    assert sha(fullmeta)==ap["input_sha256"]["analysis/bgcprophet_reference/data/mibig_12510_metadata.tsv"]
    rows=read(train,",");report=read(saved/"core_position_validation.tsv","\\t")
    assert len(rows)==len(report)==12510
    for i,r in enumerate(report):
        assert int(r["embedding_row"])==i and r["ID"]==rows[i]["ID"] and r["structural_pass"]=="True"
        for key in ("tdlabels_agree","core_interrupted","core_contiguous"):
            r[key]=r[key]=="True"
    (out/"reused_audit_provenance.json").write_text(json.dumps(dict(
        source=str(saved),audit_provenance=ap,summary=summary,
        validation_sha256=sha(saved/"core_position_validation.tsv"),structural_audit_recomputed=False),indent=2))'''
assert old in s
p.write_text(s.replace(old,new),encoding="utf-8")
import ast
ast.parse(p.read_text())

