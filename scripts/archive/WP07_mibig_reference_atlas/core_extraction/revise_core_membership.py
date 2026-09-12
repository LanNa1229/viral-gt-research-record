from pathlib import Path
p=Path("scripts/extract_core_context_pilot.py")
s=p.read_text(encoding="utf-8-sig").replace('VERSION="1.0.0"','VERSION="2.0.0"')
a=s.index('        hits=[j for j')
b=s.index('    group_errors=[]',a)
s=s[:a]+'''        mask=[p in set(core) for p in td]
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
''' + s[b:]
a=s.index('    summary=dict(');b=s.index('    (out/"core_position_summary.json")',a)
s=s[:a]+'''    summary=dict(rows=len(rows),bgcs=len(groups),structural_failures=sum(not r["structural_pass"] for r in report),
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
''' + s[b:]
s=s.replace('core_start_0based=report[i]["core_start_0based"],core_end_exclusive_0based=report[i]["core_end_exclusive_0based"]) for j,i in enumerate(pilot)])',
'''core_positions_0based=report[i]["core_positions_0based"],
         core_contiguous=report[i]["core_contiguous"]) for j,i in enumerate(pilot)])''')
s=s.replace('vectors=torch.stack([H[j,int(report[i]["core_start_0based"]):int(report[i]["core_end_exclusive_0based"])].mean(dim=0) for j,i in enumerate(indices)])',
'''masks=torch.tensor([[v=="1" for v in report[i]["derived_core_mask"].split()] for i in indices],dtype=torch.bool,device=H.device)
        counts=masks.sum(dim=1)
        expected_counts=torch.tensor([len(rows[i]["sentence"].split()) for i in indices],device=H.device)
        counts_ok=bool(torch.equal(counts,expected_counts))
        if not counts_ok:raise ValueError("Selected contextual position count differs from sentence length")
        vectors=torch.stack([H[j,masks[j]].mean(dim=0) for j in range(len(indices))])''')
s=s.replace('core_pooled_shape=str(tuple(vectors.shape)),encoder_finite=finite,core_pooled_finite=valid,',
'''core_pooled_shape=str(tuple(vectors.shape)),encoder_finite=finite,core_pooled_finite=valid,
                core_length=int(expected_counts[j].item()),selected_position_count=int(counts[j].item()),
                selected_count_matches_core=counts_ok,core_interrupted=report[i]["core_interrupted"],''')
p.write_text(s,encoding="utf-8")

