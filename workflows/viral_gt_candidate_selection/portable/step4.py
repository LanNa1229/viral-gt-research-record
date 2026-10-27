#!/usr/bin/env python3
"""Step 4 v1.0.0: selection from frozen public tables; Python standard library only."""
import argparse, collections, csv, hashlib, json, re, sys
from pathlib import Path
VERSION="1.0.0"
EXCLUDE={"KC511810.1","KX960781.1","KY189087.1"}
EXPECTED={"positive_windows15":331,"positive_windows11":344,"positive_occurrences":388,"positive_loci":288,"genome_units":466,"candidates":245,"negatives":215,"unresolved":6,"Both":190,"11-only":29,"15-only":26,"sankey":242}
def sha(p):
    with Path(p).open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()
def read(p):
    with Path(p).open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write(p,rows):
    rows=list(rows);p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
def yes(v):return str(v).lower() in {"true","1","yes"}
def ids(value):return value.split(";") if value else []
def unique(rows,key):
    d={r[key]:r for r in rows}
    if len(d)!=len(rows):raise ValueError("Duplicate "+key)
    return d
def genome_key(w):
    h=w.get("genome_sequence_sha256","").strip();t=w.get("taxid","").strip()
    return ("sequence_and_taxid",h,t) if h and t else ("accession_specific",w["genome_accession"],h+":"+t)
def stable_id(prefix,value):
    return prefix+hashlib.sha256(json.dumps(value,separators=(",",":")).encode()).hexdigest()[:24]
def genome_id(w):return stable_id("vgg_",genome_key(w))
def combine_segments(segments):
    out=[]
    for a,b in sorted(segments):
        if not out or a>out[-1][1]:out.append([a,b])
        else:out[-1][1]=max(out[-1][1],b)
    return out
def window_segments(w,genes):
    regular=[];wrapped=[]
    for g in genes:
        parts=json.loads(g.get("location_parts","[]"))
        crossing=len(parts)>1 and ((parts[0][2]==1 and parts[0][0]>parts[-1][0]) or (parts[0][2]==-1 and parts[0][0]<parts[-1][0]))
        if crossing:wrapped.extend((int(p[0]),int(p[1])) for p in parts)
        elif g.get("start_1based"):regular.append((int(g["start_1based"]),int(g["end_1based"])))
    if not wrapped:return [[int(w["start_1based"]),int(w["end_1based"])]]
    if regular:wrapped.append((min(p[0] for p in regular),max(p[1] for p in regular)))
    return combine_segments(wrapped)
def reconstruct_groups(windows,genes):
    bygenes=collections.defaultdict(list);bygen=collections.defaultdict(list)
    for g in genes:bygenes[g["window_id"]].append(g)
    segments={}
    for wid,w in windows.items():
        gg=bygenes[wid]
        if len(gg)!=int(w["n_genes"]) or len({g["gene_index"] for g in gg})!=len(gg):raise ValueError("Gene membership "+wid)
        segments[wid]=window_segments(w,gg);bygen[genome_id(w)].append(wid)
    loci={};locus_for={}
    for gid,wids in bygen.items():
        parent=list(range(len(wids)))
        def find(i):
            while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
            return i
        for i,left in enumerate(wids):
            for j in range(i):
                right=wids[j]
                if any(a<=d and c<=b for a,b in segments[left] for c,d in segments[right]):parent[find(i)]=find(j)
        groups=collections.defaultdict(list)
        for i,wid in enumerate(wids):groups[find(i)].append(wid)
        for ww in groups.values():
            ww=sorted(ww);regions=combine_segments([p for wid in ww for p in segments[wid]])
            lid=stable_id("vgl_",[genome_key(windows[ww[0]]),regions])
            loci[lid]=dict(genome_unit_id=gid,window_ids=ww,regions=regions)
            for wid in ww:locus_for[wid]=lid
    return bygen,loci,locus_for
def positive(pred,retrieval):
    # Read the frozen classifier flag; do not rethreshold its raw sigmoid score.
    return yes(pred["saccharide_predicted"]) or float(retrieval["saccharide_fraction_top10"])>=0.5
def state(wids,pred,selected):
    if any(w in selected for w in wids):return "candidate"
    return "not_selected" if all(w in pred for w in wids) else "unresolved"
def combined(a,b):
    return "candidate" if "candidate" in [a,b] else "not_selected" if a==b=="not_selected" else "unresolved"
def support(a,b):
    return "Both" if a==b=="candidate" else "11-only" if b=="candidate" else "15-only" if a=="candidate" else ""
def guarded_run(cache,dataset,out,resources):
    cache=Path(cache).resolve();dataset=Path(dataset).resolve();out=Path(out).resolve()
    if out==cache or out==dataset or cache in out.parents or dataset in out.parents:raise ValueError("Output must be separate from inputs")
    rr=read(resources);files={}
    for r in rr:
        p=cache/r["key"]
        if sha(p)!=r["sha256"]:raise ValueError("Public resource checksum mismatch "+r["key"])
        files[r["key"]]=p
    w15=unique(read(files["windows15.tsv"]),"window_id");w11=unique(read(files["windows11.tsv"]),"window_id")
    p15=unique(read(files["predictions15.tsv"]),"window_id");p11=unique(read(files["predictions11.tsv"]),"window_id")
    r15=unique(read(files["retrieval15.tsv"]),"window_id");r11=unique(read(files["retrieval11.tsv"]),"window_id")
    if set(p15)!=set(r15) or set(p11)!=set(r11):raise ValueError("Retrieval/prediction membership")
    pairs=read(files["pairs.tsv"]);pair={r["window15_id"]:r["window11_id"] for r in pairs}
    if len(pair)!=757 or set(pair)!=set(w15) or len(set(pair.values()))!=757 or set(pair.values())!=set(w11):raise ValueError("Pair identity")
    for old,new in pair.items():
        for k in ["anchor_protein_id","genome_accession","anchor_start_1based","anchor_end_1based","anchor_strand"]:
            if w15[old][k]!=w11[new][k]:raise ValueError("Pair identity "+old+" "+k)
    ps={"15":p15,"11":p11};rs={"15":r15,"11":r11};ws={"15":w15,"11":w11}
    selected={c:{wid for wid,r in ps[c].items() if positive(r,rs[c][wid])} for c in ["15","11"]}
    bygen,loci,locus_for=reconstruct_groups(w15,read(files["genes15.tsv"]))
    public_genomes=unique(read(files["genomes15.tsv"]),"genome_unit_id");public_loci=unique(read(files["loci15.tsv"]),"locus_id")
    if set(bygen)!=set(public_genomes) or set(loci)!=set(public_loci):raise ValueError("Reconstructed grouping identities")
    for gid,wids in bygen.items():
        if set(wids)!=set(ids(public_genomes[gid]["window_ids"])):raise ValueError("Genome members "+gid)
    for lid,l in loci.items():
        frozen=public_loci[lid]
        if set(l["window_ids"])!=set(ids(frozen["window_ids"])) or l["regions"]!=json.loads(frozen["genomic_segments_1based"]):raise ValueError("Overlap membership "+lid)
    def statuses(old_ids):
        a=state(old_ids,p15,selected["15"]);b=state([pair[w] for w in old_ids],p11,selected["11"])
        return a,b
    audit_loci=[]
    for lid in public_loci:
        l=loci[lid];a,b=statuses(l["window_ids"]);f=public_loci[lid]
        audit_loci.append(dict(locus_id=lid,genome_unit_id=l["genome_unit_id"],status15=a,status11=b,combined_status=combined(a,b),same_locus_both_sizes=a==b=="candidate",window15_ids=f["window_ids"],anchor_protein_ids=f["anchor_protein_ids"],genomic_segments_1based=f["genomic_segments_1based"]))
    by_locus=collections.defaultdict(list)
    for r in audit_loci:by_locus[r["genome_unit_id"]].append(r)
    inputs=collections.defaultdict(set)
    for r in read(files["anchors.tsv"]):
        if r.get("window_id"):inputs[r["window_id"]].add(r["input_accession"])
    rows=[]
    for gid,g in public_genomes.items():
        ww=ids(g["window_ids"]);nn=[pair[x] for x in ww];ll=by_locus[gid];a,b=statuses(ww);cc=combined(a,b)
        cls=any(yes(ps[c][wid]["saccharide_predicted"]) for c,ww2 in [("15",ww),("11",nn)] for wid in ww2 if wid in ps[c])
        both=a==b=="candidate";same=any(r["same_locus_both_sizes"] for r in ll)
        rows.append(dict(genome_unit_id=gid,source_genome_accessions=g["source_genome_accessions"],source_organism_labels=g["source_organism_labels"],source_ncbi_taxid=g["taxid"],genome_sequence_sha256=g["genome_sequence_sha256"],n_cohort_anchor_accessions=len(set().union(*(inputs[w] for w in ww))),n_anchor_occurrences_alias_collapsed=len({(w15[w]["anchor_start_1based"],w15[w]["anchor_end_1based"],w15[w]["anchor_strand"]) for w in ww}),n_occurrence_windows=len(ww),n_fixed_baseline_loci=len(ll),n_candidate_loci15=sum(r["status15"]=="candidate" for r in ll),n_candidate_loci11=sum(r["status11"]=="candidate" for r in ll),n_combined_candidate_loci=sum(r["combined_status"]=="candidate" for r in ll),status15=a,status11=b,combined_status=cc,window_support=support(a,b),n_same_locus_both=sum(r["same_locus_both_sizes"] for r in ll),genome_both_sizes=both,same_locus_both_sizes=same,both_sizes_different_loci_only=both and not same,classifier_positive_either_size=cls,retrieval_only_candidate=cc=="candidate" and not cls,retrieval_only_window_present=any(not yes(ps[c][wid]["saccharide_predicted"]) and wid in selected[c] for c,ww2 in [("15",ww),("11",nn)] for wid in ww2 if wid in ps[c]),n_inferred_windows15=sum(x in p15 for x in ww),n_inferred_windows11=sum(x in p11 for x in nn),fully_inferred15=all(x in p15 for x in ww),fully_inferred11=all(x in p11 for x in nn),any_inferred=any(x in p15 or pair[x] in p11 for x in ww),n_complete_windows15=sum(w15[x]["n_genes"]=="15" for x in ww),n_complete_windows11=sum(w11[x]["n_genes"]=="11" for x in nn),unresolved_other_size=cc=="candidate" and "unresolved" in [a,b],window15_ids=g["window_ids"],fixed_baseline_locus_ids=";".join(r["locus_id"] for r in ll)))
    # The expected original master is read only AFTER all scientific selection fields are reconstructed.
    frozen=unique(read(dataset/"outputs/master_genome_dataset.tsv"),"genome_unit_id")
    differences=[]
    for r in rows:
        if r["genome_unit_id"] not in frozen:raise ValueError("Unexpected genome")
        for k,v in r.items():
            if str(v)!=frozen[r["genome_unit_id"]][k]:differences.append((r["genome_unit_id"],k,str(v),frozen[r["genome_unit_id"]][k]))
    if set(frozen)!=set(bygen) or differences:raise ValueError("Master row discrepancies: "+repr(differences[:10]))
    pub_laudit=unique(read(files["locus_audit.tsv"]),"locus_id")
    for r in audit_loci:
        if any(str(v)!=pub_laudit[r["locus_id"]][k] for k,v in r.items()):raise ValueError("Locus audit "+r["locus_id"])
    registry={(r["level"],r["entity_id"]):r for r in read(files["paired_registry.tsv"]) if r["set"]=="C"}
    for level,rr,idcol in [("genome",rows,"genome_unit_id"),("locus",audit_loci,"locus_id")]:
        for r in rr:
            expected=registry[level,r[idcol]]
            if any(r[k]!=expected[k] for k in ["status15","status11"]):raise ValueError("Original registry "+r[idcol])
    candidates=[r for r in rows if r["combined_status"]=="candidate"];cand={r["genome_unit_id"]:r for r in candidates}
    pub245=unique(read(files["membership245.tsv"]),"genome_unit_id")
    if set(cand)!=set(pub245):raise ValueError("Original 245 membership")
    for g,r in cand.items():
        if r["window_support"]!=pub245[g]["historical_window_support"]:raise ValueError("Original 245 support")
    exclusion=read(dataset/"outputs/FigureB_exclusion_audit.tsv")
    excluded={g for g,r in cand.items() if set(ids(r["source_genome_accessions"]))&EXCLUDE}
    if excluded!={r["genome_unit_id"] for r in exclusion} or excluded!={r["genome_unit_id"] for r in read(files["exclusions3.tsv"])}:raise ValueError("Three exclusions")
    for r in exclusion:
        acc=r["source_genome_accessions"];gb=files[acc+".gb"];raw=gb.read_text();length=int(re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp",raw,re.M).group(1))
        if sha(gb)!=r["source_sha256"] or len(re.findall(r"^     CDS\s",raw,re.M))!=1 or length!=int(r["source_length_nt"]):raise ValueError("Single-CDS source proof "+acc)
        old=r["window15_id"];new=pair[old]
        if w15[old]["n_genes"]!="1" or w11[new]["n_genes"]!="1" or cand[r["genome_unit_id"]]["window_support"]!="Both":raise ValueError("Exclusion input "+acc)
    sankey=[r for r in candidates if r["genome_unit_id"] not in excluded]
    for g,r in pub245.items():
        if yes(r["included_in_frozen_242_sankey"])!=(g not in excluded):raise ValueError("242 membership "+g)
    counts=dict(positive_windows15=len(selected["15"]),positive_windows11=len(selected["11"]),positive_occurrences=sum(x in selected["15"] or pair[x] in selected["11"] for x in w15),positive_loci=sum(r["combined_status"]=="candidate" for r in audit_loci),genome_units=len(rows),candidates=len(candidates),negatives=sum(r["combined_status"]=="not_selected" for r in rows),unresolved=sum(r["combined_status"]=="unresolved" for r in rows),sankey=len(sankey))
    counts.update(collections.Counter(r["window_support"] for r in candidates))
    if counts!=EXPECTED:raise ValueError("Unexpected counts "+repr(counts))
    sank_support=dict(collections.Counter(r["window_support"] for r in sankey))
    if sank_support!={"Both":187,"11-only":29,"15-only":26}:raise ValueError("Sankey support")
    cat_old=next(w for w,r in w15.items() if r["genome_accession"]=="KY684083.1" and r["anchor_protein_id"]=="ARF08732.1");cat_new=pair[cat_old]
    cat=dict(accession="KY684083.1",anchor="ARF08732.1",window15=cat_old,window11=cat_new,locus_id=locus_for[cat_old],score15=float(p15[cat_old]["saccharide_probability"]),score11=float(p11[cat_new]["saccharide_probability"]),p10_15=float(r15[cat_old]["saccharide_fraction_top10"]),p10_11=float(r11[cat_new]["saccharide_fraction_top10"]),window15_candidate=cat_old in selected["15"],window11_candidate=cat_new in selected["11"],genome_support=cand[genome_id(w15[cat_old])]["window_support"])
    if abs(cat["score15"]-.024441218)>1e-9 or abs(cat["score11"]-.878088236)>1e-9 or (cat["p10_15"],cat["p10_11"],cat["window15_candidate"],cat["window11_candidate"],cat["genome_support"])!=(.2,.8,False,True,"Both"):raise ValueError("Catovirus consistency")
    # Audit-only support map; no frozen probabilities or embeddings are rewritten.
    support_rows=[]
    for old,w in sorted(w15.items()):
        for c,wid in [("15",old),("11",pair[old])]:
            if wid in selected[c]:
                a=yes(ps[c][wid]["saccharide_predicted"]);b=float(rs[c][wid]["saccharide_fraction_top10"])>=.5
                support_rows.append(dict(genome_unit_id=genome_id(w),locus_id=locus_for[old],paired_occurrence15_id=old,window_size=c,window_id=wid,anchor_protein_id=w["anchor_protein_id"],genome_accession=w["genome_accession"],evidence_type="classifier_and_retrieval" if a and b else "classifier_only" if a else "retrieval_only"))
    write(out/"selection_genomes.tsv",rows);write(out/"selection_loci.tsv",audit_loci);write(out/"positive_window_support.tsv",support_rows)
    write(out/"candidate_membership245.tsv",[dict(genome_unit_id=r["genome_unit_id"],source_genome_accessions=r["source_genome_accessions"],window_support=r["window_support"]) for r in candidates])
    write(out/"sankey_membership242.tsv",[dict(genome_unit_id=r["genome_unit_id"],source_genome_accessions=r["source_genome_accessions"],window_support=r["window_support"]) for r in sankey])
    report=dict(status="PASS",version=VERSION,counts=counts,sankey_support=sank_support,genome_fields_compared=len(rows[0]),genome_field_comparisons=len(rows)*len(rows[0]),genome_rows_compared=len(rows),locus_rows_compared=len(loci),all_genome_and_locus_identities_reconstructed=True,positive_support_rows=len(support_rows),membership245_exact=True,membership242_exact=True,discrepancies=[],catovirus=cat,public_resources_sha256_verified=len(files),python=sys.version,third_party_python_packages_required=False,frozen_metadata_note="All original taxonomy/display metadata are preserved in the dataset master; selection-only reproduction independently compares every emitted selection field, without taxonomy reconciliation.",input_hashes={r["key"]:r["sha256"] for r in read(resources)})
    (out/"reproduction_validation.json").write_text(json.dumps(report,indent=2)+"\n")
    return report
def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument("--cache",required=True);a.add_argument("--dataset",required=True);a.add_argument("--out",required=True);a.add_argument("--resources",default=str(Path(__file__).resolve().parents[1]/"references/public_resources.tsv"));o=a.parse_args()
    report=guarded_run(o.cache,o.dataset,o.out,o.resources);print(json.dumps({k:v for k,v in report.items() if k not in ["input_hashes"]},indent=2))
if __name__=="__main__":main()

