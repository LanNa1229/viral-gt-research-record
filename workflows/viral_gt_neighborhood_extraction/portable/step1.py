#!/usr/bin/env python3
"""Portable Step 1 adapter (new, 2026-10-26), with hash-pinned historical core functions.
No model imports, network access, taxonomy inference or changes to historical files.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,ast,collections,csv,hashlib,importlib.util,json,platform,time
HOME=Path(__file__).resolve().parents[1]
PIN={'extract':'9880dab746f4762ba4f7d8aacacb452feec3c0fc2384fcd71f84d42f8d3768e9','experiment':'8d6fc8848b720de15b859699bbe4685506d56bdabcce00a3fdf73201a6e182a8','analysis':'2a85da9400f6b00aa505c37e1b66b70d1ee8c0ef1e4c221eb1c4bdf92a3bddf7'}
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields=None):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    with Path(p).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fields or sorted(set().union(*(r.keys() for r in rows))),delimiter='\t');w.writeheader();w.writerows(rows)
def save(p,o):Path(p).write_text(json.dumps(o,indent=2)+'\n')
def resolved(root,value):
    path=(root/value).resolve()
    if not path.is_relative_to(root):raise ValueError('Manifest path leaves data root: '+value)
    return path

def historical():
    ep=HOME/'historical/baseline_594_v2/src/viral_gt_mining/extract.py';xp=HOME/'historical/window11/experiment.py';ap=HOME/'historical/baseline_594_v2/src/viral_gt_mining/analysis.py'
    for key,p in [('extract',ep),('experiment',xp),('analysis',ap)]:
        if sha(p)!=PIN[key]:raise ValueError('Historical source hash mismatch: '+str(p))
    spec=importlib.util.spec_from_file_location('original_step1_extraction',ep);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    # Select unchanged function ASTs only; do not execute original top-level imports/launchers.
    ns={'collections':collections,'json':json,'hashlib':hashlib,'as_bool':lambda x:str(x).lower() in {'true','1','yes'}}
    for path,names in [(xp,{'shortened'}),(ap,{'unit_status','genome_key','genome_unit_id','group_genomes','combine_segments','window_segments','merge_loci'})]:
        tree=ast.parse(path.read_text(encoding='utf-8-sig'))
        selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
        selected += [n for n in tree.body if path==ap and isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ['GENOME_DEFINITION','LOCUS_DEFINITION'] for t in n.targets)]
        if len([n for n in selected if isinstance(n,ast.FunctionDef)])!=len(names):raise ValueError('Historical function set changed')
        exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),ns)
    return module,ns

def load_inputs(root,anchors=None):
    root=root.resolve();source=read(root/'portable_inputs/genbanks.tsv');aliases=read(root/'portable_inputs/aliases.tsv')
    paths=[]
    for r in source:
        p=resolved(root,r['path'])
        if sha(p)!=r['sha256']:raise ValueError('GenBank checksum mismatch: '+r['path'])
        paths.append(str(p))
    aa=collections.defaultdict(list)
    for r in aliases:
        if anchors is not None and r['input_accession'] not in anchors:continue
        r=dict(r);p=resolved(root,r['evidence_path'])
        if sha(p)!=r['evidence_sha256']:raise ValueError('Alias evidence checksum mismatch: '+r['evidence_path'])
        r['evidence_path']=str(p);aa[r['input_accession']].append(r)
    return paths,dict(aa)

def derive11(oldrun,out,root,ex,ns):
    old=read(oldrun/'viral_gt_windows.tsv');groups=collections.defaultdict(list)
    for g in read(oldrun/'viral_gt_window_genes.tsv'):groups[g['window_id']].append(g)
    relations={r['window15_id']:r for r in read(root/'relationships/anchor_window_locus_genome.tsv') if r['window15_id']}
    rows=[];genes=[];pairs=[]
    for w in old:
        wid=w['window_id'];gg=sorted(groups[wid],key=lambda g:int(g['relative_position']));idx=int(w['anchor_gene_index']);newid='vgt11_'+ex.digest(json.dumps([w['genome_accession'],w['genome_sequence_sha256'],w['anchor_gene_id'],5,'coordinate_no_wrap_v1'],separators=(',',':')))[:24]
        row,newgenes=ns['shortened'](w,gg,newid,idx)
        # Recorded identity crosswalk supplies fixed baseline loci, not candidate labels/scores.
        rel=relations.get(wid,{})
        row.update(baseline_locus_id=rel.get('native15_locus_id',''),genome_unit_id=ns['genome_unit_id'](ns['genome_key'](w)))
        rows.append(row);genes.extend(newgenes)
        pairs.append(dict(window15_id=wid,window11_id=newid,genome_unit_id=row['genome_unit_id'],baseline_locus_id=row['baseline_locus_id'],n15=len(gg),n11=len(newgenes)))
    out.mkdir();write(out/'viral_gt_windows.tsv',rows);write(out/'viral_gt_window_genes.tsv',genes);write(out/'window15_window11_mapping.tsv',pairs)
    with (out/'viral_gt_window_proteins.faa').open('w') as f:
        for g in genes:
            if g['sequence']:f.write('>'+g['window_id']+'|'+format(int(g['relative_position']),'+d')+'|'+g['gene_id']+' protein_id='+g['protein_id']+'\n'+g['sequence']+'\n')
    return rows,genes

def grouping(windows,genes,out,ns):
    # Empty prediction mapping means no inferred classification is manufactured.
    loci,_=ns['merge_loci'](windows,genes,{});_,genomes=ns['group_genomes'](windows,{})
    drop={'classification_status','n_analyzed_windows','n_positive_windows','fully_analyzed'}
    write(out/'structural_loci.tsv',[{k:v for k,v in r.items() if k not in drop} for r in loci]);write(out/'structural_genome_units.tsv',[{k:v for k,v in r.items() if k not in drop} for r in genomes]);return len(loci),len(genomes)

def extract(args):
    root=args.data_root.resolve();out=args.output.resolve()
    if out.exists():raise ValueError('Refusing to overwrite an existing output directory')
    if out.is_relative_to(root) or out.is_relative_to(HOME):raise ValueError('Output must be outside the data and code packages')
    ex,ns=historical();anchorfile=resolved(root,args.anchors);anchors=ex.load_anchors(anchorfile)
    if args.expected_anchors is not None and len(set(anchors))!=args.expected_anchors:raise ValueError('Unexpected anchor cohort size')
    paths,aliases=load_inputs(root,set(anchors));out.mkdir(parents=True)
    result=ex.extract(anchorfile,None,paths,out/'window15',out/'regenerable_cache',aliases)
    w15=read(out/'window15/viral_gt_windows.tsv');g15=read(out/'window15/viral_gt_window_genes.tsv');w11,g11=derive11(out/'window15',out/'window11',root,ex,ns)
    sizes={}
    for name,w,g,nom in [('window15',w15,g15,15),('window11',w11,g11,11)]:
        loci,genomes=grouping(w,g,out/name,ns);sizes[name]=dict(windows=len(w),gene_rows=len(g),complete=sum(int(r['n_genes'])==nom for r in w),shortened=sum(int(r['n_genes'])<nom for r in w),eligible_for_later_inference=sum(r['inference_eligible']=='true' for r in w),native_loci=loci,genome_units=genomes)
    save(out/'portable_execution.json',dict(adapter='new portability adaptation 2026-10-26; unchanged historical function ASTs isolated from broad launchers',historical_code_pins=PIN,python=sys.version,biopython=__import__('Bio').__version__,input_anchor_sha256=sha(anchorfile),source_files=len(paths),extraction15=result,counts=sizes,inference_calls=0,network_calls=0,classification_not_recomputed=True,original_source_files_modified=False))
    print(json.dumps(sizes,indent=2))

def verify(args):
    root=args.data_root.resolve();source=read(root/'outputs/baseline_594_v2/viral_gt_windows.tsv');genes=read(root/'outputs/baseline_594_v2/viral_gt_window_genes.tsv');s11=read(root/'outputs/window11_594_v1/viral_gt_windows.tsv');g11=read(root/'outputs/window11_594_v1/viral_gt_window_genes.tsv');am=read(root/'outputs/baseline_594_v2/anchor_mapping.tsv');relations=read(root/'relationships/anchor_window_locus_genome.tsv')
    anchors=[x for x in (root/'inputs/baseline_594_anchors.txt').read_text().splitlines() if x];member=list(csv.reader((root/'inputs/upstream/cluster_membership.tsv').open(),delimiter='\t'));members=[r[1] for r in member if r[0]=='A0A3Q8Q3U2']
    assert len(anchors)==len(set(anchors))==594 and set(anchors)==set(members)
    assert len({r['input_accession'] for r in am if r['status']=='mapped'})==593 and {r['input_accession'] for r in am if r['status']=='unresolved'}=={'UPI001FBC14AE'}
    assert len(source)==len(s11)==757 and len(genes)==10688 and len(g11)==7967
    counts15=sum(r['status']=='PASS' for r in read(root/'outputs/baseline_594_v2/inference_qc.tsv'));counts11=sum(r['status']=='PASS' for r in read(root/'outputs/window11_594_v1/inference_qc.tsv'));assert (counts15,counts11)==(737,740)
    assert sum(int(r['n_genes'])==15 for r in source)==616 and sum(int(r['n_genes'])==11 for r in s11)==658
    assert len(relations)==758 and len({r['genome_unit_id'] for r in relations if r['genome_unit_id']})==466
    assert len({r['native15_locus_id'] for r in relations if r['native15_locus_id']})==565 and len({r['native11_locus_id'] for r in relations if r['native11_locus_id']})==575
    ex,ns=historical();paths,_=load_inputs(root);r=collections.defaultdict(list);gg=collections.defaultdict(list)
    for x in genes:r[x['window_id']].append(x)
    for x in g11:gg[x['window_id']].append(x)
    map11={x['baseline_window_id']:x for x in s11}
    for old in source:
        current=map11[old['window_id']];assert [x['gene_id'] for x in sorted(r[old['window_id']],key=lambda g:int(g['relative_position'])) if abs(int(x['relative_position']))<=5]==[x['gene_id'] for x in sorted(gg[current['window_id']],key=lambda g:int(g['relative_position']))]
    result=dict(passed=True,anchors=594,mapped=593,unresolved='UPI001FBC14AE',occurrences_each_size=757,genome_units=466,native_loci15=565,native_loci11=575,checked_sources=len(paths),recorded_inference_count15=counts15,recorded_inference_count11=counts11,all_frozen_757_subset_pairs_verified=True,extraction_performed=False,inference_performed=False)
    print(json.dumps(result,indent=2))

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('extract');q.add_argument('--data-root',type=Path,required=True);q.add_argument('--output',type=Path,required=True);q.add_argument('--anchors',default='inputs/baseline_594_anchors.txt');q.add_argument('--expected-anchors',type=int);q.set_defaults(func=extract)
    q=sub.add_parser('verify');q.add_argument('--data-root',type=Path,required=True);q.set_defaults(func=verify)
    a=p.parse_args();a.func(a)
if __name__=='__main__':main()
