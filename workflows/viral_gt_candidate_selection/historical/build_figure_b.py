"""Reconstruct Figure B from frozen outputs; no model calls or source writes."""
import sys
sys.dont_write_bytecode = True
import csv, json, hashlib, collections, pathlib, shutil, urllib.request, datetime, xml.etree.ElementTree as ET
ROOT = pathlib.Path(__file__).resolve().parents[3]
OUT = pathlib.Path(__file__).resolve().parents[1]
BASE = ROOT/'analysis/viral_gt_mining_20261017/runs/baseline_594_v2'
R11 = BASE.parent/'window11_594_v1'
SIZE = ROOT/'analysis/viral_gt_window_size_20261017'
def read(p):
    with p.open(encoding='utf-8-sig', newline='') as f: return list(csv.DictReader(f, delimiter='\t'))
def write(p, rows, fields=None):
    rows=list(rows); fields=fields or list(rows[0])
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f, fields, delimiter='\t'); w.writeheader(); w.writerows(rows)
def sha(p): return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
checks=[]
def check(ok, message):
    checks.append(message)
    if not ok: raise AssertionError(message)
def verify_sources():
    inventories=[(ROOT/'analysis/viral_gt_interpretation_20261017/freeze/baseline_files.tsv', BASE)]
    for name in ['viral_gt_interpretation_20261017','viral_gt_retrieval_rescue_20261017','viral_gt_window_size_20261017']:
        p=ROOT/'analysis'/name; inventories.append((p/'artifact_checksums.tsv',p))
    inventories.append((R11/'artifact_checksums.tsv',R11))
    for inv, base in inventories:
        for r in read(inv):
            p=pathlib.Path(r['path'].replace('\\','/'))
            if not p.is_absolute(): p=base/p
            check(p.exists() and sha(p)==r['sha256'], 'frozen hash: '+str(p.relative_to(ROOT)))

def build():
    verify_sources()
    genomes=read(BASE/'genome_units.tsv'); loci=read(BASE/'candidate_loci.tsv')
    windows={r['window_id']:r for r in read(BASE/'viral_gt_windows.tsv')}
    pairs=read(R11/'window15_window11_mapping.tsv')
    preds=[{r['window_id']:r for r in read(b/'viral_bgc_predictions.tsv')} for b in [BASE,R11]]
    retrieval=[{r['window_id']:r for r in read(b/'saccharide_localization.tsv')} for b in [BASE,R11]]
    pmap={r['window15_id']:r['window11_id'] for r in pairs}
    meta={r['genome_unit_id']:r for r in read(ROOT/'analysis/viral_gt_interpretation_20261017/tables/genome_relatedness_metadata.tsv')}
    tax={r['taxid']:r for r in read(BASE/'taxonomy_cache.tsv')}
    lineages={}; nodes={}; parents={}; provenance=[]
    for source in sorted({r['source_path'] for r in tax.values()}):
        p=pathlib.Path(source); dest=OUT/'sources/taxonomy'/p.name
        shutil.copyfile(p,dest); shutil.copyfile(pathlib.Path(source+'.json'),pathlib.Path(str(dest)+'.json'))
        prov=json.loads(pathlib.Path(source+'.json').read_text())
        provenance.append(dict(file=p.name,sha256=sha(p),retrieved_utc=prov.get('retrieved_utc',''),url=prov.get('url','')))
        for t in ET.parse(p).getroot().findall('Taxon'):
            path=[]
            for n in list(t.findall('LineageEx/Taxon'))+[t]:
                tid=n.findtext('TaxId'); node=dict(taxid=tid,name=n.findtext('ScientificName'),rank=n.findtext('Rank'))
                check(tid not in nodes or nodes[tid]==node,'consistent taxon '+tid)
                nodes[tid]=node
                if path:
                    check(tid not in parents or parents[tid]==path[-1], 'consistent parent '+tid)
                    parents[tid]=path[-1]
                path.append(tid)
            lineages[t.findtext('TaxId')]=path
    write(OUT/'sources/taxonomy/provenance.tsv',provenance)
    for t in tax.values():
        path=lineages[t['resolved_taxid']]
        for rank in ['family','genus','species']:
            ids=[i for i in path if nodes[i]['rank']==rank]
            check((ids[-1] if ids else '')==t[rank+'_taxid'],'cached rank matches XML '+t['taxid']+' '+rank)
    write(OUT/'sources/taxonomy/nodes.tsv',[dict(**n,parent_taxid=parents.get(i,'')) for i,n in sorted(nodes.items())])
    inputs=collections.defaultdict(set)
    for r in read(BASE/'anchor_mapping.tsv'):
        if r['window_id']: inputs[r['window_id']].add(r['input_accession'])
    def state(wids, k, evidence='union'):
        ids=[w if k==0 else pmap[w] for w in wids]
        def positive(w):
            a=preds[k][w]['saccharide_predicted'].lower()=='true'
            b=float(retrieval[k][w]['saccharide_fraction_top10'])>=.5
            return a if evidence=='classifier' else (b and not a) if evidence=='retrieval' else a or b
        if any(positive(w) for w in ids if w in preds[k]): return 'candidate'
        return 'not_selected' if all(w in preds[k] for w in ids) else 'unresolved'
    registry={(r['level'],r['entity_id'],r['set']):r for r in read(SIZE/'tables/candidate_recovery_comparison.tsv')}
    lout=[]; bygen=collections.defaultdict(list)
    for l in loci:
        ids=l['window_ids'].split(';'); a,b=state(ids,0),state(ids,1)
        r=registry['locus',l['locus_id'],'C']
        check((a,b)==(r['status15'],r['status11']),'locus registry '+l['locus_id'])
        row=dict(locus_id=l['locus_id'],genome_unit_id=l['genome_unit_id'],status15=a,status11=b,
                 combined_status='candidate' if 'candidate' in [a,b] else 'not_selected' if a==b=='not_selected' else 'unresolved',
                 same_locus_both_sizes=a==b=='candidate',window15_ids=l['window_ids'],anchor_protein_ids=l['anchor_protein_ids'],genomic_segments_1based=l['genomic_segments_1based'])
        lout.append(row); bygen[l['genome_unit_id']].append(row)
    write(OUT/'tables/locus_audit.tsv',lout)
    master=[]
    for g in genomes:
        gid=g['genome_unit_id']; ids=g['window_ids'].split(';'); ls=bygen[gid]; m=meta[gid]; t=tax.get(g['taxid'],{})
        a,b=state(ids,0),state(ids,1); reg=registry['genome',gid,'C']
        check((a,b)==(reg['status15'],reg['status11']),'genome registry '+gid)
        generic=m['generic_source_label'].lower()=='true'
        family=t.get('family_taxid',''); genus=t.get('genus_taxid','')
        combined='candidate' if 'candidate' in [a,b] else 'not_selected' if a==b=='not_selected' else 'unresolved'
        both=a==b=='candidate'; support='Both' if both else '11-only' if b=='candidate' else '15-only' if a=='candidate' else ''
        cls=any(state(ids,k,'classifier')=='candidate' for k in [0,1])
        r=dict(genome_unit_id=gid,source_genome_accessions=g['source_genome_accessions'],source_organism_labels=g['source_organism_labels'],source_ncbi_taxid=g['taxid'],genome_sequence_sha256=g['genome_sequence_sha256'],generic_source_label=generic,
               taxonomy_source='NCBI Taxonomy EFetch cached XML',taxonomy_resolved_taxid=t.get('resolved_taxid',''),taxonomy_xml_sha256=t.get('source_sha256',''),taxonomy_xml_file=pathlib.Path(t.get('source_path','')).name,
               family=t.get('family',''),family_taxid=family,genus=t.get('genus',''),genus_taxid=genus,species=t.get('species',''),species_taxid=t.get('species_taxid',''),
               conservative_family_included=bool(family and not generic),conservative_exclusion_reason='generic_source_label' if generic else 'family_rank_unresolved' if not family else '',
               taxonomy_uncertainty='generic_source_inherited_ranks' if generic else 'family_rank_unresolved' if not family else 'non_generic_rank_resolved',
               n_cohort_anchor_accessions=len(set().union(*(inputs[w] for w in ids))),n_anchor_occurrences_alias_collapsed=len({(windows[w]['anchor_start_1based'],windows[w]['anchor_end_1based'],windows[w]['anchor_strand']) for w in ids}),
               n_occurrence_windows=len(ids),n_fixed_baseline_loci=len(ls),n_candidate_loci15=sum(l['status15']=='candidate' for l in ls),n_candidate_loci11=sum(l['status11']=='candidate' for l in ls),n_combined_candidate_loci=sum(l['combined_status']=='candidate' for l in ls),
               status15=a,status11=b,combined_status=combined,window_support=support,n_same_locus_both=sum(l['same_locus_both_sizes'] for l in ls),genome_both_sizes=both,same_locus_both_sizes=any(l['same_locus_both_sizes'] for l in ls),both_sizes_different_loci_only=both and not any(l['same_locus_both_sizes'] for l in ls),
               classifier_positive_either_size=cls,retrieval_only_candidate=combined=='candidate' and not cls,retrieval_only_window_present=any(state(ids,k,'retrieval')=='candidate' for k in [0,1]),
               n_inferred_windows15=sum(w in preds[0] for w in ids),n_inferred_windows11=sum(pmap[w] in preds[1] for w in ids),fully_inferred15=all(w in preds[0] for w in ids),fully_inferred11=all(pmap[w] in preds[1] for w in ids),any_inferred=any(w in preds[0] or pmap[w] in preds[1] for w in ids),
               n_complete_windows15=sum(windows[w]['n_genes']=='15' for w in ids),n_complete_windows11=sum(int(q['n11'])==11 for q in pairs if q['window15_id'] in ids),unresolved_other_size=combined=='candidate' and 'unresolved' in [a,b],window15_ids=g['window_ids'],fixed_baseline_locus_ids=';'.join(l['locus_id'] for l in ls))
        master.append(r)
    write(OUT/'tables/master_genome_dataset.tsv',master)
    check(len(master)==466 and len({r['genome_unit_id'] for r in master})==466,'466 unique units')
    check(sum(r['combined_status']=='candidate' for r in master)==245,'245 positive genomes')
    check(sum(l['combined_status']=='candidate' for l in lout)==288,'288 positive fixed loci')
    check(collections.Counter(r['window_support'] for r in master if r['combined_status']=='candidate')=={'Both':190,'11-only':29,'15-only':26},'190/29/26 support')
    def summarize(rows,rank):
        groups=collections.defaultdict(list)
        for r in rows: groups[r[rank+'_taxid']].append(r)
        out=[]
        for tid,rs in sorted(groups.items(),key=lambda p:nodes[p[0]]['name']):
            positive=sum(r['combined_status']=='candidate' for r in rs); negative=sum(r['combined_status']=='not_selected' for r in rs); unresolved=len(rs)-positive-negative
            supports=collections.Counter(r['window_support'] for r in rs)
            conf='generic_source_only' if all(r['generic_source_label'] for r in rs) else 'mixed_generic_and_non_generic' if any(r['generic_source_label'] for r in rs) else 'non_generic_rank_resolved'
            row=dict(tip_id=rank+'_'+tid,taxid=tid,scientific_name=nodes[tid]['name'],taxonomy_confidence=conf,total_selected=len(rs),any_inferred=sum(r['any_inferred'] for r in rs),outcome_resolved=positive+negative,candidate_positive=positive,fully_analyzed_negative=negative,unresolved=unresolved,positive_fraction_pct=100*positive/(positive+negative) if positive+negative else '',both=supports['Both'],only11=supports['11-only'],only15=supports['15-only'],same_locus_both_genomes=sum(r['same_locus_both_sizes'] for r in rs),same_locus_both_loci=sum(r['n_same_locus_both'] for r in rs),classifier_positive_genomes=sum(r['classifier_positive_either_size'] for r in rs),retrieval_only_genomes=sum(r['retrieval_only_candidate'] for r in rs),combined_candidate_loci=sum(r['n_combined_candidate_loci'] for r in rs),generic_source_units=sum(r['generic_source_label'] for r in rs),genome_unit_ids=';'.join(r['genome_unit_id'] for r in rs))
            check(row['both']+row['only11']+row['only15']==positive,'support partition '+rank+tid)
            check(positive+negative+unresolved==len(rs),'resolved partition '+rank+tid)
            out.append(row)
        return out
    allstats={}; trees={}
    for variant,rank,selected,folder in [
        ('conservative','family',[r for r in master if r['conservative_family_included']],OUT),
        ('extended','family',[r for r in master if r['family_taxid']],OUT/'extended'),
        ('reliable_genus','genus',[r for r in master if not r['generic_source_label'] and r['genus_taxid']],OUT/'supplementary_genus')]:
        folder.mkdir(exist_ok=True); stats=summarize(selected,rank); allstats[variant]=stats
        write(folder/(rank+'_taxonomy_figure_data.tsv'),stats)
        write(OUT/'tables'/(variant+'_'+rank+'_summary.tsv'),stats)
        tree={'id':'10239','children':{}}
        for row in stats:
            path=[]; tid=row['taxid']
            while tid:
                path.append(tid); tid=parents.get(tid)
            path.reverse(); path=path[path.index('10239'):]
            branch=tree
            for tid in path[1:]: branch=branch['children'].setdefault(tid,{'id':tid,'children':{}})
            branch['tip_id']=row['tip_id']
        def newick(n):
            if not n['children']: return n['tip_id']
            children=sorted(n['children'].values(),key=lambda c:nodes[c['id']]['name'])
            return '('+','.join(newick(c) for c in children)+")'"+nodes[n['id']]['name'].replace("'","''")+' (NCBI '+n['id']+")'"
        nw=newick(tree)+';\n'; (folder/('viral_'+rank+'_taxonomy.nwk')).write_text(nw)
        check(':' not in nw,'no branch lengths '+variant)
        colors=['#0072B2','#E69F00','#CC79A7']; confidence={'non_generic_rank_resolved':'#009E73','mixed_generic_and_non_generic':'#E69F00','generic_source_only':'#D55E00'}
        datasets={
          'itol_evidence_multibar.txt':(['DATASET_MULTIBAR','SEPARATOR TAB','DATASET_LABEL\tEvidence among union candidate genomes','COLOR\t#0072B2','FIELD_COLORS\t#0072B2\t#D55E00','FIELD_LABELS\tClassifier-positive\tRetrieval-only','ALIGN_FIELDS\t0','WIDTH\t160','MARGIN\t15'],[[r['tip_id'],r['classifier_positive_genomes'],r['retrieval_only_genomes']] for r in stats]),
          'itol_window_support_multibar.txt':(['DATASET_MULTIBAR','SEPARATOR TAB','DATASET_LABEL\tCandidate genome support','COLOR\t#0072B2','FIELD_COLORS\t'+'\t'.join(colors),'FIELD_LABELS\tBoth\t11-only\t15-only','ALIGN_FIELDS\t0','WIDTH\t250','MARGIN\t15','LEGEND_TITLE\tWindow-size support','LEGEND_SHAPES\t1\t1\t1','LEGEND_COLORS\t'+'\t'.join(colors),'LEGEND_LABELS\tBoth\t11-only\t15-only'], [[r['tip_id'],r['both'],r['only11'],r['only15']] for r in stats]),
          'itol_positive_fraction.txt':(['DATASET_SIMPLEBAR','SEPARATOR TAB','DATASET_LABEL\tCandidate % among outcome-resolved units','COLOR\t#595959','WIDTH\t140','MARGIN\t20','DATASET_SCALE\t0\t25\t50\t75\t100','SHOW_VALUE\t0'],[[r['tip_id'],round(r['positive_fraction_pct'],6)] for r in stats if r['outcome_resolved']]),
          'itol_positive_counts.txt':(['DATASET_TEXT','SEPARATOR TAB','DATASET_LABEL\tPositive / outcome-resolved','COLOR\t#303030','MARGIN\t15'],[[r['tip_id'],str(r['candidate_positive'])+'/'+str(r['outcome_resolved'])+(' (NA)' if not r['outcome_resolved'] else ''),-1,'#303030','normal',1,0] for r in stats]),
          'itol_taxonomy_confidence.txt':(['DATASET_COLORSTRIP','SEPARATOR TAB','DATASET_LABEL\tSource taxonomy status','COLOR\t#009E73','COLOR_BRANCHES\t0','STRIP_WIDTH\t12','SHOW_STRIP_LABELS\t0','LEGEND_TITLE\tTaxonomy status','LEGEND_SHAPES\t1\t1\t1','LEGEND_COLORS\t#009E73\t#E69F00\t#D55E00','LEGEND_LABELS\tNon-generic rank resolved\tMixed generic/non-generic\tGeneric source only'],[[r['tip_id'],confidence[r['taxonomy_confidence']],r['taxonomy_confidence']] for r in stats]),
          'itol_'+rank+'_labels.txt':(['LABELS','SEPARATOR TAB'],[[r['tip_id'],r['scientific_name']] for r in stats])}
        for filename,(header,data) in datasets.items():
            ids=[str(r[0]) for r in data]; tips={r['tip_id'] for r in stats}
            check(len(ids)==len(set(ids)) and set(ids)<=tips,'annotation IDs '+variant+filename)
            if filename!='itol_positive_fraction.txt': check(set(ids)==tips,'complete annotations '+variant+filename)
            (folder/filename).write_text('\n'.join(header+['DATA']+['\t'.join(map(str,row)) for row in data])+'\n',encoding='utf-8')
        check(sum(r['total_selected'] for r in stats)==len(selected),'single taxon assignment '+variant)
        trees[variant]=tree
        excluded=[dict(**r,tree_exclusion_reason='family_rank_unresolved' if variant=='extended' else ('generic_source_label' if r['generic_source_label'] else 'genus_rank_unresolved') if rank=='genus' else r['conservative_exclusion_reason'] or 'family_rank_unresolved') for r in master if r not in selected]
        write(OUT/'tables'/(variant+'_excluded_genomes.tsv'),excluded)
    buckets=collections.defaultdict(list)
    for r in master:
        if not r['conservative_family_included']: buckets[r['conservative_exclusion_reason']].append(r)
    coverage=[]
    for name,rs in [('all_selected',master),('conservative_included',[r for r in master if r['conservative_family_included']])]+list(buckets.items()):
        coverage.append(dict(group=name,total_selected=len(rs),any_inferred=sum(r['any_inferred'] for r in rs),candidate_positive=sum(r['combined_status']=='candidate' for r in rs),fully_analyzed_negative=sum(r['combined_status']=='not_selected' for r in rs),unresolved=sum(r['combined_status']=='unresolved' for r in rs)))
    write(OUT/'tables/taxonomy_coverage.tsv',coverage)
    brief=dict(selected=466,candidate_genomes=245,candidate_loci=288,status_counts=dict(collections.Counter(r['combined_status'] for r in master)),support_counts=dict(collections.Counter(r['window_support'] for r in master if r['combined_status']=='candidate')),same_locus_both_genomes=sum(r['same_locus_both_sizes'] for r in master),same_locus_both_loci=sum(r['n_same_locus_both'] for r in master),both_sizes_different_loci_only=sum(r['both_sizes_different_loci_only'] for r in master),classifier_positive_genomes=sum(r['classifier_positive_either_size'] for r in master),strict_retrieval_only_genomes=sum(r['retrieval_only_candidate'] for r in master),generic_taxid2100421=dict(total=sum(r['source_ncbi_taxid']=='2100421' for r in master),positive=sum(r['source_ncbi_taxid']=='2100421' and r['combined_status']=='candidate' for r in master)),coverage=coverage,breadth={k:dict(total_taxa=len(v),positive_taxa=sum(r['candidate_positive']>0 for r in v),selected=sum(r['total_selected'] for r in v),positive=sum(r['candidate_positive'] for r in v)) for k,v in allstats.items()},taxonomy_provenance=provenance)
    (OUT/'summary.json').write_text(json.dumps(brief,indent=2))
    (OUT/'sources/taxonomy/trees.json').write_text(json.dumps(dict(nodes=nodes,trees=trees),indent=2))
    verify_sources()
    (OUT/'validation.json').write_text(json.dumps(dict(passed=True,n_checks=len(checks),checks=checks),indent=2))
    print(json.dumps(brief,indent=2))

if __name__=='__main__': build()
