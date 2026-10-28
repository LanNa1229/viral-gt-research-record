#!/usr/bin/env python3
"""Independent read-only validation of packaged original results. No graph optimization."""
import argparse,csv,hashlib,json,math
from pathlib import Path
from collections import Counter,defaultdict
def table(p,sep=','):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter=sep))
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def verify(data):
    wf=Path(__file__).resolve().parents[1];targets=json.loads((wf/'configs/verification_targets.json').read_text());checks=[]
    def check(name,condition,observed):checks.append(dict(test=name,status='PASS' if condition else 'FAIL',observed=observed))
    n=table(data/'outputs/primary/leiden/cytoscape_nodes_analyzed.csv');ni={r['name']:r for r in n};e=table(data/'outputs/primary/leiden/edges_clean_undirected.csv');p=table(data/'outputs/primary/placement/viral_structural_placement.tsv','\t');alln=table(data/'outputs/primary/placement/cytoscape_nodes_with_figure2b_classes.tsv','\t');families={r['name']:r['cazy_family_normalized'] for r in alln}
    counts=Counter(r['community'] for r in n);incident=Counter();adj=defaultdict(list)
    for r in e:
        a,b=r['source'],r['target'];w=float(r['weight_global_min_tm']);incident.update([a,b]);adj[a].append((b,w));adj[b].append((a,w))
    structure=table(data/'inputs/structures/structure_source_manifest.tsv','\t');display=table(data/'figures/figure2a/historical_source/figure2a_edges_display_backbone.csv')
    metrics=dict(input_structures=len(structure),unique_structural_pairs=len({(r['source'],r['target']) for r in e}),multi_protein_communities=sum(c>1 for c in counts.values()),network_isolates=sum(incident[x]==0 for x in ni),viral_proteins=sum(r['is_virus'].lower()=='true' for r in n),cazy_unassigned_viral=sum(r['cazy_status']=='CAZy_unassigned' for r in p),direct_family_supported=sum(bool(r['directly_supported_cazy_family']) for r in p),dominant_family_concordant=sum(r['unassigned_family_support_call']=='direct_support_concordant_with_community' for r in p),alternative_family_support=sum(r['unassigned_family_support_call']=='direct_support_discordant_with_community' for r in p),figure2a_display_edges=len(display))
    for k,value in metrics.items():check(k,value==targets[k],value)
    check('structure_manifest_and_node_ID_sets',len(ni)==len(n) and {r['protein_id'] for r in structure}==set(ni),len(ni))
    check('viral_placement_ID_set',len(p)==len({r['name'] for r in p}) and {r['name'] for r in p}=={r['name'] for r in n if r['is_virus'].lower()=='true'},len(p))
    pair=defaultdict(list);widths=Counter();selfhits=0;orient=Counter();rawids=set()
    with (data/'inputs/foldseek/005_GTA_all_vs_all_TM065_cov050_rich.id_corrected.m8').open() as f:
        for line in f:
            x=line.rstrip('\r\n').split('\t');widths[len(x)]+=1
            if len(x)!=17:continue
            a,b=x[:2];rawids.update([a,b])
            if a==b:selfhits+=1;continue
            orient[a,b]+=1;pair[tuple(sorted((a,b)))].append([float(v) for v in x[2:11]])
    check('corrected_m8_schema',set(widths)=={17},dict(widths));check('corrected_alignment_IDs',rawids==set(ni),len(rawids));check('raw_pair_IDs',set(pair)=={(r['source'],r['target']) for r in e},len(pair))
    bad=Counter()
    for r in e:
        rows=pair[r['source'],r['target']];den=len(rows);expected={'weight_global_min_tm':sum(min(x[5],x[6]) for x in rows)/den,'weight_alntmscore':sum(x[4] for x in rows)/den,'evalue':min(x[2] for x in rows),'bits_mean':sum(x[3] for x in rows)/den,'bits_max':max(x[3] for x in rows),'mutual_coverage_mean':sum(min(x[7],x[8]) for x in rows)/den,'mutual_coverage_min':min(min(x[7],x[8]) for x in rows),'fident_mean':sum(x[0] for x in rows)/den,'alnlen_mean':sum(x[1] for x in rows)/den,'directional_rows':den}
        for key,value in expected.items():
            if not math.isclose(float(r[key]),value,rel_tol=1e-10,abs_tol=1e-12):bad[key]+=1
    check('all_edge_aggregates_from_alignments',not bad,dict(bad));support_bad=[]
    for r in p:
        fs=defaultdict(float);fc=Counter()
        for b,w in adj[r['name']]:
            if ni[b]['community']==r['community'] and families[b]:fs[families[b]]+=w;fc[families[b]]+=1
        ordered=sorted(fs,key=lambda f:(-fs[f],f));best=ordered[0] if ordered else '';total=sum(fs.values());next_strength=fs[ordered[1]] if len(ordered)>1 else 0
        passes=bool(best) and fc[best]>=2 and fs[best]/total>=.67 and (fs[best]-next_strength)/total>=.20
        expected=best if passes and r['cazy_status']=='CAZy_unassigned' else ''
        if expected!=r['directly_supported_cazy_family']:support_bad.append(r['name'])
    check('direct_support_independently_from_full_graph',not support_bad,support_bad)
    check('classifier_original_node_hash',digest(data/'outputs/primary/leiden/cytoscape_nodes_analyzed.csv')==json.loads((data/'outputs/primary/placement/classification_run.json').read_text())['inputs']['nodes_sha256'],digest(data/'outputs/primary/leiden/cytoscape_nodes_analyzed.csv'))
    check('classifier_original_edge_hash',digest(data/'outputs/primary/leiden/edges_clean_undirected.csv')==json.loads((data/'outputs/primary/placement/classification_run.json').read_text())['inputs']['edges_sha256'],digest(data/'outputs/primary/leiden/edges_clean_undirected.csv'))
    check('Figure2B_original_input_byte_identity',digest(data/'figures/figure2b/viral_structural_placement_input.tsv')==digest(data/'outputs/primary/placement/viral_structural_placement.tsv'),digest(data/'figures/figure2b/viral_structural_placement_input.tsv'))
    labels={'established_family_context':'Established-family context','peripheral_established_family_context':'Family-peripheral context','genuine_mixed_family_community':'Mixed-family community','intercommunity_bridge':'Intercommunity connector','sparsely_annotated_community':'Sparsely annotated community','annotation_depleted_community':'Annotation-depleted community','virus_enriched_dark_community':'Virus-enriched dark community','singleton_or_isolate':'Singleton or isolate'}
    flows=Counter()
    for r in p:
        if r['cazy_status']=='CAZy_assigned':outcome='Concordant known-family placement' if r['cazy_relation_to_community']=='concordant' else 'Discordant assigned-family context'
        else:outcome={'direct_support_concordant_with_community':'Directly supported unassigned protein','no_direct_annotated_neighbor_support':'No direct annotated-neighbour support'}.get(r['unassigned_family_support_call'],'Limited, ambiguous, or alternative support')
        flows['CAZy assigned' if r['cazy_status']=='CAZy_assigned' else 'CAZy unassigned',labels[r['placement_major']],outcome]+=1
    ref=table(data/'figures/figure2b/figure2b_alluvial_flow_counts.tsv','\t');check('all_Figure2B_three_stage_flows',dict(flows)=={(r['annotation'],r['placement'],r['outcome']):int(r['n']) for r in ref},dict(groups=len(flows),total=sum(flows.values())) )
    return dict(scope='packaged original outputs independently verified; not a new Leiden run',metrics=metrics,alignment_observed=dict(rows=sum(widths.values()),self_hits=selfhits,reciprocal_pairs=sum((b,a) in orient for a,b in orient)//2,one_way_pairs=sum((b,a) not in orient for a,b in orient),duplicate_same_direction_rows=sum(c-1 for c in orient.values())),checks=checks,passed=all(c['status']=='PASS' for c in checks))
def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--data-root',required=True,type=Path);a.add_argument('--report',required=True,type=Path);x=a.parse_args();r=verify(x.data_root);x.report.parent.mkdir(parents=True,exist_ok=True);x.report.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['metrics'],indent=2));print('PASS' if r['passed'] else 'FAIL');raise SystemExit(0 if r['passed'] else 1)
if __name__=='__main__':main()
