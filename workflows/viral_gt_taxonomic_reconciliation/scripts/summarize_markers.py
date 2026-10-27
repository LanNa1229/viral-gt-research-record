"""Report raw geNomad marker evidence and remove selected GT hits analytically."""
import sys
sys.dont_write_bytecode=True
import collections,json,pathlib
from fetch_metadata import OUT,read,write
from taxonomy_mapping import TaxonomyMap
def main():
    db=OUT/'db/genomad_db';parents={};ranks={};names={}
    for line in (db/'nodes.dmp').open():
        f=[x.strip() for x in line.split('|')];parents[f[0]]=f[1];ranks[f[0]]=f[2]
    for line in (db/'names.dmp').open():
        f=[x.strip() for x in line.split('|')]
        if f[3]=='scientific name':names[f[0]]=f[1]
    lineages={}
    def lineage(t):
        if t not in lineages:
            original=t;result={};seen=set()
            while t in parents and t not in seen:
                seen.add(t);result[ranks[t]]=names.get(t,'');t=parents[t]
            lineages[original]=result
            return result
        return lineages[t]
    mapping=TaxonomyMap();old_to_current_families=collections.defaultdict(set)
    for taxid in parents:
        lin=lineage(taxid)
        if lin.get('family') and (lin.get('genus') or lin.get('species')):
            current,_=mapping.from_lineage({k:v for k,v in lin.items() if k in ['species','genus']})
            if current.get('family'):old_to_current_families[lin['family']].add(current['family'])
    write(OUT/'tables/reference_family_MSL41_crosswalk.tsv',[dict(reference_family=k,current_families=';'.join(sorted(v)),ambiguous_split=len(v)>1) for k,v in sorted(old_to_current_families.items())])
    quality={r['genome_unit_id']:r for r in read(OUT/'tables/full_genome_quality.tsv')}
    gt=collections.defaultdict(list)
    for r in read(OUT/'tables/cached_full_genome_proteins.tsv'):
        if r['selected_gt_anchor']=='True':gt[r['genome_unit_id']].append((int(r['start']),int(r['end'])))
    path=OUT/'runs/genomad/unresolved_genomes_annotate/unresolved_genomes_genes.tsv'
    genes=read(path);groups=collections.defaultdict(list);gene_evidence=[]
    for r in genes:
        gid=r['gene'].rsplit('_',1)[0];start,end=int(r['start']),int(r['end'])
        overlap=max([max(0,min(end,b)-max(start,a)+1)/min(end-start+1,b-a+1) for a,b in gt[gid]],default=0)
        isgt=overlap>=.5;lin=lineage(r['taxid']) if r['taxid']!='1' else {}
        current,reason=mapping.from_lineage(lin)
        if reason.startswith('mapped_family') and len(old_to_current_families[lin.get('family','')])>1:
            current={};reason='legacy_family_split_without_lower_rank_support'
        r.update(genome_unit_id=gid,selected_gt_overlap=isgt,gt_overlap_fraction=overlap,lineage_json=json.dumps(lin,sort_keys=True),current_lineage_json=json.dumps(current,sort_keys=True),current_mapping_reason=reason);groups[gid].append(r)
        gene_evidence.append(r)
    write(OUT/'tables/genomad_gene_taxonomy_evidence.tsv',gene_evidence)
    rows=[]
    for gid,q in quality.items():
        if q['classification_pending']!='True':continue
        rs=groups.get(gid,[]);row=dict(genome_unit_id=gid,n_predicted_genes=len(rs),n_marker_hits=sum(r['marker']!='NA' for r in rs),n_taxonomic_markers=sum(r['taxid']!='1' for r in rs),n_virus_hallmarks=sum(int(r['virus_hallmark']) for r in rs),n_selected_gt_genes=sum(r['selected_gt_overlap'] for r in rs),marker_gene_coverage=sum(r['marker']!='NA' for r in rs)/len(rs) if rs else 0,taxonomic_gene_coverage=sum(r['taxid']!='1' for r in rs)/len(rs) if rs else 0,genomad_version='1.13.0',reference_db_version='1.9',reference_taxonomy='ICTV MSL39')
        for mode in ['full','no_gt']:
            sel=[r for r in rs if r['taxid']!='1' and (mode=='full' or not r['selected_gt_overlap'])];total=sum(float(r['bitscore']) for r in sel)
            row[mode+'_n_taxonomic_markers']=len(sel);row[mode+'_total_taxonomic_bitscore']=total
            for rank in ['realm','phylum','class','order','family','genus','current_family','current_genus']:
                weights=collections.Counter();counts=collections.Counter()
                for r in sel:
                    name=json.loads(r['current_lineage_json'] if rank.startswith('current_') else r['lineage_json']).get(rank.replace('current_',''),'')
                    if name:weights[name]+=float(r['bitscore']);counts[name]+=1
                winner=weights.most_common(1)[0] if weights else ('',0)
                row[mode+'_'+rank]=winner[0];row[mode+'_'+rank+'_weight_fraction']=winner[1]/total if total else 0;row[mode+'_'+rank+'_marker_count']=counts[winner[0]]
                row[mode+'_'+rank+'_all_weights_json']=json.dumps(weights,sort_keys=True)
                supporting=[r for r in sel if winner[0] and json.loads(r['current_lineage_json'] if rank.startswith('current_') else r['lineage_json']).get(rank.replace('current_',''),'')==winner[0]]
                row[mode+'_'+rank+'_unique_marker_count']=len({r['marker'] for r in supporting})
                row[mode+'_'+rank+'_marker_span_fraction']=(max(int(r['end']) for r in supporting)-min(int(r['start']) for r in supporting)+1)/int(q['genome_length']) if supporting else 0
            row[mode+'_marker_span_fraction']=(max(int(r['end']) for r in sel)-min(int(r['start']) for r in sel)+1)/int(q['genome_length']) if sel else 0
        row['family_stable_without_gt']=row['full_family']==row['no_gt_family'] and bool(row['full_family'])
        row['gt_taxonomic_bitscore_fraction']=1-row['no_gt_total_taxonomic_bitscore']/row['full_total_taxonomic_bitscore'] if row['full_total_taxonomic_bitscore'] else 0
        rows.append(row)
    write(OUT/'tables/genomad_genome_marker_evidence.tsv',rows)
    print('Summarized',len(rows),'genomes;',sum(bool(r['full_family']) for r in rows),'with family-associated markers (not automatically accepted assignments).')
if __name__=='__main__':main()
