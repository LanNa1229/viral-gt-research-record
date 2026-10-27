"""Final 242-genome Sankeys. Reads bundled frozen tables; never runs inference."""
import sys
sys.dont_write_bytecode=True
import pathlib,csv,json,hashlib,collections,importlib.util,datetime,logging,platform
D=pathlib.Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('archived_geometry',D/'scripts/legacy_geometry.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
read,write,sha=v.read,v.write,v.sha
EXCLUDE={'KC511810.1','KX960781.1','KY189087.1'}
SCOPE={'caudo_other_families':('class','Caudoviricetes'),'caudo_genera_no_family':('class','Caudoviricetes'),'imit_other_families':('order','Imitervirales'),'pima_other_families':('order','Pimascovirales')}

def sink(r,version):
    m=r['major_lineage'];f=r['reconciled_family'];g=r['reconciled_genus'];rank=r['display_rank'];name=r['display_taxon']
    if m=='Unplaced':return 'unplaced','Unplaced\n(display only)','unplaced','display_only',.75
    if version=='full':
        if f:return 'family_'+f.lower(),f+('\npredicted placement' if f=='Schitoviridae' else ''),'family','formal_assignment',5
        if g:
            assert r['genus_assignment_status']=='Authoritative' and m=='Duplodnaviria'
            return 'genus_'+g.lower(),g+'‡','genus','formal_assignment',5
        assert name in ['Caudoviricetes','Imitervirales','Algavirales','Megaviricetes']
        return name.lower()+'_unresolved',name+'\nfamily unresolved†',rank,'direct_unresolved_assignment',5
    if m=='Duplodnaviria':
        if not f and not g:return 'caudo_unresolved','Caudoviricetes\nfamily unresolved†','class','direct_unresolved_assignment',5
        if f in ['Kyanoviridae','Darmviridae']:return f.lower(),f,'family','formal_assignment',5
        if f=='Schitoviridae':return 'schitoviridae','Schitoviridae\npredicted placement','family','formal_assignment',5
        if f:
            assert f in ['Kleczkowskaviridae','Stanwilliamsviridae','Straboviridae']
            return 'caudo_other_families','Other families of\nCaudoviricetes (group)','family','display_aggregate',5
        assert g
        return 'caudo_genera_no_family','Genera of Caudoviricetes\nfamily unassigned (group)','genus','display_aggregate',5
    if m=='Varidnaviria':
        if f in ['Mimiviridae','Phycodnaviridae','Pithoviridae']:return f.lower(),f,'family','formal_assignment',5
        if f in ['Allomimiviridae','Mesomimiviridae']:return 'imit_other_families','Imitervirales\nother families (group)','family','display_aggregate',5
        if f in ['Marseilleviridae','Orpheoviridae']:return 'pima_other_families','Pimascovirales\nother families (group)','family','display_aggregate',5
        assert not f and name in ['Imitervirales','Algavirales','Megaviricetes']
        return name.lower()+'_unresolved',name+'\nfamily unresolved†',rank,'direct_unresolved_assignment',5
    if m=='Adnaviria':
        assert f in ['Tristromaviridae','Lipothrixviridae']
        return f.lower(),f,'family','formal_assignment',5
    assert m=='Realm-unassigned Naldaviricetes' and f=='Baculoviridae'
    return 'baculoviridae','Baculoviridae','family','formal_assignment',5

def render(version,nodes,links,height):
    plt=v.plt;Path=v.Path
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    for f in ['/mnt/c/Windows/Fonts/arial.ttf','/mnt/c/Windows/Fonts/arialbd.ttf']:
        if pathlib.Path(f).exists():v.font_manager.fontManager.addfont(f)
    plt.rcParams.update({'font.family':'Arial','font.size':7.2,'svg.fonttype':'none','pdf.use14corefonts':True,'pdf.compression':6,'axes.unicode_minus':False,'savefig.facecolor':'white'})
    fig=plt.figure(figsize=(v.WIDTH/72,height/72));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,v.WIDTH);ax.set_ylim(height,0);ax.axis('off');texts=[]
    for e in sorted(links,key=lambda x:-x['n']):
        a,b=nodes[e['source']],nodes[e['target']];x0=a['x']+v.NODE_W;x1=b['x'];y0=e['source_y'];y1=e['target_y'];t=e['thickness'];d=(x1-x0)*.48
        pts=[(x0,y0),(x0+d,y0),(x1-d,y1),(x1,y1),(x1,y1+t),(x1-d,y1+t),(x0+d,y0+t),(x0,y0+t),(x0,y0)];codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY]
        ax.add_patch(v.PathPatch(Path(pts,codes),facecolor=b['ribbon_color'] if a['rank']=='root' else a['ribbon_color'],edgecolor='none',alpha=.78,zorder=1))
    def text(nid,kind,x,y,label,**kw):
        t=ax.text(x,y,label,**kw);texts.append((nid,kind,t))
    for n in nodes.values():
        ax.add_patch(v.FancyBboxPatch((n['x'],n['y']),v.NODE_W,n['height'],boxstyle='round,pad=0,rounding_size='+str(min(1.5,n['height']/2)),linewidth=0,facecolor=n['node_color'],zorder=3))
        big=n['height']>=18
        text(n['id'],'label',n['x']+v.NODE_W+3,n['y']+3 if big else n['center'],n['label'],fontsize=v.FONTS['labels'],color='#111111',va='top' if big else 'center',ha='left',linespacing=1.1,zorder=5)
        text(n['id'],'count',n['x']-.5 if big else n['x']-3,n['y']-3 if big else n['center'],str(n['n']),fontsize=v.FONTS['counts'],color='#111111',weight='bold' if n['n']>=8 else 'normal',va='bottom' if big else 'center',ha='left' if big else 'right',zorder=5)
    text('panel','panel',3,6,'B',weight='bold',fontsize=v.FONTS['panel'],va='top',color='#111111');baseline=height-48
    ax.plot([v.XS[0],v.XS[5]],[baseline,baseline],color='#111111',lw=.35)
    for i,label in enumerate(['Viruses','Realm','Phylum','Class','Order','Family / deepest\nassignment']):
        ax.plot([v.XS[i],v.XS[i]],[baseline,baseline+3],color='#111111',lw=.35)
        text('rank_'+str(i),'rank',v.XS[i],baseline+5,label,fontsize=v.FONTS['rank'],va='top',ha='center',color='#111111',linespacing=1.1)
    notes=['† Provisional higher-rank placement. Unplaced is a display category, not a clade.','Counts are candidate genome units; Both denotes screening support, not independent replication.'] if version=='compact' else ['† Provisional higher-rank placement. ‡ Accepted genus without assigned family.','Counts are candidate genome units. Unplaced is a display category, not a clade.']
    for i,line in enumerate(notes):text('note_'+str(i),'note',8,height-18+i*8,line,fontsize=v.FONTS['notes'],color='#333333',va='top')
    fig.canvas.draw();rr=fig.canvas.get_renderer();px=fig.dpi/72;boxes=[]
    for nid,kind,t in texts:
        b=t.get_window_extent(rr);boxes.append(dict(node_id=nid,text_kind=kind,text=t.get_text(),x0=b.x0/px,y0=(fig.bbox.height-b.y1)/px,x1=b.x1/px,y1=(fig.bbox.height-b.y0)/px))
    write(D/'qc'/('text_boxes_'+version+'.tsv'),boxes);stem=D/('FigureB_242_'+version)
    fig.savefig(str(stem)+'.svg',metadata={'Creator':'Viral GT 242-genome taxonomic Sankey; established Nature-inspired style'})
    fig.savefig(str(stem)+'.pdf',metadata={'Title':'Viral GT Figure B - 242-genome '+version,'Subject':'Computationally prioritized putative glycan-associated GT neighborhoods; taxonomic assignment, not gene movement'})
    fig.savefig(str(stem)+'.png',dpi=600);fig.savefig(D/'qc'/('preview_'+version+'.png'),dpi=180);plt.close(fig)

def main():
    (D/'qc').mkdir(exist_ok=True)
    master={r['genome_unit_id']:r for r in read(D/'inputs/reconciliation/final_genome_taxonomy_master.tsv')};original=read(D/'inputs/candidate_only/candidate_only_genome_membership.tsv')
    assert len(master)==466 and len(original)==len({r['genome_unit_id'] for r in original})==245
    assert {r['genome_unit_id'] for r in original}=={g for g,r in master.items() if r['combined_status']=='candidate'}
    excluded=[r for r in original if set(r['source_genome_accessions'].split(';'))&EXCLUDE];ids={r['genome_unit_id'] for r in excluded};rows=[dict(r) for r in original if r['genome_unit_id'] not in ids]
    assert len(excluded)==3 and len(rows)==242 and set(';'.join(r['source_genome_accessions'] for r in excluded).split(';'))==EXCLUDE
    audit=read(D/'inputs/incomplete_35_window_pair_audit.tsv');exaudit=[]
    for r in excluded:
        rr=[a for a in audit if a['genome_unit_id']==r['genome_unit_id']];assert len(rr)==1
        a=rr[0];assert a['source_annotated_CDS_count']==a['11_model_input_n']==a['15_model_input_n']=='1' and r['reconciled_family']=='Poxviridae'
        exaudit.append(dict(genome_unit_id=r['genome_unit_id'],source_genome_accessions=r['source_genome_accessions'],organism=r['source_organism_labels'],family=r['reconciled_family'],original_candidate_status=r['combined_status'],original_window_support=r['window_support'],source_length_nt=a['source_record_length_nt'],source_annotated_CDS_count=1,model_input_genes11=1,model_input_genes15=1,anchor_protein=a['anchor_protein_id'],source_product=a['anchor_product'],window15_id=a['window15_id'],source_genbank=a['source_path'],source_sha256=a['source_sha256'],reason='User-authorized figure-cohort exclusion: single-CDS source; multigene neighborhood cannot be evaluated',frozen_status_changed=False))
    write(D/'FigureB_exclusion_audit.tsv',exaudit)
    cm={r['genome_unit_id']:r for r in read(D/'inputs/candidate_only/compact_family_oriented/compact_genome_membership.tsv')};ss={r['terminal_id']:r for r in read(D/'inputs/candidate_only/compact_family_oriented/candidate_only_terminal_summary.tsv')}
    for r in rows:
        r['major_lineage']=v.major(r);r['detailed_terminal_id']=cm[r['genome_unit_id']]['terminal_id'];r['detailed_terminal_label']=ss[r['detailed_terminal_id']]['display_label']
    top=dict(collections.Counter(r['major_lineage'] for r in rows));support=dict(collections.Counter(r['window_support'] for r in rows))
    assert top=={'Duplodnaviria':154,'Varidnaviria':67,'Adnaviria':3,'Realm-unassigned Naldaviricetes':3,'Unplaced':15},top
    assert support=={'Both':187,'11-only':29,'15-only':26},support
    allmembers={};stats={};groups_out=[]
    for version,expected in [('compact',18),('full',31)]:
        v.sink_spec=lambda r,unused:sink(r,version)
        nodes,links,members,paths,height=v.graph(rows,'detailed')
        for n in nodes.values():n['label']=n['label'].replace('§','')
        terminal={n:r for n,r in nodes.items() if n.startswith('assignment_')};assert len(terminal)==expected and sum(n['n'] for n in terminal.values())==242
        aa=[]
        for nid,n in terminal.items():
            rs=[r for r in rows if r['genome_unit_id'] in n['ids']];code=nid[len('assignment_'):]
            scope_rank,scope_taxon=SCOPE.get(code,(n['rank'],rs[0]['reconciled_family'] if n['rank']=='family' else rs[0]['reconciled_genus'] if n['rank']=='genus' else rs[0]['display_taxon']))
            if n['rank']=='unplaced':scope_rank,scope_taxon='unplaced',''
            scope_paths=[]
            for r in rs:
                p=json.loads(r['path_json']);pd={x['rank']:x['name'] for x in p}
                if scope_taxon:
                    assert pd.get(scope_rank)==scope_taxon,(nid,pd,scope_rank,scope_taxon)
                    scope_paths.append(p[:next(i for i,t in enumerate(p) if t['rank']==scope_rank)+1])
                else:scope_paths.append([])
            assert all(p==scope_paths[0] for p in scope_paths),nid
            ev=collections.Counter(r['evidence_category'] for r in rs)
            aa.append(dict(version=version,display_group_id=nid,display_name=n['label'].replace('\n',' / '),candidate_count=n['n'],display_rank=n['rank'],group_kind=n['kind'],scope_rank=scope_rank,scope_taxon=scope_taxon,verified_scope_lineage_json=json.dumps(scope_paths[0]),constituent_family_counts_json=json.dumps(collections.Counter(r['reconciled_family'] or 'UNASSIGNED' for r in rs),sort_keys=True),constituent_genus_counts_json=json.dumps(collections.Counter(r['reconciled_genus'] or 'UNASSIGNED' for r in rs),sort_keys=True),evidence_counts_json=json.dumps(ev,sort_keys=True),genome_unit_ids=';'.join(sorted(n['ids']))))
        write(D/('FigureB_'+version+'_aggregation_audit.tsv'),aa);groups_out+=aa
        for m in members:m['version']=version;m['display_label']=nodes[m['display_group_id']]['label'].replace('\n',' / ')
        write(D/('FigureB_'+version+'_group_membership.tsv'),members);allmembers[version]={r['genome_unit_id']:r for r in members}
        write(D/('FigureB_'+version+'_nodes.tsv'),[dict(node_id=n['id'],label=n['label'].replace('\n',' / '),taxonomic_rank=n['rank'],node_kind=n['kind'],stage=n['stage'],major_lineage=n['major'],candidate_count=n['n'],x_pt=n['x'],y_pt=n['y'],width_pt=v.NODE_W,height_pt=n['height'],node_hex=n['node_color'],ribbon_hex=n['ribbon_color'],genome_unit_ids=';'.join(sorted(n['ids']))) for n in nodes.values()])
        write(D/('FigureB_'+version+'_links.tsv'),[dict(source_node=e['source'],target_node=e['target'],candidate_count=e['n'],thickness_pt=e['thickness'],source_y_pt=e['source_y'],target_y_pt=e['target_y'],genome_unit_ids=';'.join(sorted(e['ids']))) for e in links])
        write(D/('FigureB_'+version+'_genome_flow_paths.tsv'),[dict(genome_unit_id=g,ordered_node_ids=';'.join(p)) for g,p in paths.items()])
        stats[version]=dict(width_mm=180,height_mm=height/72*25.4,width_pt=v.WIDTH,height_pt=height,terminal_categories=expected,terminal_ranks=dict(collections.Counter(n['rank'] for n in terminal.values())))
        # Explicit ancestry assertions and exact membership audit gate precede rendering.
        render(version,nodes,links,height)
    final=[]
    for r in rows:
        out=dict(master[r['genome_unit_id']]);out.update(plot_validated_path_json=r['path_json'],plot_major_lineage=r['major_lineage'],plot_evidence_category=r['evidence_category'])
        for version in stats:out[version+'_terminal_id']=allmembers[version][r['genome_unit_id']]['display_group_id'];out[version+'_terminal_label']=allmembers[version][r['genome_unit_id']]['display_label']
        final.append(out)
    write(D/'FigureB_242_candidate_membership.tsv',final)
    style=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_candidates=245,excluded_genome_units=3,final_genome_units=242,versions=stats,top_level_counts=top,screening_support=support,palette=v.P,font_sizes_pt=v.FONTS,node_width_pt=v.NODE_W,points_per_genome=v.K,ribbon_opacity=.78,png_dpi=600,svg_text_editable=True,pdf_fonts='Standard Helvetica / Helvetica-Bold, vector text',source_typography='Installed Arial for SVG/PNG',no_inference=True,no_frozen_records_modified=True,exclusion_scope='Figure delivery cohort only; original 245 candidate labels untouched',reference='Established local v1/v2 style, originally based on Barrio-Hernandez et al. Nature 2023 Figure 3a',runtime=dict(python=sys.version,platform=platform.platform(),matplotlib=v.matplotlib.__version__))
    (D/'FigureB_style_specification.json').write_text(json.dumps(style,indent=2)+'\n');print(json.dumps(style,indent=2))
if __name__=='__main__':main()
