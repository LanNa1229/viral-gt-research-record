"""Custom vector-first taxonomic Sankey styled after Nature 2023 Fig.3a.
No inference; all widths are candidate-genome counts; no proprietary font embedding.
"""
import sys
sys.dont_write_bytecode=True
import pathlib,csv,json,hashlib,collections,re,math,logging,datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch,FancyBboxPatch
from matplotlib import font_manager
D=pathlib.Path(__file__).resolve().parents[1];B=D.parents[1];C=B/'itol/candidate_only_v1';M=B/'itol/mixed_rank_v1'
# D is .../figures/nature_style_sankey_v1, so B is the reconciliation directory.
RANKS=['realm','phylum','class','order'];STAGE={'root':0,'realm':1,'phylum':2,'class':3,'order':4,'assignment':5}
WIDTH=180/25.4*72;XS={0:19,1:94,2:170,3:246,4:322,5:398,.75:84};NODE_W=3.8;K=.65
P={'yellow':'#F7D25F','pale_yellow':'#FAF5CA','pink':'#D58DB9','purple':'#AB6DAB','lavender':'#E9D7EA','green':'#88C22C','darkgreen':'#519936','palegreen':'#D2E8BA','midgreen':'#AAD7A2','teal':'#058E6C','cyan':'#2ABBEC','paleblue':'#B3D7F0','ochre':'#995022','paleochre':'#E8D2B9','gray':'#A7A7A7','palegray':'#E4E4E4'}
FONTS={'labels':7.2,'counts':7.0,'rank':7.4,'panel':11.5,'notes':6.5}
MAJORS=['Duplodnaviria','Varidnaviria','Adnaviria','Realm-unassigned Naldaviricetes','Unplaced']
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def gid(prefix,text):return prefix+'_'+hashlib.sha256(text.encode()).hexdigest()[:12]
def major(g):
    p={x['rank']:x['name'] for x in json.loads(g['path_json'])}
    if not p:return 'Unplaced'
    if p.get('realm') in MAJORS:return p['realm']
    assert p.get('class')=='Naldaviricetes' and not p.get('realm'),p
    return 'Realm-unassigned Naldaviricetes'
def load():
    master={r['genome_unit_id']:r for r in read(B/'tables/final_genome_taxonomy_master.tsv')};rows=read(C/'candidate_only_genome_membership.tsv');cm={r['genome_unit_id']:r for r in read(C/'compact_family_oriented/compact_genome_membership.tsv')};ss={r['terminal_id']:r for r in read(C/'compact_family_oriented/candidate_only_terminal_summary.tsv')}
    assert len(rows)==len({r['genome_unit_id'] for r in rows})==245
    for r in rows:
        r['major_lineage']=major(r);r['detailed_terminal_id']=cm[r['genome_unit_id']]['terminal_id'];r['detailed_terminal_label']=ss[r['detailed_terminal_id']]['display_label'];assert master[r['genome_unit_id']]['combined_status']=='candidate'
    return master,rows

def sink_spec(r,version):
    m=r['major_lineage'];family=r['reconciled_family'];genus=r['reconciled_genus'];rank=r['display_rank'];name=r['display_taxon']
    if m=='Unplaced':return 'unplaced','Unplaced\n(display only)','unplaced','display_only',.75
    if version=='detailed':
        if family:return r['detailed_terminal_id'],family,'family','formal_assignment',5
        if genus:return r['detailed_terminal_id'],genus+'‡','genus','formal_assignment',5
        return r['detailed_terminal_id'],name+'\nfamily unresolved†',rank,'direct_unresolved_assignment',5
    if m=='Duplodnaviria':
        if not family and not genus:return 'duplo_class_unresolved','Caudoviricetes\nfamily unresolved†','class','direct_unresolved_assignment',5
        if family=='Kyanoviridae':return 'kyanoviridae','Kyanoviridae','family','formal_assignment',5
        return 'duplo_other','Other accepted taxa*','mixed','display_aggregate',5
    if m=='Varidnaviria':
        if family in ['Mimiviridae','Phycodnaviridae','Pithoviridae']:return family.lower(),family,'family','formal_assignment',5
        if family:return 'nuc_other_families','Other accepted families*','family','display_aggregate',5
        return name.lower()+'_unresolved',name+'\nfamily unresolved†',rank,'direct_unresolved_assignment',5
    if m=='Adnaviria':return 'adnaviria_families','Adnaviria families*','family','display_aggregate',5
    assert m=='Realm-unassigned Naldaviricetes'
    return 'baculoviridae','Baculoviridae','family','formal_assignment',5

def colors(m,rank):
    if m=='Duplodnaviria':return P['pink'] if rank=='realm' else P['purple'],P['lavender']
    if m=='Varidnaviria':return P['green'] if rank=='realm' else P['darkgreen'] if rank in ['phylum','class'] else P['teal'],P['palegreen'] if rank in ['realm','phylum','class'] else P['midgreen']
    if m=='Adnaviria':return P['cyan'],P['paleblue']
    if m=='Realm-unassigned Naldaviricetes':return P['ochre'],P['paleochre']
    return P['gray'],P['palegray']

def graph(rows,version):
    nodes={};edges=collections.defaultdict(list);members=[];paths={}
    def add(key,label,rank,kind,stage,m):
        if key not in nodes:nodes[key]=dict(id=key,label=label,rank=rank,kind=kind,stage=stage,major=m,ids=set())
        return nodes[key]
    add('root_viruses','Viruses','root','cohort_root',0,'root')
    for r in rows:
        code,label,srank,skind,sstage=sink_spec(r,version);sink='assignment_'+code
        path=json.loads(r['path_json']);p={t['rank']:t['name'] for t in path};m=r['major_lineage'];seq=['root_viruses']
        if m!='Unplaced':
            if m=='Realm-unassigned Naldaviricetes':
                key='display_realm_unassigned';add(key,'Realm\nunassigned','realm','missing_rank_display',1,m);seq.append(key)
            for rank in RANKS:
                if p.get(rank):
                    key=gid('tax_'+rank,p[rank]);add(key,p[rank],rank,'formal_taxon',STAGE[rank],m);seq.append(key)
        add(sink,label,srank,skind,sstage,m);seq.append(sink);paths[r['genome_unit_id']]=seq
        for key in seq:nodes[key]['ids'].add(r['genome_unit_id'])
        for a,b in zip(seq,seq[1:]):edges[a,b].append(r['genome_unit_id'])
        members.append(dict(version=version,display_group_id=sink,display_label=label.replace('\n',' / '),display_group_kind=skind,genome_unit_id=r['genome_unit_id'],major_lineage=m,accepted_family=r['reconciled_family'],accepted_genus=r['reconciled_genus'],original_deepest_rank=r['display_rank'],original_deepest_taxon=r['display_taxon'],evidence_category=r['evidence_category'],window_support=r['window_support']))
    for n in nodes.values():n['n']=len(n['ids']);n['source_paths']=[{t['rank']:t['name'] for t in json.loads(r['path_json'])} for r in rows if r['genome_unit_id'] in n['ids']];n['height']=n['n']*K;n['x']=XS[n['stage']];n['node_color'],n['ribbon_color']=colors(n['major'],n['rank'])
    nodes['root_viruses']['node_color']=P['yellow']
    # Family summaries containing one of the accepted predicted families remain marked.
    for n in nodes.values():
        if n['id'].startswith('assignment_') and any(r['genome_unit_id'] in n['ids'] and r['evidence_category']=='High-confidence predicted family' for r in rows):n['label']+='§'
    terminals=[n for n in nodes.values() if n['id'].startswith('assignment_') and n['major']!='Unplaced']
    def termkey(n):
        m=MAJORS.index(n['major']);label=n['label']
        if n['major']=='Duplodnaviria':order=0 if label.startswith('Caudoviricetes') else 1 if label.startswith('Kyanoviridae') else 2 if label.startswith('Other') else 3
        elif n['major']=='Varidnaviria':order=0 if label.startswith('Mimiviridae') else 1 if label.startswith('Other') else 2 if label.startswith('Imitervirales') else 3 if label.startswith('Phycodnaviridae') else 4 if label.startswith('Algavirales') else 5 if label.startswith('Pithoviridae') else 7 if label.startswith('Megaviricetes') else 6
        else:order=0
        if version=='detailed':
            p=n['source_paths'][0];o=p.get('order','');c=p.get('class','')
            if n['major']=='Duplodnaviria':order=0 if label.startswith('Caudoviricetes') else 1 if o=='Pantevenvirales' else 2 if o=='Crassvirales' else 3 if o else 4 if n['rank']=='family' else 5
            elif n['major']=='Varidnaviria':order=0 if o=='Imitervirales' else 1 if o=='Algavirales' else 2 if o=='Pimascovirales' else 3 if c=='Megaviricetes' else 4
        return m,order,-n['n'],label
    terminals.sort(key=termkey);y=23.;lastmajor=None
    for n in terminals:
        if lastmajor is not None and lastmajor!=n['major']:y+=12
        hlabel=FONTS['labels']*1.15*(n['label'].count('\n')+1);space=max(n['height'],hlabel)
        n['center']=y+space/2;n['y']=n['center']-n['height']/2;y+=space+(8 if version=='compact' else 2.5);lastmajor=n['major']
    up=nodes['assignment_unplaced'];up['y']=y+10;up['center']=up['y']+up['height']/2
    targetcenter={g:n['center'] for n in terminals+[up] for g in n['ids']}
    for stage in [4,3,2,1,0]:
        col=[n for n in nodes.values() if n['stage']==stage];col.sort(key=lambda n:sum(targetcenter[g] for g in n['ids'])/n['n']);previous=19
        for n in col:
            desired=sum(targetcenter[g] for g in n['ids'])/n['n'];lab=FONTS['labels']*1.15*(n['label'].count('\n')+1);space=max(n['height'],lab)
            center=max(desired,previous+space/2);n['center']=center;n['y']=center-n['height']/2;previous=center+space/2+12
    bottom=max(n['y']+max(n['height'],(n['label'].count('\n')+1)*FONTS['labels']*1.15) for n in nodes.values());height=bottom+69
    # Port packing: sum of exact count thicknesses equals the source/target bar height.
    out=collections.defaultdict(list);inc=collections.defaultdict(list)
    for (a,b),gids in edges.items():out[a].append(b);inc[b].append(a)
    links=[]
    for (a,b),gids in edges.items():
        na,nb=nodes[a],nodes[b];sortout=sorted(out[a],key=lambda q:nodes[q]['center']);sortin=sorted(inc[b],key=lambda q:nodes[q]['center'])
        ya=na['y']+sum(len(edges[a,q])*K for q in sortout[:sortout.index(b)]);yb=nb['y']+sum(len(edges[q,b])*K for q in sortin[:sortin.index(a)])
        links.append(dict(source=a,target=b,n=len(gids),ids=gids,source_y=ya,target_y=yb,thickness=len(gids)*K))
    return nodes,links,members,paths,height

def render(version,nodes,links,height):
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    for file in ['/mnt/c/Windows/Fonts/arial.ttf','/mnt/c/Windows/Fonts/arialbd.ttf']:
        if pathlib.Path(file).exists():font_manager.fontManager.addfont(file)
    plt.rcParams.update({'font.family':'Arial','font.size':7.2,'svg.fonttype':'none','pdf.use14corefonts':True,'pdf.compression':6,'axes.unicode_minus':False,'savefig.facecolor':'white'})
    fig=plt.figure(figsize=(WIDTH/72,height/72));ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,WIDTH);ax.set_ylim(height,0);ax.axis('off');fig.patch.set_facecolor('white')
    for e in sorted(links,key=lambda e:-e['n']):
        a,b=nodes[e['source']],nodes[e['target']];x0=a['x']+NODE_W;x1=b['x'];y0=e['source_y'];y1=e['target_y'];t=e['thickness'];d=(x1-x0)*.48
        verts=[(x0,y0),(x0+d,y0),(x1-d,y1),(x1,y1),(x1,y1+t),(x1-d,y1+t),(x0+d,y0+t),(x0,y0+t),(x0,y0)];codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY]
        color=b['ribbon_color'] if a['rank']=='root' else a['ribbon_color'];ax.add_patch(PathPatch(Path(verts,codes),facecolor=color,edgecolor='none',alpha=.78,zorder=1))
    texts=[]
    for n in nodes.values():
        rad=min(1.5,n['height']/2);ax.add_patch(FancyBboxPatch((n['x'],n['y']),NODE_W,n['height'],boxstyle='round,pad=0,rounding_size='+str(rad),linewidth=0,facecolor=n['node_color'],zorder=3))
        big=n['height']>=18;y=n['y']+3 if big else n['center'];va='top' if big else 'center'
        t=ax.text(n['x']+NODE_W+3,y,n['label'],fontsize=FONTS['labels'],color='#111111',va=va,ha='left',linespacing=1.1,zorder=5);texts.append((n['id'],'label',t))
        if big:xnum,ynum,hanum,vanum=n['x']-.5,n['y']-3,'left','bottom'
        else:xnum,ynum,hanum,vanum=n['x']-3,n['center'],'right','center'
        t=ax.text(xnum,ynum,str(n['n']),fontsize=FONTS['counts'],color='#111111',weight='bold' if n['n']>=8 else 'normal',va=vanum,ha=hanum,zorder=5);texts.append((n['id'],'count',t))
    ax.text(3,6,'B',weight='bold',fontsize=FONTS['panel'],va='top',color='#111111')
    baseline=height-48;ax.plot([XS[0],XS[5]],[baseline,baseline],color='#111111',lw=.35,zorder=5)
    heads=['Viruses','Realm','Phylum','Class','Order','Family / deepest\nassignment']
    for i,label in enumerate(heads):
        ax.plot([XS[i],XS[i]],[baseline,baseline+3],color='#111111',lw=.35);ax.text(XS[i],baseline+5,label,fontsize=FONTS['rank'],va='top',ha='center',color='#111111',linespacing=1.1)
    notes=['* Display aggregation; † provisional higher-rank placement; ‡ accepted genus without family.', '§ Includes accepted predicted-family assignment. Unplaced is a display category, not a clade.']
    for i,line in enumerate(notes):ax.text(8,height-18+i*8,line,fontsize=FONTS['notes'],color='#333333',va='top')
    fig.canvas.draw();renderer=fig.canvas.get_renderer();boxes=[]
    for nid,kind,t in texts:
        bb=t.get_window_extent(renderer);px=fig.dpi/72;boxes.append(dict(version=version,node_id=nid,text_kind=kind,text=t.get_text(),x0=bb.x0/px,y0=(fig.bbox.height-bb.y1)/px,x1=bb.x1/px,y1=(fig.bbox.height-bb.y0)/px))
    stem=D/('FigureB_Nature_style_'+version)
    fig.savefig(str(stem)+'.svg',metadata={'Creator':'Viral GT custom vector Sankey; design reference Barrio-Hernandez et al. Nature 2023 Fig.3a'})
    fig.savefig(str(stem)+'.pdf',metadata={'Title':'Viral GT Figure B - '+version,'Creator':'Reproducible custom Matplotlib Bezier Sankey','Subject':'Taxonomic assignment of 245 candidate genome units; not gene movement'})
    fig.savefig(str(stem)+'.png',dpi=600);fig.savefig(D/'qc'/('preview_'+version+'.png'),dpi=180);plt.close(fig)
    return boxes

def main():
    (D/'qc').mkdir(exist_ok=True);master,rows=load();top=collections.Counter(r['major_lineage'] for r in rows);assert top==dict(zip(MAJORS,[154,70,3,3,15])),top
    allnodes=[];alllinks=[];allmembers=[];bboxes=[];versionstats={};gidpaths=[]
    for version in ['compact','detailed']:
        nodes,links,members,paths,height=graph(rows,version);assert sum(n['n'] for n in nodes.values() if n['id'].startswith('assignment_'))==245
        assert sum(n['id'].startswith('assignment_') for n in nodes.values())==(13 if version=='compact' else 32)
        bboxes+=render(version,nodes,links,height);allmembers+=members
        for n in nodes.values():allnodes.append(dict(version=version,node_id=n['id'],label=n['label'].replace('\n',' / '),taxonomic_rank=n['rank'],node_kind=n['kind'],stage=n['stage'],major_lineage=n['major'],candidate_count=n['n'],x_pt=n['x'],y_pt=n['y'],width_pt=NODE_W,height_pt=n['height'],node_hex=n['node_color'],ribbon_hex=n['ribbon_color'],genome_unit_ids=';'.join(sorted(n['ids']))))
        for e in links:alllinks.append(dict(version=version,source_node=e['source'],target_node=e['target'],candidate_count=e['n'],thickness_pt=e['thickness'],source_y_pt=e['source_y'],target_y_pt=e['target_y'],genome_unit_ids=';'.join(sorted(e['ids']))))
        gidpaths.extend(dict(version=version,genome_unit_id=g,ordered_node_ids=';'.join(seq)) for g,seq in paths.items());versionstats[version]=dict(width_mm=180,height_mm=height/72*25.4,assignment_groups=sum(n['id'].startswith('assignment_') for n in nodes.values()),root_n=245,unplaced_separate_n=15)
    full=[]
    for r in rows:
        g=dict(master[r['genome_unit_id']]);g.update(plot_validated_path_json=r['path_json'],plot_major_lineage=r['major_lineage'],plot_original_deepest_rank=r['display_rank'],plot_evidence_category=r['evidence_category']);full.append(g)
    write(D/'FigureB_245_candidate_taxonomy_master.tsv',full);write(D/'FigureB_Sankey_nodes.tsv',allnodes);write(D/'FigureB_Sankey_links.tsv',alllinks);write(D/'FigureB_display_group_membership.tsv',allmembers);write(D/'FigureB_genome_flow_paths.tsv',gidpaths);write(D/'qc/text_bounding_boxes.tsv',bboxes)
    style=dict(reference=dict(article='Barrio-Hernandez et al., Nature 2023',doi='10.1038/s41586-023-06510-w',panel='Figure3a, PDF page5',local_pdf='/mnt/c/Users/rahnn/Zotero/storage/4DF72GFH/2023_Nature_Hernandes et al_Clustering predicted structures at the scale of the known protein universe.pdf',artwork_reused_in_output=False),palette=P,font_sizes_pt=FONTS,font_svg_png='Arial, legally installed Windows system font; not distributed',font_pdf='Standard PDF Helvetica / Helvetica-Bold core fonts; no proprietary font program embedded',editable_svg_text=True,node_width_pt=NODE_W,points_per_genome=K,ribbon_opacity=.78,node_corner_radius_pt=1.5,versions=versionstats,top_level_counts=dict(top),flow_semantics='taxonomic assignment; not movement or horizontal gene transfer',missing_ranks='skipped; realm-unassigned carrier is explicitly display-only',unplaced='15-count gray short branch outside named hierarchy; counted as terminal at stage0.75',no_inference_rerun=True)
    (D/'FigureB_style_specification.json').write_text(json.dumps(style,indent=2)+'\n');print(json.dumps(versionstats,indent=2))
if __name__=='__main__':main()


