"""Map named reference taxa to MSL41 without guessing across taxon splits."""
from collections import defaultdict
from fetch_metadata import OUT,read
RANKS=['realm','subrealm','kingdom','subkingdom','phylum','subphylum','class','subclass','order','suborder','family','subfamily','genus','subgenus','species']
class TaxonomyMap:
    def __init__(self):
        self.rows=read(OUT/'sources/ictv/ICTV_Master_Species_List_2025_MSL41.v1__MSL.tsv')
        self.index={r:{} for r in RANKS};self.aliases={r:{} for r in RANKS}
        self.byid={x['ICTV_ID']:x for x in self.rows if x.get('ICTV_ID')}
        for row in self.rows:
            for rank in RANKS:
                name=row.get(rank.capitalize(),'')
                if name and name not in self.index[rank]:
                    self.index[rank][name]={k:row.get(k.capitalize(),'') or '' for k in RANKS[:RANKS.index(rank)+1]}
        for row in read(OUT/'sources/ictv/ICTV_Master_Species_List_2025_MSL41.v1__Taxa_Renamed_or_Abolished.tsv'):
            rank=(row.get('Rank') or '').lower();old=(row.get('Old Name') or '').split(';')[-1].strip();new=(row.get('New Name') or '').split(';')[-1].strip()
            if rank not in self.index:continue
            if rank=='species' and row.get('ICVD_ID','') in self.byid:
                new=self.byid[row['ICVD_ID']]['Species']
            if 'Renamed' in (row.get('Change') or '') and new in self.index[rank]:self.aliases[rank][old]=new
    def current(self,rank,name):
        name=str(name or '').strip()
        return name if name in self.index[rank] else self.aliases[rank].get(name,'')
    def from_lineage(self,lineage):
        for rank in reversed(RANKS):
            name=self.current(rank,lineage.get(rank,''))
            if name:return dict(self.index[rank][name]),'mapped_'+rank+('_renamed' if name!=lineage.get(rank) else '')
        return {},'no_current_MSL41_rank'
