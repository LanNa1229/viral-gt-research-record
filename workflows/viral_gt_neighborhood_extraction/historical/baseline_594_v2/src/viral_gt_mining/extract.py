"""Extract +/-7 CDSs in increasing genomic coordinate order, regardless of strand.

Coordinates are 1-based inclusive. Every annotated CDS counts, including missing
proteins and pseudogenes. Contig ends are never wrapped or compensated by shifting.
"""
import argparse
import csv
import hashlib
import json
import platform
import re
import sqlite3
import time
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
from Bio.Data import CodonTable

PARSER_POLICY='all_CDS_partial_translation_v2'


def digest(text):
    return hashlib.sha256(text.encode()).hexdigest()


def sha_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path):
    with open(path, newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    path, rows = Path(path), list(rows)
    fields = fields or sorted(set().union(*(r.keys() for r in rows)))
    temp = path.with_suffix(path.suffix+'.tmp')
    with temp.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    temp.replace(path)


def json_write(path, value):
    path = Path(path)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2)+'\n')
    temp.replace(path)


def source_translation(feature, record):
    """Verify partial CDSs without demanding absent start/stop codons.

    Only explicitly fuzzy biological ends permit missing end codons. Incomplete
    trailing nucleotides can be discarded only at a documented partial 3' end.
    No source amino-acid residues are added, deleted or replaced.
    """
    q=feature.qualifiers
    offset=int(q.get('codon_start',['1'])[0])-1
    table=int(q.get('transl_table',['1'])[0])
    nt=feature.extract(record.seq)[offset:]
    first,last=feature.location.parts[0],feature.location.parts[-1]
    five=first.end if first.strand==-1 else first.start
    three=last.start if last.strand==-1 else last.end
    partial5=not isinstance(five,ExactPosition)
    partial3=not isinstance(three,ExactPosition)
    if not partial5 and not partial3:
        return str(nt.translate(table=table,cds=True)),'full_CDS'
    if len(nt)%3:
        if not partial3:
            raise ValueError('Incomplete codon at nonpartial biological 3prime end')
        nt=nt[:len(nt)//3*3]
    aa=str(nt.translate(table=table))
    if not partial3 and not aa.endswith('*'):
        raise ValueError('Missing terminal stop at nonpartial biological 3prime end')
    if aa.endswith('*'):aa=aa[:-1]
    if not partial5 and offset==0:
        if str(nt[:3]).upper() not in CodonTable.unambiguous_dna_by_id[table].start_codons:
            raise ValueError('Invalid start at nonpartial biological 5prime end')
        aa='M'+aa[1:]
    return aa,f'partial_CDS_5prime={partial5}_3prime={partial3}_codon_start={offset+1}'


def parse_record(record):
    sources = [f for f in record.features if f.type == 'source']
    taxids = sorted({x.split(':', 1)[1] for f in sources
                    for x in f.qualifiers.get('db_xref', []) if x.startswith('taxon:')})
    lineage = record.annotations.get('taxonomy', [])
    info = dict(genome_accession=record.id, genome_sequence_sha256=digest(str(record.seq).upper()),
        organism=record.annotations.get('organism', ''), taxid=';'.join(taxids), lineage=';'.join(lineage),
        taxonomy_rank_status='ranks_unresolved_source_lineage_only', family='', genus='', species='',
        family_taxid='', genus_taxid='', species_taxid='',
        viral_status='source_lineage_viruses' if 'Viruses' in lineage else 'unverified',
        genome_length=len(record), topology=record.annotations.get('topology', 'unknown'))
    genes = []
    for fi, f in enumerate(record.features):
        if f.type != 'CDS':
            continue
        if f.location is None:
            raise ValueError(f'{record.id}: CDS {fi} has no location')
        q = f.qualifiers
        pids, translations = q.get('protein_id', []), q.get('translation', [])
        seq = ''.join(translations[0].split()).upper() if len(translations) == 1 else ''
        status, detail = 'missing_translation', ''
        if len(translations)>1 or len(pids)>1:
            status = 'ambiguous_qualifiers'
        elif seq:
            status = 'invalid_amino_acids'
            if set(seq) <= set('ACDEFGHIKLMNPQRSTVWYBXZJUO'):
                try:
                    aa,method=source_translation(f,record)
                    status = ('verified_partial' if method!='full_CDS' else 'verified') if aa==seq else 'translation_mismatch'
                    detail=method
                except Exception as exc:
                    status, detail = 'translation_unverified', str(exc)
        if not pids and status in {'verified','verified_partial'}:
            status = 'missing_protein_id'
        if 'pseudo' in q or 'pseudogene' in q:
            status = 'pseudogene'
        parts = [[int(p.start)+1, int(p.end), p.strand] for p in f.location.parts]
        genes.append(dict(protein_id=pids[0] if len(pids)==1 else '', protein_accessions=';'.join(pids),
            source_feature_index=fi, start_1based=int(f.location.start)+1, end_1based=int(f.location.end),
            strand=f.location.strand or 0, location=str(f.location), location_parts=json.dumps(parts),
            partial=any(not isinstance(p.start, ExactPosition) or not isinstance(p.end, ExactPosition) for p in f.location.parts),
            locus_tag=';'.join(q.get('locus_tag', [])), product=';'.join(q.get('product', [])),
            sequence=seq, sequence_sha256=digest(seq) if seq else '', aa_length=len(seq),
            translation_status=status, translation_detail=detail))
    genes.sort(key=lambda g:(g['start_1based'], g['end_1based'], g['source_feature_index']))
    for i, g in enumerate(genes):
        g['gene_index'] = i
        g['gene_id'] = digest(json.dumps([record.id, info['genome_sequence_sha256'], g['location'],
            g['protein_id'], g['source_feature_index']], separators=(',', ':')))[:24]
    info['n_cds'] = len(genes)
    return info, genes


class AnnotationIndex:
    def __init__(self, path):
        self.db = sqlite3.connect(path)
        self.db.executescript('''
        CREATE TABLE IF NOT EXISTS files(path TEXT PRIMARY KEY, sha256 TEXT);
        CREATE TABLE IF NOT EXISTS genomes(key TEXT PRIMARY KEY, path TEXT, data TEXT);
        CREATE TABLE IF NOT EXISTS genes(genome_key TEXT, gene_index INTEGER, protein_id TEXT,
          protein_base TEXT, data TEXT, PRIMARY KEY(genome_key,gene_index));
        CREATE INDEX IF NOT EXISTS protein_exact ON genes(protein_id);
        CREATE INDEX IF NOT EXISTS protein_base ON genes(protein_base);
        ''')

    def add(self, path):
        path, checksum = str(Path(path).resolve()), sha_file(path)
        old = self.db.execute('SELECT sha256 FROM files WHERE path=?', (path,)).fetchone()
        cache_identity=checksum+'|'+PARSER_POLICY
        if old and old[0]==cache_identity:
            return dict(path=path, sha256=checksum, cache_reused=True)
        with open(path) as handle:
            records = list(SeqIO.parse(handle, 'genbank'))
        if not records:
            raise ValueError('No GenBank records: '+path)
        parsed = [parse_record(r) for r in records]
        with self.db:
            self.db.execute('DELETE FROM genes WHERE genome_key IN (SELECT key FROM genomes WHERE path=?)', (path,))
            self.db.execute('DELETE FROM genomes WHERE path=?', (path,))
            for info, genes in parsed:
                key = digest(path+'|'+info['genome_accession'])
                info.update(source_path=path, source_sha256=checksum)
                self.db.execute('INSERT INTO genomes VALUES(?,?,?)', (key, path, json.dumps(info)))
                self.db.executemany('INSERT INTO genes VALUES(?,?,?,?,?)', [
                    (key,g['gene_index'],g['protein_id'],g['protein_id'].split('.')[0],json.dumps(g)) for g in genes])
            self.db.execute('INSERT OR REPLACE INTO files VALUES(?,?)', (path,cache_identity))
        return dict(path=path, sha256=checksum, cache_reused=False)

    def occurrences(self, accession, active_paths):
        field = 'protein_id' if re.search(r'\.\d+$', accession) else 'protein_base'
        matches = self.db.execute(f'SELECT genome_key,gene_index FROM genes WHERE {field}=?', (accession,)).fetchall()
        found, seen = [], set()
        for key, index in matches:
            info = json.loads(self.db.execute('SELECT data FROM genomes WHERE key=?', (key,)).fetchone()[0])
            if info['source_path'] not in active_paths:
                continue
            genes = [json.loads(r[0]) for r in self.db.execute('SELECT data FROM genes WHERE genome_key=? ORDER BY gene_index', (key,))]
            ident = (info['genome_accession'], info['source_sha256'], genes[index]['gene_id'])
            if ident not in seen:
                found.append((info,genes,index))
                seen.add(ident)
        return sorted(found,key=lambda x:(x[0]['genome_accession'],x[2],x[0]['source_path']))


def make_window(info, genes, anchor_index):
    anchor = genes[anchor_index]
    selected = genes[max(0,anchor_index-7):anchor_index+8]
    wid = 'vgt_'+digest(json.dumps([info['genome_accession'],info['genome_sequence_sha256'],
        anchor['gene_id'],7,'coordinate_no_wrap_v1'],separators=(',', ':')))[:24]
    missing = [i for i in range(-7,8) if not 0<=anchor_index+i<len(genes)]
    reasons = sorted({g['translation_status'] for g in selected if g['translation_status'] not in {'verified','verified_partial'}})
    if info['viral_status']=='unverified':
        reasons.append('viral_lineage_unverified')
    row = dict(info, window_id=wid, anchor_protein_id=anchor['protein_id'], anchor_gene_id=anchor['gene_id'],
        anchor_gene_index=anchor_index, anchor_start_1based=anchor['start_1based'],
        anchor_end_1based=anchor['end_1based'], anchor_strand=anchor['strand'], n_genes=len(selected),
        missing_positions=','.join(map(str,missing)), missing_upstream=sum(i<0 for i in missing),
        missing_downstream=sum(i>0 for i in missing), start_1based=min(g['start_1based'] for g in selected),
        end_1based=max(g['end_1based'] for g in selected), inference_eligible=str(not reasons).lower(),
        inference_exclusion_reason=';'.join(reasons))
    members = [dict(g,window_id=wid,genome_accession=info['genome_accession'],
        relative_position=g['gene_index']-anchor_index,is_anchor=g['gene_index']==anchor_index) for g in selected]
    return row,members


def load_anchors(path,column=None):
    if column:
        return [r[column].strip() for r in read_tsv(path)]
    return [s.strip() for s in Path(path).read_text(encoding='utf-8-sig').splitlines() if s.strip() and not s.startswith('#')]


def extract(anchors_path,column,genbanks,run,cache,aliases=None):
    started = time.perf_counter()
    run,cache = Path(run),Path(cache)
    run.mkdir(parents=True,exist_ok=True)
    cache.mkdir(parents=True,exist_ok=True)
    anchors = load_anchors(anchors_path,column)
    if not anchors:
        raise ValueError('Empty accession cohort')
    aliases = aliases or {}
    index = AnnotationIndex(cache/'annotations.sqlite')
    sources,errors,active = [],[],set()
    for path in sorted({str(Path(p).resolve()) for p in genbanks}):
        try:
            sources.append(index.add(path))
            active.add(path)
        except Exception as exc:
            errors.append(dict(path=path,error=str(exc)))
    counts,first = Counter(anchors),{}
    mappings,windows,members = [],{},{}
    for i,accession in enumerate(anchors,1):
        first.setdefault(accession,i)
        targets = aliases.get(accession,[accession])
        found, seen = [],set()
        for target in targets:
            target_id = target if isinstance(target,str) else target['protein_id']
            for info,genes,pos in index.occurrences(target_id,active):
                if isinstance(target,dict):
                    if target.get('sequence_sha256') and genes[pos]['sequence_sha256']!=target['sequence_sha256']:
                        continue
                    requested_genome=target.get('genome_accession','')
                    if requested_genome:
                        genome_matches=(info['genome_accession']==requested_genome if re.search(r'\.\d+$',requested_genome)
                            else info['genome_accession'].split('.')[0]==requested_genome)
                        if not genome_matches:continue
                key = (info['genome_accession'],genes[pos]['gene_id'])
                if key not in seen:
                    found.append((info,genes,pos));seen.add(key)
        base = dict(input_row=i,input_accession=accession,duplicate_count=counts[accession],
            duplicate_of_input_row=first[accession] if first[accession]!=i else '',match_count=len(found),
            version_policy='exact_version_or_all_unversioned_occurrences',
            mapping_scope='supplied_annotations_and_explicit_crosswalks; not exhaustive_public_genomes')
        if not found:
            mappings.append(dict(base,status='unresolved',reason='no_matching_CDS_accession_or_explicit_alias'))
        for info,genes,pos in found:
            row,selected = make_window(info,genes,pos)
            windows[row['window_id']],members[row['window_id']] = row,selected
            mappings.append(dict(base,status='mapped',reason='',window_id=row['window_id'],
                genome_accession=info['genome_accession'],protein_id=genes[pos]['protein_id'],gene_index=pos,
                start_1based=genes[pos]['start_1based'],end_1based=genes[pos]['end_1based'],strand=genes[pos]['strand'],
                source_path=info['source_path'],source_sha256=info['source_sha256'],
                mapping_method='explicit_alias' if accession in aliases else 'protein_accession'))
    rows = sorted(windows.values(),key=lambda r:r['window_id'])
    gene_rows = [g for row in rows for g in members[row['window_id']]]
    qc = [dict(window_id=r['window_id'],status='PASS' if r['inference_eligible']=='true' else 'EXCLUDED_INFERENCE',
        n_genes=r['n_genes'],expected_genes=15,missing_positions=r['missing_positions'],all_CDS_retained=True,
        ordered=True,inference_eligible=r['inference_eligible'],reason=r['inference_exclusion_reason']) for r in rows]
    for name,data,empty in [('anchor_mapping',mappings,['input_row','input_accession','status','reason','window_id']),
        ('viral_gt_windows',rows,['window_id','genome_accession','anchor_protein_id','n_genes','inference_eligible']),
        ('viral_gt_window_genes',gene_rows,['window_id','relative_position','gene_index','protein_id','sequence','sequence_sha256','translation_status']),
        ('extraction_qc',qc,['window_id','status','reason'])]:
        write_tsv(run/(name+'.tsv'),data,None if data else empty)
    write_tsv(run/'source_errors.tsv',errors,['path','error'])
    write_tsv(run/'source_manifest.tsv',sources,['path','sha256','cache_reused'])
    with (run/'viral_gt_window_proteins.faa').open('w') as f:
        for g in gene_rows:
            if g['sequence']:
                f.write(f">{g['window_id']}|{int(g['relative_position']):+d}|{g['gene_id']} protein_id={g['protein_id']}\n{g['sequence']}\n")
    mapped = {r['input_row'] for r in mappings if r['status']=='mapped'}
    result = dict(status='PASS' if not errors else 'SOURCE_ERRORS',input_path=str(Path(anchors_path).resolve()),
        input_sha256=sha_file(anchors_path),input_rows=len(anchors),unique_accessions=len(counts),
        duplicate_rows=len(anchors)-len(counts),successful_input_rows=len(mapped),unresolved_input_rows=len(anchors)-len(mapped),
        mapping_occurrence_rows=sum(r['status']=='mapped' for r in mappings),windows=len(rows),gene_rows=len(gene_rows),
        complete_windows=sum(r['n_genes']==15 for r in rows),eligible_windows=sum(r['inference_eligible']=='true' for r in rows),
        source_files=len(sources),source_errors=len(errors),wall_seconds=time.perf_counter()-started,
        python=platform.python_version(),script_sha256=sha_file(__file__),software_version='0.1.0',parser_policy=PARSER_POLICY,
        parameters=dict(upstream=7,downstream=7,order='ascending_coordinates',circular_policy='no_wrap',
            CDS_policy='all_annotated_CDS',duplicates='retain_input_rows_share_identical_window'))
    json_write(run/'extraction_provenance.json',result)
    index.db.close()
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--anchors',required=True)
    p.add_argument('--accession-column')
    p.add_argument('--genbank',action='append',default=[])
    p.add_argument('--genbank-manifest',help='TSV path column, relative to manifest directory')
    p.add_argument('--alias-manifest',help='TSV input_accession,protein_id,evidence_path,evidence_sha256')
    p.add_argument('--run',required=True)
    p.add_argument('--cache',required=True)
    a = p.parse_args()
    paths = a.genbank
    if a.genbank_manifest:
        paths += [str(Path(a.genbank_manifest).resolve().parent/r['path']) for r in read_tsv(a.genbank_manifest)]
    aliases = {}
    if a.alias_manifest:
        for r in read_tsv(a.alias_manifest):
            evidence = Path(a.alias_manifest).resolve().parent/r['evidence_path']
            if sha_file(evidence)!=r['evidence_sha256']:
                raise ValueError('Alias evidence checksum mismatch: '+str(evidence))
            aliases.setdefault(r['input_accession'],[]).append(r)
    print(json.dumps(extract(a.anchors,a.accession_column,paths,a.run,a.cache,aliases),indent=2))


if __name__=='__main__':
    main()
