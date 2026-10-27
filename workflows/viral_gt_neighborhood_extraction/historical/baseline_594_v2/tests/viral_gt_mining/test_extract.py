import tempfile
import unittest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, BeforePosition, AfterPosition
from Bio import SeqIO
from viral_gt_mining.extract import parse_record, make_window, extract, read_tsv, AnnotationIndex


def fixture(n=20):
    r = SeqRecord(Seq('ATGAAATAA'+'N'*11)*n,id='TEST001.1',description='test')
    r.annotations.update(molecule_type='DNA',taxonomy=['Viruses'],organism='test virus',topology='linear')
    r.features = [SeqFeature(SimpleLocation(0,len(r)),type='source',qualifiers={'db_xref':['taxon:1']})]
    for i in range(n):
        r.features.append(SeqFeature(SimpleLocation(i*20,i*20+9,strand=1),type='CDS',
            qualifiers={'protein_id':[f'AAA{i:05d}.1'],'translation':['MK']}))
    return r


class ExtractionTests(unittest.TestCase):
    def test_symmetric_window(self):
        info,genes = parse_record(fixture())
        w,g = make_window(info,genes,10)
        self.assertEqual([x['relative_position'] for x in g],list(range(-7,8)))
        self.assertEqual(w['n_genes'],15)
        self.assertEqual(w['inference_eligible'],'true')

    def test_boundaries_no_compensation(self):
        info,genes = parse_record(fixture())
        for i,expected,missing in [(0,8,'-7,-6,-5,-4,-3,-2,-1'),(19,8,'1,2,3,4,5,6,7')]:
            w,g = make_window(info,genes,i)
            self.assertEqual(w['n_genes'],expected)
            self.assertEqual(w['missing_positions'],missing)

    def test_negative_strand_keeps_coordinate_order(self):
        r = fixture()
        r.seq = r.seq[:200]+Seq('TTATTTCAT')+r.seq[209:]
        r.features[11].location = SimpleLocation(200,209,strand=-1)
        info,genes = parse_record(r)
        w,g = make_window(info,genes,10)
        self.assertEqual(w['anchor_strand'],-1)
        self.assertEqual(genes[10]['translation_status'],'verified')
        self.assertEqual([x['gene_index'] for x in g],list(range(3,18)))

    def test_overlapping_anchors_separate(self):
        info,genes = parse_record(fixture())
        a,ag = make_window(info,genes,9)
        b,bg = make_window(info,genes,10)
        self.assertNotEqual(a['window_id'],b['window_id'])
        self.assertEqual(len({g['gene_id'] for g in ag}&{g['gene_id'] for g in bg}),14)

    def test_missing_translation_not_skipped(self):
        r = fixture()
        del r.features[9].qualifiers['translation']
        info,genes = parse_record(r)
        w,g = make_window(info,genes,10)
        self.assertEqual(len(genes),20)
        self.assertEqual(len(g),15)
        self.assertEqual(w['inference_eligible'],'false')
        self.assertEqual(g[5]['protein_id'],'AAA00008.1')

    def test_circular_no_wrap(self):
        r = fixture();r.annotations['topology']='circular'
        info,genes = parse_record(r)
        self.assertEqual(make_window(info,genes,0)[0]['n_genes'],8)

    def test_duplicate_accessions_versions_multiple_occurrences(self):
        with tempfile.TemporaryDirectory() as temp:
            t = Path(temp);r = fixture()
            r.features[3].qualifiers['protein_id']=['AAA00010.1']
            r.features[4].qualifiers['protein_id']=['AAA00010.2']
            SeqIO.write(r,t/'g.gb','genbank')
            (t/'a.txt').write_text('AAA00010.1\nAAA00010.1\nAAA00010\nAAA00010.9\n')
            p = extract(t/'a.txt',None,[t/'g.gb'],t/'out',t/'cache')
            rows = read_tsv(t/'out/anchor_mapping.tsv')
            self.assertEqual(p['successful_input_rows'],3)
            self.assertEqual(p['unresolved_input_rows'],1)
            self.assertEqual(p['windows'],3)
            self.assertEqual(sum(r['input_row']=='1' for r in rows),2)
            self.assertEqual(sum(r['input_row']=='3' for r in rows),3)
            ids = {r['window_id'] for r in rows if r['status']=='mapped'}
            extract(t/'a.txt',None,[t/'g.gb'],t/'out',t/'cache')
            self.assertEqual(ids,{r['window_id'] for r in read_tsv(t/'out/anchor_mapping.tsv') if r['status']=='mapped'})
            self.assertEqual(read_tsv(t/'out/source_manifest.tsv')[0]['cache_reused'],'True')

    def test_multiple_contigs_never_cross(self):
        with tempfile.TemporaryDirectory() as temp:
            t = Path(temp);a=fixture(4);b=fixture(4);b.id='TEST002.1'
            SeqIO.write([a,b],t/'g.gb','genbank')
            (t/'a.txt').write_text('AAA00000.1\n')
            p=extract(t/'a.txt',None,[t/'g.gb'],t/'out',t/'cache')
            self.assertEqual(p['windows'],2)
            self.assertEqual({r['n_genes'] for r in read_tsv(t/'out/viral_gt_windows.tsv')},{'4'})

    def test_mismatch_ineligible(self):
        r=fixture();r.features[11].qualifiers['translation']=['MM']
        info,genes=parse_record(r)
        self.assertEqual(make_window(info,genes,10)[0]['inference_eligible'],'false')

    def test_partial_3prime_uses_only_complete_codons(self):
        r=fixture();r.seq=Seq('ATGAAACC')+r.seq[8:]
        r.features[1].location=SimpleLocation(0,AfterPosition(8),strand=1)
        info,genes=parse_record(r)
        self.assertEqual(genes[0]['translation_status'],'verified_partial')
        self.assertEqual(genes[0]['sequence'],'MK')
        self.assertEqual(make_window(info,genes,0)[0]['inference_eligible'],'true')

    def test_partial_5prime_does_not_invent_methionine(self):
        r=fixture();r.seq=Seq('GAGAAATAA')+r.seq[9:]
        r.features[1].location=SimpleLocation(BeforePosition(0),9,strand=1)
        r.features[1].qualifiers['translation']=['EK']
        info,genes=parse_record(r)
        self.assertEqual(genes[0]['translation_status'],'verified_partial')
        self.assertEqual(genes[0]['sequence'],'EK')

    def test_partial_reverse_strand_biological_end(self):
        r=fixture();r.seq=Seq('TTATTTC')+r.seq[7:]
        r.features[1].location=SimpleLocation(0,AfterPosition(7),strand=-1)
        r.features[1].qualifiers.update(translation=['K'],codon_start=['2'])
        info,genes=parse_record(r)
        self.assertEqual(genes[0]['translation_status'],'verified_partial')

    def test_partial_translation_mismatch_still_excluded(self):
        r=fixture();r.seq=Seq('ATGAAACC')+r.seq[8:]
        r.features[1].location=SimpleLocation(0,AfterPosition(8),strand=1)
        r.features[1].qualifiers['translation']=['MM']
        info,genes=parse_record(r)
        self.assertEqual(genes[0]['translation_status'],'translation_mismatch')

    def test_alias_genome_version_exact(self):
        with tempfile.TemporaryDirectory() as temp:
            t=Path(temp);a=fixture();b=fixture();b.id='TEST001.2'
            SeqIO.write([a,b],t/'g.gb','genbank');(t/'a.txt').write_text('UPI_TEST\n')
            alias={'UPI_TEST':[dict(protein_id='AAA00010.1',genome_accession='TEST001.1')]}
            result=extract(t/'a.txt',None,[t/'g.gb'],t/'out',t/'cache',alias)
            self.assertEqual(result['windows'],1)
            self.assertEqual(read_tsv(t/'out/viral_gt_windows.tsv')[0]['genome_accession'],'TEST001.1')


if __name__=='__main__':
    unittest.main()
