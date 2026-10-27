import json
import tempfile
import unittest
from pathlib import Path
from viral_gt_mining.retrieve import resolve_anchor


class RetrievalTests(unittest.TestCase):
    def test_archived_SH68_real_record(self):
        root=Path(__file__).resolve().parents[5]
        archived=root/'analysis/viral_glycan_benchmark_A0A0K0N728_20261016/sources/A0A482ICV6_unisave.txt'
        if not archived.exists():self.skipTest('Local archived benchmark record unavailable')
        class NoNetwork:
            def get(self,*args):raise ValueError('Inactive accession')
        meta,links=resolve_anchor('A0A482ICV6',NoNetwork(),{'A0A482ICV6':{'unisave':archived}})
        self.assertEqual(meta['status'],'crosslinks_available')
        self.assertEqual(meta['sequence_sha256'],'7282978941e3e1122c35eb794142f5c27d6097936660d73bcaa7e40c031ff8f9')
        self.assertIn('QBP05883.1',{r['protein_id'] for r in links})
        self.assertIn('MK016663',{r['genome_accession'].split('.')[0] for r in links})

    def test_current_mapping_preserves_all_genomic_crosslinks(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'A.json'
            path.write_text(json.dumps(dict(sequence={'value':'MK'},organism={'taxonId':1},
                uniProtKBCrossReferences=[dict(database='EMBL',id=g,properties=[{'key':'ProteinId','value':p}])
                    for g,p in [('AAA001','AAA123.1'),('BBB002','BBB124.2')]])))
            meta,links=resolve_anchor('A',None,{'A':{'json':path}})
            self.assertEqual(len(links),2)
            self.assertEqual(meta['n_genome_links'],2)


if __name__=='__main__':unittest.main()
