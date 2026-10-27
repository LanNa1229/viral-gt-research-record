"""Representative frozen-data regressions; outputs live only in fresh temporary directories."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import unittest,tempfile,subprocess,json,csv,hashlib,shutil
HOME=Path(__file__).resolve().parents[1];FIX=HOME/'tests/fixture';CLI=HOME/'portable/step1.py'
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class PortableRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory(prefix='viral gt step1 regression ');cls.root=Path(cls.tmp.name);cls.data=cls.root/'relocated dataset with spaces';shutil.copytree(FIX,cls.data);cls.out=cls.root/'fresh outputs'
        result=subprocess.run([sys.executable,'-B',str(CLI),'extract','--data-root',str(cls.data),'--output',str(cls.out)],capture_output=True,text=True)
        if result.returncode:raise RuntimeError(result.stdout+result.stderr)
        cls.execution=json.loads((cls.out/'portable_execution.json').read_text())
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def test_all_window_fields_match_frozen(self):
        for size in ['15','11']:
            observed={r['window_id']:r for r in read(self.out/('window'+size)/'viral_gt_windows.tsv')};expected={r['window_id']:r for r in read(self.data/'expected'/(size+'_viral_gt_windows.tsv'))};self.assertEqual(set(observed),set(expected))
            for wid,r in expected.items():
                for k,v in r.items():
                    if k=='source_path':continue # relocation is explicit and checked independently
                    self.assertEqual(observed[wid].get(k),v,(size,wid,k))
                self.assertTrue(Path(observed[wid]['source_path']).is_relative_to(self.data))
    def test_all_gene_fields_match_frozen(self):
        for size in ['15','11']:
            expected={r['gene_id']+'|'+r['window_id']:r for r in read(self.data/'expected'/(size+'_viral_gt_window_genes.tsv'))};observed={r['gene_id']+'|'+r['window_id']:r for r in read(self.out/('window'+size)/'viral_gt_window_genes.tsv')};self.assertEqual(expected,observed)
    def test_original_mapping_and_unresolved_preserved(self):
        a=read(self.data/'expected/anchor_mapping.tsv');b=read(self.out/'window15/anchor_mapping.tsv');key=lambda r:(r['input_accession'],r['status'],r.get('window_id',''))
        self.assertEqual({key(r) for r in a},{key(r) for r in b});self.assertEqual({r['input_accession'] for r in b if r['status']=='unresolved'},{'UPI001FBC14AE'})
    def test_source_fixture_remains_byte_identical(self):
        for r in read(FIX/'SHA256SUMS.tsv'):self.assertEqual(sha(FIX/r['path']),r['sha256']);self.assertEqual(sha(self.data/r['path']),r['sha256'])
    def test_boundaries_and_input_sizes(self):
        for size,nominal in [('15',15),('11',11)]:
            for r in read(self.out/('window'+size)/'viral_gt_windows.tsv'):
                missing=[x for x in r['missing_positions'].split(',') if x];self.assertEqual(int(r['n_genes'])+len(missing),nominal)
        self.assertTrue(any(int(r['n_genes'])==1 for r in read(self.out/'window15/viral_gt_windows.tsv')))
    def test_structural_grouping_preserves_aliases_and_overlap(self):
        for size in ['15','11']:
            units=read(self.out/('window'+size)/'structural_genome_units.tsv');loci=read(self.out/('window'+size)/'structural_loci.tsv');self.assertEqual(len(units),7);self.assertEqual(len(loci),7)
            self.assertTrue(any(int(r['counterpart_alias_count'])==1 for r in units));self.assertTrue(any(len(json.loads(r['genomic_segments_1based']))>1 for r in loci));self.assertTrue(all('classification_status' not in r for r in loci))
    def test_no_model_or_network_work(self):
        self.assertEqual(self.execution['inference_calls'],0);self.assertEqual(self.execution['network_calls'],0);self.assertFalse(any(self.out.rglob('viral_bgc_predictions.tsv')))
    def test_refuse_existing_output(self):
        marker=self.out/'keep.txt';marker.write_text('preserve');r=subprocess.run([sys.executable,'-B',str(CLI),'extract','--data-root',str(self.data),'--output',str(self.out)],capture_output=True,text=True);self.assertNotEqual(r.returncode,0);self.assertEqual(marker.read_text(),'preserve')
    def test_detect_source_checksum_tamper(self):
        bad=self.root/'damaged fixture';shutil.copytree(self.data,bad);source=read(bad/'portable_inputs/genbanks.tsv')[0]['path'];p=bad/source;p.write_bytes(p.read_bytes()+b'\n');r=subprocess.run([sys.executable,'-B',str(CLI),'extract','--data-root',str(bad),'--output',str(self.root/'not generated')],capture_output=True,text=True);self.assertNotEqual(r.returncode,0);self.assertIn('checksum mismatch',r.stderr);self.assertFalse((self.root/'not generated').exists())
if __name__=='__main__':unittest.main(verbosity=2)
