"""Checks frozen scientific correspondence, safeguards and recorded numerical regression."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,unittest,tempfile,importlib.util,json
p=argparse.ArgumentParser();p.add_argument('--cache',type=Path,required=True);p.add_argument('--regression',type=Path,required=True);args,remaining=p.parse_known_args();sys.argv=[sys.argv[0]]+remaining
PACKAGE=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('step2',PACKAGE/'portable/step2.py');step2=importlib.util.module_from_spec(spec);spec.loader.exec_module(step2)
class ReleaseTests(unittest.TestCase):
    def test_complete_frozen_accounting_and_correspondence(self):
        result=step2.verify(args.cache);self.assertEqual([(x['success'],x['saccharide_positive']) for x in result['runs']],[(737,246),(740,265)])
    def test_classifier_rule_against_all_frozen_rows(self):
        labels=['Alkaloid','Terpene','NRP','Polyketide','RiPP','Saccharide','Other']
        for size in ['15','11']:
            for row in step2.rows(args.cache/'runs'/size/'viral_bgc_predictions.tsv'):
                scores=[float(row['probability_'+x]) for x in labels];threshold=max(.5,scores[6]);chosen=[x for x,s in zip(labels,scores) if s>=threshold]
                if 'Other' in chosen and len(chosen)>1:chosen.remove('Other')
                self.assertEqual(row['predicted_labels'],' '.join(chosen) or 'Unknown');self.assertEqual(float(row['threshold']),threshold)
    def test_tampered_resources_fail_closed(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'resource';p.write_bytes(b'tampered')
            with self.assertRaises(ValueError):step2.checked(p,'0'*64)
    def test_existing_mismatched_cache_not_overwritten(self):
        import hashlib
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'resource';p.write_bytes(b'preserve')
            with self.assertRaises(ValueError):step2.put(p,b'new',hashlib.sha256(b'new').hexdigest())
            self.assertEqual(p.read_bytes(),b'preserve')
    def test_existing_regression_output_refused(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(ValueError):step2.regress(args.cache,Path(td),PACKAGE/'tests/representative_windows.json')
    def test_full_inference_requires_explicit_authorization_flag(self):
        with self.assertRaises(ValueError):step2.infer_full(args.cache,Path('/unused_step2_output'),'15',False)
    def test_regression_covers_context_and_cold_encoder(self):
        r=json.loads((args.regression/'regression_receipt.json').read_text());self.assertEqual(r['status'],'PASS');self.assertEqual(r['windows_inferred'],12);self.assertFalse(r['full_cohort_inference']);self.assertFalse(r['historical_pilot_cache_or_LMDB_used']);self.assertTrue(r['all_actual_position_pooling']);self.assertTrue(r['model_weights_unchanged']);self.assertLessEqual(r['max_score_abs_error'],1e-6);self.assertLessEqual(r['max_embedding_abs_error'],1e-6)
if __name__=='__main__':unittest.main()
