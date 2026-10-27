"""Contract tests supplementing full frozen-data regression; no inference or network."""
import hashlib, importlib.util, json, tempfile, unittest
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('step3_adapter', ROOT/'portable/step3.py')
adapter = importlib.util.module_from_spec(spec); spec.loader.exec_module(adapter)

class RetrievalContractTests(unittest.TestCase):
    def test_exact_ties_use_accession_and_preserve_duplicate_entries(self):
        references = np.array([[1.,0.],[1.,0.],[0.,1.]])
        identifiers = np.array(['BGCb','BGCa','BGCc'])
        order, similarity = adapter.rank_one(references,identifiers,np.array([1.,0.]))
        self.assertEqual(identifiers[order].tolist(),['BGCa','BGCb','BGCc'])
        self.assertEqual(len(order),3)
        saccharide = np.array([True,True,False])
        self.assertEqual(float(saccharide[order[:2]].mean()),1.0)

    def test_observed_positive_precedes_missing_window(self):
        self.assertEqual(adapter.status(['a','b'],{'a'},{'a'}),'candidate')

    def test_uninferred_is_not_negative(self):
        self.assertEqual(adapter.status(['a','b'],{'a'},set()),'unresolved')
        self.assertEqual(adapter.status(['a','b'],{'a','b'},set()),'not_selected')

    def test_zero_or_nonfinite_query_rejected(self):
        for query in [np.zeros(2),np.array([np.nan,1.])]:
            with self.assertRaises(ValueError):
                adapter.rank_one(np.eye(2),np.array(['a','b']),query)

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)
            sentinel = path/'preserve.txt'; sentinel.write_text('unchanged')
            with self.assertRaisesRegex(ValueError,'completely new output'):
                adapter.reproduce(path/'missing-cache',path/'missing-data',path)
            self.assertEqual(sentinel.read_text(),'unchanged')

    def test_tampered_cache_fails_before_computation(self):
        with tempfile.TemporaryDirectory() as folder:
            first = adapter.spec()['resources'][0]; path=Path(folder)/first['key']
            path.parent.mkdir(parents=True);path.write_bytes(b'wrong input')
            with self.assertRaisesRegex(ValueError,'changed public input'):
                adapter.verify(Path(folder),None)

    def test_ast_loading_never_executes_historical_module_body(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);path=root/'historical/window_size/paired.py'
            path.parent.mkdir(parents=True)
            source='raise RuntimeError("module body must not execute")\ndef reference_metrics(vectors, metadata, condition):\n    return "function-only"\n'
            path.write_text(source)
            (root/'resources.json').write_text(json.dumps({'historical_sources':[{'file':'historical/window_size/paired.py','sha256':hashlib.sha256(source.encode()).hexdigest()}]}))
            previous=adapter.PACKAGE
            try:
                adapter.PACKAGE=root
                self.assertEqual(adapter.paired_function(root)(None,None,None),'function-only')
            finally:
                adapter.PACKAGE=previous

if __name__=='__main__':
    unittest.main()

