"""Focused edge cases for the portable implementation; no model dependencies."""
import importlib.util,unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location("step4",Path(__file__).resolve().parents[1]/"portable/step4.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class SelectionTests(unittest.TestCase):
    def test_retrieval_boundary_inclusive(self):
        self.assertTrue(m.positive({"saccharide_predicted":"False"},{"saccharide_fraction_top10":"0.5"}))
        self.assertFalse(m.positive({"saccharide_predicted":"False"},{"saccharide_fraction_top10":"0.4"}))
    def test_read_frozen_flag_not_raw_probability(self):
        self.assertFalse(m.positive({"saccharide_predicted":"False","saccharide_probability":"0.9"},{"saccharide_fraction_top10":"0.2"}))
    def test_classifier_alone_sufficient(self):
        self.assertTrue(m.positive({"saccharide_predicted":"True"},{"saccharide_fraction_top10":"0"}))
    def test_unknown_never_silent_negative(self):
        self.assertEqual(m.state(["a","b"],{"a":{}},set()),"unresolved")
        self.assertEqual(m.state(["a","b"],{"a":{},"b":{}},set()),"not_selected")
    def test_positive_wins_over_uninferred_member(self):
        self.assertEqual(m.state(["a","missing"],{"a":{}},{"a"}),"candidate")
    def test_cross_size_union_and_support(self):
        self.assertEqual(m.combined("candidate","unresolved"),"candidate")
        self.assertEqual(m.support("candidate","unresolved"),"15-only")
        self.assertEqual(m.support("candidate","candidate"),"Both")
        self.assertEqual(m.combined("not_selected","unresolved"),"unresolved")
    def test_same_sequence_taxid_aliases(self):
        a={"genome_sequence_sha256":"hash","taxid":"1","genome_accession":"A"}
        b=dict(a,genome_accession="B")
        self.assertEqual(m.genome_id(a),m.genome_id(b))
        self.assertNotEqual(m.genome_id(a),m.genome_id(dict(b,taxid="2")))
    def test_missing_key_stays_accession_specific(self):
        a={"genome_sequence_sha256":"hash","taxid":"","genome_accession":"A"}
        self.assertNotEqual(m.genome_id(a),m.genome_id(dict(a,genome_accession="B")))
    def test_inclusive_overlap_not_adjacency(self):
        self.assertEqual(m.combine_segments([(1,10),(10,20),(21,30)]),[[1,20],[21,30]])
    def test_origin_parts_do_not_fill_gap(self):
        w={"start_1based":"1","end_1based":"1000"}
        g={"location_parts":"[[900,1000,1],[1,50,1]]","start_1based":"1","end_1based":"1000"}
        self.assertEqual(m.window_segments(w,[g]),[[1,50],[900,1000]])
    def test_transitive_locus_and_multiple_anchors(self):
        windows={};genes=[]
        for wid,a,b in [("a",1,10),("b",8,18),("c",16,26),("d",27,35)]:
            windows[wid]={"genome_sequence_sha256":"hash","taxid":"1","genome_accession":"A","n_genes":"1","start_1based":str(a),"end_1based":str(b)}
            genes.append({"window_id":wid,"gene_index":"0","location_parts":"[]","start_1based":str(a),"end_1based":str(b)})
        genomes,loci,mapping=m.reconstruct_groups(windows,genes)
        self.assertEqual(len(genomes),1);self.assertEqual(len(loci),2)
        self.assertEqual(mapping["a"],mapping["c"]);self.assertNotEqual(mapping["c"],mapping["d"])
    def test_one_gene_inputs_not_filtered_by_selection(self):
        self.assertEqual(m.state(["one_gene"],{"one_gene":{"n_genes":"1"}},{"one_gene"}),"candidate")
if __name__=="__main__":unittest.main(verbosity=2)

