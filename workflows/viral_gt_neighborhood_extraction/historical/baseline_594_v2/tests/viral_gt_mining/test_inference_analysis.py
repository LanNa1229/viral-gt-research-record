"""Boundary tests for resumability and denominator decisions, without model loading."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from viral_gt_mining.inference import load_window_cache, save_window_cache
from viral_gt_mining.analysis import merge_loci, summarize, unit_status, make_taxonomy, group_genomes


class ResumeTests(unittest.TestCase):
    def test_resume_preserves_vectors_and_rejects_wrong_identity(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as directory:
            path = Path(directory) / "window.npz"
            self.assertIsNone(load_window_cache(path, "key"))
            values = (np.linspace(0, 1, 7, dtype=np.float32), np.ones(320, dtype=np.float32), np.ones(320, dtype=np.float32))
            save_window_cache(path, "key", values)
            for expected, actual in zip(values, load_window_cache(path, "key")):
                np.testing.assert_array_equal(expected, actual)
            with self.assertRaisesRegex(ValueError, "identity mismatch"):
                load_window_cache(path, "other-key")

    def test_resume_rejects_nonfinite_and_wrong_pooling(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as directory:
            path = Path(directory) / "window.npz"
            scores = np.zeros(7, dtype=np.float32)
            vector = np.ones(320, dtype=np.float32)
            scores[0] = np.nan
            save_window_cache(path, "key", (scores, vector, vector))
            with self.assertRaisesRegex(ValueError, "Invalid cached"):
                load_window_cache(path, "key")
            scores[0] = 0
            save_window_cache(path, "key", (scores, vector, vector * 2))
            with self.assertRaisesRegex(ValueError, "pooling mismatch"):
                load_window_cache(path, "key")


class AnalysisTests(unittest.TestCase):
    def test_origin_spanning_CDS_does_not_merge_middle_genome(self):
        import json
        windows=[dict(window_id=wid,genome_accession='CIRC.1',genome_sequence_sha256='same',taxid='1',
            n_genes=1,start_1based=start,end_1based=end,anchor_protein_id=wid,organism='virus')
            for wid,start,end in [('origin',1,1000),('middle',400,500)]]
        genes=[dict(window_id='origin',gene_index=0,start_1based=1,end_1based=1000,
            location_parts=json.dumps([[950,1000,1],[1,40,1]])),
            dict(window_id='middle',gene_index=10,start_1based=400,end_1based=500,location_parts='[[400,500,1]]')]
        loci,units=merge_loci(windows,genes,{})
        self.assertEqual(len(loci),2)
        self.assertIn('[[1, 40], [950, 1000]]',{r['genomic_segments_1based'] for r in loci})

    def test_partial_inference_is_unknown_unless_observed_positive(self):
        predictions = {"a": {"saccharide_predicted": "false"}}
        partial = unit_status(["a", "b"], predictions)
        self.assertEqual(partial["status"], "unknown")
        summary = summarize("genome", [partial])
        self.assertEqual(summary["n_resolved"], 0)
        self.assertEqual(summary["positive_fraction_resolved"], "")
        self.assertEqual(summary["positive_fraction_total_upper_bound"], 1)
        predictions["a"]["saccharide_predicted"] = "true"
        self.assertEqual(unit_status(["a", "b"], predictions)["status"], "positive")

    def test_loci_merge_only_connected_overlap_on_same_genome_identity(self):
        windows, genes = [], []
        for wid, start, end, sequence_hash in [("a", 0, 2, "same"), ("b", 2, 4, "same"), ("c", 5, 6, "same"), ("d", 1, 3, "different")]:
            windows.append(dict(window_id=wid, genome_accession="TEST.1", genome_sequence_sha256=sequence_hash, n_genes=end-start+1, taxid="1", organism="test virus", start_1based=start*100+1, end_1based=end*100+90, anchor_protein_id=wid))
            genes.extend(dict(window_id=wid, gene_index=str(i)) for i in range(start, end+1))
        loci, units = merge_loci(windows, genes, {})
        self.assertEqual(len(loci), 3)
        self.assertIn("a;b", {row["window_ids"] for row in loci})
        self.assertTrue(all(unit["status"] == "unknown" for unit in units))

    def test_species_is_unresolved_without_rank_cache(self):
        window = dict(window_id="a", taxid="1", organism="A virus species-looking name", species="Do not infer this")
        unit = unit_status(["a"], {})
        rows = make_taxonomy([window], [unit], [unit], {}, {})
        self.assertTrue(all(row["taxon"] == "Unresolved" for row in rows))
        self.assertEqual(len(rows), 9)

    def test_identical_genbank_refseq_counterparts_collapse_despite_annotation_index_changes(self):
        windows, genes = [], []
        for wid, accession, indices in [("gb", "GB0001.1", [3, 4]), ("rs", "NC_000001.1", [8, 9])]:
            windows.append(dict(window_id=wid, genome_accession=accession, genome_sequence_sha256="same_whole_genome", taxid="42", organism="test virus", n_genes=2, start_1based=100, end_1based=300, anchor_protein_id="SAME_PROTEIN.1"))
            genes.extend(dict(window_id=wid, gene_index=str(i)) for i in indices)
        genome_units, genome_rows = group_genomes(windows, {})
        self.assertEqual(len(genome_units), 1)
        self.assertEqual(genome_rows[0]["source_genome_accessions"], "GB0001.1;NC_000001.1")
        self.assertEqual(genome_rows[0]["counterpart_alias_count"], 1)
        loci, _ = merge_loci(windows, genes, {})
        self.assertEqual(len(loci), 1)
        self.assertEqual(loci[0]["window_ids"], "gb;rs")
        self.assertEqual(loci[0]["start_gene_index"], "")
        self.assertEqual(loci[0]["start_1based"], 100)

    def test_same_protein_different_whole_genomes_or_taxids_remain_distinct(self):
        windows = [dict(window_id=wid, genome_accession=wid, genome_sequence_sha256=genome_hash, taxid=taxid, organism="test virus", n_genes=1, start_1based=100, end_1based=300, anchor_protein_id="SAME_PROTEIN.1") for wid, genome_hash, taxid in [("a", "genome_a", "42"), ("b", "genome_b", "42"), ("c", "genome_a", "43")]]
        genes = [dict(window_id=w["window_id"], gene_index="0") for w in windows]
        self.assertEqual(len(group_genomes(windows, {})[0]), 3)
        self.assertEqual(len(merge_loci(windows, genes, {})[0]), 3)

    def test_missing_taxid_retains_accession_specific_units(self):
        windows = [dict(window_id=wid, genome_accession=wid, genome_sequence_sha256="same_whole_genome", taxid="") for wid in ["gb", "refseq"]]
        units, rows = group_genomes(windows, {})
        self.assertEqual(len(units), 2)
        self.assertTrue(all(row["counterpart_alias_count"] == 0 for row in rows))


if __name__ == "__main__":
    unittest.main()
