import unittest
from unittest.mock import patch, MagicMock
import os
import tempfile
import pandas as pd
from src.ResPathExplorer.CARDAnalysis import CARDAnalysis


class TestCARDAnalysis(unittest.TestCase):

    def setUp(self):
        self.mock_genes = ["geneA", "geneB"]

    def test_invalid_genes_list_raises(self):
        with self.assertRaises(ValueError):
            CARDAnalysis(genes_list="notalist")

    @patch.object(CARDAnalysis, 'download_CARD_file')
    @patch.object(CARDAnalysis, 'finding_ARG', return_value=([{"Gene Name": "geneA"}], []))
    def test_init_with_valid_input(self, mock_find, mock_download):
        obj = CARDAnalysis(genes_list=self.mock_genes, has_CARDdata=True)
        self.assertIsInstance(obj.ARGdf, pd.DataFrame)
        self.assertEqual(len(obj.ARGdf), 1)
        self.assertEqual(obj.ARGdf.iloc[0]["Gene Name"], "geneA")

    def test_find_gene_ids_matches_synonym(self):
        obo_content = """
        [Term]
        id: ARO:1234567
        name: beta-lactamase
        synonym: "blaTEM, blaSHV"
        def: "A beta-lactamase enzyme."
        relationship: confers_resistance_to_antibiotic ARO:3000001 ! Penicillin

        [Term]
        """
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, encoding='utf-8') as f:
            f.write(obo_content)
            tmp_file = f.name

        obj = CARDAnalysis.__new__(CARDAnalysis)
        result = obj.find_gene_ids(tmp_file, "blaTEM")

        self.assertIsNotNone(result)
        self.assertEqual(result["Gene Name"], "beta-lactamase")
        self.assertEqual(result["Matched Name"], "blaTEM")  # <- NOVO TESTE
        self.assertIn("blaTEM", result["All Synonyms"])  # <- NOVO TESTE
        self.assertIn("Penicillin", result["Antibiotics"])

        os.remove(tmp_file)

    def test_find_gene_ids_file_not_found(self):
        obj = CARDAnalysis.__new__(CARDAnalysis)
        with self.assertRaises(FileNotFoundError):
            obj.find_gene_ids("nonexistent.obo", "geneX")

    def test_finding_ARG_separates_found_and_not_found(self):
        obj = CARDAnalysis.__new__(CARDAnalysis)
        obj.find_gene_ids = MagicMock(side_effect=[
            {"Gene Name": "geneA", "Gene ID": "ARO:001", "Description": "desc", "Antibiotics": "Penicillin"},
            None
        ])

        found, not_found = obj.finding_ARG(["geneA", "geneB"], "fake.obo")
        self.assertEqual(len(found), 1)
        self.assertEqual(found[0]["Gene Name"], "geneA")
        self.assertEqual(not_found, ["geneB"])

    def test_add_antibiotic_existing_gene(self):
        obj = CARDAnalysis.__new__(CARDAnalysis)
        obj.ARGdf = pd.DataFrame([{
            "Gene Name": "beta-lactamase",
            "Gene ID": "ARO:1234567",
            "Description": "desc",
            "Antibiotics": pd.NA
        }])
        obj.add_antibiotic_by_id("ARO:1234567", "Penicillin")
        self.assertEqual(obj.ARGdf.loc[0, "Antibiotics"], "Penicillin")

    def test_add_antibiotic_nonexistent_gene(self):
        obj = CARDAnalysis.__new__(CARDAnalysis)
        obj.ARGdf = pd.DataFrame([{
            "Gene Name": "beta-lactamase",
            "Gene ID": "ARO:1234567",
            "Description": "desc",
            "Antibiotics": pd.NA
        }])
        # Captura print como fallback, ou apenas garante que nada muda
        obj.add_antibiotic_by_id("ARO:9999999", "Unknown")
        self.assertTrue(pd.isna(obj.ARGdf.loc[0, "Antibiotics"]))

    def test_plot_antibiotic_frequencies_valid_df(self):
        obj = CARDAnalysis.__new__(CARDAnalysis)
        df = pd.DataFrame([
            {"Gene Name": "geneA", "Gene ID": "ARO:001", "Description": "", "Antibiotics": "Penicillin"},
            {"Gene Name": "geneB", "Gene ID": "ARO:002", "Description": "", "Antibiotics": "Penicillin, Tetracycline"},
            {"Gene Name": "geneC", "Gene ID": "ARO:003", "Description": "", "Antibiotics": "Tetracycline"}
        ])

        with patch("matplotlib.pyplot.show") as mock_show:
            obj.plot_antibiotic_frequencies(df)
            mock_show.assert_called_once()

    def test_plot_antibiotic_frequencies_missing_column_raises(self):
        obj = CARDAnalysis.__new__(CARDAnalysis)
        df_invalid = pd.DataFrame([
            {"Gene Name": "geneX", "Gene ID": "ARO:999", "Description": ""}
        ])

        with self.assertRaises(ValueError):
            obj.plot_antibiotic_frequencies(df_invalid)
