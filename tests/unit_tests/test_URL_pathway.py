import unittest
from unittest.mock import patch
from src.ResPathExplorer.URL_pathway import get_url_pathway

class TestGetUrlPathway(unittest.TestCase):

    def test_valid_input_without_gene_colors(self):
        with patch("src.ResPathExplorer.URL_pathway.KEGG.show_pathway") as mock_show:
            mock_show.return_value = "http://dummy-url"
            url = get_url_pathway("hsa04110")
            self.assertEqual(url, "http://dummy-url")
            mock_show.assert_called_once()

    def test_valid_input_with_gene_colors(self):
        gene_colors = {
            "TP53": "#FF0000,#000000",
            "BRCA1": "blue,#00FF00"
        }
        with patch("src.ResPathExplorer.URL_pathway.KEGG.show_pathway") as mock_show:
            mock_show.return_value = "http://dummy-colored-url"
            url = get_url_pathway("hsa04110", gene_colors)
            self.assertEqual(url, "http://dummy-colored-url")
            mock_show.assert_called_once()

    def test_invalid_target_path_type(self):
        with self.assertRaises(ValueError):
            get_url_pathway(123)  # not a string

    def test_empty_target_path(self):
        with self.assertRaises(ValueError):
            get_url_pathway("")

    def test_invalid_gene_color_dict_type(self):
        with self.assertRaises(ValueError):
            get_url_pathway("hsa04110", ["TP53", "BRCA1"])

    def test_invalid_gene_color_format_missing_comma(self):
        with self.assertRaises(ValueError):
            get_url_pathway("hsa04110", {"TP53": "#FF0000"})  # missing comma

    def test_invalid_color_in_dict(self):
        with self.assertRaises(ValueError):
            get_url_pathway("hsa04110", {"TP53": "#GGGGGG,#000000"})  # invalid hex code

