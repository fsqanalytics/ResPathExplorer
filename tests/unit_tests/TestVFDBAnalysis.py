import unittest
from unittest.mock import patch, mock_open, MagicMock
import pandas as pd
from io import StringIO
from src.ResPathExplorer.VFDBAnalysis import VFDBAnalysis


class TestVFDBAnalysis(unittest.TestCase):

    def setUp(self):
        self.vfdb = VFDBAnalysis(db_dir="test_db")
        self.valid_df = pd.DataFrame({
            "Functional category": [
                "Adherence", "Adherence", "Toxin", "Adherence",
                "Invasion", "Biofilm", "Biofilm", "Regulation"
            ]
        })

    def test_invalid_db_dir_raises(self):
        with self.assertRaises(ValueError):
            VFDBAnalysis(db_dir="")

    @patch("os.makedirs")
    def test_init_creates_dir(self, mock_makedirs):
        vfdb = VFDBAnalysis(db_dir="some_dir")
        mock_makedirs.assert_called_once_with("some_dir", exist_ok=True)

    @patch("os.path.exists", return_value=False)
    @patch("requests.get")
    def test_download_data_downloads_files(self, mock_get, mock_exists):
        mock_response = MagicMock()
        mock_response.content = b"dummy fasta data"
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        with patch("builtins.open", mock_open()) as mock_file, \
             patch("gzip.open", mock_open(read_data=b"mock xls data")), \
             patch("shutil.copyfileobj"):
            self.vfdb.download_data()

        self.assertTrue(mock_get.called)

    @patch("gzip.open")
    @patch("os.path.exists", return_value=True)
    def test_parse_fasta_entries_valid(self, mock_exists, mock_gzip_open):
        mock_fasta_content = """>VFG037176(gb|WP_001081735) (plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii ACICU]
MNRREFLLNSTKTMFGTAALASFPLSIQKALAIDAKVESGTIQDVKHIVILTQENRSFDN
>VFG037177(gb|WP_000632986) (plc2) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii ACICU]
MITRRKFLNYSLNMGFGAAALAAFPSSIQKALAIPANNKTGTIQDVEHVIILMQENRSFD
HYFGTLKGVRGFADRFTIPLPNGRRVWEQLRSNGQVLTPFHLDGTANNAQRADGTPHTWD
DSQLAWDNGRMANWPTHKTDISMGYFKEKEIPYQFALANAFTICDAYHCSMHTGTDANRS
FHLTGTNGATPTKRSFVNNEWDWIDGNPATADRGYTWKTYAERLEEAGISWICYQNMPDE"""

        mock_gzip_open.return_value.__enter__.return_value = StringIO(mock_fasta_content)

        entries = self.vfdb.parse_fasta_entries()
        self.assertEqual(len(entries), 2)
        self.assertEqual(entries[0]["Gene_Name"], "plc1")
        self.assertEqual(entries[1]["VFID"], "VF0470")

    @patch.object(VFDBAnalysis, "download_data")
    @patch.object(VFDBAnalysis, "parse_fasta_entries")
    @patch("pandas.read_excel")
    @patch("os.path.exists", return_value=True)
    def test_load_and_process_merges_correctly(self, mock_exists, mock_read_excel, mock_parse, mock_download):
        # Mocked FASTA and Excel data
        fasta_entries = [
            {'Gene_Name': 'GeneA', 'VFID': '001', 'Description': '', 'Functional category': '', 'Bacteria': 'E.coli'},
            {'Gene_Name': 'GeneB', 'VFID': 'BADID', 'Description': '', 'Functional category': '', 'Bacteria': 'E.coli'}
        ]
        df_xls = pd.DataFrame({
            'VFID': ['001', '002'],
            'VF_Name': ['A', 'B'],
            'Function': ['FuncA', 'FuncB']
        })

        mock_parse.return_value = fasta_entries
        mock_read_excel.return_value = df_xls

        df = self.vfdb.load_and_process()
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 2)

    def test_search_virulence_genes_filters(self):
        self.vfdb.df_genes = pd.DataFrame({
            'Gene_Name': ['GeneA', 'GeneB', 'GeneC'],
            'Bacteria': ['E.coli', 'E.coli', 'Salmonella'],
            'VFID': ['001', '002', '003']
        })

        result = self.vfdb.search_virulence_genes(['GeneA', 'GeneC'], 'E.coli')
        self.assertEqual(len(result), 1)
        self.assertIn('GeneA', result['Gene_Name'].values)

    def test_search_virulence_genes_validations(self):
        self.vfdb.df_genes = pd.DataFrame()

        with self.assertRaises(ValueError):
            self.vfdb.search_virulence_genes([], "E.coli")

        with self.assertRaises(ValueError):
            self.vfdb.search_virulence_genes(['gene'], "")

        self.vfdb.df_genes = None
        with self.assertRaises(ValueError):
            self.vfdb.search_virulence_genes(['gene'], "E.coli")

    def test_plot_raises_on_empty_df(self):
        with self.assertRaises(ValueError):
            self.vfdb.plot_virulence_factors_percentage(pd.DataFrame(), "E.coli")

    def test_plot_raises_on_missing_column(self):
        df = pd.DataFrame({"Invalid": ["x"]})
        with self.assertRaises(ValueError):
            self.vfdb.plot_virulence_factors_percentage(df, "E.coli")

    def test_plot_raises_on_invalid_bacteria_name(self):
        with self.assertRaises(ValueError):
            self.vfdb.plot_virulence_factors_percentage(self.valid_df, "")

    @patch("matplotlib.pyplot.show")
    def test_plot_executes_successfully(self, mock_show):
        try:
            self.vfdb.plot_virulence_factors_percentage(self.valid_df, "E.coli")
        except Exception as e:
            self.fail(f"Plot raised unexpected error: {e}")
        mock_show.assert_called_once()

    @patch("matplotlib.pyplot.show")
    def test_plot_show_all_categories_false(self, mock_show):
        try:
            self.vfdb.plot_virulence_factors_percentage(
                self.valid_df, "E.coli", show_all_categories=False
            )
        except Exception as e:
            self.fail(f"Plot with filtered categories raised error: {e}")
        mock_show.assert_called_once()
