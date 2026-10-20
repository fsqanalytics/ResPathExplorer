import unittest
import os
import pandas as pd
from tempfile import TemporaryDirectory

from src.ResPathExplorer.save_df_as_html import save_df_as_html


class TestSaveDfAsHtml(unittest.TestCase):

    def setUp(self):
        self.df = pd.DataFrame({
            "Gene": ["TP53", "BRCA1"],
            "Expression": [2.3, 5.6]
        })
        self.title = "My Gene Data"

    def test_html_file_created_with_title_and_table(self):
        with TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_output.html")

            save_df_as_html(self.df, filepath, self.title)

            self.assertTrue(os.path.exists(filepath), "HTML file was not created")

            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()

                self.assertIn(f"<title>{self.title}</title>", content)
                self.assertIn(f"<h1>{self.title}</h1>", content)

                self.assertIn("TP53", content)
                self.assertIn("2.3", content)

    def test_empty_dataframe(self):
        with TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "empty.html")

            empty_df = pd.DataFrame(columns=["Gene", "Expression"])
            save_df_as_html(empty_df, filepath, "Empty Table")

            self.assertTrue(os.path.exists(filepath))

            with open(filepath, "r", encoding="utf-8") as f:
                content = f.read()
                self.assertIn("<title>Empty Table</title>", content)

