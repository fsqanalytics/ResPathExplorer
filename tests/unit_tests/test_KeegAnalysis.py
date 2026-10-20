import pytest
from unittest.mock import patch, mock_open, MagicMock
from src.ResPathExplorer.KeggAnalysis import KeggAnalysis
import pandas as pd
import matplotlib
matplotlib.use('Agg')

class TestKeggAnalysis:

    def test_init_raises_value_error_on_empty_args(self):
        with pytest.raises(ValueError):
            KeggAnalysis("", "file.gmt")
        with pytest.raises(ValueError):
            KeggAnalysis("hsa", "")

    @patch("os.path.exists", return_value=True)
    @patch("src.ResPathExplorer.KeggAnalysis.KeggAnalysis._load_gmt_file")
    def test_init_loads_existing_gmt(self, mock_load_gmt, mock_exists):
        ka = KeggAnalysis("hsa", "file.gmt", use_existing_gmt=True)
        mock_load_gmt.assert_called_once_with("file.gmt")
        assert ka.file_name_gmt == "file.gmt"
        assert ka.org == "hsa"
        assert ka.organism is not None

    @patch("src.ResPathExplorer.KeggAnalysis.REST.kegg_list")
    def test_get_organism_prefix_valid(self, mock_kegg_list):
        mock_kegg_list.return_value.read.return_value = (
            "T01001\thsa\tHomo sapiens (human)\n"
            "T01002\tmmu\tMus musculus (mouse)\n"
        )
        ka = KeggAnalysis.__new__(KeggAnalysis)
        prefix = ka._get_organism_prefix("Homo sapiens")
        assert prefix == "hsa"

    @patch("src.ResPathExplorer.KeggAnalysis.REST.kegg_list")
    def test_get_organism_prefix_invalid(self, mock_kegg_list):
        mock_kegg_list.return_value.read.return_value = "T01001\thsa\tHomo sapiens (human)\n"
        ka = KeggAnalysis.__new__(KeggAnalysis)
        with pytest.raises(ValueError):
            ka._get_organism_prefix("Unknown organism")

    @patch("src.ResPathExplorer.KeggAnalysis.REST.kegg_list")
    def test_get_organism_name_valid(self, mock_kegg_list):
        mock_kegg_list.return_value.read.return_value = (
            "T01001\thsa\tHomo sapiens (human)\n"
            "T01002\tmmu\tMus musculus (mouse)\n"
        )
        ka = KeggAnalysis.__new__(KeggAnalysis)
        name = ka._get_organism_name("hsa")
        assert name == "Homo sapiens (human)"

    @patch("src.ResPathExplorer.KeggAnalysis.REST.kegg_list")
    def test_get_organism_name_invalid(self, mock_kegg_list):
        mock_kegg_list.return_value.read.return_value = "T01001\thsa\tHomo sapiens (human)\n"
        ka = KeggAnalysis.__new__(KeggAnalysis)
        with pytest.raises(ValueError):
            ka._get_organism_name("xyz")

    def test_load_gmt_file_loads_data_correctly(self):
        fake_gmt_content = "pathway1\tDescription1\tgene1\tgene2\npathway2\tDescription2\tgene3"
        ka = KeggAnalysis.__new__(KeggAnalysis)
        ka.gene_set = {}

        with patch("builtins.open", mock_open(read_data=fake_gmt_content)):
            ka._load_gmt_file("fake_path.gmt")

        assert ka.gene_set == {
            ("pathway1", "Description1"): ["gene1", "gene2"],
            ("pathway2", "Description2"): ["gene3"]
        }

    def test_get_kgml_returns_valid_data(self):
        ka = KeggAnalysis.__new__(KeggAnalysis)

        class FakeService:
            def get(self, pathway_id, file_type):
                return "valid KGML content"

        fake_service = FakeService()
        result = ka.get_kgml("hsa00010", service=fake_service)
        assert "valid KGML content" in result

    def test_get_kgml_raises_on_empty(self):
        ka = KeggAnalysis.__new__(KeggAnalysis)

        class FakeService:
            def get(self, pathway_id, file_type):
                return None

        fake_service = FakeService()
        import pytest
        with pytest.raises(ValueError):
            ka.get_kgml("hsa00010", service=fake_service)

    def test_extract_path_name_from_kgml(self):
        ka = KeggAnalysis.__new__(KeggAnalysis)
        kgml_string = "<pathway name='path:eco01100' org='eco' number='01100' title='Metabolic pathways'></pathway>"

        path_name = ka.extract_path_name_from_kgml(kgml_string)
        assert path_name == "Metabolic pathways"

    def test_extract_genes_from_kgml_string(self):
        ka = KeggAnalysis.__new__(KeggAnalysis)
        kgml_string = """
        <pathway>
            <entry type="gene" name="eco:gene1 eco:gene2"/>
            <entry type="gene" name="eco:gene3"></entry>
        </pathway>
        """
        genes = ka.extract_genes_from_kgml_string(kgml_string)
        assert set(genes) == {"gene1", "gene2", "gene3"}

    def test_create_GMT_file(self, monkeypatch, tmp_path):
        ka = KeggAnalysis.__new__(KeggAnalysis)
        ka.org = "eco"
        ka.gene_set = {}

        # Mock da resposta do requests.get para retornar uma lista de pathways simulada
        fake_response = MagicMock()
        fake_response.raise_for_status = lambda: None
        fake_response.text = "path1\tDescription1\npath2\tDescription2\n"

        monkeypatch.setattr("requests.get", lambda url: fake_response)

        # Mock dos métodos do KeggAnalysis usados dentro de _create_GMT_file
        monkeypatch.setattr(ka, "get_kgml", lambda path, service=None: """
            <pathway>
                <entry type="gene" name="gene1 gene2"/>
            </pathway>
        """)

        monkeypatch.setattr(ka, "extract_genes_from_kgml_string", lambda kgml_string: {"gene1", "gene2"})

        monkeypatch.setattr(ka, "extract_path_name_from_kgml", lambda kgml_string: "Fake Pathway Name")

        # Mock para save_GeneSet_GMT para só validar que foi chamada e com o conteúdo correto
        called = {}

        def fake_save(file, gene_set):
            called['file'] = file
            called['gene_set'] = gene_set

        monkeypatch.setattr(ka, "save_GeneSet_GMT", fake_save)

        # Criar arquivo temporário para output
        output_file = tmp_path / "output.gmt"

        # Chamar o método que será testado
        ka._create_GMT_file(str(output_file), org_code="eco")

        # Assertivas para garantir que o método funcionou como esperado
        assert called['file'] == str(output_file)
        assert ("path1", "Fake Pathway Name") in called['gene_set']
        assert ("path2", "Fake Pathway Name") in called['gene_set']
        assert called['gene_set'][("path1", "Fake Pathway Name")] == {"gene1", "gene2"}

    def test_save_GeneSet_GMT(self, tmp_path):
        ka = KeggAnalysis.__new__(KeggAnalysis)
        gene_set = {
            ("pathway1", "desc1"): ["geneA", "geneB"],
            ("pathway2", "desc2"): ["geneC"]
        }
        file_path = tmp_path / "test.gmt"

        ka.save_GeneSet_GMT(str(file_path), gene_set)

        assert file_path.exists()

        with open(file_path, "r") as f:
            lines = f.read().strip().split("\n")
        expected_lines = {
            "pathway1\tdesc1\tgeneA\tgeneB",
            "pathway2\tdesc2\tgeneC"
        }
        assert set(lines) == expected_lines

        with pytest.raises(FileExistsError):
            ka.save_GeneSet_GMT(str(file_path), gene_set)

    def test_enrichment_analysis(self,monkeypatch, tmp_path):

        ka = KeggAnalysis.__new__(KeggAnalysis)
        ka.file_name_gmt = "test_pathways"
        ka.get_pathway_name = lambda term: f"{term} - Description"

        mock_df = pd.DataFrame({
            'Term': ['path1', 'path2'],
            'Adjusted P-value': [0.01, 0.04],
            'Genes': ['gene1;gene2', 'gene3;gene4']
        })

        mock_enrich = MagicMock()
        mock_enrich.res2d = mock_df

        monkeypatch.setattr("src.ResPathExplorer.KeggAnalysis.gp.enrich", lambda **kwargs: mock_enrich)

        monkeypatch.setattr("src.ResPathExplorer.KeggAnalysis.rename_file", lambda a, b, c: None)

        ka.enrichment_analysis(
            gene_list=["gene1", "gene2", "gene3"],
            cutoff=0.05,
            name_outdir=str(tmp_path),
            number_path=2,
            name_results_file="results",
            genes_background=["gene1", "gene2", "gene3", "gene4"]
        )

        assert not ka.enrichment_results.empty
        assert "Pathway name" in ka.enrichment_results.columns
        assert len(ka.limited_enrichment_results) == 2
        assert "path1" in ka.paths_genes_dict
        assert ka.paths_genes_dict["path1"] == ['gene1', 'gene2']

    def test_get_pathway_name(self, monkeypatch):
        ka = KeggAnalysis.__new__(KeggAnalysis)

        response_text = (
            "ENTRY eco00010 Pathway\nNAME Glycolysis / Gluconeogenesis - Escherichia coli K-12 MG1655\nDESCRIPTION Glycolysis is the process...")

        mock_read = MagicMock(return_value=response_text)
        monkeypatch.setattr("Bio.KEGG.REST.kegg_get", lambda id_pathway: MagicMock(read=mock_read))

        pathway_name = ka.get_pathway_name("eco00010")

        assert pathway_name == "Glycolysis / Gluconeogenesis - Escherichia coli K-12 MG1655"

    def test_visualize_enrichment_results_barplot(self, tmp_path):
        ka = KeggAnalysis.__new__(KeggAnalysis)

        ka.limited_enrichment_results = pd.DataFrame({
            "Pathway name": ["Pathway A", "Pathway B"],
            "Adjusted P-value": [0.01, 0.05],
            "Overlap": [10, 5]
        })

        output_dir = tmp_path
        file_name = "test_plot"

        ka.visualize_enrichment_results(
            name_outdir=str(output_dir),
            outplot_file_name=file_name,
            plot_title="Test Plot",
            plot_type="barplot"
        )

        expected_file = output_dir / f"{file_name}_barplot.png"
        assert expected_file.exists(), "Plot file not created correctly."

    def test_visualize_enrichment_results_raises_on_none(self):
        ka = KeggAnalysis.__new__(KeggAnalysis)
        ka.limited_enrichment_results = None

        with pytest.raises(RuntimeError, match="Run enrichment_analysis\\(\\) first\\."):
            ka.visualize_enrichment_results("some_dir", "plot", "title", "barplot")

    def test_visualize_enrichment_results_invalid_plot_type(self):
        ka = KeggAnalysis.__new__(KeggAnalysis)
        ka.limited_enrichment_results = pd.DataFrame({
            "Pathway name": ["Pathway A"],
            "Adjusted P-value": [0.01],
            "Overlap": [5]
        })

        with pytest.raises(ValueError, match="Invalid plot_type"):
            ka.visualize_enrichment_results("some_dir", "plot", "title", "invalid")

    def test_search_gene_path_in_gene_set(self):
        ka = KeggAnalysis.__new__(KeggAnalysis)
        ka.gene_set = {
            ("path1",): ["geneA", "geneB"],
            ("path2",): ["geneC"],
            ("path3",): ["geneA", "geneD"]
        }

        result = ka.search_gene_path("geneA", search_in_gene_set=True)
        assert sorted(result) == ["path1", "path3"]

    def test_search_gene_path_in_dataframe(self):
        ka = KeggAnalysis.__new__(KeggAnalysis)
        df = pd.DataFrame({
            "Term": ["path1", "path2"],
            "Genes": ["geneA;geneB", "geneC;geneD"]
        })

        result = ka.search_gene_path("geneC", data=df)
        assert result == ["path2"]

    def test_search_gene_path_invalid_usage(self):
        ka = KeggAnalysis.__new__(KeggAnalysis)

        with pytest.raises(ValueError, match="Either 'data' must be provided or 'search_in_gene_set' must be True."):
            ka.search_gene_path("geneX")







