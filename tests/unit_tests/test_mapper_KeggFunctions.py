import pytest
from unittest.mock import patch, mock_open, MagicMock
from src.ResPathExplorer.mapper_KeggFunctions import search_gene_id_kegg, get_gene_name_by_kegg_id


class Testmapper_KeggFunctions:

    def test_search_gene_id_kegg(self, monkeypatch):
        org = "eco"

        def fake_kegg_find(db, name):
            class FakeResponse:
                def read(self_inner):
                    if name == "gene_name":
                        return "eco:b0002\tgene_name\neco:b0003\tgene_name2\n"
                    else:
                        return ""

            return FakeResponse()

        monkeypatch.setattr("Bio.KEGG.REST.kegg_find", fake_kegg_find)

        gene_id = search_gene_id_kegg("gene_name", org_code= org)
        assert gene_id == "eco:b0002"

        gene_id_none = search_gene_id_kegg("nonexistent")
        assert gene_id_none is None


    def test_get_gene_name_by_kegg_id(self, monkeypatch):

        response_text = """ENTRY ecob0002 Pathway\nNAME thrL\nDESCRIPTION thr leader peptide\nSYMBOL thrL"""

        mock_read = MagicMock(return_value=response_text)
        monkeypatch.setattr("Bio.KEGG.REST.kegg_get", lambda kegg_id: MagicMock(read=mock_read))

        symbol = get_gene_name_by_kegg_id("eco:b0002")
        assert symbol == "thrL"

        with pytest.raises(ValueError):
            get_gene_name_by_kegg_id("invalid_id")

        monkeypatch.setattr("Bio.KEGG.REST.kegg_get",
                            lambda kegg_id: MagicMock(read=MagicMock(return_value="NO SYMBOL HERE")))
        symbol_none = get_gene_name_by_kegg_id("eco:b0002")
        assert symbol_none is None