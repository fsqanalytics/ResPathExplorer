"""Regression checks for annotation integrity and complete local analysis flows."""

from io import BytesIO
import tarfile
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import seaborn as sns

from src.ResPathExplorer.CARDAnalysis import CARDAnalysis
from src.ResPathExplorer.KeggAnalysis import KeggAnalysis
from src.ResPathExplorer.VFDBAnalysis import VFDBAnalysis
from src.ResPathExplorer.mapper_KeggFunctions import search_gene_id_kegg
from src.ResPathExplorer.save_df_as_html import save_df_as_html


@pytest.mark.parametrize("absolute_gmt", [False, True])
def test_real_local_enrichment_writes_named_reports_without_network(tmp_path, absolute_gmt):
    gmt = tmp_path / "custom pathways.gmt"
    gmt.write_text("eco00010\tSynthetic pathway\tg1\tg2\tg3\tg4\n")
    with patch.object(KeggAnalysis, "_get_organism_name", return_value="Escherichia coli"):
        ka = KeggAnalysis("eco", str(gmt) if absolute_gmt else gmt.name)

    outdir = tmp_path / "nested" / "reports"
    ka.enrichment_analysis(
        ["g1", "g2", "g3", "g4"], 0.05, str(outdir), 10, "analysis",
        [f"g{i}" for i in range(1, 101)]
    )

    table = pd.read_csv(outdir / "analysis.txt", sep="\t")
    assert table["Term"].tolist() == ["eco00010"]
    assert table["Pathway name"].tolist() == ["Synthetic pathway"]
    assert table["Adjusted P-value"].iloc[0] < 0.05
    assert set(ka.paths_genes_dict["eco00010"]) == {"g1", "g2", "g3", "g4"}
    assert (outdir / "analysis.pdf").read_bytes().startswith(b"%PDF")


@pytest.mark.parametrize("result_kind", ["not_significant", "no_overlap"])
def test_enrichment_without_significant_pathways_still_exports_reports(tmp_path, result_kind):
    gmt = tmp_path / "background.gmt"
    gmt.write_text("eco00010\tCommon pathway\tg1\tg2\tg3\tg4\n")
    with patch.object(KeggAnalysis, "_get_organism_name", return_value="Escherichia coli"):
        ka = KeggAnalysis("eco", str(gmt))
    query = ["g1"] if result_kind == "not_significant" else ["g5"]
    outdir = tmp_path / "reports"
    ka.enrichment_analysis(query, 0.05, str(outdir), 10, "empty", [f"g{i}" for i in range(1, 6)])

    assert ka.enrichment_results.empty
    assert ka.limited_enrichment_results.empty
    assert ka.paths_genes_dict == {}
    exported = pd.read_csv(outdir / "empty.txt", sep="\t")
    assert len(exported) == (1 if result_kind == "not_significant" else 0)
    assert (outdir / "empty.pdf").read_bytes().startswith(b"%PDF")
    ka.visualize_enrichment_results(str(outdir), "empty", "No enriched pathways")
    assert (outdir / "empty_barplot.png").exists()


def test_plot_before_enrichment_reports_clear_error(tmp_path):
    gmt = tmp_path / "minimal.gmt"
    gmt.write_text("eco00010\tPathway\tg1\n")
    with patch.object(KeggAnalysis, "_get_organism_name", return_value="Escherichia coli"):
        ka = KeggAnalysis("eco", str(gmt))
    with pytest.raises(RuntimeError, match="Run enrichment_analysis"):
        ka.visualize_enrichment_results(str(tmp_path), "plot", "title")


def test_dotplot_sizes_represent_matched_gene_counts_and_zero_pvalues():
    ka = KeggAnalysis.__new__(KeggAnalysis)
    ka.limited_enrichment_results = pd.DataFrame({
        "Pathway name": ["A", "B", "C"],
        "Adjusted P-value": [0, 0.01, 0.02],
        "Overlap": ["2/100", "10/100", "2/20"]
    })
    fig = ka._make_enrichment_figure("Enrichment", "dotplot")
    points = fig.axes[0].collections[0]
    sizes = points.get_sizes()
    assert sizes[0] == sizes[2]
    assert sizes[1] > sizes[0]
    assert np.isfinite(points.get_offsets()).all()
    assert fig.axes[0].yaxis_inverted()


def test_search_pathway_uses_complete_literal_gene_identifiers():
    ka = KeggAnalysis.__new__(KeggAnalysis)
    data = pd.DataFrame({
        "Term": ["exact", "substring", "literal", "regex", "missing"],
        "Genes": ["geneA; geneB", "geneAB;geneC", "gene.1", "geneX1", None]
    })
    assert ka.search_gene_path("geneA", data=data) == ["exact"]
    assert ka.search_gene_path("gene.1", data=data) == ["literal"]


@pytest.mark.parametrize("following_stanza", ["", "\n[Typedef]\nname: other\n"])
def test_card_processes_final_term_and_stanza_boundaries(tmp_path, following_stanza):
    obo = tmp_path / "aro.obo"
    obo.write_text(
        "[Term]\nid: ARO:1\nname: first\n"
        "relationship: confers_resistance_to_antibiotic ARO:100 ! Penicillin\n"
        "\n[Term]\nid: ARO:2\nname: geneA\nsynonym: \"aliasA\" EXACT []\n"
        + following_stanza
    )
    card = CARDAnalysis.__new__(CARDAnalysis)
    match = card.find_gene_ids(str(obo), "aliasA")
    assert match["Gene ID"] == "ARO:2"
    assert match["Matched Name"] == "aliasA"
    assert pd.isna(match["Antibiotics"])
    assert card.find_gene_ids(str(obo), "other") is None


def test_card_empty_matches_keep_table_schema(tmp_path):
    (tmp_path / "aro.obo").write_text("[Term]\nid: ARO:1\nname: geneA\n")
    card = CARDAnalysis(["absent"], has_CARDdata=True)
    assert card.ARGdf.empty
    assert {"Matched Name", "Gene ID", "Antibiotics"}.issubset(card.ARGdf.columns)
    assert card.not_ARG_list == ["absent"]


def test_card_download_copies_nested_ontology_to_expected_path(tmp_path):
    archive = BytesIO()
    content = b"[Term]\nid: ARO:1\nname: geneA\n"
    with tarfile.open(fileobj=archive, mode="w:bz2") as tar:
        entry = tarfile.TarInfo("nested/aro.obo")
        entry.size = len(content)
        tar.addfile(entry, BytesIO(content))
    with patch("src.ResPathExplorer.CARDAnalysis.requests.get") as get:
        get.return_value.status_code = 200
        get.return_value.iter_content.return_value = [archive.getvalue()]
        CARDAnalysis.__new__(CARDAnalysis).download_CARD_file()
    assert (tmp_path / "aro.obo").read_bytes() == content
    assert not (tmp_path / "card-data.tar.bz2").exists()


def test_vfdb_search_preserves_strains_without_repeating_queries(tmp_path):
    vfdb = VFDBAnalysis(str(tmp_path / "db"))
    vfdb.df_genes = pd.DataFrame({
        "Gene_Name": ["GeneA", "geneA", "geneA", "geneB"],
        "Bacteria": ["E.coli A", "E.coli B", "EXcoli C", "Salmonella"],
        "VFID": ["VF1", "VF1", "VF1", "VF2"]
    })
    result = vfdb.search_virulence_genes(["geneA", " GeneA "], "E.coli")
    assert result["Bacteria"].tolist() == ["E.coli A", "E.coli B"]
    missing = vfdb.search_virulence_genes(["absent"], "E.coli")
    assert missing.empty
    assert list(missing.columns) == list(vfdb.df_genes.columns)


def test_vfdb_distribution_deduplicates_genes_and_includes_all_categories(tmp_path):
    vfdb = VFDBAnalysis(str(tmp_path / "db"))
    data = pd.DataFrame({
        "Gene_Name": ["geneA", "GeneA", "geneB", "geneC"],
        "Functional category": ["Adherence", "Adherence", "Novel category", None]
    })
    with patch("src.ResPathExplorer.VFDBAnalysis.sns.barplot", wraps=sns.barplot) as plot:
        vfdb.plot_virulence_factors_percentage(data, "Species", show_all_categories=False)
    frequencies = plot.call_args.kwargs["data"].set_index("Functional category")["Percentage"]
    assert set(frequencies.index) == {"Adherence", "Novel category", "Unclassified"}
    assert frequencies.sum() == pytest.approx(100)
    assert frequencies["Adherence"] == pytest.approx(100 / 3)


def test_vfdb_excel_http_error_is_reported_before_decompression(tmp_path):
    import requests

    vfdb = VFDBAnalysis(str(tmp_path / "db"))
    (tmp_path / "db" / "VFDB_setA_pro.fas.gz").write_bytes(b"cached")
    with patch("src.ResPathExplorer.VFDBAnalysis.requests.get") as get:
        response = get.return_value.__enter__.return_value
        response.raise_for_status.side_effect = requests.HTTPError("503")
        with pytest.raises(requests.HTTPError):
            vfdb.download_data()
    assert not (tmp_path / "db" / "VFs.xls").exists()


def test_mapping_without_organism_and_with_literal_organism_filter():
    with patch("src.ResPathExplorer.mapper_KeggFunctions.REST.kegg_find") as find:
        find.return_value.read.return_value = "eco:b0002\tgeneA\nmmu:123\tgeneA\n"
        assert search_gene_id_kegg("geneA") == "eco:b0002"
        assert search_gene_id_kegg("geneA", "mmu") == "mmu:123"
        assert search_gene_id_kegg("geneA", "e.*") is None


def test_html_write_errors_propagate_and_titles_are_escaped(tmp_path):
    data = pd.DataFrame({"Gene": ["geneA"]})
    with pytest.raises(OSError):
        save_df_as_html(data, str(tmp_path / "missing" / "table.html"))
    output = tmp_path / "table.html"
    save_df_as_html(data, str(output), title="<Gene & pathway>")
    assert "&lt;Gene &amp; pathway&gt;" in output.read_text()
    assert "<Gene & pathway>" not in output.read_text()
