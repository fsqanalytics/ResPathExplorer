"""Execute the example cells with synthetic study inputs and mocked databases."""

import gzip
import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd


EXAMPLES = Path(__file__).resolve().parents[2] / "Examples" / "Foodborne bacteria"


def run_notebook_cells(name):
    notebook = json.loads((EXAMPLES / name).read_text(encoding="utf-8"))
    namespace = {}
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] == "code":
            code = "".join(cell["source"])
            exec(compile(code, f"{name}:cell{index}", "exec"), namespace)
    return namespace


def prepare_study_inputs(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(EXAMPLES.parents[1] / "src"))
    rows = [{
        "Gene_id": f"lmo{i:04d}",
        "gene_name": f"gene{i}",
        "C18:2n-6_fold_change": 3 if i == 1 else 0.3 if i == 2 else 1,
        "C18:2n-6_pValue": 0.01 if i <= 2 else 1,
        "C18:2n-6_fdr": 0.01 if i <= 2 else 1
    } for i in range(1, 101)]
    rows.append({**rows[-1], "Gene_id": "unknown", "gene_name": "unknown"})
    pd.DataFrame(rows).to_excel(tmp_path / "aem.03027-20-s0002.xlsx", index=False)

    with patch("Bio.KEGG.REST.kegg_find") as find, patch("Bio.KEGG.REST.kegg_get") as get:
        find.return_value.read.return_value = ""
        get.side_effect = lambda identifier: type("Response", (), {
            "read": lambda self: "SYMBOL geneA, aliasA" if identifier.endswith("0001") else "SYMBOL geneB"
        })()
        return run_notebook_cells("prepare_data.ipynb")


def test_preparation_notebook_generates_analysis_inputs_from_source(tmp_path, monkeypatch):
    result = prepare_study_inputs(tmp_path, monkeypatch)
    assert len(result["background"]) == 100
    assert (tmp_path / "upregulated.txt").read_text().splitlines() == ["lmo0001"]
    assert (tmp_path / "downregulted.txt").read_text().splitlines() == ["lmo0002"]
    assert (tmp_path / "upregulatedGNames.txt").read_text().splitlines() == ["geneA"]
    assert len(pd.read_excel(tmp_path / "unmapped_genes.xlsx")) == 1
    assert (tmp_path / "df_all.xlsx").exists()


def test_analysis_notebook_uses_independent_colors_for_expression_groups(tmp_path, monkeypatch):
    prepare_study_inputs(tmp_path, monkeypatch)
    (tmp_path / "exemple_fb.gmt").write_text(
        "lmo00010\tUp pathway\tlmo0001\n"
        "lmo00020\tDown pathway\tlmo0002\n"
    )
    (tmp_path / "aro.obo").write_text(
        "[Term]\nid: ARO:1\nname: geneA\n"
        "relationship: confers_resistance_to_antibiotic ARO:100 ! Penicillin\n"
        "[Term]\nid: ARO:2\nname: geneB\n"
    )
    database = tmp_path / "vfdb_data"
    database.mkdir()
    with gzip.open(database / "VFDB_setA_pro.fas.gz", "wt") as fasta:
        fasta.write(
            ">VFG1(gb|WP_1) (geneA) description [Factor (VF0001) - Adherence (VFC1)] "
            "[Listeria monocytogenes strain A]\nMNR\n"
        )
    (database / "VFs.xls").write_bytes(b"fixture replaced by read_excel mock")
    read_excel = pd.read_excel

    def read_input(path, *args, **kwargs):
        if Path(path).suffix == ".xls":
            return pd.DataFrame({"VFID": ["VF0001"], "VF_Name": ["Factor"], "Function": ["Adherence"]})
        return read_excel(path, *args, **kwargs)

    with patch("Bio.KEGG.REST.kegg_list") as organisms, \
            patch("pandas.read_excel", side_effect=read_input), \
            patch("ResPathExplorer.URL_pathway.get_url_pathway", return_value="https://example.invalid/pathway") as url:
        organisms.return_value.read.return_value = (
            "T1\tlmo\tListeria monocytogenes EGD-e (serotype 1/2a)\n"
        )
        result = run_notebook_cells("Analysis.ipynb")

    assert result["gene_color_dict_up"] == {"lmo0001": "skyblue,purple"}
    assert result["gene_color_dict_down"] == {"lmo0002": "indianred,blue"}
    assert result["profile_dict"]["downregulated_genes"]["resistance_genes"] == ["geneB"]
    url.assert_any_call("lmo00010", {"lmo0001": "skyblue,purple"})
    url.assert_any_call("lmo00020", {"lmo0002": "indianred,blue"})
