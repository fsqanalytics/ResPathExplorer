# 🧬 ResPathExplorer

A modular Python library for functional analysis of microbial genes using curated bioinformatics databases. This tool supports gene mapping to KEGG pathways, resistance profiles via CARD, and virulence factors through VFDB — providing a comprehensive view of microbial functionality in contexts such as immunology, food safety, and antimicrobial resistance.

## Key Features

- Functional Mapping: Annotate genes using KEGG, CARD, and VFDB.

- Enrichment Analysis: Perform pathway enrichment from gene lists.

- Pathway-Level Analysis: Visualize gene distribution across biological pathways.

- Resistance & Virulence Profiling: Contextualize resistance and virulence genes functionally.

- Custom Data Input: Supports gene lists and annotation tables in standard formats.

- Extensible Pipeline: Modular design for easy integration into broader omics workflows.

- Visualization Tools: Built-in plots for resistance profiles, virulence categories, enrichment results, and KEGG pathways.

## Examples
You can find full working examples in the `Examples/` folder:

#### `Foodborne bacteria/`:

Transcriptome Analysis of *Listeria monocytogenes* Exposed to Beef Fat Reveals Antimicrobial and Pathogenicity Attenuation Mechanisms

doi: https://doi.org/10.1128/AEM.03027-20

- The example was done just for C18:2n-6

`prepare_data.ipynb`: Preparing data for analysis.

`Analysis.ipynb`: ResPathExplorer application. Run the preparation notebook first;
see [example setup and inputs](Examples/Foodborne%20bacteria/README.md).

## Structure
```text
📁 ResPathExplorer/
├── 📁 src/ResPathExplorer/
│   ├── __init__.py
│   ├── CARDAnalysis.py
│   ├── KeggAnalysis.py
│   ├── URL_pathway.py
│   ├── VFDBAnalysis.py
│   ├── mapper_KeggFunctions.py
│   ├── rename_file.py
│   ├── save_df_as_html.py
│   └── validate_color_code.py
├── 📁 Examples/
│   ├── 📁 Foodborne bacteria/
├── 📁 tests/
```

## 🛠 Installation

To install the library directly from GitHub:

```bash
pip install git+https://github.com/lais-carvalho/ResPathExplorer.git
```

For local development and notebooks:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev,notebooks]"
python -m pytest
cd "Examples/Foodborne bacteria"
jupyter lab
```

The Excel readers `xlrd` (`.xls`) and `openpyxl` (`.xlsx`) are included in the
library dependencies. Tests use local fixtures and mocked database responses;
they do not require access to KEGG, CARD, or VFDB. Both `test_*.py` and `Test*.py`
test modules are discovered automatically.

## Analysis outputs and interpretation

- `KeggAnalysis.enrichment_analysis()` writes `<name_results_file>.txt` (all
  pathways returned by GSEApy, tab-separated) and `<name_results_file>.pdf` (the top
  significant pathways). Reports are written even when no pathways pass the
  adjusted P-value cutoff. `all_enrichment_results` contains the full table;
  `enrichment_results` contains significant rows and `limited_enrichment_results`
  contains the top N. GMT descriptions supply pathway names locally.
- Supply experiment-specific `genes_background` and use the same identifier
  format in the query, background, and GMT. This library performs pathway
  over-representation analysis, rather than ranked-list GSEA.
- CARD and VFDB annotations match gene names and synonyms. They do not align
  sequences, identify resistance-conferring variants, or establish phenotypes.
- VFDB searches retain distinct strain records. Category plots count each
  gene-category association once across strains. Unknown categories remain
  visible; missing categories appear as `Unclassified`. If no gene-name column
  is supplied, plots count input rows instead. Percentages sum to 100% and a gene
  assigned to two categories contributes one association to each.
- `search_gene_id_kegg()` returns the first match across organisms when
  `org_code` is omitted; specify the organism for biological interpretation.
  `save_df_as_html()` propagates write errors so callers can handle failed exports.

## Acknowledgements

- European Food Safety Authority (EFSA) – support via the “Pathogens-in-Foods Database” project.

- Centro de Investigação da Montanha (CIMO), Portugal.

# Authors

Developed by Laís Carvalho, Ursula Gonzales-Barron, and Vasco Cadavez at CIMO, Instituto Politécnico de Bragança, Portugal.
