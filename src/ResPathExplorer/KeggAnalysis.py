import os
import re
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import xml.etree.ElementTree as ET
from typing import List, Dict, Tuple, Optional
from bioservices import KEGG
from Bio.KEGG import REST
import gseapy as gp


class KeggAnalysis:
    """
    A class for managing KEGG pathway analysis, including GMT file generation,
    enrichment analysis, and visualization.

    Attributes:
        organism (str): Full name of the organism.
        org (str): KEGG organism code.
        gene_set (Dict[Tuple[str, str], List[str]]): Dictionary of pathways and their associated genes.
        file_name_gmt (str): Path to GMT file.
        all_enrichment_results (pd.DataFrame): Full enrichment results.
        enrichment_results (pd.DataFrame): Results passing the adjusted P-value cutoff.
        limited_enrichment_results (pd.DataFrame): Top N filtered enrichment results.
        paths_genes_dict (Dict[str, List[str]]): Dictionary of pathways and enriched genes.
    """

    def __init__(self, organism_name: str, file_name_gmt: str, use_existing_gmt: bool = True):
        if not organism_name:
            raise ValueError("Organism name cannot be empty.")
        if not file_name_gmt:
            raise ValueError("A GMT file name must be provided.")

        # Determine KEGG organism code and full name
        if " " in organism_name or len(organism_name) > 4:
            self.org = self._get_organism_prefix(organism_name)
            self.organism = organism_name
        else:
            self.org = organism_name
            self.organism = self._get_organism_name(organism_name)

        self.file_name_gmt = file_name_gmt
        self.gene_set = {}
        self.all_enrichment_results = None
        self.enrichment_results = None
        self.limited_enrichment_results = None
        self.paths_genes_dict = {}

        if use_existing_gmt:
            if not os.path.exists(file_name_gmt):
                raise FileNotFoundError(f"GMT file '{file_name_gmt}' not found.")
            self._load_gmt_file(file_name_gmt)
        else:
            self._create_GMT_file(file_name_gmt)

    def _get_organism_prefix(self, organism_name: str) -> str:
        """Find the KEGG code for a given full organism name."""
        organisms_list = REST.kegg_list("organism").read()
        organism_prefix = None
        for line in organisms_list.strip().split("\n"):
            entry = line.split("\t")
            org_code, org_id, org_name = entry[1], entry[0], entry[2]
            if organism_name.lower() in org_name.lower():
                organism_prefix = org_code
                break
        if organism_prefix:
            return organism_prefix
        else:
            raise ValueError(f"Organism name '{organism_name}' not found in KEGG database.")

    def _get_organism_name(self, organism_code: str) -> str:
        """Find the full organism name from a KEGG code."""
        organisms_list = REST.kegg_list("organism").read()
        organism_name = None
        for line in organisms_list.strip().split("\n"):
            entry = line.split("\t")
            org_code, org_id, org_name = entry[1], entry[0], entry[2]
            if organism_code.lower() == org_code.lower():
                organism_name = org_name
                break
        if organism_name:
            return organism_name
        else:
            raise ValueError(f"Organism code '{organism_code}' not found in KEGG database.")

    def _load_gmt_file(self, file_name: str) -> None:
        """Load pathways and genes from a GMT file into `self.gene_set`."""
        with open(file_name, 'r') as file:
            for line_number, line in enumerate(file, 1):
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 3 or not parts[0]:
                    raise ValueError(f"Invalid GMT record on line {line_number} of '{file_name}'.")
                pathway = parts[0]
                pathway_name = parts[1]
                genes = parts[2:]
                self.gene_set[(pathway, pathway_name)] = genes

    def get_kgml(self, pathway_id, service):
        """
        Get the kgml file of the desired pathway from the KEGG database
        """
        kgml_str = service.get(pathway_id, "kgml")
        if kgml_str is None:
            raise ValueError(f"Failed to get KGML for {pathway_id}")
        return kgml_str

    def extract_path_name_from_kgml(self, kgml_string):
        """
        Extract the pathway name from a kgml structure string
        """
        root = ET.fromstring(kgml_string)
        return root.attrib.get("title")

    def extract_genes_from_kgml_string(self, kgml_string):
        """
        Extracts genes from a KGML string, as returned by the KEGG API.
        """
        root = ET.fromstring(kgml_string)

        genes = set()
        for entry in root.findall('entry'):
            if entry.attrib.get('type') == 'gene':
                raw_names = entry.attrib.get('name', '')
                gene_ids = raw_names.strip().split()
                for gene in gene_ids:
                    genes.add(gene.split(":")[-1])
        return genes

    def _create_GMT_file(self, output_file: str, org_code: Optional[str] = None, service=None) -> None:
        """Fetch KEGG pathways and their genes, and save to a GMT file."""
        if service is None:
            service = KEGG()
        org = org_code or self.org
        url = f"https://rest.kegg.jp/list/pathway/{org}"
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        resp = response.text

        paths_loaded = []

        for line in resp.split("\n"):
            path = line.split("\t")[0]
            if len(path) > 0:
                try:
                    kgml_string = self.get_kgml(path, service)
                    genes = self.extract_genes_from_kgml_string(kgml_string)
                    path_name = self.extract_path_name_from_kgml(kgml_string)

                    if len(genes) != 0:
                        self.gene_set[(path, path_name)] = genes
                        paths_loaded.append(path)
                    else:
                        print(f"For pathway {path} no genes were found")

                except Exception as e:
                    print(f"Pathway {path} has no information in KEGG, it was not added to the gene set: {e}")
                    continue

        if not self.gene_set:
            raise ValueError(f"No pathways with genes could be loaded for '{org}'.")
        self.save_GeneSet_GMT(output_file, self.gene_set)
        print(f"GMT {output_file} saved with {len(paths_loaded)} pathways")

    def save_GeneSet_GMT(self, gmt_file_name: str, gene_set_dict: Dict[Tuple[str, str], List[str]]) -> None:
        """
        Save a gene set dictionary to a .gmt file if it doesn't already exist.

        Args:
            gmt_file_name (str): Output GMT file name.
            gene_set_dict (dict): Keys as (pathway_name, description), values as gene lists.

        Raises:
            FileExistsError: If the file already exists.
            Exception: For other I/O errors.
        """
        if os.path.exists(gmt_file_name):
            raise FileExistsError(f"File '{gmt_file_name}' already exists.")

        try:
            with open(gmt_file_name, "w") as f:
                for (name, description), genes in gene_set_dict.items():
                    line = f"{name}\t{description}\t" + "\t".join(sorted(genes)) + "\n"
                    f.write(line)
            print(f"GMT file saved: {gmt_file_name}")
        except Exception as e:
            raise Exception(f"Error saving file '{gmt_file_name}': {e}")

    def enrichment_analysis(
            self,
            gene_list: List[str],
            cutoff: float,
            name_outdir: str,
            number_path: int,
            name_results_file: str,
            genes_background: Optional[List[str]] = None) -> None:
        """Store enrichment results and write named TSV and PDF reports.

        The TSV contains all pathways returned by GSEApy. ``enrichment_results`` contains
        significant pathways; ``limited_enrichment_results`` contains the top N.
        Pathway names come from the GMT descriptions, without extra API calls.
        """
        if not isinstance(gene_list, list) or not gene_list or not all(
                isinstance(gene, str) and gene.strip() for gene in gene_list):
            raise ValueError("gene_list must be a non-empty list of gene identifiers.")
        if not 0 < cutoff <= 1:
            raise ValueError("cutoff must be greater than 0 and at most 1.")
        if not isinstance(number_path, int) or isinstance(number_path, bool) or number_path < 1:
            raise ValueError("number_path must be a positive integer.")
        os.makedirs(name_outdir, exist_ok=True)
        enr = gp.enrich(
            gene_list=gene_list,
            background=genes_background,
            gene_sets=self.file_name_gmt,
            outdir=None,
            cutoff=cutoff,
            no_plot=True
        )
        results = enr.res2d
        if results is None or (isinstance(results, list) and not results):
            results = pd.DataFrame(columns=[
                "Term", "Overlap", "P-value", "Adjusted P-value", "Genes"
            ])
        results = results.copy()
        names = {key[0]: key[1] for key in getattr(self, "gene_set", {})}
        results["Pathway name"] = results["Term"].map(names).fillna(results["Term"])
        self.all_enrichment_results = results
        filtered_results = results[results['Adjusted P-value'] <= cutoff].copy()
        self.enrichment_results = filtered_results

        data_frame = self.enrichment_results.sort_values(by='Adjusted P-value', ascending=True, inplace=False)
        data_final = data_frame.head(number_path)

        self.limited_enrichment_results = data_final

        term_genes_dict = data_final.set_index('Term')['Genes'].to_dict()
        for k, vs in term_genes_dict.items():
            term_genes_dict[k] = [v.strip() for v in vs.split(";") if v.strip()]
        self.paths_genes_dict = term_genes_dict

        results.to_csv(os.path.join(name_outdir, f"{name_results_file}.txt"), sep="\t", index=False)
        fig = self._make_enrichment_figure("Pathway enrichment", "barplot")
        try:
            fig.savefig(os.path.join(name_outdir, f"{name_results_file}.pdf"), bbox_inches="tight")
        finally:
            plt.close(fig)

    def get_pathway_name(self, id_pathway: str) -> str:
        """Fetch the pathway name given a KEGG pathway ID."""
        dic = {}
        result = REST.kegg_get(id_pathway).read()
        lin = result.split("\n")
        res = ""
        for l in lin:
            if l.startswith("NAME"):
                res = l
        import re
        match = re.search(r'NAME\s+(.+)', res)
        if match is None:
            raise ValueError(f"No pathway name found for '{id_pathway}'.")
        return match.group(1)

    def visualize_enrichment_results(
            self,
            name_outdir: str,
            outplot_file_name: str,
            plot_title: str,
            plot_type: str = "barplot"
    ) -> None:
        """Visualize enrichment results using a barplot or dotplot."""
        fig = self._make_enrichment_figure(plot_title, plot_type)
        try:
            os.makedirs(name_outdir, exist_ok=True)
            plot_path = os.path.join(name_outdir, f"{outplot_file_name}_{plot_type}.png")
            fig.savefig(plot_path, bbox_inches='tight')
            plt.show()
        finally:
            plt.close(fig)

    def _make_enrichment_figure(self, plot_title: str, plot_type: str):
        """Build a figure shared by PDF reports and interactive visualizations."""
        if getattr(self, "limited_enrichment_results", None) is None:
            raise RuntimeError("Run enrichment_analysis() first.")
        if plot_type not in {"barplot", "dotplot"}:
            raise ValueError("Invalid plot_type. Use 'barplot' or 'dotplot'.")

        data = self.limited_enrichment_results.copy()
        fig, ax = plt.subplots(figsize=(12, 8))
        if data.empty:
            ax.text(0.5, 0.5, "No pathways pass the adjusted P-value cutoff.",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(plot_title)
            ax.set_axis_off()
            return fig

        data["-log10(AdjP)"] = -np.log10(
            data["Adjusted P-value"].astype(float).clip(lower=np.finfo(float).tiny)
        )

        if plot_type == "barplot":
            norm = plt.Normalize(data["-log10(AdjP)"].min(), data["-log10(AdjP)"].max())
            cmap = plt.cm.viridis

            bars = ax.barh(data["Pathway name"], data["-log10(AdjP)"], color=cmap(norm(data["-log10(AdjP)"])))
            ax.set_xlabel('-log10(Adjusted P-value)')
            ax.set_title(f'{plot_title} - Barplot')
            ax.yaxis.set_inverted(True)

            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            plt.colorbar(sm, ax=ax, label='-log10(Adjusted P-value)')

        elif plot_type == "dotplot":
            data["Matched genes"] = pd.to_numeric(
                data["Overlap"].astype(str).str.split("/", n=1).str[0]
            )
            sns.scatterplot(
                x="-log10(AdjP)",
                y="Pathway name",
                size="Matched genes",
                hue="-log10(AdjP)",
                data=data,
                sizes=(50, 300),
                palette="viridis",
                legend="brief",
                ax=ax
            )
            ax.set_title(f'{plot_title} - Dotplot')
            ax.set_xlabel('-log10(Adjusted P-value)')
            ax.yaxis.set_inverted(True)

        fig.tight_layout()
        return fig

    def search_gene_path(
            self,
            gene_id: str,
            data: Optional[pd.DataFrame] = None,
            search_in_gene_set: bool = False
    ) -> List[str]:
        """
        Find KEGG pathways containing a specific gene, either from enrichment data or the full gene set.
        """
        if search_in_gene_set:
            return [k[0] for k, genes in self.gene_set.items() if gene_id in genes]
        elif data is not None:
            mask = data["Genes"].map(
                lambda value: isinstance(value, str)
                and gene_id in {gene.strip() for gene in value.split(";")}
            )
            return data.loc[mask, "Term"].tolist()
        else:
            raise ValueError("Either 'data' must be provided or 'search_in_gene_set' must be True.")


