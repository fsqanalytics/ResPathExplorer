import os
import re
import requests
from .rename_file import rename_file
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
        enrichment_results (pd.DataFrame): Full enrichment results.
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

        if use_existing_gmt:
            if not os.path.exists(file_name_gmt):
                raise FileNotFoundError(f"GMT file '{file_name_gmt}' not found.")
            self._load_gmt_file(file_name_gmt)
        else:
            self._create_GMT_file(file_name_gmt)

        self.enrichment_results = None

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
            for line in file:
                parts = line.strip().split("\t")
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
        url = f"http://rest.kegg.jp/list/pathway/{org}"
        response = requests.get(url)
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
                    line = f"{name}\t{description}\t" + "\t".join(genes) + "\n"
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
        """Perform pathway enrichment analysis and store results."""
        enr = gp.enrich(
            gene_list=gene_list,
            background=genes_background,
            gene_sets=self.file_name_gmt,
            outdir=name_outdir,
            cutoff=cutoff
        )
        results = enr.res2d
        filtered_results = results[results['Adjusted P-value'] <= cutoff].copy()

        dict_names = {}
        for t in filtered_results["Term"]:
            dict_names[t] = self.get_pathway_name(t).split(" -")[0]
        filtered_res = filtered_results.copy()
        filtered_res['Pathway name'] = filtered_res['Term'].map(dict_names)

        self.enrichment_results = filtered_res

        data_frame = self.enrichment_results.sort_values(by='Adjusted P-value', ascending=True, inplace=False)
        data_final = data_frame.head(number_path)

        self.limited_enrichment_results = data_final

        term_genes_dict = data_final.set_index('Term')['Genes'].to_dict()
        for k, vs in term_genes_dict.items():
            term_genes_dict[k] = [v for v in vs.split(";")]
        self.paths_genes_dict = term_genes_dict

        rename_file(name_outdir, f"{self.file_name_gmt}.human.enrichr.reports.txt", f"{name_results_file}.txt")
        rename_file(name_outdir, f"{self.file_name_gmt}.human.enrichr.reports.pdf", f"{name_results_file}.pdf")

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
        name = re.search(r'NAME\s+(.+)', res).group(1)
        return name

    def visualize_enrichment_results(
            self,
            name_outdir: str,
            outplot_file_name: str,
            plot_title: str,
            plot_type: str = "barplot"
    ) -> None:
        """Visualize enrichment results using a barplot or dotplot."""
        if self.limited_enrichment_results is None:
            raise RuntimeError("Run enrichment_analysis() first.")

        data = self.limited_enrichment_results.copy()
        fig, ax = plt.subplots(figsize=(12, 8))

        if plot_type == "barplot":
            data["-log10(AdjP)"] = -np.log10(data["Adjusted P-value"])
            norm = plt.Normalize(data["-log10(AdjP)"].min(), data["-log10(AdjP)"].max())
            cmap = plt.cm.viridis

            bars = ax.barh(data["Pathway name"], data["-log10(AdjP)"], color=cmap(norm(data["-log10(AdjP)"])))
            ax.set_xlabel('-log10(Adjusted P-value)')
            ax.set_title(f'{plot_title} - Barplot')
            ax.invert_yaxis()

            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            plt.colorbar(sm, ax=ax, label='-log10(Adjusted P-value)')

        elif plot_type == "dotplot":
            data["-log10(AdjP)"] = -np.log10(data["Adjusted P-value"])
            sns.scatterplot(
                x="-log10(AdjP)",
                y="Pathway name",
                size="Overlap",
                hue="-log10(AdjP)",
                data=data,
                sizes=(50, 300),
                palette="viridis",
                legend=None,
                ax=ax
            )
            ax.set_title(f'{plot_title} - Dotplot')
            ax.set_xlabel('-log10(Adjusted P-value)')
            ax.invert_yaxis()

        else:
            raise ValueError("Invalid plot_type. Use 'barplot' or 'dotplot'.")

        plt.tight_layout()
        plot_path = os.path.join(name_outdir, f"{outplot_file_name}_{plot_type}.png")
        plt.savefig(plot_path, bbox_inches='tight')
        plt.show()

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
            return data[data["Genes"].str.contains(gene_id)]["Term"].tolist()
        else:
            raise ValueError("Either 'data' must be provided or 'search_in_gene_set' must be True.")


