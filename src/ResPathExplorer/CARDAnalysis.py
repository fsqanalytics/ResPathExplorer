import requests
import os
import re
import tarfile
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Tuple, Optional


class CARDAnalysis:
    """
    A class for identifying antibiotic resistance genes (ARGs) using the CARD ontology (ARO).

    Attributes:
        genes_list (List[str]): List of gene names to analyze.
        ARG_list (List[dict]): List of successfully matched ARG entries.
        not_ARG_list (List[str]): List of genes with no match in the CARD ontology.
        ARGdf (pd.DataFrame): DataFrame containing the ARG_list.
    """

    def __init__(self, genes_list: List[str], has_CARDdata: bool = False):
        """
        Initializes the CARDAnalysis class, optionally downloading CARD data if not available.

        Args:
            genes_list (List[str]): Gene names to be checked against CARD ontology.
            has_CARDdata (bool): Set to True if `aro.obo` is already available locally.

        Raises:
            ValueError: If `genes_list` is not a list of strings.
        """
        if not isinstance(genes_list, list) or not all(isinstance(g, str) for g in genes_list):
            raise ValueError("'genes_list' must be a list of strings representing gene names.")

        if not has_CARDdata:
            self.download_CARD_file()

        self.genes_list = genes_list
        self.ARG_list, self.not_ARG_list = self.finding_ARG(self.genes_list, "aro.obo")
        self.ARGdf = pd.DataFrame(self.ARG_list)

    def download_CARD_file(self) -> None:
        """
        Downloads the CARD ontology file (`aro.obo`) from the official site and extracts it.
        """
        CARD_URL = "https://card.mcmaster.ca/latest/ontology"
        zip_filename = "card-data.tar.bz2"
        obo_filename = "aro.obo"

        response = requests.get(CARD_URL, stream=True)

        if response.status_code != 200:
            raise ConnectionError(f"Failed to download CARD data (status code {response.status_code})")

        with open(zip_filename, "wb") as f:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)

        with tarfile.open(zip_filename, "r:bz2") as tar:
            for member in tar.getmembers():
                if obo_filename in member.name:
                    tar.extract(member, path=".")
                    print(f"Download and extraction complete: {member.name}")
                    break

        os.remove(zip_filename)

    def find_gene_ids(self, obo_file: str, gene_name: str) -> Optional[dict]:
        """
        Searches the `aro.obo` file for a gene and returns its annotation metadata
        """
        if not os.path.exists(obo_file):
            raise FileNotFoundError(f"File not found: {obo_file}")

        with open(obo_file, "r", encoding="utf-8") as f:
            lines = f.readlines()

        current_id = None
        current_name = None
        current_description = None
        current_antibiotics = []
        current_synonyms = []
        gene_occurrences = []

        for line in lines:
            line = line.strip()

            if line.startswith("id: ARO:"):
                current_id = line.split(": ")[1]

            elif line.startswith("name: "):
                current_name = line.split(": ", 1)[1]
                current_synonyms = []


            elif line.startswith("synonym: "):
                match = re.search(r'^synonym: "([^"]+)"', line)
                if match:
                    synonyms_split = [s.strip() for s in match.group(1).split(",")]
                    current_synonyms.extend(synonyms_split)

            elif line.startswith("def: "):
                current_description = line.split(": ", 1)[1] if ": " in line else ""

            elif line.startswith("relationship: confers_resistance_to_antibiotic"):
                antibiotic = line.split(": ")[1].split("! ")[-1]
                current_antibiotics.append(antibiotic)

            elif line == "[Term]" and current_name:
                all_names = [current_name] + current_synonyms

                for name in all_names:
                    if gene_name.lower() == name.lower():
                        gene_occurrences.append({
                            "Gene Name": current_name,
                            "Matched Name": name,
                            "Gene ID": current_id,
                            "Description": current_description or " ",
                            "Antibiotics": ", ".join(current_antibiotics) if current_antibiotics else pd.NA,
                            "All Synonyms": ", ".join(current_synonyms) if current_synonyms else pd.NA
                        })
                        break

                # Reset
                current_id = None
                current_name = None
                current_description = None
                current_antibiotics = []
                current_synonyms = []

        return gene_occurrences[0] if gene_occurrences else None

    def finding_ARG(self, genes_list: List[str], obo_file: str) -> Tuple[List[dict], List[str]]:
        """
        Identifies genes in the list that are associated with antibiotic resistance.
        """
        res = []
        genes_not_found = []

        for g in genes_list:
            result = self.find_gene_ids(obo_file, g)
            if result:
                res.append(result)
            else:
                genes_not_found.append(g)

        return res, genes_not_found

    def add_antibiotic_by_id(self, gene_id: str, antibiotic: str) -> None:
        """
        Updates the antibiotic information for a gene in the ARGdf DataFrame.
        """
        if hasattr(self, "ARGdf") and not self.ARGdf.empty:
            index = self.ARGdf[self.ARGdf['Gene ID'] == gene_id].index
            if not index.empty:
                self.ARGdf.at[index[0], 'Antibiotics'] = antibiotic
            else:
                print(f"Gene ID '{gene_id}' not found.")
        else:
            raise AttributeError("ARGdf is not initialized.")

    def plot_antibiotic_frequencies(self, df: pd.DataFrame, label_fontsize: int = 10,
                                bar_color: str = 'green', bar_width: float = 0.8) -> None:
        """
        Plots the frequency of antibiotics found in genes as a bar chart.
        """
        if "Antibiotics" not in df.columns:
            raise ValueError("DataFrame must contain an 'Antibiotics' column.")

        df_exploded = df[df['Antibiotics'].notna()]
        df_exploded = df_exploded['Antibiotics'].str.split(',').explode().str.strip()

        # Count occurrences of each antibiotic
        antibiotic_counts = df_exploded.value_counts()

        # Calculate the relative frequency (as a percentage)
        relative_frequencies = (antibiotic_counts / len(df)) * 100

        # Create the bar plot with custom bar width
        plt.figure(figsize=(10, 8))
        sns.barplot(x=relative_frequencies.values, y=relative_frequencies.index,
                    color=bar_color, width=bar_width)

        # Add style and labels
        plt.title('Frequency of Antibiotics Found in Genes', fontsize=16, fontweight='bold')
        plt.xlabel('Frequency of Resistant Genes (%)', fontsize=14)
        plt.ylabel('Antibiotic', fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=label_fontsize)

        # Adjust layout
        plt.tight_layout()

        # Show the plot
        plt.show()
