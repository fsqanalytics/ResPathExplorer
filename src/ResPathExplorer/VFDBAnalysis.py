import os
import gzip
import shutil
import requests
import re
import pandas as pd
from typing import List, Optional
import seaborn as sns
import matplotlib.pyplot as plt


class VFDBAnalysis:
    """
    A class to download, parse, and analyze data from the VFDB (Virulence Factors Database).
    """

    def __init__(self, db_dir: str = "db"):
        """
        Initializes the VFDBAnalyzer.

        Args:
            db_dir (str): Directory where VFDB files will be stored.

        Raises:
            ValueError: If db_dir is not a valid string path.
        """
        if not isinstance(db_dir, str) or not db_dir.strip():
            raise ValueError("db_dir must be a non-empty string.")

        self.db_dir = db_dir
        self.fasta_url = "https://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz"
        self.xls_url = "https://www.mgc.ac.cn/VFs/Down/VFs.xls.gz"
        self.fasta_path = os.path.join(db_dir, "VFDB_setA_pro.fas.gz")
        self.xls_gz_path = os.path.join(db_dir, "VFs.xls.gz")
        self.xls_path = os.path.join(db_dir, "VFs.xls")
        self.df_genes: Optional[pd.DataFrame] = None

        os.makedirs(self.db_dir, exist_ok=True)

    def download_data(self) -> None:
        """
        Downloads and decompresses VFDB FASTA and XLS files if not already present.
        """
        # Download FASTA
        if not os.path.exists(self.fasta_path):
            print("Downloading FASTA file...")
            response = requests.get(self.fasta_url, timeout=30)
            response.raise_for_status()
            with open(self.fasta_path, 'wb') as f_out:
                f_out.write(response.content)
            print("FASTA file downloaded.")
        else:
            print("FASTA file already exists.")

        # Download and decompress XLS
        if not os.path.exists(self.xls_path):
            print("Downloading XLS file...")
            with requests.get(self.xls_url, stream=True, timeout=30) as r:
                r.raise_for_status()
                with open(self.xls_gz_path, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)

            print("🗜️ Decompressing XLS file...")
            with gzip.open(self.xls_gz_path, 'rb') as f_in:
                with open(self.xls_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print("XLS file decompressed.")
        else:
            print("XLS file already exists.")

    def parse_fasta_entries(self) -> List[dict]:
        """
        Parses FASTA headers to extract gene metadata.
        """
        if not os.path.exists(self.fasta_path):
            raise FileNotFoundError(f"FASTA file not found at {self.fasta_path}")

        entries = []
        with gzip.open(self.fasta_path, 'rt') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()
                    gene_matches = re.findall(r'\(([^()]+)\)', header)
                    gene = gene_matches[1] if len(gene_matches) > 1 else None
                    desc_part = header.split(')', 2)[-1].split('[')[0].strip()
                    description = desc_part
                    bracket_split = header.split('[', 1)
                    after_bracket = bracket_split[1] if len(bracket_split) > 1 else ""
                    vf_id_match = re.search(r'\(([^()]+)\)', after_bracket)
                    vf_id = vf_id_match.group(1) if vf_id_match else None
                    m = re.search(r'\[.*?\((?:[^()]*)\)\s*-\s*(.*?)\s*\(.*?\)\]', header)
                    category = m.group(1).strip() if m else None
                    org_match = re.findall(r'\[([^\[\]]+)\]', header)
                    organism = org_match[-1] if org_match else None

                    entries.append({
                        'Gene_Name': gene,
                        'Description': description,
                        'Functional category': category,
                        'Bacteria': organism,
                        'VFID': vf_id
                    })

        if not entries:
            raise ValueError("No entries parsed from the FASTA file.")

        return entries

    def load_and_process(self) -> pd.DataFrame:
        """
        Loads and merges FASTA and XLS data to create a complete gene DataFrame.
        """
        self.download_data()

        fasta_entries = self.parse_fasta_entries()
        df_fasta = pd.DataFrame(fasta_entries)

        if not os.path.exists(self.xls_path):
            raise FileNotFoundError(f"XLS file not found at {self.xls_path}")

        df_xls = pd.read_excel(self.xls_path, header=1)

        required_cols = {'VFID', 'VF_Name', 'Function'}
        if not required_cols.issubset(df_xls.columns):
            raise ValueError(f"Missing required columns in XLS file: {required_cols - set(df_xls.columns)}")

        # Fix missing VFID references
        filter_vf = df_fasta['VFID'].str.startswith("VF", na=False)
        df_filtered = df_fasta[~filter_vf]
        unique_ids = df_filtered["VFID"].unique()

        vf_map = {}
        for uid in unique_ids:
            pattern = f"({uid})"
            match = df_xls[df_xls['VF_Name'].apply(lambda x: pattern in str(x))]
            if not match.empty:
                vf_map[uid] = match["VFID"].iloc[0]

        for k, v in vf_map.items():
            df_fasta.loc[df_fasta['VFID'] == k, 'VFID'] = v

        df_xls_reduced = df_xls[['VFID', 'VF_Name', 'Function']]
        self.df_genes = df_fasta.merge(df_xls_reduced, on='VFID', how='left')

        if self.df_genes.empty:
            raise ValueError("The merged DataFrame is empty.")

        return self.df_genes

    def search_virulence_genes(self, important_genes: List[str], bacteria: str) -> pd.DataFrame:
        """
        Filters virulence genes based on gene name and bacterial species.
        """
        if self.df_genes is None:
            raise ValueError("Gene data not loaded. Call load_and_process() first.")
        if not isinstance(important_genes, list) or not important_genes or not all(
                isinstance(gene, str) and gene.strip() for gene in important_genes):
            raise ValueError("important_genes must be a non-empty list of gene names.")
        if not isinstance(bacteria, str) or not bacteria.strip():
            raise ValueError("bacteria must be a non-empty string.")

        gene_names = {gene.strip().casefold() for gene in important_genes}
        mask = (
            self.df_genes['Bacteria'].str.contains(bacteria.strip(), case=False, na=False, regex=False)
            & self.df_genes['Gene_Name'].str.strip().str.casefold().isin(gene_names)
        )
        # Keep distinct strain records for provenance, without repeating input genes.
        return self.df_genes.loc[mask].drop_duplicates().reset_index(drop=True)

    def plot_virulence_factors_percentage(self, df: pd.DataFrame,
                                          bacteria_name: str,
                                          show_all_categories: bool = True) -> None:
        """
        Plot the distribution of unique gene-category associations across strains.

        Each gene counts once per category when ``Gene_Name`` is available.
        Without gene names, each input row is an observation. Missing categories
        are shown as Unclassified; categories outside the standard list remain visible.
        """
        if df is None or df.empty:
            raise ValueError("DataFrame cannot be empty.")
        if "Functional category" not in df.columns:
            raise ValueError("Column 'Functional category' is missing from df.")
        if not isinstance(bacteria_name, str) or not bacteria_name.strip():
            raise ValueError("bacteria_name must be a non-empty string.")

        categories = [
            'Exotoxin', 'Nutritional/Metabolic factor', 'Biofilm',
            'Immune modulation', 'Regulation', 'Adherence',
            'Effector delivery system', 'Exoenzyme', 'Motility', 'Invasion',
            'Post-translational modification', 'Others', 'Stress survival',
            'Antimicrobial activity/Competitive advantage'
        ]

        observations = df.copy()
        observations["Functional category"] = (
            observations["Functional category"].fillna("").str.strip().replace("", "Unclassified")
        )
        if "Gene_Name" in observations:
            observations["Gene_Name"] = observations["Gene_Name"].str.strip().str.casefold()
            named = observations["Gene_Name"].notna() & observations["Gene_Name"].ne("")
            observations = pd.concat([
                observations.loc[named].drop_duplicates(["Gene_Name", "Functional category"]),
                observations.loc[~named]
            ], ignore_index=True)
        actual_counts = observations["Functional category"].value_counts(normalize=True) * 100
        categories += sorted(set(actual_counts.index) - set(categories))
        freq = pd.DataFrame({'Functional category': categories})
        freq["Percentage"] = freq["Functional category"].map(actual_counts).fillna(0)

        if not show_all_categories:
            freq = freq[freq["Percentage"] > 0].reset_index(drop=True)

        plt.figure(figsize=(14, 6))
        ax = sns.barplot(data=freq, x="Functional category", y="Percentage", color="steelblue")
        plt.title(f"Virulence Gene-Category Distribution - {bacteria_name}")
        plt.ylabel("Gene-category associations (%)" if "Gene_Name" in df else "Observations (%)")
        plt.xticks(rotation=45, ha='right')
        plt.ylim(0, max(freq["Percentage"].max() + 5, 10))
        plt.tight_layout()

        for p in ax.patches:
            height = p.get_height()
            ax.annotate(f'{height:.1f}%',
                        (p.get_x() + p.get_width() / 2., height / 2 if height > 0 else 0.5),
                        ha='center', va='center', fontsize=9,
                        color='white' if height > 0 else 'black', fontweight='bold')
        plt.show()
