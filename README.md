# ResPathExplorer

ResPathExplorer is a Python-based library designed to streamline transcriptomic data analysis by integrating the Kyoto Encyclopedia of Genes and Genomes (KEGG) and the Comprehensive Antibiotic Resistance Database (CARD). It enables the identification and visualization of enriched metabolic pathways and antibiotic resistance gene mapping.

# Features

- Interpret transcriptomic data.

- Identify metabolic pathway enrichment.

- Map prevalence of antibiotic resistance genes.

- Integrated visualization with KEGG and CARD databases.

# Installation

pip install respathexplorer

# Usage

import respathexplorer as rpe

## Analyze transcriptomic data
data = rpe.load_data('your_data.csv')

## Perform pathway analysis
pathways = rpe.analyze_pathways(data)

## Visualize results
rpe.visualize_pathways(pathways)

# Applications

- Antimicrobial resistance monitoring

- Food safety and sustainability research

- Pathogenic bacteria characterization in agro-food systems

# Authors

Developed by Laís Carvalho, Ursula Gonzales-Barron, and Vasco Cadavez at CIMO, Instituto Politécnico de Bragança, Portugal.
