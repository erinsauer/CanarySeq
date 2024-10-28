# Overview

This repository contains the code, data, and analysis accompanying the paper:

"Diet driven differences in host tolerance are linked to shifts in global gene expression in a common avian host-pathogen system"

Authors: Erin L. Sauer, Carson Stacy, Weston Perrine, Ashley C. Love, Jeffrey A. Lewis, Sarah E. DuRant

Paper currently under review.
Published on bioRxiv, 2024.

## Description
In this study, we explore how diet influences the immune response of canaries following infection with Mycoplasma gallisepticum (MG), a significant avian pathogen. By conducting RNA-Seq analysis, we aim to understand the transcriptomic changes associated with diet composition alongside MG infection.

## Table of Contents
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Data](#data)
- [Results](#results)
- [Contact](#contact)

## Project Structure
```
CanarySeq/
├── data/                   # Raw and processed data
│   ├── raw/                # Raw data files
│   └── processed/          # Processed data files
├── docs/                   # Documentation files
├── scripts/                # Scripts for analysis
├── results/                # Output results (figures, tables)
├── figures/                # Figures generated from the analysis
├── LICENSE                 # License information
├── README.md               # Project overview (this file)
└── .gitignore              # Git ignore file
```

- docs/: Contains supplementary documentation and detailed explanations.
- data/: Includes raw and processed datasets used in the analysis.
- scripts/: R scripts and functions used for data processing and analysis.
- results/: Output figures, tables, and any result files generated by the scripts.
- figures/: Figures generated for the publication.

A description of each file in the repository is available in `docs/`

## Installation
Open Rmarkdown file `docs/CodeSupplement.Rmd` in RStudio and run initial code blocks. Detailed instructions available in `docs/Installation.md`

## Usage
All analysis is available in CodeSupplement.Rmd. Helper scripts used in the main markdown document are available in `scripts`. 

## Data
### Raw Data
- Location: data/raw/
- Description: Contains the original datasets used in the study, including sample metadata and reference genomes.

### Processed Data
- Location: data/processed/
- Description: Data files that have been cleaned and processed for analysis, such as normalized expression matrices.

## Results
### Figures and Tables
- Location: results/ and figures/
- Description: Output figures and tables generated from the analysis scripts.
- Correspondence: Each figure/table filename corresponds to those in the paper for easy reference.
### Supplementary Materials
- Location: results/supplementary/
- Description: Additional data files, extended results, and supplementary analyses supporting the main findings.


## Contact
For questions, please contact:
Erin Sauer
- Email: erinsauer10@gmail.com
Carson Stacy
- Email: clstacy.stat@gmail.com
