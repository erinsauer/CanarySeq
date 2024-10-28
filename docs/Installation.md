# Installation Guide

## Prerequisites
Before you begin, ensure you have met the following requirements:
- Operating System: Linux, macOS, or Windows
- R version: ≥ 4.3.0
- RStudio (optional but recommended)

## Install Required R Packages
The analysis relies on several R packages. Install them using the instructions below.

### Install Bioconductor Packages
Open R or RStudio and run:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "AnnotationHub",
    "clusterProfiler",
    "edgeR",
    "ComplexHeatmap",
    "biomartr",
    "VGAM",
    "AnnotationDbi"
))
```

### Install CRAN Packages
Install the following CRAN packages:

```
install.packages(c(
    "tidyverse",
    "ggplot2",
    "cowplot",
    "ggpubr",
    "eulerr",
    "ComplexUpset",
    "flextable",
    "gt",
    "kableExtra",
    "GGally",
    "data.table",
    "ggrepel",
    "ggbeeswarm",
    "factoextra",
    "here"
))
```

### Optional Packages
Optional but recommended packages:

```
install.packages(c(
    "knitr",
    "rmarkdown"
))
```

## Clone the Repository
Open a terminal and run:

```{bash}
git clone https://github.com/erinsauer/CanarySeq.git
```

Navigate to the project directory:

```{bash}
cd CanarySeq
```

## Set Up the Project Directory
Ensure the project directory has the following structure:

```
CanarySeq/
├── data/
│   ├── raw/                # Raw data files
│   └── processed/          # Processed data files
├── scripts/                # Analysis scripts
├── results/                # Output results
├── figures/                # Generated figures
├── docs/                   # Documentation files
└── ...
```

