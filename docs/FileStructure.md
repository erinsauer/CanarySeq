# File Structure

Below are descriptions for each file in this repository

```
CanarySeq/
├── CanarySeq.Rproj                # R project file
├── README.md                      # Project overview (this file)
├── Reference/                     # Reference genome and annotation files
│   ├── GCFcanaryGenotype.gtf
│   └── GCFcanaryGenotype.gtf.gz
├── data/                          # Raw and processed data
│   ├── raw/                       # Raw data files
│   │   └── sample_metadata/       # Sample metadata files
│   │       ├── canary_antibody.csv
│   │       ├── canary_eyescores.csv
│   │       ├── canary_eyescores_subset.csv
│   │       ├── canary_loads.csv
│   │       ├── metadata.tsv
│   │       ├── metadata_phenotypes.tsv
│   │       └── metadata_tolerance_scores.tsv
│   └── processed/                 # Processed data files
│       ├── annotation/            # Annotation files
│       │   ├── _db_downloads/     # Downloaded databases
│       │   │   └── collections/
│       │   │       └── refseq/
│       │   │           └── Serinus_canaria/
│       │   │               ├── Serinus_canaria_assembly_stats_refseq.txt
│       │   │               ├── Serinus_canaria_cds_from_genomic_refseq.fna.gz
│       │   │               ├── Serinus_canaria_genomic_refseq.fna.gz
│       │   │               ├── Serinus_canaria_genomic_refseq.gff.gz
│       │   │               ├── Serinus_canaria_protein_refseq.faa.gz
│       │   │               ├── Serinus_canaria_rna_from_genomic_refseq.fna.gz
│       │   │               └── doc/
│       │   │                   ├── doc_Serinus canaria_db_refseq.txt
│       │   │                   └── doc_Serinus_canaria_db_refseq.tsv
│       │   └── results_BM_merged.tsv   # Merged annotation results
│       ├── rsem.merged.gene_counts.tsv # Merged gene counts from RSEM
│       └── rsem.merged.gene_tpm.tsv    # Merged TPM values from RSEM
├── docs/                          # Documentation files
│   ├── CodeSupplement.Rmd         # Code supplement for the paper
│   ├── Installation.md            # Installation guide
│   ├── FileStructure.md           # File directory guide
│   └── Usage.md                   # Usage instructions
├── reports/                       # Reports and analysis outputs
│   └── mapping_statistics/        # Mapping statistics and MultiQC reports
│       ├── multiqc_data/          # Data files generated by MultiQC
│       └── multiqc_report.html    # MultiQC summary report
├── results/                       # Output results (figures, tables)
│   ├── figures/                   # Figures generated from the analysis
│   │   ├── DEGenesBothCorrInfDietLogFC.png
│   │   ├── DEGenesBothCorrlipInfDietLogFC.png
│   │   ├── DEGenesBothCorrprotInfDietLogFC.png
│   │   ├── GO_term_enrichment_comparison.png
│   │   ├── fig.s1.png
│   │   ├── fig.s2.png
│   │   ├── fig.s3.png
│   │   ├── fig.s4.DE_genes_multinode_ES_combined_upsetplot.png
│   │   ├── fig.s5.png
│   │   ├── fig.s6.png
│   │   ├── fig.s7.png
│   │   ├── fig2.png
│   │   ├── fig3.png
│   │   ├── fig4_phenotype.png
│   │   ├── fig4a_venn_plot_phenotype_fig.png
│   │   ├── logFC_density_DE_in_diet_fig.png
│   │   ├── logFC_heatmap.png
│   │   ├── top_genes_ES_inf.pdf
│   │   ├── top_genes_ES_inf.png
│   │   ├── top_genes_all.png
│   │   ├── top_genes_ixndietinf.png
│   │   ├── top_genes_ixndietinf_colors.png
│   │   └── venn_plot_interaction_fig.png
│   └── tables/                    # Tables of results
│       ├── DE_genes/              # Differential expression results
│       │   ├── DE_genes_EyeScore.csv       # DE genes related to Eye Score
│       │   ├── DE_genes_diet.csv           # DE genes for Diet contrast
│       │   ├── DE_genes_inf_all.csv        # DE genes for Infection (all)
│       │   ├── DE_genes_inf_lipid.csv      # DE genes for Infection (Lipid diet)
│       │   ├── DE_genes_inf_protein.csv    # DE genes for Infection (Protein diet)
│       │   ├── DE_genes_infxdiet.csv       # DE genes for Infection × Diet interaction
│       │   └── DE_genes_pathogen_load.csv  # DE genes related to Pathogen Load
│       ├── all_genes/             # All genes with estimates
│       │   ├── all_genes_combined_estimates.csv  # Combined estimates for all genes
│       │   └── tbl_ordinal_regression_results.tsv # Ordinal regression results
│       └── functional_enrichments/ # Functional enrichment results
│           ├── GO_ES_enrichments_DE.csv    # GO enrichments for Eye Score DE genes
│           ├── GO_enrichments.csv          # GO enrichments for various contrasts
│           └── KEGG_enrichments.csv        # KEGG pathway enrichments
└── scripts/                       # Scripts for analysis
    ├── canary_data_processing_functions.R  # Custom functions for data processing
    ├── generate_metadata.R                 # Script to generate metadata files
    └── process_gene_annotations.R          # Script to process gene annotations
```

## Description of Key Directories and Files:

- CanarySeq.Rproj: R project file for organizing the workspace.
- README.md: Project overview and instructions.
- Reference/: Contains reference genome and annotation files used in the analysis.
- data/: Directory containing all data files.
  - raw/: Raw data files.
    - sample_metadata/: Original sample metadata and phenotype data.
  - processed/: Data files that have been processed and are ready for analysis.
    - annotation/: Processed gene annotation files and associated downloads.
    - rsem.merged.gene_counts.tsv: Merged gene counts from RSEM.
    - rsem.merged.gene_tpm.tsv: Merged TPM values from RSEM.
- docs/: Documentation files for installation and usage.
  - CodeSupplement.Rmd: R Markdown file containing the code supplement for the paper.
  - Installation.md: Installation instructions.
  - Usage.md: Guide on how to use the scripts and data.
- reports/: Contains reports and outputs from data quality assessments.
  - mapping_statistics/: Mapping statistics and MultiQC reports generated from sequencing data.
    - multiqc_data/: Data files generated by MultiQC.
    - multiqc_report.html: Interactive HTML report summarizing mapping statistics.
- results/: Output results from analyses.
  - figures/: Figures generated from the analysis, corresponding to those in the paper.
  - tables/: Tables of results.
    - DE_genes/: Differential expression results for various contrasts.
    - all_genes/: Combined estimates and regression results for all genes.
    - functional_enrichments/: Results from functional enrichment analyses (GO and KEGG).
- scripts/: R scripts used for data processing, analysis, and figure generation.
