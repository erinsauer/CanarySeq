library(biomartr)

filepath_rsem_gene_counts <- here::here("data", "rsem.merged.gene_counts.tsv")

# read in rsem.merged.gene_counts.tsv
g.counts.l <- data.table::fread(filepath_rsem_gene_counts, sep="\t") %>% as_tibble()

# biomartr::getCollection( db = "refseq", organism = "Serinus canaria")
here()
scan_gff <- biomartr::read_gff(here("_db_downloads/collections/refseq/Serinus_canaria/Serinus_canaria_genomic_refseq.gff.gz")) %>%
  filter(type=="gene", source == "Gnomon") %>%
  # create new column pulling out gene name
  mutate(gene_common_name = str_extract(attribute, "ID=gene-([^;]+)"),
         gene_ID = str_extract(attribute, "Dbxref=GeneID:([^;]+)")) %>%
  mutate(gene_common_name = str_remove(gene_common_name, "ID=gene-") ,
         gene_ID = str_remove(gene_ID, "Dbxref=GeneID:"))

# get list of gene IDs
gene_set <- unique(scan_gff$gene_ID)

# get result BM
result_BM <- biomartr::biomart( genes      = gene_set, # genes were retrieved using biomartr::getGenome()
                                mart       = "ENSEMBL_MART_ENSEMBL", # marts were selected with biomartr::getMarts()
                                dataset    = "scanaria_gene_ensembl", # datasets were selected with biomartr::getDatasets()
                                attributes = c(#'ensembl_transcript_id',
                                  # 'ensembl_peptide_id',
                                  'ensembl_gene_id',
                                  # 'uniprot_gn_symbol',
                                  'wikigene_name',
                                  #'entrezgene_id',
                                  'entrezgene_description',
                                  # '',
                                  'external_gene_name'
                                  # 'go_id',
                                  # 'description'
                                ), # attributes were selected with biomartr::getAttributes()
                                filters    = "entrezgene_id"
) 


# include gff annotation data
gene_anno <- scan_gff %>%
  dplyr::select(gene_ID, gene_common_name, attribute) %>%
  dplyr::rename(entrezgene_id = gene_ID,
                gene_gff = gene_common_name) %>%
  # split attribute column based on ; separator, column name comes before =
  separate_wider_delim(attribute, delim = ";", names= c("ID", "LOC_id", "common_name", "description",
                                                        "endrange", "genbank_key", "gene_name", "gene_type"),
                       too_many = "merge",
                       too_few = "align_start") %>%
  dplyr::select(entrezgene_id, gene_gff, description) %>%
  mutate(description = str_remove(description, "description="))



# Generate count matrix
gene_expression_matrix <- g.counts.l %>% #View()
  # remove globin genes
  dplyr::filter(gene_id != "LOC103824466" ) %>%
  dplyr::filter(gene_id != "LOC103817748" ) %>%
  dplyr::filter(gene_id != "LOC103817749" ) %>%
  dplyr::filter(gene_id != "LOC115485052" ) %>%
  dplyr::filter(gene_id != "LOC103817748" ) %>%
  # Add rownames
  column_to_rownames("gene_id") %>%
  # remove unneeded column
  dplyr::select(-`transcript_id(s)`) %>% 
  # convert to matrix for edgeR
  as.matrix()

# create merged biomart annotation
result_BM_merged <- result_BM %>%
  mutate(entrezgene_id = as.character(entrezgene_id)) %>%
  full_join(gene_anno, by = "entrezgene_id") %>%
  mutate(gene_name = case_when(
    str_detect(gene_gff, "LOC1") & !str_detect(external_gene_name, "LOC1") & !is.na(external_gene_name) & external_gene_name != "" ~ external_gene_name,
    str_detect(gene_gff, "LOC1") & !str_detect(wikigene_name, "LOC1") & !is.na(wikigene_name) ~ wikigene_name,
    #str_detect(gene_gff, "LOC10") ~ gene_gff,
    TRUE ~ gene_gff#paste0("LOC",entrezgene_id)
  ),
  gene_tag = case_when(
    paste0("LOC",entrezgene_id) %in% rownames(gene_expression_matrix) ~ paste0("LOC",entrezgene_id),
    gene_name %in% rownames(gene_expression_matrix) ~ gene_name,
    TRUE ~ NA),
  gene = case_when(
    paste0("LOC",entrezgene_id) %in% rownames(gene_expression_matrix) ~ entrezgene_id,
    gene_name %in% rownames(gene_expression_matrix) ~ gene_name,
    TRUE ~ NA)
  ) %>%
  # remove formatting issue for commas in description.
  mutate(description = str_replace_all(description, "%2", ",")) %>%
  mutate(description = str_replace_all(description, ",C", ",")) 

write_tsv(result_BM_merged, file = here::here("data/annotation/results_BM_merged.tsv"))

