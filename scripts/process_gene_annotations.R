library(biomartr)

filepath_rsem_gene_counts <- here::here("data", "rsem.merged.gene_counts.tsv")

# read in rsem.merged.gene_counts.tsv
g.counts.l <- data.table::fread(filepath_rsem_gene_counts, sep="\t") %>% as_tibble()

# biomartr::getCollection( db = "refseq", organism = "Serinus canaria")

scan_gff <- biomartr::read_gff("./_db_downloads/collections/refseq/Serinus_canaria/Serinus_canaria_genomic_refseq.gff.gz") %>%
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

