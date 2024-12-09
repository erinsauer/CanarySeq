# A function to tidy gene expression matrix
get_tidy_expression <- function(data, meta, from_file=T, type = "tpm") {
  if (type == "tpm") {
    # tpm_df <- readr::read_tsv(data_path,
    #                           show_col_types = FALSE) %>%
    if (from_file) {
      tpm_df <- data.table::fread(data_path, sep="\t") %>%
        tibble() %>%
        dplyr::select(
          -c(
            "ES89_Blue-049_Control-protein_0",
            "ES31_Blue-055_NA_0",
            "ES32_Blue-055_NA_14",
            "ES65_106_NA_0",
            'transcript_id(s)'#, contains("Blue-050")
          )
        )
      meta_subset <- meta %>%
        dplyr::select(c(1,3,4,5,6,7))
      
      tpm_df <<- tpm_df %>%
        as_tibble() %>%
        mutate(gene = gene_id) %>%
        dplyr::select(-gene_id) %>%
        pivot_longer(
          !gene,
          names_to = c("Lane", "Bird", "MG.Diet", "Day"),
          names_sep = "_",
          values_to = type
        ) %>%
        separate(col = MG.Diet, c("MG", "Diet"), remove=FALSE) %>%
        left_join(meta_subset, join_by(Lane==SampID,
                                       Bird==Bird.ID,
                                       Diet==Diet,
                                       Day==Day,
                                       MG.Diet==MG.Diet)) %>%
        left_join(result_BM_merged,
                  #gene_GO_tbl |>
                  #mutate(gene = as.character(entrezgene_id)),
                  by = "gene")
    } else {
      # when data is already in R.
      tpm_df <- data
      
      meta_subset <- meta %>%
        dplyr::select(c(1,3,4,5,6,7))
      
      tpm_df <<- tpm_df %>%
        as_tibble(rownames = "gene_gff") %>%
        pivot_longer(
          !gene_gff,
          names_to = c("Lane", "Bird", "MG.Diet", "Day"),
          names_sep = "_",
          values_to = "tpm"
        ) %>%
        separate(col = MG.Diet, c("MG", "Diet"), remove=FALSE) %>%
        left_join(#entrez_IDs |> dplyr::select(-go_id),
          result_BM_merged,
          by = "gene_gff", multiple = "first") %>%
        left_join(meta_subset, join_by(Lane==SampID,
                                       Bird==Bird.ID,
                                       Diet==Diet,
                                       Day==Day,
                                       MG.Diet==MG.Diet)) 
    }
    
    
    
    
  } else if (type == "count") {
    # count_df <- readr::read_tsv(data_path,
    #                             show_col_types = FALSE) %>%
    if (from_file) {
      count_df <- data.table::fread(data_path, sep="\t") %>%
        tibble() %>%
        dplyr::select(
          -c(
            "ES89_Blue-049_Control-protein_0",
            "ES31_Blue-055_NA_0",
            "ES32_Blue-055_NA_14",
            "ES65_106_NA_0"
          )
        )
      count_df <<-count_df %>%
        as_tibble(rownames = "gene") %>%
        pivot_longer(
          !gene,
          names_to = c("Lane", "Bird", "MG.Diet", "Day"),
          names_sep = "_",
          values_to = type
        ) %>%
        separate(col = MG.Diet, c("MG", "Diet"), remove = FALSE)
    } else {
      error("Error, input of count data object not supported at this time. Convert to TPM.")
    }
    
    
  } else if (length(type)==2){
    # tpm_df <- readr::read_tsv(filepath_rsem_gene_tpm,
    #                           show_col_types = FALSE) %>%
    tpm_df <- data.table::fread(filepath_rsem_gene_tpm, sep="\t") %>%
      tibble() %>%
      dplyr::select(
        -c(
          "ES89_Blue-049_Control-protein_0",
          "ES31_Blue-055_NA_0",
          "ES32_Blue-055_NA_14",
          "ES65_106_NA_0"
        )
      )
    tpm_df <<- tpm_df %>%
      as_tibble(rownames = "gene") %>%
      pivot_longer(
        !gene,
        names_to = c("Lane", "Bird", "MG.Diet", "Day"),
        names_sep = "_",
        values_to = "tpm"
      ) %>%
      separate(col = MG.Diet, c("MG", "Diet"))
    
    # count_df <- readr::read_tsv(filepath_rsem_gene_counts,
    # show_col_types = FALSE) %>%
    count_df <- data.table::fread(filepath_rsem_gene_counts, sep="\t") %>%
      tibble() %>%
      dplyr::select(
        -c(
          "ES89_Blue-049_Control-protein_0",
          "ES31_Blue-055_NA_0",
          "ES32_Blue-055_NA_14",
          "ES65_106_NA_0"
        )
      )
    count_df <<-count_df %>%
      as_tibble(rownames = "gene") %>%
      pivot_longer(
        !gene,
        names_to = c("Lane", "Bird", "MG.Diet", "Day"),
        names_sep = "_",
        values_to = "count"
      ) %>%
      separate(col = MG.Diet, c("MG", "Diet"),remove=FALSE)
    count_df <<- count_df
  } else {
    stop("Error in inputs")
  }
}


# negate %in%
`%ni%` <- Negate(`%in%`)


# A function to take gene name and plot ggplot scatterboxplot
geneviz <- function(gene, tpm_df, subset) {
  ggdat <- 
    if (subset=="CTLvsINF") {
      tpm_df %>%
        filter(gene_id == gene) %>%
        dplyr::filter(MG == "Control" | Day != 0)
    } else if (subset=="DAY0") {
      tpm_df %>%
        filter(gene_id == gene) %>%
        dplyr::filter(Day == 0)
    } else if (subset=="ALL_UNINF") {
      tpm_df %>%
        filter(gene_id == gene) %>%
        dplyr::filter(MG == "Control" | Day == 0)
    } else {
      tpm_df %>%
        filter(gene_id == gene)
    }
  
  # identify where anova stat should go
  anova_y_low <- ggdat %>%
    dplyr::mutate(log2tpm = log2(tpm+1)) %>%
    dplyr::arrange(log2tpm) %>%
    head(1) %>%
    pull(log2tpm) %>% as.numeric()
  
  y_high <- ggdat %>%
    dplyr::mutate(log2tpm = log2(tpm+1)) %>%
    dplyr::arrange(desc(log2tpm)) %>%
    head(1) %>%
    pull(log2tpm) %>% as.numeric()
  
  
  ggdat %>%
    ggplot(aes(
      x = interaction(Diet, MG),
      y = log2(tpm+1),
      #color=Diet,
      shape = as.factor(Day),
      group = interaction(Diet, MG)
    )) +
    # ylim(0,3) +
    geom_boxplot(width = 0.3, outlier.shape = NA) +
    ggbeeswarm::geom_quasirandom(cex = 1.5) +
    ggpubr::stat_compare_means(
      comparisons =
        list(
          c("lipid.Control", "protein.Control"),
          c("lipid.MG", "protein.MG"),
          c("lipid.Control", "lipid.MG"),
          c("protein.Control", "protein.MG")
        ),
      method = "t.test",
      # label = "p.signif",
      label = "p.format",
      # Add pairwise comparisons p-value
      label.y = c(y_high+0.2,y_high+0.2,y_high+0.7,y_high+1.3)
    ) +
    # ylim(-1,8) +
    ylim(anova_y_low-1, y_high+2) +
    ggpubr::stat_compare_means(label.y = anova_y_low-0.5, method = "anova") +
    ggtitle(paste0(gene, " expression - ", subset)) +
    theme_classic()
}


# visualize gene expression over time
timevizMG <- function(gene, tpm_df, include_control=T, diff=T) {
  if (include_control==T) {
    if (diff==T) {
      tpm_df %>%
        filter(gene_id == gene) %>%
        group_by(Bird) %>%
        dplyr::select(c('gene_id', "Bird", "MG", "Diet", 
                        "Day", "tpm")) %>%
        dplyr::mutate(Day = paste0("Day", Day)) %>% 
        # dplyr::filter(tpm>1) %>% # remove samples where gene has very low expression
        pivot_wider(names_from=Day, values_from=tpm) %>%
        dplyr::mutate(diff14v0 = log(Day14)-log(Day0),
                      diff21v0 = log(Day21)-log(Day0)) %>%
        dplyr::select("gene_id", "Bird", "MG", "Diet", 
                      "diff14v0", "diff21v0") %>%
        pivot_longer(c("diff14v0", "diff21v0"),names_to="comparison", values_to="LogDiff") %>%
        ggplot(aes(
          x = interaction(Diet,MG),
          y = LogDiff,
          shape = Diet
        )) +
        # geom_line() +
        geom_boxplot(outlier.shape = NA) +
        geom_point() +
        ggtitle(paste0(gene, " expression changes")) +
        theme_classic() +
        facet_wrap(~comparison)
    } else {
      tpm_df %>%
        filter(gene_id == gene) %>%
        ggplot(aes(
          x = Day,
          y = log2(tpm+1),
          shape = Diet,
          group = Bird
        )) +
        geom_point() +
        geom_line() +
        # geom_boxplot() +
        ggtitle(paste0(gene, " expression over time")) +
        theme_classic() +
        facet_wrap(~MG)
    }
  } else {
    tpm_df %>%
      filter(gene_id == gene,
             MG == "MG") %>% # filter only the MG group
      ggplot(aes(
        x = Day,
        y = log2(tpm+1),
        shape = Diet,
        group = Bird
      )) +
      geom_point() +
      geom_line() +
      # geom_boxplot() +
      ggtitle(paste0(gene, " expression over time in MG infected birds")) +
      theme_classic()
  }
}


# function to convert eulerr to ggplot (https://gist.github.com/danlooo/d23d8bcf8856c7dd8e86266097404ded.js")
# prefer to move to a helper function file.
ggeulerr <- function(combinations, show_quantities = TRUE, show_labels = TRUE, text_size = 2, alpha=0.5, ...) {
  data <-
    eulerr::euler(combinations = combinations) %>%
    plot(quantities = show_quantities) %>%
    pluck("data")
  
  tibble() %>%
    ggplot() +
    ggforce::geom_ellipse(
      data = data$ellipses %>% as_tibble(rownames = "Set"),
      mapping = aes(x0 = h, y0 = k, a = a, b = b, angle = 0, fill = Set),
      alpha = alpha
    ) +
    geom_text(
      data = {
        data$centers %>%
          mutate(
            label = labels %>% map2(quantities, ~ {
              if (!is.na(.x) && !is.na(.y) && show_labels) {
                paste0(.x, "\n", sprintf(.y, fmt = "%.4g"))
              } else if (!is.na(.x) && show_labels) {
                .x
              } else if (!is.na(.y)) {
                .y
              } else {
                ""
              }
            })
          )
      },
      mapping = aes(x = x, y = y, label = label), size = text_size
    ) +
    theme(panel.grid = element_blank()) +
    coord_fixed() +
    scale_fill_hue()
}



# Function to clean and mutate the data
clean_and_mutate_data <- function(data, sig.cutoff=0.05, logFC.cutoff=0.575) {
  data %>%
    mutate(delabel = case_when(
      !is.na(gene_name) ~ gene_name,
      TRUE ~ paste0("LOC",gene)
    )) %>%
    mutate(delabel = case_when(
      FDR < sig.cutoff ~ delabel,
      abs(logFC) > logFC.cutoff ~ delabel,
      TRUE ~ NA
    )) %>%
    mutate(diffexpressed = case_when(
      FDR > sig.cutoff ~ "Neither",
      logFC < 0 ~ "Down",
      logFC > 0 ~ "Up"
    )) %>%
    mutate(diffexpressed = factor(diffexpressed, levels = c("Neither", "Down", "Up"))) %>%
    arrange(diffexpressed)
}

# Function to create the plot
create_gene_plot <- function(data, title, sig.cutoff=0.05, logFC.cutoff=0.575, show.names=T) {
  ggplot(data, aes(x = logCPM, 
                   y = logFC, 
                   col = diffexpressed, 
                   label = delabel#, size=AveExpr
  )) +
    geom_point(size=0.5) +
    theme_classic() +
    scale_color_manual("Differentially\n Expressed",values = c("black", "blue", "red")) +
    #geom_vline(xintercept = c(-logFC.cutoff, logFC.cutoff), col = "red") +
    #geom_hline(yintercept = -log10(sig.cutoff), col = "red") +
    {if (show.names) ggrepel::geom_text_repel(size=2, na.rm=TRUE)} +
    scale_size_area() +
    # ggsci::scale_color_aaas() +
    ggtitle(title) + 
    theme(legend.position="none",
          text = element_text(size=10),
          axis.text = element_text(size=10))
}