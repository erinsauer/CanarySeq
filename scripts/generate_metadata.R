# assign file path locaiton of rsem.merged.gene_counts.tsv file
filepath_rsem_gene_counts <- here::here("data", "rsem.merged.gene_counts.tsv")

# read in rsem.merged.gene_counts.tsv
g.counts.l <- data.table::fread(filepath_rsem_gene_counts, sep="\t") %>% as_tibble()

# remove probelmatic samples
g.counts.l <- g.counts.l %>%
  dplyr::select(-c(
    "ES89_Blue-049_Control-protein_0",
    "ES31_Blue-055_NA_0",
    "ES32_Blue-055_NA_14",
    "ES65_106_NA_0" )) %>% 
  # dplyr::select(-contains("Blue-050")) %>%
  # change diet group of bird SEM19-1481 from MG-lipid to MG-protein
  dplyr::rename(all_of(
    # c(`ES64_SEM19-1481_MG-lipid_21`="ES64_SEM19-1481_MG-protein_21",
    #   `ES52_SEM19-1481_MG-lipid_14`="ES52_SEM19-1481_MG-protein_14",
    #   `ES83_SEM19-1481_MG-lipid_0`="ES83_SEM19-1481_MG-protein_0"))
    # )
    c(`ES64_SEM19-1481_MG-protein_21`="ES64_SEM19-1481_MG-lipid_21",
      `ES52_SEM19-1481_MG-protein_14`="ES52_SEM19-1481_MG-lipid_14",
      `ES83_SEM19-1481_MG-protein_0`="ES83_SEM19-1481_MG-lipid_0"))
  )

# create tibble metadata object, convert column names to sample ID column
metadata <-
  tibble(SampID = colnames(g.counts.l[, c(3:ncol(g.counts.l))])) %>%
  separate(
    SampID,
    c("Sample.ID", "Bird.ID",
      "MG.Diet", "Day"),
    "_",
    extra = "warn",
    fill = "warn",
    remove = FALSE
  ) %>%
  separate(
    MG.Diet,
    c("MG.Status", "Diet"),
    sep = "[-]",
    extra = "warn",
    fill = "warn",
    remove = FALSE
  ) %>%
  mutate(
    MG.Diet.Day = paste0(MG.Diet, ".", Day),
    MG.Day = paste0(MG.Status, ".", Day),
    Diet.Day = paste0(Diet, ".", Day)
  ) %>%
  mutate(MG.active = if_else(Day < 1, "NO",
                             if_else(MG.Status == "Control",
                                     "NO", "YES"))) %>%
  mutate(
    Exposure.MG = paste0(MG.active, ".", MG.Status),
    Exposure.Diet = paste0(MG.active, ".", Diet),
    Exposure.Diet.Day = paste0(Exposure.Diet, ".", Day)
  )


# read in phenotype data

## pathogen loads
canary_loads <-
  # only birds with both phenotype and genotype data
  # readr::read_csv(here("Geno.x.Pheno/canary_loads.csv")) %>% 
  data.table::fread(here("sample_metadata/canary_loads.csv"),sep=",") %>% 
  tibble() %>%
  # mutate(measure = "eye_score", value = ES) %>% # allow data to be merged with tidygeneExpr
  mutate(Bird = stringr::str_replace(ID, " ", ".")) %>%
  mutate(Bird = stringr::str_replace(Bird, "Orange", "orange")) %>%
  mutate(Bird = stringr::str_replace(Bird, "- ", "")) %>%
  # mutate(Bird = stringr::str_replace(Bird, "-", ".")) %>%
  mutate(Bird = stringr::str_replace(Bird, "LD.19 239", "LD19.239")) %>%
  mutate(Bird = stringr::str_replace(Bird, "Blue.045", "Blue.45")) %>%
  dplyr::filter(Bird %in% unique(metadata$Bird.ID)) %>%
  dplyr::select(-ID)

## semi-quantitative eyescore measures
eyescores <-
  # readr::read_csv(
  data.table::fread(
    here("sample_metadata/canary_eyescores.csv"),
    sep=","
    # col_types = cols(Day = col_integer(),
    #                  Date = col_character())
  ) %>%
  tibble() %>%
  mutate(Day = as.integer(Day),
         Date = as.character(Date)) %>%
  # correcting for errors in manual data entry
  mutate(Bird = stringr::str_replace(ID, " ", "-")) %>%
  mutate(Bird = stringr::str_replace(Bird, "Orange", "orange")) %>%
  mutate(Bird = stringr::str_replace(Bird, "- ", "")) %>%
  mutate(Bird = stringr::str_replace(Bird, "\\.", "-")) %>%
  mutate(Bird = stringr::str_replace(Bird, "Nlue", "Blue")) %>%
  mutate(Bird = stringr::str_replace(Bird, "LD-19 239", "LD19-239")) %>%
  mutate(Bird = stringr::str_replace(Bird, "Blue-045", "Blue-45")) %>%
  dplyr::filter(Bird %in% unique(metadata$Bird.ID)) %>%
  dplyr::mutate(day = replace(Day, which(Day<0), 0)) %>% # adjust day -1 to day 0
  dplyr::select(Bird, day, ES)

# inspired from https://support.bioconductor.org/p/p132992/

# generate columns of phenotypic metadata I'd like to have
# pathogen loads
day3_loads <- canary_loads %>%
  filter(day==3) %>%
  mutate(logLoad_day3 = logLoad, r_logLoad_day3=logLoad/max(logLoad)) %>%
  dplyr::select(c(Bird, logLoad_day3, r_logLoad_day3))

day7_loads <- canary_loads %>%
  filter(day==7) %>%
  mutate(logLoad_day7 = logLoad, r_logLoad_day7=logLoad/max(logLoad)) %>%
  dplyr::select(c(Bird, logLoad_day7, r_logLoad_day7))

day14_loads <- canary_loads %>%
  filter(day==14) %>%
  mutate(logLoad_day14 = logLoad, r_logLoad_day14=logLoad/max(logLoad)) %>%
  dplyr::select(c(Bird, logLoad_day14, r_logLoad_day14))

day21_loads <- canary_loads %>%
  filter(day==21) %>%
  mutate(logLoad_day21 = logLoad, r_logLoad_day21=logLoad/max(logLoad)) %>%
  dplyr::select(c(Bird, logLoad_day21, r_logLoad_day21))

#AUC estimates
day14_AUC_loads <- canary_loads %>%
  filter(day<=14) %>%
  group_by(Bird) %>%
  summarise(
    auc = MESS::auc(day,logLoad, type = "spline")
  ) %>%
  mutate(logLoad_AUC_14 = auc, r_logLoad_AUC_14=auc/max(auc)) %>%
  dplyr::select(c(Bird, logLoad_AUC_14, r_logLoad_AUC_14))

day21_AUC_loads <- canary_loads %>%
  filter(day<=21) %>%
  group_by(Bird) %>%
  summarise(
    auc = MESS::auc(day,logLoad, type = "spline")
  ) %>%
  mutate(logLoad_AUC_21 = auc, r_logLoad_AUC_21=auc/max(auc)) %>%
  dplyr::select(c(Bird, logLoad_AUC_21, r_logLoad_AUC_21))


load_phenotypes <- purrr::reduce(list(day3_loads,day7_loads,
                                      day14_loads, day21_loads,
                                      day14_AUC_loads, day21_AUC_loads), 
                                 dplyr::left_join, by = 'Bird')

## now do the same for eyescore data
day_eyes <- eyescores %>%
  filter(., day %in% c(0,14,21)) %>%
  mutate(day=as.character(day)) %>%
  dplyr::select(c(Bird, day, ES))

day3_eyes <- eyescores %>%
  filter(.,day==3) %>%
  mutate(ES_day3 = ES, r_ES_day3=ES/max(ES)) %>%
  dplyr::select(c(Bird, ES_day3, r_ES_day3))

day7_eyes <- eyescores %>%
  filter(.,day==7) %>%
  mutate(ES_day7 = ES, r_ES_day7=ES/max(ES)) %>%
  dplyr::select(c(Bird, ES_day7, r_ES_day7))

day14_eyes <- eyescores %>%
  filter(.,day==14) %>%
  mutate(ES_day14 = ES, r_ES_day14=ES/max(ES)) %>%
  dplyr::select(c(Bird, ES_day14, r_ES_day14))

day21_eyes <- eyescores %>%
  filter(.,day==21) %>%
  mutate(ES_day21 = ES, r_ES_day21=ES/max(ES)) %>%
  dplyr::select(c(Bird, ES_day21, r_ES_day21))

#AUC estimates
day7_AUC_eyes <- eyescores %>%
  filter(.,day<=7) %>%
  group_by(Bird) %>%
  summarise(
    auc = MESS::auc(day,ES, type = "spline")
  ) %>%
  mutate(ES_AUC_7 = auc, r_ES_AUC_7=auc/max(auc)) %>%
  dplyr::select(c(Bird, ES_AUC_7, r_ES_AUC_7))

day14_AUC_eyes <- eyescores %>%
  filter(.,day<=14) %>%
  group_by(Bird) %>%
  summarise(
    auc = MESS::auc(day,ES, type = "spline")
  ) %>%
  mutate(ES_AUC_14 = auc, r_ES_AUC_14=auc/max(auc)) %>%
  dplyr::select(c(Bird, ES_AUC_14, r_ES_AUC_14))

day21_AUC_eyes <- eyescores %>%
  filter(.,day<=21) %>%
  group_by(Bird) %>%
  summarise(
    auc = MESS::auc(day,ES, type = "spline")
  ) %>%
  mutate(ES_AUC_21 = auc, r_ES_AUC_21=auc/max(auc)) %>%
  dplyr::select(c(Bird, ES_AUC_21, r_ES_AUC_21))



eye_phenotypes <- purrr::reduce(list(day3_eyes,day7_eyes,
                                     day14_eyes, day21_eyes,
                                     day7_AUC_eyes,
                                     day14_AUC_eyes, day21_AUC_eyes), 
                                dplyr::left_join, by = 'Bird')

## add phenotype data to existing metadata df
metadata_phenotypes <- metadata %>%
  left_join(load_phenotypes, by = join_by('Bird.ID' == "Bird")) %>%
  left_join(eye_phenotypes, by = join_by('Bird.ID' == "Bird")) %>%
  left_join(day_eyes, by = join_by('Bird.ID' == "Bird", "Day" == "day")) %>%
  mutate_all(., ~replace_na(.,0)) %>%
  distinct() %>% # to get rid of duplicated row
  mutate(ES_minMax = ES/6) %>% # normalized to 0 to 1
  # View()
  mutate(logLoad = case_when(
    Day == 0 ~ logLoad_day3,
    Day == 14 ~ logLoad_AUC_14/14,
    Day == 21 ~ logLoad_AUC_21/21,
    TRUE ~ NA
  ),
  logLoad_alt = case_when(
    Day == 0 ~ 0,
    Day == 14 ~ logLoad_day14,
    Day == 21 ~ logLoad_day21,
    TRUE ~ NA
  ),
  ES_day0_analysis = case_when(
    Day == 0 ~ ES_day3/6,
    TRUE ~ ES_minMax
  )
  ) %>%
  mutate(tolerance_loadoverES = (1+logLoad_alt)/(ES_day0_analysis+1)-1, 
         tolerance_erin = 1-((ES)/(logLoad_alt)),
         tolerance = 1-(ES/logLoad_alt))

# create new columnts, create rownames
metadata_phenotypes_ixn <- metadata_phenotypes |>
  group_by(Bird.ID) |>
  mutate(tolerance_median = median(tolerance, na.rm = TRUE)) |>
  ungroup()|>
  mutate(
    tolerance_cat = case_when(
      tolerance_median > 0 ~ "tolerant",
      tolerance_median < 0 ~ "susceptible",
      TRUE ~ "intermediate"
    )) %>%
  mutate(tolerance_cat = factor(tolerance_cat, levels = c("susceptible", "intermediate",#"intermediate1",
                                                          "tolerant"))) %>%
  mutate(diet.tolerance = paste0(Diet, ".",tolerance_cat)) %>%
  # mutate(SampID = str_replace_all(SampID,"\\.","-")) %>%
  column_to_rownames("SampID")

metadata_phenotypes_ixn  %>%
  rownames_to_column("Sample.Name") %>%
  dplyr::select("Sample.Name", "Sample.ID", "Bird.ID", "MG.Status", "MG.active", "Diet", "Day", "ES_day0_analysis", "logLoad", "tolerance", "tolerance_median", "tolerance_cat") %>%
  rename("ES_day0_analysis"= "ES" ,  "tolerance"="tolerance.day" ,  "tolerance_median" = "tolerance.bird",
         "tolerance_cat" = "tolerance.category" 
  ) %>%
  mutate(Day = as.numeric(Day)) %>% 
  # left_join(eyescores[,-3], by = c("Bird.ID"="Bird", "Day" = "day"))
  write_tsv(., here("sample_metadata/metadata_tolerance_scores.tsv"))

write_tsv(metadata, here::here("sample_metadata/metadata.tsv"))
write_tsv(metadata_phenotypes, here::here("sample_metadata/metadata_phenotypes.tsv"))
