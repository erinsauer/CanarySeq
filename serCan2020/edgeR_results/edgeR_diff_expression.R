#setRepositories(ind = c(1:5))
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(edgeR, tidyverse, variancePartition,
               BiocParallel, limma, pander
               )

# library(edgeR)
# library(tidyverse)
# library(variancePartition)
# library(BiocParallel)
# library(limma)
# library(pander)

# setwd("~/CanarySeq/globinModifiedGTFResults")
setwd("~/Documents/GitHub/CanarySeq/serCan2020")

g.counts.l <- read.delim("rsem.merged.gene_counts.tsv", 
                         header=TRUE)
# remove globulin genes before determing cpm values
g.counts.l <- g.counts.l %>%
  filter(str_detect(gene_id, c("LOC103824466|LOC103817749|LOC103817748|LOC115485052"), negate=TRUE))

row.names(g.counts.l) <- (g.counts.l$gene_id)
# View(g.tpm.l)
g.counts.l$ES89_Blue.049_Control.protein_0 <- NULL #remove the bad samples
g.counts.l$ES31_Blue.055_NA_0 <- NULL
g.counts.l$ES32_Blue.055_NA_14 <- NULL
g.counts.l$ES65_106_NA_0 <- NULL



meta <- data.frame(SampID = colnames(g.counts.l[,c(3:ncol(g.counts.l))]))
meta <- meta %>% separate(SampID, 
                          c("Sample.ID", "Bird.ID",
                            "MG.Diet", "Day"), "_", 
                          extra = "warn", fill = "warn", 
                          remove = FALSE) 
# View(meta)

meta <- meta %>% separate(MG.Diet, 
                          c("MG.Status", "Diet"), 
                          sep = "[.]", extra = "warn", 
                          fill = "warn", 
                          remove = FALSE) 

meta$MG.Diet.Day <- paste0(meta$MG.Diet, ".", meta$Day)
meta$MG.Day <- paste0(meta$MG.Status, ".", meta$Day)
meta$Diet.Day <- paste0(meta$Diet, ".", meta$Day)

meta <- meta %>% mutate(MG.active = if_else(Day < 1, "NO", 
                                            if_else(MG.Status == "Control",
                                                    "NO", "YES")))
meta$Exposure.MG <- paste0(meta$MG.active, ".", meta$MG.Status)
meta$Exposure.Diet <- paste0(meta$MG.active, ".", meta$Diet)
meta$Exposure.Diet.Day <- paste0(meta$Exposure.Diet, ".", meta$Day)

# View(meta)
row.names(meta) <- meta$SampID
######## diff exp ########################
y<-DGEList(counts=g.counts.l[,c(3:ncol(g.counts.l))], group=meta$Exposure.Diet)

Bird.ID <- meta$Bird.ID

keep <- filterByExpr(y)
keep <- rowSums(cpm(y)>1) >= 10
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
geneExpr <- y$counts



Bird <- factor(meta$Bird.ID)
Disease <- factor(meta$MG.active, levels=c("NO","YES"))
Diet <- factor(meta$Diet, levels=c("protein","lipid"))
Treatment <- factor(meta$MG.Status, levels=c("Control","MG"))
design <- model.matrix(~Bird)

MG.protein <- Disease== "YES" & Diet=="protein"
MG.lipid <- Disease=="YES" & Diet=="lipid"
Trt.protein <- Treatment=="MG" & Diet == "protein"
Trt.lipid <- Treatment=="MG" & Diet == "lipid"
# MG <- Disease== "YES"
protein <- Diet=="protein"
Control <- Disease=="NO"
# Control.lipid <- Disease=="NO" & Diet=="lipid"
design <- cbind(design, MG.protein, MG.lipid)#, Trt.protein)#, Trt.lipid)#, Control.protein, Control.lipid)
# design <- cbind(design, MG)


y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# fit$coefficients %>% colnames()

control_birds <- meta %>%
  filter(MG.Status=="Control") %>%
  mutate(bird=paste0("Bird",Bird.ID)) %>%
  pull(bird) %>% unique()

protein_control_birds <- meta %>%
  filter(Diet=="protein") %>%
  filter(MG.Status=="Control") %>%
  mutate(bird=paste0("Bird",Bird.ID)) %>%
  pull(bird) %>% unique()

protein_birds <- meta %>%
  filter(Diet=="protein") %>%
  mutate(bird=paste0("Bird",Bird.ID)) %>%
  pull(bird) %>% unique()


control_birds_positions <- colnames(fit$coefficients) %in% control_birds
protein_control_birds_positions <- colnames(fit$coefficients) %in% protein_control_birds

protein_birds_positions <- colnames(fit$coefficients) %in% protein_birds


diet_control_birds_contrast <- c(control_birds_positions-2*protein_control_birds_positions) %>%
  data.frame() %>%
  mutate(vector = case_when(
    . == -1 ~ -0.1,
    . == 1 ~ 1/9,
    .default = 0
  )) %>% 
  pull(vector)

diet_birds_contrast <- as.numeric(protein_birds_positions) %>%
  data.frame() %>%
  mutate(vector = case_when(
    . == 1 ~ -1/sum(protein_birds_positions),
    . == 0 ~ 1/(sum(!protein_birds_positions)-3),
    .default = 0
  )) %>% 
  # mutate(vector=sample(vector)) %>%
  pull(vector)
  
diet_birds_contrast[38:39] <- 0 #mg.prot and mg.lipid contrasts
diet_birds_contrast[1] <- 0 #yintercept

#find DE genes for protein infected
qlf.prot.inf <- glmQLFTest(fit, coef="MG.protein")
# MGP.top <- topTags(qlf, n=length(row.names(qlf$table)),
        # adjust.method="BH")

#find DE genes for lipid infected
qlf.lipid.inf <- glmQLFTest(fit, coef="MG.lipid")
# MGL.top <- topTags(qlf, n = 50) 

#look at effect of diet on uninfected birds
# qlf.diet <- glmQLFTest(fit, contrast=c(control_birds_positions-2*protein_birds_positions))
qlf.diet.ctl <- glmQLFTest(fit, contrast=diet_control_birds_contrast)

#diet effect on all birds
qlf.diet <- glmQLFTest(fit, contrast=diet_birds_contrast)


# find DE genes for all infected birds
# qlf.inf <- glmQLFTest(fit, coef=c("MG.protein","MG.lipid"))
qlf.inf <- glmQLFTest(fit, contrast=c(rep(0,37),0.5,0.5))
# topTags(qlf, n = 50) 
# MG.top <- topTags(qlf, n=length(row.names(qlf$table)),
                  # adjust.method="BH", p.value=0.05)

# qlf <- glmQLFTest(fit, coef=c("MG"))
# topTags(qlf, n = 50) 
# MG.top <- topTags(qlf, n=length(row.names(qlf$table)),
                  # adjust.method="BH", p.value=0.05)

# MG.top$table %>%
  # filter(FDR < 0.05) %>%
  # arrange(desc(abs(logFC))) %>% head()

#compare difference in response to infection between diets (lipid-protein)
qlf.diet.inf.diff <- glmQLFTest(fit, contrast=c(rep(0,37),1,-1))
# MG.diet.top <- topTags(qlf) # no significant difference identified


# visualize MDS plot corrected for Bird batch effects
logCPM <- cpm(y, log=TRUE)
design.Status <- model.matrix(~Disease)
logCPM.corrected <- removeBatchEffect(logCPM, batch=Bird, design=design.Status)
plotMDS(logCPM.corrected, label=Disease, top = 1000)


### INSERTED CODE ###
getTopTags <- function(x, contrast) {
  x %>%
    topTags(n=Inf) %>%
    data.frame() %>%
    rownames_to_column("genes") %>%
    mutate(comparison=as.character(contrast))
}


qlf_each <- bind_rows(
  getTopTags(qlf.lipid.inf, "lipid.infected"),
  getTopTags(qlf.prot.inf, "protein.infected"),
  getTopTags(qlf.diet.inf.diff, "diet.inf.difference"),
  getTopTags(qlf.inf, "infected"),
  getTopTags(qlf.diet.ctl, "diet_control"),
  getTopTags(qlf.diet, "diet")
)

qlf_each %>%
  # pivot_longer(-c(genes,logCPM, F, PValue, FDR, comparison), names_prefix = "logFC.", names_to = "contrast",
               # values_to = "logFC") %>%
  mutate(`FDR<0.05` = case_when(
    FDR < 0.05 ~ "Significant (by contrast)",
    .default = "Not Significant"
  )) %>%
  filter(FDR<=0.1) %>%
  ggplot(data = ., aes(x=logCPM, y = logFC, color = `FDR<0.05`, shape=`FDR<0.05`)) + 
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.5, size =0.8) +
  # ylim(-10,10) +
  facet_wrap(~comparison) +
  # facet_grid(rows=vars(strain), cols=vars(stress)) +
  ggtitle("QLF: Effect of diet and infection on stress response") +
  theme_classic()




effectSize <- 1.5

qlf_each %>%
  mutate(`FDR<0.05` = case_when(
    FDR < 0.05 ~ "Significant (by contrast)",
    .default = "Not Significant"
  ),
  # HigherIn = case_when(
  #   logFC <= -effectSize ~ "WT",
  #   logFC >= effectSize ~ "HOG1Δ",
  #   .default = "Neither"
  # ),
  delabel = case_when(
    FDR < 0.05 & abs(logFC)>0.5 ~ genes,
    .default = NA
  )
  ) %>% 
  # filter(ESR != "neither") %>%
  # filter(ESR == "iESR") %>%
  # filter(ESR == "rESR") %>%
  ggplot(data = ., aes(x=logFC, y = -log10(PValue), color = `FDR<0.05`,
                       label=delabel)) + 
  geom_hline(yintercept = 0) +
  geom_point(alpha=0.5, size =0.8) +
  # ylim(-10,10) +
  facet_wrap(~comparison) +
  # ggrepel::geom_text_repel(size=1.5, color = "black") +
  # ggtitle("QLF: Effect of HOG1Δ on stress response") +
  # labs(caption = "postive = higher in HOG1Δ, negative = higher in WT") +
  theme_classic() +
  theme(legend.position="bottom")




# look at raw expression real quick
logCPM <- edgeR::cpm(y, prior.count=1, log=TRUE)


tidyFC <- logCPM %>%
  as_tibble(rownames = "gene") %>%
  pivot_longer(., cols = -gene,
               names_to = "sample",
               values_to = "logCPM") %>% 
  separate(sample, 
           c("Sample.ID", "Bird.ID",
             "MG.Diet", "Day"), "_", 
           extra = "warn", fill = "warn", 
           remove = FALSE) %>%
  separate(MG.Diet, 
               c("MG.Status", "Diet"), 
               sep = "[.]", extra = "warn", 
               fill = "warn", 
               remove = FALSE) %>%
  mutate(MG.Diet.Day = paste0(MG.Diet, ".", Day),
         MG.Day = paste0(MG.Status, ".", Day),
         Diet.Day <- paste0(Diet, ".", Day)
  ) %>%
  mutate(MG.active = if_else(Day < 1, "NO", 
                             if_else(MG.Status == "Control",
                                     "NO", "YES")))
  

# meta$MG.Diet.Day <- paste0(meta$MG.Diet, ".", meta$Day)
# meta$MG.Day <- paste0(meta$MG.Status, ".", meta$Day)
# meta$Diet.Day <- paste0(meta$Diet, ".", meta$Day)
# 
# meta <- meta %>% mutate(MG.active = if_else(Day < 1, "NO", 
#                                             if_else(MG.Status == "Control",
#                                                     "NO", "YES")))
# meta$Exposure.MG <- paste0(meta$MG.active, ".", meta$MG.Status)
# meta$Exposure.Diet <- paste0(meta$MG.active, ".", meta$Diet)
# meta$Exposure.Diet.Day <- paste0(meta$Exposure.Diet, ".", meta$Day)

gene_of_interest = "LOC103820730" #T-cell receptor-associated transmembrane adapter 1
gene_of_interest = "LOC103824071" #BCAN weirdly has bimodal distribution of expression.
gene_of_interest = "TENT5C"
gene_of_interest = "LOC103820912" #VWA1 gene
gene_of_interest = "LOC103820351" #RSAD2 gene
gene_of_interest = "LOC103820477" #BPGM
gene_of_interest = "LOC103812956" #POMGNT2
gene_of_interest = "LOC103816332" #LNX2 like, related to ubiquination in humans
gene_of_interest = "LOC103815956" #tln2
gene_of_interest = "LOC127060867" # ankyrin repeat and SOCS box containing 16
gene_of_interest = "LOC103822412" # RAB44, smallest FDR for positive logFC in diet_ctl contrast

# the contrast that makes the most sense is just all the control birds with diet.

tidyFC %>%
  filter(gene == gene_of_interest) %>%
  # filter(genotype == "WT") %>%
  # ggplot(aes(x=Bird, y = logCPM, shape=MG.Status, color=Diet, alpha=Day)) +
  ggplot(aes(x=Bird, y = logCPM, shape=MG.Status, color=as.numeric(Day))) +
  # ggbeeswarm::geom_beeswarm(cex=2) +
  geom_point() +
  # geom_boxplot()+
  # facet_grid(rows = vars(stress), cols = vars(genotype)) +
  # facet_wrap(~factor(genotype, levels=c("WT", "HOG1"))) +
  theme_classic() +
  ggtitle(paste0(gene_of_interest," expression")) +
  coord_flip() +
  facet_grid(rows=vars(MG.Status), cols=vars(Diet), scales = 'free_y')
  # ggtitle("expression of MOT3 in response to NaCl stress") +
  # ggpubr::stat_compare_means(label = "p.format", method = "t.test", 
  #                    paired=T, tip.length = 0, label.y = 4.5) 



median_CPM <- tidyFC %>%
  filter(gene == gene_of_interest) %>%
  group_by(Diet, MG.active) %>%
  summarize(mean = mean(logCPM),
            median=median(logCPM))
# density plot
tidyFC %>%
  filter(gene == gene_of_interest) %>%
  ggplot(aes(x=logCPM, fill=Diet)) +
  # geom_histogram(alpha=0.5,binwidth = 0.1, position = "identity") +
  geom_density(alpha=0.1,aes(y=after_stat(count), color=Diet)) +
  geom_dotplot(binwidth=0.1, dotsize = 1,method='histodot', alpha=0.7,stackgroups=TRUE) +
  # geom_density(alpha=0.5) +
  # geom_histogram(alpha=0.5,position = "identity") +
  geom_vline(data=median_CPM, aes(xintercept=mean, color=Diet))+
  facet_wrap(~MG.active, nrow=2) +
  ggtitle(paste0(gene_of_interest," expression")) +
  theme_classic()


tidyFC %>%
  filter(gene == gene_of_interest) %>%
  # filter(genotype == "WT") %>%
  # ggplot(aes(x=Bird, y = logCPM, shape=MG.Status, color=Diet, alpha=Day)) +
  ggplot(aes(x=Bird, y = logCPM, shape=MG.Status, color=as.numeric(Day))) +
  # ggbeeswarm::geom_beeswarm(cex=2) +
  geom_point() +
  # geom_boxplot()+
  # facet_grid(rows = vars(stress), cols = vars(genotype)) +
  # facet_wrap(~factor(genotype, levels=c("WT", "HOG1"))) +
  theme_classic() +
  ggtitle(paste0(gene_of_interest," expression")) +
  coord_flip() +
  facet_grid(rows=vars(Diet), cols=vars(MG.Status), scales = 'free_y')

### END INSERTED CODE ###










# volcano plot code
increased <- MG.top$table$logFC > 0 & MG.top$table$FDR < 0.05
plot(MG.top$table$logFC, -10*log10(MG.top$table$PValue), 
     main="Volcano plot", xlab="M", ylab="-10*log(P-val)",type = "p",
     pch=19
     )

# highlight our DE genes
points(MG.top$table$logFC[increased], 
       -10*log10(MG.top$table$PValue[increased]), col="red", pch=19)

# identify genes enriched in the control
decreased <- MG.top$table$logFC < 0 & MG.top$table$FDR < 0.05
points(MG.top$table$logFC[decreased], 
       -10*log10(MG.top$table$PValue[decreased]), col="blue", pch=19)

# extracts the log-fold change and p-values 
# (can be FDR-corrected), ranked by p-value
View(data.frame(d.MGLP1421))














g.tpm.l <- read.delim("rsem.merged.gene_tpm.tsv", 
                                          header=TRUE)
row.names(g.tpm.l) <- (g.tpm.l$gene_id)
View(g.tpm.l)
g.tpm.l$ES89_Blue.049_Control.protein_0 <- NULL #remove the bad samples
g.tpm.l$ES31_Blue.055_NA_0 <- NULL
g.tpm.l$ES32_Blue.055_NA_14 <- NULL
g.tpm.l$ES65_106_NA_0 <- NULL


############ make the metadata table ######################
meta <- data.frame(SampID = colnames(g.tpm.l[,c(3:ncol(g.tpm.l))]))
meta <- meta %>% separate(SampID, 
                          c("Sample.ID", "Bird.ID",
                            "MG.Diet", "Day"), "_", 
                          extra = "warn", fill = "warn", 
                          remove = FALSE) 
# View(meta)

meta <- meta %>% separate(MG.Diet, 
                          c("MG.Status", "Diet"), 
                          sep = "[.]", extra = "warn", 
                          fill = "warn", 
                          remove = FALSE) 

meta$MG.Diet.Day <- paste0(meta$MG.Diet, ".", meta$Day)
meta$MG.Day <- paste0(meta$MG.Status, ".", meta$Day)
meta$Diet.Day <- paste0(meta$Diet, ".", meta$Day)

meta <- meta %>% mutate(MG.active = if_else(Day < 1, "NO", 
                                            if_else(MG.Status == "Control",
                                                    "NO", "YES")))
meta$Exposure.MG <- paste0(meta$MG.active, ".", meta$MG.Status)
meta$Exposure.Diet <- paste0(meta$MG.active, ".", meta$Diet)
meta$Exposure.Diet.Day <- paste0(meta$Exposure.Diet, ".", meta$Day)

# View(meta)
row.names(meta) <- (meta$SampID)
############ diff exp #################################
My.dge<-DGEList(counts=g.tpm.l[,c(3:ncol(g.tpm.l))], group=meta$Exposure.Diet.Day)# was group MG.Diet.Day
summary(My.dge)
#barplot(My.dge$samples$lib.size*1e-6, names=1:ncol(My.dge), ylab="Library size (millions)")

## creates object from table of counts, 
## with a data frame counting information for each sample
View(My.dge$samples)
#keep <- filterByExpr(My.dge)
keep <- rowSums(cpm(My.dge)>1) >= 6
My.dge <- My.dge[keep, , keep.lib.sizes=FALSE]
My.dge <- calcNormFactors(My.dge)
geneExpr <- My.dge$counts
## uses TMM (weighted trimmed mean of M-values) 
## to calculate normalization factors (i.e. effective library sizes) 
## for each sample.

## rule of thumb, require that genes have at least 10-15 reads 
## in at least one condition to be considered expressed
## removing unexpressed genes increases our statistical power 
## to detect differential expression (because fewer tests)

My.dge$samples

plotMD(My.dge, column=1)
abline(h=0, col="red", lty=2, lwd=2)

plotMDS(My.dge,col=as.numeric(My.dge$samples$group), gene.selection="common",
        labels = NULL,pch = "*"
        )

Plot.dge<-DGEList(counts=g.tpm.l[,c(3:106)], group=meta$Bird.ID)
plotMDS(Plot.dge,col=as.numeric(Plot.dge$samples$group, gene.selection="pairwise"))
# legend("bottomleft", as.character(unique(My.dge$samples$group)),
#        col=(My.dge$samples$group), pch=20)

## mean-difference (MD) plot: looks at fold-changes in 
## one library vs. others (should cluster near zero)

design <- model.matrix(~0+group, data=My.dge$samples)  

# setting up our design matrix; 
# ~0 = do not include intercept (i.e. common reference)

My.yglm <- estimateDisp(My.dge, design)
My.yglm$common.dispersion
## uses negative binomial likelihood to estimate common, 
## trended, and tagwise dispersions 

plotBCV(My.yglm)

## plots dispersions

My.fit <- glmQLFit(My.yglm, design)
plotQLDisp(My.fit)
# Fit a quasi-likelihood negative binomial 
# generalized log-linear model to count data

My.contrasts <- makeContrasts(
  Diet0 = ((groupMG.protein.0 + groupControl.protein.0)/2) -
          ((groupMG.lipid.0 + groupControl.lipid.0)/2),
  MGpost0=((groupMG.protein.14 + groupMG.protein.21 +
            groupMG.lipid.14 + groupMG.lipid.21)/4) - 
          ((groupControl.protein.14 + groupControl.protein.21 +
            groupControl.lipid.14 + groupControl.lipid.21)/4),
  MGd140= ((((groupMG.lipid.14 + groupMG.protein.14)/2) -
          ((groupMG.lipid.0 + groupMG.protein.0)/2))-
            (((groupControl.lipid.14 + groupControl.protein.14)/2) -
               ((groupControl.lipid.0 + groupControl.protein.0)/2))),
  MGd210= ((((groupMG.lipid.21 + groupMG.protein.21)/2) -
              ((groupMG.lipid.0 + groupMG.protein.0)/2))-
             (((groupControl.lipid.21 + groupControl.protein.21)/2) -
                ((groupControl.lipid.0 + groupControl.protein.0)/2))),
  MGd2114=((((groupMG.lipid.21 + groupMG.protein.21)/2) -
              ((groupMG.lipid.14 + groupMG.protein.14)/2))-
             (((groupControl.lipid.21+ groupControl.protein.21)/2) -
                ((groupControl.lipid.14 + groupControl.protein.14)/2))),
  MG14=   (((groupMG.lipid.14 + groupMG.protein.14)/2) -
          ((groupControl.lipid.14 + groupControl.protein.14)/2)),
  MG21=   (((groupMG.lipid.21 + groupMG.protein.21)/2) -
             ((groupControl.lipid.21 + groupControl.protein.21)/2)),

  MGLpost0=((groupMG.lipid.14 + groupMG.lipid.21)/2) -
           ((groupControl.lipid.14 + groupControl.lipid.21)/2),
  MGLd140= ((groupMG.lipid.14 - groupMG.lipid.0) -
           (groupControl.lipid.14 - groupControl.lipid.0)),
  MGLd210= ((groupMG.lipid.21 - groupMG.lipid.0) -
           (groupControl.lipid.21 - groupControl.lipid.0)),
  MGLd2114= ((groupMG.lipid.21 - groupMG.lipid.14) -
            (groupControl.lipid.21 - groupControl.lipid.14)),
  MGL14=   (groupMG.lipid.14 - groupControl.lipid.14),
  MGL21=   (groupMG.lipid.21 - groupControl.lipid.21),

  MGPpost0=((groupMG.protein.14 + groupMG.protein.21)/2) -
           ((groupControl.protein.14 + groupControl.protein.21)/2),
  MGPd140= ((groupMG.protein.14 - groupMG.protein.0) -
              (groupControl.protein.14 - groupControl.protein.0)),
  MGPd210= ((groupMG.protein.21 - groupMG.protein.0) -
              (groupControl.protein.21 - groupControl.protein.0)),
  MGPd2114= ((groupMG.protein.21 - groupMG.protein.14) -
               (groupControl.protein.21 - groupControl.protein.14)),
  MGP14=   (groupMG.protein.14 - groupControl.protein.14),
  MGP21=   (groupMG.protein.21 - groupControl.protein.21),
  
  DietIfpost0=((groupMG.lipid.14 + groupMG.lipid.21)/2) - 
              ((groupMG.protein.14 + groupMG.protein.21)/2),
  DietIfd140= ((groupMG.lipid.14 - groupMG.lipid.0) - 
              (groupMG.protein.14 - groupMG.protein.0)),
  DietIfd210= ((groupMG.lipid.21 - groupMG.lipid.0) - 
              (groupMG.protein.21 - groupMG.protein.0)),
  DietIfd2114=((groupMG.lipid.21 - groupMG.lipid.14) - 
               (groupMG.protein.21 - groupMG.protein.14)),
  DietIf14=   (groupMG.lipid.14 - groupMG.protein.14),
  DietIf21=   (groupMG.lipid.21 - groupMG.protein.21),
  
  MGDietpost0=((((groupMG.lipid.14 + groupMG.lipid.21)/2) - 
                 ((groupControl.lipid.14 + groupControl.lipid.21)/2))-
                 (((groupMG.protein.14 + groupMG.protein.21)/2) - 
              ((groupControl.protein.14 + groupControl.protein.21)/2))),
  MGDietd140= (((groupMG.lipid.14 - groupMG.lipid.0) - 
                  (groupControl.lipid.14 - groupControl.lipid.0))-
                 ((groupMG.protein.14 - groupMG.protein.0) - 
              (groupControl.protein.14 - groupControl.protein.0))),
  MGDietd210= (((groupMG.lipid.21 - groupMG.lipid.0) - 
                  (groupControl.lipid.21 - groupControl.lipid.0))-
                 ((groupMG.protein.21 - groupMG.protein.0) - 
              (groupControl.protein.21 - groupControl.protein.0))),
  MGDietd2114= (((groupMG.lipid.21 - groupMG.lipid.14) - 
                   (groupControl.lipid.21 - groupControl.lipid.14))-
                  ((groupMG.protein.21 - groupMG.protein.14) - 
               (groupControl.protein.21 - groupControl.protein.14))),
  MGDiet14=   ((groupMG.lipid.14 - groupControl.lipid.14)-
                 (groupMG.protein.14 - groupControl.protein.14)),
  MGDiet21=   ((groupMG.lipid.21 - groupControl.lipid.21)-
                 (groupMG.protein.21 - groupControl.protein.21)),
 levels=design)

# 
names <- c("Diet0", 
        "MGpost0", "MGd140", "MGd210", "MGd2114", "MG14", "MG21", 
        "MGLpost0", "MGLd140", "MGLd210", "MGLd2114", "MGL14", "MGL21", 
        "MGPpost0", "MGPd140", "MGPd210", "MGPd2114", "MGP14", "MGP21",
        "DietIfpost0", "DietIfd140", "DietIfd210", "DietIfd2114", "DietIf14", "DietIf21",
        "MGDietpost0", "MGDietd140", "MGDietd210", "MGDietd2114", "MGDiet14", "MGDiet21")

setwd("~/Documents/GitHub/CanarySeq/serCan2020/edgeR_results")

for (i in 1:31){
  contrast <- My.contrasts[,i] 
  model <- glmQLFTest(My.fit, contrast=contrast)
  # Conduct genewise statistical tests for the given contrast
  TT <- topTags(model, n=length(row.names(model$table)),
                       adjust.method="BH")
  #save table and MD plot
  write.csv(TT, paste0(names[i],".csv"))
  jpeg(paste0(names[i],".jpeg"))
  plotMD(model)
  abline(h=c(-1, 1), col="blue")
  dev.off()
}

# volcano plot code
increased <- d.Diet0$table$logFC > 0 & d.Diet0$table$FDR < 0.01
plot(d.Diet0$table$logFC, -10*log10(d.Diet0$table$PValue), 
     main="Volcano plot", xlab="M", ylab="-10*log(P-val)")

# highlight our DE genes
points(d.Diet0$table$logFC[increased], 
       -10*log10(d.Diet0$table$PValue[increased]), col="red")

# identify genes enriched in the control
decreased <- d.Diet0$table$logFC < 0 & d.Diet0$table$FDR < 0.01
points(d.Diet0$table$logFC[decreased], 
       -10*log10(d.Diet0$table$PValue[decreased]), col="blue")

# extracts the log-fold change and p-values 
# (can be FDR-corrected), ranked by p-value
View(data.frame(d.MGLP1421))
write.csv(MGLP1421, "d.MGLP1421.csv")


gene.ann <- read.delim("rsem.merged.gene_tpm_labels.tsv", 
                      header=TRUE)
View(gene.ann)
