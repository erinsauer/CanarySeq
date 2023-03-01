#setRepositories(ind = c(1:5))
if (!require("pacman")) install.packages("pacman")
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

g.tpm.l <- read.delim("rsem.merged.gene_tpm.tsv", 
                                          header=TRUE)
row.names(g.tpm.l) <- (g.tpm.l$gene_id)
View(g.tpm.l)
g.tpm.l$ES89_Blue.049_Control.protein_0 <- NULL #remove the bad samples


############ make the metadata table ######################
meta <- data.frame(SampID = colnames(g.tpm.l[,c(3:106)]))
meta <- meta %>% separate(SampID, 
                          c("Sample.ID", "Bird.ID",
                            "MG.Diet", "Day"), "_", 
                          extra = "warn", fill = "warn", 
                          remove = FALSE) 
View(meta)

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

View(meta)
row.names(meta) <- (meta$SampID)
############ diff exp #################################
My.dge<-DGEList(counts=g.tpm.l[,c(3:106)], group=meta$MG.Diet.Day)
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

plotMDS(My.dge,col=as.numeric(My.dge$samples$group), gene.selection="common")

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
