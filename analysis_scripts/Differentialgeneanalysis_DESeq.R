# This script was used to perform pairwise differential gene expression analysis between RNAseq examples. 

#####SessionInfo#####

# > sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Catalina 10.15.7
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] GGally_2.1.2                ggplot2_3.3.6               dplyr_1.0.9                 stringr_1.4.0               edgeR_3.36.0               
# [6] DESeq2_1.34.0               SummarizedExperiment_1.24.0 Biobase_2.54.0              MatrixGenerics_1.6.0        matrixStats_0.62.0         
# [11] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4            BiocGenerics_0.40.0        
# [16] limma_3.50.3               
# 
# loaded via a namespace (and not attached):
#   [1] locfit_1.5-9.5         Rcpp_1.0.8.3           lattice_0.20-45        png_0.1-7              Biostrings_2.62.0      assertthat_0.2.1      
# [7] utf8_1.2.2             plyr_1.8.7             R6_2.5.1               RSQLite_2.2.14         httr_1.4.3             pillar_1.7.0          
# [13] zlibbioc_1.40.0        rlang_1.0.3            rstudioapi_0.13        annotate_1.72.0        blob_1.2.3             Matrix_1.4-1          
# [19] splines_4.1.2          BiocParallel_1.28.3    geneplotter_1.72.0     RCurl_1.98-1.7         bit_4.0.4              munsell_0.5.0         
# [25] DelayedArray_0.20.0    compiler_4.1.2         pkgconfig_2.0.3        tidyselect_1.1.2       KEGGREST_1.34.0        tibble_3.1.7          
# [31] GenomeInfoDbData_1.2.7 XML_3.99-0.10          reshape_0.8.9          fansi_1.0.3            withr_2.5.0            crayon_1.5.1          
# [37] bitops_1.0-7           grid_4.1.2             xtable_1.8-4           gtable_0.3.0           lifecycle_1.0.1        DBI_1.1.3             
# [43] magrittr_2.0.3         scales_1.2.0           stringi_1.7.6          cli_3.3.0              cachem_1.0.6           XVector_0.34.0        
# [49] genefilter_1.76.0      ellipsis_0.3.2         generics_0.1.3         vctrs_0.4.1            RColorBrewer_1.1-3     tools_4.1.2           
# [55] bit64_4.0.5            glue_1.6.2             purrr_0.3.4            parallel_4.1.2         fastmap_1.1.0          survival_3.3-1        
# [61] AnnotationDbi_1.56.2   colorspace_2.0-3       memoise_2.0.1 
#####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install("DESeq2")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# BiocManager::install("stringr")
# BiocManager::install("dplyr")
# BiocManager::install("GGally")

library('limma')
library("DESeq2")
library('edgeR')
library('stringr')
library('dplyr')
library('ggplot2')
library('GGally')

convertensemblgenetosymbol <- function(table, tablecolumn) {
  #detach("package:tidyverse", unload=TRUE)
  library(dplyr)
  detach("package:dplyr", unload=TRUE)
  library(DBI)
  library(org.Mm.eg.db)
  x <- org.Mm.egENSEMBL
  genesofinterest <- as.data.frame(table)
  ensembllist <- as.character(tablecolumn)
  #ensembllist <- c(rownames(genesofinterest))
  
  #select(org.Mm.eg.db, keys = ensembllist, columns =c("SYMBOL", "GOALL"), keytype = "ENSEMBL" )
  ensembllistanno <- mapIds(org.Mm.eg.db, keys = ensembllist, column = c("SYMBOL"), keytype = "ENSEMBL")
  ensembllistannodf <- as.data.frame(ensembllistanno)
  
  ensembllist <- cbind(ensembllist, ensembllistannodf)
  genesofinterestanno <- cbind(genesofinterest, ensembllistannodf)
  output <- as.data.frame(genesofinterestanno)
  #write.csv(genesofinterestanno, output)
  detach("package:org.Mm.eg.db", unload=TRUE)
  library(dplyr)
  output <- as.data.frame(genesofinterestanno)
}

# read data
counts <- read.table("~/Desktop/exoncounttables/exonCounts-RBRB01_04-2020-07-07.txt", header=T, sep="\t", row.names = 1, check.name=FALSE,
                                     stringsAsFactor=FALSE)
#counts <- counts[6:29]

##Rearrange columns for proper DGE comparison 
construct <- c("FCKO","FCKO","FCKO", "FdCD","FdCD", "FdCD", "FDKO","FDKO","FDKO", "FWT1", "FWT1", "FWT1", "NCKO","NCKO","NCKO", "NdCD","NdCD","NdCD", "NDKO","NDKO","NDKO", "NWT","NWT","NWT" )
replicate <- c("1", "2", "3")

annot <- cbind(construct, replicate)
annot <- as.data.frame(annot)
rownames(annot) <- c("FCKO1","FCKO2","FCKO3", "FdCD1","FdCD2", "FdCD3", "FDKO1","FDKO2","FDKO3", "FWT1", "FWT2", "FWT3", "NCKO1","NCKO2","NCKO3", "NdCD1","NdCD2","NdCD3", "NDKO1","NDKO2","NDKO3", "NWT1","NWT2","NWT3" )
#rownames(annot) <- c("WT1", "WT2", "DKO1", "DKO2", "dCD1", "dCD2")
#colnames(df_lcpm) <- rownames(annot)
#read in sample names from samplelist file
#samples <- read.table("/Users/Ryan/Desktop/samplelist_RB.txt", header=F, sep="\t", row.names = 1, check.name=FALSE,
#stringsAsFactor=FALSE)
#samples.df <- as.data.frame(samples)

#columns of all possible comparisons
FCKO <- counts[, 6:8]
FdCD <- counts[, 9:11]
FDKO <- counts[, 12:14]
FWT <- counts[, 15:17]
NCKO <- counts[, 18:20]
NdCD <- counts[, 21:23]
NDKO <- counts[, 24:26]
NWT <- counts[, 27:29]
# 
counts <- data.frame(FCKO, FDKO)
output <- "~/Desktop/Diffgenes/DiffgenesNWTvsNDKOfull-07.07.20.anno.csv"
output <- "~/Desktop/data/Diffgenes/DiffgenesFCKOvsFDKO-11.20.22.anno.csv"

# colnames(counts) <- c("A1", "A2", "A3", "B1","B2", "B3")
# #colnames(counts) <- c("WT1", "WT2", "DKO1", "DKO2", "dCD1", "dCD2")
#
construct <- c("A", "A", "A", "B", "B", "B")
replicate <- c("1", "2", "3")
annot <- cbind(construct, replicate)
annot <- as.data.frame(annot)
# rownames(annot) <- c("A1","A2","A3","B1","B2","B3")

colnames(counts) <- c(rownames(annot))

y <- DGEList(counts=counts,samples=annot)
samplenames <- colnames(y) 
samplenames

construct <- as.factor(y$samples$construct)
#time <- as.factor(y$samples$time)
replicate <- as.factor(y$samples$replicate)
y$geneid <- rownames(y)
lib.size <- as.factor(y$samples$lib.size)

# Subset samples:
x <- y[,,keep.lib.sizes=TRUE] ##subset by what?
dim(x)

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

table(rowSums(x$counts==0)==24)

keep.exprs <- rowSums(cpm>1)>=2 ### changed from at least 1 cpm in 2 samples to
x <- x[keep.exprs,, keep.lib.sizes=FALSE] # (=F) recounts lib size after trim

dim(x)

## Plot Log-CPM RAW and Filtered data
library(RColorBrewer)

nsamples <- ncol(x)
# col <- brewer.pal(nsamples, "Paired") #original
col <- brewer.pal(8, "Paired") # find pallete to accommodate >12 samples if necessary
par(mfrow=c(1,2))

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)

plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

lcpm_x <- cpm(x, log=TRUE)

x_norm <- calcNormFactors(x, method = "TMM")
#x_norm <- calcNormFactors(x, method = "upperquartile")
x_norm$samples$norm.factors

lcpm_x_norm <- cpm(x_norm, log=TRUE)
x_norm <- x
x_norm$samples$norm.factors <- 1
#x_norm$counts[,1] <- ceiling(x_norm$counts[,1]*0.05)
#x_norm$counts[,2] <- x_norm$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(x_norm, log=TRUE)

## Boxplots of Data pre and post-normalization
boxplot(lcpm_x, las=2, col=col, main="")
title(main="Unnormalised data",ylab="Log-cpm")
x_norm <- calcNormFactors(x_norm)
x_norm$samples$norm.factors

lcpm <- cpm(x, log=TRUE)
boxplot(lcpm_x_norm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- x$samples$group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.phase <- x$samples$phase
levels(col.phase) <- brewer.pal(nlevels(col.phase), "Set2")
col.phase <- as.character(col.phase)
col.lib_prep <- x$samples$lib_prep
levels(col.lib_prep) <- brewer.pal(nlevels(col.lib_prep), "Set3")
col.lib_prep <- as.character(col.lib_prep)



## Plot MDS
# dim 1,2 by time
par(mfrow = c(1,2))
plotMDS(lcpm, labels=x_norm$samples$time, pch=20, col=col)
title(main="Sample time")
par(mfrow = c(1,2))
legend("topright", samplenames, text.col=col.group, bty="n", ncol=1)

# # dim 1,2 by construct
par(mfrow = c(1,2))
plotMDS(lcpm, labels=x_norm$samples$construct, pch=20, col=col)
title(main="Sample Number")
par(mfrow = c(1,2))
#legend("topright", samplenames, text.col=col, bty="n", ncol=1)

# dim 1,2 by replicate
par(mfrow = c(1,2))
plotMDS(lcpm, labels=x_norm$samples$replicate, pch=20, col=col)
title(main="replicate")

par(mfrow = c(1,1))
plotMDS(lcpm, labels=x_norm$samples$litter, pch=20, col=col)
title(main="replicate")

par(mfrow = c(1,1))
plotMDS(lcpm, labels=x$samples$litter, col=col, dim=c(1,2))
title(main="MDS plot of MLL mutants")
legend("topright", samplenames, text.col=col, bty="n", ncol=1)

df_lcpm <- as.data.frame(lcpm)
df_lcpmannot <- convertensemblgenetosymbol(df_lcpm, rownames(df_lcpm))

library('DESeq2')
# counts matrix: genes as rows and samples as columns; annotation samples as rows and info in columns
DESeqOurcts <- as.matrix(counts)
DESeqOurcoldata <- as.matrix(annot)

all(rownames(DESeqOurcoldata) %in% colnames(DESeqOurcts))
all(rownames(DESeqOurcoldata) == colnames(DESeqOurcts))

##This is where the DESeq2 vignette starts

# ensure 'control level' is the 1st level of annotated input
DESeqOurdds <- DESeqDataSetFromMatrix(countData = DESeqOurcts,
                                      colData = DESeqOurcoldata,
                                      design = ~ construct)
DESeqOurdds

# add additional metadata
DESeqOurfeatureData <- data.frame(gene=rownames(DESeqOurcts))
mcols(DESeqOurdds) <- DataFrame(mcols(DESeqOurdds), DESeqOurfeatureData)
mcols(DESeqOurdds)

# prefiltering
DESeqOurdds <- DESeqOurdds[ rowSums(counts(DESeqOurdds)) > 1, ]

# explicitly set reference (default is 1st by alphabetic); not necessary if using contrasts
##DESeqOurdds$condition <- relevel(DESeqOurdds$group, ref="e75") ##obvs 7.5 is not in this but what would be good ref?

### Differential Expression
OurddsMF <- DESeqOurdds

design(OurddsMF) <- formula(~ construct)
OurddsMF <- DESeq(OurddsMF)
OurresMF <- results(OurddsMF)
head(OurresMF)




 # MA plot
plotMA(OurresMF, alpha = 0.05, ylim=c(-8,8))
 
resMFOrdered <- OurresMF[order(OurresMF$padj),] # reordered by padj
summary(resMFOrdered)
sum(OurresMF$padj < 0.05, na.rm=TRUE) # count sig genes with padj < 0.05

#filter reads for padj and log2FC and generate .csv file named by output above

resMFOrderedfiltered <-  subset(resMFOrdered, padj < 0.05 & (log2FoldChange < -1 | log2FoldChange > 1))
resMFOrderedfiltered <-  subset(resMFOrdered, padj < 0.05)

write.csv(resMFOrderedfiltered, output)
sum(OurresMF$padj < 0.05, na.rm=TRUE) # count sig genes with padj < 0.05
converted <- convertensemblgenetosymbol(resMFOrderedfiltered, rownames(resMFOrderedfiltered))
write.csv(converted, output)

#plot logcpm x logcpm of samples starting with lcpm + normalized dataframe

#define sample groups 
# FCKO <- df_lcpm[,1:3]
# FdCD <- df_lcpm[,4:6]
# FDKO <- df_lcpm[,6:9]
# FWT <- df_lcpm[,10:12]
# NCKO <- df_lcpm[,13:15]
# NdCD <- df_lcpm[,16:18]
# NDKO <- df_lcpm[,19:21]
# NWT <-df_lcpm[,22:23]
# 
# x <- cbind(FCKO, FdCD,FDKO,FWT,NCKO,NdCD,NDKO,NWT)
# 
# FWT$ave <- rowMeans(FWT, dims = 1)
# FCKO$ave <- rowMeans(FCKO, dims = 1)
# FDKO$ave <- rowMeans(FDKO, dims = 1)
# FdCD$ave <- rowMeans(FdCD, dims = 1)
# NWT$ave <- rowMeans(NWT, dims = 1)
# NCKO$ave <- rowMeans(NCKO, dims = 1)
# NDKO$ave <- rowMeans(NDKO, dims = 1)
# NdCD$ave <- rowMeans(NdCD, dims = 1)
# 
# averagedf <- cbind(FCKO$ave, FdCD$ave,FDKO$ave,FWT$ave,NCKO$ave,NdCD$ave,NDKO$ave,NWT$ave)
# colnames(averagedf) <- c("FCKO", "FdCD","FDKO","FWT","NCKO","NdCD","NDKO","NWT")
# rownames(averagedf) <- c(row.names(FWT))
# averagedf <- as.data.frame(averagedf) 
# averagedfanno <- convertensemblgenetosymbol(averagedf)
# 
# 
# averagedf2 <- as.data.frame(averagedf)
# averagedf2$ensembl_gene_id <- row.names(averagedf)
#  
# library(gridExtra)
# 
# averagedf2_k7 <- merge(k7,averagedf2, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')
# 
# p1 <- ggplot(annotlcpm, aes(x=NWT, y=NCKO)) +
#   geom_point(size=1, shape=23) +
#   geom_point(data=averagedf2_k1, aes(x=NWT, y=NCKO), shape = 19, color ='red', size =1) +
#   theme_classic()
#   #geom_text(averagedf2_k1, aes(y=0, stat=nrow(averagedf2_k1), vjust= 0))
# 
# averagedf2_k2 <- merge(k2,averagedf2, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')
# p2 <- ggplot(averagedf, aes(x=NWT, y=NDKO)) +
#   geom_point(size=1, shape=23) +
#   geom_point(data=averagedf2_k2, aes(x=NWT, y=NDKO), shape = 19, color ='red', size =1) +
#   theme_classic()
# 
# averagedf2_k3 <- merge(k3,averagedf2, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')
# p3 <- ggplot(averagedf, aes(x=NWT, y=NdCD)) +
#   geom_point(size=1, shape=23) +
#   geom_point(data=averagedf2_k3, aes(x=NWT, y=NdCD), shape = 19, color ='red', size =1) +
#   theme_classic()
# 
# averagedf2_k4 <- merge(k4,averagedf2, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')
# p4 <- ggplot(averagedf, aes(x=FWT, y=FCKO)) +
#   geom_point(size=1, shape=23) +
#   geom_point(data=averagedf2_k4, aes(x=FWT, y=FCKO), shape = 19, color ='red', size =1) +
#   theme_classic()
# 
# averagedf2_k5 <- merge(k5,averagedf2, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')
# p5 <- ggplot(averagedf, aes(x=FWT, y=FDKO)) +
#   geom_point(size=1, shape=23) +
#   geom_point(data=averagedf2_k5, aes(x=FWT, y=FDKO), shape = 19, color ='red', size =1) +
#   theme_classic()
# 
# averagedf2_k6 <- merge(k6,averagedf2, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')
# p6 <- ggplot(averagedf, aes(x=FWT, y=FdCD)) +
#   geom_point(size=1, shape=23) +
#   geom_point(data=averagedf2_k6, aes(x=FWT, y=FdCD), shape = 19, color ='red', size =1) +
#   theme_classic()
# 
# averagedf2_k7 <- merge(k7,averagedf2, by.x = 'ensembl_gene_id', by.y = 'ensembl_gene_id')
# p7 <- ggplot(averagedf, aes(x=NWT, y=FWT)) +
#   geom_point(size=1, shape=23) +
#   geom_point(data=averagedf2_k7, aes(x=NWT, y=FWT), shape = 19, color ='red', size =1) +
#   theme_classic()
# p7
# 
# grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)
# 
# # comparison <- cbind(averagedf[,8], averagedf[,6])
# # colnames(comparison) <- c("NWT", "NdCD")
# # rownames(comparison) <- c(row.names(FWT))
# # comparison <- as.data.frame(comparison)
# df_lcpm <- annotlcpm
# annotlcpmwt <- df_lcpm %>% mutate(across(NCKO:NWT) - NWT) %>% mutate(across(FCKO:FWT) - FWT)
# annotlcpmwt <- df_lcpm %>% mutate(across(FCKO:FWT) - FWT)
# rownames(annotlcpmwt) <- annotlcpmwt$ensembllistanno
# annotlcpmwt <- cbind(annotlcpmwt[2:4], annotlcpmwt[6:8])
# 
# 
# NFWT <- FWT$ave - NWT$ave
# NFCKO <- FCKO$ave - NCKO$ave
# NFDKO <- FDKO$ave - NDKO$ave
# NFdCD <-  FdCD$ave - NdCD$ave
# 
# FCKO_WT <- FCKO$ave - FWT$ave
# FDKO_WT <- FDKO$ave - FWT$ave
# FdCD_WT <- FdCD$ave - FWT$ave
# 
# NCKO_WT <- NCKO$ave - NWT$ave
# NDKO_WT <- NDKO$ave - NWT$ave
# NdCD_WT <- NdCD$ave - NWT$ave
# 
# genotypenorm <- cbind(NFWT, NFCKO, NFDKO, NFdCD)
# colnames(genotypenorm) <- c("NFWT", "NFCKO", "NFDKO", "NFdCD")
# rownames(genotypenorm) <- c(row.names(FWT))
# genotypenorm <- as.data.frame(genotypenorm)
# WTnorm <- genotypenorm
# 
# WTnorm <- cbind(FCKO_WT,FDKO_WT,FdCD_WT,NCKO_WT,NDKO_WT,NdCD_WT)
# colnames(WTnorm) <- c("FCKO","FDKO", "FdCD","NCKO","NDKO","NdCD")
# rownames(WTnorm) <- c(row.names(FWT))
# WTnorm <- as.data.frame(WTnorm)
# 
# #Import list of genes of interst, sig genes etc.
# 
# wtsiggenes <- read.table("~/Desktop/Diffgenes/DiffgeneschenNWTvsFWT-05.20.20.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
#            stringsAsFactor=FALSE)
# #wtsiggenes <- read.table("~/Desktop/Diffgenes/DiffgenesNWTvsFWT-02.23.20.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
# #                         stringsAsFactor=FALSE)
# mergedsiggenes <- merge(genotypenorm, wtsiggenes, by = 0)
# 
# wtsiggenes <- read.table("~/Desktop/Diffgenes/masterexpressiontable.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
#                          stringsAsFactor=FALSE)
# 
# subsetgenes <- read.table("~/Desktop/Subs1.csv", header=F, sep=",", row.names = 0, check.name=FALSE,
#            stringsAsFactor=FALSE)
# 
# #filter genes for significant genes or LogFC cut off highlight with geom point overlay
# sample1 <- NDKO_WT
# sample2 <- NdCD_WT
# WTnormfiltered <- WTnorm %>% filter(sample1 > 1 | sample1 < -1 | sample2 > 1 | sample2 < -1)
# #WTnormfiltered$FDKO <- WTnorm %>% na_if(WTnorm$FDKO > -1 & WTnorm$FDKO < 1)
# 
# 
# pg1 <- ggplot(WTnorm, aes(x=NFWT, y=NFCKO)) +
#   ggtitle("WT normalized Log2cpm, CKO vs WT") +
#   labs(x='FWT/NWT Log2cpm', y = "FCKO/NCKO Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   xlim(-5,5) +
#   ylim(-5,5) +
#   theme_classic()
# 
# geom_point(data=mergedsiggenes, aes(x=NFWT, y=NFCKO), shape = 19, color ='black', size =1) +
# 
# sample1 <- NFWT
# sample2 <- NFDKO
# WTnormfiltered <- WTnorm %>% filter(sample1 > 1 | sample1 < -1 | sample2 > 1 | sample2 < -1)
# 
# 
# pg2 <- ggplot(WTnorm, aes(x=NFWT, y=NFDKO)) +
#   ggtitle("WT normalized Log2cpm, DKO vs WT") +
#   labs(x='FWT/NWT Log2cpm', y = "FDKO/NDKO Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   xlim(-5,5) +
#   ylim(-5,5) +
#   theme_classic()
# 
# pg2
# geom_point(data=mergedsiggenes, aes(x=NFWT, y=NFDKO), shape = 19, color ='black', size =1) +
# 
# sample1 <- NFWT
# sample2 <- NFdCD
# WTnormfiltered <- WTnorm %>% filter(sample1 > 1 | sample1 < -1 | sample2 > 1 | sample2 < -1)
# 
# pg3 <- ggplot(WTnorm, aes(x=NFWT, y=NFdCD)) +
#   ggtitle("WT normalized Log2cpm, dCD vs WT") +
#   labs(x='FWT/NWT Log2cpm', y = "FdCD/NdCD Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   xlim(-5,5) +
#   ylim(-5,5) +
#   theme_classic()
# 
# geom_point(data=mergedsiggenes, aes(x=NFWT, y=NFdCD), shape = 19, color ='black', size =1) +
# 
# 
# grid.arrange(pg1, pg2, pg3, nrow = 1)
# 
# 
# 
# 
# 
# t1 <- ggplot(WTnorm, aes(x=FDKO, y=FCKO)) +
#   ggtitle("WT normalized Log2cpm, FCKO vs FDKO") +
#   labs(x='FDKO/WT Log2cpm', y = "FCKO/WT Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   xlim(-5,5) +
#   ylim(-5,5) +
#   theme_classic()
# 
# 
# t2 <- ggplot(WTnorm, aes(x=FDKO, y=FdCD)) +
#   ggtitle("WT normalized Log2cpm, FDKO vs FdCD") +
#   labs(x='FDKO/WT Log2cpm', y = "FdCD/WT Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_point(data=WTnormfiltered, aes(x=FDKO, y=FdCD), shape = 19, color ='black', size =1) +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   xlim(-5,5) +
#   ylim(-5,5) +
#   theme_classic()
# 
# t3 <- ggplot(WTnorm, aes(x=FDKO, y=FdCD)) +
#   ggtitle("WT normalized Log2cpm, FDKO vs FdCD") +
#   labs(x='FDKO/WT Log2cpm', y = "FdCD/WT Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_point(data=WTnormfiltered, aes(x=FDKO, y=FdCD), shape = 19, color ='black', size =1) +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   theme_classic()
# 
# t4 <- ggplot(WTnorm, aes(x=NDKO, y=NCKO)) +
#   ggtitle("WT normalized Log2cpm, NDKO vs NCKO") +
#   labs(x='NDKO/WT Log2cpm', y = " NCKO/WT Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_point(data=WTnormfiltered, aes(x=NDKO, y=NCKO), shape = 19, color ='black', size =1) +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   theme_classic()
# 
# t5 <- ggplot(WTnorm, aes(x=NDKO, y=NdCD)) +
#   ggtitle("WT normalized Log2cpm, NDKO vs NdCD") +
#   labs(x='NDKO/WT Log2cpm', y = "NdCD/WT Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_point(data=WTnormfiltered, aes(x=NDKO, y=NdCD), shape = 19, color ='black', size =1) +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   xlim(-5,5) +
#   ylim(-5,5) +
#   theme_classic()
# 
# t6 <- ggplot(WTnorm, aes(x=FDKO, y=FdCD)) +
#   ggtitle("WT normalized Log2cpm, FDKO vs FdCD") +
#   labs(x='FDKO/WT Log2cpm', y = "FdCD/WT Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_point(data=WTnormfiltered, aes(x=FDKO, y=FdCD), shape = 19, color ='black', size =1) +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   theme_classic()
# 
# grid.arrange(t5, t2, nrow =1)
# 
# t7 <- ggplot(WTnorm, aes(x=NDKO, y=FDKO)) +
#   ggtitle("WT normalized Log2cpm, NDKO vs FDKO") +
#   labs(x='NDKO/NWT Log2cpm', y = "FDKO/FWT Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_point(data=WTnormfiltered, aes(x=NDKO, y=FDKO), shape = 19, color ='black', size =1) +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   xlim(-5,5) +
#   ylim(-5,5) +
#   geom_smooth(data=WTnormfiltered, method ="lm") +
#   theme_classic()
# 
# cor.test(WTnormfiltered$NDKO, WTnormfiltered$FDKO, method = "pearson", conf.level = 0.95)
# cor.test(WTnorm$NDKO, WTnorm$FDKO, method = "pearson", conf.level = 0.95)
# 
# sample1 <- NdCD_WT
# sample2 <- FdCD_WT
# WTnormfiltered <- WTnorm %>% filter(sample1 > 1 | sample1 < -1 | sample2 > 1 | sample2 < -1)
# 
# t8 <- ggplot(WTnorm, aes(x=NdCD, y=FdCD)) +
#   ggtitle("WT normalized Log2cpm, NdCD vs FdCD") +
#   labs(x='NdCD/NWT Log2cpm', y = "FdCD/FWT Log2cpm") +
#   geom_point(size=1, shape=19, color = 'gray') +
#   geom_point(data=WTnormfiltered, aes(x=NdCD, y=FdCD), shape = 19, color ='black', size =1) +
#   geom_hline(yintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_hline(yintercept = -1, color = "blue", linetype = 'dashed' ) +   
#   geom_vline(xintercept = 1, color = "blue", linetype = 'dashed' ) +
#   geom_vline(xintercept = -1, color = "blue", linetype = 'dashed' ) +
#   xlim(-5,5) +
#   ylim(-5,5) +
#   geom_smooth(data=WTnormfiltered, method ="lm") +
#   theme_classic()
# 
# grid.arrange(t7, t8, nrow =1)
# 
# cor.test(WTnormfiltered$NdCD, WTnormfiltered$FdCD, method = "pearson", conf.level = 0.95)
# cor.test(WTnorm$NdCD, WTnorm$FdCD, method = "pearson", conf.level = 0.95)
# 
# k1 <- read.table("~/Desktop/Diffgenes/DiffgenesNWTvsNCKOv2-05.21.2020.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
#            stringsAsFactor=FALSE)
# k2 <- read.table("~/Desktop/Diffgenes/DiffgenesNWTvsNDKOv2-05.21.2020.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
#                  stringsAsFactor=FALSE)
# k3 <- read.table("~/Desktop/Diffgenes/DiffgenesNWTvsNdCDv2-05.21.2020.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
#                  stringsAsFactor=FALSE)
# k4 <- read.table("~/Desktop/Diffgenes/DiffgenesFWTvsFCKOv2-05.21.2020.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
#                  stringsAsFactor=FALSE)
# k5 <- read.table("~/Desktop/Diffgenes/DiffgenesFWTvsFDKOv2-05.21.2020.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
#                  stringsAsFactor=FALSE)
# k6 <- read.table("~/Desktop/Diffgenes/DiffgenesFWTvsFdCDv2-05.21.2020.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
#                  stringsAsFactor=FALSE)
# k7 <- read.table("~/Desktop/Diffgenes/DiffgenesNWTvsFWTv2-05.21.2020.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
#                  stringsAsFactor=FALSE)
# 
# #set up initial df
# k1 <- cbind(rownames(k1), k1[2], k1[6])
# k2 <- cbind(rownames(k2), k2[2], k2[6])
# colnames(k1) <- c("ensembl_gene_id", "log2FoldChange", "padj")
# colnames(k2) <- c("ensembl_gene_id", "log2FoldChange", "padj")
# masterexptable2 <- merge(k1, k2, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = TRUE)
# 
# #change kx to add all k dataframes to single master exp table
# kx <- k7
# kx <- cbind(rownames(kx), kx[2], kx[6])
# colnames(kx) <- c("ensembl_gene_id", "log2FoldChange", "padj")
# masterexptable2 <- merge(masterexptable2, kx, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all = TRUE)
# 
# masterexptable2 <- read.csv("~/Desktop/Diffgenes/masterexptable2.anno.csv", header = T)
# masterexptable2 <- subset(masterexptable2, (!is.na(masterexptable2[,17]))) #filter NA gene names
# masterexptable3 <- cbind(masterexptable2[3:16])
# rownames(masterexptable3) <- c(as.character(masterexptable2$ensembllistanno)) 
# 
# colnames(masterexptable3) <- c("NCKOlog2FoldChange", "NCKOpadj","NDKOlog2FoldChange", "NDKOpadj","NdCDlog2FoldChange", "NdCDpadj","FCKOlog2FoldChange", "FCKOpadj","FDKOlog2FoldChange", "FDKOpadj","FdCDlog2FoldChange", "FdCDpadj", "FWTlog2FoldChange", "FWTpadj")
# 
# masterexptablesubset <- cbind(masterexptable3[3],masterexptable3[5], masterexptable3[9], masterexptable3[11])
# colnames(masterexptable2) <- c("ensembl_gene_id", "NDKOlog2FoldChange", "NDKOpadj","NdCDlog2FoldChange", "NdCDpadj","FDKOlog2FoldChange", "FDKOpadj","FdCDlog2FoldChange", "FdCDpadj")
# 
# Subs1 <- subset(masterexptablesubset, (!is.na(masterexptablesubset[,1])) & (!is.na(masterexptablesubset[,2])))
# Subs2 <- subset(masterexptablesubset, (!is.na(masterexptablesubset[,1])) & (!is.na(masterexptablesubset[,2])) & (!is.na(masterexptablesubset[,3])) & (!is.na(masterexptablesubset[,4])))
# Subs3 <- subset(masterexptablesubset, (!is.na(masterexptablesubset[,3])) & (!is.na(masterexptablesubset[,4])))
# Subs4 <- subset(masterexptablesubset, (!is.na(masterexptablesubset[,1])) & (!is.na(masterexptablesubset[,3])))
# 
# 
# 
# write.csv(Subs4, "~/Desktop/Subs4.csv")
# 
# library(ComplexHeatmap)
# 
# Subs4 <- read.csv("~/Desktop/Subs4.csv", header = T)
# 
# 
# masterlogtable <- cbind(masterexptable2[2], masterexptable2[4])
# masterannottable <- as.data.frame(cbind(Subs1[10]))
# row.names(masterannottable) <- masterannottable$ensembllistanno
# 
# masterlogtable <- as.matrix(masterlogtable)
# masterlogtable[is.na(masterlogtable)] <- 0
# 
# masterlogtable <- as.matrix(cbind(Subs1[2], Subs1[4]))
# 
# Subs1 <- cbind(Subs1[1], Subs1[2])
# Subs3 <- cbind(Subs3[3], Subs3[4])
# Subs4 <- cbind(Subs4[1], Subs4[3])
# 
# rowha <- HeatmapAnnotation(df = masterannottable, 
#                            which = "row",
#                            show_legend = FALSE,
#                            show_annotation_name = FALSE
# )
# 
# Heatmap(Subs4, 
#         name = "Log2foldchange", #title of legend
#         #column_title = "Samples", 
#         row_title = "Gene name",
#         row_names_gp = gpar(fontsize = 7),
#         column_labels = c("1", "2"),
#         #right_annotation = rowha,
#         column_dend_reorder = FALSE
# )
# 
# annotlcpm <- read.csv("~/Desktop/Diffgenes/avedf_lcpmv2.anno.csv", header =T)
# annotlcpm <- distinct(annotlcpm, annotlcpm$ensembllistanno, .keep_all = TRUE)
# annotlcpm <- subset(annotlcpm, (!is.na(annotlcpm[,10])))
# filteredlcpm <- as.data.frame(annotlcpm[2:9])
# rownames(filteredlcpm) <- c(as.character(annotlcpm$ensembllistanno))
# 
# filteredlcpmzscore <- sapply(filteredlcpm, scale)
# filteredlcpmzscore <- as.data.frame(filteredlcpmzscore)
# rownames(filteredlcpmzscore) <- rownames(filteredlcpm)
# 
# 
# genesofinterest <- c("Porcn", "Smad2", "Smad3", "Axin2", "Axin1", "Cdx1", "Cdx2", "Cdx4", "Bmp4", "Dvl1", "Dvl2", "Dvl3", "Ctnnb1", "Sp5", "Tcf3", "Lef1", "Tcf1", "Tcf7l1",
#                      "Stat3", "Wnt3", "Dkk1", "Dkk2", "Dkk3", "Dkk4", "Wnt5a", "Wnt5b", "Wnt8a", "Wnt10a", "Wnt10b", "Wnt11", "Wnt4", "Wnt6", "Wnt7b", "Wnt9a")
# genesofinterest <- c("Tcf7l1","Sall1", "Sall4", "Pou5f1", "Esrrb", "Klf2", "Klf4", "Sox2", "Nanog", "Tcf3", "Lef1", "Runx1", "Fgf5", "Dnmt3a", "Dnmt3b", "Dnmt3l", "Otx2", "Pou3f1", "T", "Twist1", "Tbx3", "Sox9", "Sox17")
# genesofinterest <- c("Hdac1", "Hdac2", "Mta1", "Mta2", "Mta3", "Rbbp7", "Rbbp4", "Chd1", "Chd2", "Chd3", "Mbd2", "Mbd3")
# genesofinterest <- c("Kmt2a","Kmt2b", "Kmt2c", "Kmt2d", "Men1","Ash2l","Dpy30","Rbbp5","Kdm6a", "Pagr1a", "Paxip1", "Wdr5", "Kdm6a" )
# genesofinterest <- c("Fzd1", "Fzd2","Fzd3","Fzd4","Fzd5","Fzd6","Fzd7","Fzd8","Fzd9","Fzd10", "Lrp5", "Lrp6")
# 
# 
# filteredlcpmzscoreselect <- subset(filteredlcpmzscore, rownames(filteredlcpmzscore) %in% genesofinterest)
# Heatmap(filteredlcpmzscoreselect)
# 
# annotlcpmselect <- subset(filteredlcpm, rownames(filteredlcpm) %in% genesofinterest)
# pheatmap(annotlcpmselect)
# # rownames(annotlcpmselect) <- c(as.character(annotlcpm$ensembllistanno))
# 
# Heatmap(annotlcpmselect,
#         row_title = "COMPASS complex")
# 
# 
# 
# annotlcpm <- read.csv("~/Desktop/Diffgenes/avedf_lcpmv2.anno.csv", header =T)
# annotlcpm <- distinct(annotlcpm, annotlcpm$ensembllistanno, .keep_all = TRUE)
# annotlcpm <- subset(annotlcpm, (!is.na(annotlcpm[,10])))
# filteredlcpm <- as.data.frame(annotlcpm[2:9])
# rownames(filteredlcpm) <- c(as.character(annotlcpm$ensembllistanno))
# 
# wtlcpm <- read.csv("~/Desktop/Diffgenes/DiffgenesNWTvsFWTv3-07.07.20.anno.csv", header =T)
# genesofinterest <- c(as.character(wtlcpm$ensembllistanno))
# 
# annotlcpmselect <- subset(filteredlcpm, rownames(filteredlcpm) %in% genesofinterest)
# pheatmap(annotlcpmselect, show_rownames = FALSE, scale = "row", main = "Significant Naive to Formative genes")
# 
# pheatmap(filteredlcpm, show_rownames = FALSE, scale = "row", main = "All Detected Genes")
# 
# 
# filteredlcpmzscore <- sapply(filteredlcpm, scale)
# filteredlcpmzscore <- as.data.frame(filteredlcpmzscore)
# colnames(filteredlcpmzscore) <- c("FCKO", "FdCD","FDKO","FWT","NCKO","NdCD","NDKO","NWT")
# rownames(filteredlcpmzscore) <- rownames(filteredlcpm)
# filteredlcpmzscoreselect <- subset(filteredlcpmzscore, rownames(filteredlcpmzscore) %in% genesofinterest)
# pheatmap(filteredlcpmzscoreselect, show_rownames = FALSE, main = "Significant Naive to Formative genes", scale = "row")
# pheatmap(filteredlcpmzscore, show_rownames = FALSE, main = "All Detected Genes", scale = "row")
# 
# output <- "~/Desktop/Diffgenes/DiffgenesNWTvsNCKOv3-07.07.20.anno.csv"
# genetable <- read.csv("~/Desktop/Diffgenes/DiffgenesNWTvsNCKOv3-07.07.20.csv", header =T, row.names = 1)
# converted <- convertensemblgenetosymbol(genetable)
# write.csv(converted, output)
# 
# 
# 
# uniqueNdCD <- as.data.frame(setdiff(k3$ensembllistanno, k2$ensembllistanno))
# write.csv(uniqueNdCD, "~/Desktop/naiveNdCDunique.csv")
# 
# 
# uniqueNDKO <- as.data.frame(setdiff(k2$ensembllistanno, k3$ensembllistanno))
# write.csv(uniqueNDKO, "~/Desktop/NDKOunique.csv")
# 
# uniquedCD <- as.data.frame(setdiff(k6$ensembllistanno, k5$ensembllistanno))
# write.csv(sharedfiltered, "~/Desktop/NDKOdCDshared.csv")






