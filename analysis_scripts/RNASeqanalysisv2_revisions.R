#Code used to further analyze and visualize RNAseq data from differentialgeneanalysis

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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] patchwork_1.1.1 umap_0.2.8.0    ggrepel_0.9.1   gridExtra_2.3   dplyr_1.0.9     ggplot2_3.3.6  
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.8.3        RSpectra_0.16-1     pillar_1.7.0        compiler_4.1.2      tools_4.1.2         jsonlite_1.8.0      lifecycle_1.0.1    
# [8] tibble_3.1.7        gtable_0.3.0        lattice_0.20-45     png_0.1-7           pkgconfig_2.0.3     rlang_1.0.3         Matrix_1.4-1       
# [15] cli_3.3.0           DBI_1.1.3           rstudioapi_0.13     withr_2.5.0         generics_0.1.3      GlobalOptions_0.1.2 vctrs_0.4.1        
# [22] askpass_1.1         tidyselect_1.1.2    reticulate_1.25     glue_1.6.2          R6_2.5.1            fansi_1.0.3         purrr_0.3.4        
# [29] magrittr_2.0.3      scales_1.2.0        ellipsis_0.3.2      assertthat_0.2.1    shape_1.4.6         circlize_0.4.15     colorspace_2.0-3   
# [36] utf8_1.2.2          openssl_2.0.2       munsell_0.5.0       crayon_1.5.1   
#####

library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggrepel)
library(umap)
library(patchwork)

#required functions
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

#####
####Required for analyses:Import replicate lcpm dataframe#####

replicatelcpm <- read.csv("~/Desktop/data/RBRB01replicatelcpm.anno.051721.csv", header =T)
replicatelcpm <- distinct(replicatelcpm, replicatelcpm$ensembllistanno, .keep_all = TRUE)
replicatelcpm <- subset(replicatelcpm, (!is.na(replicatelcpm[,26])))
filteredreplicatelcpm <- as.data.frame(replicatelcpm[2:25])
rownames(filteredreplicatelcpm) <- c(as.character(replicatelcpm$ensembllistanno))

#####
####Required for analyses:Import average lcpm dataframe#####

annotlcpm <- read.csv("~/Desktop/data/Diffgenes/avedf_lcpmv2.anno.csv", header =T)
annotlcpm <- distinct(annotlcpm, annotlcpm$ensembllistanno, .keep_all = TRUE)
annotlcpm <- subset(annotlcpm, (!is.na(annotlcpm[,10])))
filteredlcpm <- as.data.frame(annotlcpm[2:10])
rownames(filteredlcpm) <- c(as.character(annotlcpm$ensembllistanno))

filteredlcpm$rnafold_geno_wt <- filteredlcpm$FWT - filteredlcpm$NWT
filteredlcpm$rnafold_geno_cko <- filteredlcpm$FCKO - filteredlcpm$NCKO
filteredlcpm$rnafold_geno_dko <- filteredlcpm$FDKO - filteredlcpm$NDKO
filteredlcpm$rnafold_geno_dcd <- filteredlcpm$FdCD - filteredlcpm$NdCD

filteredlcpm$FCKO_WT <- filteredlcpm$FCKO - filteredlcpm$FWT
filteredlcpm$FDKO_WT <- filteredlcpm$FDKO - filteredlcpm$FWT
filteredlcpm$FdCD_WT <- filteredlcpm$FdCD - filteredlcpm$FWT

filteredlcpm$NCKO_WT <- filteredlcpm$NCKO - filteredlcpm$NWT
filteredlcpm$NDKO_WT <- filteredlcpm$NDKO - filteredlcpm$NWT
filteredlcpm$NdCD_WT <- filteredlcpm$NdCD - filteredlcpm$NWT

filteredlcpm$NDKO_NCKO <- filteredlcpm$NDKO - filteredlcpm$NCKO
filteredlcpm$FDKO_FCKO <- filteredlcpm$FDKO - filteredlcpm$FCKO

#####
####Required for analyses:Import DESeq2 results####

siggenetables=list()
k1 <- read.table("~/Desktop/data/Diffgenes/DiffgenesNWTvsNCKOv3-07.07.20.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                 stringsAsFactor=FALSE)
k2 <- read.table("~/Desktop/data/Diffgenes/DiffgenesNWTvsNDKOv3-07.07.20.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                 stringsAsFactor=FALSE)
k3 <- read.table("~/Desktop/data/Diffgenes/DiffgenesNWTvsNdCDv3-07.07.20.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                 stringsAsFactor=FALSE)
k4 <- read.table("~/Desktop/data/Diffgenes/DiffgenesFWTvsFCKOv3-07.07.20.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                 stringsAsFactor=FALSE)
k5 <- read.table("~/Desktop/data/Diffgenes/DiffgenesFWTvsFDKOv3-07.07.20.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                 stringsAsFactor=FALSE)
k6 <- read.table("~/Desktop/data/Diffgenes/DiffgenesFWTvsFdCDv3-07.07.20.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                 stringsAsFactor=FALSE)
k7 <- read.table("~/Desktop/data/Diffgenes/DiffgenesNWTvsFWTv3-07.07.20.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                 stringsAsFactor=FALSE)
k8 <- read.table("~/Desktop/data/Diffgenes/DiffgenesNDKOvsFDKOfull-07.07.20.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                 stringsAsFactor=FALSE)
k9 <- read.table("~/Desktop/data/Diffgenes/DiffgenesNCKOvsNDKO-11.20.22.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
              stringsAsFactor=FALSE)
k10 <- read.table("~/Desktop/data/Diffgenes/DiffgenesFCKOvsFDKO-11.20.22.anno.csv", header=T, sep=",", row.names = 1, check.name=FALSE,
                  stringsAsFactor=FALSE)


siggene_group <- list(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10)

filteredlcpm_wtsig <- right_join(filteredlcpm, k7, by = "ensembllistanno", keep=TRUE)
filteredlcpm_ndkosig <- right_join(filteredlcpm, k2, by = "ensembllistanno", keep=TRUE)
filteredlcpm_fdkosig <- right_join(filteredlcpm, k5, by = "ensembllistanno", keep=TRUE)

filteredlcpm_ndkockosig <- right_join(filteredlcpm, k9, by = "ensembllistanno", keep=TRUE)
filteredlcpm_fdkockosig <- right_join(filteredlcpm, k10, by = "ensembllistanno", keep=TRUE)



#filteredlcpm_wtsig.up <- filteredlcpm_wtsig %>% filter(rnafold_geno_wt > 1)
#filteredlcpm_wtsig.down <- filteredlcpm_wtsig %>% filter(rnafold_geno_wt < -1)

#ggplot() +
#  geom_point(data=filteredlcpm_wtsig.up, aes(x=rnafold_geno_wt, y=log2FoldChange))


#####
####Generate plots for WT vs genotype LCPM##########
#lcpmxlcpm plots#
#genesofinterest <- c("Grhl2", "Nanog", "Esrrb", "Nr0b1", "Tbx3", "Klf4", "Zfp42", "Otx2", "Pou3f1", "Fgf5")
genesofinterest <- c("Grhl2", "Esrrb", "Klf4", "Zfp42", "Otx2", "Pou3f1", "Fgf5")
genesofinterest <- c("Nanog", "Sox2")

filteredlcpmhighlight <- subset(filteredlcpm, rownames(filteredlcpm) %in% genesofinterest)
filteredlcpmhighlight$genesymbol <- rownames(filteredlcpmhighlight)
  
# compare_group_first=c("NWT","NWT","NWT","FWT", "FWT", "FWT", "NWT")
# compare_group_firstname=c("Naive WT","Naive WT","Naive WT", "Form WT", "Form WT","Form WT", "Naive WT")
# compare_group_second=c("NCKO","NDKO","NdCD","FCKO", "FDKO", "FdCD", "FWT") 
# compare_group_secondname=c("Naive MLL3KO","Naive DKO","Naive dCD","Form MLL3KO", "Form DKO", "Form dCD", "Form WT")
# compare_group_z <- c("NCKO_WT","NDKO_WT", "NdCD_WT", "FCKO_WT", "FDKO_WT", "FdCD_WT", "rnafold_geno_wt")

###saved as 3x9 pdf for later edit 

siggene_group <- list(k1, k2, k3, k9)
compare_group_first=c("NWT","NWT","NWT", "NCKO")
compare_group_firstname=c("Naive WT","Naive WT","Naive WT", "Naive MLL3KO")
compare_group_second=c("NCKO","NDKO","NdCD", "NDKO")
compare_group_secondname=c("Naive MLL3KO","Naive DKO","Naive dCD", "Naive DKO")
compare_group_z <- c("NCKO_WT","NDKO_WT", "NdCD_WT", "NDKO_NCKO")

comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname, compare_group_z))
siggene_group_counter = 0
plotlist=list()
combinedplotlist=list()

plotarray <- for(i in 1:nrow(comparematrix)) {
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  zgroup <- as.character(row[1,5])
  siggene_group_counter = siggene_group_counter + 1
  
  
  
  siggenes <- merge(siggene_group[siggene_group_counter],filteredlcpm, by = 'ensembllistanno')
  #zgenes <- eval(parse(text = zgroup))
  filteredlcpm_wtsig.up <- siggenes %>% filter(eval(parse(text = zgroup)) > 1)
  filteredlcpm_wtsig.down <- siggenes %>% filter(eval(parse(text = zgroup)) < -1)
  
  #countup <- filteredlcpm_wtsig.up %>% count(rnafold_geno_wt > 1)
  
  countup <- filteredlcpm_wtsig.up %>% count(eval(parse(text = zgroup)) > 1)
  countup <- paste(countup[1,2], "Up", sep=' ')
  grobup <- grobTree(textGrob(countup, x=0.05,  y=0.93, just="left",
                              gp=gpar(col="deepskyblue3", fontsize=13)))
  
  #countdown <- filteredlcpm_wtsig.down %>% count(rnafold_geno_wt < -1)
  
  countdown <- filteredlcpm_wtsig.down %>% count(eval(parse(text = zgroup)) < -1)
  countdown <- paste(countdown[1,2], "Down", sep=' ')
  grobdown <- grobTree(textGrob(countdown, x=0.95,  y=0.07, just="right",
                                gp=gpar(col="deepskyblue3", fontsize=13)))
  
  
  
  # siggenes <- merge(siggene_group[siggene_group_counter],annotlcpm, by = 'ensembllistanno')
  # siggenes2 <- subset(siggenes, log2FoldChange > 1 | log2FoldChange < -1)
  # countup <- siggenes %>% count(log2FoldChange > 1)
  # countup <- paste(countup[2,2], "Up", sep=' ')
  # grobup <- grobTree(textGrob(countup, x=0.05,  y=0.93, just="left",
  #                             gp=gpar(col="deepskyblue3", fontsize=13)))
  # countdown <- siggenes %>% count(log2FoldChange < -1)
  # countdown <- paste(countdown[2,2], "Down", sep=' ')
  # grobdown <- grobTree(textGrob(countdown, x=0.95,  y=0.07, just="right",
  #                               gp=gpar(col="deepskyblue3", fontsize=13)))
  
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot(data=filteredlcpm, aes_string(x=xsample, y=ysample)) +
                                                  geom_point(size=0.2, shape=23, color = "gray") +
                                                  geom_point(data=filteredlcpm_wtsig.up, aes_string(x=xsample, y=ysample), shape = 19, color ='deepskyblue3', size =0.2) +
                                                  geom_point(data=filteredlcpm_wtsig.down, aes_string(x=xsample, y=ysample), shape = 19, color ='deepskyblue3', size =0.2) +
                                                  geom_text_repel(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample, label="genesymbol", size = 12), color = "black", force = 75, box.padding = 1)+ #, color = "black", force = 30) +
                                                  geom_point(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample), color = 'black') +
                                                  annotation_custom(grobup) +
                                                  annotation_custom(grobdown) +
                                                  xlab(paste(xname, "(Log2 CPM)")) +
                                                  ylab(paste(yname, "(Log2 CPM)")) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none"))
  combinedplotlist[[siggene_group_counter]] <- plotlist[[siggene_group_counter]]
  
}
wrap_plots(combinedplotlist)

###form lcpm plot
compare_group_first=c("FWT", "FWT", "FWT", "FCKO")
compare_group_firstname=c("Form WT", "Form WT","Form WT", "Form MLL3KO")
compare_group_second=c("FCKO", "FDKO", "FdCD", "FDKO")
compare_group_secondname=c("Form MLL3KO", "Form DKO", "Form dCD", "Form MLL3KO")
compare_group_z <- c("FCKO_WT", "FDKO_WT", "FdCD_WT", "FDKO_FCKO")

siggene_group <- list(k4, k5, k6, k10)
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname, compare_group_z))
siggene_group_counter = 0
plotlist=list()
combinedplotlist=list()

plotarray <- for(i in 1:nrow(comparematrix)) {
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  zgroup <- as.character(row[1,5])
  siggene_group_counter = siggene_group_counter + 1
  
  
  
  siggenes <- merge(siggene_group[siggene_group_counter],filteredlcpm, by = 'ensembllistanno')
  #zgenes <- eval(parse(text = zgroup))
  filteredlcpm_wtsig.up <- siggenes %>% filter(eval(parse(text = zgroup)) > 1)
  filteredlcpm_wtsig.down <- siggenes %>% filter(eval(parse(text = zgroup)) < -1)
  
  countup <- filteredlcpm_wtsig.up %>% count(eval(parse(text = zgroup)) > 1)
  countup <- paste(countup[1,2], "Up", sep=' ')
  grobup <- grobTree(textGrob(countup, x=0.05,  y=0.93, just="left",
                              gp=gpar(col="magenta3", fontsize=13)))
  
  #countdown <- filteredlcpm_wtsig.down %>% count(rnafold_geno_wt < -1)
  
  countdown <- filteredlcpm_wtsig.down %>% count(eval(parse(text = zgroup)) < -1)
  countdown <- paste(countdown[1,2], "Down", sep=' ')
  grobdown <- grobTree(textGrob(countdown, x=0.95,  y=0.07, just="right",
                                gp=gpar(col="magenta3", fontsize=13)))
  
  
  
  # siggenes <- merge(siggene_group[siggene_group_counter],annotlcpm, by = 'ensembllistanno')
  # siggenes2 <- subset(siggenes, log2FoldChange > 1 | log2FoldChange < -1)
  # countup <- siggenes %>% count(log2FoldChange > 1)
  # countup <- paste(countup[2,2], "Up", sep=' ')
  # grobup <- grobTree(textGrob(countup, x=0.05,  y=0.93, just="left",
  #                             gp=gpar(col="deepskyblue3", fontsize=13)))
  # countdown <- siggenes %>% count(log2FoldChange < -1)
  # countdown <- paste(countdown[2,2], "Down", sep=' ')
  # grobdown <- grobTree(textGrob(countdown, x=0.95,  y=0.07, just="right",
  #                               gp=gpar(col="deepskyblue3", fontsize=13)))
  
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot(data=filteredlcpm, aes_string(x=xsample, y=ysample)) +
                                                  geom_point(size=0.2, shape=23, color = "gray") +
                                                  geom_point(data=filteredlcpm_wtsig.up, aes_string(x=xsample, y=ysample), shape = 19, color ='magenta3', size =0.2) +
                                                  geom_point(data=filteredlcpm_wtsig.down, aes_string(x=xsample, y=ysample), shape = 19, color ='magenta3', size =0.2) +
                                                  geom_text_repel(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample, label="genesymbol", size = 12), color = "black", force = 175, box.padding = 1, max.overlaps = 100)+ #, color = "black", force = 30) +
                                                  geom_point(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample), color = 'black') +
                                                  annotation_custom(grobup) +
                                                  annotation_custom(grobdown) +
                                                  xlab(paste(xname, "(Log2 CPM)")) +
                                                  ylab(paste(yname, "(Log2 CPM)")) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none"))
  combinedplotlist[[siggene_group_counter]] <- plotlist[[siggene_group_counter]]
  
}
wrap_plots(combinedplotlist)




##naive form lcpm plot
compare_group_first=c("NWT", "NDKO")
compare_group_firstname=c("Naive WT", "Naive DKO")
compare_group_second=c("FWT", "FDKO")
compare_group_secondname=c("Form WT", "Form DKO")
compare_group_z <- c("rnafold_geno_wt", "rnafold_geno_dko")

siggene_group <- list(k7, k8)
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname, compare_group_z))
siggene_group_counter = 0
plotlist=list()
combinedplotlist=list()
#for proper wt counts
plotarray <- for(i in 1:nrow(comparematrix)) {
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  zgroup <- as.character(row[1,5])
  siggene_group_counter = siggene_group_counter + 1
  
  
  
  siggenes <- merge(siggene_group[siggene_group_counter],filteredlcpm, by = 'ensembllistanno')
  #zgenes <- eval(parse(text = zgroup))
  filteredlcpm_wtsig.up <- siggenes %>% filter(eval(parse(text = zgroup)) > 1)
  filteredlcpm_wtsig.down <- siggenes %>% filter(eval(parse(text = zgroup)) < -1)
  
  countup <- filteredlcpm_wtsig.up %>% count(eval(parse(text = zgroup)) > 1)
  countup <- paste(countup[1,2], "Up", sep=' ')
  grobup <- grobTree(textGrob(countup, x=0.05,  y=0.93, just="left",
                              gp=gpar(col="magenta3", fontsize=13)))
  
  countdown <- filteredlcpm_wtsig.down %>% count(eval(parse(text = zgroup)) < -1)
  countdown <- paste(countdown[1,2], "Down", sep=' ')
  grobdown <- grobTree(textGrob(countdown, x=0.95,  y=0.07, just="right",
                                gp=gpar(col="deepskyblue3", fontsize=13)))
  
  
  
  # siggenes <- merge(siggene_group[siggene_group_counter],annotlcpm, by = 'ensembllistanno')
  # siggenes2 <- subset(siggenes, log2FoldChange > 1 | log2FoldChange < -1)
  # countup <- siggenes %>% count(log2FoldChange > 1)
  # countup <- paste(countup[2,2], "Up", sep=' ')
  # grobup <- grobTree(textGrob(countup, x=0.05,  y=0.93, just="left",
  #                             gp=gpar(col="deepskyblue3", fontsize=13)))
  # countdown <- siggenes %>% count(log2FoldChange < -1)
  # countdown <- paste(countdown[2,2], "Down", sep=' ')
  # grobdown <- grobTree(textGrob(countdown, x=0.95,  y=0.07, just="right",
  #                               gp=gpar(col="deepskyblue3", fontsize=13)))
  
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot(data=filteredlcpm, aes_string(x=xsample, y=ysample)) +
                                                  geom_point(size=0.2, shape=23, color = "gray") +
                                                  geom_point(data=filteredlcpm_wtsig.up, aes_string(x=xsample, y=ysample), shape = 19, color ='magenta3', size =0.2) +
                                                  geom_point(data=filteredlcpm_wtsig.down, aes_string(x=xsample, y=ysample), shape = 19, color ='deepskyblue3', size =0.2) +
                                                  geom_text_repel(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample, label="genesymbol", size = 12), color = "black", force = 175, box.padding = 1, max.overlaps = 100)+ #, color = "black", force = 30) +
                                                  geom_point(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample), color = 'black') +
                                                  annotation_custom(grobup) +
                                                  annotation_custom(grobdown) +
                                                  xlab(paste(xname, "(Log2 CPM)")) +
                                                  ylab(paste(yname, "(Log2 CPM)")) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none"))
  combinedplotlist[[siggene_group_counter]] <- plotlist[[siggene_group_counter]]
  
}
wrap_plots(combinedplotlist)



#####
####Generate UMAP plots of samples####

#Filter all samples for only genes sig changed in WT n to form. 

# filteredreplicatelcpm <- filteredlcpm[2:5]
# rownames(filteredreplicatelcpm) <- filteredlcpm$ensembllistanno
# filteredlcpmWTsiggenes <- subset(filteredreplicatelcpm, rownames(filteredreplicatelcpm) %in% k7$ensembllistanno)

filteredlcpmWTsiggenes <- subset(filteredreplicatelcpm, rownames(filteredreplicatelcpm) %in% k7$ensembllistanno)
filteredlcpmWTsiggenes2 <- as.data.frame(c(filteredlcpmWTsiggenes[1:3], filteredlcpmWTsiggenes[7:15], filteredlcpmWTsiggenes[19:24]))
filteredreplicatelcpmt <- as.data.frame(t(filteredlcpmWTsiggenes2))



#filteredreplicatelcpmt <- as.data.frame(t(filteredreplicatelcpm))

#run UMAP package on df_lcpm, batchcorrected if necessary
counts.umap <- umap(filteredreplicatelcpmt, n_neighbors = 4) #n_neighbors = 5
head(counts.umap$layout)
#sample <- c("NWT2", "FWT2", "NDKO2", "FDKO2")
sample <- c("F MLL3KO","F MLL3KO","F MLL3KO", "F dCD","F dCD", "F dCD", "F DKO","F DKO","F DKO", "F WT", "F WT", "F WT", "N MLL3KO","N MLL3KO","N MLL3KO", "N dCD","N dCD","N dCD", "N DKO","N DKO","N DKO", "N WT","N WT","N WT" )
replicatename <- sample

sample <- c("F MLL3KO","F MLL3KO","F MLL3KO", "F DKO","F DKO","F DKO", "F WT", "F WT", "F WT", "N MLL3KO","N MLL3KO","N MLL3KO", "N DKO","N DKO","N DKO", "N WT","N WT","N WT" )
replicatename <- sample


layoutdf <- as.data.frame(counts.umap$layout)
colnames(layoutdf) <- c("one", "two")
layoutdf$samplename <- replicatename

# layoutdf$samplename <- factor(layoutdf$samplename, 
#                                        levels = c("N WT", "N MLL3KO", "N DKO", "N dCD", "F WT", "F MLL3KO", "F DKO", "F dCD" ))  
layoutdf$samplename <- factor(layoutdf$samplename, 
                              levels = c("N WT", "N MLL3KO", "N DKO", "F WT", "F MLL3KO", "F DKO"))  


#layoutdf <- layoutdf %>% filter(samplename == "N WT" | samplename == "F WT" | samplename == "N DKO" | samplename == "F DKO")


#300x800pixels
umap_naiveform <- ggplot(layoutdf, mapping = aes(x=one , y=two,fill = samplename)) +
  geom_point(size=3, aes(shape=samplename)) +
  scale_fill_manual(values = c("N WT" = "deepskyblue3", "F WT" = "magenta3", "N DKO" = "deepskyblue3", "F DKO" = "magenta3", "N MLL3KO" = "deepskyblue3", "F MLL3KO" = "magenta3")) +
  #  scale_fill_manual(values = c("N WT" = "deepskyblue3", "F WT" = "magenta3", "N DKO" = "deepskyblue3", "F DKO" = "magenta3", "N MLL3KO" = "deepskyblue3", "F MLL3KO" = "magenta3", "N dCD" = "deepskyblue3", "F dCD" = "magenta3" )) +
  #scale_shape_manual(values = c(21, 22, 23, 24, 21, 22, 23, 24)) +
  scale_shape_manual(values = c(21, 22, 23, 21, 22, 23)) +
    # #geom_text_repel(aes(label = samplename),
    #   box.padding = 0.1, 
    #   point.padding = 0.1,
    theme_bw() +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(aspect.ratio = 1/1.5, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2))
    # xlim(c(-10.5, 10.5)) +
    # ylim(c(-3.5,3.5)) 

umap_naiveform

#####
####Generate plots for genotypes normalized by pluripotent state####
# df_lcpm <- annotlcpm
# WT <- df_lcpm$FWT - df_lcpm$NWT 
# CKO <- df_lcpm$FCKO - df_lcpm$NCKO 
# DKO <- df_lcpm$FDKO - df_lcpm$NDKO 
# dCD <- df_lcpm$FdCD - df_lcpm$NdCD 
# df_annotlcpmstate <- as.data.frame(cbind(WT, CKO, DKO, dCD))
# rownames(df_annotlcpmstate) <- annotlcpm$ensembllistanno
# df_annotlcpmstate$genesymbol <- rownames(df_annotlcpmstate)
# df_annotlcpmstate$ensembllistanno <- rownames(df_annotlcpmstate)

#genesofinterest <- c("Dnmt3a", "Dnmt3b", "Dnmt3l", "Tet1", "Tet2", "Tet3")
genesofinterest <- c("Grhl2", "Esrrb", "Klf4", "Zfp42", "Otx2", "Pou3f1", "Fgf5")
#genesofinterest <- c("Suz12", "Ezh1", "Ezh2", "Eed", "Bmi1", "Rnf1", "Ring1", "RbAp48", "Jarid2", "Mel1r", "Rybp", "Cbx", "Usp7", "Bcor", "Pcgf1", "Pcgf2", "Pcgf4", "Smcl", "Kdm6a", "Dnmt3l", "Dnmt3b", "Dnmt3a")
#genesofinterest <- c("Grhl2", "Nanog", "Esrrb", "Nr0b1", "Nr5a2", "Tfcp2l1", "Tbx3", "Klf2", "Klf4", "Zfp42", "Otx2", "Pou3f1", "Sox11", "Lin28a", "Mir302b", "Fgf5")
filteredlcpmhighlight <- subset(filteredlcpm, ensembllistanno %in% genesofinterest)

# compare_group_first=c("WT", "WT", "WT")
# compare_group_firstname=c("Form.WT/NaiveWT", "Form.WT/NaiveWT", "Form.WT/NaiveWT")
# compare_group_second=c("CKO","DKO", "dCD") 
# compare_group_secondname=c("Form.CKO/NaiveCKO", "Form.DKO/NaiveDKO", "Form.dCD/NaivedCD")

compare_group_first=c("rnafold_geno_wt", "rnafold_geno_wt", "rnafold_geno_wt")
compare_group_firstname=c("Form.WT/NaiveWT", "Form.WT/NaiveWT", "Form.WT/NaiveWT")
compare_group_second=c("rnafold_geno_cko","rnafold_geno_dko", "rnafold_geno_dcd") 
compare_group_secondname=c("Form.CKO/NaiveCKO", "Form.DKO/NaiveDKO", "Form.dCD/NaivedCD")

compare_group_first=c("rnafold_geno_wt", "rnafold_geno_wt", "rnafold_geno_wt", "rnafold_geno_cko")
compare_group_firstname=c("Form.WT/NaiveWT", "Form.WT/NaiveWT", "Form.WT/NaiveWT", "Form.CKO/NaiveCKO")
compare_group_second=c("rnafold_geno_cko","rnafold_geno_dko", "rnafold_geno_dcd", "rnafold_geno_dko") 
compare_group_secondname=c("Form.CKO/NaiveCKO", "Form.DKO/NaiveDKO", "Form.dCD/NaivedCD", "Form.DKO/NaiveDKO")


###additional for reviews##
compare_group_first=c("rnafold_geno_cko", "rnafold_geno_dcd")
compare_group_firstname=c("Form.CKO/NaiveCKO", "Form.dCD/NaivedCD")
compare_group_second=c("rnafold_geno_dko","rnafold_geno_dko") 
compare_group_secondname=c("Form.DKO/NaiveDKO", "Form.DKO/NaiveDKO")

###


comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))

siggene_group_counter = 0
plotlist=list()

plotarray <- for(i in 1:nrow(comparematrix)) {
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  #siggene_group_counter = siggene_group_counter + 1
  siggene_group_counter = 7
  siggenes <- merge(siggene_group[siggene_group_counter],filteredlcpm, by = 'ensembllistanno')
  filteredlcpm_wtsig.up <- siggenes %>% filter(rnafold_geno_wt > 1)
  filteredlcpm_wtsig.down <- siggenes %>% filter(rnafold_geno_wt < -1)
  
  ct <- cor.test(filteredlcpm[,xsample], filteredlcpm[,ysample], method = "pearson") #, exact = FALSE
  ct_pearson <- round(as.numeric(ct[4]), digits = 2)
  #ct_pval <- round(as.numeric(ct[3]), digits = 2)
  
  # x <- xsample
  # y <- ysample
  # eq <- lm_eqn(filteredlcpm)
  # 
  grobcorr <- grobTree(textGrob((paste("r =", ct_pearson)), x=0.65,  y=0.10, just="left", gp=gpar(col="black", fontsize=12)))
  
  #grobpval <- grobTree(textGrob(paste("p.val", ct_pval, sep=' = '), x=0.75,  y=0.05, just="left",
  #    gp=gpar(col="black", fontsize=10)))
  #grobspearman <- grobTree(textGrob("Spearman", x=0.65,  y=0.17, just="left",
   #                                 gp=gpar(col="black", fontsize=12)))
  
  
  
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot(data=filteredlcpm, aes_string(x=xsample, y=ysample)) +
                                                  #geom_abline(slope = 1, color = "black", linetype = 'dashed' ) +
                                                  geom_point(size=0.5, shape=23, color = "gray", alpha = 0.2) +
                                                  geom_point(data=filteredlcpm_wtsig.up, aes_string(x=xsample, y=ysample), shape = 19, color ='magenta3', size =0.5, alpha = 0.2) +
                                                  geom_point(data=filteredlcpm_wtsig.down, aes_string(x=xsample, y=ysample), shape = 19, color ='deepskyblue3', size =0.5, alpha = 0.2) +
                                                  geom_text_repel(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample, label="ensembllistanno", size = 12), color = "black", force = 75, box.padding = 1) +
                                                  geom_point(data = filteredlcpmhighlight, aes_string(x=xsample, y=ysample), color = 'black') +
                                                  geom_smooth(method="lm", se = FALSE) +
                                                  #geom_text(x=6, y=1, label = eq, parse = TRUE) +
                                                  # annotation_custom(grobup) +
                                                  # annotation_custom(grobdown) +
                                                  xlab(paste(xname, "(Log2 CPM)")) +
                                                  ylab(paste(yname, "(Log2 CPM)")) +
                                                  # geom_hline(yintercept = 0, color = "black", linetype = 'dashed' ) +
                                                  #geom_vline(xintercept = 0, color = "black", linetype = 'dashed' ) +  
                                                  #geom_smooth(method='lm', formula= y~x, colour = "black", size =0.75) +
                                                  xlim(c(-7,7)) +
                                                  ylim(c(-7,7)) +
                                                  annotation_custom(grobcorr) +
                                                  #annotation_custom(grobspearman) +
                                                  #annotation_custom(grobpval) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none"))
  
}
wrap_plots(plotlist)





xyz <- as.data.frame(filteredlcpm[,c("rnafold_geno_wt","rnafold_geno_dcd")])
x <- filteredlcpm$rnafold_geno_wt
y <- filteredlcpm$rnafold_geno_dcd
lm_eqn(filteredlcpm)

xyz <- as.data.frame(filteredlcpm[,c("rnafold_geno_cko","rnafold_geno_dko")])
x <- filteredlcpm$rnafold_geno_cko
y <- filteredlcpm$rnafold_geno_dko
lm_eqn(filteredlcpm)

xyz <- as.data.frame(filteredlcpm[,c("rnafold_geno_wt","rnafold_geno_dcd")])
x <- filteredlcpm$rnafold_geno_wt
y <- filteredlcpm$rnafold_geno_dcd
lm_eqn(filteredlcpm)

#wtcko -> y = 0.046 + 0.76x
#wtdko -> y = 0.2 + 0.52x
#wtdcd -> y = 0.0018 + 0.67x

# subsetplotlist <- plotlist[seq(1,3)]
# grid.arrange(grobs=subsetplotlist, nrow=1)
# 
# p1 <- plotlist[1]
# p2 <- plotlist[2]
# 
# p1 + p2



#####
#####Analysis of N DKOWT and F DKOWT for premature and non-activated genes#####

naive_wtgenes <- filteredlcpm %>% filter((ensembllistanno %in% k7$ensembllistanno) & rnafold_geno_wt < -1)
form_wtgenes <- filteredlcpm %>% filter((ensembllistanno %in% k7$ensembllistanno) & rnafold_geno_wt > 1)
##naive diff gene cat plots
naive_dkogenesdown <- filteredlcpm %>% filter((ensembllistanno %in% k2$ensembllistanno) & NDKO_WT < -1)
naive_dkogenesup <- filteredlcpm %>% filter((ensembllistanno %in% k2$ensembllistanno) & NDKO_WT > 1)

a1 <- naive_dkogenesdown %>% filter(ensembllistanno %in% naive_wtgenes$ensembllistanno)
length(a1$ensembllistanno) #368
a1$phenotype <- "down" 
a1$state <- "naive" 

a2 <- naive_dkogenesdown %>% filter(ensembllistanno %in% form_wtgenes$ensembllistanno)
length(a2$ensembllistanno) #58
a2$phenotype <- "down" 
a2$state <- "form" 

a12 <- rbind(a1,a2) 
a12not <- naive_dkogenesdown %>% filter(!ensembllistanno %in% a12$ensembllistanno)
a12not$phenotype <- "down" 
a12not$state <- "other" 


a3 <- naive_dkogenesup %>% filter(ensembllistanno %in% naive_wtgenes$ensembllistanno)
length(a3$ensembllistanno) #39
a3$phenotype <- "up" 
a3$state <- "naive" 

a4 <- naive_dkogenesup %>% filter(ensembllistanno %in% form_wtgenes$ensembllistanno)
length(a4$ensembllistanno) #275
a4$phenotype <- "up" 
a4$state <- "form" 

a34 <- rbind(a3,a4)
a34not <- naive_dkogenesup %>% filter(!ensembllistanno %in% a34$ensembllistanno)
a34not$phenotype <- "up" 
a34not$state <- "other" 


aa <- rbind(a12,a12not,a34,a34not)
aa$phenotype <- factor(aa$phenotype, levels = c("up", "down"))
aa$state <- factor(aa$state, levels = c("naive", "form", "other"))

customcolors <- c("deepskyblue3", "magenta3","grey")
bp1 <- ggplot() +
  geom_bar(data=aa, aes(x=phenotype, fill=state)) +
  ylim(c(0,1000)) +
  scale_fill_manual(values=customcolors) + 
  theme_classic()


##form sites now
form_dkogenesup <- filteredlcpm %>% filter((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT > 1)
form_dkogenesdown <- filteredlcpm %>% filter((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT < -1)

a1 <- form_dkogenesdown %>% filter(ensembllistanno %in% naive_wtgenes$ensembllistanno)
length(a1$ensembllistanno) #368
a1$phenotype <- "down" 
a1$state <- "naive" 

a2 <- form_dkogenesdown %>% filter(ensembllistanno %in% form_wtgenes$ensembllistanno)
length(a2$ensembllistanno) #58
a2$phenotype <- "down" 
a2$state <- "form" 

a12 <- rbind(a1,a2) 
a12not <- form_dkogenesdown %>% filter(!ensembllistanno %in% a12$ensembllistanno)
a12not$phenotype <- "down" 
a12not$state <- "other" 


a3 <- form_dkogenesup %>% filter(ensembllistanno %in% naive_wtgenes$ensembllistanno)
length(a3$ensembllistanno) #39
a3$phenotype <- "up" 
a3$state <- "naive" 

a4 <- form_dkogenesup %>% filter(ensembllistanno %in% form_wtgenes$ensembllistanno)
length(a4$ensembllistanno) #275
a4$phenotype <- "up" 
a4$state <- "form" 

a34 <- rbind(a3,a4)
a34not <- form_dkogenesup %>% filter(!ensembllistanno %in% a34$ensembllistanno)
a34not$phenotype <- "up" 
a34not$state <- "other" 


aa <- rbind(a12,a12not,a34,a34not)
aa$phenotype <- factor(aa$phenotype, levels = c("up", "down"))
aa$state <- factor(aa$state, levels = c("naive", "form", "other"))

customcolors <- c("deepskyblue3", "magenta3","grey")
bp2 <- ggplot() +
  geom_bar(data=aa, aes(x=phenotype, fill=state)) +
  ylim(c(0,1000)) +
  scale_fill_manual(values=customcolors) + 
  theme_classic()

bp1 + bp2




naiveform_consistentloss <- naive_dkogenesdown %>% filter(((ensembllistanno %in% k2$ensembllistanno) & NDKO_WT < -1) & ((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT < -1))
naive_dkoloss <- naive_dkogenesdown %>% filter((ensembllistanno %in% naive_wtgenes$ensembllistanno) & FDKO_WT > -1)
naive_dkogenesdown_other <- naive_dkogenesdown %>% filter((!ensembllistanno %in% naiveform_consistentloss$ensembllistanno) & (!ensembllistanno %in% naive_dkoloss$ensembllistanno))

naiveform_consistentloss$designation <- c("consistentloss")
naive_dkoloss$designation <- c("naive_gene_loss")
naive_dkogenesdown_other$designation <- c("other")
x1 <- rbind(naive_dkoloss, naiveform_consistentloss,naive_dkogenesdown_other)
x1$groupname <- "naivedown"

naive_dkogenesup <- filteredlcpm %>% filter((ensembllistanno %in% k2$ensembllistanno) & NDKO_WT > 1)

naive_dkopremature <- naive_dkogenesup %>% filter((ensembllistanno %in% form_wtgenes$ensembllistanno) & NDKO_WT > 1)
naive_dkogenesup_other <- naive_dkogenesup %>% filter((!ensembllistanno %in% naive_dkopremature$ensembllistanno))
naive_dkopremature$designation <- "premature"
naive_dkogenesup_other$designation <- "other"

x2 <- rbind(naive_dkopremature,naive_dkogenesup_other)
x2$groupname <- "naiveup"


x3 <- rbind(x1,x2)
x3$designation <- factor(x3$designation, levels = c("other","naive_gene_loss","premature","consistentloss"))
x3$groupname <- factor(x3$groupname, levels = c("naiveup", "naivedown"))


customcolors <- c("grey", "deepskyblue3", "magenta3", "black")
x4 <- ggplot() +
  geom_bar(data=x3, aes(x=groupname, fill=designation)) +
  ylim(c(0,1000)) +
  scale_fill_manual(values=customcolors) + 
  theme_classic()


form_dkogenesdown <- filteredlcpm %>% filter((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT < -1)
form_dkoloss <- form_dkogenesdown %>% filter((ensembllistanno %in% form_wtgenes$ensembllistanno) & NDKO_WT > -1)
form_dkogenesdown_other <- form_dkogenesdown %>% filter((!ensembllistanno %in% naiveform_consistentloss$ensembllistanno) & (!ensembllistanno %in% form_dkoloss$ensembllistanno))

naiveform_consistentloss <- naive_dkogenesdown %>% filter(((ensembllistanno %in% k2$ensembllistanno) & NDKO_WT < -1) & ((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT < -1))

naiveform_consistentloss$designation <- c("consistentloss")
form_dkoloss$designation <- c("form_gene_loss")
form_dkogenesdown_other$designation <- c("other")
q1 <- rbind(form_dkoloss, naiveform_consistentloss, form_dkogenesdown_other)
q1$groupname <- "formdown"


form_dkogenesup <- filteredlcpm %>% filter((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT > 1)

form_dkonorepress <- naive_wtgenes %>% filter((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT > 1)
form_dkogenesup_other <- form_dkogenesup %>% filter((!ensembllistanno %in% form_dkonorepress$ensembllistanno))
form_dkonorepress$designation <- "not_repressed"
form_dkogenesup_other$designation <- "other"

q2 <- rbind(form_dkonorepress,form_dkogenesup_other )
q2$groupname <- "formup"

q3 <- rbind(q1, q2)
q3$designation <- factor(q3$designation, levels = c("other","form_gene_loss","not_repressed","consistentloss"))
q3$groupname <- factor(q3$groupname, levels = c("formup", "formdown"))


customcolors <- c("grey", "magenta3", "deepskyblue3", "black")
q4 <- ggplot() +
  geom_bar(data=q3, aes(x=groupname, fill=designation)) +
  ylim(c(0,1000)) +
  scale_fill_manual(values=customcolors) +
  theme_classic()

wrap_plots(list(x4, q4))
#####
#####Heatmaps of form genes DKO#####

annotlcpmselect <- subset(filteredlcpm, rownames(filteredlcpm) %in% form_dkogenesdown$ensembllistanno)
#genesofinterestdf  <- as.data.frame(genename.state)
#rownames(genesofinterestdf) <- genesofinterest
#Used resolution of 300x400
x <- annotlcpmselect[c("NWT", "FWT", "NDKO", "FDKO")]
#x <- x[order(FWT),]

#library(circlize)
col_fun2 = colorRamp2(c(-5, 0), c("blue", "white"))
col_fun1 = colorRamp2(c(0, 5), c("white", "red"))
ha <- rowAnnotation(df1 = annotlcpmselect$rnafold_geno_wt, df2 = annotlcpmselect$FDKO_WT,
                    col= list(df1 = col_fun1, df2 = col_fun2))

Heatmap(x, row_km = 5, column_order = order(as.numeric(gsub("column", "", colnames(x)))),
        row_order = order(-x$FWT),
        right_annotation = ha
        )
      

####

annotlcpmselect <- subset(filteredlcpm, rownames(filteredlcpm) %in% form_dkoloss$ensembllistanno)

# annotlcpmselecta <- subset(annotlcpmselect, rownames(annotlcpmselect) %in% form_dkoloss$ensembllistanno)
# annotlcpmselecta$desig <- "formgenes"
# annotlcpmselectb <- subset(annotlcpmselect, !rownames(annotlcpmselect) %in% form_dkoloss$ensembllistanno)
# annotlcpmselectb$desig <- "other"
#annotlcpmselect2 <- rbind(annotlcpmselecta, annotlcpmselectb)
annotlcpmselect <- annotlcpmselect

#genesofinterestdf  <- as.data.frame(genename.state)
#rownames(genesofinterestdf) <- genesofinterest
#Used resolution of 300x400
x <- annotlcpmselect[c("NWT", "FWT", "NDKO", "FDKO")]
#x <- x[order(FWT),]

#library(circlize)
col_fun2 = colorRamp2(c(-5, 0), c("blue", "white"))
col_fun1 = colorRamp2(c(0, 5), c("white", "red"))
# ha <- rowAnnotation(df1 = annotlcpmselect$rnafold_geno_wt, df2 = annotlcpmselect$FDKO_WT, df3 = annotlcpmselect$desig,
#                     col= list(df1 = col_fun1, df2 = col_fun2, df3 = c("formgenes" = "black" , "other" = "grey" )))


ha <- rowAnnotation(df1 = annotlcpmselect$rnafold_geno_wt, df2 = annotlcpmselect$FDKO_WT,
                    col= list(df1 = col_fun1, df2 = col_fun2))

Heatmap(x, row_km = 5, column_order = order(as.numeric(gsub("column", "", colnames(x)))),
        row_order = order(-x$FWT),
        right_annotation = ha,
        show_row_names = FALSE 
)



# kclus <- kmeans(x[1:2], 5)
# kclus$cluster
# split <- paste0("Cluster\n", kclus$cluster)
# Heatmap(x, split=split, cluster_row_slices = FALSE,
#         row_names_gp = grid::gpar(fontsize = 4), 
#         column_order = order(as.numeric(gsub("column", "", colnames(x)))),
#         right_annotation = ha
# )

###naive
annotlcpmselect <- subset(filteredlcpm, rownames(filteredlcpm) %in% naive_dkoloss$ensembllistanno)

# annotlcpmselecta <- subset(annotlcpmselect, rownames(annotlcpmselect) %in% naive_dkoloss$ensembllistanno)
# annotlcpmselecta$desig <- "naivegenes"
# annotlcpmselectb <- subset(annotlcpmselect, !rownames(annotlcpmselect) %in% naive_dkoloss$ensembllistanno)
# annotlcpmselectb$desig <- "other"

#annotlcpmselect2 <- rbind(annotlcpmselecta, annotlcpmselectb)
#annotlcpmselect <- annotlcpmselect2

#genesofinterestdf  <- as.data.frame(genename.state)
#rownames(genesofinterestdf) <- genesofinterest
#Used resolution of 300x400
x <- annotlcpmselect[c("NWT", "FWT", "NDKO", "FDKO")]
#x <- x[order(FWT),]

#library(circlize)
col_fun2 = colorRamp2(c(-5, 0), c("blue", "white"))
col_fun1 = colorRamp2(c(-5, 0), c("blue", "white"))
ha <- rowAnnotation(df1 = annotlcpmselect$rnafold_geno_wt, df2 = annotlcpmselect$NDKO_WT,
                    col= list(df1 = col_fun1, df2 = col_fun2))
Heatmap(x, row_km = 5, column_order = order(as.numeric(gsub("column", "", colnames(x)))),
        row_order = order(-x$NWT),
        right_annotation = ha,
        show_row_names = FALSE
)







#####
#####Upset plots for differential genes#####

library(UpSetR)
#requires R>3.6

library(VennDiagram)
#k2,k5,k7 are ndko,fdko,nfwt respectively





# set_up1 <- k1 %>% filter(log2FoldChange > 1)
# set_up1_names <- set_up1$ensembllistanno
# 

naive_wtgenes <- filteredlcpm %>% filter((ensembllistanno %in% k7$ensembllistanno) & rnafold_geno_wt < -1)
form_wtgenes <- filteredlcpm %>% filter((ensembllistanno %in% k7$ensembllistanno) & rnafold_geno_wt > 1)

naive_dkoloss <- naive_wtgenes %>% filter((ensembllistanno %in% k2$ensembllistanno) & NDKO_WT < -1)
naive_dkopremature <- form_wtgenes %>% filter((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT > 1)

form_dkoloss <- form_wtgenes %>% filter((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT < -1)
form_dkonorepress <- naive_wtgenes %>% filter((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT > 1)

form_wtgenes_nodko <- form_wtgenes %>% filter((!ensembllistanno %in% form_dkoloss$ensembllistanno))


naiveform_constistentloss <- filteredlcpm %>% filter(((ensembllistanno %in% k2$ensembllistanno) & NDKO_WT < -1) & ((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT < -1))

form_wtgenes_sort <- form_wtgenes[order(-form_wtgenes$FWT),]
form_wtgenes_sort$rank <- seq(1,nrow(form_wtgenes_sort))

form_wtgenes_nodko_sort <- form_wtgenes_nodko[order(-form_wtgenes_nodko$FWT),]
form_wtgenes_nodko_sort$rank <- seq(1,nrow(form_wtgenes_nodko_sort))


p1 <- ggplot() +
  geom_point(data=form_wtgenes_nodko_sort, aes(x=rank, y=FDKO), size = 0.2, color="red") +
  geom_point(data=form_wtgenes_nodko_sort, aes(x=rank, y=FWT), size = 0.2) +
  ylim(c(-2, 10)) +
  theme_classic() +
  theme(text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

form_dkoloss_sort <- form_dkoloss[order(-form_dkoloss$FWT),]
form_dkoloss_sort$rank <- seq(1,nrow(form_dkoloss_sort))

p2 <- ggplot() +
  geom_point(data=form_dkoloss_sort, aes(x=rank, y=FDKO), size = 0.1, color="red") +
  geom_point(data=form_dkoloss_sort, aes(x=rank, y=FWT), size = 0.1) +
  theme_classic() +
  ylab("") +
  ylim(c(-2, 10)) +
  theme(text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

wrap_plots(list(p1, p2), widths = c(13, 2))

p3 <- ggplot() +
  geom_point(data=form_wtgenes_nodko_sort, aes(x=rank, y=FCKO), size = 0.2, color="blue") +
  geom_point(data=form_wtgenes_nodko_sort, aes(x=rank, y=FWT), size = 0.2) +
  ylim(c(-2, 10)) +
  theme_classic() +
  theme(text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")



p4 <- ggplot() +
  geom_point(data=form_dkoloss_sort, aes(x=rank, y=FCKO), size = 0.1, color="blue") +
  geom_point(data=form_dkoloss_sort, aes(x=rank, y=FWT), size = 0.1) +
  theme_classic() +
  ylab("") +
  ylim(c(-2, 10)) +
  theme(text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2), legend.position = "none")

wrap_plots(list(p3, p4), widths = c(13, 2))

# set_up1 <- filteredlcpm %>% filter((ensembllistanno %in% k2$ensembllistanno) & NDKO_WT > 1)# %>% filter(NDKO_WT > 1)
# set_up2_names <- set_up1$ensembllistanno
# 
# set_up1 <- filteredlcpm %>% filter((ensembllistanno %in% k2$ensembllistanno) & NDKO_WT > 1)# %>% filter(NDKO_WT > 1)
# set_up2_names <- set_up1$ensembllistanno

# 
# set_up1 <- k3 %>% filter(log2FoldChange > 1)
# set_up3_names <- set_up1$ensembllistanno
# 
# set_up1 <- k4 %>% filter(log2FoldChange > 1)
# set_up4_names <- set_up1$ensembllistanno
# 
# set_up1 <- k5 %>% filter(log2FoldChange > 1)
# set_up5_names <- set_up1$ensembllistanno
# 
# set_up1 <- k6 %>% filter(log2FoldChange > 1)
# set_up6_names <- set_up1$ensembllistanno
# 
# set_up1 <- k7 %>% filter(log2FoldChange > 1)
# set_up7_names <- set_up1$ensembllistanno
# 
# set_up1 <- k1 %>% filter(log2FoldChange < -1)
# set_down1_names <- set_up1$ensembllistanno
# 
# set_up1 <- k2 %>% filter(log2FoldChange < -1)
# set_down2_names <- set_up1$ensembllistanno
# 
# set_up1 <- k3 %>% filter(log2FoldChange < -1)
# set_down3_names <- set_up1$ensembllistanno
# 
# set_up1 <- k4 %>% filter(log2FoldChange < -1)
# set_down4_names <- set_up1$ensembllistanno
# 
# set_up1 <- k5 %>% filter(log2FoldChange < -1)
# set_down5_names <- set_up1$ensembllistanno
# 
# set_up1 <- k6 %>% filter(log2FoldChange < -1)
# set_down6_names <- set_up1$ensembllistanno
# 
# set_up1 <- k7 %>% filter(log2FoldChange < -1)
# set_down7_names <- c(set_up1$ensembllistanno)

listInput <- list(NCKOup = set_up1_names, NDKOup = set_up2_names, NdCDup = set_up3_names, 
                  FCKOup = set_up4_names, FDKOup = set_up5_names, FdCDup = set_up6_names, Formative = set_up7_names,
                  NCKOdown = set_down1_names, NDKOdown = set_down2_names, NdCDdown = set_down3_names, 
                  FCKOdown = set_down4_names, FDKOdown = set_down5_names, FdCDdown= set_down6_names, Naive = set_down7_names)

#upset(fromList(listInput), nsets = 14, group.by = "sets")
upset(fromList(listInput), order.by = "freq", nsets = 14, nintersects = 25)

listInput <- list(NDKOup = set_up2_names, FDKOup = set_up5_names, Formative = set_up7_names, NDKOdown = set_down2_names, 
                  FDKOdown = set_down5_names, Naive = set_down7_names)

listInput <- list(NDKOup = set_up2_names, FDKOup = set_up5_names, Formative = set_up7_names, NDKOdown = set_down2_names, 
                  FDKOdown = set_down5_names, Naive = set_down7_names, NCKOdown = set_down1_names, FCKOdown = set_down4_names,
                  NCKOup = set_up1_names, FCKOup = set_up4_names)

listInput <- list(Formative = set_up7_names,
                 Naive = set_down7_names, NCKOdown = set_down1_names, FCKOdown = set_down4_names,
                  NCKOup = set_up1_names, FCKOup = set_up4_names)

listInput <- list(NDKOup = set_up2_names, FDKOup = set_up5_names, NDKOdown = set_down2_names, FDKOdown = set_down5_names,
                  NCKOdown = set_down1_names, FCKOdown = set_down4_names,
                  NCKOup = set_up1_names, FCKOup = set_up4_names)


listInput2 <- fromList(listInput)

upset(fromList(listInput), order.by = "freq", nsets = 10, nintersects = 15)

##use calculate.overlap from VennDiagram library to extract overlapping names, 3rd object of assigned object list
geneset_premature <- calculate.overlap(list(listInput$NDKOup, listInput$Formative))
geneset_notactivated <- calculate.overlap(list(listInput$FDKOdown, listInput$Formative))
geneset_consistentectopic <- calculate.overlap(list(listInput$NDKOup, listInput$FDKOup))
geneset_consistentloss <- calculate.overlap(list(listInput$NDKOdown, listInput$FDKOdown))
geneset_notrepressed <- calculate.overlap(list(listInput$FDKOup, listInput$Naive))
geneset_naiveloss <- calculate.overlap(list(listInput$NDKOdown, listInput$Naive))

# geneset_allbyall <- calculate.overlap(list(listInput$FDKOdown, listInput$FDKOup, listInput$NDKOup, listInput$NDKOdown, listInput$Naive, listInput$Formative))                
# threeway comparison order "a123", "a12", "a13", "a23", "a1", "a2", "a3")
geneset_premature <- calculate.overlap(list(listInput$NDKOup, listInput$Formative, listInput$FDKOup))
geneset_notactivated <- calculate.overlap(list(listInput$FDKOdown, listInput$Formative, listInput$NDKOdown))
geneset_consistentectopic <- calculate.overlap(list(listInput$NDKOup, listInput$FDKOup))
geneset_consistentloss <- calculate.overlap(list(listInput$NDKOdown, listInput$FDKOdown))
geneset_notrepressed <- calculate.overlap(list(listInput$FDKOup, listInput$Naive, listInput$NDKOup))

geneset_premature <- geneset_premature[[2]]
geneset_notactivated <- geneset_notactivated[[2]]
geneset_consistentectopic <- geneset_consistentectopic[[3]]
geneset_consistentloss <- geneset_consistentloss[[3]]
geneset_notrepressed <- geneset_notrepressed[[2]]




#####
