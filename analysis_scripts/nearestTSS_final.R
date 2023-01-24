#Code used to generate all nearest neighbor analysis data. Draws from dataframes created with rnaseq and k4m1/k27ac related scripts
#Also used to generate enhancer compensation figures in Fig.5

library(grid)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)

#####sessionInfo#####

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
#   [1] ggpubr_0.4.0    patchwork_1.1.1 forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4     readr_2.1.2     tidyr_1.2.0    
# [9] tibble_3.1.7    tidyverse_1.3.1 ggplot2_3.3.6  
# 
# loaded via a namespace (and not attached):
#   [1] cellranger_1.1.0 pillar_1.7.0     compiler_4.1.2   dbplyr_2.2.1     tools_4.1.2      lubridate_1.8.0  jsonlite_1.8.0   lifecycle_1.0.1 
# [9] gtable_0.3.0     pkgconfig_2.0.3  rlang_1.0.3      reprex_2.0.1     DBI_1.1.3        cli_3.3.0        rstudioapi_0.13  haven_2.5.0     
# [17] xml2_1.3.3       withr_2.5.0      httr_1.4.3       fs_1.5.2         generics_0.1.3   vctrs_0.4.1      hms_1.1.1        tidyselect_1.1.2
# [25] glue_1.6.2       R6_2.5.1         rstatix_0.7.0    fansi_1.0.3      readxl_1.4.0     carData_3.0-5    car_3.1-0        tzdb_0.3.0      
# [33] modelr_0.1.8     magrittr_2.0.3   scales_1.2.0     backports_1.4.1  ellipsis_0.3.2   abind_1.4-5      rvest_1.0.2      assertthat_0.2.1
# [41] colorspace_2.0-3 ggsignif_0.6.3   utf8_1.2.2       stringi_1.7.6    munsell_0.5.0    broom_1.0.0      crayon_1.5.1


#####


#####
#####Import original lcpm dataframe with TMM using all samples#####

##all genes first in geno types
#original lcpm df using all samples to TMM
#annotlcpm <- read.csv("~/Desktop/data/Diffgenes/avedf_lcpmv2.anno.csv", header =T)
# annotlcpm <- subset(annotlcpm, (!is.na(annotlcpm[,10])))
# annotlcpm <- distinct(annotlcpm, annotlcpm$ensembllistanno, .keep_all = TRUE)
# x <- annotlcpm[,c("NWT", "FWT", "NDKO", "FDKO")]
# x2 <- gather(x)
# x2$key <- factor(x2$key, levels = c("NWT", "FWT", "NDKO", "FDKO"))
#####
#####Required:Import new lcpm dataframe with TMM on only WT and DKO#####
##performing TMM with just WT and DKO samples
#df_lcpmannot_WTDKOnorm_20220427.txt
annotlcpm <- read.table("~/Desktop/data/processed_data/df_lcpmannot_WTDKOnorm_20220427.txt", header =T)
NWT <- rowMeans(annotlcpm[1:3])
FWT <- rowMeans(annotlcpm[4:6])
NDKO <- rowMeans(annotlcpm[7:9])
FDKO <- rowMeans(annotlcpm[10:12])
rna_fold_wt2 <- FWT - NWT
rna_fold_dko2 <- FDKO - NDKO
rna_fold_nwt_dko2 <- NDKO - NWT
rna_fold_fwt_dko2 <- FDKO - FWT

ensembllistanno <- annotlcpm$ensembllistanno

annotlcpm <- data.frame(ensembllistanno, NWT, FWT, NDKO, FDKO, rna_fold_wt2, rna_fold_dko2,rna_fold_nwt_dko2, rna_fold_fwt_dko2)
annotlcpm <- subset(annotlcpm, (!is.na(annotlcpm[,1])))
annotlcpm <- distinct(annotlcpm, annotlcpm$ensembllistanno, .keep_all = TRUE)

colnames(annotlcpm) <- c("neargene.name", "NWT2", "FWT2", "NDKO2", "FDKO2", "rna_fold_wt2", "rna_fold_dko2","rna_fold_nwt_dko2", "rna_fold_fwt_dko2")
annotlcpm <- annotlcpm[1:9]
filteredlcpm <- annotlcpm

#####
#####generate group categories for ntss forloop#####


##specify k4m1 alone groups
supermatrix.k4m1 <- read.delim("~/Desktop/data/processed_data/supermatrix_k4m1.20220426.txt", sep = '\t', header = TRUE)

#all k4m1 subset by atac, not for nearest tss based on figure progression :-( 
#k4m1_atac_filter <- supermatrix.k4m1 
k4m1_atac_filter <- read.table("~/Desktop/k4m1wt_atacfilter.bed", sep = '\t')
colnames(k4m1_atac_filter) <- c("chr", "start", "end", "width")
k4m1_atac_filter$start <- k4m1_atac_filter$start + 1
supermatrix.k4m1 <- inner_join(supermatrix.k4m1, k4m1_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

supermatrix.k4m1.up <- supermatrix.k4m1 %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.k4m1.down <- supermatrix.k4m1 %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

supermatrix.k4m1.non  <- supermatrix.k4m1 %>% filter(FDR > 0.1 & pk_fold_nf_wt_k4m1 > -0.7 & pk_fold_nf_wt_k4m1 < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.k4m1.non.loss <- supermatrix.k4m1.non %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1)
supermatrix.k4m1.non.non <- supermatrix.k4m1.non %>% filter(pk_fold_f_dkowt_k4m1 > -0.7 & pk_fold_n_dkowt_k4m1 > -0.7)


#k4m1 only, no k27a 
k27a_atac_filter <- read.table("~/Desktop/data/k4m1k27a_k4m1peaksnotk27a_at.bed", sep = '\t')
colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
supermatrix.k4m1only <- inner_join(supermatrix.k4m1, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

supermatrix.k4m1only.up <- supermatrix.k4m1only %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.k4m1only.down <- supermatrix.k4m1only %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

supermatrix.k4m1only.non  <- supermatrix.k4m1only %>% filter(FDR > 0.1 & pk_fold_nf_wt_k4m1 > -0.7 & pk_fold_nf_wt_k4m1 < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.k4m1only.non.loss <- supermatrix.k4m1only.non %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1)
supermatrix.k4m1only.non.non <- supermatrix.k4m1only.non %>% filter(pk_fold_f_dkowt_k4m1 > -0.7 & pk_fold_n_dkowt_k4m1 > -0.7)

#nonloss 2458
#nonnon 18135



##specify k27a alone groups
supermatrix.k27a <- read.delim("~/Desktop/data/processed_data/supermatrix_k27a.20220427.txt", sep = '\t', header = TRUE)

k27a_atac_filter <- read.table("~/Desktop/k27a_at_peaksfiltered.bed", sep = '\t')
colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
supermatrix.k27a <- inner_join(supermatrix.k27a, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

supermatrix.k27a.up <- supermatrix.k27a %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.k27a.down <- supermatrix.k27a %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.k27a.non  <- supermatrix.k27a %>% filter(FDR > 0.1 & pk_fold_nf_wt_k27a > -0.7 & pk_fold_nf_wt_k27a < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

supermatrix.up_k27a_dep <- supermatrix.k27a.up %>% filter(pk_fold_f_dkowt_k27a < -1)
supermatrix.up_k27a_ind <- supermatrix.k27a.up %>% filter(pk_fold_f_dkowt_k27a < 0.7 & pk_fold_f_dkowt_k27a > -0.7 )
#write.table(supermatrix.up_k27a_dep[1:3], "~/Desktop/xx.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
supermatrix.down_k27a_dep <- supermatrix.k27a.down %>% filter(pk_fold_n_dkowt_k27a < -1)
supermatrix.down_k27a_ind <- supermatrix.k27a.down %>% filter(pk_fold_n_dkowt_k27a < 0.7 & pk_fold_n_dkowt_k27a > -0.7 )

supermatrix.non_k27a_dep <- supermatrix.k27a.non %>% filter(pk_fold_n_dkowt_k27a < -1 & pk_fold_f_dkowt_k27a < -1)
supermatrix.non_k27a_ind <- supermatrix.k27a.non %>% filter((pk_fold_n_dkowt_k27a < 0.7 & pk_fold_n_dkowt_k27a > -0.7) | (pk_fold_f_dkowt_k27a < 0.7 & pk_fold_f_dkowt_k27a > -0.7))



##specify k4m1/k27a groups

k27a_atac_filter <- read.table("~/Desktop/data/k4m1k27a_k27apeaksbyk4m1at.bed", sep = '\t')
colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
supermatrix.k4m1k27a <- inner_join(supermatrix.k27a, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))
supermatrix.ii.k4m1k27a <- supermatrix.k4m1k27a %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

##
x <- read.table("~/Desktop/data/k4m1k27asubset/y1.bed") 
colnames(x) <- c("chr", "start", "end")
x <- inner_join(supermatrix.k27a, x, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))


# y1ind <- x %>% filter((pk_fold_n_dkowt_k27a < 0.7 & pk_fold_n_dkowt_k27a > -0.7) | (pk_fold_f_dkowt_k27a < 0.7 & pk_fold_f_dkowt_k27a > -0.7))
# y1dep <- x %>% filter(pk_fold_n_dkowt_k27a < -1 & pk_fold_f_dkowt_k27a < -1)

y1ind <- x %>% filter(pk_fold_f_dkowt_k27a > -1) #910
y1dep <- x %>% filter(pk_fold_f_dkowt_k27a < -1) ##1005

# write.table(y1ind[1:3], "~/Desktop/y1ind.bed", col.names = FALSE, quote = FALSE, row.names = FALSE)
# write.table(y1dep[1:3], "~/Desktop/y1dep.bed", col.names = FALSE, quote = FALSE, row.names = FALSE)


##
x <- read.table("~/Desktop/data/k4m1k27asubset/z3.bed") 
colnames(x) <- c("chr", "start", "end")
x <- inner_join(supermatrix.k27a, x, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

z3ind <- x %>% filter(pk_fold_n_dkowt_k27a > -1)
z3dep <- x %>% filter(pk_fold_n_dkowt_k27a < -1)

write.table(z3ind[1:3], "~/Desktop/z3ind.bed", col.names = FALSE, quote = FALSE, row.names = FALSE)
write.table(z3dep[1:3], "~/Desktop/z3dep.bed", col.names = FALSE, quote = FALSE, row.names = FALSE)


x <- read.table("~/Desktop/data/k4m1k27asubset/z1.bed")  #2164 sites
colnames(x) <- c("chr", "start", "end")
x <- inner_join(supermatrix.k27a, x, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

z1ind <- x %>% filter(pk_fold_n_dkowt_k27a > -1) #1320
z1dep <- x %>% filter(pk_fold_n_dkowt_k27a < -1) #844

write.table(z1ind[1:3], "~/Desktop/z1ind.bed", col.names = FALSE, quote = FALSE, row.names = FALSE)
write.table(z1dep[1:3], "~/Desktop/z1dep.bed", col.names = FALSE, quote = FALSE, row.names = FALSE)



generate_deepheatmapgrouped(list(z1ind,z1dep), "~/Desktop/zz1.bed")
##
x <- read.table("~/Desktop/data/k4m1k27asubset/y4.bed") 
colnames(x) <- c("chr", "start", "end")
x <- inner_join(supermatrix.k27a, x, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

y4ind <- x %>% filter(pk_fold_f_dkowt_k27a > -1)
y4dep <- x %>% filter(pk_fold_f_dkowt_k27a < -1)

# write.table(y4ind[1:3], "~/Desktop/y4ind.bed", col.names = FALSE, quote = FALSE, row.names = FALSE)
# write.table(y4dep[1:3], "~/Desktop/y4dep.bed", col.names = FALSE, quote = FALSE, row.names = FALSE)

x <- read.table("~/Desktop/data/k4m1k27asubset/x2.bed") 
colnames(x) <- c("chr", "start", "end")
x <- inner_join(supermatrix.k27a, x, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

x2ind <- x %>% filter(pk_fold_f_dkowt_k27a > -1)
x2dep <- x %>% filter(pk_fold_f_dkowt_k27a < -1)

x <- read.table("~/Desktop/data/k4m1k27asubset/y2.bed") 
colnames(x) <- c("chr", "start", "end")
x <- inner_join(supermatrix.k27a, x, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

y2ind <- x %>% filter(pk_fold_f_dkowt_k27a > -1)
y2dep <- x %>% filter(pk_fold_f_dkowt_k27a < -1)

x <- read.table("~/Desktop/data/k4m1k27asubset/z2.bed") 
colnames(x) <- c("chr", "start", "end")
x <- inner_join(supermatrix.k27a, x, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

z2ind <- x %>% filter(pk_fold_n_dkowt_k27a > -1) #168
z2dep <- x %>% filter(pk_fold_n_dkowt_k27a < -1) #973

generate_deepheatmapgrouped(list(x2ind,x2dep,z2ind,z2dep,y2ind,y2dep ), "~/Desktop/xzy2.bed")



allgenes <- annotlcpm
allk4m1peaks <- supermatrix.k4m1
allk27apeaks <- supermatrix.k27a
allk4m1k27apeaks <- supermatrix.ii.k4m1k27a

#####
#####lcpm all genes control#####

rnalist <- list(annotlcpm)
##save as 4x2" pdf per panel
plotlist=list()
statslist=list()
meanlist=list()
counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  compare1 <- c("NWT2", "FWT2", "NDKO2", "FDKO2")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather$key <- factor(matrix_oi_gather$key,levels = c("NWT2", "FWT2", "NDKO2", "FDKO2"))
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  statslist[[counter]] <- compare_means(value ~ key, data = matrix_oi_gather2, method="wilcox.test", p.adjust.method="BH" )
  
  meanlist[[counter]] <- matrix_oi_gather2 %>% group_by(key) %>% summarise_at(vars(value), list(name = mean)) 
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 geom_violin(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 stat_summary(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red") +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-3,17) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)




#####
#####lcpm forloop for peak categories#####

rnalist <- list(supermatrix.k4m1.non.non, supermatrix.k4m1.non.loss, supermatrix.k4m1.up, supermatrix.k4m1.down)
rnalist <- list(supermatrix.k27a.non, supermatrix.k27a.down, supermatrix.k27a.up)
rnalist <- list(supermatrix.non_k27a_ind, supermatrix.non_k27a_dep, supermatrix.down_k27a_ind, supermatrix.down_k27a_dep, supermatrix.up_k27a_ind, supermatrix.up_k27a_dep)
rnalist <- list(z3ind, z3dep)
rnalist <- list(y4ind, y4dep)

##k4m1 alone peaks

#saved as 4x8in
rnalist <- list(supermatrix.k4m1only.non.non, supermatrix.k4m1only.non.loss, supermatrix.k4m1only.down, supermatrix.k4m1only.up)


rnalist <- list(supermatrix.k4m1.non.non, supermatrix.k4m1.non.loss, supermatrix.k4m1.down, supermatrix.k4m1.up)
tablename <- c("k4m1alone")
statsnames <- c("k4m1nonnon", "k4m1nonloss", "k4m1down", "k4m1up")

plotlist=list()
statslist=list()
meanlist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  matrix_oi <- inner_join(matrix_oi, annotlcpm, by = c("neargene.name" = "neargene.name"))
  
  
  compare1 <- c("NWT2", "FWT2", "NDKO2", "FDKO2")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather$key <- factor(matrix_oi_gather$key,levels = c("NWT2", "FWT2", "NDKO2", "FDKO2"))
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  statslist[[counter]] <- compare_means(value ~ key, data = matrix_oi_gather2, method="wilcox.test", p.adjust.method="BH" )
  
  meanlist[[counter]] <- matrix_oi_gather2 %>% group_by(key) %>% summarise_at(vars(value), list(name = mean)) 
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-3,17) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)

counter = 0
stats_grouped <- c() 
statsarray <- for(i in 1:length(statslist)) {
  counter = counter + 1
  stats_oi <- statslist[[i]]
  stats_oi$group <- statsnames[[counter]]
  
  stats_grouped <- rbind(stats_grouped, stats_oi)
  
}

counter = 0
means_grouped <- c()
meanarray <- for(i in 1:length(meanlist)) {
  counter = counter + 1
  means_oi <- meanlist[[i]]
  means_oi$group <- statsnames[[counter]]
  
  means_grouped <- rbind(means_grouped, means_oi)
  
}

##k27a alone peaks
rnalist <- list(supermatrix.non_k27a_ind, supermatrix.non_k27a_dep, supermatrix.down_k27a_ind, supermatrix.down_k27a_dep, supermatrix.up_k27a_ind, supermatrix.up_k27a_dep)
tablename <- c("k27aalone")
statsnames <- c("k27anonind", "k27anondep", "k27adownind", "k27adowndep", "k27aupind", "k27aupdep")

#saved as 6x8in
plotlist=list()
statslist=list()
meanlist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  matrix_oi <- inner_join(matrix_oi, annotlcpm, by = c("neargene.name" = "neargene.name"))
  
  
  compare1 <- c("NWT2", "FWT2", "NDKO2", "FDKO2")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather$key <- factor(matrix_oi_gather$key,levels = c("NWT2", "FWT2", "NDKO2", "FDKO2"))
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  statslist[[counter]] <- compare_means(value ~ key, data = matrix_oi_gather2, method="wilcox.test", p.adjust.method="BH" )
  
  meanlist[[counter]] <- matrix_oi_gather2 %>% group_by(key) %>% summarise_at(vars(value), list(name = mean)) 
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-2,15) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)

counter = 0
stats_grouped <- c() 
statsarray <- for(i in 1:length(statslist)) {
  counter = counter + 1
  stats_oi <- statslist[[i]]
  stats_oi$group <- statsnames[[counter]]
  
  stats_grouped <- rbind(stats_grouped, stats_oi)
  
}

counter = 0
means_grouped <- c()
meanarray <- for(i in 1:length(meanlist)) {
  counter = counter + 1
  means_oi <- meanlist[[i]]
  means_oi$group <- statsnames[[counter]]
  
  means_grouped <- rbind(means_grouped, means_oi)
  
}



##k4m1/k27a form peaks

rnalist <- list(z3ind, z3dep)
tablename <- c("z3")
statsnames <- c("z3ind", "z3dep")

plotlist=list()
statslist=list()
meanlist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  matrix_oi <- inner_join(matrix_oi, annotlcpm, by = c("neargene.name" = "neargene.name"))
  
  
  compare1 <- c("NWT2", "FWT2", "NDKO2", "FDKO2")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather$key <- factor(matrix_oi_gather$key,levels = c("NWT2", "FWT2", "NDKO2", "FDKO2"))
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  statslist[[counter]] <- compare_means(value ~ key, data = matrix_oi_gather2, method="wilcox.test", p.adjust.method="BH" )
  
  meanlist[[counter]] <- matrix_oi_gather2 %>% group_by(key) %>% summarise_at(vars(value), list(name = mean)) 
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-3,17) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)

counter = 0
stats_grouped <- c() 
statsarray <- for(i in 1:length(statslist)) {
  counter = counter + 1
  stats_oi <- statslist[[i]]
  stats_oi$group <- statsnames[[counter]]
  
  stats_grouped <- rbind(stats_grouped, stats_oi)
  
}

counter = 0
means_grouped <- c()
meanarray <- for(i in 1:length(meanlist)) {
  counter = counter + 1
  means_oi <- meanlist[[i]]
  means_oi$group <- statsnames[[counter]]
  
  means_grouped <- rbind(means_grouped, means_oi)
  
}


##save as 2x4" pdf per panel
rnalist <- list(y4ind, y4dep)
tablename <- c("y4")
statsnames <- c("y4ind", "y4dep")

plotlist=list()
statslist=list()
meanlist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  matrix_oi <- inner_join(matrix_oi, annotlcpm, by = c("neargene.name" = "neargene.name"))
  
  
  compare1 <- c("NWT2", "FWT2", "NDKO2", "FDKO2")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather$key <- factor(matrix_oi_gather$key,levels = c("NWT2", "FWT2", "NDKO2", "FDKO2"))
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  statslist[[counter]] <- compare_means(value ~ key, data = matrix_oi_gather2, method="wilcox.test", p.adjust.method="BH" )
  
  meanlist[[counter]] <- matrix_oi_gather2 %>% group_by(key) %>% summarise_at(vars(value), list(name = mean)) 
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-3,17) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)

counter = 0
stats_grouped <- c() 
statsarray <- for(i in 1:length(statslist)) {
  counter = counter + 1
  stats_oi <- statslist[[i]]
  stats_oi$group <- statsnames[[counter]]
  
  stats_grouped <- rbind(stats_grouped, stats_oi)
  
  }

counter = 0
means_grouped <- c()
meanarray <- for(i in 1:length(meanlist)) {
  counter = counter + 1
  means_oi <- meanlist[[i]]
  means_oi$group <- statsnames[[counter]]
  
  means_grouped <- rbind(means_grouped, means_oi)
  
}

##K4me1 subset by atac
rnalist <- list(supermatrix.k4m1non.non, supermatrix.k4m1non.loss, supermatrix.down, supermatrix.up)
tablename <- c("k4m1alone")
statsnames <- c("k4m1nonnon", "k4m1nonloss", "k4m1down", "k4m1up")

plotlist=list()
statslist=list()
meanlist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  matrix_oi <- inner_join(matrix_oi, annotlcpm, by = c("neargene.name" = "neargene.name"))
  
  
  compare1 <- c("NWT2", "FWT2", "NDKO2", "FDKO2")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather$key <- factor(matrix_oi_gather$key,levels = c("NWT2", "FWT2", "NDKO2", "FDKO2"))
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  statslist[[counter]] <- compare_means(value ~ key, data = matrix_oi_gather2, method="wilcox.test", p.adjust.method="BH" )
  
  meanlist[[counter]] <- matrix_oi_gather2 %>% group_by(key) %>% summarise_at(vars(value), list(name = mean)) 
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-3,17) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)

counter = 0
stats_grouped <- c() 
statsarray <- for(i in 1:length(statslist)) {
  counter = counter + 1
  stats_oi <- statslist[[i]]
  stats_oi$group <- statsnames[[counter]]
  
  stats_grouped <- rbind(stats_grouped, stats_oi)
  
}

counter = 0
means_grouped <- c()
meanarray <- for(i in 1:length(meanlist)) {
  counter = counter + 1
  means_oi <- meanlist[[i]]
  means_oi$group <- statsnames[[counter]]
  
  means_grouped <- rbind(means_grouped, means_oi)
  
}




#####

#####fold change forloop for peak categories#####

##all phenos
rnalist <- list(annotlcpm)

plotlist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  compare1 <- c("rna_fold_nwt_dko2", "rna_fold_fwt_dko2")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-4,4) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)




##k4m1/k27a naive phenos
rnalist <- list(z3ind, z3dep)
tablename <- c("z3")
statsnames <- c("z3ind", "z3dep")

plotlist=list()
statslist=list()
meanlist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  matrix_oi <- inner_join(matrix_oi, annotlcpm, by = c("neargene.name" = "neargene.name"))
  
  #compare1 <- c("rna_fold_nwt_dko2", "rna_fold_fwt_dko2")
  compare1 <- c("rna_fold_nwt_dko", "rna_fold_fwt_dko")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-4,4) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)

##formative 
rnalist <- list(y4ind, y4dep)
tablename <- c("y4")
statsnames <- c("y4ind", "y4dep")

plotlist=list()
statslist=list()
meanlist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  matrix_oi <- inner_join(matrix_oi, annotlcpm, by = c("neargene.name" = "neargene.name"))
  
  #compare1 <- c("rna_fold_nwt_dko2", "rna_fold_fwt_dko2")
  compare1 <- c("rna_fold_nwt_dko", "rna_fold_fwt_dko")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-4,4) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)


###naive all k27a
rnalist <- list(supermatrix.down_k27a_ind, supermatrix.down_k27a_dep)
tablename <- c("naivek27a")
statsnames <- c("k27aind", "k27adep")

plotlist=list()
statslist=list()
meanlist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  
  matrix_oi <- inner_join(matrix_oi, annotlcpm, by = c("neargene.name" = "neargene.name"))
  
  
  compare1 <- c("rna_fold_nwt_dko2", "rna_fold_fwt_dko2")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  
  matrix_oi2 <- matrix_oi %>% filter(!is.na(NWT2))
  #count total peaks used
  countx <- count(matrix_oi2)
  countx <- paste("n", "=", countx, sep=' ')
  grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  #count unique genes used
  county <- length(unique(matrix_oi2$neargene.name))
  county <- paste("genes", "=", county, sep=' ')
  groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                             gp=gpar(col="black", fontsize=10)))
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
                                 stat_boxplot(data=matrix_oi_gather2, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
                                 geom_boxplot(data=matrix_oi_gather2, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
                                 scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
                                 ylim(-4,4) +
                                 ylab("") +
                                 xlab("") +
                                 annotation_custom(grobx) +
                                 annotation_custom(groby) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
                               
  )
  
}
wrap_plots(plotlist)




#####
#####code?#####
x1z1 <- read.table("~/Desktop/data/k4m1k27asubset/k27apeak_k4m1sub_x1z1.bed") 
colnames(x1z1) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
x1z1 <- inner_join(supermatrix.k27a, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))
#####
#####enhancer compensation mainfigs#####
##naive
matrix_oi <- supermatrix.ii.k4m1k27a
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a > -1))
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y2 <- y[order(-y$totalpeaks, -y$prop_lost),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

q0 <- y3 %>% filter(prop_lost == 1)
q0$genepeakpheno <- "q0"
q1 <- y3 %>% filter(prop_lost >= 0.75 & prop_lost < 1)
q1$genepeakpheno <- "q1"
q2 <- y3 %>% filter(prop_lost >= 0.50 & prop_lost < 0.75)
q2$genepeakpheno <- "q2"
q3 <- y3 %>% filter(prop_lost >= 0.25 & prop_lost < 0.50)
q3$genepeakpheno <- "q3"
q4 <- y3 %>% filter(prop_lost < 0.25 & prop_lost > 0)
q4$genepeakpheno <- "q4"
q5 <- y3 %>% filter(prop_lost == 0)
q5$genepeakpheno <- "q5"

y4 <- rbind(q0, q1, q2, q3, q4, q5)
y4v2 <- inner_join(y4, annotlcpm, by = c("gene.name" = "neargene.name"))

comp1 <- ggplot(data=y4v2, aes(x=genepeakpheno, y=rna_fold_nwt_dko)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(label=..count..), y=2, stat='count', colour="black", size=4) +
  theme_classic() +
  ylim(-2.5,2.5) +
  geom_hline(yintercept=0)


##formative
matrix_oi <- supermatrix.ii.k4m1k27a
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y2 <- y[order(-y$totalpeaks, -y$prop_lost),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

##create quintiles for enhancer proportion lost
q0 <- y3 %>% filter(prop_lost == 1)
q0$genepeakpheno <- "q0"
q1 <- y3 %>% filter(prop_lost >= 0.75 & prop_lost < 1)
q1$genepeakpheno <- "q1"
q2 <- y3 %>% filter(prop_lost >= 0.50 & prop_lost < 0.75)
q2$genepeakpheno <- "q2"
q3 <- y3 %>% filter(prop_lost >= 0.25 & prop_lost < 0.50)
q3$genepeakpheno <- "q3"
q4 <- y3 %>% filter(prop_lost < 0.25 & prop_lost > 0)
q4$genepeakpheno <- "q4"
q5 <- y3 %>% filter(prop_lost == 0)
q5$genepeakpheno <- "q5"

y4 <- rbind(q0, q1, q2, q3, q4, q5)
y4v2 <- inner_join(y4, annotlcpm, by = c("gene.name" = "neargene.name"))

comp2 <- ggplot(data=y4v2, aes(x=genepeakpheno, y=rna_fold_fwt_dko)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(label=..count..), y=2, stat='count', colour="black", size=4) +
  theme_classic() +
  ylim(-2.5,2.5) +
  geom_hline(yintercept=0)

wrap_plots(comp1 + comp2)


#####
#####enhancer compensation, greater loss mod#####
##naive
matrix_oi <- supermatrix.ii.k4m1k27a
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a < -2))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a > -2))
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y2 <- y[order(-y$totalpeaks, -y$prop_lost),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

q0 <- y3 %>% filter(prop_lost == 1)
q0$genepeakpheno <- "q0"
q1 <- y3 %>% filter(prop_lost >= 0.75 & prop_lost < 1)
q1$genepeakpheno <- "q1"
q2 <- y3 %>% filter(prop_lost >= 0.50 & prop_lost < 0.75)
q2$genepeakpheno <- "q2"
q3 <- y3 %>% filter(prop_lost >= 0.25 & prop_lost < 0.50)
q3$genepeakpheno <- "q3"
q4 <- y3 %>% filter(prop_lost < 0.25 & prop_lost > 0)
q4$genepeakpheno <- "q4"
q5 <- y3 %>% filter(prop_lost == 0)
q5$genepeakpheno <- "q5"

y4 <- rbind(q0, q1, q2, q3, q4, q5)
y4v2 <- inner_join(y4, annotlcpm, by = c("gene.name" = "neargene.name"))

comp1 <- ggplot(data=y4v2, aes(x=genepeakpheno, y=rna_fold_nwt_dko)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(label=..count..), y=2, stat='count', colour="black", size=4) +
  theme_classic() +
  ylim(-2.5,2.5) +
  geom_hline(yintercept=0)

##formative
matrix_oi <- supermatrix.ii.k4m1k27a
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a < -2))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a > -2))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y2 <- y[order(-y$totalpeaks, -y$prop_lost),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

##create quintiles for enhancer proportion lost
q0 <- y3 %>% filter(prop_lost == 1)
q0$genepeakpheno <- "q0"
q1 <- y3 %>% filter(prop_lost >= 0.75 & prop_lost < 1)
q1$genepeakpheno <- "q1"
q2 <- y3 %>% filter(prop_lost >= 0.50 & prop_lost < 0.75)
q2$genepeakpheno <- "q2"
q3 <- y3 %>% filter(prop_lost >= 0.25 & prop_lost < 0.50)
q3$genepeakpheno <- "q3"
q4 <- y3 %>% filter(prop_lost < 0.25 & prop_lost > 0)
q4$genepeakpheno <- "q4"
q5 <- y3 %>% filter(prop_lost == 0)
q5$genepeakpheno <- "q5"

y4 <- rbind(q0, q1, q2, q3, q4, q5)
y4v2 <- inner_join(y4, annotlcpm, by = c("gene.name" = "neargene.name"))

comp2 <- ggplot(data=y4v2, aes(x=genepeakpheno, y=rna_fold_fwt_dko)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(label=..count..), y=2, stat='count', colour="black", size=4) +
  theme_classic() +
  ylim(-2.5,2.5) +
  geom_hline(yintercept=0)

wrap_plots(comp1 + comp2)

#####
#####enhancer compensation, distance to tss max 10kb #####
##naive
matrix_oi <- supermatrix.ii.k4m1k27a
matrix_oi$Distance.to.TSS <- abs(matrix_oi$Distance.to.TSS)
matrix_oi <- matrix_oi %>% filter(Distance.to.TSS < 10000)

x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a > -1))
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y2 <- y[order(-y$totalpeaks, -y$prop_lost),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

q0 <- y3 %>% filter(prop_lost == 1)
q0$genepeakpheno <- "q0"
q1 <- y3 %>% filter(prop_lost >= 0.75 & prop_lost < 1)
q1$genepeakpheno <- "q1"
q2 <- y3 %>% filter(prop_lost >= 0.50 & prop_lost < 0.75)
q2$genepeakpheno <- "q2"
q3 <- y3 %>% filter(prop_lost >= 0.25 & prop_lost < 0.50)
q3$genepeakpheno <- "q3"
q4 <- y3 %>% filter(prop_lost < 0.25 & prop_lost > 0)
q4$genepeakpheno <- "q4"
q5 <- y3 %>% filter(prop_lost == 0)
q5$genepeakpheno <- "q5"

y4 <- rbind(q0, q1, q2, q3, q4, q5)
y4v2 <- inner_join(y4, annotlcpm, by = c("gene.name" = "neargene.name"))

comp3 <- ggplot(data=y4v2, aes(x=genepeakpheno, y=rna_fold_nwt_dko)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(label=..count..), y=2, stat='count', colour="black", size=4) +
  theme_classic() +
  ylim(-2.5,2.5) +
  geom_hline(yintercept=0)


##formative
matrix_oi <- supermatrix.ii.k4m1k27a
matrix_oi$Distance.to.TSS <- abs(matrix_oi$Distance.to.TSS)
matrix_oi <- matrix_oi %>% filter(Distance.to.TSS < 10000)
x <- matrix_oi %>% filter(!is.na(NWT))

##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y2 <- y[order(-y$totalpeaks, -y$prop_lost),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

##create quintiles for enhancer proportion lost
q0 <- y3 %>% filter(prop_lost == 1)
q0$genepeakpheno <- "q0"
q1 <- y3 %>% filter(prop_lost >= 0.75 & prop_lost < 1)
q1$genepeakpheno <- "q1"
q2 <- y3 %>% filter(prop_lost >= 0.50 & prop_lost < 0.75)
q2$genepeakpheno <- "q2"
q3 <- y3 %>% filter(prop_lost >= 0.25 & prop_lost < 0.50)
q3$genepeakpheno <- "q3"
q4 <- y3 %>% filter(prop_lost < 0.25 & prop_lost > 0)
q4$genepeakpheno <- "q4"
q5 <- y3 %>% filter(prop_lost == 0)
q5$genepeakpheno <- "q5"

y4 <- rbind(q0, q1, q2, q3, q4, q5)
y4v2 <- inner_join(y4, annotlcpm, by = c("gene.name" = "neargene.name"))

comp4 <- ggplot(data=y4v2, aes(x=genepeakpheno, y=rna_fold_fwt_dko)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(label=..count..), y=2, stat='count', colour="black", size=4) +
  theme_classic() +
  ylim(-2.5,2.5) +
  geom_hline(yintercept=0)

wrap_plots(comp1 | comp2 | comp3 | comp4)

#####
#####enhancer compensation, subset by high peak intensity #####
##naive
matrix_oi <- supermatrix.ii.k4m1k27a
matrix_oi$Distance.to.TSS <- abs(matrix_oi$Distance.to.TSS)
matrix_oi <- matrix_oi %>% filter(N_WT_k27a > 9)

x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a > -1))
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y2 <- y[order(-y$totalpeaks, -y$prop_lost),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

q0 <- y3 %>% filter(prop_lost == 1)
q0$genepeakpheno <- "q0"
q1 <- y3 %>% filter(prop_lost >= 0.75 & prop_lost < 1)
q1$genepeakpheno <- "q1"
q2 <- y3 %>% filter(prop_lost >= 0.50 & prop_lost < 0.75)
q2$genepeakpheno <- "q2"
q3 <- y3 %>% filter(prop_lost >= 0.25 & prop_lost < 0.50)
q3$genepeakpheno <- "q3"
q4 <- y3 %>% filter(prop_lost < 0.25 & prop_lost > 0)
q4$genepeakpheno <- "q4"
q5 <- y3 %>% filter(prop_lost == 0)
q5$genepeakpheno <- "q5"

y4 <- rbind(q0, q1, q2, q3, q4, q5)
y4v2 <- inner_join(y4, annotlcpm, by = c("gene.name" = "neargene.name"))

comp1 <- ggplot(data=y4v2, aes(x=genepeakpheno, y=rna_fold_nwt_dko)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(label=..count..), y=2, stat='count', colour="black", size=4) +
  theme_classic() +
  ylim(-2.5,2.5) +
  geom_hline(yintercept=0)


##formative
matrix_oi <- supermatrix.ii.k4m1k27a
x <- matrix_oi %>% filter(!is.na(NWT))
matrix_oi$Distance.to.TSS <- abs(matrix_oi$Distance.to.TSS)
matrix_oi <- matrix_oi %>% filter(F_WT_k27a > 9)

##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y2 <- y[order(-y$totalpeaks, -y$prop_lost),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

##create quintiles for enhancer proportion lost
q0 <- y3 %>% filter(prop_lost == 1)
q0$genepeakpheno <- "q0"
q1 <- y3 %>% filter(prop_lost >= 0.75 & prop_lost < 1)
q1$genepeakpheno <- "q1"
q2 <- y3 %>% filter(prop_lost >= 0.50 & prop_lost < 0.75)
q2$genepeakpheno <- "q2"
q3 <- y3 %>% filter(prop_lost >= 0.25 & prop_lost < 0.50)
q3$genepeakpheno <- "q3"
q4 <- y3 %>% filter(prop_lost < 0.25 & prop_lost > 0)
q4$genepeakpheno <- "q4"
q5 <- y3 %>% filter(prop_lost == 0)
q5$genepeakpheno <- "q5"

y4 <- rbind(q0, q1, q2, q3, q4, q5)
y4v2 <- inner_join(y4, annotlcpm, by = c("gene.name" = "neargene.name"))

comp2 <- ggplot(data=y4v2, aes(x=genepeakpheno, y=rna_fold_fwt_dko)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(aes(label=..count..), y=2, stat='count', colour="black", size=4) +
  theme_classic() +
  ylim(-2.5,2.5) +
  geom_hline(yintercept=0)

wrap_plots(comp1 + comp2)

#####
#####enhancer compensation weighted by number#####
##naive
matrix_oi <- supermatrix.ii.k4m1k27a
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a > -1))
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks
y$weightedloss <- log2(y$prop_lost * y$totalpeaks)
#y$weightedloss <- log2(y$deppeaks * y$totalpeaks)

y2 <- y[order(-y$weightedloss),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

y4v2 <- inner_join(y3, annotlcpm, by = c("gene.name" = "neargene.name"))

comp1 <- ggplot(data=y4v2, aes(x=weightedloss, y=rna_fold_nwt_dko)) +
  geom_point(size = 0.2) +
  #stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue") +
  ylim(-3,3) +
  xlim(0,5) +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4)) +
  theme_classic() +
  geom_hline(yintercept = 0)

##formative
matrix_oi <- supermatrix.ii.k4m1k27a
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y$weightedloss <- y$prop_lost * y$totalpeaks
#y$weightedloss <- log2(y$deppeaks * y$totalpeaks)
y$weightedloss <- log2(y$prop_lost * y$totalpeaks)

y2 <- y[order(-y$weightedloss),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

y4v2 <- inner_join(y3, annotlcpm, by = c("gene.name" = "neargene.name"))

comp2 <- ggplot(data=y4v2, aes(x=weightedloss, y=rna_fold_fwt_dko)) +
  geom_point(size = 0.2) +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4)) +
 # stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue") +
  ylim(-3,3) +
  xlim(0,5) +
  theme_classic() +
  geom_hline(yintercept = 0)


wrap_plots(comp1 + comp2)


#####
#####enhancer compensation weighted by number tss distance max 10kb#####
##naive
matrix_oi <- supermatrix.ii.k4m1k27a
matrix_oi$Distance.to.TSS <- abs(matrix_oi$Distance.to.TSS)
matrix_oi <- matrix_oi %>% filter(Distance.to.TSS < 10000)
x <- matrix_oi %>% filter(!is.na(NWT))

##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a > -1))
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y$weightedloss <- y$prop_lost * y$totalpeaks
y$weightedloss <- log2(y$prop_lost * y$totalpeaks)
#y$weightedloss <- log2(y$deppeaks * y$totalpeaks)

#y$weightedloss <- (y$deppeaks * y$totalpeaks)

y2 <- y[order(-y$weightedloss),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

y4v2 <- inner_join(y3, annotlcpm, by = c("gene.name" = "neargene.name"))

comp3 <- ggplot(data=y4v2, aes(x=weightedloss, y=rna_fold_nwt_dko)) +
  geom_point(size = 0.2) +
  #stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue") +
  ylim(-3,3) +
  xlim(0,2) +
  #geom_smooth(method="lm") +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4)) +
  theme_classic() +
  geom_hline(yintercept = 0)

##formative
matrix_oi <- supermatrix.ii.k4m1k27a
matrix_oi$Distance.to.TSS <- abs(matrix_oi$Distance.to.TSS)
matrix_oi <- matrix_oi %>% filter(Distance.to.TSS < 10000)
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y$weightedloss <- log2(y$prop_lost * y$totalpeaks)

y2 <- y[order(-y$weightedloss),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

y4v2 <- inner_join(y3, annotlcpm, by = c("gene.name" = "neargene.name"))

comp4 <- ggplot(data=y4v2, aes(x=weightedloss, y=rna_fold_fwt_dko)) +
  geom_point(size = 0.2) +
  #geom_smooth(method="lm") +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4)) +
  #stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue") +
  ylim(-3,3) +
  xlim(0,2) +
  theme_classic() +
  geom_hline(yintercept = 0)


wrap_plots(comp1 | comp2 | comp3 | comp4)


#####
#####enhancer compensation weighted by dep times total and new tmm#####
##naive
matrix_oi <- supermatrix.ii.k4m1k27a
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a > -1))
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

#y$weightedloss <- log2(y$prop_lost * y$totalpeaks)
y$weightedloss <- log2(y$deppeaks * y$totalpeaks)

y2 <- y[order(-y$weightedloss),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

y4v2 <- inner_join(y3, annotlcpm, by = c("gene.name" = "neargene.name"))

comp1 <- ggplot(data=y4v2, aes(x=weightedloss, y=rna_fold_nwt_dko2)) +
  geom_point(size = 0.2) +
  #stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue") +
  ylim(-3,3) +
  xlim(0,10) +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4)) +
  theme_classic() +
  geom_hline(yintercept = 0)

##formative
matrix_oi <- supermatrix.ii.k4m1k27a
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y$weightedloss <- y$prop_lost * y$totalpeaks
y$weightedloss <- log2(y$deppeaks * y$totalpeaks)
#y$weightedloss <- log2(y$prop_lost * y$totalpeaks)

y2 <- y[order(-y$weightedloss),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

y4v2 <- inner_join(y3, annotlcpm, by = c("gene.name" = "neargene.name"))

comp2 <- ggplot(data=y4v2, aes(x=weightedloss, y=rna_fold_fwt_dko2)) +
  geom_point(size = 0.2) +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4)) +
  # stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue") +
  ylim(-3,3) +
  xlim(0,10) +
  theme_classic() +
  geom_hline(yintercept = 0)


wrap_plots(comp1 + comp2)
#####
#####enhancer compensation weighted by dep times total and new tmm max 10kb#####
##naive
matrix_oi <- supermatrix.ii.k4m1k27a
matrix_oi$Distance.to.TSS <- abs(matrix_oi$Distance.to.TSS)
matrix_oi <- matrix_oi %>% filter(Distance.to.TSS < 10000)
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_n_dkowt_k27a > -1))
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

#y$weightedloss <- log2(y$prop_lost * y$totalpeaks)
y$weightedloss <- log2(y$deppeaks * y$totalpeaks)

y2 <- y[order(-y$weightedloss),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

y4v2 <- inner_join(y3, annotlcpm, by = c("gene.name" = "neargene.name"))

comp1 <- ggplot(data=y4v2, aes(x=weightedloss, y=rna_fold_nwt_dko2)) +
  geom_point(size = 0.2) +
  #stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue") +
  ylim(-3,3) +
  xlim(0,4) +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4)) +
  theme_classic() +
  geom_hline(yintercept = 0)

##formative
matrix_oi <- supermatrix.ii.k4m1k27a
matrix_oi$Distance.to.TSS <- abs(matrix_oi$Distance.to.TSS)
matrix_oi <- matrix_oi %>% filter(Distance.to.TSS < 10000)
x <- matrix_oi %>% filter(!is.na(NWT))
##change filter criteria to match peakset
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1))
genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a < -1))
#genes_deppeaks <- count(x %>% group_by(neargene.name) %>% filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000)))

colnames(genes_deppeaks) <- c("gene.name", "deppeaks")
#genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1))
genes_indpeaks <- count(x %>% group_by(neargene.name) %>% filter(pk_fold_f_dkowt_k27a > -1))

colnames(genes_indpeaks) <- c("gene.name", "indpeaks")

y <- full_join(genes_indpeaks, genes_deppeaks, by = "gene.name", .keep_all)
y[is.na(y)] <- 0
y$totalpeaks <- y$indpeaks + y$deppeaks
y$prop_lost <- y$deppeaks / y$totalpeaks

y$weightedloss <- y$prop_lost * y$totalpeaks
y$weightedloss <- log2(y$deppeaks * y$totalpeaks)
#y$weightedloss <- log2(y$prop_lost * y$totalpeaks)

y2 <- y[order(-y$weightedloss),] 
y2$rank <- seq(1, nrow(y2))

y3 <- left_join(y2, matrix_oi, by = c("gene.name" = "neargene.name"))
y3 <- y3 %>% distinct(gene.name, .keep_all = TRUE)

y4v2 <- inner_join(y3, annotlcpm, by = c("gene.name" = "neargene.name"))

comp2 <- ggplot(data=y4v2, aes(x=weightedloss, y=rna_fold_fwt_dko2)) +
  geom_point(size = 0.2) +
  geom_smooth(method="gam", formula = y ~ s(x, bs = "cs", k=4)) +
  # stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue") +
  ylim(-3,3) +
  xlim(0,4) +
  theme_classic() +
  geom_hline(yintercept = 0)


wrap_plots(comp1 + comp2)
#####
#####Distance cor plots?#####

matrix_oi2 <- supermatrix.down_k27a_dep
matrix_oi2$Distance.to.TSS <- abs(matrix_oi2$Distance.to.TSS)
matrix_oi2$log2Distance.to.TSS <- log2(matrix_oi2$Distance.to.TSS)

matrix_oi2[matrix_oi2['Distance.to.TSS'] > 250000, 'Distance.to.TSS'] = 250000
matrix_oi2[matrix_oi2['Distance.to.TSS'] > 100000, 'Distance.to.TSS'] = 100000


matrix_oi3 <- supermatrix.up_k27a_dep
matrix_oi3$Distance.to.TSS <- abs(matrix_oi3$Distance.to.TSS)
matrix_oi3$log2Distance.to.TSS <- log2(matrix_oi3$Distance.to.TSS)
matrix_oi3[matrix_oi3['Distance.to.TSS'] > 250000, 'Distance.to.TSS'] = 250000
matrix_oi3[matrix_oi3['Distance.to.TSS'] > 100000, 'Distance.to.TSS'] = 100000

matrix_oi4 <- supermatrix.k27a.non
matrix_oi4$Distance.to.TSS <- abs(matrix_oi4$Distance.to.TSS)
matrix_oi4$log2Distance.to.TSS <- log2(matrix_oi4$Distance.to.TSS)
matrix_oi4[matrix_oi4['Distance.to.TSS'] > 250000, 'Distance.to.TSS'] = 250000
matrix_oi4[matrix_oi4['Distance.to.TSS'] > 100000, 'Distance.to.TSS'] = 100000


ggplot() +
  geom_histogram(data=matrix_oi4, aes(x=Distance.to.TSS, y=stat(count / sum(count))), bins = 50, color="grey50", alpha=0.5) +
  geom_histogram(data=matrix_oi3, aes(x=Distance.to.TSS, y=stat(count / sum(count))), bins = 50, color="magenta3", alpha=0.2) +
  geom_histogram(data=matrix_oi2, aes(x=Distance.to.TSS, y=stat(count / sum(count))), bins = 50, color="deepskyblue3", alpha=0.2) +
  xlim(0,102000) +
  theme_classic()


matrix_oi2v2 <- matrix_oi2 %>% filter(Distance.to.TSS < 100000)
matrix_oi2v2 <-  matrix_oi2v2 %>% filter(!is.na(rna_fold_nwt_dko))
matrix_oi2v3 <- inner_join(matrix_oi2v2, annotlcpm, by = c("neargene.name" = "neargene.name"))



nwt_mean <- median(matrix_oi2v3$rna_fold_nwt_dko2)
nwt_distancemedian <- median(matrix_oi2v3$Distance.to.TSS)
nwt_log2distancemedian <- median(matrix_oi2v3$log2Distance.to.TSS)

ggplot() +
  geom_point(data=matrix_oi2v3, aes(x=rna_fold_nwt_dko2, y=Distance.to.TSS), size = 0.2, alpha = 0.5) +
  xlim(-5,5) +
  ylim(0, 50000) + 
  geom_vline(xintercept = nwt_mean) +
  geom_hline(yintercept = nwt_distancemedian) +
  theme_classic()

ggplot() +
  geom_point(data=matrix_oi2v3, aes(x=rna_fold_nwt_dko2, y=log2Distance.to.TSS), size = 0.2, alpha = 0.5) +
  xlim(-5,5) +
  ylim(11, 17) + 
  geom_vline(xintercept = nwt_mean) +
  geom_hline(yintercept = nwt_log2distancemedian) +
  theme_classic()


matrix_oi3v2 <- matrix_oi3 %>% filter(Distance.to.TSS < 100000)
matrix_oi3v2 <-  matrix_oi3v2 %>% filter(!is.na(rna_fold_fwt_dko))
matrix_oi3v3 <- inner_join(matrix_oi3v2, annotlcpm, by = c("neargene.name" = "neargene.name"))

fwt_mean <- median(matrix_oi3v3$rna_fold_fwt_dko2)
fwt_distancemedian <- median(matrix_oi3v3$Distance.to.TSS)
fwt_log2distancemedian <- mean(matrix_oi3v3$log2Distance.to.TSS)


ggplot() +
  geom_point(data=matrix_oi3v3, aes(x=rna_fold_fwt_dko2, y=Distance.to.TSS), size = 0.2, alpha = 0.5) +
  xlim(-5,5) +
  ylim(0, 50000) + 
  geom_vline(xintercept = fwt_mean) +
  geom_hline(yintercept = fwt_distancemedian) +
  theme_classic() 
  

ggplot() +
  geom_point(data=matrix_oi3v3, aes(x=rna_fold_fwt_dko2, y=log2Distance.to.TSS), size = 0.2, alpha = 0.5) +
  xlim(-5,5) +
  ylim(11,17) + 
  geom_vline(xintercept = fwt_mean) +
  geom_hline(yintercept = fwt_log2distancemedian) +
  theme_classic()


#####
#####Creation of gimme compatible matrix of y1y4 inddep#####
##form peaks
y_one <- read.csv("~/Desktop/data/y1ind_at.bed", sep='\t', header = FALSE)
colnames(y_one) <- c("chr", "start", "end")
y_one$peaktype <- "y1ind"
y_two <- read.csv("~/Desktop/data/y1dep_at.bed", sep='\t', header = FALSE)
colnames(y_two) <- c("chr", "start", "end")
y_two$peaktype <- "y1dep"
y_three <- read.csv("~/Desktop/data/y4ind_at.bed", sep='\t', header = FALSE)
colnames(y_three) <- c("chr", "start", "end")
y_three$peaktype <- "y4ind"
y_four <- read.csv("~/Desktop/data/y4dep_at.bed", sep='\t', header = FALSE)
colnames(y_four) <- c("chr", "start", "end")
y_four$peaktype <- "y4dep"

x <- rbind(y_one, y_two, y_three, y_four) 
gimme <- as.data.frame(paste0(x$chr, ":", x$start, "-", x$end))
gimme$cluster <- x$peaktype
colnames(gimme) <- c("loc", "cluster")
write.table(gimme, "~/Desktop/gimme2022ypeaks.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep='\t')


##naive peaks
y_one <- read.csv("~/Desktop/data/z1ind_at.bed", sep='\t', header = FALSE)
colnames(y_one) <- c("chr", "start", "end")
y_one$peaktype <- "z1ind"
y_two <- read.csv("~/Desktop/data/z1dep_at.bed", sep='\t', header = FALSE)
colnames(y_two) <- c("chr", "start", "end")
y_two$peaktype <- "z1dep"
y_three <- read.csv("~/Desktop/data/z3ind_at.bed", sep='\t', header = FALSE)
colnames(y_three) <- c("chr", "start", "end")
y_three$peaktype <- "z4ind"
y_four <- read.csv("~/Desktop/data/z3dep_at.bed", sep='\t', header = FALSE)
colnames(y_four) <- c("chr", "start", "end")
y_four$peaktype <- "z4dep"

x <- rbind(y_one, y_two, y_three, y_four) 
gimme <- as.data.frame(paste0(x$chr, ":", x$start, "-", x$end))
gimme$cluster <- x$peaktype
colnames(gimme) <- c("loc", "cluster")
write.table(gimme, "~/Desktop/gimme2022zpeaks.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep='\t')




#####
#####
