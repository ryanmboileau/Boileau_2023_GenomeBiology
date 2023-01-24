#Similar code to epigenome_master_k4m1 but focused on H3K27ac peaks. Underlies primarily Fig.3

library(ggplot2)
library(tidyverse)
library(patchwork)
library(gridExtra)
library(grid)
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
#   [1] ggpubr_0.4.0    gridExtra_2.3   patchwork_1.1.1 forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4     readr_2.1.2    
# [9] tidyr_1.2.0     tibble_3.1.7    tidyverse_1.3.1 ggplot2_3.3.6  
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.1.2 haven_2.5.0      carData_3.0-5    colorspace_2.0-3 vctrs_0.4.1      generics_0.1.3   utf8_1.2.2       rlang_1.0.3     
# [9] pillar_1.7.0     glue_1.6.2       withr_2.5.0      DBI_1.1.3        dbplyr_2.2.1     modelr_0.1.8     readxl_1.4.0     lifecycle_1.0.1 
# [17] plyr_1.8.7       munsell_0.5.0    ggsignif_0.6.3   gtable_0.3.0     cellranger_1.1.0 rvest_1.0.2      tzdb_0.3.0       fansi_1.0.3     
# [25] broom_1.0.0      Rcpp_1.0.8.3     backports_1.4.1  scales_1.2.0     jsonlite_1.8.0   abind_1.4-5      fs_1.5.2         hms_1.1.1       
# [33] stringi_1.7.6    rstatix_0.7.0    cli_3.3.0        tools_4.1.2      magrittr_2.0.3   car_3.1-0        crayon_1.5.1     pkgconfig_2.0.3 
# [41] ellipsis_0.3.2   xml2_1.3.3       reprex_2.0.1     lubridate_1.8.0  assertthat_0.2.1 reshape_0.8.9    httr_1.4.3       rstudioapi_0.13 
# [49] R6_2.5.1         compiler_4.1.2 

#####
#####load required functions#####

generate_deepheatmapgrouped <- function(listofinputs, outputfilename){
  "input a list of bed file coordinates chr,start,end and output a single txt file with # separating groups for use with deeptools plot heatmap"
  stuffer <- "#"
  write.table(listofinputs[[1]][1:3], outputfilename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  for(i in 2:length(listofinputs)){
    nextgroup <- listofinputs[[i]][1:3]
    write.table(stuffer, outputfilename, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
    write.table(nextgroup, outputfilename, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t') 
    
  }
}
homerpeakannotationcleanup <- function(mergedannottable) {
  
  # df_annot <- as.data.frame(lapply(mergedannottable, function(x) gsub("intron.*", "Intron", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("exon.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Exon.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Intron.*", "Intron", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("promoter.*", "Promoter-TSS", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("TTS.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3' UTR.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5' UTR.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3 UTR.*", "Other", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5 UTR.*", "Other", x)))
  
  df_annot <- as.data.frame(lapply(mergedannottable, function(x) gsub("intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("exon.*", "Exon", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Exon.*", "Exon", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("promoter.*", "Promoter-TSS", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("TTS.*", "TTS", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3' UTR.*", "3'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5' UTR.*", "5'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3 UTR.*", "3'UTR", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5 UTR.*", "5'UTR", x)))
  
  #df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("p.*", "Promoter-TSS", x)))
  # 
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_1$", "cluster_01", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_2", "cluster_02", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_3", "cluster_03", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_4", "cluster_04", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_5", "cluster_05", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_6", "cluster_06", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_7", "cluster_07", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_8", "cluster_08", x)))
  # df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("cluster_9", "cluster_09", x)))
}

#####
#####Import data sets ######
# 
# deseq_wt_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_WT.txt", sep = '\t', header=TRUE)
# deseq_cko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_CKO.txt", sep = '\t', header=TRUE)
# deseq_dko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_DKO.txt", sep = '\t', header=TRUE)
# deseq_dcd_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_dCD.txt", sep = '\t', header=TRUE)
# 
# anno_wt_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_WT.anno.txt", sep = '\t', header=TRUE)
# anno_cko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_CKO.anno.txt", sep = '\t', header=TRUE)
# anno_dko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_DKO.anno.txt", sep = '\t', header=TRUE)
# anno_dcd_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_dCD.anno.txt", sep = '\t', header=TRUE)
# 
# readcounts_k4m1_peaks_wt <- read.table('~/Desktop/readcounts_deseq_RBRB11_k4m1.tab', sep = '\t', header=TRUE)

#generate homer compatible bed files for annotation + RNAseq data appending
# fileoi <- cbind(x[1:3], seq(1,nrow(x)))
# filename <- "~/Desktop/data/k27a_diffbind_results/cell_state_internal/diffbindresults_k27a_fDKO_WT.homer.txt"
# write.table(fileoi, filename, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')

anno_wt_k27a <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k27a_nf_WT.anno.txt", sep = '\t', header=TRUE)
anno_dko_k27a <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k27a_nf_DKO.anno.txt", sep = '\t', header=TRUE)

# anno_wt_k27m3 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k27m3_nf_WT.anno.txt", sep = '\t', header=TRUE)
# anno_wt_k4m3 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k4m3_nf_WT.anno.txt", sep = '\t', header=TRUE)
# anno_wt_rad21 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_Rad21_nf_WT.anno.txt", sep = '\t', header=TRUE)

# anno_dko_k27m3 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k27m3_nf_DKO.anno.txt", sep = '\t', header=TRUE)
# anno_dko_k4m3 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k4m3_nf_DKO.anno.txt", sep = '\t', header=TRUE)
# anno_dko_rad21 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_Rad21_nf_DKO.anno.txt", sep = '\t', header=TRUE)

# deseq_wt_k27m3 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k27m3_nf_WT.txt", sep = '\t', header=TRUE)
# deseq_wt_k4m3 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k4m3_nf_WT.txt", sep = '\t', header=TRUE)
# deseq_wt_rad21 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_Rad21_nf_WT.txt", sep = '\t', header=TRUE)


deseq_wt_k27a <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k27a_nf_WT.txt", sep = '\t', header=TRUE)
deseq_dko_k27a <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k27a_nf_DKO.txt", sep = '\t', header=TRUE)
deseq_wt_k27a <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k27a_nfDKOatWT.txt_full.txt", sep = '\t', header=TRUE)

# deseq_dko_k27m3 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k27m3_nf_DKO.txt", sep = '\t', header=TRUE)
# deseq_dko_k4m3 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_k4m3_nf_DKO.txt", sep = '\t', header=TRUE)
# deseq_dko_rad21 <- read.csv("~/Desktop/data/k27a_diffbind_results/diffbindresults_Rad21_nf_DKO.txt", sep = '\t', header=TRUE)

# readcounts_k27a_peaks_wt <- read.csv("~/Desktop/data/k27a_diffbind_results/readcounts_k27apeaks_allmarks_WT.tab", sep = '\t', header=TRUE, quote = "")
# readcounts_k27m3_peaks_wt <- read.csv("~/Desktop/data/k27a_diffbind_results/readcounts_k27m3peaks_allmarks_WT.tab", sep = '\t', header=TRUE, quote = "")
# readcounts_k4m3_peaks_wt <- read.csv("~/Desktop/data/k27a_diffbind_results/readcounts_k4m3peaks_allmarks_WT.tab", sep = '\t', header=TRUE, quote = "")
# readcounts_rad21_peaks_wt <- read.csv("~/Desktop/data/k27a_diffbind_results/readcounts_Rad21peaks_allmarks_WT.tab", sep = '\t', header=TRUE, quote = "")

readcounts_k27a_peaks_wt_at <- read.csv("~/Desktop/data/readcounts_k27apeaks_allmarks_WT_v2.tab", sep = '\t', header=TRUE, quote = "")
readcounts_k27a_peaks_dko_at <- read.csv("~/Desktop/data/readcounts_k27a_nf_dko_allmarks.tab", sep = '\t', header=TRUE, quote = "")

# 
# readcounts_k27a_peaks_dko <- read.csv("~/Desktop/data/k27a_diffbind_results/readcounts_k27a_nf_DKO.tab", sep = '\t', header=TRUE, quote = "")
# readcounts_k27m3_peaks_dko <- read.csv("~/Desktop/data/k27a_diffbind_results/readcounts_k27m3_nf_DKO.tab", sep = '\t', header=TRUE, quote = "")
# readcounts_k4m3_peaks_dko <- read.csv("~/Desktop/data/k27a_diffbind_results/readcounts_k4m3_nf_DKO.tab", sep = '\t', header=TRUE, quote = "")
# readcounts_rad21_peaks_dko <- read.csv("~/Desktop/data/k27a_diffbind_results/readcounts_rad21_nf_DKO.tab", sep = '\t', header=TRUE, quote = "")


#####
######Integrate data frames######

#choose peak set and corresponding data sets on which to run script
anno_oi <- anno_wt_k27a
deseq_oi <- deseq_wt_k27a
readcounts_oi <- readcounts_k27a_peaks_wt_at

anno_oi <- anno_dko_k27a
deseq_oi <- deseq_dko_k27a
readcounts_oi <- readcounts_k27a_peaks_dko_at



# filename.up.df <- "~/Desktop/data/sigpeaks_k27a_wt_up.bed"
# filename.down.df <- "~/Desktop/data/sigpeaks_k27m3_wt_down.bed"
# filename.non.df <- "~/Desktop/data/sigpeaks_k27m3_wt_non.bed"

##massage data into mergeable/joinable dfs, start by making uniform labels and first layer processing of primary data
anno_oi <- anno_oi[2:ncol(anno_oi)]
anno_oi2 <- cbind(anno_oi[1:3], anno_oi[9], anno_oi[20:27])
anno_oi2$simple_annot <- t(homerpeakannotationcleanup(anno_oi$Annotation))
anno_oi2$neargene.name <- anno_oi$Gene.Name

deseq_oi$chr <- deseq_oi$seqnames
deseq_oi <- cbind(deseq_oi[12], deseq_oi[2:11])
deseq_oi$start <- deseq_oi$start + 1

colnames(readcounts_oi) <- c("chr", "start", "end", "N_WT_k4m1", "F_WT_k4m1","N_DKO_k4m1", "F_DKO_k4m1", "N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a", "N_WT_k4m3", "F_WT_k4m3", "N_DKO_k4m3", "F_DKO_k4m3",
                             "N_WT_rad21", "F_WT_rad21", "N_DKO_rad21", "F_DKO_rad21", "N_WT_k27m3", "F_WT_k27m3", "N_DKO_k27m3", "F_DKO_k27m3", "N_WT_AT", "F_WT_AT", "N_DKO_AT", "F_DKO_AT")


# readcounts_oi$width <- readcounts_oi$end - readcounts_oi$start
# readcounts_oi_density <- readcounts_oi[4:27] / readcounts_oi$width
# readcounts_oi_ldensity <- log2(readcounts_oi_density * 1000000 + 1) #multiply by constant, add 1 to avoid divide by zero and log2 scale
# readcounts_oi2 <- cbind(readcounts_oi[1:3], readcounts_oi_ldensity)
# readcounts_oi2$start <- readcounts_oi2$start + 1

readcounts_oi$width <- readcounts_oi$end - readcounts_oi$start
readcounts_oi_density <- readcounts_oi[4:31] / readcounts_oi$width
readcounts_oi_ldensity <- log2(readcounts_oi_density * 1000000 + 1) #multiply by constant, add 1 to avoid divide by zero and log2 scale
readcounts_oi2 <- cbind(readcounts_oi[1:3], readcounts_oi_ldensity)
readcounts_oi2$start <- readcounts_oi2$start + 1
######
#####second layer processing for folds#####
##perform additional processing to achieve second layer data by conducting log fold changes etc. 
#RNAseq fold changes in anno sheet
rna_fold_nf_wt <- anno_oi2$FWT - anno_oi2$NWT
rna_fold_nf_cko <- anno_oi2$FCKO - anno_oi2$NCKO
rna_fold_nf_dko <- anno_oi2$FDKO - anno_oi2$NDKO
rna_fold_nf_dcd <- anno_oi2$FdCD - anno_oi2$NdCD

rna_fold_nwt_cko <- anno_oi2$NCKO - anno_oi2$NWT
rna_fold_nwt_dko <- anno_oi2$NDKO - anno_oi2$NWT
rna_fold_nwt_dcd <- anno_oi2$NdCD - anno_oi2$NWT

rna_fold_fwt_cko <- anno_oi2$FCKO - anno_oi2$FWT
rna_fold_fwt_dko <- anno_oi2$FDKO - anno_oi2$FWT
rna_fold_fwt_dcd <- anno_oi2$FdCD - anno_oi2$FWT

anno_oi3 <- cbind(anno_oi2, rna_fold_nf_wt, rna_fold_nf_cko,rna_fold_nf_dko,rna_fold_nf_dcd,rna_fold_nwt_cko,rna_fold_nwt_dko,
                  rna_fold_nwt_dcd,rna_fold_fwt_cko,rna_fold_fwt_dko, rna_fold_fwt_dcd)
#anno_oi3 <- cbind(anno_oi2, rna_fold_nf_wt, rna_fold_nf_dko,rna_fold_nwt_dko, rna_fold_fwt_dko)

#readcount fold changes in multibigwigsummary
pk_fold_nf_wt_k4m1 <- readcounts_oi2$F_WT_k4m1 - readcounts_oi2$N_WT_k4m1
pk_fold_nf_wt_k27a <- readcounts_oi2$F_WT_k27a - readcounts_oi2$N_WT_k27a
pk_fold_nf_wt_k4m3 <- readcounts_oi2$F_WT_k4m3 - readcounts_oi2$N_WT_k4m3
pk_fold_nf_wt_k27m3 <- readcounts_oi2$F_WT_k27m3 - readcounts_oi2$N_WT_k27m3
pk_fold_nf_wt_rad21 <- readcounts_oi2$F_WT_rad21 - readcounts_oi2$N_WT_rad21

pk_fold_nf_dko_k4m1 <- readcounts_oi2$F_DKO_k4m1 - readcounts_oi2$N_DKO_k4m1
pk_fold_nf_dko_k27a <- readcounts_oi2$F_DKO_k27a - readcounts_oi2$N_DKO_k27a
pk_fold_nf_dko_k4m3 <- readcounts_oi2$F_DKO_k4m3 - readcounts_oi2$N_DKO_k4m3
pk_fold_nf_dko_k27m3 <- readcounts_oi2$F_DKO_k27m3 - readcounts_oi2$N_DKO_k27m3
pk_fold_nf_dko_rad21 <- readcounts_oi2$F_DKO_rad21 - readcounts_oi2$N_DKO_rad21

pk_fold_f_dkowt_k4m1 <- readcounts_oi2$F_DKO_k4m1 - readcounts_oi2$F_WT_k4m1
pk_fold_f_dkowt_k27a <- readcounts_oi2$F_DKO_k27a - readcounts_oi2$F_WT_k27a
pk_fold_f_dkowt_k4m3 <- readcounts_oi2$F_DKO_k4m3 - readcounts_oi2$F_WT_k4m3
pk_fold_f_dkowt_k27m3 <- readcounts_oi2$F_DKO_k27m3 - readcounts_oi2$F_WT_k27m3
pk_fold_f_dkowt_rad21 <- readcounts_oi2$F_DKO_rad21 - readcounts_oi2$F_WT_rad21

pk_fold_n_dkowt_k4m1 <- readcounts_oi2$N_DKO_k4m1 - readcounts_oi2$N_WT_k4m1
pk_fold_n_dkowt_k27a <- readcounts_oi2$N_DKO_k27a - readcounts_oi2$N_WT_k27a
pk_fold_n_dkowt_k4m3 <- readcounts_oi2$N_DKO_k4m3 - readcounts_oi2$N_WT_k4m3 
pk_fold_n_dkowt_k27m3 <-readcounts_oi2$N_DKO_k27m3 - readcounts_oi2$N_WT_k27m3
pk_fold_n_dkowt_rad21 <- readcounts_oi2$N_DKO_rad21 - readcounts_oi2$N_WT_rad21

pk_fold_nf_wt_at <- readcounts_oi2$F_WT_AT - readcounts_oi2$N_WT_AT 
pk_fold_nf_dko_at <- readcounts_oi2$F_DKO_AT - readcounts_oi2$N_DKO_AT

pk_fold_f_dkowt_at <- readcounts_oi2$F_DKO_AT - readcounts_oi2$F_WT_AT 
pk_fold_n_dkowt_at <- readcounts_oi2$N_DKO_AT - readcounts_oi2$N_WT_AT

##
# readcounts_oi3 <- cbind(readcounts_oi2, pk_fold_nf_wt_k4m1, pk_fold_nf_wt_k27a,pk_fold_nf_wt_k4m3,pk_fold_nf_wt_k27m3,pk_fold_nf_wt_rad21,
#                         pk_fold_nf_dko_k4m1,pk_fold_nf_dko_k27a,pk_fold_nf_dko_k4m3,pk_fold_nf_dko_k27m3,pk_fold_nf_dko_rad21,
#                         pk_fold_n_dkowt_k4m1,pk_fold_n_dkowt_k27a,pk_fold_n_dkowt_k4m3,pk_fold_n_dkowt_k27m3,pk_fold_n_dkowt_rad21,
#                         pk_fold_f_dkowt_k4m1,pk_fold_f_dkowt_k27a,pk_fold_f_dkowt_k4m3,pk_fold_f_dkowt_k27m3, pk_fold_f_dkowt_rad21,
#                         pk_fold_nf_cko_k4m1,pk_fold_nf_dcd_k4m1)

readcounts_oi3 <- cbind(readcounts_oi2, pk_fold_nf_wt_k4m1, pk_fold_nf_wt_k27a,pk_fold_nf_wt_k4m3,pk_fold_nf_wt_k27m3,pk_fold_nf_wt_rad21,
                        pk_fold_nf_dko_k4m1,pk_fold_nf_dko_k27a,pk_fold_nf_dko_k4m3,pk_fold_nf_dko_k27m3,pk_fold_nf_dko_rad21,
                        pk_fold_n_dkowt_k4m1,pk_fold_n_dkowt_k27a,pk_fold_n_dkowt_k4m3,pk_fold_n_dkowt_k27m3,pk_fold_n_dkowt_rad21,
                        pk_fold_f_dkowt_k4m1,pk_fold_f_dkowt_k27a,pk_fold_f_dkowt_k4m3,pk_fold_f_dkowt_k27m3, pk_fold_f_dkowt_rad21,
                        pk_fold_nf_wt_at, pk_fold_nf_dko_at, pk_fold_f_dkowt_at, pk_fold_n_dkowt_at)


##combine into SUPER matrix SUPERRRR
supermatrix1 <- full_join(anno_oi3, deseq_oi, by = c("Chr" = "chr", "Start" = "start", "End" = "end"), keep=TRUE)
supermatrix <- full_join(supermatrix1, readcounts_oi3, by = c("Chr" = "chr", "Start" = "start", "End" = "end"), keep=TRUE)

#write.table(supermatrix, "~/Desktop/data/processed_data/supermatrix_k27a.20220427.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

#####
#####skip processing steps#####
supermatrix <- read.delim("~/Desktop/data/processed_data/supermatrix_k27a.20220427.txt", sep = '\t', header = TRUE)

## use with wt data
k27a_atac_filter <- read.table("~/Desktop/k27a_at_peaksfiltered.bed", sep = '\t')
colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
supermatrix <- inner_join(supermatrix, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))



#####

#####filter all peaks by atac data#####

k27a_atac_filter <- read.table("~/Desktop/data/k4m1k27a_k27apeaksbyk4m1at.bed", sep = '\t')
colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
supermatrix <- inner_join(supermatrix, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

k27a_atac_filter <- read.table("~/Desktop/data/k4m1k27a_k27apeaksnotk4m1_at.bed", sep = '\t')
colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
supermatrix <- inner_join(supermatrix, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))




## use with dko data
# k27a_atac_filter <- read.table("~/Desktop/k27apeaks_at_peaksfiltered_DKO.bed", sep = '\t')
# colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
# k27a_atac_filter$start <- k27a_atac_filter$start + 1
# supermatrix <- inner_join(supermatrix, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

# k27a_k4m1sub_filter <- read.table("~/Desktop/data/diffbindresults_k27awt_subbyk4m1wt.txt", sep = '\t')
# colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
# k27a_atac_filter$start <- k27a_atac_filter$start + 1
# supermatrix <- inner_join(supermatrix, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))



#####


#####fig3a foldchange vs FDR plots#####

foldchangedf <- select(supermatrix, pk_fold_nf_wt_k27a)

foldchangedf$FDR <- supermatrix$FDR
foldchangedf$log10fdr <- -1 * log10(supermatrix$FDR) 
foldchangedf2 <- foldchangedf %>% drop_na()

foldchangedf2sig <- foldchangedf2 %>% filter(log10fdr > 1.3 & (pk_fold_nf_wt_k27a > 1 | pk_fold_nf_wt_k27a < -1))

count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_wt_k27a > 1 | pk_fold_nf_wt_k27a < -1))) #21726 sites
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_wt_k27a > 1))) #8924 sites
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_wt_k27a < -1))) #12802 sites


ggplot() +
  geom_point(data=foldchangedf2, aes(x=pk_fold_nf_wt_k27a, y=log10fdr), color = "gray", size = 0.3) +
  geom_point(data=foldchangedf2sig, aes(x=pk_fold_nf_wt_k27a, y=log10fdr), color = "black", size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=1.30103, linetype = "dashed") +
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0,80)

# ggplot() +
#   geom_point(data=foldchangedf2, aes(x=pk_fold_nf_cko_k4m1, y=log10fdr), color = "gray", size = 0.3) +
#   geom_point(data=foldchangedf2sig, aes(x=pk_fold_nf_cko_k4m1, y=log10fdr), color = "black", size = 0.3) +
#   theme_bw() +
#   geom_hline(yintercept=1.30103, linetype = "dashed") +
#   theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   ylim(0,15)
# 
# ggplot() +
#   geom_point(data=foldchangedf2, aes(x=pk_fold_nf_dko_k4m1, y=log10fdr), color = "gray", size = 0.3) +
#   #geom_point(data=foldchangedf2sig, aes(x=pk_fold_nf_wt_k4m1, y=log10fdr), color = "black", size = 0.1, alpha = 0.2) +
#   theme_bw() +
#   geom_hline(yintercept=1.30103, linetype = "dashed") +
#   theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   ylim(0,4)
# 
# ggplot() +
#   geom_point(data=foldchangedf2, aes(x=pk_fold_nf_dcd_k4m1, y=log10fdr), color = "gray", size = 0.3) +
#   geom_point(data=foldchangedf2sig, aes(x=pk_fold_nf_dcd_k4m1, y=log10fdr), color = "black", size = 0.3) +
#   theme_bw() +
#   geom_hline(yintercept=1.30103, linetype = "dashed") +
#   theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   ylim(0,6)




###DKO sig results
foldchangedf <- select(supermatrix, pk_fold_nf_dko_k27a)

foldchangedf$FDR <- supermatrix$FDR
foldchangedf$log10fdr <- -1 * log10(supermatrix$FDR) 
foldchangedf2 <- foldchangedf %>% drop_na()

foldchangedf2sig <- foldchangedf2 %>% filter(log10fdr > 1.3 & (pk_fold_nf_dko_k27a > 1 | pk_fold_nf_dko_k27a < -1))

count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_dko_k27a > 1 | pk_fold_nf_dko_k27a < -1))) #6230sites
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_dko_k27a > 1))) #1935sites
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_dko_k27a < -1))) #4295 sites


ggplot() +
  geom_point(data=foldchangedf2, aes(x=pk_fold_nf_dko_k27a, y=log10fdr), color = "gray", size = 0.3) +
  geom_point(data=foldchangedf2sig, aes(x=pk_fold_nf_dko_k27a, y=log10fdr), color = "black", size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=1.30103, linetype = "dashed") +
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0,50)

#####

######peak filter parameters: up down nonsig, standard <0.5 1FC, >0.1 FC <0.7#####
# all wt
# supermatrix.up <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a > 1)
# supermatrix.down <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a < -1)
# supermatrix.non <- supermatrix %>% filter(FDR > 0.1 & pk_fold_nf_wt_k27a > -0.7 & pk_fold_nf_wt_k27a < 0.7)

# all dko
# supermatrix.up <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_dko_k27a > 1)
# supermatrix.down <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_dko_k27a < -1)
# supermatrix.non <- supermatrix %>% filter(FDR > 0.1 & pk_fold_nf_dko_k27a > -0.7 & pk_fold_nf_dko_k27a < 0.7)


#distal wt only
supermatrix.up <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.down <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.non <- supermatrix %>% filter(FDR > 0.1 & pk_fold_nf_wt_k27a > -0.7 & pk_fold_nf_wt_k27a < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

#sharedsites 12205
#downsites 9963
#upsites 7205

# write.table(supermatrix.up[1:3], "~/Desktop/supermatrix.up.k27a.ii.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE )
# write.table(supermatrix.non[1:3], "~/Desktop/supermatrix.non.k27a.ii.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE )
# write.table(supermatrix.down[1:3], "~/Desktop/supermatrix.down .k27a.ii.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE )



supermatrix.down.k27adown <- supermatrix.down %>% filter(pk_fold_n_dkowt_k27a < -1 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.down.k27anon <- supermatrix.down %>% filter(pk_fold_n_dkowt_k27a > -0.7 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

supermatrix.up.k27adown <- supermatrix.up %>% filter(pk_fold_f_dkowt_k27a < -1 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.up.k27anon <- supermatrix.up %>% filter(pk_fold_f_dkowt_k27a > -0.7 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

supermatrix.non.k27adown <- supermatrix.non %>% filter(pk_fold_f_dkowt_k27a < -1 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.non.k27anon <- supermatrix.non %>% filter(pk_fold_f_dkowt_k27a > -0.7 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

supermatrix.ii <- supermatrix %>%  filter((simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

# generate_deepheatmapgrouped(list(supermatrix.down.k27anon, supermatrix.down.k27adown), "~/Desktop/singlek27a_wtdown.bed")
# generate_deepheatmapgrouped(list(supermatrix.up.k27anon, supermatrix.up.k27adown), "~/Desktop/singlek27a_wtup.bed")
# generate_deepheatmapgrouped(list(supermatrix.non.k27anon, supermatrix.non.k27adown), "~/Desktop/singlek27a_wtnon.bed")
# 
# 
# supermatrix.fdkodf.non <-supermatrix.fdkodf %>% filter(FDR > 0.1 & pk_fold_f_dkowt_k27a > -0.7 & pk_fold_nf_wt_k4m1 < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))
# supermatrix.fdkodf.down <- supermatrix.fdkodf %>% filter(FDR < 0.05 & pk_fold_f_dkowt_k27a < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))
# 
# x <- semi_join(supermatrix.up, supermatrix.fdkodf.down, by=c("Chr", "Start", "End"))
# xx <- semi_join(supermatrix.non, supermatrix.fdkodf.down, by=c("Chr", "Start", "End"))
# 
# y <- semi_join(supermatrix.up, supermatrix.fdkodf.non, by=c("Chr", "Start", "End"))
# yy <- semi_join(supermatrix.non, supermatrix.fdkodf.non, by=c("Chr", "Start", "End"))



# supermatrix.up.k27adown.k4m1down <- supermatrix.up.k27adown %>% filter(pk_fold_nf_wt_k4m1 > 1)
# supermatrix.up.k27adown.k4m1non <- supermatrix.up.k27adown %>% filter(pk_fold_nf_wt_k4m1 < 0.7 & pk_fold_nf_wt_k4m1 > -0.7)
# 
# supermatrix.up.k27anon.k4m1down <- supermatrix.up.k27anon %>% filter(pk_fold_nf_wt_k4m1 > 1)
# supermatrix.up.k27anon.k4m1non <- supermatrix.up.k27anon %>% filter(pk_fold_nf_wt_k4m1 < 0.7 & pk_fold_nf_wt_k4m1 > -0.7)
# 
# listx <- list(supermatrix.up.k27anon.k4m1non, supermatrix.up.k27anon.k4m1down, supermatrix.up.k27adown.k4m1non, supermatrix.up.k27adown.k4m1down )
# listxx <- list(supermatrix.up.k27anon.k4m1non, supermatrix.up.k27anon.k4m1down, supermatrix.up.k27adown.k4m1non, supermatrix.up.k27adown.k4m1down )
# 
# 
# generate_deepheatmapgrouped(list(supermatrix.up, supermatrix.down), "~/Desktop/xx.bed")
# 
# supermatrix.tss <- supermatrix %>% filter(simple_annot == "Promoter-TSS")
# 
# supermatrix.tss.up <- supermatrix.tss %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a > 1 & simple_annot == "Promoter-TSS")
# supermatrix.tss.down <- supermatrix.tss %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a < -1 & simple_annot == "Promoter-TSS")
# supermatrix.tss.non <- supermatrix.tss %>% filter(FDR > 0.1 & pk_fold_nf_wt_k27a > -0.7 & pk_fold_nf_wt_k27a < 0.7 & simple_annot == "Promoter-TSS" )
# 
# 
# 
# 
# 
# supermatrix.non.k27adown <- supermatrix.non %>% filter(pk_fold_f_dkowt_k27a < -1 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
# supermatrix.non.k27anon <- supermatrix.non %>% filter(pk_fold_f_dkowt_k27a > -0.7 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
# 
# xx <- supermatrix.up.k27anon %>% filter(pk_fold_f_dkowt_k4m1 < -1 )
# xx <- supermatrix.up.k27anon %>% filter(pk_fold_f_dkowt_k4m1 < 0.7 & pk_fold_f_dkowt_k4m1 > -0.7 )
# 
# xx <- supermatrix.up.k27adown %>% filter(pk_fold_f_dkowt_k4m1 < -1 )
# 
# xx <- supermatrix.up.k27adown %>% filter(pk_fold_f_dkowt_k4m1 < 0.7 & pk_fold_f_dkowt_k4m1 > -0.7 )



##Premature K27ac
# supermatrix.up <- supermatrix %>% filter(FDR < 0.05 & Fold > 1 & pk_fold_n_dkowt_k27a > 5 & pk_fold_nf_wt_k27a > 5 )# & simple_annot == "Promoter-TSS")
# supermatrix.down <- supermatrix %>% filter(FDR < 0.05 & Fold < -1 & pk_fold_n_dkowt_k27a > 5 & pk_fold_nf_wt_k27a > 5 )# & simple_annot == "Promoter-TSS")
# supermatrix.non <- supermatrix %>% filter(FDR > 0.1 & Fold > -0.7 & Fold < 0.7 & pk_fold_n_dkowt_k27a > 5 & pk_fold_nf_wt_k27a > 5 )# & simple_annot == "Promoter-TSS" )
# filename.up.df <- "~/Desktop/premature_k27a.bed"
# write.table(supermatrix.up[1:3], filename.up.df, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
# 
# 

#subset sig peaks by sig RNAseq
wtsigrna <- read.csv("~/Desktop/data/Diffgenes/DiffgenesNWTvsFWTv3-07.07.20.anno.csv")
wtsigrnagenes <- wtsigrna$ensembllistanno



# supermatrix.wtsigrna <- subset(supermatrix, neargene.name %in% wtsigrnagenes)
# 
# supermatrix.up <- supermatrix.wtsigrna %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m3 > 1 & simple_annot == "Promoter-TSS")
# supermatrix.down <- supermatrix.wtsigrna %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m3 < -1 & simple_annot == "Promoter-TSS")
# supermatrix.non <- supermatrix.wtsigrna %>% filter(FDR > 0.1 & pk_fold_nf_wt_k4m3 > -0.7 & pk_fold_nf_wt_k4m3 < 0.7 & simple_annot == "Promoter-TSS" )


#####
#####peak annotation post filter #####
supermatrix.up$peaktype <- c("Sig.up")
supermatrix.down$peaktype <- c("Sig.down")
supermatrix.non$peaktype <- c("Nonsig")

peaktype_df <- rbind(supermatrix.up, supermatrix.down, supermatrix.non)

####proportion of annotations underlying peaks
peaktype_df_annot <- as.data.frame(cbind(peaktype_df$peaktype, peaktype_df$simple_annot))
colnames(peaktype_df_annot) <- c("peaktype", "simpleannot")


#literally cannot figure out how to sub this non-coding with other. 
#peaktype_df_annot2 <- t(as.data.frame(lapply(peaktype_df_annot$simpleannot, function(x) gsub("non-coding*", "Other2", x))))
peaktype_df_annot$simpleannot <- factor(peaktype_df_annot$simpleannot, levels = c("TTS", "3'UTR", "5'UTR", "Exon", "Promoter-TSS", "Intergenic", "Intron"))

ggplot(peaktype_df_annot, aes(x=peaktype, fill = simpleannot)) +
  geom_bar(width = 0.9) +
  theme_classic() +
  ylim(c(0, 25000)) +
  theme(aspect.ratio = 2, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1))#, legend.position = "none")


# ggplot(supermatrix.up_k27a_dep, aes(x=simple_annot, fill = simple_annot)) +
#   geom_bar(width = 0.9) +
#   theme_bw() +
#   theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1))



#####
######Nearest TSS on up, down, non primary filter#####
#supermatrix.up.rna <- gather(supermatrix.up[15:24])
supermatrix.up.rna <- gather(supermatrix.up[15:18])
supermatrix.up.rna <- gather(supermatrix.up_k27a[15:18])
supermatrix.up.rna$peaktype <- c("Sig.Up")

#supermatrix.down.rna <- gather(supermatrix.down[15:24])
supermatrix.down.rna <- gather(supermatrix.down[15:18])
supermatrix.down.rna$peaktype <- c("Sig.Down")
supermatrix.down.rna <- gather(supermatrix.down_k27a[15:18])

#supermatrix.non.rna <- gather(supermatrix.non[15:24])
supermatrix.non.rna <- gather(supermatrix.non[15:18])
supermatrix.non.rna$peaktype <- c("Non-Sig.")
supermatrix.non.rna <- gather(supermatrix.non_k27a[15:18])

supermatrix.rna <- rbind(supermatrix.up.rna,supermatrix.down.rna,supermatrix.non.rna)
supermatrix.rna$key <- factor(supermatrix.rna$key, 
                                   levels = c("rna_fold_nf_wt", 
                                            "rna_fold_nf_cko","rna_fold_nf_dko","rna_fold_nf_dcd",
                                            "rna_fold_nwt_cko","rna_fold_nwt_dko",
                                            "rna_fold_nwt_dcd","rna_fold_fwt_cko","rna_fold_fwt_dko", 
                                            "rna_fold_fwt_dcd"))


ggplot() +
  #geom_violin(data=supermatrix.rna, aes(x=key, y=value, fill=peaktype)) +
  geom_boxplot(data=supermatrix.rna, aes(x=key, y=value, fill=peaktype)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))


##plot lcpm 
rnalist <- list(supermatrix.up_k27a, supermatrix.down_k27a, supermatrix.non_k27a)
rnalist <- list(supermatrix.up, supermatrix.down, supermatrix.non)
rnalist <- listx
plotlist=list()
counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  compare1 <- c("NWT", "FWT", "NDKO", "FDKO")
  #matrix_oi <- matrix_oi[,c("NWT", "FWT", "NDKO", "FDKO")] 
  matrix_oi_gather <- gather(matrix_oi[compare1])
  
  #sample <- samplenames
  
  xyz <- compare_means(value ~ key, data = matrix_oi_gather, method="wilcox.test", p.adjust.method="BH" )
  xyz <- xyz %>% slice(-c(3,4))
  xyz <- xyz %>% mutate(y.position= c(14, 16,18, 20))
  
  plotlist[[counter]] <- print(ggboxplot(data=matrix_oi_gather, x = "key", y="value", notch = TRUE, 
                                         add = "jitter", add.params = list(size = 0.15, color = "black",alpha=0.2), 
                                         color = "key", palette = "aaas",
                                         xlab="", ylab="Log2 CPM") +
                                 stat_pvalue_manual(xyz, label = "p.adj") +
                                 theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none")
  )
  
}
wrap_plots(plotlist)


##plot fold changes
plotlist=list()
counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  compare1 <- c("rna_fold_nf_wt", "rna_fold_nf_dko")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  
  #sample <- samplenames
  
  xyz <- compare_means(value ~ key, data = matrix_oi_gather, method="wilcox.test", p.adjust.method="BH" )
  xyz <- xyz %>% mutate(y.position= c(7))
  
  plotlist[[counter]] <- print(ggboxplot(data=matrix_oi_gather, x = "key", y="value", notch = TRUE, 
                                         add = "jitter", add.params = list(size = 0.05, color = "black",alpha=0.2), 
                                         color = "key", palette = "aaas",
                                         xlab="", ylab="Log2 CPM") +
                                 stat_pvalue_manual(xyz, label = "p.adj") +
                                 theme(legend.position="none")
  )
  
}
wrap_plots(plotlist)




#pk_fold_nf_wt_k27a,pk_fold_nf_dko_k27a,pk_fold_n_dkowt_k27a,pk_fold_f_dkowt_k27a
# test <- as.data.frame(cbind(supermatrix.up$pk_fold_nf_wt_k27a, supermatrix.up$pk_fold_nf_dko_k27a, supermatrix.up$pk_fold_n_dkowt_k27a, supermatrix.up$pk_fold_f_dkowt_k27a, supermatrix.up$rna_fold_nwt_dko))
# colnames(test) <- c("pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a", "pk_fold_n_dkowt_k27a", "pk_fold_f_dkowt_k27a", "rna_fold_nwt_dko")
# test2 <- gather(test, key, value)
# ggplot() +
#   #geom_violin(data=supermatrix.rna, aes(x=key, y=value, fill=peaktype)) +
#   geom_boxplot(data=test2, aes(x=key, y=value)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45))

#####
#####cpm scatter plots wt and dko#####
compare_group_first=c("N_WT_k4m1", "N_DKO_k4m1", "N_WT_k27a", "N_DKO_k27a", "N_WT_k4m3", "N_DKO_k4m3",
                             "N_WT_rad21", "N_DKO_rad21", "N_WT_k27m3","N_DKO_k27m3")
compare_group_firstname=c("N_WT_k4m1", "N_DKO_k4m1", "N_WT_k27a", "N_DKO_k27a", "N_WT_k4m3", "N_DKO_k4m3",
                          "N_WT_rad21", "N_DKO_rad21", "N_WT_k27m3","N_DKO_k27m3")
compare_group_second=c("F_WT_k4m1", "F_DKO_k4m1", "F_WT_k27a", "F_DKO_k27a", "F_WT_k4m3", "F_DKO_k4m3", "F_WT_rad21", "F_DKO_rad21","F_WT_k27m3","F_DKO_k27m3")
compare_group_secondname=c("F_WT_k4m1", "F_DKO_k4m1", "F_WT_k27a", "F_DKO_k27a", "F_WT_k4m3", "F_DKO_k4m3", "F_WT_rad21", "F_DKO_rad21","F_WT_k27m3","F_DKO_k27m3")

compare_group_first=c("N_WT_k4m3", "N_WT_k4m3", "F_WT_k4m3")
compare_group_firstname <- compare_group_first
compare_group_second=c("F_WT_k4m3", "N_DKO_k4m3", "F_DKO_k4m3" )
compare_group_secondname <- compare_group_second

compare_group_first=c("N_WT_k27m3", "N_WT_k27m3", "F_WT_k27m3")
compare_group_firstname <- compare_group_first
compare_group_second=c("F_WT_k27m3", "N_DKO_k27m3", "F_DKO_k27m3" )
compare_group_secondname <- compare_group_second

compare_group_first=c("N_WT_k27a", "N_DKO_k27a")
compare_group_firstname <- compare_group_first
compare_group_second=c("F_WT_k27a", "F_DKO_k27a")
compare_group_secondname <- compare_group_second

compare_group_first=c("F_WT_k4m1", "F_DKO_k4m1")
compare_group_firstname <- compare_group_first
compare_group_second=c("F_WT_k27a", "F_DKO_k27a")
compare_group_secondname <- compare_group_second

compare_group_first=c("N_WT_k4m1", "N_DKO_k4m1")
compare_group_firstname <- compare_group_first
compare_group_second=c("N_WT_k27a", "N_DKO_k27a")
compare_group_secondname <- compare_group_second

compare_group_first=c("N_WT_k4m1")
compare_group_firstname <- compare_group_first
compare_group_second=c("N_WT_k27a")
compare_group_secondname <- compare_group_second


comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
plotlist=list()
plotlist_hist_top=list()
plotlist_hist_right=list()
combinedplotlist=list()
counter = 0

plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.non, aes_string(x=xsample, y=ysample), color='gray', size = 0.1, alpha = 0.05) +
                                                  geom_point(data=supermatrix.up, aes_string(x=xsample, y=ysample), color='red', size = 0.1, alpha = 0.05) +
                                                  geom_point(data=supermatrix.down, aes_string(x=xsample, y=ysample), color='blue', size = 0.1, alpha = 0.05) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(0,15) +
                                                  ylim(0,15))

  combinedplotlist[[counter]] <- plotlist[[counter]]
  
  #combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)


comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
plotlist=list()
combinedplotlist=list()
counter = 0


supermatrix.subcat <- supermatrix.non
plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.subcat, aes_string(x=xsample, y=ysample), color='black', size = 0.1, alpha = 0.2) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(0,10) +
                                                  ylim(0,15))
  
  combinedplotlist[[counter]] <- plotlist[[counter]]
  
  #combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)

cor(supermatrix.subcat$N_WT_k4m1, supermatrix.subcat$N_WT_k27a, method = "spearman")
cor(supermatrix.subcat$N_DKO_k4m1, supermatrix.subcat$N_DKO_k27a, method = "spearman")

cor(supermatrix.subcat$F_WT_k4m1, supermatrix.subcat$F_WT_k27a, method = "spearman")
cor(supermatrix.subcat$F_DKO_k4m1, supermatrix.subcat$F_DKO_k27a, method = "spearman")


cor(x, y = NULL, use = "everything",
    method = c("pearson", "kendall", "spearman"))
#plotlist_hist_top[[5]] + plot_spacer() + plotlist[[5]] + plotlist_hist_right[[5]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))


#####
######fold change scatter plots wt and dko#####

compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_firstname=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_second=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")
compare_group_secondname=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")

compare_group_first=c("pk_fold_nf_wt_k27a", "pk_fold_nf_wt_k27a")
compare_group_firstname <- compare_group_first
compare_group_second=c("pk_fold_n_dkowt_k27a", "pk_fold_f_dkowt_k27a")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_nf_wt_k4m3", "pk_fold_nf_wt_k4m3")
compare_group_firstname <- compare_group_first
compare_group_second=c("pk_fold_n_dkowt_k4m3", "pk_fold_f_dkowt_k4m3")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_nf_wt_k4m3")
compare_group_firstname <- compare_group_first
compare_group_second=c("pk_fold_nf_dko_k4m3")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_nf_wt_k27a")
compare_group_firstname <- compare_group_first
compare_group_second=c("pk_fold_nf_dko_k27a")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_nf_wt_k27m3")
compare_group_firstname <- compare_group_first
compare_group_second=c("pk_fold_nf_dko_k27m3")
compare_group_secondname <- compare_group_second


compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_dko_k4m1")
compare_group_firstname <- compare_group_first
compare_group_second=c("pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_n_dkowt_k4m1", "pk_fold_f_dkowt_k4m1")
compare_group_firstname <- compare_group_first
compare_group_second=c("pk_fold_n_dkowt_k27a", "pk_fold_f_dkowt_k27a")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_n_dkowt_k27a")
compare_group_firstname <- compare_group_first
compare_group_second=c("pk_fold_f_dkowt_k27a")
compare_group_secondname <- compare_group_second




supermatrix.all <- rbind(supermatrix.down, supermatrix.up, supermatrix.non)

comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
plotlist=list()
plotlist_hist_top=list()
plotlist_hist_right=list()
combinedplotlist=list()
counter = 0

###fold change vs fold change
plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  ct <- cor.test(supermatrix.all[,xsample], supermatrix.all[,ysample], method = "spearman", conf.level = 0.95, exact = FALSE) 
  ct <- round(as.numeric(ct[4]), digits = 3)
  grobcorr <- grobTree(textGrob(paste(parse(text="\u03C1"), ct, sep=' = '), x=0.75,  y=0.10, just="left",
                                gp=gpar(col="black", fontsize=8)))
  grobpearson <- grobTree(textGrob("Spearman", x=0.75,  y=0.15, just="left",
                                   gp=gpar(col="black", fontsize=8)))
  # annotation_custom(grobcorr) +
  #   annotation_custom(grobpearson) +
  
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.non, aes_string(x=xsample, y=ysample), color='gray', size = 0.1, alpha = 0.2) +
                                                  geom_point(data=supermatrix.up, aes_string(x=xsample, y=ysample), color='red', size = 0.1, alpha = 0.2) +
                                                  geom_point(data=supermatrix.down, aes_string(x=xsample, y=ysample), color='blue', size = 0.1, alpha = 0.2) +
                                                  xlab(paste(xname, "(Log2 Read Density)")) +
                                                  ylab(paste(yname, "(Log2 Read Density)")) +
                                                  annotation_custom(grobcorr) +
                                                  annotation_custom(grobpearson) +
                                                  #geom_smooth(data=supermatrix.all, aes_string(x=xsample, y=ysample), method='lm', color = "black", size = 0.1) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(-4,4) +
                                                  ylim(-10,10))
  
  plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                           geom_density(data=supermatrix.up, aes_string(x=xsample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.down, aes_string(x=xsample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.non, aes_string(x=xsample), color='black', size = 0.2, alpha = 0.2) +
                                                           theme_void() +
                                                           theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                           xlim(-4,4))
                                         
  
  plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                             geom_density(data=supermatrix.up, aes_string(x=ysample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.down, aes_string(x=ysample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.non, aes_string(x=ysample), color='black', size = 0.2, alpha = 0.2) +
                                                             theme_void() +
                                                             theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                             xlim(-10,10) +
                                                             coord_flip())
  
  combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)

supermatrix.subcat <- supermatrix.up
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
plotlist=list()
combinedplotlist=list()
counter = 0
plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.subcat, aes_string(x=xsample, y=ysample), color='black', size = 0.2, alpha = 0.3) +
                                                  xlab(paste(xname, "(Log2 Read Density)")) +
                                                  ylab(paste(yname, "(Log2 Read Density)")) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(-4,4) +
                                                  ylim(-10,10))
  
  combinedplotlist[[counter]] <- plotlist[[counter]]
  
  
}
wrap_plots(combinedplotlist)




compare_group_first=c("pk_fold_nf_wt_k27a", "pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a","pk_fold_nf_dko_k27a","pk_fold_nf_wt_k27a", "pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a","pk_fold_nf_dko_k27a")
compare_group_firstname <- compare_group_first
compare_group_second=c("N_WT_k27a", "N_DKO_k27a", "N_WT_k27a", "N_DKO_k27a","F_WT_k27a","F_DKO_k27a","F_WT_k27a", "F_DKO_k27a")
compare_group_secondname <- compare_group_second
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))

compare_group_first=c("pk_fold_nf_wt_rad21", "pk_fold_nf_wt_rad21","pk_fold_nf_dko_rad21","pk_fold_nf_dko_rad21","pk_fold_nf_wt_rad21", "pk_fold_nf_wt_rad21","pk_fold_nf_dko_rad21","pk_fold_nf_dko_rad21")
compare_group_firstname <- compare_group_first
compare_group_second=c("N_WT_rad21", "N_DKO_rad21", "N_WT_rad21", "N_DKO_rad21","F_WT_rad21","F_DKO_rad21","F_WT_rad21", "F_DKO_rad21")
compare_group_secondname <- compare_group_second
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))

compare_group_first=c("pk_fold_nf_wt_k4m3", "pk_fold_nf_wt_k4m3","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k4m3","pk_fold_nf_wt_k4m3", "pk_fold_nf_wt_k4m3","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k4m3")
compare_group_firstname <- compare_group_first
compare_group_second=c("N_WT_k4m3", "N_DKO_k4m3", "N_WT_k4m3", "N_DKO_k4m3","F_WT_k4m3","F_DKO_k4m3","F_WT_k4m3", "F_DKO_k4m3")
compare_group_secondname <- compare_group_second
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))

compare_group_first=c("pk_fold_nf_wt_k27m3", "pk_fold_nf_wt_k27m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_k27m3","pk_fold_nf_wt_k27m3", "pk_fold_nf_wt_k27m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_k27m3")
compare_group_firstname <- compare_group_first
compare_group_second=c("N_WT_k27m3", "N_DKO_k27m3", "N_WT_k27m3", "N_DKO_k27m3","F_WT_k27m3","F_DKO_k27m3","F_WT_k27m3", "F_DKO_k27m3")
compare_group_secondname <- compare_group_second
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))

compare_group_first=c("pk_fold_nf_wt_k27a", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k27a", "pk_fold_nf_wt_k27a")
compare_group_firstname <- compare_group_first
compare_group_second=c("N_WT_k27a", "N_DKO_k27a", "F_WT_k27a", "F_DKO_k27a")
compare_group_secondname <- compare_group_second
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))



plotlist=list()
plotlist_hist_top=list()
plotlist_hist_right=list()
combinedplotlist=list()
counter = 0


##foldchange peaks vs readcounts
plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.non, aes_string(x=xsample, y=ysample), color='gray', size = 0.1, alpha = 0.2) +
                                                  geom_point(data=supermatrix.up, aes_string(x=xsample, y=ysample), color='red', size = 0.1, alpha = 0.2) +
                                                  geom_point(data=supermatrix.down, aes_string(x=xsample, y=ysample), color='blue', size = 0.1, alpha = 0.2) +
                                                  xlab(paste(xname, "(Log2 Read Density)")) +
                                                  ylab(paste(yname, "(Log2 Read Density)")) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(-10,10) +
                                                  ylim(0,15))
  
  plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                           geom_density(data=supermatrix.up, aes_string(x=xsample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.down, aes_string(x=xsample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.non, aes_string(x=xsample), color='black', size = 0.2, alpha = 0.2) +
                                                           theme_void() +
                                                           theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                           xlim(-10,10))
  
  
  plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                             geom_density(data=supermatrix.up, aes_string(x=ysample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.down, aes_string(x=ysample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.non, aes_string(x=ysample), color='black', size = 0.2, alpha = 0.2) +
                                                             theme_void() +
                                                             theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                             xlim(0,15) +
                                                             coord_flip())
  
  combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)

#####
######fold change plots with RNA folds#####
#WT
compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_firstname=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_second=c("rna_fold_nf_wt")
compare_group_secondname=c("rna_fold_nf_wt")
# 
# compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
# compare_group_firstname=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
# compare_group_second=c("rna_fold_nf_dko")
# compare_group_secondname=c("rna_fold_nf_dko")
# 
compare_group_first=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")
compare_group_firstname=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")
compare_group_second=c("rna_fold_nf_dko")
compare_group_secondname=c("rna_fold_nf_dko")
# rna_fold_fwt_dko
# compare_group_first=c("rna_fold_nf_wt")
# compare_group_firstname=c("rna_fold_nf_wt")
# compare_group_second=c("rna_fold_nf_dcd")
# compare_group_secondname=c("rna_fold_nf_dcd")

compare_group_first=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")
compare_group_firstname=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")
compare_group_second=c("rna_fold_fwt_dko")
compare_group_secondname=c("rna_fold_fwt_dko")

compare_group_second=c("pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a")
compare_group_secondname <- compare_group_second
#compare_group_firstname=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_first=c("rna_fold_nf_wt","rna_fold_nf_dko")
compare_group_firstname=c("rna_fold_nf_wt","rna_fold_nf_dko")

compare_group_second=c("pk_fold_nf_wt_rad21","pk_fold_nf_dko_rad21")
compare_group_secondname <- compare_group_second
#compare_group_firstname=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_first=c("rna_fold_nf_wt","rna_fold_nf_dko")
compare_group_firstname=c("rna_fold_nf_wt","rna_fold_nf_dko")


compare_group_second=c("pk_fold_nf_wt_k27m3","pk_fold_nf_dko_k27m3")
compare_group_secondname <- compare_group_second
#compare_group_firstname=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_first=c("rna_fold_nf_wt","rna_fold_nf_dko")
compare_group_firstname=c("rna_fold_nf_wt","rna_fold_nf_dko")

compare_group_second=c("pk_fold_nf_wt_k4m3","pk_fold_nf_dko_k4m3")
compare_group_secondname <- compare_group_second
#compare_group_firstname=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_first=c("rna_fold_nf_wt","rna_fold_nf_dko")
compare_group_firstname=c("rna_fold_nf_wt","rna_fold_nf_dko")



# compare_group_first=c("rna_fold_nf_dko","rna_fold_nf_dko")
# compare_group_firstname=c("rna_fold_nf_dko","rna_fold_nf_dko")
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
supermatrix.all <- rbind(supermatrix.down, supermatrix.non, supermatrix.up)

plotlist=list()
plotlist_hist_top=list()
plotlist_hist_right=list()
combinedplotlist=list()
counter = 0

plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  ct <- cor.test(supermatrix.all[,xsample], supermatrix.all[,ysample], method = "spearman", conf.level = 0.95, exact = FALSE) 
  ct <- round(as.numeric(ct[4]), digits = 3)
  grobcorr <- grobTree(textGrob(paste(parse(text="\u03C1"), ct, sep=' = '), x=0.75,  y=0.10, just="left",
                                gp=gpar(col="black", fontsize=8)))
  grobpearson <- grobTree(textGrob("Spearman", x=0.75,  y=0.15, just="left",
                                   gp=gpar(col="black", fontsize=8)))
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.non, aes_string(x=xsample, y=ysample), color='gray', size = 0.1, alpha = 0.4) +
                                                  geom_point(data=supermatrix.up, aes_string(x=xsample, y=ysample), color='red', size = 0.1, alpha = 0.4) +
                                                  geom_point(data=supermatrix.down, aes_string(x=xsample, y=ysample), color='blue', size = 0.1, alpha = 0.4) +
                                                  xlab(paste(xname, "(Log2 Read Density)")) +
                                                  ylab(paste(yname, "(Log2 CPM)")) +
                                                  annotation_custom(grobcorr) +
                                                  annotation_custom(grobpearson) +
                                                  theme_bw() +
                                                  geom_smooth(data=supermatrix.all, aes_string(x=xsample, y=ysample), method='lm', color = "black", size = 0.1) +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(-8,8) +
                                                  ylim(-8,8))
  
  plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                           geom_density(data=supermatrix.up, aes_string(x=xsample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.down, aes_string(x=xsample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.non, aes_string(x=xsample), color='black', size = 0.2, alpha = 0.2) +
                                                           theme_void() +
                                                           theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                           xlim(-8,8))
  
  
  plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                             geom_density(data=supermatrix.up, aes_string(x=ysample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.down, aes_string(x=ysample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.non, aes_string(x=ysample), color='black', size = 0.2, alpha = 0.2) +
                                                             theme_void() +
                                                             theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                             xlim(-8,8) +
                                                             coord_flip())
  
  combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)

#####rna foldchange and lcpm sites
compare_group_first=c("F_WT_k27a","F_DKO_k27a")
compare_group_firstname <- compare_group_first
compare_group_second=c("FWT","FDKO")
compare_group_secondname=c("FWT","FDKO")

compare_group_first=c("N_WT_k27a","N_DKO_k27a")
compare_group_firstname <- compare_group_first
compare_group_second=c("NWT","NDKO")
compare_group_secondname=c("NWT","NDKO")

compare_group_first=c("N_WT_k27m3","N_DKO_k27m3")
compare_group_firstname <- compare_group_first
compare_group_second=c("NWT","NDKO")
compare_group_secondname=c("NWT","NDKO")

compare_group_first=c("N_WT_k4m3","N_DKO_k4m3")
compare_group_firstname <- compare_group_first
compare_group_second=c("NWT","NDKO")
compare_group_secondname=c("NWT","NDKO")

compare_group_first=c("N_WT_rad21","N_DKO_rad21")
compare_group_firstname <- compare_group_first
compare_group_second=c("NWT","NDKO")
compare_group_secondname=c("NWT","NDKO")

# compare_group_first=c("rna_fold_nf_dko","rna_fold_nf_dko")
# compare_group_firstname=c("rna_fold_nf_dko","rna_fold_nf_dko")

comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
supermatrix.all <- rbind(supermatrix.up, supermatrix.non, supermatrix.down)

plotlist=list()
plotlist_hist_top=list()
plotlist_hist_right=list()
combinedplotlist=list()
counter = 0

plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  ct <- cor.test(supermatrix.all[,xsample], supermatrix.all[,ysample], method = "spearman", conf.level = 0.95, exact = FALSE) 
  ct <- round(as.numeric(ct[4]), digits = 3)
  grobcorr <- grobTree(textGrob(paste(parse(text="\u03C1"), ct, sep=' = '), x=0.75,  y=0.10, just="left",
                                gp=gpar(col="black", fontsize=8)))
  #grobpearson <- grobTree(textGrob("Spearman", x=0.75,  y=0.15, just="left",
  #                                 gp=gpar(col="black", fontsize=8)))
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.non, aes_string(x=xsample, y=ysample), color='gray', size = 0.1, alpha = 0.4) +
                                                  geom_point(data=supermatrix.up, aes_string(x=xsample, y=ysample), color='red', size = 0.1, alpha = 0.4) +
                                                  geom_point(data=supermatrix.down, aes_string(x=xsample, y=ysample), color='blue', size = 0.1, alpha = 0.4) +
                                                  xlab(paste(xname, "(Log2 Read Density)")) +
                                                  ylab(paste(yname, "(Log2 CPM)")) +
                                                  annotation_custom(grobcorr) +
                                                  #annotation_custom(grobpearson) +
                                                  theme_bw() +
                                                  geom_smooth(data=supermatrix.all, aes_string(x=xsample, y=ysample), method='glm', color = "black", size = 0.1) +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(0,15) +
                                                  ylim(-2,12))
  
  plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                           geom_density(data=supermatrix.up, aes_string(x=xsample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.down, aes_string(x=xsample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.non, aes_string(x=xsample), color='black', size = 0.2, alpha = 0.2) +
                                                           theme_void() +
                                                           theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                           xlim(0,15))
  
  
  plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                             geom_density(data=supermatrix.up, aes_string(x=ysample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.down, aes_string(x=ysample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.non, aes_string(x=ysample), color='black', size = 0.2, alpha = 0.2) +
                                                             theme_void() +
                                                             theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                             xlim(-2,12) +
                                                             coord_flip())
  
  combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)








####
# compare_group_second=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")
# compare_group_secondname=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")
# 






#compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21",
 #                     "pk_fold_nf_dko_k4m1","pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")

































#####
#####k27a filters including k27a only peaks no k4m1#####
#####k27a only peaks no k4m1
##use bedops --not-element-of 1 to identify k27a peaks that do not have k4m1 peaks
k27apeaks_k4m1filtered <- read.table("~/Desktop/data/k27apeaksonly_nok4m1.txt", sep = '\t', header=FALSE)
colnames(k27apeaks_k4m1filtered) <- c("chr", "start", "end", "peakid")
k27apeaks_k4m1filtered$start <- k27apeaks_k4m1filtered$start + 1
k27apeaksonly_super <- left_join(k27apeaks_k4m1filtered,supermatrix,  by = c("chr" = "Chr", "start" = "Start", "end" = "End"), keep=TRUE)

supermatrix.up_k27a <- k27apeaksonly_super %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))
supermatrix.down_k27a <- k27apeaksonly_super %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))
supermatrix.non_k27a <- k27apeaksonly_super %>% filter(FDR > 0.1 & pk_fold_nf_wt_k27a > -0.7 & pk_fold_nf_wt_k27a < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))

supermatrix.up_k27a <- k27apeaksonly_super %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a > 1 & (simple_annot == "Intergenic" ))
supermatrix.down_k27a <- k27apeaksonly_super %>% filter(FDR < 0.05 & pk_fold_nf_wt_k27a < -1 & (simple_annot == "Intergenic" ))
supermatrix.non_k27a <- k27apeaksonly_super %>% filter(FDR > 0.1 & pk_fold_nf_wt_k27a > -0.7 & pk_fold_nf_wt_k27a < 0.7 & (simple_annot == "Intergenic" ))


# oioio <- supermatrix.up_k27a %>% filter(peakid == 37228)

# k27ak4m1peaks_super <- anti_join(supermatrix, k27apeaks_k4m1filtered,  by = c("Chr" = "chr", "Start" = "start", "End" = "end"), keep=TRUE)
# k27apeaksonly_super <- k27ak4m1peaks_super

generate_deepheatmapgrouped(list(supermatrix.up_k27a, supermatrix.down_k27a, supermatrix.non_k27a), "k4m1k27a_k27aonly.bed")

supermatrix.up_k27a_dep <- supermatrix.up_k27a %>% filter(pk_fold_f_dkowt_k27a < -1)
supermatrix.up_k27a_ind <- supermatrix.up_k27a %>% filter(pk_fold_f_dkowt_k27a < 0.7 & pk_fold_f_dkowt_k27a > -0.7 )
#write.table(supermatrix.up_k27a_dep[1:3], "~/Desktop/xx.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
supermatrix.down_k27a_dep <- supermatrix.down_k27a %>% filter(pk_fold_n_dkowt_k27a < -1)
supermatrix.down_k27a_ind <- supermatrix.down_k27a %>% filter(pk_fold_n_dkowt_k27a < 0.7 & pk_fold_n_dkowt_k27a > -0.7 )

supermatrix.non_k27a_depn <- supermatrix.non_k27a %>% filter(pk_fold_n_dkowt_k27a < -1)
supermatrix.non_k27a_depf <- supermatrix.non_k27a %>% filter(pk_fold_f_dkowt_k27a < -1)
supermatrix.non_k27a_ind <- supermatrix.down_k27a %>% filter((pk_fold_n_dkowt_k27a < 0.7 & pk_fold_n_dkowt_k27a > -0.7) | (pk_fold_f_dkowt_k27a < 0.7 & pk_fold_f_dkowt_k27a > -0.7))

#subset k27a peaks by k4m1 sigup category
k27apeaks_k4m1sigup <- read.table("~/Desktop/data/k27a_at_subbyk4m1sigup.bed", sep = '\t', header=FALSE)
colnames(k27apeaks_k4m1sigup) <- c("chr", "start", "end", "peakid")
k27apeaks_k4m1sigup$start <- k27apeaks_k4m1sigup$start + 1
k27apeaks_k4m1sigup_super <- left_join(k27apeaks_k4m1sigup,supermatrix,  by = c("chr" = "Chr", "start" = "Start", "end" = "End"), keep=TRUE)

k27apeaks_k4m1sigup_k27adown <- k27apeaks_k4m1sigup_super %>% filter(pk_fold_f_dkowt_k27a < -1)
k27apeaks_k4m1sigup_k27anon <- k27apeaks_k4m1sigup_super %>% filter(pk_fold_f_dkowt_k27a < 0.7 & pk_fold_f_dkowt_k27a > -0.7 )
listyy <- list(k27apeaks_k4m1sigup_k27anon, k27apeaks_k4m1sigup_k27adown)

#generate_deepheatmapgrouped(listyy, "k27apeaks_subk4m1sigup.grouped.bed")

##subset k27a peaks by k4m1 nonsig category
k27apeaks_k4m1sigup<- read.table("~/Desktop/data/k27a_at_subbyk4m1nonsig.bed", sep = '\t', header=FALSE)
colnames(k27apeaks_k4m1sigup) <- c("chr", "start", "end", "peakid")
k27apeaks_k4m1sigup$start <- k27apeaks_k4m1sigup$start + 1
k27apeaks_k4m1nonsig_super <- left_join(k27apeaks_k4m1sigup,supermatrix,  by = c("chr" = "Chr", "start" = "Start", "end" = "End"), keep=TRUE)

k27apeaks_k4m1nonsig_k27adown <- k27apeaks_k4m1nonsig_super %>% filter(pk_fold_f_dkowt_k27a < -1 & pk_fold_nf_wt_k27a > 1)

k27apeaks_k4m1nonsig_k27anon <- k27apeaks_k4m1nonsig_super %>% filter(pk_fold_f_dkowt_k27a < 0.7 & pk_fold_f_dkowt_k27a > -0.7 & pk_fold_nf_wt_k27a > 1)
listyyy <- list(k27apeaks_k4m1nonsig_k27anon, k27apeaks_k4m1nonsig_k27adown)
listyyyy <- list(k27apeaks_k4m1nonsig_k27adown, "#")

generate_deepheatmapgrouped(listyyyy, "k27apeaks_subk4m1nonsig_k27aup_dkodown.grouped.bed")





generate_deepheatmapgrouped(rnalist, "k27ak4m1_atac_up_mll34grouped.bed")

generate_deepheatmapgrouped(rnalist, "k27ak4m1_atac_down_mll34grouped.bed")
generate_deepheatmapgrouped(rnalist, "k4m1only_atac__mll34grouped.bed")
generate_deepheatmapgrouped(rnalist, "k27ak4m1_atac_non_mll34grouped.bed")


##nearestTSS analysis
##plot lcpm 
rnalist <- list(supermatrix.up_k27a, supermatrix.down_k27a, supermatrix.non_k27a)
rnalist <- list(supermatrix.up_k27a_dep, supermatrix.up_k27a_ind)
rnalist <- list(supermatrix.down_k27a_dep, supermatrix.down_k27a_ind)
rnalist <- list(supermatrix.non_k27a_depn, supermatrix.non_k27a_depf,supermatrix.non_k27a_ind )
rnalist <- listyyy
plotlist=list()
counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  #compare1 <- c("NWT", "FWT", "NCKO", "FCKO", "NDKO", "FDKO")
  compare1 <- c("NWT", "FWT","NDKO", "FDKO")
  #matrix_oi <- matrix_oi[,c("NWT", "FWT", "NDKO", "FDKO")] 
  matrix_oi_gather <- gather(matrix_oi[compare1])
  
  #sample <- samplenames
  
   xyz <- compare_means(value ~ key, data = matrix_oi_gather, method="wilcox.test", p.adjust.method="BH" )
   xyz <- xyz %>% slice(-c(3,4))
   xyz <- xyz %>% mutate(y.position= c(14, 16,18, 20))
  
  plotlist[[counter]] <- print(ggboxplot(data=matrix_oi_gather, x = "key", y="value", notch = TRUE, 
                                         add = "jitter", add.params = list(size = 0.15, color = "black",alpha=0.2), 
                                         color = "key", palette = "aaas",
                                         xlab="", ylab="Log2 CPM") +
                                 stat_pvalue_manual(xyz, label = "p.adj") +
                                 theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none")
  )
  
}
wrap_plots(plotlist)


##plot fold changes
plotlist=list()
counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  #compare1 <- c("rna_fold_nf_wt", "rna_fold_nf_dko")
  compare1 <- c("rna_fold_nf_wt", "rna_fold_nf_dko")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  
  #sample <- samplenames
  
  xyz <- compare_means(value ~ key, data = matrix_oi_gather, method="wilcox.test", p.adjust.method="BH" )
  xyz <- xyz %>% mutate(y.position= c(7))
  
  plotlist[[counter]] <- print(ggboxplot(data=matrix_oi_gather, x = "key", y="value", notch = TRUE, 
                                         add = "jitter", add.params = list(size = 0.05, color = "black",alpha=0.2), 
                                         color = "key", palette = "aaas",
                                         xlab="", ylab="Log2 CPM") +
                                 stat_pvalue_manual(xyz, label = "p.adj") +
                                 theme(legend.position="none")
  )
  
}
wrap_plots(plotlist)

#####
#####sig gene deseq centric analysis#####

#import deseq genesets
wtsigrna <- read.csv("~/Desktop/data/Diffgenes/DiffgenesNWTvsFWTv3-07.07.20.anno.csv")
# nwtdkorna <- read.csv("~/Desktop/data/Diffgenes/DiffgenesNWTvsNDKOv3-07.07.20.anno.csv")
# wtsigrna <- nwtdkorna
# fwtdkorna <- read.csv("~/Desktop/data/Diffgenes/DiffgenesFWTvsFDKOv3-07.07.20.anno.csv")
# wtsigrna <- fwtdkorna
wtsigrnagenes <- wtsigrna$ensembllistanno


#supermatrix <- subset(supermatrix, neargene.name %in% wtsigrna$ensembllistanno)

supermatrix.wtsigrna <- subset(supermatrix, neargene.name %in% wtsigrna$ensembllistanno)
supermatrix.wtsigrna <- subset(supermatrix_atac_k27a, neargene.name %in% wtsigrna$ensembllistanno)


supermatrix.wtsigrna2 <- supermatrix.wtsigrna %>% filter(simple_annot == "Intron" | simple_annot == "Intergenic")
#supermatrix.wtsigrna2 <- supermatrix.wtsigrna %>% filter(simple_annot == "Intergenic") 
#supermatrix.wtsigrna2 <- supermatrix.wtsigrna %>% filter(simple_annot == "Promoter-TSS") 

supermatrix.wtsigrna2$transitionsrnaratio <- supermatrix.wtsigrna2$rna_fold_nf_dko / supermatrix.wtsigrna2$rna_fold_nf_wt 
supermatrix.wtsigrna3 <- supermatrix.wtsigrna2 %>% filter(transitionsrnaratio < 1.1 | transitionsrnaratio > 0.9)
supermatrix.wtsigrna2 <- supermatrix.wtsigrna3


##for n to f wt sig
supermatrix.wtsigrna.up <- supermatrix.wtsigrna2 %>% filter(rna_fold_nf_wt > 1)
supermatrix.wtsigrna.up_cands <- supermatrix.wtsigrna.up %>% filter(pk_fold_f_dkowt_k27a < -5)
supermatrix.wtsigrna.down <- supermatrix.wtsigrna2 %>% filter(rna_fold_nf_wt < 1)

##for n to n dkowt
supermatrix.wtsigrna.up <- supermatrix.wtsigrna2 %>% filter(rna_fold_nwt_dko > 1)
supermatrix.wtsigrna.down <- supermatrix.wtsigrna2 %>% filter(rna_fold_nwt_dko < 1)

##for f to f dkowt
supermatrix.wtsigrna.up <- supermatrix.wtsigrna2 %>% filter(rna_fold_fwt_dko > 1)
supermatrix.wtsigrna.down <- supermatrix.wtsigrna2 %>% filter(rna_fold_fwt_dko < 1)

xx <- gather(supermatrix.wtsigrna.up[33:44]) 
# xx <- gather(supermatrix.wtsigrna.down[33:52]) 
xx$key <- factor(xx$key, levels = c("N_WT_k4m1", "F_WT_k4m1", "N_DKO_k4m1", "F_DKO_k4m1", "N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a", "N_WT_k4m3", "F_WT_k4m3", "N_DKO_k4m3", "F_DKO_k4m3",
                                    "N_WT_rad21", "F_WT_rad21", "N_DKO_rad21", "F_DKO_rad21", "N_WT_k27m3", "F_WT_k27m3", "N_DKO_k27m3", "F_DKO_k27m3"))



ggplot() +
  #geom_violin(data=supermatrix.rna, aes(x=key, y=value, fill=peaktype)) +
  geom_boxplot(data=xx, aes(x=key, y=value)) +
  geom_jitter(data=xx, aes(x=key, y=value), size = 0.05) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))


xx <- gather(supermatrix.wtsigrna.up[37:40])
xyz <- compare_means(value ~ key, data = xx, method="wilcox.test", p.adjust.method="BH" )
xyz <- xyz %>% slice(-c(3,4))
xyz <- xyz %>% mutate(y.position= c(14, 16,18, 20))

ggboxplot(data=xx, x="key", y="value", notch = TRUE, 
  add = "jitter", add.params = list(size = 0.15, color = "black",alpha=0.2), 
  color = "key", palette = "aaas",
  xlab="", ylab="Log2 Read Density") +
  stat_pvalue_manual(xyz, label = "p.adj") +
  theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none")

xx2 <- gather(supermatrix.wtsigrna.up[,c("pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a")]) 
xx2$key <- factor(xx2$key, levels = c("pk_fold_nf_wt_k4m1","pk_fold_nf_dko_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a",
                                      "pk_fold_nf_wt_k4m3","pk_fold_nf_dko_k4m3","pk_fold_nf_wt_k27m3", "pk_fold_nf_dko_k27m3","pk_fold_nf_wt_rad21", "pk_fold_nf_dko_rad21"))


xyz <- compare_means(value ~ key, data = xx2, method="wilcox.test", p.adjust.method="BH" )
xyz <- xyz %>% mutate(y.position= c(7))

ggboxplot(data=xx2, x="key", y="value", notch = TRUE, 
          add = "jitter", add.params = list(size = 0.15, color = "black",alpha=0.2), 
          color = "key", palette = "aaas",
          xlab="", ylab="Log2 Read Density") +
  stat_pvalue_manual(xyz, label = "p.adj") +
  theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none")



                      
ggplot() +
  #geom_violin(data=supermatrix.rna, aes(x=key, y=value, fill=peaktype)) +
  geom_boxplot(data=xx2, aes(x=key, y=value)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))


# compare_group_first=c("N_WT_k4m1", "N_DKO_k4m1", "N_WT_k27a", "N_DKO_k27a", "N_WT_k4m3", "N_DKO_k4m3",
#                       "N_WT_rad21", "N_DKO_rad21", "N_WT_k27m3","N_DKO_k27m3")
# compare_group_firstname=c("N_WT_k4m1", "N_DKO_k4m1", "N_WT_k27a", "N_DKO_k27a", "N_WT_k4m3", "N_DKO_k4m3",
#                           "N_WT_rad21", "N_DKO_rad21", "N_WT_k27m3","N_DKO_k27m3")
# compare_group_second=c("F_WT_k4m1", "F_DKO_k4m1", "F_WT_k27a", "F_DKO_k27a", "F_WT_k4m3", "F_DKO_k4m3", "F_WT_rad21", "F_DKO_rad21","F_WT_k27m3","F_DKO_k27m3")
# compare_group_secondname=c("F_WT_k4m1", "F_DKO_k4m1", "F_WT_k27a", "F_DKO_k27a", "F_WT_k4m3", "F_DKO_k4m3", "F_WT_rad21", "F_DKO_rad21","F_WT_k27m3","F_DKO_k27m3")
# 
# compare_group_first=c("N_WT_k27a", "N_DKO_k27a")
# compare_group_firstname <- compare_group_first
# compare_group_second=c("F_WT_k27a", "F_DKO_k27a")
# compare_group_secondname <- compare_group_second
# 
# 
# 
# 
# comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
# plotlist=list()
# plotlist_hist_top=list()
# plotlist_hist_right=list()
# combinedplotlist=list()
# counter = 0
# 
# plotarray <- for(i in 1:nrow(comparematrix)) {
#   counter = counter + 1
#   row <- comparematrix[i,]
#   xsample <- as.character(row[1,1])
#   ysample <- as.character(row[1,2])
#   xname <- as.character(row[1,3])
#   yname <- as.character(row[1,4])
#   
#   plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
#                                                   geom_point(data=supermatrix, aes_string(x=xsample, y=ysample), color='gray', size = 0.1, alpha = 0.25) +
#                                                   geom_point(data=supermatrix.wtsigrna.up, aes_string(x=xsample, y=ysample), color='red', size = 0.1, alpha = 0.25) +
#                                                   geom_point(data=supermatrix.wtsigrna.down, aes_string(x=xsample, y=ysample), color='blue', size = 0.1, alpha = 0.25) +
#                                                   xlab(paste(xname, "(Log2 Read Density)")) +
#                                                   ylab(paste(yname, "(Log2 Read Density)")) +
#                                                   theme_bw() +
#                                                   theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#                                                   xlim(0,15) +
#                                                   ylim(0,15))
#   
#   plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
#                                                            geom_density(data=supermatrix, aes_string(x=xsample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
#                                                            geom_density(data=supermatrix.wtsigrna.up, aes_string(x=xsample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
#                                                            geom_density(data=supermatrix.wtsigrna.down, aes_string(x=xsample), color='black', size = 0.2, alpha = 0.2) +
#                                                            theme_void() +
#                                                            theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")))
#   
#   plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
#                                                              geom_density(data=supermatrix, aes_string(x=ysample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
#                                                              geom_density(data=supermatrix.wtsigrna.up, aes_string(x=ysample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
#                                                              geom_density(data=supermatrix.wtsigrna.down, aes_string(x=ysample), color='black', size = 0.2, alpha = 0.2) +
#                                                              theme_void() +
#                                                              theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
#                                                              coord_flip())
#   
#   combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
#   
#   
# }
# wrap_plots(combinedplotlist)

xx <- gather(supermatrix.wtsigrna.down[34:52]) 
xx$key <- factor(xx$key, levels = c("N_WT_k4m1", "F_WT_k4m1", "N_DKO_k4m1", "F_DKO_k4m1", "N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a", "N_WT_k4m3", "F_WT_k4m3", "N_DKO_k4m3", "F_DKO_k4m3",
                                    "N_WT_rad21", "F_WT_rad21", "N_DKO_rad21", "F_DKO_rad21", "N_WT_k27m3", "F_WT_k27m3", "N_DKO_k27m3", "F_DKO_k27m3"))
ggplot() +
  #geom_violin(data=supermatrix.rna, aes(x=key, y=value, fill=peaktype)) +
  geom_boxplot(data=xx, aes(x=key, y=value)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45))



###
yy <- read.table("~/Desktop/k27a_atac_merged.sort.bed")
yy <- yy[1:3]
colnames(yy) <- c("chr", "start", "end")
yy$start <- yy$start + 1

supermatrix_atac_k27a <- inner_join(supermatrix, yy, by = c("Chr" = "chr", "Start" = "start", "End" = "end"), keep=TRUE)




####

#####
#####peak overlap categories fig4b#####

xx <- read.table("~/Desktop/data/k4m1k27asubset/k27apeak_k4m1sub_x1z1.bed")
xx <- read.table("~/Desktop/data/k4m1k27asubset/k27apeak_k4m1sub_x2z2.bed")
xx <- read.table("~/Desktop/data/k4m1k27asubset/k27apeak_k4m1sub_x3z3.bed")
xx <- read.table("~/Desktop/data/k4m1k27asubset/k27apeak_k4m1sub_x4z4.bed")

colnames(xx) <- c("chr", "start", "end")

supermatrix.subset <- inner_join(supermatrix, xx, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))
xxx <- supermatrix.subset %>% filter(pk_fold_n_dkowt_k27a < -1)
xxx2 <- supermatrix.subset %>% filter(pk_fold_n_dkowt_k27a > -1)
#5459, 1698, 1883, 1412, 2283, 1018, 194, 94



xx <- read.table("~/Desktop/data/k4m1k27asubset/k27apeak_k4m1sub_x1y1.bed")
xx <- read.table("~/Desktop/data/k4m1k27asubset/k27apeak_k4m1sub_x2y2.bed")
xx <- read.table("~/Desktop/data/k4m1k27asubset/k27apeak_k4m1sub_x3y3.bed")
xx <- read.table("~/Desktop/data/k4m1k27asubset/k27apeak_k4m1sub_x4y4.bed")

xx <- read.table("~/Desktop/data/k4m1k27asubset/z4.bed")
colnames(xx) <- c("chr", "start", "end")

supermatrix.subset <- inner_join(supermatrix, xx, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

# xxx <- supermatrix.subset %>% filter(pk_fold_f_dkowt_k27a < -1 & pk_fold_n_dkowt_k27a < -1 )
# xxx2 <- supermatrix.subset %>% filter(pk_fold_f_dkowt_k27a > -1 & pk_fold_n_dkowt_k27a > -1)
xxx <- supermatrix.subset %>% filter(pk_fold_n_dkowt_k27a < -1)
xxx2 <- supermatrix.subset %>% filter(pk_fold_n_dkowt_k27a > -1)

nrow(xxx)
nrow(xxx2)
nrow(supermatrix.subset) - nrow(xxx) - nrow(xxx2)

xxx <- supermatrix.subset %>% filter(pk_fold_f_dkowt_k27a < -1)
xxx2 <- supermatrix.subset %>% filter(pk_fold_f_dkowt_k27a > -1)

xxx <- supermatrix.subset %>% filter(pk_fold_n_dkowt_k27a < -1)
xxx2 <- supermatrix.subset %>% filter(pk_fold_n_dkowt_k27a > -1)
#5210, 1995, 967, 607,1062, 334, 1318, 977
generate_deepheatmapgrouped(list(xxx2, xxx), "~/Desktop/x2_v2.bed")


1005, 910
183, 42
81, 68
901, 271
  
844, 1320
973, 168
809, 561
19, 29



#####
#####
#####
