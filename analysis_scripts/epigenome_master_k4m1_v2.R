#Processing and analysis scripts for h3k4me1 data and analysis using h3k4me1 peaks. Primarily Fig.2

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
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggsci_2.9       reshape_0.8.9   svglite_2.1.0   umap_0.2.8.0    forcats_0.5.1   stringr_1.4.0   purrr_0.3.4     readr_2.1.2    
# [9] tidyr_1.2.0     tibble_3.1.7    tidyverse_1.3.1 ggsignif_0.6.3  ggpubr_0.4.0    UpSetR_1.4.0    ggplot2_3.3.6   patchwork_1.1.1
# [17] gridExtra_2.3   dplyr_1.0.9    
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.8.3      lubridate_1.8.0   lattice_0.20-45   png_0.1-7         assertthat_0.2.1  utf8_1.2.2        RSpectra_0.16-1  
# [8] R6_2.5.1          cellranger_1.1.0  plyr_1.8.7        backports_1.4.1   reprex_2.0.1      httr_1.4.3        pillar_1.7.0     
# [15] rlang_1.0.3       readxl_1.4.0      rstudioapi_0.13   car_3.1-0         Matrix_1.4-1      reticulate_1.25   munsell_0.5.0    
# [22] broom_1.0.0       compiler_4.1.2    modelr_0.1.8      systemfonts_1.0.4 pkgconfig_2.0.3   askpass_1.1       openssl_2.0.2    
# [29] tidyselect_1.1.2  fansi_1.0.3       crayon_1.5.1      tzdb_0.3.0        dbplyr_2.2.1      withr_2.5.0       grid_4.1.2       
# [36] jsonlite_1.8.0    gtable_0.3.0      lifecycle_1.0.1   DBI_1.1.3         magrittr_2.0.3    scales_1.2.0      cli_3.3.0        
# [43] stringi_1.7.6     carData_3.0-5     fs_1.5.2          xml2_1.3.3        ellipsis_0.3.2    generics_0.1.3    vctrs_0.4.1      
# [50] tools_4.1.2       glue_1.6.2        hms_1.1.1         abind_1.4-5       colorspace_2.0-3  rstatix_0.7.0     rvest_1.0.2      
# [57] haven_2.5.0  

#####

library(dplyr)
library(gridExtra)
library(patchwork)
library(ggplot2)
library(UpSetR)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(umap)
library(svglite)
library(reshape)
library(ggsci)

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

#####Import deseq2 results and Homer annotation of deseq results for H3K4m1 peaks#####
#prepare deseq2 results, add 1 to start position to be mergecompatiblewith homer annotate peaks file

deseq_wt_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_WT.txt", sep = '\t', header=TRUE)

#deseq_cko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_CKO.txt", sep = '\t', header=TRUE)
# deseq_dko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_DKO.txt", sep = '\t', header=TRUE)
# deseq_dcd_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_dCD.txt", sep = '\t', header=TRUE)
deseq_wt_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_CKOatWT.txt_full.txt", sep = '\t', header=TRUE)
deseq_wt_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_DKOatWT.txt_full.txt", sep = '\t', header=TRUE)
deseq_wt_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_dCDatWT.txt_full.txt", sep = '\t', header=TRUE)


anno_wt_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_WT.anno.txt", sep = '\t', header=TRUE)

#anno_cko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_CKO.anno.txt", sep = '\t', header=TRUE)
# anno_dko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_DKO.anno.txt", sep = '\t', header=TRUE)
# anno_dcd_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_dCD.anno.txt", sep = '\t', header=TRUE)

# readcounts_k4m1_peaks_wt <- read.table('~/Desktop/data/k4m1_diffbind_results/readcounts_k4m1_nf_WT.tab', sep = '\t', header=TRUE)
# readcounts_k4m1_peaks_cko <- read.table('~/Desktop/data/k4m1_diffbind_results/readcounts_k4m1_nf_CKO.tab', sep = '\t', header=TRUE)
#readcounts_k4m1_peaks_dko <- read.table('~/Desktop/data/k4m1_diffbind_results/readcounts_k4m1_nf_DKO.tab', sep = '\t', header=TRUE)
# readcounts_k4m1_peaks_dcd <- read.table('~/Desktop/data/k4m1_diffbind_results/readcounts_k4m1_nf_dCD.tab', sep = '\t', header=TRUE)

#readcounts_k4m1_peaks_wt_atac <- read.table('~/Desktop/data/k4m1_diffbind_results/readcounts_k4m1_nf_WT_v2.tab', sep = '\t', header=TRUE)
readcounts_k4m1_peaks_wt_igg <- read.table('~/Desktop/data/readcounts_k4m1_nf_wt_atigg.tab', sep = '\t', header=TRUE)


# deseq_cko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_CKOatWT.txt_full.txt", sep = '\t', header=TRUE)
# anno_cko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_CKO.anno.txt", sep = '\t', header=TRUE)
# readcounts_k4m1_peaks_cko <- read.table('~/Desktop/data/k4m1_diffbind_results/readcounts_k4m1_nf_CKO.tab', sep = '\t', header=TRUE)
# 
# deseq_dko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_DKO.txt", sep = '\t', header=TRUE)
# anno_dko_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_DKO.anno.txt", sep = '\t', header=TRUE)
# readcounts_k4m1_peaks_dko <- read.table('~/Desktop/data/k4m1_diffbind_results/readcounts_k4m1_nf_DKO.tab', sep = '\t', header=FALSE)
# 
# deseq_dcd_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_dCD.txt", sep = '\t', header=TRUE)
# anno_dcd_k4m1 <- read.csv("~/Desktop/data/k4m1_diffbind_results/diffbindresults_nf_dCD.anno.txt", sep = '\t', header=TRUE)
# readcounts_k4m1_peaks_dcd <- read.table('~/Desktop/data/k4m1_diffbind_results/readcounts_k4m1_nf_dCD.tab', sep = '\t', header=TRUE)




#####
######Integrate data frames######

#choose peak set and corresponding data sets on which to run script
anno_oi <- anno_wt_k4m1
deseq_oi <- deseq_wt_k4m1
readcounts_oi <- readcounts_k4m1_peaks_wt_igg
# 
# anno_oi <- anno_cko_k4m1
# deseq_oi <- deseq_cko_k4m1
# readcounts_oi <- readcounts_k4m1_peaks_cko
# # 
# anno_oi <- anno_dcd_k4m1
# deseq_oi <- deseq_dcd_k4m1
# readcounts_oi <- readcounts_k4m1_peaks_dcd
# 


homerpeakannotationcleanup <- function(mergedannottable) {
  
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
##massage data into mergeable/joinable dfs, start by making uniform labels and first layer processing of primary data
anno_oi <- anno_oi[2:ncol(anno_oi)]
anno_oi2 <- cbind(anno_oi[1:3], anno_oi[9], anno_oi[20:27])
anno_oi2$simple_annot <- t(homerpeakannotationcleanup(anno_oi$Annotation))
anno_oi2$neargene.name <- anno_oi$Gene.Name

deseq_oi$chr <- paste0("chr", deseq_oi$seqnames)
deseq_oi <- cbind(deseq_oi[12], deseq_oi[2:11])
deseq_oi$start <- deseq_oi$start + 1

##for processing cko,dko,dcd genotype specific diffbind data
# colnames(readcounts_oi) <- c("chr", "start", "end", "N_WT_k4m1", "F_WT_k4m1","N_CKO_k4m1", "F_CKO_k4m1", "N_DKO_k4m1", "F_DKO_k4m1","N_dCD_k4m1", "F_dCD_k4m1", "N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a", "N_WT_k4m3", "F_WT_k4m3", "N_DKO_k4m3", "F_DKO_k4m3",
#                               "N_WT_rad21", "F_WT_rad21", "N_DKO_rad21", "F_DKO_rad21", "N_WT_k27m3", "F_WT_k27m3", "N_DKO_k27m3", "F_DKO_k27m3")


##for processing cko,dko,dcd genotype specific diffbind data
# readcounts_oi$width <- readcounts_oi$end - readcounts_oi$start
# readcounts_oi_density <- readcounts_oi[4:27] / readcounts_oi$width
# readcounts_oi_ldensity <- log2(readcounts_oi_density * 1000000 + 1) #multiply by constant, add 1 to avoid divide by zero and log2 scale
# readcounts_oi2 <- cbind(readcounts_oi[1:3], readcounts_oi_ldensity)
# readcounts_oi2$start <- readcounts_oi2$start + 1

#main wt readcounts labels
colnames(readcounts_oi) <- c("chr", "start", "end", "N_WT_k4m1", "F_WT_k4m1","N_CKO_k4m1", "F_CKO_k4m1", "N_DKO_k4m1", "F_DKO_k4m1","N_dCD_k4m1", "F_dCD_k4m1", "N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a", "N_WT_k4m3", "F_WT_k4m3", "N_DKO_k4m3", "F_DKO_k4m3",
                             "N_WT_rad21", "F_WT_rad21", "N_DKO_rad21", "F_DKO_rad21", "N_WT_k27m3", "F_WT_k27m3", "N_DKO_k27m3", "F_DKO_k27m3", "N_WT_AT", "F_WT_AT", "N_DKO_AT", "F_DKO_AT",
                             "N_jfrna_for", "N_jfrna_rev", "F_jfrna_for", "F_jfrna_rev",
                             "CR2_IgG", "N_WT_igg", "F_WT_igg")
#main wt readcounts processing
readcounts_oi$width <- readcounts_oi$end - readcounts_oi$start
readcounts_oi_density <- readcounts_oi[4:38] / readcounts_oi$width
readcounts_oi_ldensity <- log2(readcounts_oi_density * 1000000 + 1) #multiply by constant, add 1 to avoid divide by zero and log2 scale
readcounts_oi2 <- cbind(readcounts_oi[1:3], readcounts_oi_ldensity)
readcounts_oi2$start <- readcounts_oi2$start + 1


######
#####second layer processing for folds and supermatrix dataframe join##########
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

#anno_oi3 <- cbind(anno_oi2, rna_fold_nf_wt, rna_fold_nf_cko,rna_fold_nf_dko,rna_fold_nf_dcd,rna_fold_nwt_cko,rna_fold_nwt_dko,
#                  rna_fold_nwt_dcd,rna_fold_fwt_cko,rna_fold_fwt_dko, rna_fold_fwt_dcd)
anno_oi3 <- cbind(anno_oi2, rna_fold_nf_wt, rna_fold_nf_dko,rna_fold_nwt_dko, rna_fold_fwt_dko)

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

pk_fold_nf_wt_k4m1 <- readcounts_oi2$F_WT_k4m1 - readcounts_oi2$N_WT_k4m1
pk_fold_nf_cko_k4m1 <- readcounts_oi2$F_CKO_k4m1 - readcounts_oi2$N_CKO_k4m1
pk_fold_nf_dko_k4m1 <- readcounts_oi2$F_DKO_k4m1 - readcounts_oi2$N_DKO_k4m1
pk_fold_nf_dcd_k4m1 <- readcounts_oi2$F_dCD_k4m1 - readcounts_oi2$N_dCD_k4m1

##
pk_fold_nf_igg <- readcounts_oi2$F_WT_igg - readcounts_oi2$N_WT_igg

N_WT_k4m1_igg <- readcounts_oi2$N_WT_k4m1 - readcounts_oi2$N_WT_igg
N_CKO_k4m1_igg <- readcounts_oi2$N_CKO_k4m1 - readcounts_oi2$N_WT_igg
N_DKO_k4m1_igg <- readcounts_oi2$N_DKO_k4m1 - readcounts_oi2$N_WT_igg
N_dCD_k4m1_igg <- readcounts_oi2$N_dCD_k4m1 - readcounts_oi2$N_WT_igg

F_WT_k4m1_igg <- readcounts_oi2$F_WT_k4m1 - readcounts_oi2$F_WT_igg
F_CKO_k4m1_igg <- readcounts_oi2$F_CKO_k4m1 - readcounts_oi2$F_WT_igg
F_DKO_k4m1_igg <- readcounts_oi2$F_DKO_k4m1 - readcounts_oi2$F_WT_igg
F_dCD_k4m1_igg <- readcounts_oi2$F_dCD_k4m1 - readcounts_oi2$F_WT_igg
pk_fold_nf_wt_k4m1_igg <- F_WT_k4m1_igg - N_WT_k4m1_igg

N_WT_k4m1_igg2 <- readcounts_oi2$N_WT_k4m1 - readcounts_oi2$CR2_IgG
N_DKO_k4m1_igg2 <- readcounts_oi2$N_DKO_k4m1 - readcounts_oi2$CR2_IgG
F_WT_k4m1_igg2 <- readcounts_oi2$F_WT_k4m1 - readcounts_oi2$CR2_IgG
F_DKO_k4m1_igg2 <- readcounts_oi2$F_DKO_k4m1 - readcounts_oi2$CR2_IgG


pk_fold_nf_wt_at <- readcounts_oi2$F_WT_AT - readcounts_oi2$N_WT_AT 
pk_fold_nf_dko_at <- readcounts_oi2$F_DKO_AT - readcounts_oi2$N_DKO_AT

pk_fold_f_dkowt_at <- readcounts_oi2$F_DKO_AT - readcounts_oi2$F_WT_AT 
pk_fold_n_dkowt_at <- readcounts_oi2$N_DKO_AT - readcounts_oi2$N_WT_AT

##for processing cko,dko,dcd diffbind data
readcounts_oi3 <- cbind(readcounts_oi2, pk_fold_nf_wt_k4m1, pk_fold_nf_wt_k27a,pk_fold_nf_wt_k4m3,pk_fold_nf_wt_k27m3,pk_fold_nf_wt_rad21,
                        pk_fold_nf_dko_k4m1,pk_fold_nf_dko_k27a,pk_fold_nf_dko_k4m3,pk_fold_nf_dko_k27m3,pk_fold_nf_dko_rad21,
                        pk_fold_n_dkowt_k4m1,pk_fold_n_dkowt_k27a,pk_fold_n_dkowt_k4m3,pk_fold_n_dkowt_k27m3,pk_fold_n_dkowt_rad21,
                        pk_fold_f_dkowt_k4m1,pk_fold_f_dkowt_k27a,pk_fold_f_dkowt_k4m3,pk_fold_f_dkowt_k27m3, pk_fold_f_dkowt_rad21,
                        pk_fold_nf_cko_k4m1,pk_fold_nf_dcd_k4m1)


# readcounts_oi3 <- cbind(readcounts_oi2, pk_fold_nf_wt_k4m1, pk_fold_nf_wt_k27a,pk_fold_nf_wt_k4m3,pk_fold_nf_wt_k27m3,pk_fold_nf_wt_rad21,
#                         pk_fold_nf_dko_k4m1,pk_fold_nf_dko_k27a,pk_fold_nf_dko_k4m3,pk_fold_nf_dko_k27m3,pk_fold_nf_dko_rad21,
#                         pk_fold_n_dkowt_k4m1,pk_fold_n_dkowt_k27a,pk_fold_n_dkowt_k4m3,pk_fold_n_dkowt_k27m3,pk_fold_n_dkowt_rad21,
#                         pk_fold_f_dkowt_k4m1,pk_fold_f_dkowt_k27a,pk_fold_f_dkowt_k4m3,pk_fold_f_dkowt_k27m3, pk_fold_f_dkowt_rad21,
#                         pk_fold_nf_cko_k4m1,pk_fold_nf_dcd_k4m1,pk_fold_nf_wt_at, pk_fold_nf_dko_at, pk_fold_f_dkowt_at, pk_fold_n_dkowt_at)

#including igg subtractions
# readcounts_oi3 <- cbind(readcounts_oi2, pk_fold_nf_wt_k4m1, pk_fold_nf_wt_k27a,pk_fold_nf_wt_k4m3,pk_fold_nf_wt_k27m3,pk_fold_nf_wt_rad21,
#                         pk_fold_nf_dko_k4m1,pk_fold_nf_dko_k27a,pk_fold_nf_dko_k4m3,pk_fold_nf_dko_k27m3,pk_fold_nf_dko_rad21,
#                         pk_fold_n_dkowt_k4m1,pk_fold_n_dkowt_k27a,pk_fold_n_dkowt_k4m3,pk_fold_n_dkowt_k27m3,pk_fold_n_dkowt_rad21,
#                         pk_fold_f_dkowt_k4m1,pk_fold_f_dkowt_k27a,pk_fold_f_dkowt_k4m3,pk_fold_f_dkowt_k27m3, pk_fold_f_dkowt_rad21,
#                         pk_fold_nf_cko_k4m1,pk_fold_nf_dcd_k4m1, pk_fold_nf_igg, N_WT_k4m1_igg,N_CKO_k4m1_igg, N_DKO_k4m1_igg, N_dCD_k4m1_igg,
#                         F_WT_k4m1_igg,F_CKO_k4m1_igg, F_DKO_k4m1_igg, F_dCD_k4m1_igg, pk_fold_nf_wt_k4m1_igg,
#                         N_WT_k4m1_igg2, F_WT_k4m1_igg2, N_DKO_k4m1_igg2, F_DKO_k4m1_igg2)


supermatrix1 <- full_join(deseq_oi, anno_oi3, by = c("chr" = "Chr", "start" = "Start", "end" = "End"), keep=TRUE)
supermatrix <- full_join(supermatrix1, readcounts_oi3, by = c("Chr" = "chr", "Start" = "start", "End" = "end"), keep=TRUE)

# supermatrix1 <- full_join(deseq_oi, anno_oi3, by = c("Chr", "Start", "End"), keep=TRUE)
# supermatrix <- full_join(supermatrix1, readcounts_oi3, by = c("Chr" = "chr", "Start" = "start", "End" = "end"), keep=TRUE)


#write.table(supermatrix, "~/Desktop/data/processed_data/supermatrix_k4m1.20220426.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

######
#####skip processing and import super matrix #####
supermatrix <- read.delim("~/Desktop/data/processed_data/supermatrix_k4m1.20220426.txt", sep = '\t', header = TRUE)

######
######subset by co k27a peaks if needed######
k27a_atac_filter <- read.table("~/Desktop/data/k4m1k27a_k4m1peaksbyk27aat.bed", sep = '\t')
colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
supermatrix <- inner_join(supermatrix, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))

k27a_atac_filter <- read.table("~/Desktop/data/k4m1k27a_k4m1peaksnotk27a_at.bed", sep = '\t')
colnames(k27a_atac_filter) <- c("chr", "start", "end", "width")
k27a_atac_filter$start <- k27a_atac_filter$start + 1
supermatrix <- inner_join(supermatrix, k27a_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))




######
######naive or formative plotted individually#####

supermatrix_dkoloss <- supermatrix %>% filter(pk_fold_n_dkowt_k4m1 < -1)
supermatrix_dkoloss  <- supermatrix_dkoloss[order(supermatrix_dkoloss$N_WT_k4m1_igg2),]
supermatrix_dkoloss$peakid <- seq(1,nrow(supermatrix_dkoloss))


ggplot(supermatrix_dkoloss) +
  geom_point(aes(x=peakid, y=N_WT_igg), size = 0.1, color = "gray") +
  geom_point(aes(x=peakid, y=N_DKO_k4m1), color= "red",size = 0.1) +
  geom_point(aes(x=peakid, y=N_WT_k4m1), size = 0.1)

ggplot(supermatrix_dkoloss) +
  geom_point(aes(x=peakid, y=N_DKO_k4m1_igg2), color= "red",size = 0.1) +
  geom_point(aes(x=peakid, y=N_WT_k4m1_igg2), size = 0.1)


supermatrix_dkoloss <- supermatrix %>% filter(pk_fold_f_dkowt_k4m1 < -1) #37712 sites
supermatrix_dkoloss  <- supermatrix_dkoloss[order(supermatrix_dkoloss$F_WT_k4m1_igg2),]
supermatrix_dkoloss$peakid <- seq(1,nrow(supermatrix_dkoloss))
p1 <- ggplot(supermatrix_dkoloss) +
  geom_point(aes(x=peakid, y=F_DKO_k4m1_igg2), color= "magenta3",size = 0.1) +
  geom_point(aes(x=peakid, y=F_WT_k4m1_igg2), size = 0.1) +
  ylim(c(-5,10)) +
  theme_bw()
supermatrix_dkoloss %>% count(F_DKO_k4m1_igg2 < 0 & F_WT_k4m1_igg2 > 0) #18959 sites completely lost
supermatrix_dkoloss %>% count((F_DKO_k4m1_igg2 < F_WT_k4m1_igg2) & F_DKO_k4m1_igg2 > 0 & F_WT_k4m1_igg2 > 0) #15851 sites partial loss


supermatrix_dkoloss <- supermatrix %>% filter(pk_fold_n_dkowt_k4m1 < -1) ##30228 sites
supermatrix_dkoloss  <- supermatrix_dkoloss[order(supermatrix_dkoloss$N_WT_k4m1_igg2),]
supermatrix_dkoloss$peakid <- seq(1,nrow(supermatrix_dkoloss))
p2 <- ggplot(supermatrix_dkoloss) +
  geom_point(aes(x=peakid, y=N_DKO_k4m1_igg2), color= "dodgerblue2",size = 0.1) +
  geom_point(aes(x=peakid, y=N_WT_k4m1_igg2), size = 0.1) +
  ylim(c(-5,10)) +
  theme_bw()

p2 + p1
supermatrix_dkoloss %>% count(N_DKO_k4m1_igg2 < 0 & N_WT_k4m1_igg2 > 0) #17629 sites complete loss
supermatrix_dkoloss %>% count((N_DKO_k4m1_igg2 < N_WT_k4m1_igg2) & N_DKO_k4m1_igg2 > 0 & N_WT_k4m1_igg2 > 0) #8920 sites partial loss

#####
######peak filter parameters: up down nonsig, standard FDR<0.5 1FC, >0.1 FC <0.7#####
# supermatrix.up <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 > 1 )
# supermatrix.down <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 < -1)
# supermatrix.non <- supermatrix %>% filter(FDR > 0.1 & pk_fold_nf_wt_k4m1 > -0.7 & pk_fold_nf_wt_k4m1 < 0.7 )

# supermatrix.up <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 > 1 & simple_annot == "Promoter-TSS")
# supermatrix.down <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 < -1 & simple_annot == "Promoter-TSS")
# supermatrix.non <- supermatrix %>% filter(FDR > 0.1 & pk_fold_nf_wt_k4m1 > -0.7 & Fold < 0.7 & simple_annot == "Promoter-TSS" )

supermatrix.up <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.down <- supermatrix %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
supermatrix.non <- supermatrix %>% filter(FDR > 0.1 & pk_fold_nf_wt_k4m1 > -0.7 & pk_fold_nf_wt_k4m1 < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ) & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))

#generate_deepheatmapgrouped(list(supermatrix.down, supermatrix.up, supermatrix.non), "~/Desktop/k4m1.all.groupedheatmap.bed")
# write.table(supermatrix.up[1:3], "~/Desktop/supermatrix.up.k4m1.ii.subat.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE )
# write.table(supermatrix.down[1:3], "~/Desktop/supermatrix.down.k4m1.ii.subat.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE )
#####
#####filtering for nonchanging site phenotypes#####

supermatrix.k4m1non.loss <- supermatrix.non %>% filter(pk_fold_f_dkowt_k4m1 < -1)
supermatrix.k4m1non.lossv2 <- supermatrix.non %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1)
supermatrix.k4m1non.non <- supermatrix.non %>% filter(pk_fold_f_dkowt_k4m1 > -0.7 & pk_fold_f_dkowt_k4m1 < 0.7)
supermatrix.k4m1non.nonv2 <- supermatrix.non %>% filter(pk_fold_f_dkowt_k4m1 > -0.7 & pk_fold_f_dkowt_k4m1 < 0.7 & pk_fold_n_dkowt_k4m1 > -0.7 & pk_fold_n_dkowt_k4m1 < 0.7)
supermatrix.k4m1non.nonv3 <- supermatrix.non %>% filter(pk_fold_f_dkowt_k4m1 > -0.7 & pk_fold_n_dkowt_k4m1 > -0.7)

#supermatrix.k4m1non.nonv4 <- supermatrix.non %>% filter(pk_fold_f_dkowt_k4m1 > -1 & pk_fold_n_dkowt_k4m1 > -1)



#supermatrix.k4m1non.nonv3 <- supermatrix.non %>% filter(pk_fold_n_dkowt_k4m1 > -0.7 & pk_fold_f_dkowt_k4m1 > -0.7 & (Distance.to.TSS > 2000 | Distance.to.TSS < -2000))
#supermatrix.k4m1non.other <- supermatrix.non %>% filter(pk_fold_f_dkowt_k4m1 < -1 & pk_fold_n_dkowt_k4m1 < -1)
# supermatrix.k4m1non.lossv3 <- supermatrix.non %>% filter(pk_fold_f_dkowt_k4m1 < -1 | pk_fold_n_dkowt_k4m1 < -1)
# supermatrix.k4m1non.nloss <- supermatrix.non %>% filter(pk_fold_f_dkowt_k4m1 > -0.7 & pk_fold_n_dkowt_k4m1 < -1)
# supermatrix.k4m1non.floss <- supermatrix.non %>% filter(pk_fold_n_dkowt_k4m1 > -0.7 & pk_fold_f_dkowt_k4m1 < -1)

#write.table(supermatrix.k4m1.non.non[1:3], "~/Desktop/supermatrix.k4m1non.lossv3.subat.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE )
#write.table(supermatrix.k4m1non.nonv2[1:3], "~/Desktop/supermatrix.k4m1non.nonv2.subat.bed", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE )
generate_deepheatmapgrouped(list(supermatrix.k4m1non.nonv3, supermatrix.k4m1non.lossv2, supermatrix.down, supermatrix.up), "~/Desktop/k4m1.ii.all.groupedheatmapv3.bed")
#####
#####peakannotation barplots#####

##use row counts or use bedops --element-of 1 to find overlaps
df1 <- data.frame(group = c("k4m1non", "k4m1non loss", "down", "up", "totalk4m1non" ),
                  value = c(55862, 8632, 12382, 12982, 91405))

df1$group <- factor(df1$group, levels = c("k4m1non", "k4m1non loss", "down", "up", "totalk4m1non"))

bar1 <- ggplot(df1, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 0.9, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  ylim(0,100000) +
  theme(legend.position = "none")


#total k4m1: non, loss, down, up   30712, 5150, 6859, 6189
#x1 3611, 742, 913, 46 
#y1 2039, 225,149,1172
#z1 2357, 1141, 1370, 488
# 3611 + 2039 + 2357 =8007
# 742 + 225 +1141 = 2108
# 913 + 149 + 1370 = 2432
# 46 + 1172 + 488 = 1706

# df1 <- data.frame(group = c("k4m1non", "k4m1non loss", "down", "up", "x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4", "z1", "z2", "z3", "z4"),
#                   value = c(30712-8007, 5150-2108, 6859-2432, 6189-1706, 3611, 742, 913, 46, 2039, 225,149,1172, 2357, 1141, 1370, 488))
# df1$group <- factor(df1$group, levels = c("k4m1non", "k4m1non loss", "down", "up", "x1", "x2", "x3", "x4","z1", "z2", "z3", "z4", "y1", "y2", "y3", "y4"))

df1 <- data.frame(group = c("x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4", "z1", "z2", "z3", "z4"),
                  value = c(3611, 742, 913, 146, 2039, 225,149,1172, 2357, 1141, 1370, 48))
df1$group <- factor(df1$group, levels = c("x1", "x2", "x3", "x4","z1", "z2", "z3", "z4", "y1", "y2", "y3", "y4"))


bar1 <- ggplot(df1, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 0.9, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none")


###post atac filter barplot
df1 <- data.frame(group = c("k4m1non", "k4m1non loss", "down", "up", "totalk4m1non" ),
                  value = c(30712, 5150, 6859, 6189, 49173))

df1$group <- factor(df1$group, levels = c("k4m1non", "k4m1non loss", "down", "up", "totalk4m1non"))

bar1 <- ggplot(df1, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 0.9, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  ylim(0, 100000) +
  theme(legend.position = "none")

###post k4m1 atac filter barplot

###naive k4m1 k27a co sites and ind dep k27a

df1 <- data.frame(group = c("ind", "dep","ind", "dep","ind", "dep","ind", "dep" ),
                  stackon = c("nonind", "nonind", "nondep", "nondep", "down", "down", "up", "up"),
                  value = c(5459-1698, 1698, 1883-1412, 1412, 2283-1018, 1018, 194-94, 94))

df1$group <- factor(df1$group, levels = c("ind", "dep"))
df1$stackon <- factor(df1$stackon, levels = c("nonind", "nondep", "down", "up"))



bar1 <- ggplot(df1, aes(x=stackon, y=value, fill=group)) +
  geom_bar(width = 0.9, position= "stack",stat = "identity") +
  theme_classic() +
  ylab("") +
  ylim(0,6000) +
  xlab("") +
  theme(legend.position = "none")

prop1 <- ggplot(df1, aes(x=stackon, y=value, fill=group)) +
  geom_bar(width = 0.9, position= "fill",stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none")


###form k4m1 k27a co sites and ind dep k27a

df1 <- data.frame(group = c("ind", "dep","ind", "dep","ind", "dep","ind", "dep" ),
                  stackon = c("nonind", "nonind", "nondep", "nondep", "down", "down", "up", "up"),
                  value = c(5210-1995, 1995, 967-607, 607,1062-334, 334, 1318-977, 977))

df1$group <- factor(df1$group, levels = c("ind", "dep"))
df1$stackon <- factor(df1$stackon, levels = c("nonind", "nondep", "down", "up"))

bar2 <- ggplot(df1, aes(x=stackon, y=value, fill=group)) +
  geom_bar(width = 0.9, stat = "identity") +
  theme_classic() +
  ylab("") +
  ylim(0,6000) +
  xlab("") +
  theme(legend.position = "none")

prop1 <- ggplot(df1, aes(x=stackon, y=value, fill=group)) +
  geom_bar(width = 0.9, position= "fill",stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(legend.position = "none")

bar1 + bar2




#####
#######perform pairwise overlap counts for categories of sig marks and genos#####
# intervene pairwise -i k4m1peaks_wt_sigup.ii.bed k4m1peaks_wt_nonsig.ii.bed k4m1peaks_wt_sigdown.ii.bed k4m1peaks_cko_sigup.ii.bed k4m1peaks_cko_nonsig.ii.bed k4m1peaks_cko_sigdown.ii.bed k4m1peaks_dko_sigup.ii.bed k4m1peaks_dko_nonsig.ii.bed k4m1peaks_dko_sigdown.ii.bed k4m1peaks_dcd_sigup.ii.bed k4m1peaks_dcd_nonsig.ii.bed k4m1peaks_dcd_sigdown.ii.bed -o k4m1_all_pairwise_ii --compute count --triangle=upper --htype=square
# intervene pairwise -i k4m1peaks_wt_sigup.tss.bed k4m1peaks_wt_nonsig.tss.bed k4m1peaks_wt_sigdown.tss.bed k4m1peaks_cko_sigup.tss.bed k4m1peaks_cko_nonsig.tss.bed k4m1peaks_cko_sigdown.tss.bed k4m1peaks_dko_sigup.tss.bed k4m1peaks_dko_nonsig.tss.bed k4m1peaks_dko_sigdown.tss.bed k4m1peaks_dcd_sigup.tss.bed k4m1peaks_dcd_nonsig.tss.bed k4m1peaks_dcd_sigdown.tss.bed -o k4m1_all_pairwise_tss --compute count --triangle=upper --htype=square

#intersection_matrix <- as.matrix(read.table("~/Desktop/data/stratified_peaks_for_intervene_k4m1/k4m1_all_pairwise_tss/Intervene_pairwise_count_matrix.txt"))

intersection_matrix <- as.matrix(read.table("~/Desktop/data/stratified_peaks_for_intervene_k4m1/k4m1_all_pairwise_ii/Intervene_pairwise_count_matrix.txt"))
intersection_matrix2 <- intersection_matrix[4:ncol(intersection_matrix), 1:3]
intersection_matrix3 <- t(intersection_matrix2)
intersection_matrix4 <- reshape::melt.array(intersection_matrix3)
colnames(intersection_matrix4) <- c("wt", "mutant", "count")

intersection_matrix4$wt <- factor(intersection_matrix4$wt, 
                                  levels = c("k4m1peaks_wt_sigdown.ii","k4m1peaks_wt_nonsig.ii","k4m1peaks_wt_sigup.ii" ))
intersection_matrix4$mutant <- factor(intersection_matrix4$mutant, 
                                      levels = c("k4m1peaks_cko_sigdown.ii","k4m1peaks_cko_nonsig.ii","k4m1peaks_cko_sigup.ii",
                                                 "k4m1peaks_dko_sigdown.ii","k4m1peaks_dko_nonsig.ii","k4m1peaks_dko_sigup.ii",
                                                 "k4m1peaks_dcd_sigdown.ii","k4m1peaks_dcd_nonsig.ii","k4m1peaks_dcd_sigup.ii"))    

ggplot(data=intersection_matrix4, aes(x=wt, y = mutant, fill = count)) +
  scale_fill_gradientn(colors = c("white", "deepskyblue4", "deepskyblue4", "deepskyblue4")) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  geom_text(aes(label = count), color = "black", size = 4) +
  coord_flip() +
  theme_void()


##calculate fraction of peaks shared with WT, for ii wtup 
frac_up <- intersection_matrix3[1,] / 13293
frac_non <- intersection_matrix3[2,] / 98645
frac_down <- intersection_matrix3[3,] / 13293

frac_matrix <- rbind(frac_up, frac_non, frac_down)
row.names(frac_matrix) <- c("k4m1peaks_wt_sigup.ii","k4m1peaks_wt_nonsig.ii","k4m1peaks_wt_sigdown.ii" )
frac_matrix <- round(frac_matrix, 2)

intersection_matrix4 <- reshape::melt.array(frac_matrix)
colnames(intersection_matrix4) <- c("wt", "mutant", "count")

intersection_matrix4$wt <- factor(intersection_matrix4$wt, 
                                  levels = c("k4m1peaks_wt_sigdown.ii","k4m1peaks_wt_nonsig.ii","k4m1peaks_wt_sigup.ii" ))
intersection_matrix4$mutant <- factor(intersection_matrix4$mutant, 
                                      levels = c("k4m1peaks_cko_sigdown.ii","k4m1peaks_cko_nonsig.ii","k4m1peaks_cko_sigup.ii",
                                                 "k4m1peaks_dko_sigdown.ii","k4m1peaks_dko_nonsig.ii","k4m1peaks_dko_sigup.ii",
                                                 "k4m1peaks_dcd_sigdown.ii","k4m1peaks_dcd_nonsig.ii","k4m1peaks_dcd_sigup.ii"))    


ggplot(data=intersection_matrix4, aes(x=wt, y = mutant, fill = count)) +
  scale_fill_gradientn(colors = c("white", "deepskyblue4")) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  geom_text(aes(label = count), color = "black", size = 4) +
  coord_flip() +
  theme_void()

#scale_fill_gradient2(limits=c(0, 20000))+


#####
#####k4m1 cpm plots fig3 with igg #####
compare_group_first=c("N_WT_k4m1", "N_CKO_k4m1","N_DKO_k4m1","N_dCD_k4m1", "N_WT_igg")
compare_group_firstname=c()
compare_group_second=c("F_WT_k4m1", "F_CKO_k4m1","F_DKO_k4m1","F_dCD_k4m1", "F_WT_igg")
compare_group_secondname=c()

compare_group_first=c("N_WT_k4m1", "N_CKO_k4m1","N_DKO_k4m1","N_dCD_k4m1")
compare_group_firstname=c()
compare_group_second=c("F_WT_k4m1", "F_CKO_k4m1","F_DKO_k4m1","F_dCD_k4m1")
compare_group_secondname=c()

# compare_group_first=c("N_WT_k4m1_igg", "N_CKO_k4m1_igg","N_DKO_k4m1_igg","N_dCD_k4m1_igg")
# compare_group_firstname=c()
# compare_group_second=c("F_WT_k4m1_igg", "F_CKO_k4m1_igg","F_DKO_k4m1_igg","F_dCD_k4m1_igg")
# compare_group_secondname=c()


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
                                                  xlim(0,10) +
                                                  ylim(0,10))
  # xlab(paste(xname, "(Log2 Read Density)")) +
  #   ylab(paste(yname, "(Log2 Read Density)")) +
  
  # plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
  #                                                          geom_density(data=supermatrix.up, aes_string(x=xsample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
  #                                                          geom_density(data=supermatrix.down, aes_string(x=xsample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
  #                                                          geom_density(data=supermatrix.non, aes_string(x=xsample), color='black', size = 0.2, alpha = 0.2) +
  #                                                          theme_void() +
  #                                                          theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black"))) +
  #   xlim(0,10)
  # 
  
  # plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
  #                                                            geom_density(data=supermatrix.up, aes_string(x=ysample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
  #                                                            geom_density(data=supermatrix.down, aes_string(x=ysample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
  #                                                            geom_density(data=supermatrix.non, aes_string(x=ysample), color='black', size = 0.2, alpha = 0.2) +
  #                                                            theme_void() +
  #                                                            theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
  #                                                            coord_flip()) +
  #   xlim(0,10)
  combinedplotlist[[counter]] <- plotlist[[counter]]
  
  #combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)



#####
#####k4m1 cpm plots fig4 with k27a #####
compare_group_first=c("N_WT_k27a", "N_DKO_k27a","F_WT_k27a","F_DKO_k27a")
compare_group_firstname=c()
compare_group_second=c("F_WT_k27a", "F_DKO_k27a","F_WT_k27a","F_DKO_k27a")
compare_group_secondname=c()

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
                                                  xlim(0,10) +
                                                  ylim(0,10))
  # xlab(paste(xname, "(Log2 Read Density)")) +
  #   ylab(paste(yname, "(Log2 Read Density)")) +
  
  plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                           geom_density(data=supermatrix.up, aes_string(x=xsample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.down, aes_string(x=xsample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.non, aes_string(x=xsample), color='black', size = 0.2, alpha = 0.2) +
                                                           theme_void() +
                                                           theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black"))) +
    xlim(0,10)
  
  
  plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                             geom_density(data=supermatrix.up, aes_string(x=ysample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.down, aes_string(x=ysample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.non, aes_string(x=ysample), color='black', size = 0.2, alpha = 0.2) +
                                                             theme_void() +
                                                             theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                             coord_flip()) +
    xlim(0,10)
  
  combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)


##
compare_group_first=c("N_WT_k27a", "N_DKO_k27a","N_WT_k27a","F_WT_k27a")
compare_group_firstname=c()
compare_group_second=c("F_WT_k27a", "F_DKO_k27a","N_DKO_k27a", "F_DKO_k27a")
compare_group_secondname=c()

comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
plotlist=list()
plotlist_hist_top=list()
plotlist_hist_right=list()
combinedplotlist=list()
counter = 0

##up
plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.up, aes_string(x=xsample, y=ysample), color='black', size = 0.1, alpha = 0.05) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(0,10) +
                                                  ylim(0,10))
  # xlab(paste(xname, "(Log2 Read Density)")) +
  #   ylab(paste(yname, "(Log2 Read Density)")) +
  
  plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                           geom_density(data=supermatrix.up, aes_string(x=xsample), fill= "black", color='red', size = 0.1, alpha = 0.4) +
                                                           theme_void() +
                                                           theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black"))) +
    xlim(0,10)
  
  
  plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                             geom_density(data=supermatrix.up, aes_string(x=ysample), fill= "black", color='red', size = 0.1, alpha = 0.4) +
                                                             theme_void() +
                                                             theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                             coord_flip()) +
    xlim(0,10)
  
  combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)
##down
plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.down, aes_string(x=xsample, y=ysample), color='black', size = 0.1, alpha = 0.05) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(0,10) +
                                                  ylim(0,10))
  # xlab(paste(xname, "(Log2 Read Density)")) +
  #   ylab(paste(yname, "(Log2 Read Density)")) +
  
  plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                           geom_density(data=supermatrix.down, aes_string(x=xsample), fill= "black", color='red', size = 0.1, alpha = 0.4) +
                                                           theme_void() +
                                                           theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black"))) +
    xlim(0,10)
  
  
  plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                             geom_density(data=supermatrix.down, aes_string(x=ysample), fill= "black", color='red', size = 0.1, alpha = 0.4) +
                                                             theme_void() +
                                                             theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                             coord_flip()) +
    xlim(0,10)
  
  combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)





#####
#####peak annotation post filter #####
supermatrix.up$peaktype <- c("Sig up")
supermatrix.down$peaktype <- c("Sig down")
supermatrix.non$peaktype <- c("Non-sig")

peaktype_df <- rbind(supermatrix.up, supermatrix.down, supermatrix.non)

####proportion of annotations underlying peaks
peaktype_df_annot <- as.data.frame(cbind(peaktype_df$peaktype, peaktype_df$simple_annot))
colnames(peaktype_df_annot) <- c("peaktype", "simpleannot")
peaktype_df_annot$peaktype <- factor(peaktype_df_annot$peaktype, 
                                  levels = c("Non-sig", "Sig up", "Sig down"))

peaktypesimplify <- function(df_annot) {
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("Exon", "Other", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("TTS", "Other", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("3'UTR", "Other", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("5'UTR", "Other", x)))
}

peaktype_df_annot2 <- peaktypesimplify(peaktype_df_annot) 
# ggplot(peaktype_df_annot, aes(x=factor(peaktype), fill = factor(simpleannot))) +
#   geom_bar(position = "fill", width = 0.9) +
#   theme_bw() +
#   theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1))
peaktype_df_annot2$simpleannot <- factor(peaktype_df_annot2$simpleannot, levels = c("Other", "Promoter-TSS", "Intergenic", "Intron"))



ggplot(peaktype_df_annot2, aes(x=peaktype, fill = simpleannot)) +
  geom_bar(width = 0.9) +
  theme_classic() +
  scale_fill_npg() +
  theme(text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1))


#####
#####peak fold change plots#####
compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_firstname=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_second=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")
compare_group_secondname=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")


compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k4m1", "pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k4m1")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a", "pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_f_dkowt_at")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_n_dkowt_at")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_nf_wt_at")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_dko_at")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_nf_wt_at")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_f_dkowt_at")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_nf_wt_k4m1")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k27a")
compare_group_secondname <- compare_group_second



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
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.non, aes_string(x=xsample, y=ysample), color='gray', size = 0.1, alpha = 0.2) +
                                                  geom_point(data=supermatrix.up, aes_string(x=xsample, y=ysample), color='red', size = 0.1, alpha = 0.2) +
                                                  geom_point(data=supermatrix.down, aes_string(x=xsample, y=ysample), color='blue', size = 0.1, alpha = 0.2) +
                                                  xlab(paste(xname, "(Log2 Read Density)")) +
                                                  ylab(paste(yname, "(Log2 Read Density)")) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(-8,8) +
                                                  ylim(-10,10))
  
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
                                                             xlim(-10,10) +
                                                             coord_flip())
  
  combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)
#####
#####fold change vs cpm#####

compare_group_first=c("pk_fold_nf_wt_k27a", "pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a","pk_fold_nf_dko_k27a","pk_fold_nf_wt_k27a", "pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a","pk_fold_nf_dko_k27a")
compare_group_firstname <- compare_group_first
compare_group_second=c("N_WT_k27a", "N_DKO_k27a", "N_WT_k27a", "N_DKO_k27a","F_WT_k27a","F_DKO_k27a","F_WT_k27a", "F_DKO_k27a")
compare_group_secondname <- compare_group_second

compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k4m1")
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
#####IgG corrected cpm plots for supps#####

compare_group_first=c("N_WT_k4m1_igg", "N_CKO_k4m1_igg","N_DKO_k4m1_igg","N_dCD_k4m1_igg")
compare_group_firstname=c()
compare_group_second=c("F_WT_k4m1_igg", "F_CKO_k4m1_igg","F_DKO_k4m1_igg","F_dCD_k4m1_igg")
compare_group_secondname=c()


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
                                                  xlim(0,6) +
                                                  ylim(0,6))
  # xlab(paste(xname, "(Log2 Read Density)")) +
  #   ylab(paste(yname, "(Log2 Read Density)")) +
  
  plotlist_hist_top[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                           geom_density(data=supermatrix.up, aes_string(x=xsample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.down, aes_string(x=xsample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                           geom_density(data=supermatrix.non, aes_string(x=xsample), color='black', size = 0.2, alpha = 0.2) +
                                                           theme_void() +
                                                           theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black"))) +
    xlim(0,6)
  
  
  plotlist_hist_right[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                             geom_density(data=supermatrix.up, aes_string(x=ysample), fill= "red", color='red', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.down, aes_string(x=ysample), fill= "blue", color='blue', size = 0.1, alpha = 0.4) +
                                                             geom_density(data=supermatrix.non, aes_string(x=ysample), color='black', size = 0.2, alpha = 0.2) +
                                                             theme_void() +
                                                             theme(legend.position = "none", axis.line = element_line(size = 0.5, colour = "black")) +
                                                             coord_flip()) +
    xlim(0,6)
  
  combinedplotlist[[counter]] <- plotlist_hist_top[[counter]] + plot_spacer() + plotlist[[counter]] + plotlist_hist_right[[counter]] + plot_layout(ncol = 2, nrow = 2, widths = c(7, 1), heights = c(1, 7))
  
  
}
wrap_plots(combinedplotlist)

#####
#####fig2a,e,g foldchange vs FDR plots#####
foldchangedf <- dplyr::select(supermatrix, pk_fold_nf_wt_k4m1)
foldchangedf <- dplyr::select(supermatrix, pk_fold_nf_dko_k4m1)

foldchangedf <- dplyr::select(supermatrix, pk_fold_nf_cko_k4m1)
foldchangedf <- dplyr::select(supermatrix, pk_fold_nf_dcd_k4m1)


foldchangedf$FDR <- supermatrix$FDR
foldchangedf$log10fdr <- -1 * log10(supermatrix$FDR) 
foldchangedf2 <- foldchangedf %>% drop_na()

foldchangedf2sig <- foldchangedf2 %>% filter(log10fdr > 1.3 & (pk_fold_n_dkowt_k4m1 > 1 | pk_fold_n_dkowt_k4m1 < -1))


foldchangedf2sig <- foldchangedf2 %>% filter(log10fdr > 1.3 & (pk_fold_nf_wt_k4m1 > 1 | pk_fold_nf_wt_k4m1 < -1))
foldchangedf2sig <- foldchangedf2 %>% filter(log10fdr > 1.3 & (pk_fold_nf_cko_k4m1 > 1 | pk_fold_nf_cko_k4m1 < -1))
foldchangedf2sig <- foldchangedf2 %>% filter(log10fdr > 1.3 & (pk_fold_nf_dko_k4m1 > 1 | pk_fold_nf_dko_k4m1 < -1))
foldchangedf2sig <- foldchangedf2 %>% filter(log10fdr > 1.3 & (pk_fold_nf_dcd_k4m1 > 1 | pk_fold_nf_dcd_k4m1 < -1))



count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_wt_k4m1 > 1 | pk_fold_nf_wt_k4m1 < -1))) #27458 sites
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_wt_k4m1 > 1))) #13769
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_wt_k4m1 < -1))) #13689


count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_cko_k4m1 > 1 | pk_fold_nf_cko_k4m1 < -1))) #35971 sites #29088 peaks at WT peaks only
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_cko_k4m1 > 1))) #13015
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_cko_k4m1 < -1))) #16073


count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_dko_k4m1 > 1 | pk_fold_nf_dko_k4m1 < -1))) #0 sites

count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_dcd_k4m1 > 1 | pk_fold_nf_dcd_k4m1 < -1))) #429 sites #272 sites at WT peaks only
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_dcd_k4m1 > 1))) #29
count(foldchangedf2 %>% filter(FDR < 0.05 & (pk_fold_nf_dcd_k4m1 < -1))) #243



ggplot() +
  geom_point(data=foldchangedf2, aes(x=pk_fold_nf_wt_k4m1, y=log10fdr), color = "gray", size = 0.3) +
  geom_point(data=foldchangedf2sig, aes(x=pk_fold_nf_wt_k4m1, y=log10fdr), color = "black", size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=1.30103, linetype = "dashed") +
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0,15)


ggplot() +
  geom_point(data=foldchangedf2, aes(x=pk_fold_nf_cko_k4m1, y=log10fdr), color = "gray", size = 0.3) +
  geom_point(data=foldchangedf2sig, aes(x=pk_fold_nf_cko_k4m1, y=log10fdr), color = "black", size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=1.30103, linetype = "dashed") +
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0,15)

ggplot() +
  geom_point(data=foldchangedf2, aes(x=pk_fold_nf_dko_k4m1, y=log10fdr), color = "gray", size = 0.3) +
  #geom_point(data=foldchangedf2sig, aes(x=pk_fold_nf_wt_k4m1, y=log10fdr), color = "black", size = 0.1, alpha = 0.2) +
  theme_bw() +
  geom_hline(yintercept=1.30103, linetype = "dashed") +
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0,4)

ggplot() +
  geom_point(data=foldchangedf2, aes(x=pk_fold_nf_dcd_k4m1, y=log10fdr), color = "gray", size = 0.3) +
  geom_point(data=foldchangedf2sig, aes(x=pk_fold_nf_dcd_k4m1, y=log10fdr), color = "black", size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=1.30103, linetype = "dashed") +
  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  ylim(0,6)






#####
#####fig4 filtering with k27a and atac #####

###use bedops not element of to filter k4m1 peaks by those that overlap with atac/k27a peaks
##filtering by atac overlapping peaks to decrease false positives
k4m1_atac_filter <- read.table("~/Desktop/k4m1wt_atacfilter.bed", sep = '\t')
colnames(k4m1_atac_filter) <- c("chr", "start", "end", "width")
k4m1_atac_filter$start <- k4m1_atac_filter$start + 1
supermatrix <- inner_join(supermatrix, k4m1_atac_filter, by = c("Start" = "start", "Chr" = "chr", "End" = "end"))


##take atac filtered k4m1 peaks and filter by overlap w/wo k27a peaks
k4m1peaks_k27afiltered <- read.table("~/Desktop/data/super_k4m1_peaks_k27afilter.bed", sep = '\t', header=FALSE)
colnames(k4m1peaks_k27afiltered) <- c("chr", "start", "end")
supermatrix_k4m1_k27a <- inner_join(supermatrix, k4m1peaks_k27afiltered, by = c("Chr" = "chr", "Start" = "start", "End" = "end"), keep=TRUE)
supermatrix_k4m1_only <- anti_join(supermatrix, k4m1peaks_k27afiltered, by = c("Chr" = "chr", "Start" = "start", "End" = "end"), keep=TRUE)

supermatrix_k4m1_k27a2 <- supermatrix_k4m1_k27a %>% filter(Distance.to.TSS > 2000 | Distance.to.TSS < -2000)
supermatrix_k4m1_only2 <- supermatrix_k4m1_only %>% filter(Distance.to.TSS > 2000 | Distance.to.TSS < -2000)

supermatrix_k4m1_k27a <- supermatrix_k4m1_k27a2
supermatrix_k4m1_only <- supermatrix_k4m1_only2

##
supermatrix.up_k4m1only <- supermatrix_k4m1_only %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))
supermatrix.down_k4m1only <- supermatrix_k4m1_only %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))
supermatrix.non_k4m1only <- supermatrix_k4m1_only %>% filter(FDR > 0.1 & pk_fold_nf_wt_k4m1 > -0.7 & pk_fold_nf_wt_k4m1 < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))

supermatrix.k4m1up <- supermatrix_k4m1_k27a %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))
supermatrix.k4m1down <- supermatrix_k4m1_k27a %>% filter(FDR < 0.05 & pk_fold_nf_wt_k4m1 < -1 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))
supermatrix.k4m1non <- supermatrix_k4m1_k27a %>% filter(FDR > 0.1 & pk_fold_nf_wt_k4m1 > -0.7 & pk_fold_nf_wt_k4m1 < 0.7 & (simple_annot == "Intron" | simple_annot == "Intergenic" ))

supermatrix.k4m1upk27aup <- supermatrix.k4m1up %>% filter(pk_fold_nf_wt_k27a > 2 )

supermatrix.k4m1up_k27a_ect <- supermatrix.k4m1up %>% filter(pk_fold_f_dkowt_k27a > 1 )
supermatrix.k4m1up_k27a_dep <- supermatrix.k4m1up %>% filter(pk_fold_f_dkowt_k27a < -1 )
supermatrix.k4m1up_k27a_ind <- supermatrix.k4m1up %>% filter(pk_fold_f_dkowt_k27a > -0.7 & pk_fold_f_dkowt_k27a < 0.7 )

supermatrix.k4m1down_k27a_ect <- supermatrix.k4m1down %>% filter(pk_fold_n_dkowt_k27a > 1 )
supermatrix.k4m1down_k27a_dep <- supermatrix.k4m1down %>% filter(pk_fold_n_dkowt_k27a < -1 )
supermatrix.k4m1down_k27a_ind <- supermatrix.k4m1down %>% filter(pk_fold_n_dkowt_k27a > -0.7 & pk_fold_n_dkowt_k27a < 0.7 )

supermatrix.k4m1non_k27a_depn <- supermatrix.k4m1non %>% filter(pk_fold_n_dkowt_k27a < -1 & pk_fold_f_dkowt_k27a > -0.7)
supermatrix.k4m1non_k27a_depf <- supermatrix.k4m1non %>% filter(pk_fold_f_dkowt_k27a < -1 & pk_fold_n_dkowt_k27a > -0.7)
supermatrix.k4m1non_k27a_constloss <- supermatrix.k4m1non %>% filter(pk_fold_f_dkowt_k27a < -1 & pk_fold_n_dkowt_k27a < -1 )
supermatrix.k4m1non_k27a_ind <- supermatrix.k4m1non %>% filter((pk_fold_n_dkowt_k27a > -0.7 & pk_fold_n_dkowt_k27a < 0.7) | (pk_fold_f_dkowt_k27a > -0.7 & pk_fold_f_dkowt_k27a < 0.7))

supermatrix.k4m1non_k27a_indf <- supermatrix.k4m1non %>% filter((pk_fold_f_dkowt_k27a > -0.7 & pk_fold_f_dkowt_k27a < 0.7) & pk_fold_nf_wt_k27a > 1)


listx2 <- list(supermatrix.k4m1non_k27a_indf, supermatrix.k4m1non_k27a_depf, supermatrix.k4m1non_k27a_constloss)
listx3 <- list(supermatrix.k4m1up_k27a_ind, supermatrix.k4m1up_k27a_dep)

listx <- list(supermatrix.k4m1non_k27a_ind, supermatrix.k4m1non_k27a_depf, supermatrix.k4m1non_k27a_constloss, supermatrix.k4m1up_k27a_ind, supermatrix.k4m1up_k27a_dep)
generate_deepheatmapgrouped(listx3, "k4m1k27a_up_fig4grouped.bed")

#count(supermatrix.up_k4m1 %>% drop_na(rna_fold_nf_wt))

generate_deepheatmapgrouped(list(supermatrix.up_k4m1only,  supermatrix.down_k4m1only,supermatrix.non_k4m1only), "k4m1k27a_k4m1onlygrouped.bed")
generate_deepheatmapgrouped(list(supermatrix.k4m1up_k27a_dep, supermatrix.k4m1up_k27a_ind), "k4m1k27a_upgrouped.bed")
generate_deepheatmapgrouped(list(supermatrix.k4m1down_k27a_dep, supermatrix.k4m1down_k27a_ind), "k4m1k27a_downgrouped.bed")
generate_deepheatmapgrouped(list(supermatrix.k4m1non_k27a_depn, supermatrix.k4m1non_k27a_depf, supermatrix.k4m1non_k27a_constloss, supermatrix.k4m1non_k27a_ind), "k4m1k27a_nongroupedv2.bed")

supermatrix.tss <- supermatrix %>% filter(simple_annot == "Promoter-TSS")
supermatrix.ii <- supermatrix %>% filter(simple_annot == "Intron" | simple_annot == "Intergenic" )


##
compare_group_first=c("N_WT_k27a", "N_DKO_k27a","N_WT_k27a","F_WT_k27a")
compare_group_firstname=c()
compare_group_second=c("F_WT_k27a", "F_DKO_k27a","N_DKO_k27a", "F_DKO_k27a")
compare_group_secondname=c()

compare_group_first=c("N_WT_k4m1", "F_WT_k4m1","N_DKO_k4m1", "F_DKO_k4m1","N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a", "N_WT_k4m3", "F_WT_k4m3", "N_DKO_k4m3", "F_DKO_k4m3",
                      "N_WT_rad21", "F_WT_rad21", "N_DKO_rad21", "F_DKO_rad21", "N_WT_k27m3", "F_WT_k27m3", "N_DKO_k27m3", "F_DKO_k27m3", "N_WT_AT", "F_WT_AT", "N_DKO_AT", "F_DKO_AT")
compare_group_firstname=c()
compare_group_second=c("NWT", "FWT", "NDKO", "FDKO")
compare_group_secondname=c()

compare_group_first=c("N_WT_k4m1", "F_WT_k4m1","N_DKO_k4m1", "F_DKO_k4m1","N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a",
                      "N_WT_rad21", "F_WT_rad21", "N_DKO_rad21", "F_DKO_rad21", "N_WT_AT", "F_WT_AT", "N_DKO_AT", "F_DKO_AT")
compare_group_firstname=c()
compare_group_second=c("NWT", "FWT", "NDKO", "FDKO")
compare_group_secondname=c()


compare_group_first=c("F_WT_k4m1","F_DKO_k4m1","F_WT_k27a","F_DKO_k27a",
                     "F_WT_AT", "F_DKO_AT","F_WT_rad21", "F_DKO_rad21")
compare_group_firstname=c()
compare_group_second=c("FWT","FDKO")
compare_group_secondname=c()

compare_group_first=c("N_WT_k4m1","N_DKO_k4m1","N_WT_k27a","N_DKO_k27a",
                      "N_WT_AT", "N_DKO_AT","N_WT_rad21", "N_DKO_rad21")
compare_group_firstname=c()
compare_group_second=c("NWT","NDKO")
compare_group_secondname=c()




compare_group_first=c("N_WT_k4m1", "F_WT_k4m1","N_DKO_k4m1", "F_DKO_k4m1","N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a",
                      "N_WT_AT", "F_WT_AT", "N_DKO_AT", "F_DKO_AT")
compare_group_firstname=c()
compare_group_second=c("NWT", "FWT", "NDKO", "FDKO")
compare_group_secondname=c()


comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
plotlist=list()
combinedplotlist=list()
counter = 0
supermatrix_subcat <- supermatrix_k4m1_k27a
#supermatrix_subcat <- matrix_list[[4]]

plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  ct <- cor.test(supermatrix_subcat[,xsample], supermatrix_subcat[,ysample], method = "spearman", exact=FALSE) #, exact = FALSE
  ct_spearman <- round(as.numeric(ct[4]), digits = 2)
  
  grobcorr <- grobTree(textGrob(paste(parse(text="\u03C1"), ct_spearman, sep=' = '), x=0.5,  y=0.90, just="left",
                                gp=gpar(col="black", fontsize=12)))
  
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix_subcat, aes_string(x=xsample, y=ysample), color='black', size = 0.1, alpha = 0.05) +
                                                  theme_bw() +
                                                  annotation_custom(grobcorr) +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(0,15) +
                                                  ylim(0,15))
  # xlab(paste(xname, "(Log2 Read Density)")) +
  #   ylab(paste(yname, "(Log2 Read Density)")) +
  
  combinedplotlist[[counter]] <- plotlist[[counter]]
  
  
}
wrap_plots(combinedplotlist, ncol = 4, nrow = 2)




####
compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_firstname=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_wt_k4m3","pk_fold_nf_wt_k27m3","pk_fold_nf_wt_rad21")
compare_group_second=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")
compare_group_secondname=c("pk_fold_nf_dko_k4m1", "pk_fold_nf_dko_k27a","pk_fold_nf_dko_k4m3","pk_fold_nf_dko_k27m3","pk_fold_nf_dko_rad21")

compare_group_first=c("pk_fold_nf_wt_rad21", "pk_fold_nf_dko_rad21")

compare_group_first=c("rna_fold_nf_wt", "rna_fold_nf_dko")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_dko_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a", "pk_fold_nf_wt_k4m3","pk_fold_nf_dko_k4m3", "pk_fold_nf_wt_k27m3","pk_fold_nf_dko_k27m3","pk_fold_nf_wt_rad21","pk_fold_nf_dko_rad21", "pk_fold_nf_wt_at", "pk_fold_nf_dko_at")
compare_group_secondname <- compare_group_second 


compare_group_first=c("rna_fold_nf_wt", "rna_fold_nf_dko")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_dko_k4m1", "pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a", "pk_fold_nf_wt_at", "pk_fold_nf_dko_at")
compare_group_secondname <- compare_group_second 

compare_group_second=c("rna_fold_nf_wt", "rna_fold_nf_dko")
compare_group_secondname <-compare_group_first
compare_group_first=c("pk_fold_nf_wt_k27a","pk_fold_nf_dko_k27a", "pk_fold_nf_wt_at", "pk_fold_nf_dko_at", "pk_fold_nf_wt_rad21", "pk_fold_nf_dko_rad21")
compare_group_firstname <- compare_group_second 

compare_group_first=c("pk_fold_n_dkowt_k4m1", "pk_fold_n_dkowt_k27a", "pk_fold_n_dkowt_at", "pk_fold_n_dkowt_rad21")
compare_group_firstname <- compare_group_first 
compare_group_second=c("rna_fold_nwt_dko")
compare_group_secondname <-compare_group_second

compare_group_first=c("pk_fold_f_dkowt_k4m1", "pk_fold_f_dkowt_k27a", "pk_fold_f_dkowt_at", "pk_fold_f_dkowt_rad21")
compare_group_firstname <- compare_group_first 
compare_group_second=c("rna_fold_fwt_dko")
compare_group_secondname <-compare_group_second




comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
plotlist=list()
plotlist_hist_top=list()
plotlist_hist_right=list()
combinedplotlist=list()
counter = 0
supermatrix_subcat <- supermatrix_k4m1_k27a


###fold change vs fold change
plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  ct <- cor.test(supermatrix_subcat[,xsample], supermatrix_subcat[,ysample], method = "spearman", exact=FALSE) #, exact = FALSE
  ct_spearman <- round(as.numeric(ct[4]), digits = 2)
  
  grobcorr <- grobTree(textGrob(paste(parse(text="\u03C1"), ct_spearman, sep=' = '), x=0.7,  y=0.90, just="left",
                                gp=gpar(col="black", fontsize=12)))
  
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix_subcat, aes_string(x=xsample, y=ysample), color='gray', size = 0.2, alpha = 0.7) +
                                                  xlab(paste(xname, "(Log2 Read Density)")) +
                                                  ylab(paste(yname, "(Log2 Read Density)")) +
                                                  theme_bw() +
                                                  annotation_custom(grobcorr) +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
                                                  xlim(-4,4) +
                                                  ylim(-4.5,4.5))
  
  combinedplotlist[[counter]] <- plotlist[[counter]]
  
  
}
wrap_plots(combinedplotlist)


#####
#####fig4supp barplot filtering schematic#####
df1 <- data.frame(group = c("supermatrix_all", "supermatrix_atac", "supermatrix_k4m1only", "supermatrix_atac_k27a"),
           value = c(186824, 102330,48378 ,53952))

df1$group <- factor(df1$group, levels = c("supermatrix_all", "supermatrix_atac", "supermatrix_k4m1only", "supermatrix_atac_k27a"))

bar1 <- ggplot(df1, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(aspect.ratio = 1, legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_npg()

  
df2 <- data.frame(group = c("k4m1only", "k4m1non", "k4m1gain", "k4m1lose"),
                  value = c(48378, 25898 ,4198 , 3056))

df2$group <- factor(df2$group, levels = c("k4m1only", "k4m1non", "k4m1gain", "k4m1lose"))

bar2 <- ggplot(df2, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(aspect.ratio = 1, legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_npg()

df3 <- data.frame(group = c("k4m1k27a_all", "k4m1non", "k4m1gain", "k4m1lose"),
                  value = c(53952, 29117 ,2185 , 3890))

df3$group <- factor(df3$group, levels = c("k4m1k27a_all", "k4m1non", "k4m1gain", "k4m1lose"))

bar3 <- ggplot(df3, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(aspect.ratio = 1, legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_npg()

df5 <- data.frame(group = c("k4m1k27a_lose", "k27a_ect", "k27a_dep", "k27a_ind"),
               value = c(3890, 412,1024 ,1885))
df5$group <- factor(df5$group, levels = c("k4m1k27a_lose", "k27a_ect", "k27a_dep", "k27a_ind"))
bar5 <- ggplot(df5, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(aspect.ratio = 1, legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_npg()

df6 <- data.frame(group = c("k4m1k27a_gain", "k27a_ect", "k27a_dep", "k27a_ind"),
                  value = c(2185, 150, 1227 ,556))
df6$group <- factor(df6$group, levels = c("k4m1k27a_gain", "k27a_ect", "k27a_dep", "k27a_ind"))
bar6 <- ggplot(df6, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(aspect.ratio = 1, legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_npg()


df7 <- data.frame(group = c("k4m1k27a_non", "k27a_formdep", "k27a_naivedep", "k27a_ind"),
                  value = c(29117, 5979,6100 ,19131))
df7$group <- factor(df7$group, levels = c("k4m1k27a_non", "k27a_formdep", "k27a_naivedep", "k27a_ind"))
bar7 <- ggplot(df7, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(aspect.ratio = 1, legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_npg()

#k27a all subset by atac and no k4m1 peak overlap
df8 <- data.frame(group = c("k27a_all", "k27a_atac", "k27a_only"),
                  value = c(112137, 55447,28138))
df8$group <- factor(df8$group, levels = c("k27a_all", "k27a_atac", "k27a_only"))
bar8 <- ggplot(df8, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(aspect.ratio = 1, legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_npg()


#k27a only subset by sigupdownnon
df9 <- data.frame(group = c("k27a_only", "k27a_gain", "k27a_lose", "k27a_non"),
                  value = c(28138, 575, 834, 1652))
df9$group <- factor(df9$group, levels = c("k27a_only", "k27a_gain", "k27a_lose", "k27a_non"))
bar9 <- ggplot(df9, aes(x=group, y=value, fill=group)) +
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ylab("") +
  xlab("") +
  theme(aspect.ratio = 1, legend.position = "none", panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_fill_npg()


dev.new(width = 900, height = 900, unit = "px")
wrap_plots(bar1, bar2, bar3, bar5, bar6, bar7, bar8, bar9, ncol = 3, nrow = 3)
ggsave("peakset_categorization.svg")


#bar1 + bar2










#####
######Nearest TSS on up, down, non primary filter#####

##
rnalist <- list(supermatrix.down_k4m1,supermatrix.up_k4m1,supermatrix.non_k4m1)
rnalist <- list(supermatrix.up_k27a_down, supermatrix.up_k27a_non, supermatrix.up_k27a_up)
rnalist <- list(supermatrix.down_k27a_down, supermatrix.down_k27a_non)
rnalist <- list(supermatrix.non_k27a_downn,supermatrix.non_k27a_downf, supermatrix.non_k27a_non)
rnalist <- list(supermatrix.down_filter, supermatrix.non_filter)
rnalist <- list(supermatrix.up, supermatrix.non, supermatrix.down)
rnalist <- list(supermatrix.k4m1non_k27a_depn, supermatrix.k4m1non_k27a_depf, supermatrix.k4m1non_k27a_constloss, supermatrix.k4m1non_k27a_ind)
rnalist <- list(supermatrix.k4m1down_k27a_ect, supermatrix.k4m1down_k27a_dep, supermatrix.k4m1down_k27a_ind)
samplenames <- c("supermatrix.up_k27a_up", "supermatrix.up_k27a_down", "supermatrix.up_k27a_non")
rnalist <- list(supermatrix.k4m1non_k27a_depf, supermatrix.k4m1non_k27a_constloss, supermatrix.k4m1up_k27a_dep)
rnalist <- list(supermatrix.k4m1non_k27a_constloss, supermatrix.k4m1up_k27a_dep)


length(unique(supermatrix.k4m1up_k27a_dep$neargene.name))

rnalist <- list(supermatrix.up)
rnalist <- list(supermatrix.up.k27anon, supermatrix.up.k27adown)
rnalist <- list(supermatrix.non.k27adown, supermatrix.non.k27anon)
rnalist <- list(supermatrix.down.k27adown, supermatrix.down.k27anon)

##plot lcpm 
plotlist=list()
counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  compare1 <- c("NWT", "FWT", "NDKO", "FDKO")
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

##plotlcpm ggplot
plotlist=list()
counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  matrix_oi <- rnalist[[i]]
  compare1 <- c("NWT", "FWT", "NDKO", "FDKO")
  matrix_oi_gather <- gather(matrix_oi[compare1])
  matrix_oi_gather$key <- factor(matrix_oi_gather$key,levels = c("NWT", "FWT", "NDKO", "FDKO"))
  
  matrix_oi_gather %>% count(!is.na(value))
  #sample <- samplenames
  
  # xyz <- compare_means(value ~ key, data = matrix_oi_gather, method="wilcox.test", p.adjust.method="BH" )
  # xyz <- xyz %>% slice(-c(3,4))
  # xyz <- xyz %>% mutate(y.position= c(14, 16,18, 20))
  
  plotlist[[counter]] <- print(ggplot() +
                                 geom_jitter(data=matrix_oi_gather, aes(x = key, y=value), width=0.2, size = 0.15) +
                                 geom_boxplot(data=matrix_oi_gather, aes(x = key, y=value), width=0.8, outlier.shape=NA) +
                                 theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none") +
                                 theme_classic()
  
                               )
  
  # plotlist[[counter]] <- print(ggboxplot(data=matrix_oi_gather, x = "key", y="value", notch = TRUE, 
  #                                        add = "jitter", add.params = list(size = 0.15, color = "black",alpha=0.2), 
  #                                        color = "key", palette = "aaas",
  #                                        xlab="", ylab="Log2 CPM") +
  #                                stat_pvalue_manual(xyz, label = "p.adj") +
  #                                theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none")
  #)
  
}
wrap_plots(plotlist)

generate_deepheatmapgrouped(rnalist, "~/Desktop/xxx.bed")
# ##unique plot forloop
# plotarray <- for(i in 1:length(rnalist)) {
#   counter = counter + 1
#   matrix_oi <- rnalist[[i]]
#   compare1 <- c("NWT", "FWT", "NDKO", "FDKO")
#   matrix_oi <- distinct(matrix_oi, neargene.name, .keep_all)
#   matrix_oi_gather <- gather(matrix_oi[compare1])
#   matrix_oi_gather$key <- factor(matrix_oi_gather$key,levels = c("NWT", "FWT", "NDKO", "FDKO"))
#   
#   matrix_oi_gather %>% count(!is.na(value))
#   #sample <- samplenames
#   
#   # xyz <- compare_means(value ~ key, data = matrix_oi_gather, method="wilcox.test", p.adjust.method="BH" )
#   # xyz <- xyz %>% slice(-c(3,4))
#   # xyz <- xyz %>% mutate(y.position= c(14, 16,18, 20))
#   
#   plotlist[[counter]] <- print(ggplot() +
#                                  geom_violin(data=matrix_oi_gather, aes(x = key, y=value)) +
#                                  geom_boxplot(data=matrix_oi_gather, aes(x = key, y=value), width=0.5, outlier.shape=NA) +
#                                  theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none") +
#                                  theme_bw()
#                                
#   )
#   
#   # plotlist[[counter]] <- print(ggboxplot(data=matrix_oi_gather, x = "key", y="value", notch = TRUE, 
#   #                                        add = "jitter", add.params = list(size = 0.15, color = "black",alpha=0.2), 
#   #                                        color = "key", palette = "aaas",
#   #                                        xlab="", ylab="Log2 CPM") +
#   #                                stat_pvalue_manual(xyz, label = "p.adj") +
#   #                                theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none")
#   #)
#   
# }
# wrap_plots(plotlist)

supermatrix.up %>% count(is.na(NWT)) #4507 genes
supermatrix.non %>% count(is.na(NWT)) #44356 genes
supermatrix.down %>% count(is.na(NWT)) #5951 genes


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
                                         add = "jitter", add.params = list(size = 0.05, color = "gray",alpha=0.6), 
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
#####Nearest TSS subset by sig rna groups#####

wtsigrna <- read.csv("~/Desktop/data/Diffgenes/DiffgenesNWTvsFWTv3-07.07.20.anno.csv")
wtsigrnagenes <- wtsigrna$ensembllistanno

fdkosigrna <- read.csv("~/Desktop/data/Diffgenes/DiffgenesFWTvsFDKOv3-07.07.20.anno.csv")
fdkosigrna <- fdkosigrna %>% filter(log2FoldChange < -1)
fdkosigrnagenes <- fdkosigrna$ensembllistanno

supermatrix_k4m1_k27a %>% count(pk_fold_f_dkowt_k27a < -1 & pk_fold_nf_wt_k27a < -1)
supermatrix_k4m1_k27a %>% count(pk_fold_nf_wt_k27a < -1)

formloss <- subset(wtsigrna, ensembllistanno %in% fdkosigrnagenes)
supermatrix.wtsigrna <- subset(supermatrix_k4m1_k27a, neargene.name %in% formloss$ensembllistanno)
supermatrix.wtsigrna2 <- supermatrix.wtsigrna %>% filter(rna_fold_nf_wt > 1 & (simple_annot == "Intron" | simple_annot == "Intergenic"))
supermatrix.wtsigrna.up <- supermatrix.wtsigrna2 #%>% filter(rna_fold_nf_wt > 1 & pk_fold_nf_wt_k27a > 1)
# 
# 
supermatrix.wtsigrna <- subset(supermatrix_k4m1_k27a, neargene.name %in% wtsigrna$ensembllistanno)
supermatrix.wtsigrna2 <- supermatrix.wtsigrna %>% filter(simple_annot == "Intron" | simple_annot == "Intergenic")

supermatrix.wtsigrna2$transitionsrnaratio <- supermatrix.wtsigrna2$rna_fold_nf_dko / supermatrix.wtsigrna2$rna_fold_nf_wt 
supermatrix.wtsigrna3 <- supermatrix.wtsigrna2 %>% filter(transitionsrnaratio < 1.2 | transitionsrnaratio > 0.8)
supermatrix.wtsigrna3 <- supermatrix.wtsigrna2 %>% filter((transitionsrnaratio < 1.3 & transitionsrnaratio > 0.7) & rna_fold_nf_wt > 1)
supermatrix.wtsigrna3 <- supermatrix.wtsigrna2 %>% filter((transitionsrnaratio < 2 | transitionsrnaratio > 0) & rna_fold_nf_wt > 1)
supermatrix.wtsigrna3 <- supermatrix.wtsigrna2 %>% filter(rna_fold_nf_wt > 1)
supermatrix.wtsigrna3 <- supermatrix.wtsigrna2 %>% filter(rna_fold_fwt_dko < 0.3 & rna_fold_fwt_dko > -0.3 & rna_fold_nf_wt > 1)



supermatrix.wtsigrna.up <- supermatrix.wtsigrna3 #%>% filter(rna_fold_nf_wt > 1 & pk_fold_nf_wt_k27a > 1)

#supermatrix.wtsigrna3 <- supermatrix.wtsigrna2 %>% filter(rna_fold_fwt_dko < 1.1 | rna_fold_fwt_dko > 0.9 )


supermatrix.wtsigrna.up <- supermatrix.wtsigrna3 #%>% filter(rna_fold_nf_wt > 1 & pk_fold_nf_wt_k27a > 1)


supermatrix.wtsigrna.down <- supermatrix.wtsigrna2 %>% filter(rna_fold_nf_wt < 1) #& pk_fold_nf_wt_k27a < -1)
supermatrix.wtsigrna.up <- supermatrix.wtsigrna2 %>% filter(pk_fold_nf_wt_k27a > 1)
supermatrix.wtsigrna.up <- supermatrix.wtsigrna2 %>% filter(pk_fold_nf_wt_k27a > 1 & pk_fold_nf_wt_k4m1 > 1)


supermatrix.wtsigrna.floss <- supermatrix.wtsigrna2 %>% filter(rna_fold_nf_wt > 1 & pk_fold_f_dkowt_k27a < -1)
supermatrix.wtsigrna.dko <- supermatrix.wtsigrna2 %>% filter(rna_fold_nf_wt > 1 & rna_fold_fwt_dko < 0.1 & rna_fold_fwt_dko > -0.1)


supermatrix.wtsigrna.subcat <- supermatrix.wtsigrna.up

yy <- supermatrix.wtsigrna.subcat %>% filter(pk_fold_f_dkowt_k27a < 0)
length(unique(supermatrix.wtsigrna.subcat$neargene.name)) #unique genes
length(supermatrix.wtsigrna.subcat$neargene.name) #peaks

xx <- supermatrix.wtsigrna.subcat[,c("pk_fold_nf_wt_k4m1", "pk_fold_nf_dko_k4m1", "pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a", "pk_fold_n_dkowt_k27a", "pk_fold_f_dkowt_k27a")]
xxx <- gather(xx)
xxx$key <- factor(xxx$key, levels=c(c("pk_fold_nf_wt_k4m1", "pk_fold_nf_dko_k4m1", "pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a", "pk_fold_n_dkowt_k27a", "pk_fold_f_dkowt_k27a")))

xx <- supermatrix.wtsigrna.subcat[,c("pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a", "pk_fold_n_dkowt_k27a", "pk_fold_f_dkowt_k27a")]
xxx <- gather(xx)
xxx$key <- factor(xxx$key, levels=c(c("pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a", "pk_fold_n_dkowt_k27a", "pk_fold_f_dkowt_k27a")))

xyz <- compare_means(value ~ key, data = xxx, method="wilcox.test", p.adjust.method="BH" )
xyz <- xyz %>% slice(-c(3,4))
xyz <- xyz %>% mutate(y.position= c(14, 16,18, 20))



ggplot() +
  geom_jitter(data=xxx, aes(x=key, y=value), width = 0.2, size = 0.15) +
  #geom_jitter(data=xxx, aes(x=key, y=value), size = 0.05) +
  geom_boxplot(data=xxx, aes(x=key, y=value), width = 0.8, outlier.shape = NA) +
  theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none") +
  theme_classic()

xx <- supermatrix.wtsigrna.subcat[,c("N_WT_k4m1", "F_WT_k4m1", "N_DKO_k4m1", "F_DKO_k4m1","N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a")]
xxx <- gather(xx)
xxx$key <- factor(xxx$key, levels=c("N_WT_k4m1", "F_WT_k4m1", "N_DKO_k4m1", "F_DKO_k4m1","N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a"))

xx <- supermatrix.wtsigrna.subcat[,c("N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a")]
xxx <- gather(xx)
xxx$key <- factor(xxx$key, levels=c("N_WT_k4m1", "F_WT_k4m1", "N_DKO_k4m1", "F_DKO_k4m1","N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a"))


ggplot() +
  geom_jitter(data=xxx, aes(x=key, y=value), width = 0.2, size = 0.15) +
  #geom_jitter(data=xxx, aes(x=key, y=value), size = 0.05) +
  geom_boxplot(data=xxx, aes(x=key, y=value), width = 0.8, outlier.shape = NA) +
  theme(aspect.ratio = 1.5, text = element_text(size = 20), legend.position = "none") +
  theme_classic()


#####
#####UMAP plotting k4m1 peaks#####


#run UMAP package on df_lcpm, batchcorrected if necessary
umap_mark_oi <- supermatrix.non[33:40]
umap_mark_oi <- supermatrix.non[41:44]
umap_mark_oi2 <- umap_mark_oi %>% drop_na()
filteredreplicatelcpmt <- as.data.frame(t(umap_mark_oi2))
filteredreplicatelcpmt2 <- filteredreplicatelcpmt %>% drop_na()

counts.umap <- umap(filteredreplicatelcpmt, n_neighbors = 2) #n_neighbors = 6
head(counts.umap$layout)
sample <- c("N_WT", "F_WT", "N_CKO", "F_CKO", "N_DKO", "F_DKO", "N_dCD", "F_dCD")
#sample <- c("N_WT", "F_WT", "N_DKO", "F_DKO")
replicatename <- sample

layoutdf <- as.data.frame(counts.umap$layout)
colnames(layoutdf) <- c("one", "two")
layoutdf$samplename <- replicatename

layoutdf$samplename <- factor(layoutdf$samplename, 
                              levels = c("N_WT", "F_WT", "N_CKO", "F_CKO", "N_DKO", "F_DKO", "N_dCD", "F_dCD"))  

#300x800pixels
umap_naiveform <- ggplot(layoutdf, mapping = aes(x=one , y=two,color = samplename)) +
  geom_point(size=3) +
  # #geom_text_repel(aes(label = samplename),
  #   box.padding = 0.1, 
  #   point.padding = 0.1,
  #   segment.color = 'grey50') +
  theme_bw() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(aspect.ratio = 1/1.5, text = element_text(size = 20), panel.border = element_rect(color = "black", fill = NA, size = 2))
# xlim(c(-10.5, 10.5)) +
# ylim(c(-3.5,3.5)) 
umap_naiveform

##pca plotting
#center data
umap_mark_oi2$average <- rowMeans(umap_mark_oi2, dims = 1)

pca_mark_oi <- umap_mark_oi2 %>% mutate(across(1:8) - average)
pca_mark_oi <- pca_mark_oi[1:8]

pca_mark <- prcomp(pca_mark_oi, scale. = TRUE, center = TRUE)
rotations <- as.data.frame(pca_mark$rotation)
rotations$samplenames <- rownames(rotations)
rotations$samplenames <- c("N_WT", "F_WT", "N_CKO", "F_CKO", "N_DKO", "F_DKO", "N_dCD", "F_dCD")

percentage <- round(pca_mark$sdev / sum(pca_mark$sdev) * 100, 2)
percentage <- paste(colnames(rotations), "(", paste(as.character(percentage), "%", " )", sep="") )

ggplot(rotations, aes(x=PC2, y=PC3, color=samplenames))+
  geom_point() +
  ylim(c(-0.6, 0.6)) +
  xlab(percentage[2]) +
  ylab(percentage[3])




#####
##### peakfoldchange by rna fold change in peakset categories ####


#####
#####peak fold change plots#####


compare_group_first=c("rna_fold_nf_wt", "rna_fold_nf_dko")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_dko_k4m1")
compare_group_secondname <- compare_group_second
supermatrix_subcat <- supermatrix.k4m1up_k27a_dep

compare_group_first=c("rna_fold_nf_wt", "rna_fold_nf_dko")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k27a", "pk_fold_f_dkowt_k27a")
compare_group_secondname <- compare_group_second
supermatrix_subcat <- supermatrix.up_k27a_down

compare_group_first=c("rna_fold_nf_wt", "rna_fold_fwt_dko")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k27a", "pk_fold_f_dkowt_k27a")
compare_group_secondname <- compare_group_second
supermatrix_subcat <- supermatrix.up_k27a_down

compare_group_first=c("FWT")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k27a")
compare_group_secondname <- compare_group_second
supermatrix_subcat <- supermatrix.up_k27a_down



compare_group_first=c("rna_fold_nf_wt", "rna_fold_nf_dko")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a" )
compare_group_secondname <- compare_group_second
supermatrix_subcat <- supermatrix.up_k27a_down

compare_group_first=c("rna_fold_nf_wt", "rna_fold_nf_dko")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_at", "pk_fold_nf_dko_at")
compare_group_secondname <- compare_group_second
supermatrix_subcat <- supermatrix.up_k27a

compare_group_first=c("rna_fold_nf_wt", "rna_fold_nf_dko")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_at", "pk_fold_nf_dko_at")
compare_group_secondname <- compare_group_second

compare_group_first=c("NWT", "FWT","NDKO", "FDKO" )
compare_group_firstname <-compare_group_first
compare_group_second=c("N_WT_AT", "F_WT_AT")
compare_group_second=c("N_WT_k27a", "F_WT_k27a")
compare_group_second=c("N_WT_k27a", "F_WT_k27a", "N_DKO_k27a", "F_DKO_k27a")
compare_group_secondname <- compare_group_second


compare_group_first=c("pk_fold_nf_wt_k4m1", "pk_fold_nf_dko_k4m1")
compare_group_firstname <-compare_group_first
compare_group_second=c("pk_fold_nf_wt_k27a", "pk_fold_nf_dko_k27a")
compare_group_secondname <- compare_group_second

#supermatrix_subcat <- supermatrix
comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
plotlist=list()
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
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  geom_point(data=supermatrix.non, aes_string(x=xsample, y=ysample), color='gray', size = 0.1, alpha = 0.3) +
                                                  geom_point(data=supermatrix.up, aes_string(x=xsample, y=ysample), color='red', size = 0.1, alpha = 0.3) +
                                                  geom_point(data=supermatrix.down, aes_string(x=xsample, y=ysample), color='blue', size = 0.1, alpha = 0.3) +
                                                  #geom_point(data=supermatrix_subcat, aes_string(x=xsample, y=ysample), color='black', size = 0.1) +
                                                  #geom_smooth(data=supermatrix_subcat, aes_string(x=xsample, y=ysample),method='lm') +
                                                  xlab(paste(xname, "(Log2 Read Density)")) +
                                                  ylab(paste(yname, "(Log2 Read Density)")) +
                                                  theme_bw() +
                                                  xlim(-4,4) +
                                                  ylim(-10,10) +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)))
  #xlim(-8,8) +
  #ylim(-10,10))
  
  combinedplotlist[[counter]] <- plotlist[[counter]]
}
wrap_plots(combinedplotlist)

supermatrixdko.up <- supermatrix.up %>% filter(pk_fold_nf_dko_k27a > 1)



comparematrix <- as.data.frame(cbind(compare_group_first, compare_group_second, compare_group_firstname, compare_group_secondname))
plotlist=list()
combinedplotlist=list()
counter = 0
supermatrix_subcat <- supermatrix.k4m1upk27aup
supermatrix_subcat <- supermatrix.k4m1up_k27a_dep
supermatrix_subcat <- supermatrix.k4m1down_k27a_dep

###fold change vs fold change
plotarray <- for(i in 1:nrow(comparematrix)) {
  counter = counter + 1
  row <- comparematrix[i,]
  xsample <- as.character(row[1,1])
  ysample <- as.character(row[1,2])
  xname <- as.character(row[1,3])
  yname <- as.character(row[1,4])
  
  plotlist[[paste0(xsample, ysample)]] <- print(ggplot() +
                                                  # geom_point(data=supermatrix.non, aes_string(x=xsample, y=ysample), color='gray', size = 0.1, alpha = 0.2) +
                                                  # geom_point(data=supermatrix.up, aes_string(x=xsample, y=ysample), color='red', size = 0.1, alpha = 0.2) +
                                                  # geom_point(data=supermatrix.down, aes_string(x=xsample, y=ysample), color='blue', size = 0.1, alpha = 0.2) +
                                                  geom_point(data=supermatrix_subcat, aes_string(x=xsample, y=ysample), color='black', size = 0.1) +
                                                  geom_smooth(data=supermatrix_subcat, aes_string(x=xsample, y=ysample),method='lm') +
                                                  xlab(paste(xname, "(Log2 Read Density)")) +
                                                  ylab(paste(yname, "(Log2 Read Density)")) +
                                                  xlim(0,10) +
                                                  ylim(0,10) +
                                                  theme_bw() +
                                                  theme(aspect.ratio = 1, text = element_text(size = 15), panel.border = element_rect(color = "black", fill = NA, size = 1)))
  #xlim(-8,8) +
  #ylim(-10,10))
  
  combinedplotlist[[counter]] <- plotlist[[counter]]
}
wrap_plots(combinedplotlist)


