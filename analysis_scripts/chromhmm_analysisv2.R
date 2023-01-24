#chromhmm state processing and plotting for fig.5A and 5C

library(dplyr)
library(tidyverse)
library(ggpubr)
library(patchwork)

#####sessionInfo#####

> sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] patchwork_1.1.1 ggpubr_0.4.0    forcats_0.5.1   stringr_1.4.0   purrr_0.3.4     readr_2.1.2     tidyr_1.2.0     tibble_3.1.7   
[9] ggplot2_3.3.6   tidyverse_1.3.1 dplyr_1.0.9    

loaded via a namespace (and not attached):
  [1] cellranger_1.1.0 pillar_1.7.0     compiler_4.1.2   dbplyr_2.2.1     tools_4.1.2      jsonlite_1.8.0   lubridate_1.8.0  lifecycle_1.0.1 
[9] gtable_0.3.0     pkgconfig_2.0.3  rlang_1.0.3      reprex_2.0.1     DBI_1.1.3        cli_3.3.0        rstudioapi_0.13  haven_2.5.0     
[17] xml2_1.3.3       withr_2.5.0      httr_1.4.3       fs_1.5.2         generics_0.1.3   vctrs_0.4.1      hms_1.1.1        grid_4.1.2      
[25] tidyselect_1.1.2 glue_1.6.2       R6_2.5.1         rstatix_0.7.0    fansi_1.0.3      readxl_1.4.0     carData_3.0-5    car_3.1-0       
[33] tzdb_0.3.0       modelr_0.1.8     magrittr_2.0.3   backports_1.4.1  scales_1.2.0     ellipsis_0.3.2   abind_1.4-5      rvest_1.0.2     
[41] assertthat_0.2.1 colorspace_2.0-3 ggsignif_0.6.3   utf8_1.2.2       stringi_1.7.6    munsell_0.5.0    broom_1.0.0      crayon_1.5.1 

#####

homerpeakannotationcleanup_short <- function(mergedannottable) {
  
  df_annot <- as.data.frame(lapply(mergedannottable, function(x) gsub("intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("non-coding.*Intron.*", "Intron", x)))
  df_annot <- as.data.frame(lapply(df_annot, function(x) gsub("promoter.*", "Promoter-TSS", x)))
  
}


transitions_bed <- read.table('~/Desktop/chromhmm/transitional_16_dense.atacsub75p.bed', sep='\t', header = FALSE)
colnames(transitions_bed) <- c("chr", "start", "end", "state")
transitions_bed <- transitions_bed[1:4]


transitions_bed_homerformat <- transitions_bed[1:3]
transitions_bed_homerformat$peakid <- seq(1,nrow(transitions_bed_homerformat))
#write.table(transitions_bed_homerformat, "transitional_16_dense.atacsub75p.homer.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

anno_oi <- read.delim('~/Desktop/chromhmm/transitional_16_dense.atacsub75p.anno.bed', sep='\t', header = TRUE)
anno_oi <- anno_oi[2:ncol(anno_oi)]
anno_oi2 <- cbind(anno_oi[1:3], anno_oi[7], anno_oi[9], anno_oi[19:27])
anno_oi2$simple_annot <- t(homerpeakannotationcleanup_short(anno_oi$Annotation))
anno_oi2$neargene.name <- anno_oi$Gene.Name
anno_oi2$Start <- anno_oi2$Start - 1

rna_fold_nf_wt <- anno_oi2$FWT - anno_oi2$NWT
rna_fold_nf_dko <- anno_oi2$FDKO - anno_oi2$NDKO

rna_fold_nwt_dko <- anno_oi2$NDKO - anno_oi2$NWT
rna_fold_fwt_dko <- anno_oi2$FDKO - anno_oi2$FWT


anno_oi3 <- cbind(anno_oi2, rna_fold_nf_wt, rna_fold_nf_dko,rna_fold_nwt_dko, rna_fold_fwt_dko)
super_chromhmm <- full_join(transitions_bed, anno_oi3, by = c("chr" = "Chr", "start" = "Start", "end" = "End"), keep=TRUE)
super_chromhmm <- super_chromhmm %>% filter(simple_annot == "Intergenic" | simple_annot == "Intron" | simple_annot == "Promoter-TSS") 

# 
barplot(table(super_chromhmm$simple_annot))
barplot(table(super_chromhmm$state))
# barplot(table(super_chromhmm.tss$state))

###filter for intergenic, intronic vs tss then filter out state 19 "empty" 
##remove empty state 19, convert state # to state_01 if neccessito, generate grouped heatmap
generate_deepheatmapgrouped_usecolumn <- function(listofinputs, outputfilename, dataframe){
  "input a list of bed file coordinates chr,start,end and output a single txt file with # separating groups for use with deeptools plot heatmap"
  "perform function on bed entries with additional criteria from dataframe, e.g. chromhmm state"
  stuffer <- "#"
  df_bed <- dataframe %>% filter(state_name == "state_01")
  write.table(df_bed[1:3], outputfilename, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
  
  for(i in 2:length(listofinputs)){
    next_df_bed <- dataframe %>% filter(state_name == listofinputs[[i]])
    write.table(stuffer, outputfilename, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
    write.table(next_df_bed[1:3], outputfilename, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t') 
    
  }
}
super_chromhmm.ii <- super_chromhmm %>% filter(simple_annot == "Intergenic" | simple_annot == "Intron")
super_chromhmm.tss <- super_chromhmm %>% filter(simple_annot == "Promoter-TSS")


transitions_bed2 <- super_chromhmm.tss
transitions_bed2$state_name <- ifelse(nchar(transitions_bed2$state)!=2,paste0("0",transitions_bed2$state),transitions_bed2$state)
transitions_bed2$state_name <- paste0("state_", transitions_bed2$state_name)

transitions_bed3 <- filter(transitions_bed2, state_name != "state_15")
transitions_bed3_filter <- transitions_bed3 %>% group_by(state_name) %>% filter(n()>= 100) %>% ungroup()

state_array <- c(paste0("state_0", seq(1,9)), paste0("state_", seq(10,14)), "state_16")
#state_array <- c("state_02", "state_11", "state_12", "state_13", "state_14", "state_15", "state_16", "state_23", "state_24")

#generate_deepheatmapgrouped_usecolumn(state_array, "chormhmm_at3subset.formonly.heatmapgrouped.tss.bed", transitions_bed3_filter)


transitions_bed2 <- filter(super_chromhmm.ii, state != 15)
transitions_bed2$state_name <- ifelse(nchar(transitions_bed2$state)!=2,paste0("0",transitions_bed2$state),transitions_bed2$state)
transitions_bed2$state_name <- paste0("state_", transitions_bed2$state_name)
transitions_bed3 <- transitions_bed2 %>% filter(Distance.to.TSS > 2000 | Distance.to.TSS < -2000)
transitions_bed3_filter <- transitions_bed3 %>% group_by(state_name) %>% filter(n()>= 100) %>% ungroup()

state_array <- c(paste0("state_0", seq(1,9)), paste0("state_", seq(10,18)), paste0("state_", seq(20,25)))

#generate_deepheatmapgrouped_usecolumn(state_array, "chormhmm_at3subset.heatmapgrouped.ii.bed", transitions_bed3_filter)


super_chromhmm.ii2 <- filter(super_chromhmm.ii, state != 15)
super_chromhmm.ii2$state_name <- ifelse(nchar(super_chromhmm.ii2$state)!=2,paste0("0",super_chromhmm.ii2$state),super_chromhmm.ii2$state)
super_chromhmm.ii2$state_name <- paste0("state_", super_chromhmm.ii2$state_name)


super_chromhmm.ii_nn <- super_chromhmm.ii2[c("NWT", "FWT", "NDKO", "FDKO", "state_name")]

super_chromhmm.ii_nn <- transitions_bed3_filter[c("NWT", "NDKO","state_name")]
xx$name <- factor(xx$name, levels = c("NWT", "NDKO"))
xx <- super_chromhmm.ii_nn %>% pivot_longer(cols=NWT:NDKO)

super_chromhmm.ii_nn <- transitions_bed3_filter[c("FWT", "FDKO","state_name")]
xx <- super_chromhmm.ii_nn %>% pivot_longer(cols=FWT:FDKO)
xx$name <- factor(xx$name, levels = c("FWT", "FDKO"))

super_chromhmm.ii_nn <- super_chromhmm.ii2[c("CpG.","GC." ,"state_name")]
xx <- super_chromhmm.ii_nn %>% pivot_longer(cols="CpG.":"GC.")

super_chromhmm.ii_nn <- transitions_bed3_filter[c("Distance.to.TSS" ,"state_name")]
xx <- super_chromhmm.ii_nn %>% pivot_longer(cols="Distance.to.TSS")

xx$name <- factor(xx$name, levels = c("FWT", "FDKO"))

#super_chromhmm.ii_nn <- transitions_bed3_filter[c("rna_fold_nf_wt", "rna_fold_nf_dko", "state_name")]
super_chromhmm.ii_nn <- transitions_bed3_filter[c("rna_fold_fwt_dko","state_name")]
# super_chromhmm.ii_nn <- super_chromhmm.ii2[c("rna_fold_nwt_dko","state_name")]
# super_chromhmm.ii_nn <- super_chromhmm.ii2[c("rna_fold_nwt_dko","rna_fold_fwt_dko" ,"state_name")]

# xx <- super_chromhmm.ii_nn %>% pivot_longer(cols="rna_fold_nwt_dko":"rna_fold_fwt_dko")
xx <- super_chromhmm.ii_nn %>% pivot_longer(cols="rna_fold_fwt_dko")
# xx <- super_chromhmm.ii_nn %>% pivot_longer(cols="rna_fold_nwt_dko")
# xx <- super_chromhmm.ii_nn %>% pivot_longer(cols="rna_fold_nf_wt":"rna_fold_nf_dko")

xx$state_name <- factor(xx$state_name, levels =(c(state_array)))
#xx$name <- factor(xx$name, levels = (c("rna_fold_nf_wt","rna_fold_nf_dko")))
xxx <- drop_na(xx)

#xx$name <- factor(xx$name, levels = (c("rna_fold_nf_wt","rna_fold_nf_dko")))

barplot(table(super_chromhmm.ii_nn$state_name))
# for(1:length(state_array)):
  

ggplot(data=xxx, aes(x=state_name, y=value, fill=name)) +
  geom_hline(yintercept = 0, colour = "black") +
  #geom_violin(outlier.shape = NA, fill = "magenta3") +
  geom_boxplot(outlier.shape = NA, fill = "magenta3") +
  geom_text(aes(label=..count..), y=1.8, stat='count', colour="black", size=4) +
  theme_classic() +
  scale_fill_npg() +
  ylim(c(-2,2)) +
  theme(axis.text.x = element_text(angle = 90))
  

##plot emissions table
emissions <- read.delim('~/Desktop/data/chromhmmmodel2/emissions_16.txt', sep='\t', header = TRUE)

emissions_long <- pivot_longer(emissions, cols = c("form_k27a_DKO",
                               "form_k27a_WT",
                               "form_k4m1_WT",
                               "form_k4m1_DKO"))
colnames(emissions_long) <- c("state", "name", "value")
emissions_long$name <- factor(emissions_long$name, levels = c("form_k4m1_WT","form_k4m1_DKO", "form_k27a_DKO","form_k27a_WT"))
emissions_long$state <- factor(emissions_long$state, levels = rev(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)))

ggplot() +
  geom_tile(data=emissions_long, aes(x=name, y=state, fill=value)) +
  scale_fill_gradientn(limits  = c(0,1), colours = c("white", "magenta3")) +
  theme_classic()
  









