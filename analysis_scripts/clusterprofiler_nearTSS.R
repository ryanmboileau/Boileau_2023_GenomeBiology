#Code used to perform gene ontology analysis on differential gene lists using clusterProfiler

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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] org.Mm.eg.db_3.14.0   forcats_0.5.1         stringr_1.4.0         dplyr_1.0.9           purrr_0.3.4           readr_2.1.2          
# [7] tidyr_1.2.0           tibble_3.1.7          tidyverse_1.3.1       patchwork_1.1.1       ggplot2_3.3.6         AnnotationDbi_1.56.2 
# [13] IRanges_2.28.0        S4Vectors_0.32.4      Biobase_2.54.0        BiocGenerics_0.40.0   clusterProfiler_4.2.2
# 
# loaded via a namespace (and not attached):
#   [1] fgsea_1.20.0           colorspace_2.0-3       ggtree_3.2.1           ellipsis_0.3.2         qvalue_2.26.0          XVector_0.34.0        
# [7] fs_1.5.2               aplot_0.1.6            rstudioapi_0.13        farver_2.1.0           graphlayouts_0.8.0     ggrepel_0.9.1         
# [13] bit64_4.0.5            fansi_1.0.3            scatterpie_0.1.7       lubridate_1.8.0        xml2_1.3.3             splines_4.1.2         
# [19] cachem_1.0.6           GOSemSim_2.20.0        polyclip_1.10-0        jsonlite_1.8.0         broom_1.0.0            GO.db_3.14.0          
# [25] dbplyr_2.2.1           png_0.1-7              ggforce_0.3.3          BiocManager_1.30.18    compiler_4.1.2         httr_1.4.3            
# [31] backports_1.4.1        assertthat_0.2.1       Matrix_1.4-1           fastmap_1.1.0          lazyeval_0.2.2         cli_3.3.0             
# [37] tweenr_1.0.2           tools_4.1.2            igraph_1.3.2           gtable_0.3.0           glue_1.6.2             GenomeInfoDbData_1.2.7
# [43] reshape2_1.4.4         DO.db_2.9              fastmatch_1.1-3        Rcpp_1.0.8.3           enrichplot_1.14.2      cellranger_1.1.0      
# [49] vctrs_0.4.1            Biostrings_2.62.0      ape_5.6-2              nlme_3.1-158           ggraph_2.0.5           rvest_1.0.2           
# [55] lifecycle_1.0.1        DOSE_3.20.1            zlibbioc_1.40.0        MASS_7.3-57            scales_1.2.0           tidygraph_1.2.1       
# [61] hms_1.1.1              parallel_4.1.2         RColorBrewer_1.1-3     memoise_2.0.1          gridExtra_2.3          downloader_0.4        
# [67] ggfun_0.0.6            yulab.utils_0.0.5      stringi_1.7.6          RSQLite_2.2.14         tidytree_0.3.9         BiocParallel_1.28.3   
# [73] GenomeInfoDb_1.30.1    rlang_1.0.3            pkgconfig_2.0.3        bitops_1.0-7           lattice_0.20-45        treeio_1.18.1         
# [79] shadowtext_0.1.2       bit_4.0.4              tidyselect_1.1.2       plyr_1.8.7             magrittr_2.0.3         R6_2.5.1              
# [85] generics_0.1.3         DBI_1.1.3              pillar_1.7.0           haven_2.5.0            withr_2.5.0            KEGGREST_1.34.0       
# [91] RCurl_1.98-1.7         modelr_0.1.8           crayon_1.5.1           utf8_1.2.2             tzdb_0.3.0             viridis_0.6.2         
# [97] readxl_1.4.0           grid_4.1.2             data.table_1.14.2      blob_1.2.3             reprex_2.0.1           digest_0.6.29         
# [103] gridGraphics_0.5-1     munsell_0.5.0          viridisLite_0.4.0      ggplotify_0.1.0 


#####


# BiocManager::install("clusterProfiler", force =TRUE)
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("svglite")
library(clusterProfiler)
library(AnnotationDbi)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(org.Mm.eg.db)
#library(svglite)

##specify and load desired species annotations
#species_database <- "org.Hs.eg.db"
species_database <- "org.Mm.eg.db"
library(species_database, character.only = TRUE)

##load file containing genes of interest
raw_df1 <- read.csv("~/Desktop/Data/Diffgenes/DiffgenesNWTvsFWTv3-07.07.20.anno.csv")
raw_df1$comparison <- "nwt_fwt"
raw_df2 <- read.csv("~/Desktop/Data/Diffgenes/DiffgenesFWTvsFDKOv3-07.07.20.anno.csv")
raw_df2$comparison <- "fwt_fdko"
raw_df3 <- read.csv("~/Desktop/Data/Diffgenes/DiffgenesNWTvsNDKOv3-07.07.20.anno.csv")
raw_df3$comparison <- "nwt_ndko"
raw_df4 <- read.csv("~/Desktop/Data/Diffgenes/DiffgenesNWTvsNdCDv3-07.07.20.anno.csv")
raw_df4$comparison <- "nwt_ndcd"
raw_df5 <- read.csv("~/Desktop/Data/Diffgenes/DiffgenesFWTvsFdCDv3-07.07.20.anno.csv")
raw_df5$comparison <- "fwt_fdcd"
raw_df6 <- read.csv("~/Desktop/Data/Diffgenes/DiffgenesNWTvsNCKOv3-07.07.20.anno.csv")
raw_df6$comparison <- "nwt_ncko"
raw_df7 <- read.csv("~/Desktop/Data/Diffgenes/DiffgenesFWTvsFCKOv3-07.07.20.anno.csv")
raw_df7$comparison <- "fwt_fcko"

input_deseq_array <- rbind(raw_df1, raw_df2, raw_df3, raw_df4, raw_df5, raw_df6, raw_df7)
#comparison_group <- c("nwt_fwt", "fwt_fdko", "nwt_ndko", "nwt_ndcd", "fwt_fdcd","nwt_ncko","fwt_fcko")
#comparison_group <- c("nwt_fwt", "nwt_ndko", "fwt_fdko")
comparison_group <- c("fwt_fdko")

form_dkogenesdown <- filteredlcpm %>% filter((ensembllistanno %in% k5$ensembllistanno) & FDKO_WT < -1)
form_dkoloss <- form_dkogenesdown %>% filter((ensembllistanno %in% form_wtgenes$ensembllistanno) & NDKO_WT > -1)

x <- raw_df2 %>% filter((ensembllistanno %in% form_dkoloss$ensembllistanno))
x$comparison <- "dkoloss"
input_deseq_array <- x
comparison_group <- c("dkoloss")

x <- raw_df2 %>% filter((ensembllistanno %in% form_dkogenesdown$ensembllistanno))
x$comparison <- "dkoloss"
input_deseq_array <- x
comparison_group <- c("dkoloss")

input_deseq_array <- rbind(raw_df1)
#comparison_group <- c("nwt_fwt", "fwt_fdko", "nwt_ndko", "nwt_ndcd", "fwt_fdcd","nwt_ncko","fwt_fcko")
#comparison_group <- c("nwt_fwt", "nwt_ndko", "fwt_fdko")
comparison_group <- c("nwt_fwt")


setwd("~")
plotarray_loss = list()
plotarray_gain = list()
counter=0
for(i in 1:length(comparison_group)) {
  counter = counter + 1
  input_df <- comparison_group[[i]]
  
  deseq_loss <- input_deseq_array %>% filter(comparison == comparison_group[[i]] & log2FoldChange < -1)
  deseq_gain <- input_deseq_array %>% filter(comparison == comparison_group[[i]] & log2FoldChange > 1)
  
  gene_vector <- as.character(deseq_loss$ensembllistanno)
  gene_vector2 <- bitr(gene_vector, fromType = "SYMBOL",
                       toType = c("ENSEMBL", "ENTREZID"),
                       OrgDb = species_database)
  gene_vector3 <- as.character(gene_vector2$ENTREZID)
  
  ego_bp <- enrichGO(gene          = gene_vector3,
                     OrgDb         = species_database,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     #pvalueCutoff  = 0.01,
                     qvalueCutoff  = 1,
                     readable      = TRUE)
  
  p1 <- dotplot(ego_bp, showCategory = 15, color = "qvalue") + ggplot2::theme_bw() + 
    scale_y_discrete(labels = scales::wrap_format(40)) + 
    ggplot2::theme(aspect.ratio =1.5, 
                   text = element_text(size = 18),
                   panel.border = element_rect(fill=NA, colour = "black", size=3))
  plotarray_loss[[counter]] <- print(p1)
  #ggsave(paste0(comparison_group[[i]], "_loss.go_analysis.svg"))
  
  gene_vector <- as.character(deseq_gain$ensembllistanno)
  gene_vector2 <- bitr(gene_vector, fromType = "SYMBOL",
                       toType = c("ENSEMBL", "ENTREZID"),
                       OrgDb = species_database)
  gene_vector3 <- as.character(gene_vector2$ENTREZID)
  
  ego_bp <- enrichGO(gene          = gene_vector3,
                     OrgDb         = species_database,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     #pvalueCutoff  = 0.01,
                     qvalueCutoff  = 1,
                     readable      = TRUE)
  
  p2 <- dotplot(ego_bp, showCategory = 15, color = "qvalue") + 
    ggplot2::theme_bw() +
    scale_y_discrete(labels = scales::wrap_format(40)) + 
    ggplot2::theme(aspect.ratio =1.5, 
                   text = element_text(size = 18),
                   panel.border = element_rect(fill=NA, colour = "black", size=3))#+ ggplot2::theme(aspect.ratio =4)
  
  plotarray_gain[[counter]] <- print(p2)
}

#sad but I couldn't find a loopable solution to file saving for svg files. Doesn't help that svglite and svg() appear to be R ver sensitive
#so I definitely called each object and manually saved svg using 800x1200 pixel res
plotarray_loss[[1]]
plotarray_gain[[1]]

plotarray_loss[[2]]
plotarray_gain[[2]]

plotarray_loss[[3]]
plotarray_gain[[3]]





######old code######
#raw_df <- read.csv("~/Desktop/ncc_data/D15cranial_vs_D15vagal_sig.up.csv")
#####



