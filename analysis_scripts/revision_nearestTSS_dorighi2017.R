#Script used for nearest neighbor analysis of MLL3/4 chip-seq peaks from Dorighi et al 2017

xx <- read.table("~/Desktop/data/Mll34_dKO_dorighi2017vsMll34_WT_dorighi2017.deseq2.FDR0.01.results.bed")
xx <- read.table("~/Desktop/data/Mll34_dKO_dorighi2017vsMll34_WT_dorighi2017.deseq2.FDR0.05.results.bed")

xx2 <- xx %>% filter(V5 < -1)
xx3 <- xx2[1:3]
colnames(xx3) <- c("chr", "start", "end")
xx3$chr <- paste0("chr", xx3$chr)
xx3$peakid <- c(seq(1,nrow(xx3)))
nrow(xx3)
#write.table(xx3, "~/Desktop/data/naiveMLL34peaks_dkocorrect.homer.bed", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
yy <- read.delim("~/Desktop/data/naiveMLL34peaks_dkocorrect.anno.bed")
yy <- yy[2:28]
yyy <- yy[,c("NWT", "FWT", "NDKO", "FDKO")]
yyyy <- gather(yyy)
yyyy$key <- factor(yyyy$key, levels = c("NWT", "FWT", "NDKO", "FDKO"))

countx <- count(yy)
countx <- paste("peaks", "=", countx, sep=' ')
grobx <- grobTree(textGrob(countx, x=0.05,  y=0.93, just="left",
                           gp=gpar(col="black", fontsize=10)))

#count unique genes used
county <- length(unique(yy$Gene.Name))
county <- paste("genes", "=", county, sep=' ')
groby <- grobTree(textGrob(county, x=0.05,  y=0.85, just="left",
                           gp=gpar(col="black", fontsize=10)))


ggplot() +
  #geom_jitter(data=matrix_oi_gather2, aes(x = key, y=value), width=0.25, size = 0.1, color = "grey30", alpha = 0.5) +
  stat_boxplot(data=yyyy, aes(x = key, y=value), geom = "errorbar", width = 0.2) +
  geom_boxplot(data=yyyy, aes(x = key, y=value, fill = key), width=0.9, outlier.shape=NA) +
  scale_fill_manual(values=c("deepskyblue3", "magenta3","deepskyblue3", "magenta3" )) +
  ylim(-3,17) +
  ylab("") +
  xlab("") +
  annotation_custom(grobx) +
  annotation_custom(groby) +
  theme_classic() +
  theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")

statstest  <- compare_means(value ~ key, data = yyyy, method="wilcox.test", p.adjust.method="BH" )


