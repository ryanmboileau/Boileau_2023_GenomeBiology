#Generation of Fig1 supplemental plots including qPCR, histone western, and growth assays. Includes code to perform their statistic analysis

#import qPCRdata
qpcrdata <- read.csv("~/Desktop/data/MLLqPCRsummaryv2.csv", header = TRUE, row.names = 1)
qpcrdata <- as.data.frame(t(qpcrdata))

qpcrdata$replicate <- c("rep1","rep2","rep3")
qpcrdata$genotype <- c(rep("WT", 6), rep("MLL3KO", 6), rep("DKO", 6))

qpcrdata$cellstate <- c(rep("naive", 3),rep("form", 3) )
qpcrdata$geno_cellstate <- paste0(qpcrdata$cellstate, "_", qpcrdata$genotype)

rnalist <- list("Klf2", "Klf4", "Rex1", "Fgf5", "Otx2", "Dnmt3b")
statsnames <- c("Klf4", "Klf4", "Rex1", "Fgf5", "Otx2", "Dnmt3b")

plotlist=list()
statslist=list()

counter = 0
plotarray <- for(i in 1:length(rnalist)) {
  counter = counter + 1
  
  matrix_oi_gather <- pivot_longer(qpcrdata, cols = c(rnalist[[i]]))
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  matrix_oi_gather2$genotype <- factor(matrix_oi_gather2$genotype, levels = c("WT", "MLL3KO", "DKO"))
  matrix_oi_gather2$cellstate <- factor(matrix_oi_gather2$cellstate, levels = c("naive", "form"))
  matrix_oi_gather2$geno_cellstate <- factor(matrix_oi_gather2$geno_cellstate, levels = c("naive_WT", "form_WT", "naive_MLL3KO", "form_MLL3KO", "naive_DKO", "form_DKO"))
  
  #statslist[[counter]] <- compare_means(value~cellstate + genotype, data = matrix_oi_gather2, method="anova",paired=TRUE)
  two.way <- aov(value ~ genotype * cellstate, data=matrix_oi_gather2) #conduct two way anova
  statslist1 <- TukeyHSD(two.way) #post hoc analysis of levels by tukey's test
  statslist[[counter]] <- as.data.frame(statslist1[[3]], optional = TRUE)
  
  plotlist[[counter]] <- print(ggplot(data=matrix_oi_gather2, aes(x = geno_cellstate, y=value, color = cellstate)) +
                                 stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", color="black", width=1) +
                                 stat_summary(fun=mean, geom="point", color="black") +
                                 geom_point() +
                                 ylab("") +
                                 xlab("") + geom_hline(yintercept=0) +
                                 scale_color_manual(values=c("deepskyblue3", "magenta3")) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
      )
}
wrap_plots(plotlist)

#statslist[[1]]

counter = 0
stats_grouped <- c()
statsarray <- for(i in 1:length(statslist)) {
  counter = counter + 1
  stats_oi <- statslist[[i]]
  stats_oi$group <- statsnames[[counter]]

  stats_grouped <- rbind(stats_grouped, stats_oi)

}


##import histone western data, norm to naive wt
histdata <- read.csv("~/Desktop/data/histonewestsummaryv2.csv", header = TRUE, row.names = 1)
histdata <- as.data.frame(t(histdata))

histdata$replicate <- c("rep1","rep2")
histdata$genotype <- c(rep("WT", 4), rep("MLL3KO", 4), rep("DKO", 4))

histdata$cellstate <- c(rep("naive", 2),rep("form", 2) )
histdata$geno_cellstate <- paste0(histdata$cellstate, "_", histdata$genotype)

histlist <- list("H3K4m1", "H3K4m2", "H3K4m3", "H3K27a")
statsnames <- c("H3K4m1", "H3K4m2", "H3K4m3", "H3K27a")

plotlist=list()
statslist=list()

counter = 0
plotarray <- for(i in 1:length(histlist)) {
  counter = counter + 1
  
  matrix_oi_gather <- pivot_longer(histdata, cols = c(histlist[[i]]))
  matrix_oi_gather2 <-  matrix_oi_gather %>% filter(!is.na(value))
  matrix_oi_gather2$genotype <- factor(matrix_oi_gather2$genotype, levels = c("WT", "MLL3KO", "DKO"))
  matrix_oi_gather2$cellstate <- factor(matrix_oi_gather2$cellstate, levels = c("naive", "form"))
  matrix_oi_gather2$geno_cellstate <- factor(matrix_oi_gather2$geno_cellstate, levels = c("naive_WT", "form_WT", "naive_MLL3KO", "form_MLL3KO", "naive_DKO", "form_DKO"))
  
  #statslist[[counter]] <- compare_means(value~cellstate + genotype, data = matrix_oi_gather2, method="anova",paired=TRUE)
  two.way <- aov(value ~ genotype * cellstate, data=matrix_oi_gather2) #conduct two way anova
  statslist1 <- TukeyHSD(two.way) #post hoc analysis of levels by tukey's test
  statslist[[counter]] <- as.data.frame(statslist1[[3]], optional = TRUE)
  
  plotlist[[counter]] <- print(ggplot() +
                                 #geom_boxplot(data=matrix_oi_gather2, aes(x = geno_cellstate, y=value), color = "black") +
                                 geom_point(data=matrix_oi_gather2, aes(x = geno_cellstate, y=value, color = cellstate), shape="-", size=10, color="black", stat="summary", fun="mean") +
                                 geom_point(data=matrix_oi_gather2, aes(x = geno_cellstate, y=value, color = cellstate)) +
                                 #stat_summary(fun=mean, geom="point", shape=23, size=4) +
                                 ylab("") +
                                 xlab("") + geom_hline(yintercept=0) + geom_hline(yintercept=1) +
                                 scale_color_manual(values=c("deepskyblue3", "magenta3")) +
                                 theme_classic() +
                                 theme(aspect.ratio = 2, text = element_text(size = 20), legend.position = "none")
  )
}
wrap_plots(plotlist)
statslist[[1]]

counter = 0
stats_grouped <- c()
statsarray <- for(i in 1:length(statslist)) {
  counter = counter + 1
  stats_oi <- statslist[[i]]
  stats_oi$group <- statsnames[[counter]]
  
  stats_grouped <- rbind(stats_grouped, stats_oi)
  
}

##supp growth assay
naigrowthdata <- read.csv("~/Desktop/data/naive_growthassaysummary2.csv", header = TRUE)
ggplot() +
  geom_point(data=naigrowthdata, aes(x=time, y=abs, color=genotype)) +
  stat_summary(data=naigrowthdata, aes(x=time, y=abs, color=genotype), fun=mean, geom="line") +
  theme_classic() + scale_x_discrete(drop=TRUE)



formgrowthdata <- read.csv("~/Desktop/data/form_growthassaysummary2.csv", header = TRUE)
ggplot() +
  geom_point(data=formgrowthdata, aes(x=time, y=abs, color=genotype)) +
  stat_summary(data=formgrowthdata, aes(x=time, y=abs, color=genotype), fun=mean, geom="line") +
  theme_classic() + scale_x_discrete(drop=TRUE)





