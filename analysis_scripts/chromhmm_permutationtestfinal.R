#Function and code utilized to perform statistical tests on nearby transcription for chromhmm states in Fig.5B and 5D

perm <- function(exp, control) {
  x_shuff <- as.data.frame(sample(exp))
  y_shuff <- as.data.frame(sample(control, size = length(exp)))
  x_shuff$grouping <- "experimental"
  colnames(x_shuff) <- c("value", "grouping")
  y_shuff$grouping <- "background"
  colnames(y_shuff) <- c("value", "grouping")
  bb <- as.data.frame(rbind(x_shuff, y_shuff))
  bbb <- compare_means(value ~ grouping, data = bb, paired=FALSE)
  bbb$p
}

samp2 <- xxxx %>% filter(state_name == "state_10")
#test1 <- permutation.test(samp1$value2, samp2$value2, 10000)
#form state array
state_array <- c(paste0("state_0", seq(1,9)), paste0("state_", seq(10,14)), "state_16")

#naive state array
state_array <- c(paste0("state_0", seq(1,6)),paste0("state_0", seq(8,9)), paste0("state_", seq(10,15)), "state_16")

samp1 <- xxxx %>% filter(state_name == "all")
permlist <- vector()
for(i in 1:length(state_array)){
  print(state_array[[i]])
  samp1 <- xxxx %>% filter(state_name == "all")
  samp2 <- xxxx %>% filter(state_name == state_array[[i]])
  dfxx <- replicate(n=1000, expr = perm(samp2$value2, samp1$value2))
  permlist[[i]] <- (mean(dfxx))
  #print(nrow(samp2))
  }

dfy <- as.data.frame(cbind(state_array, permlist))
dfy$bhcorrect <- p.adjust(dfy$permlist, method = "BH")
print(dfy$bhcorrect)
write.table(dfy, "~/Desktop/chromhmm/hmmstats_naive_unique.tsv", quote = FALSE, row.names = F, sep = '\t')

