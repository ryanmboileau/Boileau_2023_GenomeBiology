#Generates pie charts shown in fig.4A

##figure 4 filtering
#####piechartfiltering fig4a######

#losingk4m1 naive with k27a categories

#k4m1only 3528
#nonk27a 913
#losek27a 1370
#gaink27a 149

##H3K4m1 loss changing pie chart
tt1 <- c("k4m1only", "nonk27a", "losek27a", "gaink27a")
tt2 <- c(3528, 913,1370,149)
tt <- data.frame(type = tt1,value = tt2)

tt$type <- factor(tt$type, levels = c("k4m1only", "nonk27a", "losek27a", "gaink27a"))


ttt1 <- ggplot(tt, aes(x="", y=value, fill=type))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  theme_void()

#gainingk4m1 naive with k27a categories

#k4m1only 4217
#nonk27a 146
#losek27a 48
#gaink27a 1172

##H3K4m1 loss changing pie chart
tt1 <- c("k4m1only", "nonk27a", "losek27a", "gaink27a")
tt2 <- c(4217,146,48,1172)
tt <- data.frame(type = tt1,value = tt2)
tt$type <- factor(tt$type, levels = c("k4m1only", "nonk27a", "losek27a", "gaink27a"))
                  

ttt2 <- ggplot(tt, aes(x="", y=value, fill=type))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  theme_void()

#nonchangingk4m1 with k27a categories

#k4m1only 27474
#nonk27a 4353
#losek27a 3498
#gaink27a 2264

##H3K4m1 non changing pie chart
tt1 <- c("k4m1only", "nonk27a", "losek27a", "gaink27a")
tt2 <- c(27474, 4353,3498,2264)
tt <- data.frame(type = tt1,value = tt2)
tt$type <- factor(tt$type, levels = c("k4m1only", "nonk27a", "losek27a", "gaink27a"))
                  
ttt3 <- ggplot(tt, aes(x="", y=value, fill=type))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  theme_void()

ttt1 + ttt2 + ttt3

