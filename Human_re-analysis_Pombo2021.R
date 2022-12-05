library(tidyverse)
library(plyr)
library(ggplot2)
library(ggpubr)
library(viridis)

#read data
data_ND<-read.csv("annot_Human_ND_GBM_Full.csv")
data_R<-read.csv("annot_Human_R_GBM_Full.csv")
metadata<-read.csv("metadata.csv")
data<-rbind(data_ND, data_R)

#assign newly diagnosed and recurrent
data$type<-"newly diagnosed"
data[grep("R",data$sample), "type"]<-"recurrent"

#add sex variable
data$Sex<-metadata[match(data$sample, metadata$Sample_id), "Sex"]
data[data$cluster=="TAM 1", "cluster"]<-"Macrophages"

#Stacked bar frequency per patient
data2<-data[data$cluster %in% c("Monocytes", "Macrophages"),]
data2$cluster<-factor(data2$cluster, levels=c("Monocytes", "Macrophages"))
freq<-dplyr::count(data2, cluster, Sex, type, sample)
freq<-group_by(freq, Sex, sample, type) %>% mutate(per = n/sum(n))
aggregate(freq$per, by=list(Category=freq$sample), FUN=sum)

ggplot(freq, aes(x = sample, y = per, fill=cluster)) +
  geom_bar(position="fill",stat="identity", color="black")+
  geom_text(aes(label=n), position = position_stack(vjust = 0.5))+
  facet_wrap(Sex~., scales="free_x", nrow = 1)+
  scale_fill_manual(values=c("#eb0258", "royalblue", "gold2"))+
  ylab("% of population")+
  xlab("")+
  theme_light()+
  theme(panel.grid.major = element_blank(), legend.position = "top")


