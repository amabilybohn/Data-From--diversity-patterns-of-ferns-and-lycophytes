######################################################
###SCRIPT PCoA, PERMANOVA, BETA DIVERSITY ANALYSIS####
######################################################
#this step was repeated for ferns and lycophytes separately (pam_fern.csv and pam_lyco.csv)

setwd("your_directory")

library(vegan)
library(pairwiseAdonis)
library(ggplot2)

#PCOA----
composition<-read.csv("pam_fern.csv", header=TRUE, sep=",")
ecorregions<-read.csv("ecoregions.csv", header=T, sep=",")

colnames(eco)<-c("site", "local", "value")

#FIRST STEP:Merging the spreadsheets----
data_1<-merge(composition, eco, by="site")

#selecting the species data
composition2<-data_1[,3:781]

##SECOND STEP:running distance matrix
dist <- vegdist(composition2, method = "jaccard")

#running the PCoA with the Cailliez correction
df.pcoa<-pcoa(dist, correction="cailliez")

#calculating the percentage of explanation of each axis of the PCoA
x_label<-round(df.pcoa$values$Rel_corr_eig[1]*100,2)
y_label<-round(df.pcoa$values$Rel_corr_eig[2]*100,2)
x_label
y_label

names(data_1)
region<-data_1$local
data <- cbind(df.plot, region)

library(ggplot2)

ggplot(data,aes(x=Axis.1,y=Axis.2, color = region, bg = region, fill=region))+ 
  geom_point(alpha=0.7, size=5) + 
  theme_bw() + theme(panel.grid = element_blank()) + 
  labs(x=paste0("PCoA1 (",x_label,"%)"), y=paste0("PCoA2 (",y_label,"%)")) + 
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) + 
  stat_ellipse()+
  #scale_fill_discrete(labels = c("North", "South", "Central", "Central-North", "Central-South"))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent', size = 1.2),
        legend.title=element_blank()) + 
  theme(axis.text = element_text(colour = 'black', size = 12), 
        axis.ticks = element_line(colour = 'black', size = 0.6), axis.title = element_text(size = 12), 
        legend.key.size = unit(0.9, 'cm'), legend.text = element_text(size=12), 
        axis.ticks.x = element_line(size=0.8), axis.ticks.y = element_line(size=0.8), 
        axis.ticks.length = unit(0.18, 'cm') ) + scale_shape_manual(values = c(21,21,21,21,21,21))  

#THIRD STEP: calculating the difference between compositions and ecoregions with PERMANOVA
fern_ado<-adonis2 (dist~as.factor(data_1$local))
fern_ado

##calculating the pairwise difference between ecoregions
fern_ado_p<-pairwise.adonis(dist,factors=as.factor(data_1$local))
fern_ado_p
