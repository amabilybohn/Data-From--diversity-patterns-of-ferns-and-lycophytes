###Additiva partitioning of diversity - adipart. This procedure was repeated for ferns and lycophytes separately, with the pam created at the gdm.R####

setwd("your_directory")
pam<-read.csv("pam_fern.csv", row.names = 1)

setwd("C:/Users/amabi/OneDrive/Área de Trabalho/OneDrive/Doutorado/Projeto/Cap1/ecoregions")
eco<-read.csv("eco_ferns.csv", header=T, sep=",") #this was an intersection of the grids and the shapefile with the ecoregions (eco_ferns.csv and eco_lyco.csv)

library(metacom)

test<-adipart(y=pam[,2:780], x=eco, nsimul=999) #analyzing the data

test #checking the results

###separating the results for each component

a_obs<-paste(c(26, "Observed", "Alpha"))
a_esp<-paste(c(23, "Expected", "Alpha"))
b1_obs<-paste(c(363, "Observed", "Beta1"))
b1_esp<-paste(c(471, "Expected", "Beta1"))
b2_obs<-paste(c(390, "Observed", "Beta2"))
b2_esp<-paste(c(293, "Expected", "Beta2"))

#binding the results
all<-as.data.frame(rbind(a_obs, a_esp, b1_obs, b1_esp, b2_obs, b2_esp))
colnames(all)<-c("Value", "obs_esp", "level")
all$Value<-as.numeric(all$Value)


#plotting
library(ggsci)
library(ggplot2)
library(dplyr)

ggplot(all,aes(x=obs_esp, y=Value, fill=level)) +
  geom_bar(stat="identity", width=0.3)+
  theme_minimal()+
  scale_fill_brewer(palette=8, direction = 1)+
  scale_y_reverse()


