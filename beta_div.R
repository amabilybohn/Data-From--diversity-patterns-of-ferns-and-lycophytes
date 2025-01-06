####Beta diversity AF and within each ecoregion - this procedure was repeated for ferns and lycophytes and for each ecoregion####

#creating empty dataframes for each component
sim<-data.frame()
sne<-data.frame()
sor<-data.frame()
local<-data.frame()

library(betapart)
setwd("your_directory")

pams<-list.files() #retrieving all pams

#looping to open all pams, run the multiple-site dissimilarities analysis and storing the values for the sim, sne, and sor values
for(i in 1:length(pams)){
composition<-read.csv(pams[[i]]) #opening the pam for each ecoregion and for the AF
composition<-composition[, -c(1, ncol(composition)-1, ncol(composition))] #removing the site, longitude and latitude columns

fern.AP <- beta.multi(composition)

sim<-rbind(sim, fern.AP$beta.SIM)
sne<-rbind(sne, fern.AP$beta.SNE)
sor<-rbind(sor, fern.AP$beta.SOR)
local<-rbind(local, pams[[i]])
}

#adding column names
colnames(sim)<-paste(c("beta.SIM"))
colnames(sne)<-paste(c("beta.SNE"))
colnames(sor)<-paste(c("beta.SOR"))

#joining all data in one table
all<-cbind(sim, sne, sor, local)

write.csv(all, "beta_div_ferns.csv")

#plotting the results
library(tidyverse)

setwd("your_directory")

df_fern<-read.csv("beta_div_ferns.csv", row.names = 1) #opening the results

colnames(df_fern)<-c("Turnover", "Nestedness", "Beta diversity", "Ecoregion") #arranging the column names

df_fern$Ecoregion<-c("Atlantic Dry", "Alto Paraná", "Araucaria", "Bahia Coastal", "Bahia Interior", "Mangroves", "Pernambuco Coastal", "Pernambuco Interior", "Restinga", "Serra do Mar", "Atlantic Forest") #adding the names of each ecoregion

df_fern<-pivot_longer(df_fern, cols = 1:3) #transforming into a long table
df_fern$lineage<-paste("Fern") #adding a column with the lineage name

setwd("your_directory")

df_lyco<-read.csv("beta_div_lyco.csv", row.names = 1) #opening the results

colnames(df_lyco)<-c("Turnover", "Nestedness", "Beta diversity", "Ecoregion") #arranging the column names

df_lyco$Ecoregion<-c("Atlantic Dry", "Alto Paraná", "Araucaria", "Bahia Coastal", "Bahia Interior", "Mangroves", "Pernambuco Coastal", "Pernambuco Interior", "Restinga", "Serra do Mar", "Atlantic Forest") #adding the names of each ecoregion

df_lyco<-pivot_longer(df_lyco, cols = 1:3) #transforming into a long table
df_lyco$lineage<-paste("Lycophyte")#adding a column with the lineage name

all<-rbind(df_fern, df_lyco) #joining the data

filtered <-
  all %>%  
  dplyr::filter( name=="Turnover"|name=="Nestedness") #selecting only lines of turnover and nestedness info

filtered_af<-
  filtered %>% dplyr::filter(Ecoregion=="Atlantic Forest") #filtering only the atlantic forest to plot the figure separately

library(ggplot2)
library(ggsci)
library(forcats)

#plot
ggplot(filtered_af, aes(x = lineage, y =value , fill = name)) +
  geom_bar(stat = "identity", position = "stack", width=0.5) +
  labs(x = "Lineage", y = NULL)+
  theme(panel.spacing = unit(0.1, "lines"))  +
  # facet_wrap(~lineage)+
  ylab("% of Beta diversity (βsØr)")+
  theme_bw()+
  scale_fill_d3()+
  scale_y_continuous(expand = c(0, 0), limits=c(0,1))+
  theme(text=element_text(family = "sans"),axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), legend.position = "bottom", legend.title = element_blank())+
  coord_flip()

setwd("your_directory")
ggsave("betadiv_af.tiff", plot = last_plot(), device = "tiff", width = 210, height = 150, units = "mm", dpi = 600)

#ploting for each ecoregion
filtered <-
  filtered %>%  
  dplyr::filter(!Ecoregion=="Atlantic Forest") #selecting all ecoregions but the AF

ggplot(filtered, aes(x = factor(Ecoregion, levels=rev(sort(c("Atlantic Dry", "Alto Paraná", "Araucaria", "Bahia Coastal", "Bahia Interior", "Mangroves", "Pernambuco Coastal", "Pernambuco Interior", "Restinga", "Serra do Mar")))), y = value, fill = name)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Ecoregion", y="% of Beta diversity (βsØr)")+
  theme(panel.spacing = unit(0.1, "lines"))  +
  facet_grid(lineage ~ .)+
  theme_bw()+
  scale_fill_d3()+
  scale_y_continuous(expand = c(0, 0), limits=c(0,1))+
  theme(text=element_text(family = "sans"),axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5), legend.position = "bottom", legend.title = element_blank())+
  coord_flip()

setwd("your_directory")
ggsave("betadiv_ecor.tiff", plot = last_plot(), device = "tiff", width = 210, height = 150, units = "mm", dpi = 600)

####Paired Beta-diversity. To this analysis we used the ocurrence data of each lineage intersected with the information of each ecoregion (to_pam_ecor_ferns.csv and to_pam_ecor_lyco.csv)####

setwd("your_directory")

pam<-read.csv("to_pam_ecor_ferns.csv")

library(stringr)

#changing the original names of the shapefile from the ones used here
pam$ECO_NAME<-pam$ECO_NAME %>% 
  str_replace_all(c("Alto Paraná Atlantic forests"="Alto_Parana","Araucaria moist forests"="Araucaria", "Araucaria moist forests"="Araucaria","Atlantic Coast restingas"="Restingas", "Atlantic dry forests"="Atlantic_dry", "Bahia coastal forests"="Bahia_Coastal","Bahia interior forests"="Bahia_Interior", "Pernambuco coastal forests"="Pernambuco_Coastal", "Pernambuco interior forests"="Pernambuco_Interior", "Serra do Mar coastal forests"="Serra_do_Mar","Southern Atlantic mangroves"="Mangroves")) 

pam_fern <- reshape2::dcast(pam, ECO_NAME~scientificName, length) #creating a "presence-absence matrix" with the ecoregions as lines and species as columns

#removing the first column and adding it as colnames
rownames(pam_fern)<-pam_fern[,1] 
pam_fern<-pam_fern[,2:770]

#transforming all values higher than one (repetitions) into 1 (presence)
pam_fern[pam_fern >= 2] <- 1

library(betapart)

beta_pair<- beta.pair(pam_fern) #running the paired beta diversity among ecoregions

#retrieving each component
sim<-as.matrix(beta_pair$beta.sim) 
sne<-as.matrix(beta_pair$beta.sne)
sor<-as.matrix(beta_pair$beta.sor)

#saving
write.csv(sim, "pairbeta_sim_ferns.csv")
write.csv(sne, "pairbeta_sne_ferns.csv")
write.csv(sor, "pairbeta_sor_ferns.csv")


#To plot, we manually checked the values of the combination for each ecoregion for the sim, sne and sor values and added to a unique dataframe for each lineage (paired_betadiv_all_fern.csv and paired_betadiv_all_lyco.csv)

setwd("your_directory")

fern_pair<-read.csv("paired_betadiv_all_fern.csv")
lyco_pair<-read.csv("paired_betadiv_all_lyco.csv")

library(tidyverse)
library(ggsci)

df_fern<-pivot_longer(fern_pair, cols = 2:4) #transforming into long
df_lyco<-pivot_longer(lyco_pair, cols = 2:4) #transforming into long

filtered_fern <-
  df_fern %>%  
  dplyr::filter( name=="Turnover"|name=="Nestedness") #keeping only turnover and nestedness values

filtered_lyco <-
  df_lyco %>%  
  dplyr::filter( name=="Turnover"|name=="Nestedness") #keeping only turnover and nestedness values

filtered_fern<-dplyr::rename(filtered_fern, "Component"="name") #adding the name as the component
filtered_lyco<-dplyr::rename(filtered_lyco, "Component"="name") #adding the name as the component

fern<-ggplot(filtered_fern, aes(x = reorder(Paired.Ecoregions,value), y = value, fill = Component)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Ecoregions", y = NULL)+
  theme(panel.spacing = unit(0.1, "lines"))  +
  labs(title="Ferns")+
  # facet_grid(lineage ~ .)+
  # ylab("% of Beta diversity (βsØr)")+
  theme_bw()+
  scale_fill_d3()+
  scale_y_continuous(expand = c(0, 0), limits=c(0,1), labels = c(0, 0.25, 0.5, 0.75, 1))+
  theme(text=element_text(family = "sans"), legend.position=c(.885,.07), legend.title = element_blank(),legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), plot.title = element_text(hjust = 0.5))+
  coord_flip()

lyco<-ggplot(filtered_fern, aes(x = reorder(Paired.Ecoregions,value), y = value, fill = Component)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Ecoregions", y = NULL)+
  theme(panel.spacing = unit(0.1, "lines"))  +
  labs(title="Ferns")+
  # facet_grid(lineage ~ .)+
  # ylab("% of Beta diversity (βsØr)")+
  theme_bw()+
  scale_fill_d3()+
  scale_y_continuous(expand = c(0, 0), limits=c(0,1), labels = c(0, 0.25, 0.5, 0.75, 1))+
  theme(text=element_text(family = "sans"), legend.position=c(.885,.07), legend.title = element_blank(),legend.background = element_blank(), legend.box.background = element_rect(colour = "black"), plot.title = element_text(hjust = 0.5))+
  coord_flip()

a<-ggarrange(fern, lyco, ncol=2, bottom="% of Beta diversity (βsØr)", widths = c(1,1))

ggsave("betadiv_pair_eco.tiff", plot = a, device = "tiff", width = 300, height = 150, units = "mm", dpi = 600)
