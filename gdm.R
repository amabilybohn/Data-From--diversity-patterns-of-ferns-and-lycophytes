####Generalized Dissimilarity Modelling###

####creating presence-absence matrices for all occurrences from the AF####
setwd("your_directory") #add your working directory

#load the occurrences of ferns (topam_ferns.csv) or lycophytes (topam_lyco.csv)
topam<-read.csv("topam_ferns.csv")

#renaming the columns 
topam<-dplyr::rename(topam, site=id, Long=Coordenadas.médias_MEAN_X, Lat=Coordenadas.médias_MEAN_Y, species=scientificName)

#removing empty spaces and .shp from species names, also transforming into factor
topam$species<-gsub(" ", "_", topam$species)
topam$species<-gsub(".shp", "", topam$species)
topam$species<-as.factor(topam$species)

#retrieving only one occurrence of sites and coordinates
coords<-unique(topam[,c(4:6)])

#transforming into a presence-absence matrix
pam_fern <- reshape2::dcast(topam, site~species, length)

#removing the species name columns
pam2<-pam_fern[,2:length(pam_fern)]

#changing numbers higher than one to one
pam2[pam2 >= 2] <- 1

#adding the site column to pam2
pam2$site<-pam_fern$site

#merging the coordinate data to the final pam
pam_fern<-merge(pam2, coords, by="site")

#saving the pam
write.csv(pam_fern, "pam_fern.csv")

####creating presence-absence matrices for each AF ecoregion ####

setwd("your_directory")

#selecting all .csv files for the 10 ecoregions
a<-list.files(pattern=".csv")

#loop to create the pams as we did above#
for(i in 1:length(a)){
setwd("your_directory")

topam<-read.csv(a[[i]])

topam<-dplyr::rename(topam, site=id, Long=Coordenadas.médias_MEAN_X, Lat=Coordenadas.médias_MEAN_Y, species=scientificName)
topam<-dplyr::select(topam, site, species, Long, Lat)
topam$species<-gsub(" ", "_", topam$species)
topam$species<-gsub(".shp", "", topam$species)
topam$species<-as.factor(topam$species)

coords<-unique(topam[,c(1,3:4)])

pam <- reshape2::dcast(topam, site~species, length)
pam2<-pam[,2:length(pam)]

pam2[pam2 >= 2] <- 1

pam2$site<-pam$site
              
pam2<-merge(pam2, coords, by="site")

setwd("your_directory")
write.csv(pam2, paste("pam_", a[[i]], sep=""), row.names = F)
}

####GDM for all AF - the same step by step was done for ferns and lycophytes and for each ecoregion alone####

setwd("your_directory")

fern<-read.csv("pam_fern.csv", row.names = 1) #opening the pam

setwd("C:/Users/amabi/OneDrive - ufpr.br/Doutorado Amabily/Projeto/Cap1/gdm/env/")
a<-list.files(pattern=".tif") #creating a list with all environmental variables
env<-raster::stack(a) #loading as a stack file

var1<-terra::as.data.frame(terra::rast(env[[1]])) #selecting the first variable and transforming it into a data.frame
var1$grids<-rownames(var1) #creating a new column with the row names of the data frame

#doing the same with all other variables and then merging them together in the same data frame (var1) by the grid number
for(i in 2:length(env)){
  var<-terra::as.data.frame(terra::rast(env[[i]]))
  var$grids<-rownames(var)
  var1<-merge(var1, var, by="grids")
}

#calculating the correlation among variables
library(usdm)
cor<-vifcor(var1[,2:37], th=0.7)
cor

#running the gdm with only uncorrelated variables

library(gdm)

setwd("your_directory")

a<-list.files(pattern=".tif") #creating a list with all environmental variables
a<-a[(c(1,4,7,12,17,18,19,24,25,27,28,31))] #selecting only those with no correlation
env<-raster::stack(a) #creating a stack with them

#this function is from the gdm package and is used to format our pam and joining the environmental  data. Here we also weighted the sites by richness
exFormat1a <- formatsitepair(fern, 1, XColumn="Long",
                             YColumn="Lat", predData=env, siteColumn="site", weightType = "richness")

gdmTab<-na.omit(exFormat1a) #removing NAs
gdm.1 <- gdm(data=gdmTab, geo=TRUE) #running the gdm
summary(gdm.1) #looking to the results

#calculating the porcentages of explanation for each variable
(1.253*100)/gdm.1$explained #6.122186
(0.972*100)/gdm.1$explained #4.749214
(0.874*100)/gdm.1$explained #4.270384
(0.672*100)/gdm.1$explained #3.283407
(0.636*100)/gdm.1$explained #3.10751
(0.623*100)/gdm.1$explained #3.043992
(0.615*100)/gdm.1$explained #3.004904
(0.615*100)/gdm.1$explained #3.004904
(0.588*100)/gdm.1$explained #2.872981
(0.57*100)/gdm.1$explained #2.785033
(0.504*100)/gdm.1$explained #2.462555

#creating the PCA axis for plots
transRasts <- gdm.transform(model=gdm.1, data=env) #transform the environmental data using the gdm results

rastDat <- na.omit(raster::getValues(transRasts)) #get the data as a table

pcaSamp <- prcomp(rastDat) #perform the principle components analysis

# Predict the first three principle components for every cell in the rasters
pcaRast <- raster::predict(transRasts, pcaSamp, index=1:3)

setwd("your_directory")

terra::writeRaster(pcaRast, "resu_fern_allMA.tif",overwrite=TRUE)

#scale the PCA rasters to make full use of the color spectrum
pcaRast[[1]] <- (pcaRast[[1]]-pcaRast[[1]]@data@min) /
  (pcaRast[[1]]@data@max-pcaRast[[1]]@data@min)*255
pcaRast[[2]] <- (pcaRast[[2]]-pcaRast[[2]]@data@min) /
  (pcaRast[[2]]@data@max-pcaRast[[2]]@data@min)*255
pcaRast[[3]] <- (pcaRast[[3]]-pcaRast[[3]]@data@min) /
  (pcaRast[[3]]@data@max-pcaRast[[3]]@data@min)*255

#Plot the three PCA rasters simultaneously, each representing a different color 
#(red, green, blue)
raster::plotRGB(pcaRast, r=1, g=2, b=3)


