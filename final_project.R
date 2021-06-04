#load packages
library(readr)
library(tidyr)
library(dplyr)
options(stringsAsFactors = FALSE)
#pre-processing the data
#read in Aedes aegypti data
aedes_aegypti <- read_csv("aedes_aegypti_data.csv")
nrow(aedes_aegypti)
#filter out species data so it is presence only - ie no 0 values in specimens collected
aedes_aegypti<- aedes_aegypti%>%
  filter(Specimens_collected>0)

#load sdm package
library(sdm)
installAll() #only first time after installing sdm package

#create coordinates dataframe
aegypti_coords <- aedes_aegypti%>%
  dplyr::select(Longitudes,Latitudes)
#add a column with boolean indication of presence (all will be 1 as 0 values are omitted in previous step)
aegypti_coords$Species <- 1

#drop na values
aegypti_coords <- aegypti_coords%>%
  drop_na()
#convert species data to spatial point dataframe
coordinates(aegypti_coords) <- c('Longitudes','Latitudes')

#-----------------------------------------------------------
#load packages
library(raster)
library(dismo)
library(mapview)

#download bioclimatic and elevation data from WorldClim
#read in data from WorldClim
bio <- list.files(path = "wc2.1_2.5m_bio",pattern = '.tif', full.names = TRUE)
elevation <- list.files(pattern = "wc2.1_2.5m_elev.tif", full.names = TRUE)
#append the list of files to add elevation raster layer
variables <- c(bio, elevation)
variables <- raster::stack(variables)
#change the names of variables to bio1 bio2 bio3 etc
names(variables) <- gsub(pattern = "wc2.1_2.5m_", replacement = "", x = names(variables))
#plot elevation
plot(variables[[20]], main="Elevation", xlab="Longitude",ylab="Latitude",cex.axis=1.2, cex.lab=1.2, cex.main=1.2,col= terrain.colors(12))
#plot presence points on the map
points(aegypti_coords)
#select the needed geospatial subset of data (the study area) and fit the model based on the data spread over this area by specifying a rectangle
box <- drawExtent()
#crop aegypti data
aegypti_coords <- aegypti_coords%>%
  crop(box)
#highlight the chosen points in red
points(aegypti_coords,col = "red")
#crop the map
variables_c<- crop(variables, box)
plot(variables_c[[20]])
points(aegypti_coords,col = "red")
#-----------------------------------
#check for multicollinearity
#load package
library(usdm)
#determine VIF of each variable
vif(variables_c)
#extract environmental data at coordinates of occurrences from raster variable
extracts <- raster::extract(variables_c,aegypti_coords)
head(extracts)
#use vifstep to check bias inflection factor and if it is greater than 10 then it is a sign of multicollinearity. Alternatively can use v_cor <- vifcor(extracts)
v_step <- vifstep(extracts)
#exclude variables with collinearity problem which are used as predictor variables in the next step
variables_ex<- exclude(variables_c,v_step)




#-------------------------------------
#creating data object
#specify a formula by species and predictor variables_ex, the train data is aegypti_coords
#randomly generate pseudo-absence data as the type is presence-only to use other methods available further
data_object <- sdmData(Species ~.,aegypti_coords, predictors = variables_ex, bg = list(method= "gRandom", n = 150))
#run sdm to fit the models --> generate 15 bootstrap replicates per method
model <- sdm(Species~., data = data_object, methods = c('glm','brt','rf','maxent','gam'), replication = c('boot'),test.p = 30, n = 5, parallelSettings = list(ncore= 4, method = "parallel"))#4 model methods/test.p = test percentage/ n = 5 means that the replication of bootstapping or subsampling is repeated 5 times + parallel is to set the number of cores
#get model performance statistics
info <- sdm::getModelInfo(model)
info
model_stats <- getEvaluation(model)
info <- info%>%
  inner_join(model_stats)
info
write_csv(info, "model_info.csv")

method_info<- info%>%
  group_by(method)%>%
  summarise(AUC = mean(AUC),TSS = mean(TSS), COR = mean(COR), Deviance = mean(Deviance))%>%
  as.data.frame()
method_info
write_csv(method_info, "method_info.csv")
#generates responce curve
rc <- rcurve(model)
rc
#generate roc curves
roc<- roc(model)
#see the percentage of importance for different variables only for a specified model based on correlation metric and AUC metric
var_imp<- getVarImp(model)
var_imp
glm <- getVarImp(model,method = "glm")
glm
gam <- getVarImp(model,method = "gam")
gam
brt <- getVarImp(model,method = "brt")
brt
rf<- getVarImp(model,method = "rf")
rf
maxent <- getVarImp(model,method = "maxent")
maxent
#results for all variables and models
plot(getVarImp(model))
#results for each method separately 
plot(getVarImp(model,method = "glm"), main = "GLM Relative Variable Importance")
plot(getVarImp(model,method = "gam"), main = "GAM Relative Variable Importance")
plot(getVarImp(model,method = "brt"), main = "BRT Relative Variable Importance")
plot(getVarImp(model,method = "rf"), main = "RF Relative Variable Importance")
plot(getVarImp(model,method = "maxent"), main = "Maxent Relative Variable Importance")
#graphical user interface of the package
gui(model)

#---------------------------------------------
#generate predictions based on current climate
prediction <- predict(model, variables_c, filename = "prediction.tif", overwrite = TRUE, mean = T)
#generate color pallet
clr <- colorRampPalette(c("#3E49BB","#349808","yellow","orange","red","darkred"))
#set smaller margins so that there is less empty space
par(mar=c(3,3,1,1), mgp=c(2,1,0))
#plot predictions for each method
plot(prediction[[1]],col = clr(200),main = "GLM")
plot(prediction[[2]],col = clr(200),main = "BRT")
plot(prediction[[3]],col = clr(200),main = "RF")
plot(prediction[[4]],col = clr(200),main = "Maxent")
plot(prediction[[5]],col = clr(200),main = "GAM")
#ensemble the models
ensem <- ensemble(model, variables_c, filename = "ensemble.tif", setting = list(method = "mean-weighted",stat = "TSS",opt =2 ))
#plot ensemble model's projection
plot(ensem,col = clr(200))


#--------------------------------
#2041-2160 projections
#SSP 126
ssp126 <- list.files(path = "ssp126",pattern = '.tif', full.names = TRUE, recursive = T)
ssp126<- raster::stack(ssp126)
elevation <- raster("wc2.1_2.5m_elev.tif")
#add elevation layer to bioclimatic variables
ssp126 <- addLayer(ssp126,elevation)
names(ssp126) <- names(variables)
#crop to coordinates of Florida
ssp126_c <- crop(ssp126, box)
#use ssp126 layers as predictor variables in ensemble
ensem_ssp126<- ensemble(model, ssp126_c, filename = "ssp126.tif", setting = list(mehtod = "mean-weighted", stat = "TSS",opt=2))

#ssp 245
ssp245 <- list.files(path = "ssp245",pattern = '.tif', full.names = TRUE, recursive = T)
ssp245<- raster::stack(ssp245)
elevation <- raster("wc2.1_2.5m_elev.tif")
#add elevation layer to bioclimatic variables
ssp245 <- addLayer(ssp245,elevation)
names(ssp245) <- names(variables)
#crop by coordinates of Florida
ssp245_c <- crop(ssp245, box)
#use ssp245 layers as predictor variables in ensemble
ensem_ssp245<- ensemble(model, ssp245_c, filename = "ssp245.tif", setting = list(mehtod = "mean-weighted", stat = "TSS",opt=2))

#ssp 370
ssp370 <- list.files(path = "ssp370",pattern = '.tif', full.names = TRUE, recursive = T)
ssp370<- raster::stack(ssp370)
elevation <- raster("wc2.1_2.5m_elev.tif")
#add elevation layer to bioclimatic variables
ssp370 <- addLayer(ssp370,elevation)
#change variable names to bio1, bio2, bio3 etc
names(ssp370) <- names(variables)
#crop to coordinates of FLorida
ssp370_c <- crop(ssp370, box)
#use ssp370 layers as predictor variables in ensemble using TSS threshold
ensem_ssp370<- ensemble(model, ssp370_c, filename = "ssp370.tif", setting = list(mehtod = "mean-weighted", stat = "TSS",opt=2))

#ssp 585
ssp585 <- list.files(path = "ssp585",pattern = '.tif', full.names = TRUE, recursive = T)
ssp585<- raster::stack(ssp585)
elevation <- raster("wc2.1_2.5m_elev.tif")
#add elevation layer to bioclimatic variables
ssp585 <- addLayer(ssp585,elevation)
#change variable names to bio1, bio2, bio3 etc
names(ssp585) <- names(variables)
#crop to coordinates of FLorida
ssp585_2100_c <- crop(ssp585, box)
#use ssp585 layers as predictor variables in ensemble
ensem_ssp585<- ensemble(model, ssp585_c, filename = "ssp585.tif", setting = list(mehtod = "mean-weighted", stat = "TSS",opt=2))

#2081-2100 projections
#SSP 126
ssp126_2100 <- list.files(path = "ssp126_100",pattern = '.tif', full.names = TRUE, recursive = T)
ssp126_2100<- raster::stack(ssp126_2100)
elevation <- raster("wc2.1_2.5m_elev.tif")
#add elevation layer to bioclimatic variables
ssp126_2100 <- addLayer(ssp126_2100,elevation)
names(ssp126_2100) <- names(variables)
#crop to coordinates of Florida
ssp126_2100_c <- crop(ssp126_2100, box)
#use ssp126 layers as predictor variables in ensemble
ensem_ssp126_2100<- ensemble(model, ssp126_2100_c, filename = "ssp126_100.tif", setting = list(mehtod = "mean-weighted", stat = "TSS",opt=2))

#ssp 245
ssp245_2100 <- list.files(path = "ssp245_100",pattern = '.tif', full.names = TRUE, recursive = T)
ssp245_2100<- raster::stack(ssp245_2100)
elevation <- raster("wc2.1_2.5m_elev.tif")
#add elevation layer to bioclimatic variables
ssp245_2100 <- addLayer(ssp245_2100,elevation)
names(ssp245_2100) <- names(variables)
#crop by coordinates of Florida
ssp245_2100_c <- crop(ssp245_2100, box)
#use ssp245 layers as predictor variables in ensemble
ensem_ssp245_2100<- ensemble(model, ssp245_2100_c, filename = "ssp245_100.tif", setting = list(mehtod = "mean-weighted", stat = "TSS",opt=2))

#ssp 370
ssp370_2100 <- list.files(path = "ssp370_100",pattern = '.tif', full.names = TRUE, recursive = T)
ssp370_2100<- raster::stack(ssp370_2100)
elevation <- raster("wc2.1_2.5m_elev.tif")
#add elevation layer to bioclimatic variables
ssp370_2100 <- addLayer(ssp370_2100,elevation)
#change variable names to bio1, bio2, bio3 etc
names(ssp370_2100) <- names(variables)
#crop to coordinates of FLorida
ssp370_2100_c <- crop(ssp370_2100, box)
#use ssp370 layers as predictor variables in ensemble
ensem_ssp370_2100<- ensemble(model, ssp370_2100_c, filename = "ssp370_100.tif", setting = list(mehtod = "mean-weighted", stat = "TSS",opt=2))

#ssp 585
ssp585_2100 <- list.files(path = "ssp585",pattern = '.tif', full.names = TRUE, recursive = T)
ssp585_2100<- raster::stack(ssp585_2100)
elevation <- raster("wc2.1_2.5m_elev.tif")
#add elevation layer to bioclimatic variables
ssp585_2100 <- addLayer(ssp585_2100,elevation)
#change variable names to bio1, bio2, bio3 etc
names(ssp585_2100) <- names(variables)
#crop to coordinates of FLorida
ssp585_2100_c <- crop(ssp585_2100, box)
#use ssp585 layers as predictor variables in ensemble
ensem_ssp585_2100<- ensemble(model, ssp585_2100_c, filename = "ssp585_tss_p.tif", setting = list(mehtod = "mean-weighted", stat = "TSS",opt=2))

#---------------------------------------
#suitability change 
#set suitability threshold 
eval<- getEvaluation(model, stat = c("threshold"),opt =2)
mean_threshold <- mean(eval$threshold)
now_bin <- ensem >= mean_threshold

#2041-2060 
#ssp126
ssp126_bin  <- ensem_ssp126 >= mean_threshold
# Create a map that shows gains and losses pf suitability
change1 <- 2 * now_bin + ssp126_bin
#ssp245
ssp245_bin  <- ensem_ssp245 >= mean_threshold
# Create a map that shows gains and losses pf suitability
change2 <- 2 * now_bin + ssp245_bin
#ssp370
ssp370_bin  <- ensem_ssp370 >= mean_threshold
# Create a map that shows gains and losses pf suitability
change3 <- 2 * now_bin + ssp370_bin
#ssp585
ssp585_bin  <- ensem_ssp585 >= mean_threshold
# Create a map that shows gains and losses pf suitability
change4 <- 2 * now_bin + ssp585_bin
#2081-2100
#ssp126
ssp126_bin_100 <- ensem_ssp126_2100 >= mean_threshold
# Create a map that shows gains and losses pf suitability
change5 <- 2 * now_bin + ssp126_bin_100
#ssp245
ssp245_bin_100 <- ensem_245_2100 >= mean_threshold
# Create a map that shows gains and losses pf suitability
change6 <- 2 * now_bin + ssp245_bin_100
#ssp370
ssp370_bin_100 <- ensem_ssp370_2100 >= mean_threshold
# Create a map that shows gains and losses pf suitability
change7 <- 2 * now_bin + ssp370_bin_100
#ssp585
ssp585_bin_100 <- ensem_ssp585_2100 >= mean_threshold
# Create a map that shows gains and losses pf suitability
change8 <- 2 * now_bin + ssp585_bin_100

#create figure for 2041-2060
par(mfrow=c(2,4))
plot(ensem_ssp126)
plot(change1, col=c('white', 'green', 'red', 'blue'), 
     breaks=c(-0.5,0.5,1.5,2.5, 3.5), 
     axis.args=list(at=c(0,1,2,3), label=c('Neither', 'Suitability Gain', 'Suitability Loss', 'No change')))
plot(ensem_ssp245)
plot(change2, col=c('white', 'green', 'red', 'blue'), 
     breaks=c(-0.5,0.5,1.5,2.5, 3.5), 
     axis.args=list(at=c(0,1,2,3), label=c('Neither', 'Suitability Gain', 'Suitability Loss', 'No change')))
plot(ensem_ssp370)
plot(change3, col=c('white', 'green', 'red', 'blue'), 
     breaks=c(-0.5,0.5,1.5,2.5, 3.5), 
     axis.args=list(at=c(0,1,2,3), label=c('Neither', 'Suitability Gain', 'Suitability Loss', 'No change')))
plot(ensem_ssp585)
plot(change4, col=c('white', 'green', 'red', 'blue'), 
     breaks=c(-0.5,0.5,1.5,2.5, 3.5), 
     axis.args=list(at=c(0,1,2,3), label=c('Neither', 'Suitability Gain', 'Suitability Loss', 'No change')))

#create figure for 2081-2100
par(mfrow=c(4,2))
plot(ensem_ssp126_2100,col= clr(200))
plot(change5, col=c('white', 'green', 'red', 'blue'), 
     breaks=c(-0.5,0.5,1.5,2.5, 3.5), 
     axis.args=list(at=c(0,1,2,3), label=c('Neither', 'Suitability Gain', 'Suitability Loss', 'No change')))
plot(ensem_ssp245_2100)
plot(change6, col=c('white', 'green', 'red', 'blue'), 
     breaks=c(-0.5,0.5,1.5,2.5, 3.5), 
     axis.args=list(at=c(0,1,2,3), label=c('Neither', 'Suitability Gain', 'Suitability Loss', 'No change')))
plot(ensem_ssp370_2100)
plot(change7, col=c('white', 'green', 'red', 'blue'), 
     breaks=c(-0.5,0.5,1.5,2.5, 3.5), 
     axis.args=list(at=c(0,1,2,3), label=c('Neither', 'Suitability Gain', 'Suitability Loss', 'No change')))
plot(ensem_ssp585_2100)
plot(change8, col=c('white', 'green', 'red', 'blue'), 
     breaks=c(-0.5,0.5,1.5,2.5, 3.5), 
     axis.args=list(at=c(0,1,2,3), label=c('Neither', 'Suitability Gain', 'Suitability Loss', 'No change')))




