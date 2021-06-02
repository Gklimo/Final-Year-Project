#load packages
library(readr)
library(tidyr)
library(dplyr)
plot(variables_c[[20]])

stringsAsFactors = FALSE
#read in Aedes aegypti data
aedes_aegypti <- read_csv("aedes_aegypti_data.csv")
nrow(aedes_aegypti)
#filter out species data so it is presence only - ie no 0 values in specimens collected
aedes_aegypti<- aedes_aegypti%>%
  filter(Specimens_collected>0)
locations<- aedes_aegypti%>%
  distinct(Locations)
#file with availiable locations
write_csv(locations,"locations.csv")
#load sdm package
library(sdm)
installAll() #only first time after installing sdm package
#create coordinates dataframe
aegypti_coords <- aedes_aegypti%>%
  dplyr::select(Longitudes,Latitudes)
#add a column with boolean indication of presence (all will be 1 as 0 values are ommited in previous step)
aegypti_coords$Species <- 1
#---------------------------------------------
#drop na values to avoid any errors
aegypti_coords <- aegypti_coords%>%
  drop_na()
#convert species data to spatial point dataframe
coordinates(aegypti_coords) <- c('Longitudes','Latitudes')
#load packages
library(raster)
library(dismo)
library(mapview)
#download bioclim data from WorldClim
#read in data from WorldClim
bio <- list.files(path = "wc2.1_2.5m_bio",pattern = '.tif', full.names = TRUE)
elevation <- list.files(pattern = "wc2.1_2.5m_elev.tif", full.names = TRUE)
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
#check for multicolinearity
#load package
library(usdm)

#vifstep is used for checking bias inflection factor and if it is greater than 10 then it is a sign of multicolinearity - stepwize procedure and checks vif and identifies those greater than 10 (by default it is 10) and if it is greater than the threashold it removes the highest score and goes through the process again and again until only variables below the threashold are present
#vifcor first checks correlation coefficient and if it is greater than threashold (eg 0.7 or 0.9(default))and it decides between the pair of variables which one should be excluded and excludes the one with higher VIF and then repeats the procedure until all variables scores are below the threashold
vif(variables_c)
#extract species data using raster from raster variable
extracts <- raster::extract(variables_c,aegypti_coords)
head(extracts)

#vifstep(variables_c)
v_step <- vifstep(extracts)
v_cor <- vifcor(extracts)
# VIFs of the remained variables -------- 
#Variables      VIF
#1    bio_12 6.238347
#2    bio_14 3.769333
#3    bio_17 6.212042
#4     bio_5 4.556403
#5     bio_8 4.161169
#6      elev 1.691756
v_step
v_cor#states that 14 variables from the 19 input variables have colinearity problem and says which ones should be kept and which excluded
#exclude variables with colinearity problem
variables_ex<- exclude(variables_c,v_step)#5 variables remain
variables_ex#predictor variables
#----------------------------------------
#creating data object


#specify a formula by species and predictor variables(bioc_ex), the train data is aegypti_coords

#sdmData(Species ~.,aegypti_coords, predictors = variables_ex)#number of records is only 105 because some could be duplicates or missing values in the predictor variables so they are automatically excluded 


#add pseudoabsence data as the type is presence-only to use other methods availible further
data_object <- sdmData(Species ~.,aegypti_coords, predictors = variables_ex, bg = list(method= "gRandom", n = 150))




#getmethodNames()#look at the availiable methods
#run sdm to fit the models --> generate 5 models per method (25 models for sub and 25 for boot) - do only 1: boot or sub later and change replicates number to 15
#model <- sdm(Species~., data = data_object, methods = c('glm','brt','rf','gam','maxent'), replication = c('sub','boot'),test.p = 30, n = 5, parallelSettings = list(ncore= 4, method = "parallel"))#4 model methods/test.p = test percentage/ n = 5 means that the replication of bootstapping or subsampling is repeated 5 times + parallel is to set the number of cores
model <- sdm(Species~., data = data_object, methods = c('glm','brt','rf','maxent'), replication = c('boot'),test.p = 30, n = 5, parallelSettings = list(ncore= 4, method = "parallel"))#4 model methods/test.p = test percentage/ n = 5 means that the replication of bootstapping or subsampling is repeated 5 times + parallel is to set the number of cores

#get model statistics
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


?roc
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

#results only for eg glm (can specify either id or method)
plot(getVarImp(model,method = "glm"), main = "GLM Relative Variable Importance")
plot(getVarImp(model,method = "gam"), main = "GAM Relative Variable Importance")
plot(getVarImp(model,method = "brt"), main = "BRT Relative Variable Importance")
plot(getVarImp(model,method = "rf"), main = "RF Relative Variable Importance")
plot(getVarImp(model,method = "maxent"), main = "Maxent Relative Variable Importance")





#look up the density function to generate the density plot
#extract details of the evaluation out of this model object getEvaluation()/getVarImp()/getRoc() for roc curves used

#check details inside the object to make it reproducible eg:
#uncomment below
#model@models$Species$rf$`13`@object



#graphical user interface of the package
gui(model)
model
#generates 50 outputs and generate 
prediction <- predict(model, variables_c, filename = "prediction_n.img", overwrite = TRUE, mean = T)
p <-raster::stack("prediction.img")
e <- raster("ensemble.tif")
par(mfrow=c(1,1))
plot(prediction)
plot(p[[1]])
plot(p[[2]])
plot(p[[3]])
plot(p[[4]])
plot(p[[5]])
#combines 50 outputs
ensem <- ensemble(model, variables_c, filename = "ensemble_n.tif", setting = list(method = "mean-weighted",stat = "TSS",opt =2 ))
ensem_auc <- ensemble(model, variables_c, filename = "ensemble_auc_n.img", setting = list(method = "mean-weighted",stat = "AUC" ))


par(mfrow=c(1,1))
plot(prediction)
par(mfrow=c(1,3))
plot(ensem,col=clr(200))
plot(ensem_auc,col = clr(200))
plot(ensem - ensem_auc,col = clr(200))

p_ssp126 <- predict(model, ssp126_2100_c, filename = "prediction_126.img", overwrite = TRUE, mean = T)
p_ssp245 <- predict(model, ssp245_2100_c, filename = "prediction_245.img", overwrite = TRUE, mean = T)
p_ssp370 <- predict(model, ssp370_2100_c, filename = "prediction_370.img", overwrite = TRUE, mean = T)
p_ssp585 <- predict(model, ssp585_2100_c, filename = "prediction_585.img", overwrite = TRUE, mean = T)
plot(p_ssp126)
plot(p_ssp245)
plot(p_ssp370)
plot(p_ssp585)
plot(prediction)
#####################################
#future predictions
#ssp245_2100 <- list.files(pattern = "ssp245_2100.tif", full.names = TRUE)
#ssp245_2100<- raster::stack(ssp245_2100)

ssp245_2100 <- list.files(path = "ssp245",pattern = '.tif', full.names = TRUE, recursive = T)
ssp245_2100<- raster::stack(ssp245_2100)
elevation <- list.files(pattern = "wc2.1_2.5m_elev.tif", full.names = TRUE)
elevation <- raster(elevation)
ssp245_2100 <- addLayer(ssp245_2100,elevation)

#change the names so they are the same as names in variables
names(ssp245_2100) <- names(variables)

#crop to same coordinates
ssp245_2100_c <- crop(ssp245_2100, box)
#ssp585 is the predictor variables so use is as an argument in ensemble models
ensem_ssp245_2100<- ensemble(model, ssp245_2100_c, filename = "ssp245_2100.img", setting = list(mehtod = "mean-weighted", stat = "AUC"))

#generate color pallet
clr <- colorRampPalette(c("#3E49BB","#349808","yellow","orange","red","darkred"))

plot(ensem,col = clr(6))
plot(ensem_ssp245_2100,col = clr(6))
#quantify the changes quantitatively
#calculate suitability --> negative or positive values. (-ve)  = less suitability in the future compared to the current time (suitability declined), (+ve) = more suitability in the future so suitability increased for the species
change <- ensem_ssp245_2100 - ensem

clr2 <- colorRampPalette(c("red","yellow","green","blue"))

plot(change,  col = clr2(4))


writeRaster(ensem,"ensemble_tss.tif")
writeRaster(ensem_ssp245_2100,"ensem_245_100.tif")

#####################
#calculate area
plot(ensem,col = clr(200))
eval<- getEvaluation(model, stat = c("threshold"),opt =4)

mean_threshold <- mean(eval$threshold)
mean_threshold <-0.568
now <- raster("ensemble.tif")

future <- raster("ssp126_2100.tif")

now_bin <- now >= mean_threshold
future_bin  <- future >= mean_threshold
# Create a coding that shows gains and losses
change <- 2 * now_bin + future_bin

par(mfrow=c(1,3))
plot(now_bin)
plot(future_bin)
plot(change, col=c('white', 'green', 'red', 'blue'), 
     breaks=c(-0.5,0.5,1.5,2.5, 3.5), 
     axis.args=list(at=c(0,1,2,3), label=c('Neither', 'Suitability Gain', 'Suitability Loss', 'No change')))



plot(ensem_ssp245_2100,col = clr(200))
plot(ensem_ssp126_2100)

par(mar=c(3,3,1,1), mgp=c(2,1,0))


gui(model)
s