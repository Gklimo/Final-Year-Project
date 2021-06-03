
#start the analysis here 

library(dplyr)
library(tidyr)
library(readr)
# set strings as factors to false
options(stringsAsFactors = FALSE)
#world abundance 

data <- read_csv("world_abundance_data.csv")
#check class
class(data)
#check dimensions
dim(data)
#convert data to data.frame class
#data <- data%>%
# as_data_frame()
#filter out Aedes aegypti data
aedes_aegypti <- data%>%
  select("Sample ID","Species","Collection date range","Locations","Latitudes","Longitudes","Specimens collected","Sex","Developmental stage")%>%
  rename("Sample_ID" = "Sample ID","Developmental_stage" = "Developmental stage", "Date_collected"= "Collection date range","Specimens_collected"= "Specimens collected")%>%
  filter(grepl("aegypti",Species), Sex == "female",`Developmental_stage`== "adult")%>%
  mutate(Date_collected = as.Date(Date_collected))

#create csv
write_csv(aedes_aegypti,"aedes_aegypti_data.csv")