# Load libraries
require(data.table)
require(geosphere)
require(zoo)
require(lubridate)
require(tidyverse)


##############################################################################################################################################


######################################
## Set directories and get PEP data ##
######################################


# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/Remote_sensing/Analysis")

## Paths

#Phenology
phenology_path = "Analysis_input/Phenology_data"

#Output
Drivers_path  = "Analysis_input/Drivers"


## load phenology data (Rbind files)
pheno.df = fread(paste(phenology_path, "PhenologyData_VNP_2013_2021.csv", sep="/")) 


##############################################################################################################################################


#################
## Photoperiod ##
#################


# Get all time-points
site <- unique(pheno.df$geometry)

# Initialize data frame to store results
photo.df <- data.frame()
i=1

for(id in site) {

  # Subset table according to latitude
  photo.sub <- as.data.frame(pheno.df %>% 
    filter(geometry==id) %>% 
    dplyr::select(geometry, Lat) %>%
    distinct(Lat, .keep_all = T)) #delete duplicates
  
  # Calculate daily photoperiod for the whole year
  photo <- geosphere::daylength(photo.sub$Lat,1:366)
  
  # Add daily photoperiods to the subset table
  photo.sub[as.character(1:366)] <- 0
  photo.sub[,3:368] <- photo
  
  photo.df <- rbind(photo.df,photo.sub)
  print(paste0(round(i/length(site)*100, 1),"% of photoperiods calculated!"))
  i=i+1
}
rm(photo.sub)

# Export dataset
write.table(photo.df, paste0(Drivers_path, "/Photoperiod_VNP.csv"), sep=",", row.names=FALSE)


##############################################################################################################################################
