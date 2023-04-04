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
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_analysis/Analysis")

## Paths

#Phenology
PEPcsv_path = "PEP_data/PEPcsv"

#Output
PEP_drivers_path  = "Analysis_input/Drivers"


## load phenology data
PEP.df <- fread(paste(PEPcsv_path, "pepData.csv", sep="/"))


##############################################################################################################################################


#################
## Photoperiod ##
#################


# Get all time-points
site <- unique(PEP.df$pep_id)

# Initialize data frame to store results
photo.df <- data.frame()
i=1

for(id in site) {

  # Subset table according to latitude
  photo.sub <- PEP.df %>% 
    filter(pep_id==id) %>% 
    dplyr::select(pep_id,lat)
  photo.sub <- photo.sub[!duplicated(photo.sub$pep_id),]
  
  # Calculate daily photoperiod for the whole year
  photo <- daylength(photo.sub$lat,1:366)
  
  # Add daily photoperiods to the subset table
  photo.sub[as.character(1:366)] <- 0
  photo.sub[,3:368] <- photo
  
  photo.df <- rbind(photo.df,photo.sub)
  print(paste0(round(i/length(site)*100, 1)," % of photoperiods calculated!"))
  i=i+1
}
rm(photo.sub)

# Export dataset
write.table(photo.df,paste0(PEP_drivers_path ,"/Photoperiod.csv"),sep=",",row.names=FALSE)


##############################################################################################################################################
