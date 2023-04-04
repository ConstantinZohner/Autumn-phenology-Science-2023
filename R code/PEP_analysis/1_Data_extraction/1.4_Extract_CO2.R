require(data.table)
require(tidyverse)

##########################################
## Set directories and get PEP data ##
##########################################

# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_analysis/Analysis")

## Paths

#CO2
CO2_path  = "Analysis_input/Drivers/CO2"

#Phenology
PEPcsv_path = "PEP_data/PEPcsv"


## load phenology data
PEP.df <- fread(paste(PEPcsv_path, "pepData.csv", sep="/"))

#load mean of monthly atmospheric CO2 (mole fraction)
CO2.df = fread(paste(CO2_path, "CO2_monthly.csv", sep="/"))

##############################################################################################################################################

##############
## CO2 data ##
##############

#add year 2015 (= 2014 value + 2, based on Mauna Loa observations, https://www.esrl.noaa.gov/gmd/ccgg/trends/data.html)
CO2_2015.df = CO2.df %>%
  filter(year %in% 2014) %>% #keep only 2014
  mutate(data = data+2, #add 2 ppm
         year = as.integer(2015))
CO2.df = rbind(CO2.df, CO2_2015.df)
rm(CO2_2015.df)

#filter latitude and year
CO2.df = CO2.df %>%
  filter(lat %in% 52.5 & 
           year %in% unique(PEP.df$year)) %>%
  rename(CO2=data)%>% #rename columns
  dplyr::select(-c(datetime, lat))

#PEP site x year data
PEP.df = PEP.df %>% 
  #add site x year identifier to PEP data
  mutate(site_year=paste0(PEP.df$pep_id,"_",PEP.df$year)) %>%
  #delete duplicates
  filter(!duplicated(site_year)) 

#create new dataframe
PEP_CO2.df <- data.frame()

# Loop through observation years
for(yr in min(PEP.df$year):max(PEP.df$year)) {
  
  # Subset by year
  pheno.sub  <- PEP.df %>% 
    filter(year == yr)
  CO2.sub <- CO2.df %>% 
    filter(year == yr)
  
  # Transpose CO2 monthly values to the phenological subset
  pheno.sub[as.character(1:12)] <- 0
  for(r in 1:nrow(pheno.sub)) {
    pheno.sub[r,as.character(1:12)] <- t(CO2.sub$CO2)
  }
  
  # Bind final dataset
  PEP_CO2.df <- rbind(PEP_CO2.df,pheno.sub)
  
  print(yr)
}

# Data wrangling 
PEP_CO2.df <- PEP_CO2.df %>% 
  arrange(pep_id,year) %>% 
  select(pep_id, year, lat, lon, site_year, as.character(1:12))

#Check CO2 data
CO2.df$date = as.Date(paste0(CO2.df$year,"_15_",CO2.df$month),"%Y_%d_%m")
ggplot(data=CO2.df, aes(x=date, y=CO2)) +
  geom_line()+
  xlab("Year") +
  ylab("Atmospheric [CO2] (ppm)") 

# remove stuff
rm(CO2.sub, CO2.df, PEP.df, pheno.sub)

#safe CO2 table
write.csv(PEP_CO2.df, paste(CO2_path, "CO2.csv", sep="/"))

##############################################################################################################################################
