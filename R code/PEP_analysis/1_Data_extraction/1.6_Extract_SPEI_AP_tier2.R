require(data.table)
require(lubridate)
require(SPEI)
require(tidyverse)

##############################################################################################################################################

##########################################
## Set directories and get own PEP data ##
##########################################

# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/Analysis")

## Paths

#1. Input
GLDAS_path  = "Analysis_input/Drivers/GLDAS/GLDAS_Ordered"

#2. Output
SPEI.df_path  = "Analysis_input/Drivers"

##############################################################################################################################################

#################
## Import data ##
#################

## Import daily climatic datasets from GLDAS
vn <- c('Daily_Mean_Data_Tair_f_inst',
        'Daily_Data_Rainf_f_tavg')
DataList <- replicate(length(vn),data.frame())
for(i in 1:length(vn)) {
  data <- fread(paste0(GLDAS_path, "/", vn[i],".csv"))
  DataList[[i]] <- data
}
names(DataList)=vn

## Unit conversions
# Precipitation estimates are given as rate in kg m-2 s-1. We need mm d-1.
# 1 kg of rain water spread over 1 square meter of surface is 1 mm in thickness
# there are 60 X 60 X 24 = 86400 seconds in one day. Therefore, 1 kg m-2 s-1 = 86400 mm d-1
#DataList[[4]][,as.character(1:366)]=DataList[[4]][,as.character(1:366)]*86400

##############################################################################################################################################

################################ 
## Get SPEI and water balance ##
################################

Tmean.df = DataList[[1]] %>%
  #long format
  pivot_longer(., -c(pep_id, year, site_year, lat, lon)) %>%
  #rename
  rename(Tmean=value)%>%
  #delete NAs
  filter(!is.na(Tmean)) %>%
  #create month-year identifier
  mutate(doy = as.numeric(name)-1,
         date = month(as.Date(doy, origin = paste0(year,"-01-01")))) %>%
  #delete columns
  select(pep_id, doy, date, year, Tmean, site_year, lat)

Prcp.df = DataList[[2]] %>%
  #long format
  pivot_longer(., -c(pep_id, year, site_year, lat, lon)) %>%
  #rename
  rename(Prcp=value) %>%
  #delete NAs
  filter(!is.na(Prcp)) %>%
  #create month-year identifier
  mutate(doy = as.numeric(name)-1,
         date = month(as.Date(doy, origin = paste0(year,"-01-01"))))%>%
  #delete columns
  select(pep_id, doy, date, year, Prcp)

#merge dataframes
daily.df = inner_join(Tmean.df, Prcp.df, by=c('pep_id','doy','date','year'))

#get monthly values
SPEI.df = as.data.frame(daily.df %>% 
                          group_by(site_year,pep_id,year,date,lat) %>% 
                          summarize(Tmonth=mean(Tmean),
                                    Prcp=sum(Prcp))%>% 
                          mutate(month=as.numeric(gsub(".*-","",date))) %>% 
                          ungroup())

#get SPEI, water balance and PET
SPEI.df = as.data.frame(SPEI.df %>%
                          group_by(pep_id) %>% 
                          mutate(PET = thornthwaite(Tmonth, lat[1]),
                                 WB = Prcp-PET,
                                 SPEI = spei(Prcp-PET,1)$fitted))%>%
  select(c(site_year,year,month,SPEI,PET,WB))

#Safe table
write.csv(SPEI.df, paste(SPEI.df_path, "SPEI.csv", sep="/"))

##############################################################################################################################################

