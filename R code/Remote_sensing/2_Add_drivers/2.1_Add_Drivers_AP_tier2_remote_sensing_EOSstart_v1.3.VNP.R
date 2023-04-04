


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Climate driver extraction of the remote sensing analysis (Senescence onset VNP) ###########################
#############################################################################################################



#required packages
require(data.table)
require(sf)
require(ncdf4)
require(raster)
require(tidyverse)
require(sp)
require(rpmodel)
require(purrr)
require(pbmcapply)
require(zoo)
require(chillR)
require(lubridate)
require(weathermetrics)
require(rgdal)



##############################################################################################################################################
##############################################################################################################################################



#############################
## Set directory and paths ##
#############################



# Set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2")


#########
# Paths #
#########


# 1. Input
##########

#Climate
GLDAS_path       = "Remote_sensing/Analysis/Analysis_input/Drivers/GLDAS"
GLDAS19_21_path  = "Remote_sensing/Analysis/Analysis_input/Drivers/GLDAS19_21"

#CO2
CO2_path         = "Remote_sensing/Analysis/Analysis_input/Drivers/CO2"

#Photoperiod
photo_path       = "Remote_sensing/Analysis/Analysis_input/Drivers"

# GPP and LAI
GPP_path         = "Remote_sensing/Analysis/Analysis_input/Drivers/Modis_GPP_LAI"

#Phenology
Pheno_path       = "Remote_sensing/Analysis/Analysis_input/Phenology_data"


# 2. Output
###########

Drivers_path  = "Remote_sensing/Analysis/Analysis_input/Drivers_final_onset_VNP/Individual_files"
Drivers_path2 = "Remote_sensing/Analysis/Analysis_input/Drivers_final_onset_VNP/Merged_file"
Drivers_path3 = "Remote_sensing/Analysis/Analysis_input/Drivers_final_onset_VNP/Missing_observations"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



## Phenology data
#################

Pheno.df <- fread(paste(Pheno_path, "PhenologyData_VNP_2013_2021.csv", sep="/")) %>% 
  filter(Lat>20,
         Onset_Greenness_Increase<171,
         Onset_Greenness_Decrease>180,
         Onset_Greenness_Decrease<280,
         Onset_Greenness_Decrease<Date_Mid_Senescence) %>% 
  group_by(geometry) %>% 
  #delete pixels with less than 9 years
  filter(n() >= 9) %>%
  #get autumn phenology means per pixel
  mutate(Onset_Mean   = mean(Onset_Greenness_Decrease)) %>%
  ungroup() %>%
  #delete duplicates
  distinct(geometry, Year, .keep_all = T) 
  

##############################################################################################################################################


## Elevation map
################

elev.raster = raster(paste(photo_path, "topo_elevation.asc", sep="/"))


## Biomes map
#############

biome.raster = raster(paste(photo_path, "WWF_Biomes_HalfDegree.tif", sep="/"))


## CO2 data
###########

CO2.df = fread(paste(CO2_path, "CO2_Annual.csv", sep="/")) %>% 
  add_row(Year=c(2019,2020,2021), CO2=c(411.7,414.2,416.5))


## Photoperiod
##############

photo.df = fread(paste(photo_path, "Photoperiod_VNP.csv", sep="/"))


##############################################################################################################################################


## GPP
######

GPP.df = fread(paste(GPP_path, "GppData_Forests_025.csv", sep='/')) %>%
  #keep only numbers in string
  mutate(geometry2 = gsub("POINT ","",  geometry),
         geometry2 = gsub("\\(|\\)","", geometry2)) %>%
  separate(geometry2, into = c("Lat","Lon"), sep=" ") %>%
  mutate(Lat = as.numeric(Lat),
         Lon = as.numeric(Lon)) %>%
  #add site x year identifier
  mutate(site_year = paste0(geometry, '_', Year)) %>% 
  #order table
  dplyr::select(Year, geometry, site_year, Lat, Lon, everything())
#rename columns
colnames(GPP.df)[6:ncol(GPP.df)] <- seq(1, 366, by=8)

# get only GPP values
GPP = GPP.df %>% 
  dplyr::select(-c(Year, geometry, site_year, Lat, Lon)) 
#Get only info
Info = GPP.df %>% 
  dplyr::select(c(Year, geometry, site_year, Lat, Lon)) 
#Get empty data frame to store duplicated data
GPPfinal = tibble(.rows=nrow(GPP))

#Loop to duplicate columns
for (i in 1:8) {
  GPPsub = GPP
  colnames(GPPsub)[1:ncol(GPPsub)] <- seq(i, 368, by=8)
  GPPfinal = cbind(GPPfinal, GPPsub)
}

#sort column names
GPPfinal = GPPfinal %>% select(str_sort(names(.), numeric=T)) 

#divide GPP values by 8 to obtain daily GPP
GPPfinal = GPPfinal[,c(1:366)]/8

#Add site info
GPPfinal = cbind(Info, GPPfinal)

#load 2019 to 2021 data
GPP19_21.df = fread(paste(GPP_path, "GPP_MODIS_025_2019-21.csv", sep='/')) %>% 
  mutate(site_year = paste(geometry, Year, sep="_")) %>% 
  dplyr::select(c(Year, geometry, site_year, Lat, Lon, everything())) %>% 
  filter(Lat>20)
#add column names
colnames(GPP19_21.df)[6:ncol(GPP19_21.df)] <- seq(1, 366, by=1)
#scale units to kg*C/m^2	
GPP19_21.df[,6:ncol(GPP19_21.df)] = GPP19_21.df[,6:ncol(GPP19_21.df)]/0.0001

#merge data
GPP.df = rbind(GPPfinal, GPP19_21.df)


##############################################################################################################################################


## Import daily climatic datasets from GLDAS
############################################

#define climate variables
vn <- c('GLDAS_Daily_Data_Tair_f_inst_Mean',
        'GLDAS_Daily_Data_Tair_f_inst_Min',
        'GLDAS_Daily_Data_Tair_f_inst_Max',
        'GLDAS_Daily_Data_Rainf_f_tavg',
        'GLDAS_Daily_Data_Qair_f_inst',
        'GLDAS_Daily_Data_SoilMoi0_10cm_inst',
        'GLDAS_Daily_Data_SoilMoi10_40cm_inst',
        'GLDAS_Daily_Data_Swnet_tavg')

vn2 <- c('GLDAS_Daily_Data_of_Tair_f_inst_mean_2019-21',
         'GLDAS_Daily_Data_of_Tair_f_inst_min_2019-21',
         'GLDAS_Daily_Data_of_Tair_f_inst_max_2019-21',
         'GLDAS_Daily_Data_of_Rainf_f_tavg_2019-21',
         'GLDAS_Daily_Data_of_Qair_f_inst_2019-21',
         'GLDAS_Daily_Data_of_SoilMoi0_10cm_inst_2019-21',
         'GLDAS_Daily_Data_of_SoilMoi10_40cm_inst_2019-21',
         'GLDAS_Daily_Data_of_Swnet_tavg_2019-21')

#create empty list
DataList <- replicate(length(vn),data.frame())

#loop through climate variables
for(i in 1:length(vn)) {

  #read data
  ##########
  
  data = fread(paste0(GLDAS_path, "/", vn[i],".csv")) %>%
    #keep only numbers in string
    mutate(geometry2 = gsub("POINT ","",  geometry),
           geometry2 = gsub("\\(|\\)","", geometry2)) %>%
    separate(geometry2, into = c("Lat","Lon"), sep=" ") %>%
    mutate(Lat = as.numeric(Lat),
           Lon = as.numeric(Lon)) %>%
    #add site x year identifier
    mutate(site_year = paste0(geometry, '_', Year)) %>% 
    #order table
    dplyr::select(Year, geometry, site_year, Lat, Lon, everything())
  
  #rename columns
  colnames(data)[6:ncol(data)] <- as.numeric(1:366)
  
  #delete NAs
  data = data %>% filter(!is.na(`170`))  
  
  #----------------------
  
  #read 2019 to 2021 data
  #######################
  
  data2 = fread(paste0(GLDAS19_21_path, "/", vn2[i],".csv")) %>%
    #keep only numbers in string
    mutate(geometry2 = gsub("POINT ","",  geometry),
           geometry2 = gsub("\\(|\\)","", geometry2)) %>%
    separate(geometry2, into = c("Lat","Lon"), sep=" ") %>%
    mutate(Lat = as.numeric(Lat),
           Lon = as.numeric(Lon)) %>%
    #add site x year identifier
    mutate(site_year = paste0(geometry, '_', Year)) %>% 
    #order table
    dplyr::select(Year, geometry, site_year, Lat, Lon, everything())
  
  #rename columns
  colnames(data2)[6:ncol(data2)] <- as.numeric(1:366)
  
  #delete NAs
  data2 = data2 %>% filter(!is.na(`170`))  
  
  #----------------------
  
  #merge tables
  data = rbind(data, data2)

  #add table to list
  DataList[[i]] <- data
}
#add names to list
names(DataList)=vn
# Note: Precipitation is given as rate in mm d-1.



##############################################################################################################################################
##############################################################################################################################################



######################################
## biome information and elevation  ##
######################################



#Add biome information and elevation 
Pheno.df = 
  #both tables together
  cbind(Pheno.df, 
        # intersection
        data.frame(alt     = raster::extract(elev.raster,    Pheno.df[, c("Lon", "Lat")])),
        data.frame(biome   = raster::extract(biome.raster,   Pheno.df[, c("Lon", "Lat")])) )

#remove stuff
rm(elev.raster, biome.raster)



##############################################################################################################################################
##############################################################################################################################################



#######################################################
## Calculate climatic predictors using Parallel calc ##
#######################################################



# Identifier (all site x year combinations)
Pheno.df$site_year = paste0(Pheno.df$geometry, '_', Pheno.df$Year)
timeseries_year    = unique(Pheno.df$site_year)
  
# add Pheno, CO2 and photoperiod data to list
DataList[[9]] = photo.df
DataList[[10]] = CO2.df
DataList[[11]] = Pheno.df
DataList[[12]] = GPP.df

rm(photo.df, CO2.df, data, Pheno.df, GPP.df)
names(DataList)=c(vn,"photoperiod",'CO2',"Pheno","GPP")
names(DataList)
#[1] "GLDAS_Daily_Data_Tair_f_inst_Mean"    "GLDAS_Daily_Data_Tair_f_inst_Min"     "GLDAS_Daily_Data_Tair_f_inst_Max"    
#[4] "GLDAS_Daily_Data_Rainf_f_tavg"        "GLDAS_Daily_Data_Qair_f_inst"         "GLDAS_Daily_Data_SoilMoi0_10cm_inst" 
#[7] "GLDAS_Daily_Data_SoilMoi10_40cm_inst" "GLDAS_Daily_Data_Swnet_tavg"          "photoperiod"                         
#[10] "CO2"                                  "Pheno"                                "GPP"         


##############################################################################################################################################


################################
# Loop through all time-points #
################################


parallelCalc <- function(timeseries_years){ 

  # Subset input data by time-point
  #################################

  #phenology data
  pheno.sub  <- DataList[[11]][which(DataList[[11]]$site_year==timeseries_years),]
  
  #daily mean temperature
  TMEAN <- DataList[[1]][which(DataList[[1]]$site_year==pheno.sub$site_year),]%>% 
    dplyr::select(as.character(1:366))
  
  #Skip timeseries for which there is no data
  if (nrow(TMEAN)==0) {
    write.table(pheno.sub, file=paste0(Drivers_path3, '/', timeseries_years, '.csv'), sep=',', row.names = F, col.names = T)
    } else {
      
    #daily minimum temperature
    TMIN      <- DataList[[2]][which(DataList[[2]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #daily maximum temperature
    TMAX      <- DataList[[3]][which(DataList[[3]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #precipitation
    PRCP      <- DataList[[4]][which(DataList[[4]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #soil moisture (10-40 cm)
    MOIST   <- DataList[[7]][which(DataList[[7]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #net short-wave radiation
    SWRAD     <- DataList[[8]][which(DataList[[8]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #day length
    PHOTO     <- DataList[[9]][which(DataList[[9]]$geometry==pheno.sub$geometry),][1]%>% 
      dplyr::select(as.character(1:366))
    
    #CO2
    CO2       <- DataList[[10]][which(DataList[[10]]$Year==pheno.sub$Year),]$CO2
    
    #GPP
    GPP <- DataList[[12]][which(DataList[[12]]$site_year==pheno.sub$site_year),] %>%
      dplyr::select(as.character(1:366))
    
    
    ##############################################################################################################################################
    
    
    # Create table of daily climate
    ###############################
    
    # Generate sub-dataframe to store results
    factors.sub <- pheno.sub %>% 
      dplyr::select(geometry, Lat, Lon, alt, Year, Onset_Greenness_Increase, Onset_Greenness_Decrease, Date_Mid_Senescence) %>%
      mutate(CO2 = CO2)
    
    # Define the current year in calendar units
    year      <- as.character(pheno.sub$Year)
    start_doy <- paste(year,"-01-01", sep="") 
    end_doy   <- paste(year,"-12-31", sep="")
    days      <- seq(as.Date(start_doy), as.Date(end_doy), by="days")
    
    #create table
    daily_vals <- data.frame(Year     = year,
                             Month    = 0,
                             Day      = 0,
                             Tmin     = as.numeric(TMIN), 
                             Tmean    = as.numeric(TMEAN), 
                             Tmax     = as.numeric(TMAX), 
                             SWrad    = as.numeric(SWRAD),
                             Moist    = as.numeric(MOIST),
                             Prcp     = as.numeric(PRCP), 
                             Photo    = as.numeric(PHOTO),
                             GPP      = as.numeric(GPP))

    #Add climate variables and data wrangling
    daily_vals = daily_vals %>%
      filter(!is.na(Tmean)) %>%#delete NAs
      mutate(
        #add month and day identifiers
        Month = lubridate::month(as.Date(days,origin=days[1])),
        Day   = lubridate::day(as.Date(days,origin=days[1])) )
      
    #set NAs to 0
    daily_vals[is.na(daily_vals)] <- 0.0001
    
    
    ##############################################################################################################################################
    
    
    # Get average daytime temperature (chillR package)
    ##################################################
    
    #Get hourly values
    hourly_vals = stack_hourly_temps(daily_vals, latitude=pheno.sub$Lat)$hourtemps
    
    #get daytime temperature
    daytime_temp = data.frame(hourly_vals %>%
                                 group_by_at(vars(-c(Temp,Hour))) %>%
                                 #delete night hours
                                 filter(
                                   #filter only if daylength is less than 24 hours
                                   if(daylength(latitude=pheno.sub$Lat,JDay=JDay[1])$Daylength < 24) 
                                   between(Hour,daylength(latitude=pheno.sub$Lat,JDay=JDay[1])$Sunrise,
                                                daylength(latitude=pheno.sub$Lat,JDay=JDay[1])$Sunset) 
                                   else TRUE
                                   ) %>% 
                                 #summarise daytime hours
                                 summarise(Tday = mean(Temp))%>%
                                 ungroup() %>%
                                 #order
                                 dplyr::select(Tday) )
    
    #get nighttime temperature
    nighttime_temp = data.frame(hourly_vals %>%
                                  group_by_at(vars(-c(Temp,Hour))) %>%
                                  #delete day hours
                                  filter(
                                    #filter if daylength is less than 24 hours
                                    if(daylength(latitude=pheno.sub$Lat,JDay=JDay[1])$Daylength < 24) 
                                    !between(Hour,daylength(latitude=pheno.sub$Lat,JDay=JDay[1])$Sunrise,
                                                  daylength(latitude=pheno.sub$Lat,JDay=JDay[1])$Sunset)
                                    #if no darkness select minimum Temp
                                    else Temp == min(Temp)
                                    )%>%
                                  #summarise nighttime hours
                                  summarise(Tnight = mean(Temp))%>%
                                  ungroup() %>%
                                  #order
                                  dplyr::select(Tnight) )
    
    #combine
    daily_vals = cbind(daily_vals, daytime_temp, nighttime_temp) %>%
      #order
      dplyr::select(Year, Month, Day, Tmin, Tmean, Tmax, Tday, Tnight, everything()) 


    ##############################################################################################################################################
    
    
    # Get important dates
    #####################
    
    # warmest day of year
    factors.sub$HottestDOY = mean(which(daily_vals$Tmax==max(daily_vals$Tmax)))
    
    # day of maximum radiation
    factors.sub$MaxRadDOY = mean(which(daily_vals$SWrad==max(daily_vals$SWrad)))
    
    # longest day of year (summer solstice)
    solstice    = which(daily_vals$Photo==max(daily_vals$Photo))[1] 
    
    # March equinox
    equinox.Mar = solstice - 97     
    
    # September equinox
    equinox.Sep = solstice + 97                                    
    
    # Mean leaf senescence
    DOY_off <- round(pheno.sub$Onset_Mean)
    
    # leaf-out
    DOY_out <- pheno.sub$Onset_Greenness_Increase
    
    
    ##############################################################################################################################################
    
    
    # Set GPP and LAI before greenup to zero
    ########################################
    
    
    #GPP
    daily_vals$GPPstart = daily_vals$GPP
    daily_vals$GPPstart[1:DOY_out]=0
    
      
      ##############################################################################################################################################
      
     
      ###################
      ## Store drivers ##
      ###################
      
      
      ######################
      ## Seasonal drivers ##
      ######################
      
      #define variables
      variable.names = c('GPPstart','Tday', 'Tnight','SWrad','Moist','Prcp')
      
      #---------------------------------------------------------------------------------------------------------
      
      for(i in 1:length(variable.names)) {

        #choose variable (daily values)
        variable = daily_vals[,variable.names[i]]
        
        #---------------------------------------------------------------------------------------------------------
        
        # Name variables
        ################
        
        # Seasonal
        ##########
        
        # LO...leaf-out date
        # SE...mean senescence date
        # SO...Summer solstice (~22 June)
        # SOm30...Summer solstice -30 (~22 May)
        # SOp30...Summer solstice +30 (~21 July)
        # SOp60...Summer solstice +60 (~22 August)
        varname.LO.SO       <- paste(variable.names[i], "LO.SO",    sep=".")
        varname.LO.SOm30    <- paste(variable.names[i], "LO.SOm30", sep=".")
        varname.LO.SOp30    <- paste(variable.names[i], "LO.SOp30", sep=".")
        varname.LO.SOp60    <- paste(variable.names[i], "LO.SOp60", sep=".")
        varname.LO.SE       <- paste(variable.names[i], "LO.SE",    sep=".")
        varname.SO.SE       <- paste(variable.names[i], "SO.SE",    sep=".")
        varname.SOm30.SE    <- paste(variable.names[i], "SOm30.SE", sep=".")
        varname.SOp30.SE    <- paste(variable.names[i], "SOp30.SE", sep=".")
        varname.SOp60.SE    <- paste(variable.names[i], "SOp60.SE", sep=".")
       
        # Solstice
        ##########
        
        # solstice1...sum of 40 to 10 days before solstice
        # solstice2...sum of 30 to 0 days before solstice
        # solstice3...sum of 20 days before to 10 days after solstice
        # solstice4...sum of 10 days before to 20 days after solstice
        # solstice5...sum of 0 to 30 days after solstice
        # solstice6...sum of 10 to 40 days after solstice
        varname.solstice1   <- paste(variable.names[i], "solstice1", sep=".")
        varname.solstice2   <- paste(variable.names[i], "solstice2", sep=".")
        varname.solstice3   <- paste(variable.names[i], "solstice3", sep=".")
        varname.solstice4   <- paste(variable.names[i], "solstice4", sep=".")
        varname.solstice5   <- paste(variable.names[i], "solstice5", sep=".")
        varname.solstice6   <- paste(variable.names[i], "solstice6", sep=".")
        
        #---------------------------------------------------------------------------------------------------------
        
        # Create columns
        ################
        
        if(variable.names[i] %in% c('GPPstart','LAIstart')){
          
          # Sums from leaf-out
          ####################
          
          factors.sub = factors.sub %>%
            mutate(
              #seasonal
              !!varname.LO.SO    := sum(variable[DOY_out:solstice]),
              !!varname.LO.SOm30 := sum(variable[DOY_out:(solstice-30)]),
              !!varname.LO.SOp30 := sum(variable[DOY_out:(solstice+30)]),
              !!varname.LO.SOp60 := sum(variable[DOY_out:(solstice+60)]),
              !!varname.LO.SE    := sum(variable[DOY_out:DOY_off]),
              !!varname.SO.SE    := sum(variable[solstice:DOY_off]),
              !!varname.SOm30.SE := sum(variable[(solstice-30):DOY_off]),
              !!varname.SOp30.SE := sum(variable[(solstice+30):DOY_off]),
              !!varname.SOp60.SE := sum(variable[(solstice+60):DOY_off]),
              
              #solstice
              !!varname.solstice1 := sum(variable[(solstice-39):(solstice-10)]),
              !!varname.solstice2 := sum(variable[(solstice-29):solstice]),
              !!varname.solstice3 := sum(variable[(solstice-19):(solstice+10)]),
              !!varname.solstice4 := sum(variable[(solstice-9):(solstice+20)]),
              !!varname.solstice5 := sum(variable[(solstice+1):(solstice+30)]),
              !!varname.solstice6 := sum(variable[(solstice+11):(solstice+40)]) )
        } 
        
        if(variable.names[i] %in% c('Tday','Tnight','Moist','SWrad')){
          
          # Means from fixed date
          #######################
          
          factors.sub = factors.sub %>%
            mutate(
              #seasonal
              !!varname.LO.SO    := mean(variable[equinox.Mar:solstice]),
              !!varname.LO.SOm30 := mean(variable[equinox.Mar:(solstice-30)]),
              !!varname.LO.SOp30 := mean(variable[equinox.Mar:(solstice+30)]),
              !!varname.LO.SOp60 := mean(variable[equinox.Mar:(solstice+60)]),
              !!varname.LO.SE    := mean(variable[equinox.Mar:DOY_off]),
              !!varname.SO.SE    := mean(variable[solstice:DOY_off]),
              !!varname.SOm30.SE := mean(variable[(solstice-30):DOY_off]),
              !!varname.SOp30.SE := mean(variable[(solstice+30):DOY_off]),
              !!varname.SOp60.SE := mean(variable[(solstice+60):DOY_off]),
              
              #solstice
              !!varname.solstice1 := mean(variable[(solstice-39):(solstice-10)]),
              !!varname.solstice2 := mean(variable[(solstice-29):solstice]),
              !!varname.solstice3 := mean(variable[(solstice-19):(solstice+10)]),
              !!varname.solstice4 := mean(variable[(solstice-9):(solstice+20)]),
              !!varname.solstice5 := mean(variable[(solstice+1):(solstice+30)]),
              !!varname.solstice6 := mean(variable[(solstice+11):(solstice+40)]) )
        }
         
        if(variable.names[i] %in% c('Prcp')){
          
          # Sums from fixed date
          ######################
          
          factors.sub = factors.sub %>%
            mutate(
              #seasonal
              !!varname.LO.SO    := sum(variable[equinox.Mar:solstice]),
              !!varname.LO.SOm30 := sum(variable[equinox.Mar:(solstice-30)]),
              !!varname.LO.SOp30 := sum(variable[equinox.Mar:(solstice+30)]),
              !!varname.LO.SOp60 := sum(variable[equinox.Mar:(solstice+60)]),
              !!varname.LO.SE    := sum(variable[equinox.Mar:DOY_off]),
              !!varname.SO.SE    := sum(variable[solstice:DOY_off]),
              !!varname.SOm30.SE := sum(variable[(solstice-30):DOY_off]),
              !!varname.SOp30.SE := sum(variable[(solstice+30):DOY_off]),
              !!varname.SOp60.SE := sum(variable[(solstice+60):DOY_off]),
              
              #solstice
              !!varname.solstice1 := sum(variable[(solstice-39):(solstice-10)]),
              !!varname.solstice2 := sum(variable[(solstice-29):solstice]),
              !!varname.solstice3 := sum(variable[(solstice-19):(solstice+10)]),
              !!varname.solstice4 := sum(variable[(solstice-9):(solstice+20)]),
              !!varname.solstice5 := sum(variable[(solstice+1):(solstice+30)]),
              !!varname.solstice6 := sum(variable[(solstice+11):(solstice+40)]) )
        }
      }
      
      
      #---------------------------------------------------------------------------------------------------------
      
      
      ####################################
      ## Calculate the monthly averages ##
      ####################################
      
      #create variable vectors
      VariableMeanVector = c("Tday","Tnight","Moist","SWrad")
      VariableSumVector  = c('GPPstart', 'Prcp')
      
      #get means and sums
      monthly_means = data.frame(daily_vals %>% 
                           group_by(Month) %>% 
                           summarize_at(VariableMeanVector, mean, na.rm = TRUE))
      monthly_sums = data.frame(daily_vals %>% 
                                   group_by(Month) %>% 
                                   summarize_at(VariableSumVector, sum, na.rm = TRUE))
    
      #merge
      monthly_vals = cbind(monthly_means,monthly_sums[,-c(1)])
      
      #Transform data
      monthly_vals = as.data.frame(t(monthly_vals))
      
      #Add to table
      #############
      
      #loop through variables
      for(i in 1:length(variable.names)) {
        #select  variable
        MONTHLY.DF = monthly_vals[variable.names[i],]
        #add column names
        names(MONTHLY.DF)=paste0(row.names(MONTHLY.DF), c(1:12))
        #cbind with table
        factors.sub = cbind(factors.sub, MONTHLY.DF)
      }
    
  
      #--------------------------------------------------------------------------
    
      
      ################################
      ## Get preseason temperatures ##
      ################################

      ## Calculate the average preseason temperatures prior to mean senescence date
      
      #get preseason length vector (10 to 120 days with 10-day steps)
      preseason.lengths = seq(10, 120, 10)
      
      #loop through preseasons
      for(preseason.length in preseason.lengths) {
        #name columns
        preseason.Tday   <- paste("Tday.PS",   preseason.length, sep=".")
        preseason.Tnight <- paste("Tnight.PS", preseason.length, sep=".")
        #add columns to table
        factors.sub = factors.sub %>%
          mutate(!!preseason.Tday   := mean(daily_vals$Tday[(DOY_off-preseason.length):DOY_off]),
                 !!preseason.Tnight := mean(daily_vals$Tnight[(DOY_off-preseason.length):DOY_off]) )
      }
      
    
      ##############################################################################################################################################
     
      
      # Safe the table
      write.table(factors.sub, file=paste0(Drivers_path, '/', timeseries_years, '.csv'), sep=',', row.names = F, col.names = T)
      
      }
}



##############################################################################################################################################
##############################################################################################################################################



##################
## Run the Loop ##
##################



#initialize the loop
outputlist <- pbmclapply(timeseries_year, parallelCalc, mc.cores=12, mc.preschedule=T)

#check how many files there are
length(list.files(path=Drivers_path, pattern='.csv'))
length(list.files(path=Drivers_path3, pattern='.csv'))

#Rbind files
climate.factors.table = rbindlist(lapply(list.files(path = Drivers_path), 
                                         function(n) fread(file.path(Drivers_path, n))))



##############################################################################################################################################
##############################################################################################################################################



###################
## Safe the data ##
###################



#Safe table
write.csv(climate.factors.table, paste(Drivers_path2, "Remote_sensing_drivers_data_onset_VNP.csv", sep="/"))

#Remove individual files
do.call(file.remove, list(list.files(Drivers_path, 
                                     full.names = TRUE)))
do.call(file.remove, list(list.files(Drivers_path3, 
                                     full.names = TRUE)))



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


