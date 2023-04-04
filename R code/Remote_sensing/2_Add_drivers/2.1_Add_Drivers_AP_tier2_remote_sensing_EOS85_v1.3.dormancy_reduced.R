


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Climate driver extraction of the remote sensing analysis (EOS85) ##########################################
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

Drivers_path  = "Remote_sensing/Analysis/Analysis_input/Drivers_final_EOS85/Individual_files"
Drivers_path2 = "Remote_sensing/Analysis/Analysis_input/Drivers_final_EOS85/Merged_file"
Drivers_path3 = "Remote_sensing/Analysis/Analysis_input/Drivers_final_EOS85/Missing_observations"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



## Phenology data
#################

Pheno.df <- fread(paste(Pheno_path, "PhenologyData_Forests_025_North_Filtered_New.csv", sep="/")) %>%
  filter(Dormancy_DOY>230,
         #delete observation if EOS85 date occurs before EOS50 or EOS10 dates
         Dormancy_DOY>MidGreendown_DOY,
         Dormancy_DOY>Senesc_DOY) %>% 
  group_by(geometry) %>% 
  #get autumn phenology means per pixel
  mutate(DormancyMean          = mean(Dormancy_DOY),
         MidGreendownMean      = mean(MidGreendown_DOY),
         SenescMean            = mean(Senesc_DOY),
         meanDuration_EOS10_50 = round(MidGreendownMean-SenescMean),
         meanDuration_EOS10_85 = round(DormancyMean-SenescMean),
         meanDuration_EOS50_85 = round(DormancyMean-MidGreendownMean)) %>%
  #delete pixels with less than 15 years
  filter(n() >= 15) %>%
  ungroup() %>%
  #delete duplicates
  distinct(geometry, Year, .keep_all = T) 


## Elevation map
################

elev.raster = raster(paste(photo_path, "topo_elevation.asc", sep="/"))


## Biomes map
#############

biome.raster = raster(paste(photo_path, "WWF_Biomes_HalfDegree.tif", sep="/"))


## CO2 data
###########

CO2.df = fread(paste(CO2_path, "CO2_Annual.csv", sep="/"))


## Photoperiod
##############

photo.df = fread(paste(photo_path, "Photoperiod.csv", sep="/"))


## LAI and GPP
##############

GPP.df = fread(paste(GPP_path, "GppData_Forests_025.csv", sep='/'))%>%
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

LAI.df = fread(paste(GPP_path, "LaiData_Forests_025.csv", sep='/'))%>%
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
colnames(LAI.df)[6:ncol(LAI.df)] <- seq(1, 366, by=8)


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
        'GLDAS_Daily_Data_Swnet_tavg',
        'GLDAS_Daily_Data_Lwnet_tavg',
        'GLDAS_Daily_Data_SWdown_f_tavg')

#create empty list
DataList <- replicate(length(vn),data.frame())

#loop through climate variables
for(i in 1:length(vn)) {
  #read data
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
  
  #add table to list
  DataList[[i]] <- data
}
#add names to list
names(DataList)=vn
# Note: Precipitation is given as rate in mm d-1.



##############################################################################################################################################
##############################################################################################################################################



#########################################################
## Add biome information, AET-PET ratio and elevation  ##
#########################################################



Pheno.df = 
  #cbind tables
  cbind(Pheno.df, 
        # intersection
        data.frame(alt     = raster::extract(elev.raster,    Pheno.df[, c("Lon", "Lat")])),
        data.frame(biome   = raster::extract(biome.raster,   Pheno.df[, c("Lon", "Lat")]))
        )

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
DataList[[11]] = photo.df
DataList[[12]] = CO2.df
DataList[[13]] = Pheno.df
DataList[[14]] = GPP.df
DataList[[15]] = LAI.df

rm(photo.df, CO2.df, data, Pheno.df, GPP.df, LAI.df)
names(DataList)=c(vn,"photoperiod",'CO2',"Pheno","GPP","LAI")
names(DataList)
#[1] "GLDAS_Daily_Data_Tair_f_inst_Mean"    "GLDAS_Daily_Data_Tair_f_inst_Min"     "GLDAS_Daily_Data_Tair_f_inst_Max"     "GLDAS_Daily_Data_Rainf_f_tavg"       
#[5] "GLDAS_Daily_Data_Qair_f_inst"         "GLDAS_Daily_Data_SoilMoi0_10cm_inst"  "GLDAS_Daily_Data_SoilMoi10_40cm_inst" "GLDAS_Daily_Data_Swnet_tavg"         
#[9] "GLDAS_Daily_Data_Lwnet_tavg"          "GLDAS_Daily_Data_SWdown_f_tavg"       "photoperiod"                          "CO2"                                 
#[13] "Pheno"                                "GPP"                                  "LAI"           


##############################################################################################################################################


################################
# Loop through all time-points #
################################


parallelCalc <- function(timeseries_years){ 

  # Subset input data by time-point
  #################################

  #phenology data
  pheno.sub  <- DataList[[13]][which(DataList[[13]]$site_year==timeseries_years),]
  
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
    MOIST     <- DataList[[7]][which(DataList[[7]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #net short-wave radiation
    SWRAD     <- DataList[[8]][which(DataList[[8]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #day length
    PHOTO     <- DataList[[11]][which(DataList[[11]]$geometry==pheno.sub$geometry),][1]%>% 
      dplyr::select(as.character(1:366))
    
    #CO2 (monthly)
    CO2       <- DataList[[12]][which(DataList[[12]]$Year==pheno.sub$Year),]$CO2
    
    #GPP
    GPP <- as.numeric(DataList[[14]][which(DataList[[14]]$site_year==pheno.sub$site_year),] %>%
      dplyr::select(6:ncol(DataList[[14]])) )
    GPP = rep(GPP, each=8) / 8
    GPP = GPP[1:366]
    
    #LAI
    LAI <- as.numeric(DataList[[15]][which(DataList[[15]]$site_year==pheno.sub$site_year),] %>%
                        dplyr::select(6:ncol(DataList[[15]])) )
    LAI = rep(LAI, each=8)
    LAI = LAI[1:366]
    
    
    ##############################################################################################################################################
    
    
    # Create table of daily climate
    ###############################
    
    # Generate sub-dataframe to store results
    factors.sub <- pheno.sub %>% 
      dplyr::select(geometry, Lat, Lon, alt, Year, Greenup_DOY, MidGreenup_DOY, Senesc_DOY, MidGreendown_DOY, Dormancy_DOY) %>%
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
                             GPP      = GPP,
                             LAI      = LAI)

    #Add climate variables and data wrangling
    daily_vals = daily_vals %>%
      filter(!is.na(Tmean)) %>%#delete NAs
      mutate(
        #add month and day identifiers
        Month = lubridate::month(as.Date(days,origin=days[1])),
        Day   = lubridate::day(as.Date(days,origin=days[1])))
      
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
    
    # Mean EOS85
    Mean_EOS85 <- round(pheno.sub$DormancyMean)
    
    # Mean EOS50
    Mean_EOS50 <- round(pheno.sub$MidGreendownMean)
    
    # Mean EOS10
    Mean_EOS10 <- round(pheno.sub$SenescMean)
    
    # EOS10
    DOY_EOS10 <- pheno.sub$Senesc_DOY
    
    # EOS50
    DOY_EOS50 <- pheno.sub$MidGreendown_DOY
    
    # leaf-out
    DOY_out <- pheno.sub$Greenup_DOY
    
    
    ##############################################################################################################################################
    
    
    # Set GPP and LAI before greenup to zero
    ########################################
    
    
    #GPP
    daily_vals$GPPstart = daily_vals$GPP
    daily_vals$GPPstart[1:DOY_out]=0
    
    #LAI
    daily_vals$LAIstart = daily_vals$LAI
    daily_vals$LAIstart[1:DOY_out]=0
    
      
      ##############################################################################################################################################
      
     
      ###################
      ## Store drivers ##
      ###################
      
      
      ######################
      ## Seasonal drivers ##
      ######################
      
      #define variables
      variable.names = c('GPPstart',
                         'SWrad',
                         'Tday', 'Tnight', 'Moist',
                         'Prcp')
      
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
        
        if(variable.names[i] %in% c('GPPstart')){
          
          # Sums from leaf-out
          ####################
          
          factors.sub = factors.sub %>%
            mutate(
              #seasonal
              !!varname.LO.SO    := sum(variable[DOY_out:solstice]),
              !!varname.LO.SOm30 := sum(variable[DOY_out:(solstice-30)]),
              !!varname.LO.SOp30 := sum(variable[DOY_out:(solstice+30)]),
              !!varname.LO.SOp60 := sum(variable[DOY_out:(solstice+60)]),
              !!varname.LO.SE    := sum(variable[DOY_out:Mean_EOS85]),
              !!varname.SO.SE    := sum(variable[solstice:Mean_EOS85]),
              !!varname.SOm30.SE := sum(variable[(solstice-30):Mean_EOS85]),
              !!varname.SOp30.SE := sum(variable[(solstice+30):Mean_EOS85]),
              !!varname.SOp60.SE := sum(variable[(solstice+60):Mean_EOS85]),
              
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
              !!varname.LO.SE    := mean(variable[equinox.Mar:Mean_EOS85]),
              !!varname.SO.SE    := mean(variable[solstice:Mean_EOS85]),
              !!varname.SOm30.SE := mean(variable[(solstice-30):Mean_EOS85]),
              !!varname.SOp30.SE := mean(variable[(solstice+30):Mean_EOS85]),
              !!varname.SOp60.SE := mean(variable[(solstice+60):Mean_EOS85]),
              
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
              !!varname.LO.SE    := sum(variable[equinox.Mar:Mean_EOS85]),
              !!varname.SO.SE    := sum(variable[solstice:Mean_EOS85]),
              !!varname.SOm30.SE := sum(variable[(solstice-30):Mean_EOS85]),
              !!varname.SOp30.SE := sum(variable[(solstice+30):Mean_EOS85]),
              !!varname.SOp60.SE := sum(variable[(solstice+60):Mean_EOS85]),
              
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
      VariableSumVector  = c('GPPstart','Prcp')
      
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

      
      ## Calculate the average preseason temperatures prior to mean EOS85 date
      ########################################################################
      
      #get preseason length vector (10 to 120 days with 10-day steps)
      preseason.lengths = seq(10, 120, 10)
      
      #loop through preseasons
      for(preseason.length in preseason.lengths) {
        #name columns
        preseason.Tday   <- paste("Tday.PS",   preseason.length, sep=".")
        preseason.Tnight <- paste("Tnight.PS", preseason.length, sep=".")
        #add columns to table
        factors.sub = factors.sub %>%
          mutate(!!preseason.Tday   := mean(daily_vals$Tday[(Mean_EOS85-preseason.length):Mean_EOS85], na.rm=T),
                 !!preseason.Tnight := mean(daily_vals$Tnight[(Mean_EOS85-preseason.length):Mean_EOS85], na.rm=T) )
      }
      
      #--------------------------------------------------------------------------
      
      ## Calculate the average temperature after mean EOS10 date
      ##########################################################
      
      #get preseason length vector (30 to 120 days with 30-day steps)
      preseason.lengths = seq(30, 120, 30)
      
      #loop through preseasons
      for(preseason.length in preseason.lengths) {
        #name columns
        preseason.Tday   <- paste("Tday.EOS10mean.PS",   preseason.length, sep=".")
        preseason.Tnight <- paste("Tnight.EOS10mean.PS", preseason.length, sep=".")
        #add columns to table
        factors.sub = factors.sub %>%
          mutate(!!preseason.Tday   := mean(daily_vals$Tday[Mean_EOS10:(Mean_EOS10+preseason.length)], na.rm=T),
                 !!preseason.Tnight := mean(daily_vals$Tnight[Mean_EOS10:(Mean_EOS10+preseason.length)], na.rm=T) )
      }
      
      
      ## Calculate the average temperature after EOS10 date
      #####################################################
      
      #get preseason length vector (30 to 120 days with 30-day steps)
      preseason.lengths = seq(30, 120, 30)
      
      #loop through preseasons
      for(preseason.length in preseason.lengths) {
        #name columns
        preseason.Tday   <- paste("Tday.EOS10.PS",   preseason.length, sep=".")
        preseason.Tnight <- paste("Tnight.EOS10.PS", preseason.length, sep=".")
        #add columns to table
        factors.sub = factors.sub %>%
          mutate(!!preseason.Tday   := mean(daily_vals$Tday[DOY_EOS10:(DOY_EOS10+preseason.length)], na.rm=T),
                 !!preseason.Tnight := mean(daily_vals$Tnight[DOY_EOS10:(DOY_EOS10+preseason.length)], na.rm=T) )
      }
      
      #--------------------------------------------------------------------------
      
      ## Calculate the average temperature after mean EOS50 date
      ##########################################################
      
      #get preseason length vector (30 to 120 days with 30-day steps)
      preseason.lengths = seq(30, 120, 30)
      
      #loop through preseasons
      for(preseason.length in preseason.lengths) {
        #name columns
        preseason.Tday   <- paste("Tday.EOS50mean.PS",   preseason.length, sep=".")
        preseason.Tnight <- paste("Tnight.EOS50mean.PS", preseason.length, sep=".")
        #add columns to table
        factors.sub = factors.sub %>%
          mutate(!!preseason.Tday   := mean(daily_vals$Tday[Mean_EOS50:(Mean_EOS50+preseason.length)], na.rm=T),
                 !!preseason.Tnight := mean(daily_vals$Tnight[Mean_EOS50:(Mean_EOS50+preseason.length)], na.rm=T) )
      }
      
      
      ## Calculate the average temperature after EOS50 date
      #####################################################
      
      #get preseason length vector (30 to 120 days with 30-day steps)
      preseason.lengths = seq(30, 120, 30)
      
      #loop through preseasons
      for(preseason.length in preseason.lengths) {
        #name columns
        preseason.Tday   <- paste("Tday.EOS50.PS",   preseason.length, sep=".")
        preseason.Tnight <- paste("Tnight.EOS50.PS", preseason.length, sep=".")
        #add columns to table
        factors.sub = factors.sub %>%
          mutate(!!preseason.Tday   := mean(daily_vals$Tday[DOY_EOS50:(DOY_EOS50+preseason.length)], na.rm=T),
                 !!preseason.Tnight := mean(daily_vals$Tnight[DOY_EOS50:(DOY_EOS50+preseason.length)], na.rm=T) )
      }
      
      #--------------------------------------------------------------------------
      
      ## Calculate the average temperatures starting at EOS10/50 date with a length of the average duration
      #####################################################################################################
      
      factors.sub = factors.sub %>%
        mutate(Tday.EOS10_50   = mean(daily_vals$Tday[DOY_EOS10:(DOY_EOS10+pheno.sub$meanDuration_EOS10_50)], na.rm=T),
               Tnight.EOS10_50 = mean(daily_vals$Tnight[DOY_EOS10:(DOY_EOS10+pheno.sub$meanDuration_EOS10_50)], na.rm=T),
               
               Tday.EOS10_85   = mean(daily_vals$Tday[DOY_EOS10:(DOY_EOS10+pheno.sub$meanDuration_EOS10_85)], na.rm=T),
               Tnight.EOS10_85 = mean(daily_vals$Tnight[DOY_EOS10:(DOY_EOS10+pheno.sub$meanDuration_EOS10_85)], na.rm=T),
               
               Tday.EOS50_85   = mean(daily_vals$Tday[DOY_EOS50:(DOY_EOS50+pheno.sub$meanDuration_EOS50_85)], na.rm=T),
               Tnight.EOS50_85 = mean(daily_vals$Tnight[DOY_EOS50:(DOY_EOS50+pheno.sub$meanDuration_EOS50_85)], na.rm=T)
               )
      
      
      ## Calculate the average temperatures starting at mean EOS10/50 date with a length of the average duration
      ##########################################################################################################
      
      factors.sub = factors.sub %>%
        mutate(Tday.EOS10mean_50   = mean(daily_vals$Tday[Mean_EOS10:(Mean_EOS10+pheno.sub$meanDuration_EOS10_50)], na.rm=T),
               Tnight.EOS10mean_50 = mean(daily_vals$Tnight[Mean_EOS10:(Mean_EOS10+pheno.sub$meanDuration_EOS10_50)], na.rm=T),
               
               Tday.EOS10mean_85   = mean(daily_vals$Tday[Mean_EOS10:(Mean_EOS10+pheno.sub$meanDuration_EOS10_85)], na.rm=T),
               Tnight.EOS10mean_85 = mean(daily_vals$Tnight[Mean_EOS10:(Mean_EOS10+pheno.sub$meanDuration_EOS10_85)], na.rm=T),
               
               Tday.EOS50mean_85   = mean(daily_vals$Tday[Mean_EOS50:(Mean_EOS50+pheno.sub$meanDuration_EOS50_85)], na.rm=T),
               Tnight.EOS50mean_85 = mean(daily_vals$Tnight[Mean_EOS50:(Mean_EOS50+pheno.sub$meanDuration_EOS50_85)], na.rm=T)
        )
      
      
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
                                         function(n) fread(file.path(Drivers_path, n))), fill=T)



##############################################################################################################################################
##############################################################################################################################################



###################
## Safe the data ##
###################



#Safe table
write.csv(climate.factors.table, paste(Drivers_path2, "Remote_sensing_drivers_data_EOS85.csv", sep="/"))

#Remove individual files
do.call(file.remove, list(list.files(Drivers_path, 
                                     full.names = TRUE)))
do.call(file.remove, list(list.files(Drivers_path3, 
                                     full.names = TRUE)))



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


