


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Climate driver extraction of the remote sensing analysis (EOS50) ##########################################
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
setwd("/Users/crowtherlabstation02/Desktop/Analysis")


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

# AET/PET
AET.PET_path     = "PEP_analysis/Analysis/Analysis_input/Drivers"

# GPP and LAI
GPP_path         = "Remote_sensing/Analysis/Analysis_input/Drivers/Modis_GPP_LAI"

#Phenology
Pheno_path       = "Remote_sensing/Analysis/Analysis_input/Phenology_data"


# 2. Output
###########

Drivers_path  = "Remote_sensing/Analysis/Analysis_input/Drivers_final_EOS50/Individual_files"
Drivers_path2 = "Remote_sensing/Analysis/Analysis_input/Drivers_final_EOS50/Merged_file"
Drivers_path3 = "Remote_sensing/Analysis/Analysis_input/Drivers_final_EOS50/Missing_observations"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



## Phenology data
#################

Pheno.df <- fread(paste(Pheno_path, "PhenologyData_Forests_025_North_Filtered_New.csv", sep="/")) %>%
  group_by(geometry) %>% 
  #delete pixels with less than 15 years
  filter(n() >= 15) %>%
  #get autumn phenology means per pixel
  mutate(Dormancy_DOY     = ifelse(Dormancy_DOY<170,365,Dormancy_DOY),#set dormancy that falls in next year to end of year
         MidGreendownMean = mean(MidGreendown_DOY),
         DormancyMax      = max(Dormancy_DOY)) %>%
  ungroup() %>%
  #delete duplicates
  distinct(geometry, Year, .keep_all = T) 
  
  
## AET/PET map
##############

#annual AET/PET ratio from SPLASH model
AET_PET.raster = raster(paste(AET.PET_path , "AET_PET_ratio_global_FULL_MODIS-C006_MOD15A2_v1.alpha_MEANANN.nc", sep="/"))


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



####################################################
## Add Soil texture, PFT, and AET/PET (meanalpha) ##
####################################################



#Add plant functional type info
Pheno.df$PFT <- "TBL" # T-BL-SG: Temperate broad-leaved summergreen tree
#PEP.df[PEP.df$species=='Larix',]$PFT <- "BNL" # B-NL-SG: Boreal needle-leaved summergreen tree

#Add biome information, AET-PET ratio and elevation (required for pmodel)
Pheno.df = 
  #both tables together
  cbind(Pheno.df, 
        # intersection
        data.frame(AET_PET = raster::extract(AET_PET.raster, Pheno.df[, c("Lon", "Lat")])),
        data.frame(alt     = raster::extract(elev.raster,    Pheno.df[, c("Lon", "Lat")])),
        data.frame(biome   = raster::extract(biome.raster,   Pheno.df[, c("Lon", "Lat")]))
        )

#remove stuff
rm(AET_PET.raster, elev.raster, biome.raster)



##############################################################################################################################################
##############################################################################################################################################



###############
## Constants ##
###############



## Constants in the Photosynthesis module
po2                 <- 20.9e3 #O2 partial pressure in Pa
p                   <- 1.0e5 # atmospheric pressure in Pa
bc3                 <- 0.015 # leaf respiration as fraction of Vmax for C3 plants
theta               <- 0.7 # colimitation (shape) parameter
q10ko               <- 1.2 #q10 for temperature-sensitive parameter ko
q10kc               <- 2.1 # q10 for temperature-sensitive parameter kc
q10tau              <- 0.57 # q10 for temperature-sensitive parameter tau
ko25                <- 3.0e4 # value of ko at 25 deg C
kc25                <- 30.0 # value of kc at 25 deg C
tau25               <- 2600.0 # value of tau at 25 deg C
alphaa              <- 0.5 # fraction of PAR assimilated at ecosystem level relative to leaf level
alphac3             <- 0.08 # intrinsic quantum efficiency of CO2 uptake in C3 plants
lambdamc3           <- 0.8 # optimal (maximum) lambda in C3 plants
cmass               <- 12.0107 # molecular mass of C [g mol-1]
cq                  <- 2.04e-6 # conversion factor for solar radiation from J m-2 to mol m-2
n0                  <- 7.15 # leaf N concentration (mg/g) not involved in photosynthesis
m                   <- 25.0 # corresponds to parameter p in Eqn 28, Haxeltine & Prentice 1996
t0c3                <- 250.0 # base temperature (K) in Arrhenius temperature response function for C3 plants
e0                  <- 308.56 # parameter in Arrhenius temp response function
tk25                <- 298.15 # 25 deg C in Kelvin
tmc3                <- 45.0 # maximum temperature for C3 photosynthesis
## Constants in the Water balance module
gamma               <- 65 # psychrometer constant gamma [Pa/K]
L                   <- 2.5*10^6 # latent heat of vaporization of water L [J/kg]
emissivity          <- 0.6 # emissivity for coniferous and deciduous surface type
k_sb                <- 5.670367*10^-8 # Stefan-Boltzman constant [W/m^2 K^4]
d1                  <- 0.5 # thickness of upper soil layer [m]
d2                  <- 1 # thickness of lower soil layer [m]    
a_m                 <- 1.391 # maximum Priestley-Taylor coefficient a_m 
g_m                 <- 3.26 # scaling conductance g_m [mm/s]
k_melt              <- 3 # rate of snowmelt [mm/???C d]

## Soil parameters depending on texture [Phenologoy_CO2_soil dataset]
E_max               <- 5 # maximum transpiration rate that can be sustained under well-watered conditions E_max [mm/d] --> depends on plant functional type (same for T-BD-SG and B-NL-SG)
# w_max = soil texture-dependent difference between field capacity and wilting point w_max [%]
# c_soil = soil texture-dependent maximum rate of ETA from the bare soil [mm/h]
# k_perc = soil texture-dependent conductivity cond_soil or percolation rate field capacity [mm/d]



##############################################################################################################################################
##############################################################################################################################################



######################
## Helper functions ##
######################


# Temperate inhibition function from LPJ-GUESS 
##############################################

temp_opt.fun <- function(temp) {
  x1        <- 1
  x2        <- 18
  x3        <- 25
  x4        <- 45
  k1        <- 2.*log((1/0.99)-1.)/(x1-x2)
  k2        <- (x1+x2)/2
  low       <- 1/(1+exp(k1*(k2-temp)))
  k3        <- log(0.99/0.01)/(x4-x3)
  high      <- 1-0.01*exp(k3*(temp-x3))
  tstress   <- low*high 
  if(tstress>=0) {
    tstress <- tstress
  } else {
    tstress <- 0
  }
  return(tstress)
}


# convert degC to kPa
#####################
degC_to_kPa.fun <- function(temp) {
  out       <- 0.6108*exp((17.27*temp)/(temp+237.3)) 
  return(out)
}


# Photoperiod function
######################

# photo = photoperiod 
# photo_min = minimum value during the growing season --> limited canopy development
# photo_max = maxmum value during the growing season --> allows canopies to develop unconstrained
photoperiod.fun <- function(photo, photo_min, photo_max) {
  if(photo<=photo_min) {
    photo_resp    <- 0
  }
  if(photo<photo_max & photo>photo_min) {
    photo_resp    <- (photo-photo_min)/(photo_max-photo_min)
  }
  if(photo>=photo_max) {
    photo_resp   <- 1
  }
  return(photo_resp)
}


# Vapour Pressure Deficit (VPD) function
########################################

# VPD = vapour pressure deficit [kPa]
# T_min & T_max = minimum and maximum daily temperature [C]
# VPD_min --> at low values, latent heat losses are unlikely to exceed available water 
# little effect on stomata
# VPD_max --> at high values, particularly if sustained, photosynthesis and growth are likely to be significantly limited
# complete stomatal closure

VPD.fun <- function(VPD, VPD_min, VPD_max) {
  if(VPD>=VPD_max) {
    y             <- 0
  }
  if(VPD<VPD_max & VPD>VPD_min) {
    y             <- 1-((VPD-VPD_min)/(VPD_max-VPD_min))
  }
  if(VPD<=VPD_min) {
    y             <- 1
  }
  return(y)
}


# Convert specific to relative humidity
#######################################

qair2rh <- function(qair, temp, press = 1013.25){
  es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
  e <- qair * press / (0.378 * qair + 0.622)
  rh <- e / es
  rh[rh > 1] <- 1
  rh[rh < 0] <- 0
  return(rh)
}



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
    
    #air humidity 
    QAIR      <- DataList[[5]][which(DataList[[5]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #soil moisture (<10cm)
    MOIST10   <- DataList[[6]][which(DataList[[6]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #soil moisture (10-40 cm)
    MOIST40   <- DataList[[7]][which(DataList[[7]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #net short-wave radiation
    SWRAD     <- DataList[[8]][which(DataList[[8]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #net long-wave radiation
    LWRAD     <- DataList[[9]][which(DataList[[9]]$site_year==pheno.sub$site_year),]%>% 
      dplyr::select(as.character(1:366))
    
    #short-wave radiation down
    SWRADdown <- DataList[[10]][which(DataList[[10]]$site_year==pheno.sub$site_year),]%>% 
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
      dplyr::select(geometry, Lat, Lon, alt, Year, Greenup_DOY, MidGreenup_DOY, MidGreendown_DOY) %>%
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
                             LWrad    = as.numeric(LWRAD),
                             SWradDown= as.numeric(SWRADdown),
                             Moist10  = as.numeric(MOIST10),
                             Moist40  = as.numeric(MOIST40),
                             Prcp     = as.numeric(PRCP), 
                             Qair     = as.numeric(QAIR),
                             Photo    = as.numeric(PHOTO),
                             GPP      = GPP,
                             LAI      = LAI)

    #Add climate variables and data wrangling
    daily_vals = daily_vals %>%
      filter(!is.na(Tmean)) %>%#delete NAs
      mutate(
        #add month and day identifiers
        Month = lubridate::month(as.Date(days,origin=days[1])),
        Day   = lubridate::day(as.Date(days,origin=days[1])),
        #relative humidity
        RH    = qair2rh(Qair, Tmean)*100,
        #dewpoint temperature
        Tdew  = weathermetrics::humidity.to.dewpoint(t = Tmean,
                                                     rh = RH,
                                                     temperature.metric = "celsius"))
      
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
      #growing-degree-days (>0??C)
      mutate(GDDday   = ifelse(Tday < 0 , 0, Tday),
             GDDnight = ifelse(Tnight < 0 , 0, Tnight))%>%
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
    DOY_off <- round(pheno.sub$MidGreendownMean)
    
    # Latest dormancy
    DOY_dorm <- round(pheno.sub$DormancyMax)
    
    # leaf-out
    DOY_out <- pheno.sub$Greenup_DOY
    
    # Greenup
    DOY_up  <- ifelse(pheno.sub$MidGreenup_DOY >= solstice, solstice-1, pheno.sub$MidGreenup_DOY)
    
    
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
    
    
    # Jmax limitation (photoperiod-dependency following Bauerle et al. 2012)
    #########################################################################
    
    #Number of days from solstice to mean senescence date
    PostSolsticeSpan = DOY_dorm-solstice
    
    #Spring degree-day threshold = 300
    GDD1 = daily_vals$GDDday #degree-day vector
    GDD1[1:DOY_out]=0 #set degree-days before leaf-out to zero
    GDD1 = GDD1[1:solstice] #set degree-days after solstice to zero
    GDDthreshold = ifelse(sum(GDD1) < 300, ifelse(sum(GDD1)>0, sum(GDD1), 1), 300) #set degree-day threshold
    GDD1 = cumsum(GDD1) #get cumulative degree-day vector
    GDD1[GDD1>GDDthreshold] <- GDDthreshold #cut of degree-day vector at threshold
    GDD1 = GDD1/GDDthreshold #bound between 0 and 1
    
    GDD2 = daily_vals$GDDday #degree-day vector
    GDD2[1:DOY_up]=0 #set degree-days before leaf-out to zero
    GDD2 = GDD2[1:solstice] #set degree-days after solstice to zero
    GDDthreshold = ifelse(sum(GDD2) < 300, ifelse(sum(GDD2)>0, sum(GDD2), 1), 300) #set degree-day threshold
    GDD2 = cumsum(GDD2) #get cumulative degree-day vector
    GDD2[GDD2>GDDthreshold] <- GDDthreshold #cut of degree-day vector at threshold
    GDD2 = GDD2/GDDthreshold #bound between 0 and 1
    
    #Daily Jmax vector (degree-day-based Jmax increase after leaf-out and linear Jmax decline after solstice)
    JmaxSAout = c(
      GDD1,
      rev(c(1:PostSolsticeSpan) / PostSolsticeSpan),
      rep(0, nrow(daily_vals)-DOY_dorm) )
    
    JmaxSAup = c(
      GDD2,
      rev(c(1:PostSolsticeSpan) / PostSolsticeSpan),
      rep(0, nrow(daily_vals)-DOY_dorm) )
    
    JmaxAout = c(
      rep(0,DOY_out),#Jmax=0 before leaf-out
      rep(1,solstice-DOY_out),#Jmax=1 after leaf-out
      rev(c(1:PostSolsticeSpan) / PostSolsticeSpan),
      rep(0, nrow(daily_vals)-DOY_dorm) )
    
    JmaxAup = c(
      rep(0,DOY_up),#Jmax=0 before leaf-out
      rep(1,solstice-DOY_up),#Jmax=1 after leaf-out
      rev(c(1:PostSolsticeSpan) / PostSolsticeSpan),
      rep(0, nrow(daily_vals)-DOY_dorm) )
    
    #Add Jmax vector to daily data
    daily_vals = daily_vals %>%
      mutate(JmaxA = JmaxSAout,
             JmaxB = JmaxSAup,
             JmaxC = JmaxAout,
             JmaxD = JmaxAup)
      
   
    
    ##############################################################################################################################################
    
    
    # Photosynthesis calculation
    ############################
    
    
    # GSI, Daily Net Photosynthesis rate (dA_n) and water stress factor (dw) are calculated daily 
    # and then accumulated by summation
    
    # Initialize vector to store daily values
    iGSI_year         <- vector()
    iGSIrad_year      <- vector()
    VPD_year          <- vector()
    iVPD_year         <- vector()  
    #dA_tot_year       <- vector()
    dA_totw_year      <- vector()
      
      # Loop through days of the growing season
      for(i in 1:nrow(daily_vals)) {

        ############################################
        ## Cumulative Growing Season Index (cGSI) ##
        ############################################
        
        # modified from Jolly et al. 2005
        
        # GSI...photoperiod-based growing-season index
        # GSI...irradiance-based growing-season index
        # VPD...vapor pressure deficit
        # iVPD...vapor pressure deficit function values
        
        # set VPD min and max
        #####################
        
        # Reference: White MA, Thornton PE, Running SW et al. (2000) Parameterization and sensitivity analysis of 
        # the BIOME???BGC terrestrial ecosystem model: net primary production controls. Earth Interactions, 4, 1???85.
        # mean of all evergreen needleleaf  tree species
        if(pheno.sub$biome %in% c(5,6,11,98,99)) {
          VPD_min       <- 0.61
          VPD_max       <- 3.1 
          } else {
          #mean of all deciduous broadleaf tree species
          VPD_min       <- 1.1
          VPD_max       <- 3.6 }
        
        # Estimate phoperiod thresholds based on the maximum and minimum values of the growing season
        photo_min <- min(daily_vals$Photo) 
        photo_max <- max(daily_vals$Photo) 

        # e_s: saturation vapour pressure [kPa] 
        e_s <- (degC_to_kPa.fun(temp=daily_vals$Tmax[i]) + degC_to_kPa.fun(temp=daily_vals$Tmin[i])) / 2
        
        # e_a: derived from dewpoint temperature [kPa]
        e_a <- degC_to_kPa.fun(temp=daily_vals$Tdew[i])
        
        # VPD: Vapour pressure deficit [kPa]
        VPD  <- e_s-e_a
        VPD_year <- c(VPD_year,VPD)

        # apply vapor pressure deficit funtion
        iVPD <- VPD.fun(VPD, VPD_min, VPD_max)
        iVPD_year <- c(iVPD_year, iVPD)
        
        # iOpt_temp: response to optimal temperature (Gompertz function)
        iOpt <- temp_opt.fun(daily_vals$Tday[i])
        
        # iPhoto: photoperiod response
        iPhoto <- photoperiod.fun(daily_vals$Photo[i], photo_min, photo_max)
        
        # iRadiation
        # get maximum radiation at the site (field capacity)
        max.rad = max(DataList[[8]][which(DataList[[8]]$geometry==pheno.sub$geometry), c(as.character(1:365))], na.rm=T)
        iRad <- daily_vals$SWrad[i] / max.rad
        
        # Calculate daily GSI
        iGSI    <- as.numeric(iVPD*iOpt*iPhoto)
        iGSIrad <- as.numeric(iVPD*iOpt*iRad)
        
        # Add to the cumulative cGSI
        iGSI_year    <- c(iGSI_year,iGSI)
        iGSIrad_year <- c(iGSIrad_year,iGSIrad)

        #----------------------------------------------------------------------------------------------
        
        ############################
        ## Zani et al. 2020 model ##
        ############################
        
        # Net photosynthesis rate (PHOTOSYNTHESIS-CONDUCTANCE MODEL, ref. Sitch et al. 2003)
        
        # apar: daily integral of absorbed photosynthetically active radiation (PAR), J m-2 d-1
        # Eqn 4, Haxeltine & Prentice 1996
        # alphaa: scaling factor for absorbed PAR at ecosystem, versus leaf, scale
        # nearly half of short-wave radiation is PAR --> mean annual value of 0.473 observed for the irradiance ratio 
        # in the PAR (ref. Papaioannou et al. 1993) plus 8% reflected and transmitted
        # convert in J/m^-2 day: the power in watts (W) is equal to the energy in joules (J), divided by the time period in seconds (s): 
        # --> 1 Watt = 1 Joule/second, therefore j = W*86400
        apar <- alphaa * daily_vals$SWrad[i] * 60 * 60 * 24
        
        # Calculate temperature inhibition function limiting photosynthesis at low and high temperatures (ref. Sitch et al. 2002)
        tstress <- temp_opt.fun(daily_vals$Tday[i])
        
        # Calculate catalytic capacity of rubisco, Vm, assuming optimal (non-water-stressed) value for lambda, i.e. lambdamc3
        # adjust kinetic parameters for their dependency on temperature 
        # i.e. relative change in the parameter for a 10 degC change in temperature
        # Eqn 22, Haxeltine & Prentice 1996a
        
        ko  <- ko25*q10ko**((daily_vals$Tday[i]-25.0)/10.0)   # Michaelis constant of rubisco for O2
        kc  <- kc25*q10kc**((daily_vals$Tday[i]-25.0)/10.0)   # Michaelis constant for CO2
        tau <- tau25*q10tau**((daily_vals$Tday[i]-25.0)/10.0)# CO2/O2 specificity ratio
        
        # gammastar: CO_2 compensation point [CO2 partial pressure, Pa]   
        # Eqn 8, Haxeltine & Prentice 1996
        gammastar <- po2/(2.0*tau)
        
        # Convert ambient CO2 level from mole fraction to partial pressure, Pa
        pa <- CO2*p
        
        # p_i: non-water-stressed intercellular CO2 partial pressure, Pa
        # Eqn 7, Haxeltine & Prentice 1996
        p_i <- pa*lambdamc3 
        
        # Calculate coefficients
        # Eqn 4, Haxeltine & Prentice 1996
        c1 <- tstress*alphac3*((p_i-gammastar)/(p_i+2.0*gammastar))
        
        # Eqn 6, Haxeltine & Prentice 1996
        c2 <- (p_i-gammastar)/(p_i+kc*(1.0+po2/ko)) 
        b <- bc3 # choose C3 value of b for Eqn 10, Haxeltine & Prentice 1996
        t0 <- t0c3 # base temperature for temperature response of rubisco
        
        # Eqn 13, Haxeltine & Prentice 1996
        s <- (24.0 / daily_vals$Photo[i] ) * b
        
        # Eqn 12, Haxeltine & Prentice 1996
        sigma <- sqrt(max(0.0,1.0-(c2-s)/(c2-theta*s)))
        
        # vm: optimal rubisco capacity, gC m-2 d-1 
        # Eqn 11, Haxeltine & Prentice 1996
        # cmass: the atomic weight of carbon, used in unit conversion from molC to g 
        # cq: conversion factor from apar [J m-2] to photosynthetic photon flux density [mol m-2]
        vm <- (1.0/b)*(c1/c2)*((2.0*theta-1.0)*s-(2.0*theta*s-c2)*sigma)*apar*cmass*cq
        
        # je: PAR-limited photosynthesis rate, gC m-2 h-1
        # Eqn 3, Haxeltine & Prentice 1996
        # Convert je from daytime to hourly basis
        if(daily_vals$Photo[i]==0) {
          je <- 0
        } else {
          je <- c1*apar*cmass*cq / daily_vals$Photo[i]
        }
        
        # jc: rubisco-activity-limited photosynthesis rate, gC m-2 h-1
        # Eqn 5, Haxeltine & Prentice 1996
        jc <- c2*vm/24.0
        
        # agd: daily gross photosynthesis, gC m-2 d-1
        # Eqn 2, modified with k_shape (theta)
        if(je<1e-10 | jc<=1e-10) {
          agd <- 0
        } else {
          agd <- (je+jc-sqrt((je+jc)**2.0-4.0*theta*je*jc))/(2.0*theta) * daily_vals$Photo[i]
        }
        
        # rd: daily leaf respiration, gC m-2 d-1
        # Eqn 10, Haxeltine & Prentice 1996
        rd <- b*vm
        
        # and: daily net photosynthesis (at leaf level), gC m-2 d-1
        and <- agd-rd
        
        # adt: total daytime net photosynthesis, gC m-2 d-1
        # Eqn 19, Haxeltine & Prentice 1996
        adt <- and + (1.0 - daily_vals$Photo[i] / 24.0) * rd
        
        # Convert adt from gC m-2 d-1 to mm m-2 d-1 using ideal gas equation
        #adtmm <- adt / cmass * 8.314 * (daily_vals$TMEAN[i] + 273.3) / p * 1000.0
        
        # Store the daily result in the yearly vector
        #dA_tot_year <- c(dA_tot_year,adt) #daytime net photosynthesis
        
        
        ## Water Stress Factor (ref. Gerten et al. 2004)
        ################################################
        
        # soil is treated as a simple bucket consisting of two layers with fixed thickness
        
        # Calculate potential evapotranspiration (ETA) rate, E_pot, mm d-1
        
        # delta: rate of increase of the saturation vapour pressure with temperature
        delta <- (2.503*10^6 * exp((17.269 * daily_vals$Tday[i]) / (237.3 + daily_vals$Tday[i]))) / (237.3 + daily_vals$Tday[i])^2
        
        # R_n: istantaneous net radiation, W m-2 = R_s net short-wave radiation flux + R_l net long-wave flux 
        R_n <- daily_vals$SWrad[i] + daily_vals$LWrad[i]
        
        # E_eq: equilibrium EvapoTranspiration
        # from seconds to day
        E_eq <- 24 * 3600 * (delta / (delta + gamma)) * (R_n / L) 
        
        # E_pot: potential EvapoTranspiration = equilibrium ETA * Priestley-Taylor coefficient 
        E_pot <- E_eq*a_m
        
        # ratio: stomata-controlled ratio between intercellular and ambient CO2 partial pressure in the absence of water limitation
        ratio <- p_i/pa # ca. 0.8
        
        # g_min: minimum canopy conductance, mm s-1
        # depends on PFT (broadleaf = 0.5, needleleaf = 0.3)
        if(pheno.sub$biome %in% c(1,4,8,9,10,12,13)) { 
          g_min <- 0.5*3600*24 # from seconds to day
        } else {
          g_min <- 0.3*3600*24
        }
        
        # g_pot: nonwater-stressed potential canopy conductance, mm s-1
        g_pot <- g_min + ((1.6*adt)/((pa/p)*(1-ratio)))
        
        # E_demand: atmoshperic demand 
        # unstressed transpiration which occurs when stomatal opening is not limited by reduced water potential in the plant
        E_demand <- E_pot/(1+(g_m/g_pot))
        
        # root1/2: fraction of roots present in the respective layers
        # depends on PFT (temperate = 0.7/0.3, boreal = 0.9/0.1)
        if (pheno.sub$biome %in% c(1,3,4,8,9,10,12,13)) { 
          root1 <- 0.7  
          root2 <- 0.3
        } else {
          root1 <- 0.9
          root2 <- 0.1
        }
        
        # relative soil moisture wr:
        # ratio between current soil water content and plant-available water capacity
        # wr ratio is computed for both soil layers by 
        # weighting their relative soil water contents (w1, w2) 
        # with the fraction of roots present in the respective layer
        w1  <- daily_vals$Moist10[i]
        w2  <- daily_vals$Moist40[i]
        
        # soil texture-dependent difference between field capacity and wilting point w_max [%]
        w_max <- 15
        wr <- root1*(w1/w_max) + root2*(w2/w_max)
        
        # E_supply: plant- and soil-limited supply function 
        E_supply <- as.numeric(E_max*wr)
        
        # dw: daily water stress factor
        dw <- min(1,(E_supply/E_demand))
        
        # dA_totw: daily net photosynthesis modified by water stress factor
        dA_totw <- adt*dw
        
        # Add daily result to the yearly vector
        dA_totw_year <- c(dA_totw_year, dA_totw)
        
      } # END loop through days of the growing season
      
      #set values before leaf-out to zero
      iGSI_year[1:DOY_out]    = 0
      iGSIrad_year[1:DOY_out] = 0
      dA_totw_year[1:DOY_out] = 0
      
      #set negative values to zero
      VPD_year[VPD_year<=0]        = 0.001
      dA_totw_year[dA_totw_year<0] = 0
      
      #Jmax corrected photosynthesis
      dA_totw.JmaxA_year <- dA_totw_year * daily_vals$JmaxA
      dA_totw.JmaxB_year <- dA_totw_year * daily_vals$JmaxB
      dA_totw.JmaxC_year <- dA_totw_year * daily_vals$JmaxC
      dA_totw.JmaxD_year <- dA_totw_year * daily_vals$JmaxD
      
      #add VPD to daily table
      daily_vals$VPD               = VPD_year *1000 #VPD in Pa
         
      #----------------------------------------------------------------------------------------------
      
      ##################
      ## P-model v1.0 ##
      ##################
      
      ## Benjamin D. Stocker et al. 2020
      ## optimality-based light use efficiency model 
      ## for simulating ecosystem gross primary production
      
      # constant variables
      alt = as.numeric(pheno.sub$alt) # elevation z [m a.s.l.]
      meana = pheno.sub$AET_PET # Local annual mean ratio of actual over potential evapotranspiration
      
      ## Calculate Photosynthetic Photon Flux Density, ppfd [mol m-2]
      # PAR as irradiance [W m-2] is given by incoming short-wave radiation
      ppfd = 60 * 60 * 24 * 10^-6 * 2.04 * (daily_vals$SWradDown)
      
      ## get maximum soil moisture at the site (field capacity)
      field.capacity = max(DataList[[7]][which(DataList[[7]]$geometry==pheno.sub$geometry), c(as.character(1:365))], na.rm=T)
      
      ## P-model v1.0
      pmodel.df <- tibble(
        tc             = daily_vals$Tday,
        vpd            = daily_vals$VPD, #VPD in Pa 
        co2            = CO2,
        fapar          = 1,
        ppfd           = ppfd,
        soilm          = daily_vals$Moist40 / field.capacity
      ) %>%
        mutate(out_pmodel = purrr::pmap(., rpmodel,
                                        elv            = alt,
                                        kphio          = 0.087,
                                        beta           = 146,
                                        method_optci   = "prentice14",
                                        method_jmaxlim = "wang17",
                                        do_ftemp_kphio = T,
                                        do_soilmstress = T,
                                        meanalpha=meana
        ))
      pmodel.df = do.call(rbind.data.frame, pmodel.df$out_pmodel)
      
      #set Photosynthesis before leaf-out to zero
      pmodel.df[1:DOY_out,]=0
      
      ## Dark respiration, rd [mol C m-2]
      rd = pmodel.df$rd
      rd = rd * cmass # convert (carbon mass)
      
      #get daily values of net daytime photosynthesis [g C m-2]
      Apm = (pmodel.df$gpp - rd) + (1.0-daily_vals$Photo/24.0)*rd
      
      #set negative values to zero
      Apm[Apm<0] = 0
      
      # Jmax-limited photosynthesis
      ApmJmaxA = Apm * daily_vals$JmaxA
      ApmJmaxB = Apm * daily_vals$JmaxB 
      ApmJmaxC = Apm * daily_vals$JmaxC 
      ApmJmaxD = Apm * daily_vals$JmaxD 
      
      #----------------------------------------------------------------------------------------------
      
      #Store the results
      ##################
      
      daily_vals = daily_vals %>%
        mutate(GSI        = iGSI_year,          #photoperiod-influenced GSI
               GSIrad     = iGSIrad_year,       #radiation-influenced GSI
               Azani      = dA_totw_year,       #net daytime photosynthesis (Zani et al., water-stressed)
               AzaniJmaxA = dA_totw.JmaxA_year, #net daytime photosynthesis spring and autumn Jmax-limited (Zani et al., water-stressed)
               AzaniJmaxB = dA_totw.JmaxB_year, #net daytime photosynthesis autumn Jmax-limited (Zani et al., water-stressed)
               AzaniJmaxC = dA_totw.JmaxC_year,
               AzaniJmaxD = dA_totw.JmaxD_year,
               Apm        = Apm,                #net daytime photosynthesis (p model)
               ApmJmaxA   = ApmJmaxA,           #net daytime photosynthesis Jmax-limited (p model)
               ApmJmaxB   = ApmJmaxB,
               ApmJmaxC   = ApmJmaxC,          
               ApmJmaxD   = ApmJmaxD
        )  %>%
        rename(Moist=Moist40) 
    
      
      ##############################################################################################################################################
      
     
      ###################
      ## Store drivers ##
      ###################
      
      
      ######################
      ## Seasonal drivers ##
      ######################
      
      #define variables
      variable.names = c('Azani', 'AzaniJmaxA', 'AzaniJmaxB', 'AzaniJmaxC', 'AzaniJmaxD', 
                         'Apm', 'ApmJmaxA', 'ApmJmaxB', 'ApmJmaxC', 'ApmJmaxD',
                         'GPP', 'LAI','GPPstart','LAIstart',
                         'GSI', 'GSIrad', 
                         'GDDday', 'GDDnight', 'SWrad',
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
        
        if(variable.names[i] %in% c('Azani', 'AzaniJmaxA', 'AzaniJmaxB', 'AzaniJmaxC', 'AzaniJmaxD', 
                                    'Apm', 'ApmJmaxA', 'ApmJmaxB', 'ApmJmaxC', 'ApmJmaxD', 
                                    'GPP', 'LAI','GPPstart','LAIstart',
                                    'GSI', 'GSIrad', 'GDDday', 'GDDnight')){
          
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
      VariableMeanVector = c('LAI','LAIstart',"Tday","Tnight","Moist","SWrad")
      VariableSumVector  = c('Azani', 'AzaniJmaxA', 'AzaniJmaxB', 'AzaniJmaxC', 'AzaniJmaxD', 
                             'Apm', 'ApmJmaxA', 'ApmJmaxB', 'ApmJmaxC', 'ApmJmaxD',   
                             'GPP','GPPstart',
                             'GSI', 'GSIrad', 'GDDday', 'GDDnight', 'Prcp')
      
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
outputlist <- pbmclapply(timeseries_year, parallelCalc, mc.cores=5, mc.preschedule=T)

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
write.csv(climate.factors.table, paste(Drivers_path2, "Remote_sensing_drivers_data.csv", sep="/"))

#Remove individual files
do.call(file.remove, list(list.files(Drivers_path, 
                                     full.names = TRUE)))
do.call(file.remove, list(list.files(Drivers_path3, 
                                     full.names = TRUE)))



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


