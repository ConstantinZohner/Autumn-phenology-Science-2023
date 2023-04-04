##############################################################################################################################################
############################################################# Rscript for: ###################################################################
##############################################################################################################################################
#### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ####################################
##############################################################################################################################################
# This script extracts the climate drivers for the PEP725 data ###############################################################################
##############################################################################################################################################



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



##############################################################################################################################################
##############################################################################################################################################



##########################################
## Set directories and get own PEP data ##
##########################################



# Set the working directory
setwd("/Users/crowtherlabstation02/Desktop/Analysis/PEP_analysis/Analysis")


#########
# Paths #
#########


# 1. Input
##########

#Climate
GLDAS_path       = "Analysis_input/Drivers/GLDAS"

#CO2
CO2_path         = "Analysis_input/Drivers/CO2"

#Photoperiod and AET/PET
clim_path        = "Analysis_input/Drivers"

#Soil
SoilTexture_path = "Analysis_input/Drivers/SoilTexture"

#Phenology
PEPcsv_path      = "PEP_data/PEPcsv"


# 2. Output
###########

PEP_drivers_path  = "Analysis_input/PEP_drivers_final/Individual_files"
PEP_drivers_path2 = "Analysis_input/PEP_drivers_final/Merged_file"
PEP_drivers_path3 = "Analysis_input/PEP_drivers_final/Missing_observations"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



## Phenology data
#################

PEP.df <- fread(paste(PEPcsv_path, "pepData.csv", sep="/")) %>%
  #order table
  arrange(species, pep_id, year) %>%
  #keep only Aesculus, Betula, Fagus, and Quercus 
  filter(species %in% c('Aesculus','Betula pen.','Fagus','Quercus')) %>% 
  group_by(timeseries) %>% 
  #delete groups with less than 15 rows
  filter(n() >= 15)%>%  
  #delete timeseries above 65 latitude
  filter(lat<65) %>% 
  #rename Betula pendula
  mutate(species    = dplyr::recode(species, `Betula pen.`="Betula"),
         timeseries = paste0(pep_id, '_', species),
         leaf_off_max = max(leaf_off)) %>%
  ungroup() 


## AET/PET map
##############

#annual AET/PET ratio from SPLASH model
AET_PET.raster = raster(paste(clim_path, "AET_PET_ratio_global_FULL_MODIS-C006_MOD15A2_v1.alpha_MEANANN.nc", sep="/"))


## CO2 data
###########

CO2.df = fread(paste(CO2_path, "CO2.csv", sep="/"))


## Photoperiod
##############

photo.df = fread(paste(clim_path, "Photoperiod.csv", sep="/"))


## Soil Texture
###############

SoilTexture.df = fread(paste(SoilTexture_path, "SoilTexture.csv", sep="/"))


## SPEI
#######

SPEI.df = fread(paste(clim_path, "SPEI.csv", sep="/"))


## Import daily climatic datasets from GLDAS
############################################

#define climate variables
vn <- c('Daily_Mean_Data_Tair_f_inst','Daily_Min_Data_Tair_f_inst','Daily_Max_Data_Tair_f_inst',
        'Daily_Data_Rainf_f_tavg','Daily_Data_Qair_f_inst',
        'Daily_Data_SoilMoi0_10cm_inst','Daily_Data_SoilMoi10_40cm_inst',
        'Daily_Data_Swnet_tavg','Daily_Data_Lwnet_tavg','Daily_Data_SWdown_f_tavg')
#create empty list
DataList <- replicate(length(vn),data.frame())
#loop through climate variables
for(i in 1:length(vn)) {
  #read data
  data = fread(paste0(GLDAS_path, "/", vn[i],".csv")) %>%
    #add site x year identifier
    mutate(site_year = paste0(pep_id, '_', year)) %>%
    #order table
    dplyr::select(site_year, pep_id, year, lat, lon, everything()) 
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



#Add w_max (soil texture-dependent difference between field capacity and wilting point [%])
PEP.df = merge(PEP.df, SoilTexture.df[,c('pep_id','w_max')], by='pep_id')

#Add plant functional type info
PEP.df$PFT <- "TBL" # T-BL-SG: Temperate broad-leaved summergreen tree
#PEP.df[PEP.df$species=='Larix',]$PFT <- "BNL" # B-NL-SG: Boreal needle-leaved summergreen tree

#Add AET-PET ratio (required for pmodel)
PEP.df = 
  #both tables together
  cbind(PEP.df, 
        # intersection
        data.frame(AET_PET=raster::extract(AET_PET.raster, PEP.df[, c("lon", "lat")])))

#remove stuff
rm(AET_PET.raster, SoilTexture.df)



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



# Identifier 1 (all site x year combinations)
PEP.df$site_year = paste0(PEP.df$pep_id,"_",PEP.df$year)

# Identifier 2 (all timeseries x year combinations)
PEP.df$ts_yr     = paste0(PEP.df$timeseries,"_",PEP.df$year)
timeseries_year  = unique(PEP.df$ts_yr)
  
# add PEP data (+plant functional type label) and photoperiod to list
DataList[[11]] = photo.df
DataList[[12]] = CO2.df
DataList[[13]] = SPEI.df
DataList[[14]] = PEP.df

rm(photo.df, CO2.df, data, PEP.df)
names(DataList)=c(vn,"photoperiod",'CO2','SPEI',"PEP")
names(DataList)
#[1] "Daily_Mean_Data_Tair_f_inst"    "Daily_Min_Data_Tair_f_inst"     "Daily_Max_Data_Tair_f_inst"    
#[4] "Daily_Data_Rainf_f_tavg"        "Daily_Data_Qair_f_inst"         "Daily_Data_SoilMoi0_10cm_inst" 
#[7] "Daily_Data_SoilMoi10_40cm_inst" "Daily_Data_Swnet_tavg"          "Daily_Data_Lwnet_tavg"         
#10] "Daily_Data_SWdown_f_tavg"       "photoperiod"                    "CO2"                           
#[13] "SPEI"                           "PEP"       


##############################################################################################################################################


################################
# Loop through all time-points #
################################


parallelCalc <- function(timeseries_years){ 

  # Subset input data by time-point
  #################################

  #phenology data
  pheno.sub  <- DataList[[14]][which(DataList[[14]]$ts_yr==timeseries_years),]
  
  #daily mean temperature
  TMEAN <- DataList[[1]][which(DataList[[1]]$site_year==pheno.sub$site_year),]%>% 
    dplyr::select(as.character(1:366))
  
  #Skip timeseries for which there is no data
  if (nrow(TMEAN)==0) {
    write.table(pheno.sub, file=paste0(PEP_drivers_path3, '/', timeseries_years, '.csv'), sep=',', row.names = F, col.names = T)
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
    PHOTO     <- DataList[[11]][which(DataList[[11]]$lat==pheno.sub$lat),][1]%>% 
      dplyr::select(as.character(1:366))
    
    #CO2 (monthly)
    CO2       <- as.data.frame(t(DataList[[12]][which(DataList[[12]]$site_year==pheno.sub$site_year),]%>%
      dplyr::select(as.character(1:12)))) %>%
      rename(CO2 = V1)%>%
      mutate(Month = as.numeric(1:12))
    
    
    ##############################################################################################################################################
    
    
    # Create table of daily climate
    ###############################
    
    # Generate sub-dataframe to store results
    factors.sub <- pheno.sub %>% 
      dplyr::select(pep_id,species,timeseries,year,lat,lon,alt,leaf_out,leaf_off,leaf_off_mean) %>%
      mutate(CO2 = CO2$CO2[6])
    
    # Define the current year in calendar units
    year      <- as.character(pheno.sub$year)
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
                             Photo    = as.numeric(PHOTO))

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
                                                     temperature.metric = "celsius")) %>%
      #Add CO2
      left_join(CO2, by = "Month")
      
    #set NAs to 0
    daily_vals[is.na(daily_vals)] <- 0.0001
    
    
    ##############################################################################################################################################
    
    
    # Get average daytime temperature (chillR package)
    ##################################################
    
    #Get hourly values
    hourly_vals = stack_hourly_temps(daily_vals, latitude=pheno.sub$lat)$hourtemps
    
    #get daytime temperature
    daytime_temp = data.frame(hourly_vals %>%
                                 group_by_at(vars(-c(Temp,Hour))) %>%
                                 #delete night hours
                                 filter(between(Hour,
                                                daylength(latitude=pheno.sub$lat,JDay=JDay[1])$Sunrise,
                                                daylength(latitude=pheno.sub$lat,JDay=JDay[1])$Sunset))%>%
                                 #summarise daytime hours
                                 summarise(Tday = mean(Temp))%>%
                                 ungroup() %>%
                                 #order
                                 dplyr::select(Tday) )
    
    #get nighttime temperature
    nighttime_temp = data.frame(hourly_vals %>%
                                  group_by_at(vars(-c(Temp,Hour))) %>%
                                  #delete day hours
                                  filter(!between(Hour,
                                                  daylength(latitude=pheno.sub$lat,JDay=JDay[1])$Sunrise,
                                                  daylength(latitude=pheno.sub$lat,JDay=JDay[1])$Sunset))%>%
                                  #summarise daytime hours
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
    DOY_off <- pheno.sub$leaf_off_mean
    
    # leaf-out
    DOY_out <- ifelse(pheno.sub$leaf_out >= solstice, solstice-1, pheno.sub$leaf_out)
   
    
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
        if(pheno.sub$species=="Betula") { 
          VPD_min       <- 1.0
          VPD_max       <- 4.0
        } else {
        if(pheno.sub$species=="Fagus") { 
          VPD_min       <- 0.6
          VPD_max       <- 3.0
        } else {
          #median of all broadleaf tree species in White et al. 2000
          VPD_min       <- 1.0
          VPD_max       <- 3.5
        } }
        
        # Estimate photoperiod thresholds based on the maximum and minimum values of the growing season
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
        max.rad = max(DataList[[8]][which(DataList[[8]]$pep_id==pheno.sub$pep_id), c(as.character(1:365))], na.rm=T)
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
        pa <- daily_vals$CO2[i]*p
        
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
        # depends on PFT
        if(pheno.sub$PFT=="TBL") { 
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
        # depends on PFT
        if (pheno.sub$PFT=="TBL") { 
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
        w_max <- pheno.sub$w_max
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
      
      #add VPD to daily table
      daily_vals$VPD = VPD_year *1000 #VPD in Pa
              
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
      field.capacity = max(DataList[[7]][which(DataList[[7]]$pep_id==pheno.sub$pep_id), c(as.character(1:365))], na.rm=T)
      
      ## P-model v1.0
      pmodel.df <- tibble(
        tc             = daily_vals$Tday,
        vpd            = daily_vals$VPD, #VPD in Pa 
        co2            = daily_vals$CO2,
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
      
      #----------------------------------------------------------------------------------------------
      
      #Store the results
      ##################
      
      daily_vals = daily_vals %>%
        mutate(GSI        = iGSI_year,          #photoperiod-influenced GSI
               GSIrad     = iGSIrad_year,       #radiation-influenced GSI
               Azani      = dA_totw_year,       #net daytime photosynthesis (Zani et al., water-stressed)
               Apm        = Apm,                #net daytime photosynthesis (p model)
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
      variable.names = c('Azani', 'Apm',
                         'GSI', 'GSIrad', 
                         'GDDday', 'GDDnight', 
                         'SWrad',
                         'Tday', 'Tnight', 
                         'Moist', 'Prcp', 'VPD')
      
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
        
        varname.LO.p30      <- paste(variable.names[i], "LO.p30",   sep=".")
        varname.LO.p60      <- paste(variable.names[i], "LO.p60",   sep=".")
        varname.LO.p90      <- paste(variable.names[i], "LO.p90",   sep=".")
        
        varname.LO.SOm10    <- paste(variable.names[i], "LO.SOm10", sep=".")
        varname.LO.SOm20    <- paste(variable.names[i], "LO.SOm20", sep=".")
        varname.LO.SOm30    <- paste(variable.names[i], "LO.SOm30", sep=".")
        varname.LO.SOm40    <- paste(variable.names[i], "LO.SOm40", sep=".")
        varname.LO.SOm50    <- paste(variable.names[i], "LO.SOm50", sep=".")
        varname.LO.SOm60    <- paste(variable.names[i], "LO.SOm60", sep=".")
        
        varname.LO.SOp10    <- paste(variable.names[i], "LO.SOp10", sep=".")
        varname.LO.SOp20    <- paste(variable.names[i], "LO.SOp20", sep=".")
        varname.LO.SOp30    <- paste(variable.names[i], "LO.SOp30", sep=".")
        varname.LO.SOp40    <- paste(variable.names[i], "LO.SOp40", sep=".")
        varname.LO.SOp50    <- paste(variable.names[i], "LO.SOp50", sep=".")
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
        
        if(variable.names[i] %in% c('Azani', 'Apm',
                                    'GSI', 'GSIrad', 'GDDday', 'GDDnight')){
          
          # Sums
          ######
          
          factors.sub = factors.sub %>%
            mutate(
              #seasonal
              !!varname.LO.SO    := sum(variable[DOY_out:solstice]),
              
              !!varname.LO.p30 := sum(variable[DOY_out:(DOY_out+30)]),
              !!varname.LO.p60 := sum(variable[DOY_out:(DOY_out+60)]),
              !!varname.LO.p90 := sum(variable[DOY_out:(DOY_out+90)]),
              
              !!varname.LO.SOm10 := ifelse(DOY_out<(solstice-10), sum(variable[DOY_out:(solstice-10)]), 0),
              !!varname.LO.SOm20 := ifelse(DOY_out<(solstice-20), sum(variable[DOY_out:(solstice-20)]), 0),
              !!varname.LO.SOm30 := ifelse(DOY_out<(solstice-30), sum(variable[DOY_out:(solstice-30)]), 0),
              !!varname.LO.SOm40 := ifelse(DOY_out<(solstice-40), sum(variable[DOY_out:(solstice-40)]), 0),
              !!varname.LO.SOm50 := ifelse(DOY_out<(solstice-50), sum(variable[DOY_out:(solstice-50)]), 0),
              !!varname.LO.SOm60 := ifelse(DOY_out<(solstice-60), sum(variable[DOY_out:(solstice-60)]), 0),
              
              !!varname.LO.SOp10 := sum(variable[DOY_out:(solstice+10)]),
              !!varname.LO.SOp20 := sum(variable[DOY_out:(solstice+20)]),
              !!varname.LO.SOp30 := sum(variable[DOY_out:(solstice+30)]),
              !!varname.LO.SOp40 := sum(variable[DOY_out:(solstice+40)]),
              !!varname.LO.SOp50 := sum(variable[DOY_out:(solstice+50)]),
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
        
        if(variable.names[i] %in% c("Tday","Tnight","Moist","SWrad","VPD")){
          
          # Means from fixed date
          #######################
          
          factors.sub = factors.sub %>%
            mutate(
              #seasonal
              !!varname.LO.SO    := mean(variable[equinox.Mar:solstice]),
              
              !!varname.LO.p30 := mean(variable[DOY_out:(DOY_out+30)]),
              !!varname.LO.p60 := mean(variable[DOY_out:(DOY_out+60)]),
              !!varname.LO.p90 := mean(variable[DOY_out:(DOY_out+90)]),
              
              !!varname.LO.SOm10 := mean(variable[equinox.Mar:(solstice-10)]),
              !!varname.LO.SOm20 := mean(variable[equinox.Mar:(solstice-20)]),
              !!varname.LO.SOm30 := mean(variable[equinox.Mar:(solstice-30)]),
              !!varname.LO.SOm40 := mean(variable[equinox.Mar:(solstice-40)]),
              !!varname.LO.SOm50 := mean(variable[equinox.Mar:(solstice-50)]),
              !!varname.LO.SOm60 := mean(variable[equinox.Mar:(solstice-60)]),
              
              !!varname.LO.SOp10 := mean(variable[equinox.Mar:(solstice+10)]),
              !!varname.LO.SOp20 := mean(variable[equinox.Mar:(solstice+20)]),
              !!varname.LO.SOp30 := mean(variable[equinox.Mar:(solstice+30)]),
              !!varname.LO.SOp40 := mean(variable[equinox.Mar:(solstice+40)]),
              !!varname.LO.SOp50 := mean(variable[equinox.Mar:(solstice+50)]),
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
        
        if(variable.names[i] %in% c("Prcp")){
          
          # Sums from fixed date
          #######################
          
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
      VariableMeanVector = c("Tday","Tnight","Moist","VPD","SWrad")
      VariableSumVector  = c('Azani', 'Apm',
                             'GSI', 'GSIrad', 
                             'GDDday', 'GDDnight',
                             'Prcp')
      
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
      write.table(factors.sub, file=paste0(PEP_drivers_path, '/', timeseries_years, '.csv'), sep=',', row.names =F, col.names = T)
      
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
length(list.files(path=PEP_drivers_path, pattern='.csv'))
length(list.files(path=PEP_drivers_path3, pattern='.csv'))

#Rbind files
climate.factors.table = rbindlist(lapply(list.files(path = PEP_drivers_path)[1:200000], 
                                         function(n) fread(file.path(PEP_drivers_path, n))))
climate.factors.table2 = rbindlist(lapply(list.files(path = PEP_drivers_path)[200001:length(list.files(path=PEP_drivers_path, pattern='.csv'))], 
                                          function(n) fread(file.path(PEP_drivers_path, n))))
climate.factors.table = rbind(climate.factors.table,climate.factors.table2)



##############################################################################################################################################
##############################################################################################################################################



###################
## Safe the data ##
###################



#Safe table
write.csv(climate.factors.table, paste(PEP_drivers_path2, "pep_drivers_data.csv", sep="/"))

#Remove individual files
do.call(file.remove, list(list.files(PEP_drivers_path, 
                                     full.names = TRUE)))

do.call(file.remove, list(list.files(PEP_drivers_path3, 
                                     full.names = TRUE)))



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


