require(rgdal)
require(raster)
require(soiltexture)
require(data.table)
require(tidyverse)


##############################################################################################################################################


######################################
## Set directories and get PEP data ##
######################################


# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_analysis/Analysis")

## Paths

#Soil texture
SoilTexture_path  = "Analysis_input/Drivers/SoilTexture"

#Phenology
PEPcsv_path = "PEP_data/PEPcsv"


## load phenology data
PEPsites.df <- fread(paste(PEPcsv_path, "pepData.csv", sep="/")) %>%
  filter(!duplicated(pep_id)) %>%
  dplyr::select(pep_id, lat, lon)%>%
  arrange(pep_id)


##############################################################################################################################################


###############
## Soil data ##
###############


#Texture and clay content (for Zani photosynthesis model)

# Get coordinates
xy <- PEPsites.df %>% 
  dplyr::select(lon,lat)

# Import raster files of soil texture

# Soil images provided by ISRIC (World Soil Information)
# SoilGrids
# https://maps.isric.org/
Soils_coarse <- raster(paste0(SoilTexture_path,'/','Layers/SoilTexture_0cm.tif')) # Coarse Fragments Volumetric in % at surface
Soils_fine <- raster(paste0(SoilTexture_path,'/',"Layers/ClayContent_0cm.tif")) # Fine "Clay" content Mass Fraction in % at surface
plot(Soils_coarse)
plot(Soils_fine)

# Create a spatial object
spdf <- SpatialPointsDataFrame(coords = xy, data = PEPsites.df)

# Extract fragment values from rasters using coordinates 
proj4string(spdf) <- CRS("+init=epsg:4326")
Clay <- raster::extract(Soils_fine, spdf, cellnumbers = T)[,2]
Silt <- raster::extract(Soils_coarse, spdf, cellnumbers = T)[,2]
Sand <- 100-(Clay+Silt)

fragment.df <- data.frame(
  "CLAY" = Clay,
  "SILT" = Silt,
  "SAND" = Sand
)

# Get soil texture from fragment percentages
texture.df <- TT.points.in.classes(
  tri.data = fragment.df,
  class.sys = "HYPRES.TT",
  PiC.type = "l"
)

# Create dataframe of soil parameters and texture-to-parameter conversion
# Ref. Sitch, S. et al. Evaluation of ecosystem dynamics, plant geography and terrestrial carbon cycling in the LPJ dynamic global vegetation model. Glob. Chang. Biol. 9, 161-185 (2003).
soil_pars.df <- data.frame(
  "w_max" = c(13,13,21,14.8,10),
  "k_perc" = c(2,3,4,3.5,5),
  "c_soil" = c(0.1,0.1,0.02,0.02,0.2),
  row.names = c("VF","F","M","MF","C")
)

# Find soil parameters based on soil texture for each timeseries
soil_pars_ts.df <- data.frame()
for(row in 1:nrow(PEPsites.df)) {
  sub.df <- texture.df[row,]
  index <- min(which(sub.df == TRUE))
  pars_ts <- soil_pars.df[index,]
  pars_ts$pep_id <- PEPsites.df[row,]$pep_id
  soil_pars_ts.df <- rbind(soil_pars_ts.df,pars_ts)
  print(paste0("Soil parameters for ",PEPsites.df[row,]$pep_id," executed!"))
}

#remove stuff
rm(PEP_CO2.df, fragment.df, pars_ts, soil_pars.df, 
   Soils_coarse, Soils_fine, xy, texture.df, spdf, data, PEPsites.df)

hist(soil_pars_ts.df$w_max)
table(soil_pars_ts.df$w_max)

#safe soil texture table
write.csv(soil_pars_ts.df, paste(SoilTexture_path, "SoilTexture.csv", sep="/"))


##############################################################################################################################################
