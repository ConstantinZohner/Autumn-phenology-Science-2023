


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Linear models for the satellite data (EOSstart) ###########################################################
#############################################################################################################



#required packages
require(tidyverse)
require(ggplot2)
require(data.table)
require(broom)
require(sjmisc)



#define plot themes
plotTheme1 = theme(legend.position   = "right",
                   legend.background = element_rect(fill=NA, size=0.5, linetype="solid"),
                   legend.text       = element_text(color="black"),
                   panel.grid.major  = element_blank(),
                   panel.grid.minor  = element_blank(),
                   panel.background  = element_blank(),
                   panel.border      = element_rect(colour = "black", fill=NA),
                   axis.line         = element_line(color = "black"),
                   axis.text         = element_text(colour = "black"),
                   strip.background  = element_rect(fill=NA),
                   strip.text        = element_text(colour = 'black'),
                   plot.title        = element_text(face="bold"))



##############################################################################################################################################
##############################################################################################################################################



#####################
## Set directories ##
#####################



# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/Remote_sensing/Analysis")

# paths
Drivers_path    = "Analysis_input/Drivers_final_onset_VNP/Merged_file"
output_path     = "Analysis_output_startSen_VNP/Data"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



#Phenology data frame
#####################

Pheno.df <- fread(paste(Drivers_path, "Remote_sensing_drivers_data_onset_VNP_preseason.csv", sep="/")) %>% 
  #transform GPP to gCm-2
  mutate_at(c("GPPstart.LO.SO", 
              "GPPstart.SO.SE", 
              "GPPstart1", 
              "GPPstart2", 
              "GPPstart3", 
              "GPPstart4", 
              "GPPstart5", 
              "GPPstart6", 
              "GPPstart7", 
              "GPPstart8", 
              "GPPstart9", 
              "GPPstart10"),
            function(x)(x*0.1))



##############################################################################################################################################
##############################################################################################################################################



########################
## VARIABLE EXPLANATION:
########################

#General
########

#geometry...site identifier
#Year...observation year
#Lat...site latitude (decimal degrees)
#Lon...site longitude (decomal degrees)
#alt...site altitude (m a.s.l.)

#---------------------------------------------------

#Phenology
##########

#Onset_Greenness_Decrease...senescence date (DOY)

#---------------------------------------------------

#Day-of-year maximum temperature and radiation
##############################################

#HottestDOY...DOY with maximum annual temperature
#MaxRadDOY...DOY with maximum irradiance

#---------------------------------------------------

#Seasonal drivers
#################

#SUMS (flexible start = leaf-out):
#---------------------------------
#Apm...Daily net photosynthesis p-model
#Azani...Daily net photosynthesis Zani model
#GSI...photoperiod-influenced growing-season index
#GSIrad...radiation-influenced growing-season index
#GDDday...daytime degree days
#GDDnight...nighttime degree days
#SWrad...Radiation sum

#MEANS (fixed start = March equinox):
#------------------------------------
#Tday...mean daytime temperature
#Tnight...mean nighttime temperature
#Moist...mean soil moisture (10-40cm)
#Prcp...precipitation sum

#-----------------

#PERIODS

#seasonal
#--------
#LO.SOm30...leaf-out to 30 days before summer solstice
#LO.SO...leaf-out to summer solstice
#LO.SOp30...leaf-out to 30 days after solstice
#LO.SOp60...leaf-out to 60 days after solstice
#LO.SE...leaf-out to mean senescence
#SOm30.SE...30 days before solstice to mean senescence
#SO.SE...solstice to mean senescence
#SOp30.SE...30 days after solstice to mean senescence
#SOp60.SE...60 days after solstice to mean senescence

#around summer solstice
#----------------------
#solstice1...40 to 10 days before solstice
#solstice2...30 days before solstice
#solstice3...20 days before until 10 days after solstice
#solstice4...10 days before until 20 days after solstice
#solstice5...30 days after solstice
#solstice6...10 to 40 days after solstice

#---------------------------------------------------           

#Monthly drivers
################

#Tnight1-12...mean of nighttime temperatures in January (1) to December (12)
#Tday1-12...mean of daytime temperatures in January (1) to December (12)
#GDDnight1-12...sum of nighttime degree-days in January (1) to December (12)
#GDDday1-12...sum of daytime degree-days in January (1) to December (12)
#SWrad1-12...mean of daily short-wave radiation in January (1) to December (12)
#Moist1-12...mean of daily soil moisture (10-40cm) in January (1) to December (12)
#Prcp1-12...sum of daily precipitation in January (1) to December (12)
#Apm1-12...sum of daily gross photosynthesis in January (1) to December (12)
#Azani1-12...sum of daily net photosynthesis in January (1) to December (12)
#GSI...day length-influenced growing-season index
#GSIrad...radiation-influenced growing-season index

#---------------------------------------------------           

#Preseason temperatures
#######################

#Tday...time series-specific best preseason (R2-based) for daytime temperatures
#tnight...time series-specific best preseason (R2-based) for nighttime temperatures



##############################################################################################################################################
##############################################################################################################################################



#################
## Full models ##
#################



#Define variables
variables = c('GPPstart','Tday') 

#create List object to store results
DataList = replicate(length(variables), data.frame())
names(DataList) = variables


##############################################################################################################################################


#Loop through variables
#######################

for (i in 1:length(variables)){

  #delete pixels with no photosynthesis before solstice
  if (variables[i] %in% c('GPPstart'))
  Pheno.df2 = Pheno.df %>%
    group_by(geometry) %>%
    filter(!mean(!!as.name(paste0(variables[i],".LO.SO"))) < .1) %>%
    ungroup() else Pheno.df2 = Pheno.df
  
  
  #set equations
  ##############
  
  #define variable names
  covariates = paste0(variables[i], c('.LO.SO','.SO.SE')) 

  
  #full models
  equation.scaled1 = as.formula(paste("Onset_Greenness_Decrease ~ ", paste0(covariates[1], '+', covariates[2], 
                                                     collapse="+"), 
                               '+ Prcp.LO.SO + Prcp.SO.SE + CO2', collapse=""))
  
  equation.scaled2 = as.formula(paste("Onset_Greenness_Decrease ~ ", paste0(covariates[1], '+', covariates[2],
                                                                     collapse="+")))
  
  
  ##############################################################################################################################################
  
  
  ##############
  #linear models
  ##############
  
  
  ModelResults.df = Pheno.df2 %>%
    group_by(geometry, Lat, Lon, LC_Type) %>%
    do({
      
      #run model
      ##########
      
      model1 = lm(equation.scaled1,  data=.)
      model2 = lm(equation.scaled2,  data=.)
      
      #create combined data frame
      ###########################
      
      data.frame(rbind(
        
        #Equation 1
        tidy(model1) %>% mutate(equation = 'full model 1'),  #add model name
        
        #Equation 2
        tidy(model2) %>% mutate(equation = 'full model 2') 
        
        ))
          })  %>%
    
    #add variable name and delete "scale()" from term column
    mutate(variable = variables[i]) %>%
    #delete intercept
    filter(!term %in% c("(Intercept)")) %>%
    ungroup()
    
  
  ##############################################################################################################################################
  
  
  #store data frame in variable list
  DataList[[i]] = ModelResults.df 
  
  #count
  print(paste0('...',i,' out of ',length(variables), ' (',variables[i],') done'))
}

#bind rows
FullAnalysis.df = bind_rows(DataList) 



##############################################################################################################################################
##############################################################################################################################################



##########################
## Monthly correlations ##
##########################



#Get sums for January to April photosynthesis parameters
Pheno.monthly.df     = Pheno.df %>%
  mutate(GPPstart4   = rowSums(dplyr::select(.,c("GPPstart1","GPPstart2","GPPstart3","GPPstart4"))))


#create List object to store results
DataList = replicate(length(variables), data.frame())
names(DataList) = variables


##############################################################################################################################################


#Loop through covariate groups
##############################


for (i in 1:length(variables)){

  #delete pixels with no photosynthesis in April
  ##############################################
  
  if (variables[i] %in% c("GPPstart"))
    Pheno.df2 = Pheno.monthly.df %>%
      group_by(geometry) %>%
      filter(!mean(!!as.name(paste0(variables[i],"4"))) < .1) %>%
      ungroup() else Pheno.df2 = Pheno.monthly.df
      
      #create explanatory variables
      #############################
      
      covariates.monthly  = paste0(variables[i], c(1:9))
      
      #---------------------------------------------------------
      
      #set equations
      ##############
      
      if(variables[i] %in% c('GPPstart')) {  
        equation = as.formula(paste("Onset_Greenness_Decrease ~ ", paste(covariates.monthly[4:9], collapse="+"), 
                                            collapse=""))
      } else {
        equation = as.formula(paste("Onset_Greenness_Decrease ~ ", paste(covariates.monthly[3:9], collapse="+"), 
                                            collapse=""))
      }
      
      #---------------------------------------------------------
      
      ###############
      # Linear models
      ###############
      
      ModelResults.df = Pheno.df2 %>%
        group_by(geometry, Lat, Lon, LC_Type) %>%
        do({
          
          #run models
          ###########
          
          model = lm(equation, data=.)
  
          
          #create combined dataframe
          ##########################
          
          data.frame(tidy(model) %>% mutate(equation = 'monthly') )
          }) %>% 
        
        ungroup() %>%
        #delete intercept
        filter(!term %in% c("(Intercept)")) %>% 
        #add variable name and keep only numbers in month column 
        mutate(variable = variables[i],
               term = readr::parse_number(term)) 
      
      
      #---------------------------------------------------------
      
      #store dataframe in variable list
      DataList[[i]] = ModelResults.df 
      
      #count
      print(paste0('...',i,' out of ',length(variables), ' (',variables[i],') done'))
}

#bind rows
MonthlyAnalysis.df = bind_rows(DataList) 



##############################################################################################################################################
##############################################################################################################################################



##########
## Safe ##
##########



#bind full model, monthly and seasonal analyses
Analysis.df = rbind(FullAnalysis.df, MonthlyAnalysis.df)

#Safe tables
write.csv(Analysis.df, paste(output_path, "Spatial_effect_data_onset_VNP_no_scaling.csv", sep="/"))



##############################################################################################################################################
##############################################################################################################################################



#####################
## Reproducibility ##	
#####################



## datetime
Sys.time()
#"2021-12-05 19:41:36 CET"


## session info
sessionInfo()
#R version 4.1.0 (2021-05-18)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Big Sur 11.2.3

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] sjmisc_2.8.7      broom_0.7.8       data.table_1.14.0 forcats_0.5.1     stringr_1.4.0     dplyr_1.0.7      
#[7] purrr_0.3.4       readr_1.4.0       tidyr_1.1.3       tibble_3.1.2      ggplot2_3.3.4     tidyverse_1.3.1  

#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.6       cellranger_1.1.0 pillar_1.6.1     compiler_4.1.0   dbplyr_2.1.1     tools_4.1.0      jsonlite_1.7.2  
#[8] lubridate_1.7.10 lifecycle_1.0.0  gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.11     reprex_2.0.0     cli_2.5.0       
#[15] rstudioapi_0.13  DBI_1.1.1        haven_2.4.1      xml2_1.3.2       withr_2.4.2      httr_1.4.2       fs_1.5.0        
#[22] generics_0.1.0   vctrs_0.3.8      hms_1.1.0        sjlabelled_1.1.8 grid_4.1.0       tidyselect_1.1.1 glue_1.4.2      
#[29] R6_2.5.0         fansi_0.5.0      readxl_1.3.1     modelr_0.1.8     magrittr_2.0.1   backports_1.2.1  scales_1.1.1    
#[36] ellipsis_0.3.2   insight_0.14.2   rvest_1.0.0      assertthat_0.2.1 colorspace_2.0-1 utf8_1.2.1       stringi_1.6.2   
#[43] munsell_0.5.0    crayon_1.4.1    



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


