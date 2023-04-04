


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Linear models for the satellite data (EOS50) ##############################################################
#############################################################################################################



#required packages
require(tidyverse)
require(data.table)
require(broom)
require(sjmisc)



##############################################################################################################################################
##############################################################################################################################################



#####################
## Set directories ##
#####################



# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/Remote_sensing/Analysis")

# paths
Drivers_path    = "Analysis_input/Drivers_final_EOS50/Merged_file"
output_path     = "Analysis_output/Data"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



#Phenology data frame
#####################

Pheno.df <- fread(paste(Drivers_path, "Remote_sensing_drivers_data_preseason.csv", sep="/")) 



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

#MidGreendown_DOY...senescence date (DOY)
#MidGreendownMean...mean timeseries senescence date (DOY)

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

#Tmin.PS15-90...mean daily minimum temperature 15 to 90 days prior to mean timeseries senescence date
#tmean.PS15-90...mean daily mean temperature 15 to 90 days prior to mean timeseries senescence date
#tmax.PS15-90...mean daily maximum temperature 15 to 90 days prior to mean timeseries senescence date

#Tday...time series-specific best preseason (R2-based) for daytime temperatures
#tnight...time series-specific best preseason (R2-based) for nighttime temperatures



##############################################################################################################################################
##############################################################################################################################################



#################
## Full models ##
#################



#Define variables
variables = c('GPPstart', 'LAIstart',
              'Apm', 'Azani', 
              'Tnight', 'Tday', 
              'SWrad')

#create List object to store results
DataList = replicate(length(variables), data.frame())
names(DataList) = variables


##############################################################################################################################################


#Loop through variables
#######################

for (i in 1:length(variables)){

  #delete pixels with no photosynthesis before solstice
  if (variables[i] %in% c("GPPstart","Azani","Apm"))
  Pheno.df2 = Pheno.df %>%
    group_by(geometry) %>%
    filter(!mean(!!as.name(paste0(variables[i],".LO.SO"))) < .1) %>%
    ungroup() else Pheno.df2 = Pheno.df
  
  
  #set equations
  ##############
  
  #define variable names
  covariates = paste0(variables[i], c('.LO.SO','.SO.SE')) 

  
  #full models
  equation1 = as.formula(paste("scale(MidGreendown_DOY) ~ ", 
                               paste0('scale(',covariates[1], ') + scale(', covariates[2], ')', 
                                                                             collapse="+"), 
                               '+ scale(Prcp.LO.SO) + scale(Prcp.SO.SE) + scale(CO2) + scale(Tday)', collapse=""))

  equation2 = as.formula(paste("scale(MidGreendown_DOY) ~ ", 
                               paste0('scale(',covariates[1], ') + scale(', covariates[2], ')')))
  
  equation3 = as.formula(paste("scale(MidGreendown_DOY) ~ ", 
                               paste0('scale(',covariates[1], ') + scale(Tday)')))
  
  ##############################################################################################################################################
  
  
  ##############
  #linear models
  ##############
  
  
  ModelResults.df = Pheno.df2 %>%
    group_by(geometry, Lat, Lon, LC_Type) %>%
    do({
      
      #run models
      ###########
      
      model1 = lm(equation1,  data=.)
      model2 = lm(equation2,  data=.)
      model3 = lm(equation3,  data=.)
      
      #create combined data frame
      ###########################
      
      data.frame(rbind(
        
        #Equations 1-3
        tidy(model1) %>% mutate(equation = 'full model1'), 
        tidy(model2) %>% mutate(equation = 'full model2'), 
        tidy(model3) %>% mutate(equation = 'full model3') )
        
        
        )})  %>%
    
    #add variable name and delete "scale()" from term column
    mutate(variable = variables[i],
           term     = gsub("scale","",term),
           term     = gsub("\\(|\\)","",term) ) %>%
    #delete intercept
    filter(!term %in% c("Intercept")) %>%
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
  mutate(LAIstart4   = rowSums(dplyr::select(.,c("LAIstart1","LAIstart2","LAIstart3","LAIstart4"))),
         GPPstart4   = rowSums(dplyr::select(.,c("GPPstart1","GPPstart2","GPPstart3","GPPstart4"))),
         Apm4        = rowSums(dplyr::select(.,c("Apm1","Apm2","Apm3","Apm4"))),
         Azani4      = rowSums(dplyr::select(.,c("Azani1","Azani2","Azani3","Azani4")))
         )


#create List object to store results
DataList = replicate(length(variables), data.frame())
names(DataList) = variables


##############################################################################################################################################


#Loop through covariate groups
##############################


for (i in 1:length(variables)){

  #delete pixels with no photosynthesis before solstice
  #####################################################
  
  if (variables[i] %in% c("GPPstart","LAIstart","Azani","Apm"))
    Pheno.df2 = Pheno.monthly.df %>%
      group_by(geometry) %>%
      filter(!mean(!!as.name(paste0(variables[i],"4"))) < .1) %>%
      ungroup() else Pheno.df2 = Pheno.monthly.df
      
  #create explanatory variables
  #############################
  
  covariates.monthly  = paste0(variables[i], c(1:10))
  
  #---------------------------------------------------------
  
  #set equations
  ##############
  
  if(variables[i] %in% c('GPPstart', 'LAIstart', 'SWrad', 'Apm', 'Azani')) {  
    equation = as.formula(paste("scale(MidGreendown_DOY) ~ ", paste('scale(', covariates.monthly[4:10], ')', collapse="+")))
  } else {
    equation = as.formula(paste("scale(MidGreendown_DOY) ~ ", paste('scale(', covariates.monthly, ')', collapse="+")))
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
      
      data.frame(tidy(model) %>% mutate(equation = 'monthly'))
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



######################
## Seasonal drivers ##
######################


#Covariates
###########

#Variable length (leaf-out influenced):
#--------------------------------------
#Apm...Daily net photosynthesis (p-model)
#Azani...Daily net photosynthesis (Zani model)
#GSI...photoperiod-influenced growing-season index
#GSIrad...radiation-influenced growing-season index
#GDDday...daytime degree days
#GDDnight...nighttime degree days
#SWrad...Radiation sum

#Fixed length:
#-------------
#Tday...mean daytime temperature
#Tnight...mean daytime temperature


#-------------------------------------------------------------


## Define covariate groups
seasons    = c('LO.SOm30', 'LO.SO', 'LO.SOp30', 'LO.SOp60', 'LO.SE', 'SOm30.SE', 'SO.SE', 'SOp30.SE')
solstice   = c('solstice1', 'solstice2', 'solstice3', 'solstice4', 'solstice5', 'solstice6')

covariates1 = paste(rep(variables, each=length(seasons)),  seasons,  sep = '.')
covariates2 = paste(rep(variables, each=length(solstice)), solstice, sep = '.')
covariates  = c(covariates1,covariates2)

#Check if all variables are in dataframe
table(names(Pheno.df) %in% covariates)[2]/length(covariates)==1

#-------------------------------------------------------------

## Create List object to store results
DataList = replicate(length(covariates), data.frame())
names(DataList) = covariates
i=1


##############################################################################################################################################


#Loop through covariates
########################

for (covariate in covariates){

  #get variable name
  variable = gsub("\\..*","", covariate)
  
  #delete pixels with no photosynthesis for the respective period
  if (variable %in% c("GPPstart","Azani","Apm"))
    Pheno.df2 = Pheno.df %>%
      group_by(geometry) %>%
      filter(!mean(!!as.name(covariates[i])) < .1) %>%
      ungroup() else Pheno.df2 = Pheno.df
      
  #---------------------------------------------------------
  
  #set equations
  ##############
  
  if(grepl("solstice",covariate) & variable %in% c("GPPstart","Azani","Apm","LAIstart","SWrad")) {
    #autumn-temperature controlled scaled
    equation = as.formula(paste("scale(MidGreendown_DOY) ~ ", paste('scale(', covariate, ')', collapse="+"), 
                                 '+ scale(Tday)', collapse=""))
  } else {
    #univariate scaled
    equation = as.formula(paste("scale(MidGreendown_DOY) ~ ", paste('scale(', covariate, ')', collapse="+"), collapse=""))
  }
      
  
  ##############################################################################################################################################
  
  
  ##################
  #Run linear models
  ##################
  

  
  ModelResults.df = Pheno.df2 %>% 
    
    group_by(geometry, Lat, Lon, LC_Type) %>%
    
    do({
      
      #run models
      ###########
      
      model = lm(equation, data=.)
      
      #create combined dataframe
      ##########################
      
      data.frame(
        
        #Equation 1
        tidy(model) %>% 
          add_column(equation = paste0(ifelse(str_contains(covariate, c("solstice"), logic="or"),"Solstice","Seasonal"),
                                       '.scaled') ) %>% 
          filter(term %in% paste0('scale(',covariate,')')) )#delete intercept
      
    }) %>%
    
    #add variable name
    mutate(term = covariate,
           variable = sub("\\..*", "", covariate))
  
  
  ##############################################################################################################################################
  
  
  #store dataframe in variable list
  DataList[[i]] = ModelResults.df 
  
  print(paste0('..... ',i, ' out of ', length(covariates), ' done'))
  i=i+1
}

#bind tables
SeasonalAnalysis.df = bind_rows(DataList) 



##############################################################################################################################################
##############################################################################################################################################



##########
## Safe ##
##########



#bind full model, monthly and seasonal analyses
Analysis.df = rbind(FullAnalysis.df, MonthlyAnalysis.df, SeasonalAnalysis.df)

#Safe tables
write.csv(Analysis.df, paste(output_path, "Spatial_effect_data.csv", sep="/"))



##############################################################################################################################################
##############################################################################################################################################



#####################
## Reproducibility ##	
#####################



## datetime
Sys.time()
#"2021-12-08 14:21:44 CET"


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
#  [1] xfun_0.24         tidyselect_1.1.1  sjlabelled_1.1.8  haven_2.4.1       colorspace_2.0-1  vctrs_0.3.8      
#[7] generics_0.1.0    htmltools_0.5.1.1 yaml_2.2.1        utf8_1.2.1        rlang_0.4.11      pillar_1.6.1     
#[13] glue_1.4.2        withr_2.4.2       DBI_1.1.1         dbplyr_2.1.1      modelr_0.1.8      readxl_1.3.1     
#[19] plyr_1.8.6        lifecycle_1.0.0   munsell_0.5.0     gtable_0.3.0      cellranger_1.1.0  rvest_1.0.0      
#[25] evaluate_0.14     knitr_1.33        fansi_0.5.0       Rcpp_1.0.6        scales_1.1.1      backports_1.2.1  
#[31] jsonlite_1.7.2    fs_1.5.0          hms_1.1.0         digest_0.6.27     stringi_1.6.2     insight_0.14.2   
#[37] grid_4.1.0        cli_2.5.0         tools_4.1.0       magrittr_2.0.1    crayon_1.4.1      pkgconfig_2.0.3  
#[43] ellipsis_0.3.2    xml2_1.3.2        reprex_2.0.0      lubridate_1.7.10  assertthat_0.2.1  rmarkdown_2.9    
#[49] httr_1.4.2        rstudioapi_0.13   R6_2.5.0          compiler_4.1.0   



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


