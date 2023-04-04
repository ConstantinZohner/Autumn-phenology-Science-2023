


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Linear models for the satellite data (EOS85) ##############################################################
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
Drivers_path    = "Analysis_input/Drivers_final_EOS85/Merged_file"
output_path     = "Analysis_output_EOS85/Data"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



#Phenology data frame
#####################

Pheno.df <- fread(paste(Drivers_path, "Remote_sensing_drivers_data_EOS85_preseason.csv", sep="/")) %>% 
  filter(Dormancy_DOY>MidGreendown_DOY,
         Dormancy_DOY>Senesc_DOY,
         MidGreendown_DOY>Senesc_DOY) %>% 
  mutate(Duration_EOS10_50 = MidGreendown_DOY - Senesc_DOY,
         Duration_EOS10_85 = Dormancy_DOY - Senesc_DOY,
         Duration_EOS50_85 = Dormancy_DOY - MidGreendown_DOY) %>% 
  filter(Duration_EOS10_50 < 150,
         Duration_EOS10_85 < 175,
         Duration_EOS50_85 < 125) %>% 
  group_by(geometry, Lat, Lon) %>%
  filter(n() >= 15) %>% 
  ungroup()



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

#Dormancy_DOY...senescence date (DOY)

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
#GPP...Daily net photosynthesis
#LAI...Daily LAI

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



################################
## Senescence duration models ##
################################



#Define variables
predictors = c('Tday.EOS10.PS.60', 'Tday.EOS10.PS.90', "Tday.EOS50.PS.60",
               'Tday.EOS10mean_50', 'Tday.EOS10mean_85', "Tday.EOS50mean_85") 

responses = rep(c('Duration_EOS10_50', 'Duration_EOS10_85', "Duration_EOS50_85"),2)

#create List object to store results
DataList = replicate(length(predictors), data.frame())
names(DataList) = predictors


##############################################################################################################################################


#Loop through variables
#######################

for (i in 1:length(predictors)){
 
      #set equation
      #############
      
      equation = as.formula(paste(responses[i], " ~ ", predictors[i]))
      
      
      ##############################################################################################################################################
      
      
      ##############
      #linear models
      ##############
      
      
      ModelResults.df = Pheno.df %>%
        group_by(geometry, Lat, Lon, LC_Type) %>%
        do({
          
          #run model
          ##########
          
          model = lm(equation,  data=.)
      
          
          #create combined data frame
          ###########################
          
          data.frame(tidy(model))})  %>%
        
        #add variable name
        mutate(response = responses[i],
               predictor = predictors[i]) %>%
        #delete intercept
        filter(!term %in% c("(Intercept)")) %>%
        ungroup()
      
      
      ##############################################################################################################################################
      
      
      #store data frame in variable list
      DataList[[i]] = ModelResults.df 
      
      #count
      print(paste0('...',i,' out of ',length(predictors), ' (',predictors[i],') done'))
}

#bind rows
DurationAnalysis.df = bind_rows(DataList) 

#Safe table
write.csv(DurationAnalysis.df, paste(output_path, "Duration_data.csv", sep="/"))



##############################################################################################################################################
##############################################################################################################################################



#################
## Full models ##
#################



#Define variables
variables = c('GPPstart', 'LAIstart',
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
  equation.scaled1 = as.formula(paste("scale(Dormancy_DOY) ~ ", paste0('scale(',covariates[1], ') + scale(', covariates[2], ')', 
                                                     collapse="+"), 
                               '+ scale(Prcp.LO.SO) + scale(Prcp.SO.SE) + scale(CO2)', collapse=""))
  
  equation.scaled2 = as.formula(paste("scale(Dormancy_DOY) ~ ", paste0('scale(',covariates[1], ') + scale(', covariates[2], ')', 
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
         GPPstart4   = rowSums(dplyr::select(.,c("GPPstart1","GPPstart2","GPPstart3","GPPstart4")))
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
  
  if (variables[i] %in% c("GPPstart","LAIstart"))
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
      
      if(variables[i] %in% c('GPPstart', 'LAIstart', 'SWrad')) {  
        equation = as.formula(paste("scale(Dormancy_DOY) ~ ", paste('scale(', covariates.monthly[4:9], ')', collapse="+"), 
                                            collapse=""))
      } else {
        equation = as.formula(paste("scale(Dormancy_DOY) ~ ", paste('scale(', covariates.monthly, ')', collapse="+"), 
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
  if (variable %in% c("GPPstart"))
    Pheno.df2 = Pheno.df %>%
      group_by(geometry) %>%
      filter(!mean(!!as.name(covariates[i])) < .1) %>%
      ungroup() else Pheno.df2 = Pheno.df
      
      #---------------------------------------------------------
      
      #set equations
      ##############
      
      
      #univariate scaled
      equation = as.formula(paste("scale(Dormancy_DOY) ~ ", paste('scale(', covariate, ')', collapse="+"), collapse=""))
      
      
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
          
          data.frame(tidy(model) %>% 
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
#############################################################THE END##########################################################################
##############################################################################################################################################


