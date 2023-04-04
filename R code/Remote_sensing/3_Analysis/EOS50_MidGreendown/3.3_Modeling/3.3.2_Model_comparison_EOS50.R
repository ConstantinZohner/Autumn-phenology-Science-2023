


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Univariate pre-/post-solstice models (EOS50) ##############################################################
#############################################################################################################



#required packages
require(tidyverse)
require(data.table)
require(broom)



##############################################################################################################################################
##############################################################################################################################################



#####################
## Set directories ##
#####################



# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/Remote_sensing/Analysis")

# paths
drivers_path    = "Analysis_input/Drivers_final_EOS50/Merged_file"
output_path     = "Analysis_output/Data"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



#Phenology dataframe
####################

Pheno.df <- fread(paste(drivers_path, "Remote_sensing_drivers_data_preseason.csv", sep="/"))  %>%
  #delete pixels with no photosynthesis before solstice
  group_by(geometry) %>%
  filter(!(mean(GPPstart.LO.SO)<.1)) %>%
  ungroup()



##############################################################################################################################################
##############################################################################################################################################



#######################################
## Timeseries-level model assessment ##
#######################################



#variable vector
variables=c("GPPstart","Tday","SWrad","Moist","Greenup_DOY")

#create List object to store results
DataList = replicate(length(variables), data.frame())
names(DataList) = variables


#Loop through variables
#######################

for (i in 1:length(variables)){

  #define variable names
  if (variables[i] == "Greenup_DOY") {
    covariates = c('Greenup_DOY','CO2')
  } else {covariates = paste0(variables[i], c('.LO.SO','.SO.SE')) }
  
  
  #set equations
  ##############
  
  equation.pre  = as.formula(paste("MidGreendown_DOY ~ ", paste0(covariates[1])))
  equation.post = as.formula(paste("MidGreendown_DOY ~ ", paste0(covariates[2])))
  
  #---------------------------------------------------------
  
  ###############
  #Get model info
  ###############
  
  ModelResults.df = Pheno.df %>%
    group_by(LC_Type, geometry)%>%
    do({
      
      #run models
      ###########
      
      modelPre  = lm(equation.pre,  data=.)
      modelPost = lm(equation.post, data=.)
      
      
      #create combined dataframe
      ##########################
      
      data.frame(rbind(
        
        #Equation pre-solstice
        glance(modelPre) %>% 
          mutate(model='pre'),
        
        #Equation post-solstice
        glance(modelPost) %>% 
          mutate(model='post') ) )
    })%>%
    mutate(variable = variables[i]) %>%
    ungroup()
  
  #---------------------------------------------------------
  
  #store dataframe in variable list
  DataList[[i]] = ModelResults.df
  
  #count
  print(paste0('...',i,' out of ',length(variables), ' (',variables[i],') done'))
}

#bind rows
Analysis.df = bind_rows(DataList) 



##############################################################################################################################################
##############################################################################################################################################



##########
## Safe ##
##########



write.csv(Analysis.df, paste(output_path, "Model_R2_data.csv", sep="/"))



##############################################################################################################################################
##############################################################################################################################################



#####################
## Reproducibility ##	
#####################


## datetime
Sys.time()
#"2021-12-05 20:06:47 CET"

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
#  [1] broom_0.7.8       data.table_1.14.0 forcats_0.5.1     stringr_1.4.0     dplyr_1.0.7       purrr_0.3.4      
#[7] readr_1.4.0       tidyr_1.1.3       tibble_3.1.2      ggplot2_3.3.4     tidyverse_1.3.1  

#loaded via a namespace (and not attached):
#  [1] Rcpp_1.0.6       cellranger_1.1.0 pillar_1.6.1     compiler_4.1.0   dbplyr_2.1.1     tools_4.1.0      jsonlite_1.7.2  
#[8] lubridate_1.7.10 lifecycle_1.0.0  gtable_0.3.0     pkgconfig_2.0.3  rlang_0.4.11     reprex_2.0.0     cli_2.5.0       
#[15] rstudioapi_0.13  DBI_1.1.1        haven_2.4.1      xml2_1.3.2       withr_2.4.2      httr_1.4.2       gtools_3.9.2    
#[22] fs_1.5.0         generics_0.1.0   vctrs_0.3.8      hms_1.1.0        grid_4.1.0       tidyselect_1.1.1 glue_1.4.2      
#[29] R6_2.5.0         fansi_0.5.0      readxl_1.3.1     gdata_2.18.0     modelr_0.1.8     magrittr_2.0.1   MASS_7.3-54     
#[36] gmodels_2.18.1   backports_1.2.1  scales_1.1.1     ellipsis_0.3.2   rvest_1.0.0      assertthat_0.2.1 colorspace_2.0-1
#[43] utf8_1.2.1       stringi_1.6.2    munsell_0.5.0    crayon_1.4.1 



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


