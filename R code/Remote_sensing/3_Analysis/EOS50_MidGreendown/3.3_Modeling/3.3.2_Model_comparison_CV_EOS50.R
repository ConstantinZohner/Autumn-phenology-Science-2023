


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Multivariate pre-/post-solstice models (EOS50) - Leave-one-out cross validation ###########################
#############################################################################################################



#required packages
require(tidyverse)
require(data.table)
require(broom)
require(caret)



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



#set equations
##############
  
equation.full = as.formula("MidGreendown_DOY ~ GPPstart.LO.SO + Tday.LO.SO + SWrad.LO.SO + Moist.LO.SO + Greenup_DOY +
                                               GPPstart.SO.SE + Tday +       SWrad.SO.SE + Moist.SO.SE") 
equation.pre  = as.formula("MidGreendown_DOY ~ GPPstart.LO.SO + Tday.LO.SO + SWrad.LO.SO + Moist.LO.SO + Greenup_DOY")
equation.post = as.formula("MidGreendown_DOY ~ GPPstart.SO.SE + Tday       + SWrad.SO.SE + Moist.SO.SE")
  
#---------------------------------------------------------
  
###############
#Get model info
###############
  
ModelResults.df = Pheno.df %>%
    group_by(LC_Type, geometry)%>%
    do({
      
      #run models
      ###########
      
      modelFull = lm(equation.full, data=.)
      modelPre  = lm(equation.pre,  data=.)
      modelPost = lm(equation.post, data=.)
      
      CVmodelFull <- train(
        equation.full, ., method = "lm",
        trControl = trainControl(method = "LOOCV") )
      
      CVmodelPre <- train(
        equation.pre, ., method = "lm",
        trControl = trainControl(method = "LOOCV") )
      
      CVmodelPost <- train(
        equation.post, ., method = "lm",
        trControl = trainControl(method = "LOOCV") )
      
      #create combined dataframe
      ##########################
      
      data.frame(rbind(
        
        #Equation full
        glance(modelFull) %>% 
          mutate(model = 'full', 
                 CV.R2 = as.numeric(CVmodelFull[4]$results[3])),
        
        #Equation pre-solstice
        glance(modelPre) %>% 
          mutate(model = 'pre', 
                 CV.R2 = as.numeric(CVmodelPre[4]$results[3])),
        
        #Equation post-solstice
        glance(modelPost) %>% 
          mutate(model ='post', 
                 CV.R2 = as.numeric(CVmodelPost[4]$results[3])) 
        ) )
    })%>%
    mutate(CV.R2 = ifelse(CV.R2 > r.squared, r.squared, CV.R2)) %>%
    ungroup()
  


##############################################################################################################################################
##############################################################################################################################################



##########
## Safe ##
##########



write.csv(ModelResults.df, paste(output_path, "Model_R2_CV_data.csv", sep="/"))



##############################################################################################################################################
##############################################################################################################################################



#####################
## Reproducibility ##	
#####################


## datetime
Sys.time()
#"2022-06-22 15:05:11 CEST"

## session info
sessionInfo()
#R version 4.1.0 (2021-05-18)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Big Sur 11.6.2

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] caret_6.0-92      lattice_0.20-44   broom_0.7.8       data.table_1.14.0 forcats_0.5.1     stringr_1.4.0    
#[7] dplyr_1.0.8       purrr_0.3.4       readr_1.4.0       tidyr_1.2.0       tibble_3.1.7      ggplot2_3.3.4    
#[13] tidyverse_1.3.1  

#loaded via a namespace (and not attached):
#  [1] httr_1.4.2           jsonlite_1.7.3       splines_4.1.0        foreach_1.5.2        prodlim_2019.11.13  
#[6] modelr_0.1.8         assertthat_0.2.1     stats4_4.1.0         cellranger_1.1.0     yaml_2.2.2          
#[11] globals_0.15.0       ipred_0.9-13         pillar_1.7.0         backports_1.2.1      glue_1.6.2          
#[16] pROC_1.18.0          digest_0.6.29        rvest_1.0.2          hardhat_1.1.0        colorspace_2.0-2    
#[21] recipes_0.2.0        htmltools_0.5.2      Matrix_1.3-3         plyr_1.8.6           timeDate_3043.102   
#[26] pkgconfig_2.0.3      listenv_0.8.0        haven_2.4.1          scales_1.1.1         gower_1.0.0         
#[31] lava_1.6.10          generics_0.1.2       ellipsis_0.3.2       withr_2.4.3          nnet_7.3-16         
#[36] cli_3.2.0            survival_3.2-11      magrittr_2.0.2       crayon_1.5.0         readxl_1.3.1        
#[41] evaluate_0.14        parallelly_1.32.0    fs_1.5.2             fansi_1.0.2          future_1.26.1       
#[46] nlme_3.1-152         MASS_7.3-54          xml2_1.3.3           class_7.3-19         tools_4.1.0         
#[51] hms_1.1.0            lifecycle_1.0.1      munsell_0.5.0        reprex_2.0.0         compiler_4.1.0      
#[56] rlang_1.0.2          grid_4.1.0           iterators_1.0.14     rstudioapi_0.13      rmarkdown_2.9       
#[61] ModelMetrics_1.2.2.2 gtable_0.3.0         codetools_0.2-18     DBI_1.1.2            reshape2_1.4.4      
#[66] R6_2.5.1             lubridate_1.7.10     knitr_1.33           fastmap_1.1.0        future.apply_1.9.0  
#[71] utf8_1.2.2           stringi_1.7.6        parallel_4.1.0       Rcpp_1.0.8           vctrs_0.4.1         
#[76] rpart_4.1-15         dbplyr_2.1.1         tidyselect_1.1.1     xfun_0.24   



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


