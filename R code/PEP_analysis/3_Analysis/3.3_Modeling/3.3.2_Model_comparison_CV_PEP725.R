


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Multivariate pre-/post-solstice models (PEP725 data) - Leave-one-out cross validation #####################
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
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_analysis/Analysis")

# paths
PEP_drivers_path    = "Analysis_input/PEP_drivers_final/Merged_file"
output_path         = "Analysis_output/Autumn/Data"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



#Phenology dataframe
####################

PEP.df <- fread(paste(PEP_drivers_path, "pep_drivers_data_preseason.csv", sep="/")) %>%
  mutate(SWrad.LO.SO = rowSums(.[,363:365]))



##############################################################################################################################################
##############################################################################################################################################



#######################################
## Timeseries-level model assessment ##
#######################################



#set equations
##############
  
equation.full = as.formula("leaf_off ~ Azani.LO.SO + Tday.LO.SO + SWrad.LO.SO + Moist.LO.SO + leaf_out +
                                       Azani.SO.SE + Tnight     + SWrad.SO.SE + Moist.SO.SE") 
equation.pre  = as.formula("leaf_off ~ Azani.LO.SO + Tday.LO.SO + SWrad.LO.SO + Moist.LO.SO + leaf_out")
equation.post = as.formula("leaf_off ~ Azani.SO.SE + Tnight     + SWrad.SO.SE + Moist.SO.SE")
  
#---------------------------------------------------------
  
###############
#Get model info
###############
  
ModelResults.df = PEP.df %>%
    group_by(species, timeseries)%>%
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
#############################################################THE END##########################################################################
##############################################################################################################################################


