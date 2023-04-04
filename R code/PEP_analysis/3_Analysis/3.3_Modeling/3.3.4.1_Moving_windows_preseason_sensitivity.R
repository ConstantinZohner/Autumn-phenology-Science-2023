


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Autumn temperature (preseason) moving-window analysis for the PEP725 data set #############################
#############################################################################################################



#required packages
require(tidyverse)
require(data.table)
require(broom)
require(gmodels)
require(patchwork)
require(MASS)



##############################################################################################################################################
##############################################################################################################################################



##########################################
## Set directories and get own PEP data ##
##########################################



# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_analysis/Analysis")

# paths
PEP_drivers_path = "Analysis_input/PEP_drivers_final/Merged_file"
output_path = "Analysis_output/Autumn/Data"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



PEP.df <- fread(paste(PEP_drivers_path, "pep_drivers_data.csv", sep="/"))%>% 
  dplyr::select(-V1) 

#------------------------------------------

#Strong filter dataframe
########################
PEPshort.df <- PEP.df %>%
  #delete years before 1980
  filter(!year<1980)%>%
  #delete groups with less than 30 years
  group_by(timeseries)%>%
  filter(n() >= 30)%>% 
  ungroup()

#------------------------------------------

#delete high elevation (>600 m) sites in full dataframe
PEPlong.df = PEP.df %>% filter(!alt>600)



##############################################################################################################################################
##############################################################################################################################################



############################
## Moving window analysis ##
############################



#create year vector for loop
year.vector.long  = c(1966:(max(PEPlong.df$year)-19))
year.vector.short = c(min(PEPshort.df$year):(max(PEPshort.df$year)-14))

#create List object to store results
DataList = replicate(2, data.frame())
names(DataList) = c("Long","Short")

#create List object of dataframes
PEPdataList = list(PEPlong.df, PEPshort.df)
#list of year vectors
YearList    = list(year.vector.long, year.vector.short)
#moving window length (20 / 15 years)
MovingWindowLength = c(20,15)



#########################################################################################################################
#########################################################################################################################



################################
# Run univariate linear models #
################################



for(k in 1:length(PEPdataList)){

    #create moving window dataframes
    mw.df = data.frame()
    mw.all.df = data.frame()
    
    #loop through years (moving windwows)
    for (Year in YearList[[k]]){
 
      #create table subset
      PEP.df.sub = PEPdataList[[k]] %>%
        filter(year >= Year, 
               year < Year+MovingWindowLength[k]) %>% 
        group_by(timeseries) %>%  
        filter(if (k==1) {n() >= 15} else {n() >= 12}
        ) %>% #delete time series with less than 15/12 years
        ungroup()
      
      
      #reshape table to long format
      #############################
      
      preseason.df = PEP.df.sub %>%
        #select columns
        dplyr::select(timeseries,year,species,pep_id,leaf_off,leaf_off_mean,
                      Tnight.PS.10,Tnight.PS.20,Tnight.PS.30,
                      Tnight.PS.40,Tnight.PS.50,Tnight.PS.60,
                      Tnight.PS.70,Tnight.PS.80,Tnight.PS.90,
                      Tnight.PS.100,Tnight.PS.110,Tnight.PS.120)%>%
        #long format
        pivot_longer(.,cols=starts_with(c("Td","Tn")), names_to = "preseason", values_to = "temp") %>%
        #create preseason length and temperature class columns
        mutate(preseason_length = readr::parse_number(gsub("[.]", "", preseason)) #keep only numbers in string
               #,temp_class = gsub("\\..*","", preseason)
               ) %>%
        dplyr::select(-preseason)
      
      
      #Run linear models
      ##################
      
      resultsLM = preseason.df %>% 
        group_by(timeseries, species, pep_id, preseason_length, leaf_off_mean) %>% 
        do({model = lm(scale(leaf_off) ~ scale(temp), data=.)  # linear model
        data.frame(tidy(model),         # coefficient info
                   glance(model))}) %>% # model info
        filter(!term %in% c("(Intercept)")) %>%
        dplyr::select(timeseries,species,pep_id,preseason_length,leaf_off_mean,estimate,r.squared)%>%
        ungroup()
      
      #--------------------------------------------------------------------------------------------------------------------
      
      #keep only models with best R2 predictions
      ##########################################
      
      #Species-specific
      resultsLM2 = resultsLM %>%   
        group_by(timeseries) %>% 
        top_n(1, r.squared) %>%
        ungroup() %>%
        mutate(preseason_start = leaf_off_mean - preseason_length) %>%
        #Summarize by species
        group_by(species) %>%
        summarise(length       = mean(preseason_length),
                  length.lowCI = t.test(preseason_length)$conf.int[1],
                  length.hiCI  = t.test(preseason_length)$conf.int[2],
                  start        = mean(preseason_start),
                  start.lowCI  = t.test(preseason_start)$conf.int[1],
                  start.hiCI   = t.test(preseason_start)$conf.int[2]) %>%
        ungroup() 
      
      #All species
      resultsLM2all = resultsLM %>%   
        group_by(timeseries) %>% 
        top_n(1, r.squared) %>%
        ungroup() %>%
        mutate(preseason_start = leaf_off_mean - preseason_length) %>%
        #summarize all
        summarise(length       = mean(preseason_length),
                  length.lowCI = t.test(preseason_length)$conf.int[1],
                  length.hiCI  = t.test(preseason_length)$conf.int[2],
                  start        = mean(preseason_start),
                  start.lowCI  = t.test(preseason_start)$conf.int[1],
                  start.hiCI   = t.test(preseason_start)$conf.int[2]) %>%
        mutate(species="Aall")
      
      #rbind species-specific and all species results
      resultsLM2 = rbind(resultsLM2, resultsLM2all) %>%
        mutate(variable = "R2",
               year=Year)
      
      #--------------------------------------------------------------------------------------------------------------------
      
      #keep only models with highest coefficients
      ###########################################
      
      #Species-specific
      resultsLM3 = resultsLM %>%   
        group_by(timeseries) %>% 
        top_n(1, estimate) %>%
        ungroup() %>%
        #Summarize by species
        mutate(preseason_start = leaf_off_mean - preseason_length) %>%
        group_by(species) %>%
        summarise(length       = mean(preseason_length),
                  length.lowCI = t.test(preseason_length)$conf.int[1],
                  length.hiCI  = t.test(preseason_length)$conf.int[2],
                  start        = mean(preseason_start),
                  start.lowCI  = t.test(preseason_start)$conf.int[1],
                  start.hiCI   = t.test(preseason_start)$conf.int[2]) %>%
        ungroup() 
      
      #All species
      resultsLM3all = resultsLM %>%   
        group_by(timeseries) %>% 
        top_n(1, estimate) %>%
        ungroup() %>%
        mutate(preseason_start = leaf_off_mean - preseason_length) %>%
        #summarize all
        summarise(length       = mean(preseason_length),
                  length.lowCI = t.test(preseason_length)$conf.int[1],
                  length.hiCI  = t.test(preseason_length)$conf.int[2],
                  start        = mean(preseason_start),
                  start.lowCI  = t.test(preseason_start)$conf.int[1],
                  start.hiCI   = t.test(preseason_start)$conf.int[2]) %>%
        mutate(species="Aall")
      
      #rbind species-specific and all species results
      resultsLM3 = rbind(resultsLM3, resultsLM3all) %>%
        mutate(variable = "coefficient",
               year=Year)
      
      #--------------------------------------------------------------------------------------------------------------------
      
      #rbind moving window subsets  
      mw.df = rbind(mw.df, resultsLM2, resultsLM3)
      
      #--------------------------------------------------------------------------------------------------------------------
      
      print(paste0('year ', Year, ' done (', min(PEPdataList[[k]]$year),'-',(max(PEPdataList[[k]]$year)-19),') [Dataset ', k, ' of ', length(DataList),']'))
    }
      
      #rbind all data
      mw.all.df = bind_rows(mw.all.df, mw.df) %>% 
        #rename factors
        mutate(dataset = ifelse(k==1, 'Long', 'Short'))
      
      #---------------------------------------------------------
      
      #store dataframe in variable list
      DataList[[k]]  = mw.all.df
  }
      
      
#bind rows
MovingWindowAnalysis.df = bind_rows(DataList) 

#Safe table
write.csv(MovingWindowAnalysis.df, paste(output_path, "Moving_window_data_preseason.csv", sep="/"))
      


##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################

      
      