


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Pre-solstice effect moving-window analysis for the PEP725 data set ########################################
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



PEP.df <- fread(paste(PEP_drivers_path, "pep_drivers_data_new.csv", sep="/"))%>% 
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

#Define covariate groups
variables = c('Azani','Tday','Tnight','SWrad')

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

  mw.all.df = data.frame()
  mw.variable.df = data.frame()
  
  #Loop through covariate groups
  for (i in 1:length(variables)){
    
    #create moving window dataframes
    mw.year.df = data.frame()
    
    variable.names=c(paste0(variables[i],".LO.SOm60"), 
                     paste0(variables[i],".LO.SOm50"),
                     paste0(variables[i],".LO.SOm40"),
                     paste0(variables[i],".LO.SOm30"),
                     paste0(variables[i],".LO.SOm20"),
                     paste0(variables[i],".LO.SOm10"),
                     paste0(variables[i],".LO.SO"),
                     paste0(variables[i],".LO.SOp10"),
                     paste0(variables[i],".LO.SOp20"),
                     paste0(variables[i],".LO.SOp30"), 
                     paste0(variables[i],".LO.SOp40"), 
                     paste0(variables[i],".LO.SOp50"), 
                     paste0(variables[i],".LO.SOp60"))
    
    
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
      
      presolstice.df = PEP.df.sub %>%
        #select columns
        dplyr::select(timeseries,year,species,pep_id,leaf_out,leaf_off,leaf_off_mean,
                      all_of(variable.names)) %>%
        #long format
        pivot_longer(.,cols=starts_with(variables[i]), names_to = "period", values_to = "value")%>%
        #create preseason length and temperature class columns
        mutate(period_end = ifelse(period == paste0(variables[i],".LO.SOm60"), 172-60,
                                   ifelse(period == paste0(variables[i],".LO.SOm50"), 172-50, 
                                          ifelse(period == paste0(variables[i],".LO.SOm40"), 172-40, 
                                                 ifelse(period == paste0(variables[i],".LO.SOm30"), 172-30, 
                                                        ifelse(period == paste0(variables[i],".LO.SOm20"), 172-20, 
                                                               ifelse(period == paste0(variables[i],".LO.SOm10"), 172-10, 
                                                                      ifelse(period == paste0(variables[i],".LO.SO"), 172, 
                                                                             ifelse(period == paste0(variables[i],".LO.SOp10"), 172+10, 
                                                                                    ifelse(period == paste0(variables[i],".LO.SOp20"), 172+20, 
                                                                                           ifelse(period == paste0(variables[i],".LO.SOp30"), 172+30, 
                                                                                                  ifelse(period == paste0(variables[i],".LO.SOp40"), 172+40, 
                                                                                                         ifelse(period == paste0(variables[i],".LO.SOp50"), 172+50,172+60)))))))))))))
      
      
      #Run models
      ###########
      
      resultsLM = presolstice.df %>% 
        group_by(species, timeseries, period_end) %>% 
        do({coeff = cor(.$leaf_off, .$value)
        leafout = mean(.$leaf_out)
        data.frame(estimate = coeff,
                   leaf_out = leafout)}) %>% # model info
        dplyr::select(species,leaf_out,period_end,timeseries,estimate)%>%
        ungroup()

      
      #--------------------------------------------------------------------------------------------------------------------
      
      #keep only models with lowest estimates
      #######################################
      
      #Species-specific
      resultsLM2 = resultsLM %>%   
        group_by(timeseries) %>% 
        top_n(-1, estimate) %>%
        ungroup() %>%
        mutate(length = period_end - leaf_out) %>%
        #Summarize by species
        group_by(species) %>%
        summarise(start          = mean(leaf_out),
                  start.lowCI    = t.test(leaf_out)$conf.int[1],
                  start.hiCI     = t.test(leaf_out)$conf.int[2],
                  
                  end            = mean(period_end),
                  end.lowCI      = t.test(period_end)$conf.int[1],
                  end.hiCI       = t.test(period_end)$conf.int[2],
                  
                  duration       = mean(length),
                  duration.lowCI = t.test(length)$conf.int[1],
                  duration.hiCI  = t.test(length)$conf.int[2]) %>%
        ungroup() 
      
      #All species
      resultsLM2all = resultsLM %>%   
        group_by(timeseries) %>% 
        top_n(-1, estimate) %>%
        ungroup() %>%
        mutate(length = period_end - leaf_out) %>%
        #summarize all
        summarise(start          = mean(leaf_out),
                  start.lowCI    = t.test(leaf_out)$conf.int[1],
                  start.hiCI     = t.test(leaf_out)$conf.int[2],
                  
                  end            = mean(period_end),
                  end.lowCI      = t.test(period_end)$conf.int[1],
                  end.hiCI       = t.test(period_end)$conf.int[2],
                  
                  duration       = mean(length),
                  duration.lowCI = t.test(length)$conf.int[1],
                  duration.hiCI  = t.test(length)$conf.int[2]) %>%
        mutate(species="Aall")
      
      #rbind species-specific and all species results
      resultsLM2 = rbind(resultsLM2, resultsLM2all) %>%
        mutate(year     = Year,
               variable = variables[i],
               )
      
      #--------------------------------------------------------------------------------------------------------------------
      
      #rbind moving window subsets  
      mw.year.df = rbind(mw.year.df, resultsLM2)
      
      #--------------------------------------------------------------------------------------------------------------------
      
      print(paste0('year ', Year, ' done (', min(PEPdataList[[k]]$year),'-',(max(PEPdataList[[k]]$year)-19),') ', variables[i], ' [Dataset ', k, ' of ', length(DataList),']'))
    }
    
    #rbind all data
    mw.variable.df = bind_rows(mw.variable.df, mw.year.df) 
  }
  
  #rbind all data
  mw.all.df = bind_rows(mw.all.df, mw.variable.df) %>% 
    #rename factors
    mutate(dataset = ifelse(k==1, 'Long', 'Short'))
  
  #---------------------------------------------------------
  
  #store dataframe in variable list
  DataList[[k]]  = mw.all.df
}


#bind rows
MovingWindowAnalysis.df = bind_rows(DataList) 

#Safe table
write.csv(MovingWindowAnalysis.df, paste(output_path, "Moving_window_data_pre_solstice.csv", sep="/"))
      


##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################

      
      