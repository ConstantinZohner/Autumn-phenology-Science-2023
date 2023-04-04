require(data.table)
require(tidyverse)
require(dplyr)

## Set directories and get own PEP data

# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/Analysis")
# paths
GLDAS_path = "Analysis_input/Drivers/GLDAS/GLDAS_Extracted_Ver2"
#PEPcsv_path = "PEP_data/PEPcsv"
output_path = "Analysis_input/Drivers/GLDAS/GLDAS_Ordered_Ver2"

##----------------------------------------

## Phenology data
#PEP.df <- fread(paste(PEPcsv_path, "pepData.csv", sep="/")) %>%
#  mutate(site_year = paste0(pep_id, '_', year))

#identifiers
vn1 <- c('Tair_f_inst.csv','Tair_f_inst.csv','Tair_f_inst.csv',
         'Qair_f_inst.csv','Rainf_f_tavg.csv',
         'SWdown_f_tavg.csv','Swnet_tavg.csv','Lwnet_tavg.csv',
         'SoilMoi0_10cm_inst.csv','SoilMoi10_40cm_inst.csv')
vn2 <- c('Daily_Min_Data','Daily_Mean_Data','Daily_Max_Data',
         'Daily_Data','Daily_Data',
         'Daily_Data','Daily_Data','Daily_Data',
         'Daily_Data','Daily_Data')

#function to subset last three characters in string
subsetFunc = function(n) str_sub(n, -3,-1)
#year vector
year.vec   = c(1948:2019)

## Loop to safe tables (one for each variable)
for(i in 1:length(vn1)) {

  #create list of files to import
  list=intersect(list.files(pattern = vn1[i], path=GLDAS_path), 
                 list.files(pattern = vn2[i], path=GLDAS_path))
  
  #list all tables
  data=lapply(list, function(n) fread(file.path(GLDAS_path, n)))
  
  #add year to each list element
  for(x in 1:length(data)){
    data[[x]]$year=year.vec[x]
    }
  
  #rbind them
  data=rbindlist(data, fill=T)
  
  #data wrangling
  data = data %>% 
    dplyr::select(pep_id, year, lat, lon, everything(), -.geo, -`system:index`) %>% 
    rename_at(vars(-(1:4)), subsetFunc) %>% #keep only last three characters in column names
    rename_at(vars(-(1:4)), readr::parse_number) %>% #keep only numbers in column names
    #filter(site_year %in% PEP.df$site_year) %>% #match with PEP df
    filter(!is.na(`180`)) %>%
    arrange(pep_id,year)
  
  #safe table
  write.table(data, paste0(output_path,'/',vn2[i],"_",vn1[i]),sep=",",row.names=FALSE)
}
