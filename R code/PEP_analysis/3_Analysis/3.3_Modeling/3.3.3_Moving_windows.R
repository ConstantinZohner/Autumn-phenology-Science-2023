


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Moving window analysis for the PEP725 data set ############################################################
#############################################################################################################



#required packages
require(tidyverse)
require(data.table)
require(broom)
require(broom)
require(broom.mixed)
require(gmodels)
require(sjmisc)
require(pbmcapply)
require(pracma)
require(raster)
require(lme4)
require(car)
require(ggplot2)
require(wesanderson)
require(patchwork)


#plot theme
plotTheme1 = theme(
  legend.position   = "none",
  legend.background = element_rect(fill=NA, size=0.5, linetype="solid"),
  legend.text       = element_text(color="black"),
  panel.grid.major  = element_line(colour = "lightgrey"), 
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



# set the working dirctory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/PEP_analysis/Analysis")

# paths
PEP_drivers_path    = "Analysis_input/PEP_drivers_final/Merged_file"
GLDAS_path          = "Analysis_input/Drivers/GLDAS"
CO2_path            = "Analysis_input/Drivers/CO2"
output_path         = "Analysis_output/Autumn/Data"
output_path_figure  = "Analysis_output/Autumn/Moving_window"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



#Phenology dataframe
####################

PEP.df <- as.data.frame(fread(paste(PEP_drivers_path, "pep_drivers_data_preseason.csv", sep="/"))) %>%
  #Get sums for photosynthesis parameters
  mutate(Apm5        = rowSums(.[c("Apm1","Apm2","Apm3","Apm4",'Apm5')]),
         ApmJmaxA5   = rowSums(.[c("ApmJmaxA1","ApmJmaxA2","ApmJmaxA3","ApmJmaxA4","ApmJmaxA5")]),
         ApmJmaxB5   = rowSums(.[c("ApmJmaxB1","ApmJmaxB2","ApmJmaxB3","ApmJmaxB4","ApmJmaxB5")]),
         Azani5      = rowSums(.[c("Azani1","Azani2","Azani3","Azani4","Azani5")]),
         AzaniJmaxA5 = rowSums(.[c("AzaniJmaxA1","AzaniJmaxA2","AzaniJmaxA3","AzaniJmaxA4","AzaniJmaxA5")]),
         AzaniJmaxB5 = rowSums(.[c("AzaniJmaxB1","AzaniJmaxB2","AzaniJmaxB3","AzaniJmaxB4","AzaniJmaxB5")]),
         GSI5        = rowSums(.[c("GSI1","GSI2","GSI3","GSI4","GSI5")]),
         GSIrad5     = rowSums(.[c("GSIrad1","GSIrad2","GSIrad3","GSIrad4","GSIrad5")]),
         
         Apm8        = rowSums(.[c("Apm8","Apm9")]),
         ApmJmaxA8   = rowSums(.[c("ApmJmaxA8","ApmJmaxA9")]),
         ApmJmaxB8   = rowSums(.[c("ApmJmaxB8","ApmJmaxB9")]),
         Azani8      = rowSums(.[c("Azani8","Azani9")]),
         AzaniJmaxA8 = rowSums(.[c("AzaniJmaxA8","AzaniJmaxA9")]),
         AzaniJmaxB8 = rowSums(.[c("AzaniJmaxB8","AzaniJmaxB9")]),
         GSI8        = rowSums(.[c("GSI8","GSI9")]),
         GSIrad8     = rowSums(.[c("GSIrad8","GSIrad9")])
         )

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


## Daily climatic data from GLDAS averaged to site-level annual mean over all years (1948-2015)
###############################################################################################

#list of climate variables (Daily air temp, rainfall and short-wave radiation)
vn <- c('Daily_Mean_Data_Tair_f_inst','Daily_Data_Rainf_f_tavg','Daily_Data_Swnet_tavg')

#create empty dataframe
Climate.df <- data.frame()

#loop
for(i in 1:length(vn)) {

  data <- fread(paste0(GLDAS_path, "/", vn[i],".csv"))%>%
    #get annual means
    mutate(Climate = rowMeans(dplyr::select(., -c(pep_id,year,lat,lon)),na.rm=T))%>%
    dplyr::select(pep_id,year,lat,lon,Climate)%>%
    #average across all years (1948-2015)
    group_by(pep_id)%>%
    summarise(Climate = ci(Climate)[1])%>%
    #rename variable
    mutate(variable = vn[i])
    
  #Rbind climate variables  
  Climate.df = rbind(Climate.df, data)
}

#Wide format
Climate.df = pivot_wider(Climate.df, names_from = variable, values_from = Climate)%>%
  rename(MAT = Daily_Mean_Data_Tair_f_inst, MAP = Daily_Data_Rainf_f_tavg, RAD = Daily_Data_Swnet_tavg)%>%
  mutate(MAP = MAP*365)

#merge with PEP.df
PEPlong.df = PEPlong.df %>%
  left_join(., Climate.df, by='pep_id')

#remove stuff
rm(Climate.df, data)



##############################################################################################################################################
##############################################################################################################################################



############################
## Moving window analysis ##
############################



#create year vector for loop
year.vector.long  = c(1966:(max(PEPlong.df$year)-19))
year.vector.short = c(min(PEPshort.df$year):(max(PEPshort.df$year)-14))

#Define covariate groups
variables = c('Apm', 'ApmJmaxA', 'ApmJmaxB',
              'Azani', 'AzaniJmaxA', 'AzaniJmaxB', 
              'GSI', 'GSIrad', 
              'GDDday', 'GDDnight', 
              'Tday', 'Tnight', 'SWrad')

#create List object to store results
DataList = replicate(2*length(variables), data.frame())
names(DataList) = rep(variables,2)

#create List object of dataframes
PEPdataList = list(PEPlong.df, PEPshort.df)
#list of year vectors
YearList    = list(year.vector.long, year.vector.short)
#moving window length (20 / 15 years)
MovingWindowLength = c(20,15)



#########################################################################################################################
#########################################################################################################################



###########################################
# get sample sizes for each moving window #
###########################################


#get full data sample sizes per year and species
################################################

FullSampleData.df = PEP.df %>%
  group_by(species, year)%>%
  summarise(count = n())


#get sample sizes per year and species for moving windows
#########################################################

#empty dataframe
SampleData.df = data.frame()

#loop
for(k in 1:length(PEPdataList)){

  #loop through years (moving windwows)
  for (year in YearList[[k]]){

    #create table subset
    PEP.sub = PEPdataList[[k]][PEPdataList[[k]]$year>=year & PEPdataList[[k]]$year< year+MovingWindowLength[k],] %>% 
      #delete time series with less than 15/12 years
      group_by(timeseries) %>%  
      filter(if (k==1) {n() >= 15} else {n() >= 12} ) %>% 
      ungroup()%>% 
      #get sample size info
      group_by(species)%>%
      summarise(count = n())%>%
      #add identifiers
      mutate(year=year,
             dataset = ifelse(k==1, 'Long', 'Short'))
      
    #rbind
    SampleData.df = rbind(SampleData.df, PEP.sub)

    #count
    print(paste0('year ', year,', dataset ', k))
  }
}

#------------------------------------------------------------------------------------------------

#Plots
######

#Map of the observations
mp <- NULL
mapWorld <- borders("world", colour="gray60", fill="gray60") # create a layer of borders
mp <- ggplot() + mapWorld + plotTheme1
#Layer the stations on top
mp <- mp + geom_point(data=PEP.df[!duplicated(PEP.df[ , c("lat", "lon")]), ], 
                      aes(x=lon, y=lat) ,color="blue", size=.3) +
  coord_cartesian(ylim = c(43, 60), xlim = c(-8, 28))+
  xlab("") + ylab('') 

#Full data
SampleSizePlot = ggplot(FullSampleData.df, 
                            aes(x = year, y = count,
                                group=species, color=species)) + 
  geom_line(size = 1) +
  scale_color_manual(values = rev(wes_palette("Darjeeling2", n = 4))) +
  xlab("") + ylab('Sample size')+
  plotTheme1+
  coord_cartesian(xlim=c(1954,2012), ylim=c(115,2500))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#Long data
SampleSizeLongPlot = ggplot(SampleData.df[SampleData.df$dataset=='Long',], 
                        aes(x = year, y = count/20,
                            group=species, color=species)) + 
  geom_line(size = 1) +
  scale_color_manual(values = rev(wes_palette("Darjeeling2", n = 4))) +
  xlab("") + ylab('Sample size')+
  plotTheme1+
  scale_x_continuous(breaks = seq(1966,1996,by=5),
                     labels = c('1966-1985','1971-1990','1976-1995','1981-2000','1986-2005','1991-2010','1996-2015'))+
  coord_cartesian(xlim=c(1967.35,1994.65), ylim=c(115,2500))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#Short data
SampleSizeShortPlot = ggplot(SampleData.df[SampleData.df$dataset=='Short',], 
                            aes(x = year, y = count/15,
                                group=species, color=species)) + 
  geom_line(size = 1) + 
  scale_color_manual(values = rev(wes_palette("Darjeeling2", n = 4))) +
  xlab("") + ylab('')+
  plotTheme1+
  theme(legend.position = 'right')+
  scale_x_continuous(breaks = seq(1981,2001,by=5),
                     labels = c('1981-1995','1986-2000','1991-2005','1996-2010','2001-2015'))+
  coord_cartesian(xlim=c(1981.5,2000.08), ylim=c(115,2500))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#define plot layout
layout <- "AB
CD"

#Merge plots
SampleSize_Plot = mp + SampleSizePlot + SampleSizeLongPlot + SampleSizeShortPlot +
  plot_layout(design = layout) + plot_annotation(tag_levels = 'a')&
  theme(plot.tag = element_text(face = 'bold'))

#save plots as .pdf
ggsave(SampleSize_Plot, file="SampleSizes.pdf", path=output_path_figure,
       width=8, height=6)

#remove stuff
rm(PEP.df, SampleData.df, PEP.sub, SampleSize_Plot, SampleSizePlot, SampleSizeLongPlot, SampleSizeShortPlot, mp, FullSampleData.df, mapWorld)



#########################################################################################################################
#########################################################################################################################



######################################################
# CHECK bias in site-level average climate over time #
######################################################



#create data frame to store results
site.mw.df = data.frame()

#choose long moving window analysis
k=1

#loop through years (moving windwows)
for (year in YearList[[k]]){
  
  #create table subset
  PEP.df.sub = PEPdataList[[k]][PEPdataList[[k]]$year >= year & PEPdataList[[k]]$year < year+MovingWindowLength[k],] %>% 
    #delete time series with less than 15/12 years
    group_by(timeseries) %>%  
    filter(if (k==1) {n() >= 15} else {n() >= 12}
    ) %>% 
    ungroup() %>%
    #delete duplicates (unique site/species per moving window)
    distinct(pep_id, species, .keep_all = T) 
  
  #get means
  ResultsMean.df = PEP.df.sub %>%
    summarise(MAT = ci(MAT)[1], 
              MAP = ci(MAP)[1], 
              RAD = ci(RAD)[1],
              ELE = ci(alt)[1])%>%
    mutate(year=year)%>% 
    pivot_longer(., cols = -year)%>% 
    rename(mean.climate = value)
  
  #get SDs
  ResultsSD.df = PEP.df.sub %>%
    summarise(MAT = sd(MAT), 
              MAP = sd(MAP), 
              RAD = sd(RAD),
              ELE = sd(alt))%>%
    mutate(year=year)%>% 
    pivot_longer(., cols = -year)%>% 
    rename(sd.climate = value)
  
  #Merge
  Results.df = inner_join(ResultsMean.df, ResultsSD.df, by=c('year','name'))
  
  #Rbind loop subsets
  site.mw.df = rbind(site.mw.df, Results.df)
  
  #count
  print(paste0('year ', year, ' done (', min(year.vector.long),'-',max(year.vector.long),')'))
}

#order variables
site.mw.df = site.mw.df%>%
  mutate(name = factor(name, levels=c("MAT", "MAP", 'RAD',"ELE"), ordered=T))

#create unique site dataframe
PEP.sites.df = PEPlong.df%>%
  filter(!duplicated(pep_id))

#Run linear model (in response to year)
SiteModelResults.df = site.mw.df %>%
  group_by(name)%>%
  do({
    model = lm(mean.climate ~ year, data=.)
    #create combined dataframe
    data.frame(tidy(model))}) %>%  
  ungroup() %>%
  #delete intercept
  filter(!term %in% c("(Intercept)"))%>%  
  mutate(estimate = estimate*10,
         mean.climate = NA,
         sd.climate = NA)%>%  
  #order factors
  mutate(name = factor(name, levels=c("MAT", "MAP", 'RAD',"ELE"), ordered=T))

#create dummy dataset to control y axis
dummy <- data.frame(year = c(1965,2015), 
                    mean.climate = c(mean(PEP.sites.df$MAT)+2*sd(PEP.sites.df$MAT),
                                     mean(PEP.sites.df$MAT)-2*sd(PEP.sites.df$MAT),
                                     mean(PEP.sites.df$MAP)+2*sd(PEP.sites.df$MAP),
                                     mean(PEP.sites.df$MAP)-2*sd(PEP.sites.df$MAP),
                                     mean(PEP.sites.df$RAD)+2*sd(PEP.sites.df$RAD),
                                     mean(PEP.sites.df$RAD)-2*sd(PEP.sites.df$RAD),
                                     mean(PEP.sites.df$alt)+2*sd(PEP.sites.df$alt),
                                     mean(PEP.sites.df$alt)-2*sd(PEP.sites.df$alt)),
                    sd.climate=0,
                    name = rep(c("MAT",'MAP','RAD','ELE'),each=2), stringsAsFactors=FALSE)%>%
  mutate(name = factor(name, levels=c("MAT", "MAP", 'RAD',"ELE"), ordered=T))

#plot
SiteClimatePlot = ggplot(site.mw.df, aes(x = year, y = mean.climate, ymin = mean.climate-sd.climate, ymax = mean.climate+sd.climate,
                                         group=name, color=name, fill=name)) + 
  geom_ribbon(color=NA) + 
  geom_line(size = 1) + 
  scale_color_manual(values = c('#F21A00','#3B9AB2','#E1AF00','#78B7C5'))+
  scale_fill_manual(values = alpha(c('#F21A00','#3B9AB2','#E1AF00','#78B7C5'),0.3))+
  xlab("Year") + ylab('Value')+
  plotTheme1+
  geom_text(data    = SiteModelResults.df,
            mapping = aes(x = Inf, y = Inf, hjust = 1.5, vjust = 2.5, 
                          label = paste('Trend: ', round(estimate,1), " per decade", sep="")), 
            size=3.5, color="black")+
  xlab("") + ylab('')+
  scale_x_continuous(breaks = seq(1966,1996,by=5),
                     labels = c('1966-1985','1971-1990','1976-1995','1981-2000','1986-2005','1991-2010','1996-2015'))+
  coord_cartesian(xlim=c(1967.3,1994.7))+
  facet_wrap(~name,ncol=2,scale='free_y', strip.position="top")+
  geom_blank(data=dummy) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

#Save PDFs
pdf(paste(output_path_figure,"SiteClimateOverTime.pdf",sep="/"), width=6, height=6, useDingbats=FALSE)
SiteClimatePlot
dev.off()

#remove stuff
rm(dummy,PEP.df.sub,PEP.sites.df,SiteClimatePlot,Results.df,ResultsMean.df,ResultsSD.df,site.mw.df,SiteModelResults.df)



#########################################################################################################################
#########################################################################################################################



#######################################
# Monthly and full model correlations #
#######################################



# Run Mixed effects models
##########################

for(k in 1:length(PEPdataList)){

  #Loop through covariate groups
  for (i in 1:length(variables)){

    #create explanatory variables
    #############################

    covariates.monthly  = paste0(variables[i], c(1:10))
    covariates.seasonal = paste0(variables[i], c('.LO.SO','.SO.SE')) 
  
    #---------------------------------------------------------
  
    #set equations
    ##############
  
  
    #################
    # Type 1: monthly
    #################
  
    if(variables[i] %in% c('Apm', 'Azani', 'GSI', 'GSIrad')) {  
    equation1 = as.formula(paste("scale(leaf_off) ~ ", paste('scale(', covariates.monthly[5:8], ')', collapse="+"), 
                                 '+ (1|timeseries) + (1|species)', collapse=""))
    equation1.species = as.formula(paste("scale(leaf_off) ~ ", paste('scale(', covariates.monthly[5:8], ')', collapse="+"), 
                                  '+ (1|timeseries)', collapse="")) 
    } else {
      equation1 = as.formula(paste("scale(leaf_off) ~ ", paste('scale(', covariates.monthly[3:9], ')', collapse="+"), 
                                       '+ (1|timeseries) + (1|species)', collapse=""))
      equation1.species = as.formula(paste("scale(leaf_off) ~ ", paste('scale(', covariates.monthly[3:9], ')', collapse="+"), 
                                       '+ (1|timeseries)', collapse="")) }
  
    #####################
    # Type 2: Full models
    #####################
    
    equation2 = as.formula(paste("scale(leaf_off) ~ ", paste0('scale(',covariates.seasonal[1], ') + scale(',covariates.seasonal[2],')', collapse="+"), 
                                       '+ scale(Prcp.LO.SO) + scale(Prcp.SO.SE) + scale(CO2) + scale(Tnight) + (1|timeseries) + (1|species)', collapse=""))
    equation2.species = as.formula(paste("scale(leaf_off) ~ ", paste0('scale(',covariates.seasonal[1], ') + scale(',covariates.seasonal[2],')', collapse="+"), 
                                 '+ scale(Prcp.LO.SO) + scale(Prcp.SO.SE) + scale(CO2) + scale(Tnight) + (1|timeseries)', collapse=""))
    
    #---------------------------------------------------------
  
    #create moving window dataframes
    mw.df = data.frame()
    mw.species.df = data.frame()
  
    #loop through years (moving windwows)
    for (year in YearList[[k]]){

      #create table subset
      PEP.df.sub = PEPdataList[[k]][PEPdataList[[k]]$year>=year & PEPdataList[[k]]$year< year+MovingWindowLength[k],] %>% 
        group_by(timeseries) %>%  
        filter(if (k==1) {n() >= 15} else {n() >= 12}
               ) %>% #delete time series with less than 15/12 years
        ungroup()
  
      #---------------------------------------------------------
  
      #################################
      #mixed effects models all species
      #################################

      ModelResults.df = PEP.df.sub %>%
       
         do({
      
          #run models
          ###########
      
          #Equation 1
          modelEq1 = lmer(equation1, data=., control = lmerControl(optimizer ="Nelder_Mead"))
          
          #Equation 2
          modelEq2 = lmer(equation2, data=.,control = lmerControl(optimizer ="Nelder_Mead"))
          
      
          #create combined dataframe
          ##########################
      
          data.frame(rbind(
            
            #Equation 1
            tidy(modelEq1, effects="fixed") %>% 
              filter(!term %in% c("(Intercept)")) %>%
              mutate(equation = 'monthly',
                     term = readr::parse_number(term)),
             
            #Equation 2
            tidy(modelEq2, effects="fixed") %>% 
              mutate(equation = 'full model') ) ) 
          
      }) %>% 
        
        #add and edit information
        mutate(species = 'Aall',
               variable = variables[i],
               year     = year) %>%
        filter(!term %in% c("(Intercept)"))
  
      #rbind moving window subsets  
      mw.df = rbind(mw.df, ModelResults.df)
      
      #----------------------------------------------------
  
      ######################################
      #species-specific mixed effects models 
      ######################################
  
      ModelResults.df = PEP.df.sub %>%
        group_by(species)%>%
        do({
          
          #run models
          ###########
      
          #Equation 1
          modelEq1 = lmer(equation1.species, data=.,control = lmerControl(optimizer ="Nelder_Mead"))
          
          #Equation 2
          modelEq2 = lmer(equation2.species, data=.,control = lmerControl(optimizer ="Nelder_Mead"))
      
          
          #create combined dataframe
          ##########################

          data.frame(rbind(
        
            #Equation 1
            tidy(modelEq1, effects="fixed") %>%
              filter(!term %in% c("(Intercept)")) %>%
              mutate(equation = 'monthly',
                     term = readr::parse_number(term)),#rename factors
            
            #Equation 2
            tidy(modelEq2, effects="fixed") %>% 
              filter(!term %in% c("(Intercept)")) %>%
              mutate(equation = 'full model') ) )
          
          }) %>% 
        mutate(variable = variables[i],
               year = year) %>%
        ungroup()
    
      #rbind moving window subsets  
      mw.species.df = rbind(mw.species.df, ModelResults.df)
  
      print(paste0('year ', year, ' done (', min(PEPdataList[[k]]$year),'-',(max(PEPdataList[[k]]$year)-19),')'))
      }
  
    #rbind all species and species-specific results
    mw.df = bind_rows(mw.df, mw.species.df) %>% 
      #rename factors
      mutate(dataset = ifelse(k==1, 'Long', 'Short'),
             term = gsub("scale","",term),
             term = gsub("\\(|\\)","",term)) #remove brackets

     #---------------------------------------------------------
  
    #store dataframe in variable list
    DataList[[i+(k-1)*length(variables)]] = mw.df
  
    #count
    print(paste0(i+(k-1)*length(variables),' out of ',length(DataList), ' (',variables[i],') done'))
    }
  }

#bind rows
MovingWindowAnalysis.df = bind_rows(DataList) 

#Safe table
write.csv(MovingWindowAnalysis.df, paste(output_path, "Moving_window_data.csv", sep="/"))



##############################################################################################################################################
##############################################################################################################################################



#####################
## Reproducibility ##	
#####################


## date time
Sys.time()
#"2021-07-10 08:42:19 CEST"


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
#  [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] broom.mixed_0.2.6 patchwork_1.1.1   wesanderson_0.3.6 car_3.0-11        carData_3.0-4     lme4_1.1-27.1    
#[7] Matrix_1.3-3      raster_3.4-13     sp_1.4-5          pbmcapply_1.5.0   sjmisc_2.8.7      gmodels_2.18.1   
#[13] broom_0.7.8       data.table_1.14.0 forcats_0.5.1     stringr_1.4.0     dplyr_1.0.7       purrr_0.3.4      
#[19] readr_1.4.0       tidyr_1.1.3       tibble_3.1.2      ggplot2_3.3.4     tidyverse_1.3.1  

#loaded via a namespace (and not attached):
#  [1] nlme_3.1-152      fs_1.5.0          lubridate_1.7.10  insight_0.14.2    httr_1.4.2        TMB_1.7.20       
#[7] tools_4.1.0       backports_1.2.1   utf8_1.2.1        R6_2.5.0          sjlabelled_1.1.8  DBI_1.1.1        
#[13] colorspace_2.0-1  withr_2.4.2       tidyselect_1.1.1  curl_4.3.2        compiler_4.1.0    cli_2.5.0        
#[19] rvest_1.0.0       xml2_1.3.2        labeling_0.4.2    scales_1.1.1      digest_0.6.27     foreign_0.8-81   
#[25] minqa_1.2.4       rmarkdown_2.9     rio_0.5.27        pkgconfig_2.0.3   htmltools_0.5.1.1 maps_3.3.0       
#[31] dbplyr_2.1.1      rlang_0.4.11      readxl_1.3.1      rstudioapi_0.13   farver_2.1.0      generics_0.1.0   
#[37] jsonlite_1.7.2    gtools_3.9.2      zip_2.2.0         magrittr_2.0.1    Rcpp_1.0.6        munsell_0.5.0    
#[43] fansi_0.5.0       abind_1.4-5       lifecycle_1.0.0   stringi_1.6.2     yaml_2.2.1        MASS_7.3-54      
#[49] plyr_1.8.6        grid_4.1.0        gdata_2.18.0      crayon_1.4.1      lattice_0.20-44   haven_2.4.1      
#[55] splines_4.1.0     hms_1.1.0         knitr_1.33        pillar_1.6.1      boot_1.3-28       reshape2_1.4.4   
#[61] codetools_0.2-18  reprex_2.0.0      glue_1.4.2        evaluate_0.14     modelr_0.1.8      vctrs_0.3.8      
#[67] nloptr_1.2.2.2    cellranger_1.1.0  gtable_0.3.0      assertthat_0.2.1  xfun_0.24         openxlsx_4.2.4   
#[73] coda_0.19-4       ellipsis_0.3.2   



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################



