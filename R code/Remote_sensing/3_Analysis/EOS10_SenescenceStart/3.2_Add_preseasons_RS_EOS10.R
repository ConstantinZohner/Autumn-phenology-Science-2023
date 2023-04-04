


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Run autumn temperature (preseason) models for the satellite data (EOS10) ##################################
#############################################################################################################



#required packages
require(tidyverse)
require(data.table)
require(broom)
require(gmodels)
require(patchwork)
require(MASS)


#define plot themes
plotTheme1 = theme(legend.position   = "none",
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



##########################################
## Set directories and get own PEP data ##
##########################################



# set the working directory
setwd("/Users/consti/Desktop/PhD/Publication_material/17_Autumn_phenology_tier2/Remote_sensing/Analysis")

# paths
Drivers_path     = "Analysis_input/Drivers_final_startSen/Merged_file"
Land_cover_path  = "Analysis_input/Drivers"
output_path      = "Analysis_output_startSen/Preseason_temperatures"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



#Phenology data
Pheno.df <- fread(paste(Drivers_path, "Remote_sensing_drivers_data_startSen.csv", sep="/"))%>% 
  dplyr::select(-V1) %>%
  filter(
    #delete senescence dates before DOY 140 and after DOY 290
    Senesc_DOY>140,Senesc_DOY<290,
    #delete observation if senescence date occurs before MidGreenup date
    Senesc_DOY>MidGreenup_DOY)%>%
  group_by(geometry) %>%
  #delete pixels with less than 15 years
  filter(n() >= 15) %>%
  ungroup()

#Land cover info
LandCover.df <- fread(paste(Land_cover_path, "Land_Cover_Type_025.csv", sep="/")) %>%
  mutate(geometry = gsub("POINT ","",  geometry),
         geometry = gsub("\\(|\\)","", geometry)) %>%
  separate(geometry, into = c("Lon","Lat"), sep=" ") %>%
  mutate(Lat = round(as.numeric(Lat),3),
         Lon = round(as.numeric(Lon),3) ) 
          
#Merge tables
Pheno.df = Pheno.df %>%
  mutate(Year = as.numeric(Year)) %>%
  inner_join(LandCover.df,  by=c("Lat","Lon")) %>%
  #delete evergreen broadleaf pixels
  filter(!LC_Type == "EvgB")

rm(LandCover.df)



##############################################################################################################################################
##############################################################################################################################################



########################
## Get best preseason ##
########################



#reshape table to long format
#############################

preseason.df = Pheno.df %>%
  #select columns
  dplyr::select(geometry,Year,LC_Type,Senesc_DOY,
         Tday.PS.10,Tday.PS.20,Tday.PS.30,
         Tday.PS.40,Tday.PS.50,Tday.PS.60,
         Tday.PS.70,Tday.PS.80,Tday.PS.90,
         Tday.PS.100,Tday.PS.110,Tday.PS.120,
         Tnight.PS.10,Tnight.PS.20,Tnight.PS.30,
         Tnight.PS.40,Tnight.PS.50,Tnight.PS.60,
         Tnight.PS.70,Tnight.PS.80,Tnight.PS.90,
         Tnight.PS.100,Tnight.PS.110,Tnight.PS.120)%>%
  #long format
  pivot_longer(.,cols=starts_with(c("Td","Tn")), names_to = "preseason", values_to = "temp") %>%
  #create preseason length and temperature class columns
  mutate(preseason_length = readr::parse_number(gsub("[.]", "", preseason)), #keep only numbers in string
         temp_class = gsub("\\..*","", preseason)) %>%
  dplyr::select(-preseason)


#Run linear models
##################

resultsLM = preseason.df %>% 
  group_by(geometry, LC_Type, temp_class, preseason_length) %>% 
  do({
    
    model = lm(scale(Senesc_DOY) ~ scale(temp), data=.)  # linear model
    data.frame(tidy(model), glance(model) )}) %>% # model info
  
  filter(!term %in% c("(Intercept)")) %>%
  dplyr::select(geometry,LC_Type,temp_class,preseason_length,estimate,r.squared)%>%
  ungroup()



##############################################################################################################################################
##############################################################################################################################################



############################################
## Plot preseason-senescence correlations ##
############################################



#R2
###

resultsLM = resultsLM %>%
  mutate(LC_Type = factor(LC_Type, levels=c("Mixed","DecB","EvgN","DecN"), ordered=T)) 

plot.R2 = resultsLM %>%
  
  ggplot()+
  aes(x=preseason_length, y=r.squared, 
      colour=temp_class) +
  
  stat_summary(fun=mean, geom="line", size = .7) +
  scale_colour_manual(values = c('#F21A00', '#3B9AB2')) +
  
  stat_summary(fun.data = "mean_cl_normal", geom="errorbar", size = 0.5, width=0) +
  stat_summary(fun.data = "mean_cl_normal", geom="point", size = 1) +
  scale_fill_manual(values = c('#F21A00', '#3B9AB2'))+
  
  xlab("Preseason length (days)") +
  ylab("Coefficient of determination (R2)") +
  coord_cartesian(ylim = c(0.01, 0.15))+
  facet_wrap(~LC_Type, ncol=1) +
  plotTheme1 +
  theme(strip.text.x = element_blank())


#Correlation coefficient
########################

plot.estimate = resultsLM %>%
  
  ggplot()+
  aes(x=preseason_length, y=estimate, 
      colour=temp_class) +
  
  geom_hline(yintercept = 0)+
  stat_summary(fun = mean, geom="line", size = .7) +
  scale_colour_manual(values = c('#F21A00', '#3B9AB2')) +
  
  stat_summary(fun.data = "mean_cl_normal", geom="errorbar", size = 0.5, width=0) +
  stat_summary(fun.data = "mean_cl_normal", geom="point", size = 1) +
  scale_fill_manual(values = c('#F21A00', '#3B9AB2'))+
  
  xlab("Preseason length (days)") +
  ylab("Standardized coefficient") +
  #coord_cartesian(ylim = c(0.01, 0.28))+
  facet_wrap(~LC_Type, ncol=1) +
  plotTheme1+
  theme(strip.text.x = element_blank())




##############################################################################################################################################
##############################################################################################################################################



#####################################################
## Plot best preseason length for each temperature ##
#####################################################



#keep only models with best predictions
resultsLM2 = resultsLM %>%   
  group_by(geometry,temp_class) %>% 
  top_n(1, r.squared) %>%
  ungroup()

#plot
plot.length = resultsLM2 %>%
  mutate(LC_Type      = factor(LC_Type, levels=c("Mixed","DecB","EvgN","DecN"), ordered=T)) %>%
  dplyr::select(LC_Type,temp_class,preseason_length)%>%
  
  ggplot() + aes(x=temp_class, y=preseason_length) +
  
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x), 
               fun.max = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               size=0.5,
               aes(colour   = temp_class)) +
  
  scale_colour_manual(values = c('#F21A00', '#3B9AB2')) +
  
  coord_cartesian(ylim = c(5, 120))+
  xlab("Daily temperature") +
  ylab("Best preseason length (days)") +
  facet_wrap(~LC_Type, ncol=1) +
  plotTheme1+
  theme(strip.text.x     = element_blank(),
        axis.text.x      = element_text(angle = 45, hjust=1))



##############################################################################################################################################
##############################################################################################################################################



#####################################
#Add best preseason temps to PEP data
#####################################



Pheno.df = Pheno.df %>%
  inner_join(., preseason.df %>%
               #filter by model data
               semi_join(resultsLM2, by=c('geometry','temp_class','preseason_length')) %>%
               dplyr::select(c(geometry,Year,temp_class,temp))%>%
               pivot_wider(.,names_from = temp_class, values_from = temp),
             by = c("Year", "geometry"))%>%
  dplyr::select(-(cols=starts_with(c("Tday.PS","Tnight.PS"))))

#Safe table
write.csv(Pheno.df, paste(Drivers_path, "Remote_sensing_drivers_data_startSen_preseason.csv", sep="/"))



##############################################################################################################################################
##############################################################################################################################################



##################################
#Run linear ridge regression model
##################################



resultsLM3 = Pheno.df %>% 
  group_by(geometry,LC_Type) %>% 
  do({model = lm.ridge(scale(Senesc_DOY) ~ scale(Tday)+scale(Tnight), data=.)  # linear model
  data.frame(tidy(model),         # coefficient info
             glance(model))}) %>% # model info
  ungroup() %>%
  #rename temperature class
  mutate(term=dplyr::recode(term, `scale(Tday)`="Tday", `scale(Tnight)`="Tnight"))

#plot preseason-senescence correlations
plot.ridge = resultsLM3 %>%
  mutate(LC_Type = factor(LC_Type, levels=c("Mixed","DecB","EvgN","DecN"), ordered=T) ) %>%
  
  ggplot()+
  aes(x=term, y=estimate, 
      colour=term, fill = term) +
  scale_colour_manual(values = c('#F21A00','#3B9AB2')) +
  stat_summary(fun.data = "mean_cl_normal", geom="errorbar", size = 0.9, width=0) +
  stat_summary(fun.data = "mean_cl_normal", geom="point", size = 1) +
  
  geom_hline(yintercept = 0)+
  xlab("Daily temperature") +
  ylab("Standardized coefficient (ridge regression)") +
  coord_cartesian(ylim = c(-0.3, 0.3))+
  facet_wrap(~LC_Type, ncol=1,strip.position = "right") +
  plotTheme1+
  theme(axis.text.x = element_text(angle = 45, hjust=1))



##############################################################################################################################################
##############################################################################################################################################



##########################
# Arrange and safe plots #
##########################



#define plot layout
layout <- "AABBCD"

#Merge plots
PreseasonPlot = plot.R2 + plot.estimate + plot.length + plot.ridge +
  plot_layout(design = layout) + plot_annotation(tag_levels = 'a')&
  theme(plot.tag = element_text(face = 'bold'))

#Safe plot
pdf(paste(output_path,"Preseason_sensitivity_RS_startSen.pdf",sep="/"), width=8, height=7, useDingbats=FALSE)
PreseasonPlot
dev.off()



##############################################################################################################################################
##############################################################################################################################################



#####################
## Reproducibility ##	
#####################


## datetime
Sys.time()
#"2021-12-05 10:24:06 CET"


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
#  [1] MASS_7.3-54       patchwork_1.1.1   gmodels_2.18.1    broom_0.7.8       data.table_1.14.0 forcats_0.5.1    
#[7] stringr_1.4.0     dplyr_1.0.7       purrr_0.3.4       readr_1.4.0       tidyr_1.1.3       tibble_3.1.2     
#[13] ggplot2_3.3.4     tidyverse_1.3.1  

#loaded via a namespace (and not attached):
#  [1] httr_1.4.2          jsonlite_1.7.2      splines_4.1.0       modelr_0.1.8        gtools_3.9.2        Formula_1.2-4      
#[7] assertthat_0.2.1    latticeExtra_0.6-29 cellranger_1.1.0    pillar_1.6.1        backports_1.2.1     lattice_0.20-44    
#[13] glue_1.4.2          digest_0.6.27       RColorBrewer_1.1-2  checkmate_2.0.0     rvest_1.0.0         colorspace_2.0-1   
#[19] htmltools_0.5.1.1   Matrix_1.3-3        pkgconfig_2.0.3     haven_2.4.1         scales_1.1.1        gdata_2.18.0       
#[25] jpeg_0.1-8.1        htmlTable_2.2.1     generics_0.1.0      farver_2.1.0        ellipsis_0.3.2      withr_2.4.2        
#[31] nnet_7.3-16         cli_2.5.0           survival_3.2-11     magrittr_2.0.1      crayon_1.4.1        readxl_1.3.1       
#[37] fs_1.5.0            fansi_0.5.0         xml2_1.3.2          foreign_0.8-81      tools_4.1.0         hms_1.1.0          
#[43] lifecycle_1.0.0     munsell_0.5.0       reprex_2.0.0        cluster_2.1.2       compiler_4.1.0      rlang_0.4.11       
#[49] grid_4.1.0          rstudioapi_0.13     htmlwidgets_1.5.3   base64enc_0.1-3     labeling_0.4.2      gtable_0.3.0       
#[55] DBI_1.1.1           R6_2.5.0            gridExtra_2.3       lubridate_1.7.10    knitr_1.33          utf8_1.2.1         
#[61] Hmisc_4.5-0         stringi_1.6.2       Rcpp_1.0.6          vctrs_0.3.8         rpart_4.1-15        png_0.1-7          
#[67] dbplyr_2.1.1        tidyselect_1.1.1    xfun_0.24  



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


