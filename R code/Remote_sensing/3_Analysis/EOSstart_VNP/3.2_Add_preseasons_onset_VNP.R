


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Run autumn temperature (preseason) models for the satellite data (EOSstart) ###############################
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
Drivers_path     = "Analysis_input/Drivers_final_onset_VNP/Merged_file"
Land_cover_path  = "Analysis_input/Drivers"
output_path      = "Analysis_output_startSen_VNP/Preseason_temperatures"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



#Phenology data
Pheno.df <- fread(paste(Drivers_path, "Remote_sensing_drivers_data_onset_VNP.csv", sep="/"))%>% 
  dplyr::select(-V1) %>%
  group_by(geometry) %>%
  #delete pixels with less than 9 years
  filter(n() >= 9) %>%
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
  dplyr::select(geometry,Year,LC_Type,Onset_Greenness_Decrease,
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
    
    model = lm(scale(Onset_Greenness_Decrease) ~ scale(temp), data=.)  # linear model
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
  coord_cartesian(ylim = c(0.01, 0.25))+
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
  plotTheme1 +
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
write.csv(Pheno.df, paste(Drivers_path, "Remote_sensing_drivers_data_onset_VNP_preseason.csv", sep="/"))



##############################################################################################################################################
##############################################################################################################################################



##################################
#Run linear ridge regression model
##################################



resultsLM3 = Pheno.df %>% 
  group_by(geometry,LC_Type) %>% 
  do({model = lm.ridge(scale(Onset_Greenness_Decrease) ~ scale(Tday)+scale(Tnight), data=.)  # linear model
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
pdf(paste(output_path,"Preseason_sensitivity_startSen_VNP.pdf",sep="/"), width=8, height=7, useDingbats=FALSE)
PreseasonPlot
dev.off()



##############################################################################################################################################
##############################################################################################################################################



#####################
## Reproducibility ##	
#####################


## datetime
Sys.time()
#"2023-01-26 08:50:51 CET"


## session info
sessionInfo()
#R version 4.1.0 (2021-05-18)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS 12.5.1

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] MASS_7.3-54          patchwork_1.1.1      gmodels_2.18.1       broom_0.7.8          viridis_0.6.1       
#[6] viridisLite_0.4.0    rgdal_1.5-27         weathermetrics_1.2.2 lubridate_1.7.10     chillR_0.72.4       
#[11] zoo_1.8-9            pbmcapply_1.5.0      rpmodel_1.2.0        forcats_0.5.1        stringr_1.4.0       
#[16] dplyr_1.0.10         purrr_0.3.4          readr_1.4.0          tidyr_1.2.0          tibble_3.1.8        
#[21] ggplot2_3.3.6        tidyverse_1.3.1      raster_3.5-15        sp_1.4-6             ncdf4_1.19          
#[26] sf_1.0-6             data.table_1.14.0   

#loaded via a namespace (and not attached):
#  [1] minqa_1.2.4         colorspace_2.0-3    ellipsis_0.3.2      class_7.3-19        htmlTable_2.2.1    
#[6] base64enc_0.1-3     pls_2.8-0           fs_1.5.2            rstudioapi_0.13     proxy_0.4-26       
#[11] farver_2.1.1        fansi_1.0.3         xml2_1.3.3          codetools_0.2-18    splines_4.1.0      
#[16] R.methodsS3_1.8.1   knitr_1.33          Formula_1.2-4       spam_2.7-0          jsonlite_1.8.0     
#[21] nloptr_2.0.3        cluster_2.1.2       dbplyr_2.1.1        png_0.1-7           R.oo_1.24.0        
#[26] Kendall_2.2         compiler_4.1.0      httr_1.4.2          backports_1.2.1     assertthat_0.2.1   
#[31] Matrix_1.3-3        fastmap_1.1.0       cli_3.3.0           htmltools_0.5.2     tools_4.1.0        
#[36] dotCall64_1.0-1     gtable_0.3.0        glue_1.6.2          maps_3.3.0          Rcpp_1.0.9         
#[41] cellranger_1.1.0    vctrs_0.4.1         gdata_2.18.0        nlme_3.1-152        xfun_0.24          
#[46] lme4_1.1-30         rvest_1.0.2         lifecycle_1.0.1     gtools_3.9.2        XML_3.99-0.8       
#[51] terra_1.5-17        scales_1.2.0        hms_1.1.0           RColorBrewer_1.1-3  fields_12.5        
#[56] yaml_2.2.2          gridExtra_2.3       rpart_4.1-15        latticeExtra_0.6-29 stringi_1.7.6      
#[61] checkmate_2.0.0     e1071_1.7-9         GenSA_1.1.7         boot_1.3-28         rlang_1.0.4        
#[66] pkgconfig_2.0.3     bitops_1.0-7        pracma_2.3.3        evaluate_0.15       lattice_0.20-44    
#[71] labeling_0.4.2      htmlwidgets_1.5.3   tidyselect_1.1.2    remef_1.0.7         plyr_1.8.6         
#[76] magrittr_2.0.3      R6_2.5.1            Hmisc_4.5-0         generics_0.1.3      DBI_1.1.2          
#[81] foreign_0.8-81      pillar_1.8.0        haven_2.4.1         withr_2.5.0         units_0.8-0        
#[86] nnet_7.3-16         survival_3.2-11     RCurl_1.98-1.5      modelr_0.1.8        crayon_1.5.1       
#[91] KernSmooth_2.23-20  utf8_1.2.2          rmarkdown_2.9       jpeg_0.1-8.1        grid_4.1.0         
#[96] readxl_1.3.1        reprex_2.0.0        digest_0.6.29       classInt_0.4-3      R.utils_2.11.0     
#[101] munsell_0.5.0  



##############################################################################################################################################
#############################################################THE END##########################################################################
##############################################################################################################################################


