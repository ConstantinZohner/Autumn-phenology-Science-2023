


#############################################################################################################
############################################## R script for: ################################################
#############################################################################################################
##### Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice ##
#############################################################################################################


#############################################################################################################
# Full models excluding CO2 for the satellite data (EOS10) ##################################################
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
Drivers_path    = "Analysis_input/Drivers_final_EOS10/Merged_file"
output_path     = "Analysis_output_startSen/Reviewer5"



##############################################################################################################################################
##############################################################################################################################################



#################
## Import data ##
#################



#Phenology data frame
#####################

Pheno.df <- fread(paste(Drivers_path, "Remote_sensing_drivers_data_startSen_preseason.csv", sep="/")) 



##############################################################################################################################################
##############################################################################################################################################



#################
## Full models ##
#################



#Define variables
variables = c('GPPstart', 'Tday') 

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
      equation.scaled1 = as.formula(paste("scale(Senesc_DOY) ~ ", paste0('scale(',covariates[1], ') + scale(', covariates[2], ')', 
                                                                         collapse="+"), 
                                          '+ scale(Prcp.LO.SO) + scale(Prcp.SO.SE)', collapse=""))
      
      equation.2 = as.formula(paste("Senesc_DOY ~ ", paste0(covariates[1], '+', covariates[2], 
                                                                         collapse="+"), 
                                          '+ Prcp.LO.SO + Prcp.SO.SE', collapse=""))
      

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
          model2 = lm(equation.2,  data=.)
          
          #create combined data frame
          ###########################
          
          data.frame(rbind(
            
            #Equation 1
            tidy(model1) %>% mutate(equation = 'full model 1'),
            tidy(model2) %>% mutate(equation = 'full model 2'))
          
          )
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



###########
## Plots ##
###########



variable.name="Tday"
plotTday = FullAnalysis.df %>%
  filter(variable == variable.name,
         equation=="full model 1") %>%
  mutate(term = factor(term, 
                       levels=c(paste0(variable.name, ".LO.SO"),
                                "Prcp.LO.SO", 'Prcp.SO.SE', 
                                paste0(variable.name, ".SO.SE")), ordered=T) ) %>% 
  ggplot(aes(x = term, y = estimate, fill=term)) + 
  geom_boxplot(outlier.shape = NA, notch=T)+
  geom_hline(yintercept=0)+
  xlab("") + ylab("Standardized effect") +
  coord_cartesian(ylim = c(-.9,.9)) +
  scale_fill_manual(values = c('#F21A00','grey60','grey35','#3B9AB2'))+
  scale_x_discrete(labels = c(paste0(variable.name, " pre"),
                              'Prcp pre','Prcp post',
                              paste0(variable.name, " post")))+
  plotTheme1 +
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        legend.position = "none")

variable.name="GPPstart"
plotGPP = FullAnalysis.df %>%
  filter(variable == variable.name,
         equation=="full model 1") %>%
  mutate(term = factor(term, 
                       levels=c(paste0(variable.name, ".LO.SO"),
                                "Prcp.LO.SO", 'Prcp.SO.SE', 
                                paste0(variable.name, ".SO.SE")), ordered=T) ) %>% 
  ggplot(aes(x = term, y = estimate, fill=term)) + 
  geom_boxplot(outlier.shape = NA, notch=T)+
  geom_hline(yintercept=0)+
  xlab("") + ylab("") +
  coord_cartesian(ylim = c(-.9,.9)) +
  scale_fill_manual(values = c('#F21A00','grey60','grey35','#3B9AB2'))+
  scale_x_discrete(labels = c(paste0(variable.name, " pre"),
                              'Prcp pre','Prcp post',
                              paste0(variable.name, " post")))+
  plotTheme1 +
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        legend.position = "none")



##############################################################################################################################################
##############################################################################################################################################



################
## Safe plots ##
################



layout <- "AB"

#Merge plots
Fig_Plot = plotTday + plotGPP +
  plot_layout(design = layout) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

#save plots as .pdf
ggsave(Fig_Plot, file='Fig2B_noCo2.pdf', 
       path=output_path,
       width=6, height=3.5)



##############################################################################################################################################
##############################################################################################################################################


