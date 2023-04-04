## This repository contains all R code included in the publication "Effect of climate warming on the timing of autumn leaf senescence reverses after the summer solstice" (Zohner et al. 2023, Science)

### Workflow

1. Experiment 1
    - Seasonal temperature and light manipulation experiment
    - Run `Experiment1_chlorophyll.Rmd` for the leaf senescence analysis
    - Run `Experiment1_photosynthesis.Rmd` for the photosynthesis analysis
2. Experiment 2
    - The effects of pre- and post-solstice warming on bud set timing
    - Run `Experiment2_budset.Rmd` for the bud set analysis 
3. FluxNet analysis
    - Analyse autumn photosynthetic declines using flux tower measurements 
    - Run `FluxNet_analysis.Rmd` for the analysis and figures
4. Harvard analysis
    - Ground-sourced phenology observations for Eastern North America (Harvard forest)
    - Run `Harvard_models.Rmd` to create models and figures
5. PEP725 analysis
    - Ground-sourced phenology observations for Europe (PEP725 data)
    - Step 1: `1_Data_extraction` - download `www.PEP725.eu` data, GLDAS climate data and other info
    - Step 2: `2_Add_drivers` - generate seasonal climate and photosynthesis drivers
    - Step 3: `3_Analysis` - run models and generate figures
6. Remote sensing analysis
    - Analysis of satellite-derived phenology observations
    - Step 1: `1_Data_extraction` - extract photoperiod info
    - Step 2: `2_Add_drivers` - generate seasonal climate and photosynthesis drivers
        - EOS10 = Date when EVI last crossed 90% threshold
        - EOSstart = Date of maximum EVI curvature change rate (breakpoint)
        - EOS50 = Date when EVI last crossed 50% threshold
        - EOS85 = Date when EVI last crossed 15% threshold
    - Step 3: `3_Analysis` - run models and generate figures
        - `3.1` Check sample sizes, `3.2` add preseason temperatures, `3.3` run models and `3.4` create figures
        - EOS10 folder = `EOS10_SenescenceStart`
        - EOSstart folder = `EOSstart_VNP`
        - EOS50 foler = `EOS50_MidGreendown`
        - EOS85 folder = `EOS85_Dormancy`
