setwd("~/EcoEvolution/spatial_ecoevo_v1.0.0/ecoevo_original")

suppressPackageStartupMessages({
  suppressWarnings({
    library(gridExtra)
    require(ggpmisc) # adding statistics to plots
    require(Rcpp) # importing C functions
    library(ggplot2)
    library(dplyr)
    library(readr)
    source("./plotting_functions.R") # various functions for plotting final data
  })
})

# param_folder <-"/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/Parameters/"
# name <-"large_time_v3e-05_d1e-07idPeriodicCC11093_FullSin_23_11"
# #name <- "small_time_v3e-05_d1e-07idPeriodicCC11093_FullSin_23_11_FAILED"
# file <- paste(param_folder,name,sep = "")
# load(file)
outfolder <-"/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/outputs/"
file <-"v3.225165e-05_d0.001idKozaiPreciseDesign2_ForPaper"
dat<-read_csv(paste(outfolder,file,sep = ""))
plotCombined <- grid.arrange(p1, p2, nrow=2)

