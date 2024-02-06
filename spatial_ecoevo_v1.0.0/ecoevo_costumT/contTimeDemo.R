setwd("~/EcoEvolution/spatial_ecoevo_v1.0.0/ecoevo_costumT")
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
#--- large time plot----
folder <-"/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/ecoevo_costumT/parameters/"
file <- "v3.225165e-05_d0.001idKozaiPreciseDesign2_ForPaper"
load(paste(folder,file,sep = ""))
large_dat <-suppressMessages(read_csv(outfile))
plot1 <- plot_landscape(large_dat %>% filter(patch %in% c(1,11,20)))+ggtitle("High Genetic Variance")
#--- small time plot----
# if small has failed (no way of knowing before reading the correct file)
file <- "v1.612582e-05_d0.001idKozaiPreciseDesign2_ForPaper_FAILED"
load(paste(folder,file,sep = ""))
small_dat <-read_csv(outfile)
plot2 <- plot_landscape(small_dat %>% filter(patch %in% c(1,11,20)))+ggtitle("Low Genetic Variance")

#---Combined plot----
plotCombined <- grid.arrange(plot1, plot2, nrow=2)
print(plotCombined)
png_fold <- "/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/Plots/"
ggsave(filename =  paste(png_fold,"KozaiPreciseDesign_contTime_Paper.png",sep="")
       , plot = plotCombined, width = 8, height = 5, units = "in")