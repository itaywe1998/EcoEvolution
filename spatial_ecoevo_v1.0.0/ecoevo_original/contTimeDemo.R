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


clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs)>0) { # command-line arguments
  dual <- as.logical(clargs[1]) # if dual display
  id <- clargs[2] # current run name
  vbar <-as.numeric(clargs[3])
  dbar <-as.numeric(clargs[4])
  cycles<-as.numeric(clargs[5])
} else { # sample input parameters, if no command line arguments are given
  id <-"PeriodicCC11093_FullSin_23_11"
  vbar <-3e-5
  dbar <-1e-7
}
# Assuming always dual display
#--- large time plot----
large_dat <-suppressMessages(read_csv(paste("outputs/large_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep ="")))
plot1 <- plot_landscape(large_dat %>% filter(patch %in% c(1,11,20)))+ggtitle("Long Adaptation Period")
#--- small time plot----
suppressMessages(small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),"_FAILED",sep ="")))
plot2 <- plot_landscape(small_dat %>% filter(patch %in% c(1,11,20)))+ggtitle("Short Adaptation Period")

#---Combined plot----
plotCombined <- grid.arrange(plot1, plot2, nrow=2)
print(plotCombined)
png_fold <- "/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/Plots/"
# ggsave(filename =  paste(png_fold,"v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),"contTime_Paper.png",sep ="")
#        , plot = plotCombined, width = 20, height = 8, units = "in")