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

times <- (c(0,5e6,5e7,1e8)+5e6*0)
#--- large time plot----
# folder <-"/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/ecoevo_costumT/parameters/"
# file <- "v3.225165e-05_d0.001idKozaiPreciseDesign2_ForPaper"
# load(paste(folder,file,sep = ""))
folder <- "/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/outputs/"
file <- "large_time_v3e-05_d1e-07idPeriodicCC11093_FullSin_23_11"
#file <- "large_time_v3e-05_d1e-07idBiggerCompetition_Seed7234"
outfile <- paste(folder,file,sep="")
large_dat <-suppressMessages(read_csv(outfile))
plot1 <- plot_landscape(large_dat %>% filter(patch %in% c(1,11,20)))+ggtitle("Long Adaptaion Time")

#--- small time plot----
file <- "small_time_v3e-05_d1e-07idPeriodicCC11093_FullSin_23_11_FAILED"
#file <- "small_time_v3e-05_d1e-07idBiggerCompetition_Seed7234_FAILED"
outfile <- paste(folder,file,sep="")
small_dat <-read_csv(outfile)
plot2 <- plot_landscape(small_dat %>% filter(patch %in% c(1,11,20)))+ggtitle("Short Adaptaion Time")

# plot1 <- plot_traitLag(large_dat%>%filter(time %in% times),nmin)+ggtitle("High Genetic Variance")
# plot2 <- plot_timeseries(large_dat%>%filter(time %in% times)%>%filter(n>nmin))
# plot3 <- plot_traitLag(small_dat%>%filter(time %in% times),nmin)+ggtitle("Low Genetic Variance")
# plot4 <- plot_timeseries(small_dat%>%filter(time %in% times)%>%filter(n>nmin))
#---Combined plot----
plotCombined <- grid.arrange(plot1,plot2, nrow=2)
print(plotCombined)
png_fold <-"/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/Plots/"
ggsave(filename =  paste(png_fold,"PeriodicCC_contTime_Paper.png",sep="")
         , plot = plotCombined, width = 16, height = 10, units = "in")