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
folder <-"/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/Parameters/"
#file <- "large_time_v3e-05_d1e-07idPeriodicCC11093_FullSin_23_11"
file <- "v3.225165e-05_d0.001idKozaiPreciseDesign2_ForPaper"
load(paste(folder,file,sep = ""))
ofolder<-"/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/outputs/"
ofile <-paste(ofolder,file,sep="")
large_dat <-suppressMessages(read_csv(ofile))
tstart <- min(large_dat$time)
final_time <-max(large_dat$time)
during_step <- final_time/200
before_step <- 0#-tstart/1000

#ctimes <- c(tstart,seq(0,tE,tE/10))
ctimes<-c(0,5.5e7,1.35e8,2.65e8,3.45e8,4.8e8,5.6e8,6.9e8,7.7e8,9.05e8,9.85e8)
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
plot1 <- plot_timeseries(large_dat %>% filter(time %in% ctimes))+ggtitle("High Genetic Variance")
#--- small time plot----
# if small has failed (no way of knowing before reading the correct file)
#file <- "small_time_v3e-05_d1e-07idPeriodicCC11093_FullSin_23_11_FAILED"
#file <- "v3.225165e-05_d1000idKozaiPreciseDesign2_ForPaper_FAILED"
file <- "v1.612582e-05_d0.001idKozaiPreciseDesign2_ForPaper_FAILED"
load(paste(folder,file,sep = ""))
#file <- "small_time_v3e-05_d1e-07idPeriodicCC11093_FullSin_23_11_FAILED"
file <- "v1.612582e-05_d0.001idKozaiPreciseDesign2_ForPaper_FAILED"
ofile <-paste(ofolder,file,sep="")
small_dat <-read_csv(ofile)
tstart <- min(small_dat$time)
final_time<- max(small_dat$time)
before_step <- 0#-tstart/1000
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
ctimes <- c(tstart,seq(0,final_time,final_time/2))
plot2 <- plot_timeseries(small_dat %>% filter(time %in% ctimes))+ggtitle("Low Genetic Variance")

#---Combined plot----
plotCombined <- grid.arrange(plot1, plot2, nrow=2)
print(plotCombined)
png_fold <- "/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/Plots/"


s_dat <- small_dat %>% filter(n > 1e-5)
l_dat <- large_dat %>% filter(n > 1e-5)
# ggsave(filename =  paste(png_fold,"KozaiPreciseDesign_Paper.png",sep="")
#        , plot = plotCombined, width = 15, height = 8, units = "in")