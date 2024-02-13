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
#file <- "v3.225165e-05_d0.01idKozaiPreciseDesign2_ForPaper"
load(paste(folder,file,sep = ""))
large_dat <-suppressMessages(read_csv(outfile))
tstart <- min(large_dat$time)
final_time <-max(large_dat$time)
during_step <- final_time/200
before_step <- 0#-tstart/1000
ctimes <- c(min_times[c(TRUE, FALSE)],max_times[c(TRUE, FALSE)])
plot1 <- plot_timeseries(large_dat %>% filter(time %in% c(0,ctimes[c(TRUE, FALSE)])))+ggtitle("High Genetic Variance")
#--- small time plot----
# if small has failed (no way of knowing before reading the correct file)
file <- "v1.612582e-05_d0.001idKozaiPreciseDesign2_ForPaper_FAILED"
#file <- "v3.225165e-05_d1000idKozaiPreciseDesign2_ForPaper_FAILED"
load(paste(folder,file,sep = ""))
small_dat <-read_csv(outfile)
tstart <- min(small_dat$time)
final_time<- max(small_dat$time)
before_step <- 0#-tstart/1000
plot2 <- plot_timeseries(small_dat %>% filter(time %in% seq(from=0,to=final_time,by=final_time/2)))+ggtitle("Low Genetic Variance")

#---Combined plot----
plotCombined <- grid.arrange(plot1, plot2, nrow=2)
print(plotCombined)
png_fold <- "/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/Plots/"


s_dat <- small_dat %>% filter(n > 1e-5)
l_dat <- large_dat %>% filter(n > 1e-5)
# ggsave(filename =  paste(png_fold,"KozaiPreciseDesign_Paper.png",sep="")
#        , plot = plotCombined, width = 15, height = 8, units = "in")