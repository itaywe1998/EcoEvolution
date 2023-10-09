
library(gridExtra)
require(ggpmisc) # adding statistics to plots
require(Rcpp) # importing C functions
library(ggplot2)
library(dplyr)
library(readr)
source("./plotting_functions.R") # various functions for plotting final data

clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs)>0) { # command-line arguments
  dual <- as.logical(clargs[1]) # if dual display
  id <- clargs[2] # current run name
  y <-as.numeric(clargs[3])
} else { # sample input parameters, if no command line arguments are given
  dual <-TRUE
  id <-"DualDisp&Run"
  y <- 50
}

print(clargs)
# Currently coordinated with ecoevo_main
vbar <- 3e-3 /y  
dbar <- 1e-5 / y 
final_time <- 2e4 * y

#--- large time plot----
large_dat <-read_csv(paste("outputs/large_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep =""))
large_tstart <- -1e6 * y
print("problem in plot 1")
plot1 <- plot_timeseries(large_dat %>% filter(time %in% c(large_tstart, large_tstart-large_tstart/200, final_time)))+ggtitle("Long Adaptation Period")

#--- small time plot----
# if small has failed (no way of knowing before reading the correct file)
if (dual){
  small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),"_FAILED",sep =""))
  final_time<<- max(small_dat$time)
}else{ small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep =""))}
small_tstart <- -1e3 * y
print("problem in plot 2")
plot2 <- plot_timeseries(small_dat %>% filter(time %in% c(small_tstart, 0, final_time)))+ggtitle("Short Adaptation Period")


plotCombined <- grid.arrange(plot1, plot2, ncol=2)
ggsave(filename =  paste("plots/v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),".pdf",sep ="")
      , plot = plotCombined)