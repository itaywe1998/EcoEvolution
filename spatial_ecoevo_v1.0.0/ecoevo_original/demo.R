library(gridExtra)
require(ggpmisc) # adding statistics to plots
require(Rcpp) # importing C functions
library(ggplot2)
library(readr)
source("./plotting_functions.R") # various functions for plotting final data
vbar <- 6e-5  
dbar <- 6e-5
id <-6
final_time <- 2e8
large_dat <-read_csv(paste("outputs/large_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep =""))
large_tstart <--1e8
plot1 <- plot_timeseries(large_dat %>% filter(time %in% c(large_tstart, 0, final_time)))+ggtitle("Long Adaptation Period")
small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep =""))
small_tstart <--1e6
plot2 <- plot_timeseries(small_dat %>% filter(time %in% c(small_tstart, 0, final_time)))+ggtitle("Short Adaptation Period")
grid.arrange(plot1, plot2, ncol=2)
