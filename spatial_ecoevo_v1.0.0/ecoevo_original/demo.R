library(gridExtra)
require(ggpmisc) # adding statistics to plots
require(Rcpp) # importing C functions
library(ggplot2)
library(readr)
source("./plotting_functions.R") # various functions for plotting final data
dual <-TRUE
vbar <- 3e-5  
dbar <- 1e-7
id <-"DualDisp"
final_time <- 2e6
large_dat <-read_csv(paste("outputs/large_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep =""))
large_tstart <--1e8
plot1 <- plot_timeseries(large_dat %>% filter(time %in% c(large_tstart, large_tstart-large_tstart/200, final_time)))+ggtitle("Long Adaptation Period")
# if small has failed (no way of knowing before reading the correct file)
if (dual){
  small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),"_FAILED",sep =""))
  final_time<<- max(small_dat$time)
}else{
  small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep =""))
}
small_tstart <--1e6
plot2 <- plot_timeseries(small_dat %>% filter(time %in% c(small_tstart, 0, final_time)))+ggtitle("Short Adaptation Period")
grid.arrange(plot1, plot2, ncol=2)