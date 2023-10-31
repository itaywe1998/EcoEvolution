

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
  dual <-TRUE
  id <-"PeriodicCC_Factory"
  vbar <-3e-5
  dbar <-1e-7
  cycles<-3
}

#--- large time plot----
large_dat <-suppressMessages(read_csv(paste("outputs/large_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep ="")))
tstart <- min(large_dat$time)
final_time <-max(large_dat$time)
during_step <- final_time/200
before_step <- -tstart/1000

if(cycles>0){
  req_times <- seq(from=0,to=final_time,l=2*cycles+1)
  obs_times <- seq(from=0,to=2*cycles)
  for (i in seq(from=1,to=2*cycles+1)) {
    obs_times[i] <-large_dat$time[which.min(abs(large_dat$time - req_times[i]))]
  }
  print(obs_times)
  plot1 <- plot_timeseries(large_dat %>% filter(time %in% c(tstart,tstart+before_step, obs_times)))
}else {
  plot1 <- plot_timeseries(large_dat %>% filter(time %in% c(tstart, tstart+before_step, seq(from=0,to=final_time,by=25*during_step))))+ggtitle("Long Adaptation Period")
}



#--- small time plot----
# if small has failed (no way of knowing before reading the correct file)
if (dual){
  suppressMessages(small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),"_FAILED",sep ="")))
}else{ small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep =""))}
tstart <- min(small_dat$time)
final_time<- max(small_dat$time)
before_step <- -tstart/1000
if(cycles>0){
  plot2 <- plot_timeseries(small_dat %>% filter(time %in% c(tstart,tstart+before_step, seq(from=0,to=final_time,l=2*cycles+1))))
}else{
  plot2 <- plot_timeseries(small_dat %>% filter(time %in% c(tstart, tstart+before_step,seq(from=0,to=final_time,by=25*during_step))))+ggtitle("Short Adaptation Period")
}

#---Combined plot----
plotCombined <- grid.arrange(plot1, plot2, nrow=2)
ggsave(filename =  paste("plots/v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),".png",sep ="")
       , plot = plotCombined)