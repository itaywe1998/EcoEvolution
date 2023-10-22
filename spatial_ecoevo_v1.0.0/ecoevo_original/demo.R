

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
  y <-as.numeric(clargs[3])
  x <-as.numeric(clargs[4])
  cycles<-as.numeric(clargs[5])
  periodic<-as.logical(clargs[6])
} else { # sample input parameters, if no command line arguments are given
  dual <-TRUE
  id <-"PeriodicCC_Factory"
  y <- 100
  x <- 1
  cycles<-5
  periodic<-TRUE
}

#print(clargs)


##### ALWAYS COORDINATE THIS WITH ECOEVO_MAIN.R !!!!!!   ########
vbar <- 3e-3 /y  
dbar <- (1e-5 / y) / x *(cycles)
#final_time <- 2e4 * y * (2*cycles)
final_time <- 2e4 * y 
half_cycle_duration <- final_time/(2*cycles)
during_step <- final_time/200
##### ALWAYS COORDINATE THIS WITH ECOEVO_MAIN.R !!!!!!   ########



#--- large time plot----
large_dat <-suppressMessages(read_csv(paste("outputs/large_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep ="")))
large_tstart <- -1e6 * y
before_step <- -large_tstart/1000
plot1 <- plot_timeseries(large_dat %>% filter(time %in% c(large_tstart, large_tstart+before_step, 0 ,25*during_step, 50*during_step,
                                                          75*during_step,100*during_step,125*during_step, 150*during_step, final_time)))+ggtitle("Long Adaptation Period")

#--- small time plot----
# if small has failed (no way of knowing before reading the correct file)
if (dual){
  suppressMessages(small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),"_FAILED",sep ="")))
  #suppressMessages(small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep ="")))
  final_time<<- max(small_dat$time)
}else{ small_dat <-read_csv(paste("outputs/small_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep =""))}
small_tstart <- -1e3 * y
before_step <- -small_tstart/1000
plot2 <- plot_timeseries(small_dat %>% filter(time %in% c(small_tstart, small_tstart+before_step,0, 25*during_step, 50*during_step,75*during_step,
                                                          100*during_step,125*during_step, 150*during_step,175*during_step,final_time)))+ggtitle("Short Adaptation Period")


plotCombined <- grid.arrange(plot1, plot2, nrow=2)
ggsave(filename =  paste("plots/v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),".png",sep ="")
      , plot = plotCombined)