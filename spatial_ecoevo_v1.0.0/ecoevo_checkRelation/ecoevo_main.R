# Copyright (C) 2021 György Barabás ----
# This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
# you are welcome to redistribute it under certain conditions. for details,
# see the GNU General Public License Agreement (in the file COPYING.txt).

# To run, either execute within R or enter the following at the command prompt:
# Rscript ecoevo.R [vbar] [dbar] [model] [replicate] [outfile]---- 
setwd("~/EcoEvolution/spatial_ecoevo_v1.0.0/ecoevo_checkRelation")
suppressPackageStartupMessages({
  suppressWarnings({
    rm(list = ls())
    start <- Sys.time()
    require(gridExtra)
    require(deSolve) # solving ordinary differential equations (ODEs)
    require(tidyverse) # manipulating and visualizing data
    require(ggpmisc) # adding statistics to plots
    require(Rcpp) # importing C functions
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(dplyr)
    sourceCpp("rhs_eval.cpp") # compile external C functions
    source("./plotting_functions.R") # various functions for plotting final data
  })
})

# ---------------------------- input parameters --------------------------------
arg <- commandArgs(trailingOnly=TRUE)
clargs = unlist(strsplit(arg[1], "#"))
print(clargs)
if (!is.na(clargs)) { # command-line arguments
  small <- as.logical(clargs[1]) # true for short adaptation time, false for long
  id <- clargs[2] # current run name
  C <- as.numeric(clargs[3]) # projected temperature increase at poles
  tE <-as.numeric(clargs[4])
  cycles <- as.numeric(clargs[5])
  updown <- as.logical(clargs[6])
} else { # sample input parameters, if no command line arguments are given
  small <-TRUE
  id <-"Relation"
  cycles <- -1
  updown <- TRUE
}
nmin <- 1e-5 # below this threshold density, genetic variances are reduced
tE <- 1e7
C <- 50 # projected temperature increase at poles
magnitude <-4
kappa <- 1*10^(ceiling(log10(-log(nmin)/tE))+magnitude) # intrinsic mortality parameter


str <- if (small) "small" else "large"
periodic <- if (cycles>0) TRUE else FALSE # Temporary Convention
replicate <- 1 # replicate number = 1
rho <- 1 # resource growth-tolerance tradeoff parameter
b <- 4 
a <- 0.1
#crit_v <- max((rho/kappa)^2 ,1.875*C/tE)
#crit_v <- (rho/kappa)^2
v <- 1000
file <- paste(str,"_time_v",toString(format(v, scientific = TRUE)),"id",toString(id),sep ="")
outfile <- paste("outputs/",file, sep = "") 
workspace <-paste("parameters/",file, sep="")
# --------------------------------functions ------------------------------------

# put the results of the numerical integration into a tidy table
organize_data <- function(dat, times, pars) {
  dat <- dat %>%
    as.data.frame() %>% # convert to data frame (needed for next step)
    as_tibble() %>% # convert to tibble (tidyverse's improved data frame)
    filter(time %in% times) # only keep specified time points
  names(dat)[1] <- "time" # name the first column "time"
  names(dat)[2] <- "n" # naming convention:
  names(dat)[3] <- "m" # (same naming convention)
  dat %>%
    return()
}


# ------------------------------- parameters -----------------------------------

# number of species and number of patches----
T0 <- 25.0 # initial mean temperature at equator
save.image(file = workspace)

# matrices----
ninit <- 1 # reserve memory for initial densities
muinit <- T0
ic <- c(ninit, muinit) # merge initial conditions into a vector

# coerce parameters into a list----
pars <- list(rho=rho, kappa=kappa, b=b,a=a,v=v, nmin=nmin,T0=T0, C=C, tE=tE,periodic=periodic, cycles=cycles, updown=updown)


# --------------------------- integrate ODEs -----------------------------------
#consider changing rtol and atol
at <-1e-2
rt <-1e-2
during_step <- tE/1000
fail_time <- 0
original_tE <- tE
tryCatch({during_cc <-ode(y=ic, times=seq(0, tE, by=during_step), func=eqs, parms=pars,
                          method = "bdf",atol  = at, rtol = rt, maxsteps = 1000)},
         error=function(fail_time){
           message("All Species Extinct")
           fail_time<<-as.numeric(fail_time$message)},
         finally = {
           if (fail_time > 0) {
             outfile <<- paste(outfile,"_FAILED",sep="")
             unlink(workspace) # Deleting old name workspace
             workspace <<- paste(workspace,"_FAILED",sep="")
             save.image(file = workspace)
             tE <<-fail_time*0.9 #alternative for round_any
             during_step <<- tE/200
             # if needed in another place will move to a function
             during_cc <-ode(y=ic, times=seq(0, tE, by=during_step), func=eqs, parms=pars,
                             method = "bdf",atol  = at, rtol = rt, maxsteps = 1000) 
           }
           diagnostics(during_cc)
           during_cc <- during_cc %>% # put during-climate-change solution into tidy tibble:
             organize_data(times=seq(from=0, to=tE, by=during_step), pars = pars) #%>%
           
           # merge data from before, during, and after climate change
           dat <- during_cc
         }) 
# --------------------------- generate output ----------------------------------
print(original_tE-max(during_cc$time))
temp <-(during_cc %>% filter(time %in% c(max(during_cc$time))))
print(mean(temp$n))
#print(min(dat$time[dat$n < 0]))
# req_times <- seq(from=0,to=tE,l=2*cycles+1)
# obs_times <- seq(from=0,to=2*cycles)
# for (i in seq(from=1,to=2*cycles+1)) {
#   obs_times[i] <-during_cc$time[which.min(abs(during_cc$time - req_times[i]))]
# }
# if(mean(temp$n) > 0){ # if ode converged till final time and no significant negative n
#   if (outfile!="") { # if data file to save to was not specified as empty (""):
#     suppressWarnings(write_csv(dat, path=outfile)) }# save data to specified file
#   plot_timeseries(dat %>% filter(time %in% c(tstart,tstart+before_step, obs_times)))
#   #plot_timeseries(dat %>% filter(time %in% c(tstart,tstart+200*before_step,
#   #                                          tstart+400*before_step,tstart+600*before_step,tstart+800*before_step,0)))
# }
print("Final Runtime")
print(Sys.time()-start)

