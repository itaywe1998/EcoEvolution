# Copyright (C) 2021 György Barabás ----
# This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
# you are welcome to redistribute it under certain conditions. for details,
# see the GNU General Public License Agreement (in the file COPYING.txt).

# To run, either execute within R or enter the following at the command prompt:
# Rscript ecoevo.R [vbar] [dbar] [model] [replicate] [outfile]---- 
rm(list = ls())

setwd("~/EcoEvolution/spatial_ecoevo_v1.0.0/ecoevo_checkRelation")
suppressPackageStartupMessages({
  suppressWarnings({
    start <- Sys.time()
    require(gridExtra)
    require(deSolve) # solving ordinary differential equations (ODEs)
    # require(tidyverse) # manipulating and visualizing data
    require(ggpmisc) # adding statistics to plots
    require(Rcpp) # importing C functions
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(dplyr)
    source("~/EcoEvolution/kozai.R")
    source("./plotting_functions.R") # various functions for plotting final data
  })
})
sourceCpp("./rhs_eval.cpp") # compile external C functions

# ---------------------------- input parameters --------------------------------
arg <- commandArgs(trailingOnly=TRUE)
clargs = unlist(strsplit(arg[1], "#"))
print(clargs)
if (!is.na(clargs)) { # command-line arguments
  model <- clargs[1] # "baseline", "trophic", "Tdep", or "Tdep_trophic"
  seed <- as.numeric(clargs[2]) # for seeding random number generator
  id <- clargs[3] # current run name
  vbar <- as.numeric(clargs[4]) 
  dbar <- as.numeric(clargs[5]) 
} else { # sample input parameters, if no command line arguments are given
  model <- "baseline" # 2 trophic levels & temperature-dependent competition
  id <-"tryForMaxDiff"
  dbar <- 0 
}
C <- 20.0
tE <-1e6

crit_diff <- 1.875 * C / tE
factor <- 5
vbar <- crit_diff * factor

S <- 1 # fifty species per trophic level
file <- paste("v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep ="")
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
  names(dat)[2] <- paste0("n_", 1, "_", 1) # naming convention:
  names(dat)[3] <- paste0("m_", 1, "_", 1) # naming convention:
  
  dat %>%
    # normalize table by collapsing columns into a key-value column pair
    pivot_longer(cols=2:ncol(.), names_to="variable", values_to="v") %>%
    # split "variable" into value type (density or trait), species, and patch
    separate(variable, c("type", "species", "patch"), sep="_") %>%
    # convert species & patch from string ("1","2",...) to integer (1,2,...)
    mutate(species=as.integer(species), patch=as.integer(patch)) %>%
    # split trait and abundance values into two columns
    pivot_wider(names_from="type", values_from="v") %>%
    # return tidy table
    return()
}


# ------------------------------- parameters -----------------------------------

# scalars----

kappa <- 0.1 # intrinsic mortality parameter
nmin <- 1e-5 # below this threshold density, genetic variances are reduced
aw <- 0.1 # (negative) slope of trait-dependence of tolerance width
bw <- 4 # intercept of trait-dependence of tolerance width
# aw <- 0 # (negative) slope of trait-dependence of tolerance width
# bw <- 0 # intercept of trait-dependence of tolerance width
Tmin <- 15

# matrices----
rho <- 1 # resource growth-tolerance tradeoff parameter

# initial conditions----
ninit <- matrix(1, 1, 1) # reserve memory for initial densities
muinit <- matrix(seq(Tmin, Tmin, l=1), 1, 1) # initial trait means
ic <- c(ninit, muinit) # merge initial conditions into a vector
pars <- list(rho=rho, kappa=kappa, C=C, s=vbar, nmin=nmin, aw=aw, bw=bw,Tmin=Tmin,tE=tE)


# --------------------------- integrate ODEs -----------------------------------
#consider changing rtol and atol
at <-1e-8
rt <-1e-8
maxsteps <- 1000
step <- tE/200
fail_time <- 0
original_tE <- tE

tryCatch({results <-ode(y=ic, times=seq(0, tE-step, by=step), func=eqs, parms=pars,
                        method = "bdf",atol  = at, rtol = rt, maxsteps = maxsteps)
},
error=function(this){
  message(this$message)
  fail_time<<-as.numeric(this$message)
},
finally = {
  if (fail_time > 0) {
    outfile <<- paste(outfile,"_FAILED",sep="")
    unlink(workspace) # Deleting old name workspace
    workspace <<- paste(workspace,"_FAILED",sep="")
    save.image(file = workspace)
    # lets try without the fail_time - step, to see better what happens
    tE <<-floor((fail_time-step)/(step)) * (step) #alternative for round_any
    step <<-tE/200
    # if needed in another place will move to a function
    results <-ode(y=ic, times=seq(0, tE, by=step), func=eqs, parms=pars,
                  method = "bdf",atol  = at, rtol = rt, maxsteps = maxsteps) 
  }
  print(Sys.time()-start)
  start <- Sys.time()
  
  diagnostics(results)
  dat <- results %>% # put during-climate-change solution into tidy tibble:
    organize_data(times=seq(from=0, to=tE, by=step), pars = pars) 
}) 
# --------------------------- generate output ----------------------------------
print(original_tE-max(dat$time))
temp <-(dat %>% filter(time %in% c(max(dat$time))))
print(mean(temp$n))
# # if data file to save to was not specified as empty (""):
# suppressWarnings(write_csv(dat, path=outfile)) # save data to specified file
# #plot_timeseries(dat %>% filter(time %in% seq(from=0,to=tE,by=40*step)))
# plot_landscape(dat %>% filter(patch %in% c(1)))
# 
# toSave <- FALSE
# if (toSave){
#   plt <- plot_timeseries(dat %>% filter(time %in% seq(from=0,to=tE,by=41*step)))
#   #plot_timeseries(dat %>% filter(time %in% seq(from=9e8,to=max(dat$time),by=2*step)))
#   ggsave(filename =  paste("plots/v",toString(format(vbar, scientific = TRUE)),"_d",
#                            toString(dbar),"id",toString(id),".png",sep =""), plot = plt,
#          dpi=300, height = 7, width = 10, units = "in")
#   
#   
# }
# print("R total Runtime")
# print(Sys.time()-start + add)
# #print(Sys.getenv("loop"))

