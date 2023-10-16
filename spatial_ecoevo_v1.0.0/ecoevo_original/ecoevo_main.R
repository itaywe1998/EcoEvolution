# Copyright (C) 2021 György Barabás ----
# This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
# you are welcome to redistribute it under certain conditions. for details,
# see the GNU General Public License Agreement (in the file COPYING.txt).

# To run, either execute within R or enter the following at the command prompt:
# Rscript ecoevo.R [vbar] [dbar] [model] [replicate] [outfile]---- 

rm(list = ls())
start <- Sys.time()
require(gridExtra)
require(deSolve) # solving ordinary differential equations (ODEs)
require(tidyverse) # manipulating and visualizing data
require(ggpmisc) # adding statistics to plots
require(Rcpp) # importing C functions
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
sourceCpp("rhs_eval.cpp") # compile external C functions
source("./plotting_functions.R") # various functions for plotting final data

# ---------------------------- input parameters --------------------------------
clargs <- commandArgs(trailingOnly=TRUE)
if (length(clargs)>0) { # command-line arguments
  model <- clargs[1] # "baseline", "trophic", "Tdep", or "Tdep_trophic"
  small <- as.logical(clargs[2]) # true for short adaptation time, false for long
  seed <- as.numeric(clargs[3]) # for seeding random number generator
  id <- clargs[4] # current run name
  y <- as.numeric(clargs[5]) # scaling factor
  x <- as.numeric(clargs[6]) # scaling factor
} else { # sample input parameters, if no command line arguments are given
  model <- "Tdep" # 2 trophic levels & temperature-dependent competition
  small <-FALSE
  id <-"BiggerPeriodicCC"
  seed <- 3695
  y <- 100
  x <- 1
  }
S <- 4 # fifty species per trophic level
vbar <- 3e-3 /y  # average genetic variance in Celsius squared
dbar <- (1e-5 / y) / x  # average dispersal (1e-7 <=> 1 meter per year)
# more precisely, in units of pole to equator distance , which is ~100,000 km (1e7 meter)

if (small){
  ts<--1e3 * y
  str<-"small"
}else {
  ts<--1e6 * y
  str<-"large"
}
replicate <- 1 # replicate number = 1
file <- paste(str,"_time_v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),sep ="")
outfile <- paste("outputs/",file, sep = "") 
workspace <-paste("parameters/",file, sep="")
# --------------------------------functions ------------------------------------

# return matrix W[i,j], which is nonzero if consumer i eats resource j;
# SR is the number of resource, SC the number of consumer species
generate_network <- function(SR, SC) {
  w <- matrix(0, SR+SC, SR+SC) # initialize adjacency matrix
  for (i in 1:SR) { # determine which resources each consumer eats: it must eat
    indices <- sort(c(i, sample((1:SC)[-i], 1))) # the one with
    w[i+SR,indices] <- 1 # matching trait value, plus a fixed number of
  } # randomly assigned ones (in this case 4 more, for 5 resources per consumer)
  omega <- numeric(0) # initialize matrix of consumption efforts
  rsum <- rowSums(w) # omega[i,j] is the proportion of i's consumption rate,
  for (i in 1:(SR+SC)) omega <- cbind(omega, rsum) # targeted at consuming j
  omega[omega!=0] <- 1/omega[omega!=0] # if not 0, set to proportion
  W <- unname(w*omega) # only the product of w and omega is used
  return(W)
}

# put the results of the numerical integration into a tidy table
organize_data <- function(dat, times, pars) {
  dat <- dat %>%
    as.data.frame() %>% # convert to data frame (needed for next step)
    as_tibble() %>% # convert to tibble (tidyverse's improved data frame)
    filter(time %in% times) # only keep specified time points
  names(dat)[1] <- "time" # name the first column "time"
  index <- 2 # keep track of which column we are naming with this counter
  for (k in 1:pars$L) {
    for (i in 1:pars$S) { # name columns for densities
      names(dat)[index] <- paste0("n_", i, "_", k) # naming convention:
      index <- index + 1 # "type_species_patch" - type is either m (trait),
    } # or n (density)
  }
  for (k in 1:pars$L) {
    for (i in 1:pars$S) { # name columns for trait values
      names(dat)[index] <- paste0("m_", i, "_", k) # (same naming convention)
      index <- index + 1
    }
  }
  dat %>%
    # normalize table by collapsing columns into a key-value column pair
    pivot_longer(cols=2:ncol(.), names_to="variable", values_to="v") %>%
    # split "variable" into value type (density or trait), species, and patch
    separate(variable, c("type", "species", "patch"), sep="_") %>%
    # convert species & patch from string ("1","2",...) to integer (1,2,...)
    mutate(species=as.integer(species), patch=as.integer(patch)) %>%
    # split trait and abundance values into two columns
    pivot_wider(names_from="type", values_from="v") %>%
    # trophic level (tl): species with index greater than SR are consumers ("C"),
    # the rest are resources ("R")
    mutate(tl=ifelse(species>SR, "C", "R")) %>%
    # return tidy table
    return()
}


# ------------------------------- parameters -----------------------------------

# number of species and number of patches----
SR <- S # number of resource species
SC <- 0 # number of consumer species: 0, unless we have...
if (model %in% c("trophic", "Tdep_trophic")) SC <- S # ...consumer species
S <- SR + SC # set S to be the total number of species
L <- 20 # number of patches

# scalars----
set.seed(seed) # set random seed for reproducibility
v <- runif(SR, 1.0*vbar, 2.0*vbar) # resource genetic variances
d <- runif(SR, 1.0*dbar, 2.0*dbar) # resource dispersal rates

periodic <- TRUE
cycles <- 5

kappa <- 0.1 # intrinsic mortality parameter
venv <- vbar # environmental variance
vmat <- matrix(rep(v, L), S, L) # genetic variances at each patch
s <- v + venv # species' total phenotypic variances
eta <- 1 # competition width (centigrade; only for Tdep and Tdep_trophic)
eps <- c(rep(0, SR), rep(0.3, SC)) # feeding efficiency of consumers
nmin <- 1e-5 # below this threshold density, genetic variances are reduced
aw <- 0.1 # (negative) slope of trait-dependence of tolerance width
bw <- 4 # intercept of trait-dependence of tolerance width
Tmax <- 25.0 # initial mean temperature at equator
Tmin <- Tmax-40 # initial mean temperature at poles
Cmax <- 15 # projected temperature increase at poles
Cmin <- 7 # projected temperature increase at equator
#Cmax <- 0
#Cmin <- 0 
tstart <- ts # starting time (relative to start of climate change at t = 0)
tE <- 2e4*y*10 # time at which climate change stops (assuming it starts at t = 0)
save.image(file = workspace)

# matrices----
rho <- runif(SR, 0.1, 11) # resource growth-tolerance tradeoff parameter
a <- matrix(0, S, S) # initialize full competition matrix (resources+consumers)
# assigned 0.7 & 0.9 instead of 0.5 & 1.5 as margins in aP, to lower competition
aP <- matrix(runif(SR*SR, 0.15*0.4, 0.15*0.9), SR, SR) # resource comp coeffs
diag(aP) <- runif(SR, 0.2*0.4, 0.2*0.9) # resource intraspecific comp coeffs
a[1:SR,1:SR] <- aP # top left block: resources
W <- matrix(0, S, S) # create feeding network: nothing if no consumers
Th <- rep(1, S) # handling times in type II f.r. (dummy value if no consumers)
arate <- rep(1, S) # attack rates in type II f.r. (dummy value if no consumers)
if (model %in% c("trophic", "Tdep_trophic")) {
  v <- c(v, runif(SC, 0.5*vbar, 1.5*vbar)) # add consumer genetic variances
  d <- c(d, runif(SC, 0.1*dbar, 10.0*dbar)) # add consumer dispersal rates
  rho <- c(rho, runif(SC, 0.9*0.1, 1.1*0.1)) # add consumer tradeoff parameters
  aH <- matrix(0, SC, SC) # initialize competition matrix (consumers)
  a[(SR+1):S,(SR+1):S] <- aH # bottom right: consumers
  W <- generate_network(SR, SC) # trophic feeding network
  Th[(SR+1):S] <- runif(S-SR, 0.5, 1) # handling times in type II f.r.
  arate[(SR+1):S] <- runif(S-SR, 1, 10) # attack rates in type II f.r.
}



# dispersal matrix----
mig <- matrix(0, L, L) # initialize dispersal matrix
for (k in 2:L) mig[k-1,k] <- 1 # each species can only migrate to the two
mig <- mig + t(mig) # nearest-neighbor patches

# initial conditions----
ninit <- matrix(0, S, L) # reserve memory for initial densities
muinit <- matrix(seq(Tmin, Tmin, l=SR), SR, L) # initial trait means
# Edit ! all initial species start with same location controlled de-facto by muninit 
# initial temperatures
Tempinit <- Temp(seq(from=0, to=1, l=L), 0, tE, Cmax, Cmin, Tmax, Tmin, periodic, cycles)
for (i in 1:SR) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
# initial traits and densities for consumers
if (model %in% c("trophic", "Tdep_trophic")) {
  muinit <- rbind(muinit, matrix(seq(Tmin, Tmax, l=SC), SC, L))
  for (i in (SR+1):S) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
}
ic <- c(ninit, muinit) # merge initial conditions into a vector

# coerce parameters into a list----
pars <- list(SR=SR, SC=SC, S=S, L=L, rho=rho, kappa=kappa, a=a, eta=eta,
             eps=eps, W=W, venv=venv, vmat=vmat, s=s, nmin=nmin, aw=aw, bw=bw,
             Tmax=Tmax, Tmin=Tmin, Th=Th, arate=arate, Cmax=Cmax, Cmin=Cmin,
             tE=tE, d=d, mig=mig, model=model, periodic=periodic, cycles=cycles)


# --------------------------- integrate ODEs -----------------------------------
#consider changing rtol and atol
at <-1e-14
rt <-1e-14
before_step <- -tstart/1000
tryCatch({before_cc <-ode(y=ic, times=seq(tstart, 0, by=before_step), func=eqs, parms=pars,
       method="bdf", atol  = at, rtol = rt, maxsteps = 10000)},
      error=function(e){message("All Species Extinct")
                        return(NA)}) # integrate ODEs before climate change starts
diagnostics(before_cc)
ic <- as.numeric(before_cc[nrow(before_cc),-1]) # final state -> new initial cond.
before_cc <- before_cc %>% # put before-climate-change solution into tidy tibble:
  organize_data(times=seq(from=tstart, to=0, by=before_step), pars = pars) %>%
  filter(time!=0) # remove time point 0 (will be starting point of during_cc)

print(Sys.time()-start)

during_step <- tE/200
at <-1e-14
rt <-1e-14
fail_time <- 0
tryCatch({during_cc <-ode(y=ic, times=seq(0, tE, by=during_step), func=eqs, parms=pars,
      method = "bdf",atol  = at, rtol = rt, maxsteps = 10000)},
      error=function(fail_time){
          message("All Species Extinct")
          fail_time<<-as.numeric(fail_time$message)},
      finally = {
        if (fail_time > 0) {
          outfile <<- paste(outfile,"_FAILED",sep="")
          unlink(workspace) # Deleting old name workspace
          workspace <<- paste(workspace,"_FAILED",sep="")
          save.image(file = workspace)
          during_cc <-ode(y=ic, times=seq(0, fail_time-during_step, by=during_step), func=eqs, parms=pars,
                          method = "bdf",atol  = at, rtol = rt, maxsteps = 10000) 
        }
        diagnostics(during_cc)
        during_cc <- during_cc %>% # put during-climate-change solution into tidy tibble:
          organize_data(times=seq(from=0, to=tE, by=during_step), pars = pars) #%>%
        
        # merge data from before, during, and after climate change
        dat <- bind_rows(before_cc, during_cc) %>%
          # add replicate, genetic var., dispersal rate, and structure as new columns
          mutate(replicate=replicate, vbar=vbar, dbar=dbar, model=model) %>%
          # merge average genetic variance and dispersal into a single column
          mutate(parameterization=paste0("V=", vbar, " d=", dbar)) %>%
          # create regions
          mutate(region=case_when(
            (patch<=round(max(patch)/3))   ~ "polar", # top third of patches are "polar"
            (patch>=round(2*max(patch)/3)) ~ "tropical", # bottom third are "tropical"
            TRUE                           ~ "temperate")) # the rest are "temperate"
        
      })  # integrate from start to end of climate change
# --------------------------- generate output ----------------------------------
print(tE-max(during_cc$time))
temp <-(during_cc %>% filter(time %in% c(max(during_cc$time))))
print(mean(temp$n))
#print(min(dat$time[dat$n < 0]))
if(mean(temp$n) > 0){ # if ode converged till final time and no significant negative n
  if (outfile!="") { # if data file to save to was not specified as empty (""):
    write_csv(dat, path=outfile) }# save data to specified file
  plot_timeseries(dat %>% filter(time %in% c(tstart,tstart+before_step, 0,25*during_step, 50*during_step,75*during_step,
                                             100*during_step,125*during_step, 150*during_step,175*during_step,tE)))
  }
print(Sys.time()-start)

