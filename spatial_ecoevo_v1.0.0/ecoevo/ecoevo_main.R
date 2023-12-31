# Copyright (C) 2021 György Barabás
# This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
# you are welcome to redistribute it under certain conditions. for details,
# see the GNU General Public License Agreement (in the file COPYING.txt).

# To run, either execute within R or enter the following at the command prompt:
# Rscript ecoevo.R [vbar] [dbar] [model] [replicate] [outfile]
start <- Sys.time()
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
  S <- as.numeric(clargs[1]) # number of species per trophic level
  vbar <- as.numeric(clargs[2]) # mean genetic variance
  dbar <- as.numeric(clargs[3]) # mean dispersal rate
  model <- clargs[4] # "baseline", "trophic", "Tdep", or "Tdep_trophic"
  replicate <- as.numeric(clargs[5]) # for seeding random number generator
  outfile <- clargs[6] # name of file to save data in (w/ path & extension)
} else { # sample input parameters, if no command line arguments are given
  S <-  2# fifty species per trophic level
  vbar <- 1e-1 # average genetic variance = 0.1 celsius squared
  dbar <- 1e-5 # average dispersal = 1e-5 (100 meters per year)
  model <- "Tdep" # 2 trophic levels & temperature-dependent competition
  replicate <- 1 # replicate number = 1
  outfile <- "myOut" # no output file; make plot instead
}


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

# number of species and number of patches
SR <- S # number of resource species
SC <- 0 # number of consumer species: 0, unless we have...
if (model %in% c("trophic", "Tdep_trophic")) SC <- S # ...consumer species
S <- SR + SC # set S to be the total number of species
L <- 20 # number of patches

# random- and trophic-dependent quantities
set.seed(1000*replicate+321) # set random seed for reproducibility
v <- runif(SR, 0.5*vbar, 1.5*vbar) # resource genetic variances
d <- runif(SR, 0.1*dbar, 10.0*dbar) # resource dispersal rates
rho <- runif(SR, 0.9, 1.1) # resource growth-tolerance tradeoff parameter
a <- matrix(0, S, S) # initialize full competition matrix (resources+consumers)
aP <- matrix(runif(SR*SR, 0.15*0.5, 0.15*1.5), SR, SR) # resource comp coeffs
diag(aP) <- runif(SR, 0.2*0.5, 0.2*1.5) # resource intraspecific comp coeffs
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

# all other parameters
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
Tmin <- -10.0 # initial mean temperature at poles
Cmax <- 30 # projected temperature increase at poles
Cmin <- 27 # projected temperature increase at equator
y <- 100
tstart <- -40*y # starting time (relative to start of climate change at t = 0)
tE <- 200*y # time at which climate change stops (assuming it starts at t = 0)
tend <- tE+250*y # time at which integration ends

# dispersal matrix
mig <- matrix(0, L, L) # initialize dispersal matrix
for (k in 2:L) mig[k-1,k] <- 1 # each species can only migrate to the two
mig <- mig + t(mig) # nearest-neighbor patches

# initial conditions
ninit <- matrix(0, S, L) # reserve memory for initial densities
muinit <- matrix(seq(Tmin, Tmax, l=SR), SR, L) # initial trait means
# initial temperatures
Tempinit <- Temp(seq(from=0, to=1, l=L), 0, tE, Cmax, Cmin, Tmax, Tmin)
for (i in 1:SR) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
# initial traits and densities for consumers
if (model %in% c("trophic", "Tdep_trophic")) {
  muinit <- rbind(muinit, matrix(seq(Tmin, Tmax, l=SC), SC, L))
  for (i in (SR+1):S) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
}
ic <- c(ninit, muinit) # merge initial conditions into a vector

# coerce parameters into a list
pars <- list(SR=SR, SC=SC, S=S, L=L, rho=rho, kappa=kappa, a=a, eta=eta,
             eps=eps, W=W, venv=venv, vmat=vmat, s=s, nmin=nmin, aw=aw, bw=bw,
             Tmax=Tmax, Tmin=Tmin, Th=Th, arate=arate, Cmax=Cmax, Cmin=Cmin,
             tE=tE, d=d, mig=mig, model=model)


# --------------------------- integrate ODEs -----------------------------------
before_step <- 2 * y
before_cc <- ode(y=ic, times=seq(tstart, 0, by=before_step), func=eqs, parms=pars,
                 method="bdf_d") # integrate ODEs before climate change starts
ic <- as.numeric(before_cc[nrow(before_cc),-1]) # final state -> new initial cond.
before_cc <- before_cc %>% # put before-climate-change solution into tidy tibble:
  organize_data(times=seq(from=tstart, to=0, by=before_step), pars=pars) %>%
  filter(time!=0) # remove time point 0 (will be starting point of during_cc)

during_step <- 0.5 * y 
during_cc <- ode(y=ic, times=seq(0, tE, by=during_step), func=eqs, parms=pars,
                 method="bdf_d") # integrate from start to end of climate change
ic <- as.numeric(during_cc[nrow(during_cc),-1]) # final state -> new initial cond.
during_cc <- during_cc %>% # put during-climate-change solution into tidy tibble:
  organize_data(times=seq(from=0, to=tE, by=during_step), pars=pars) %>%
  filter(time!=tE) # remove time point tE (will be starting point of after_cc)

after_step <- 20 * y
after_cc <- ode(y=ic, times=seq(tE, tend, by=after_step), func=eqs, parms=pars,
                method="bdf_d") %>% # integrate from end of climate change to end
  # put after-climate-change solution into tidy tibble:
  organize_data(times=seq(from=tE, to=tend, by=after_step), pars=pars)

# merge data from before, during, and after climate change
dat <- bind_rows(before_cc, during_cc, after_cc) %>%
  # add replicate, genetic var., dispersal rate, and structure as new columns
  mutate(replicate=replicate, vbar=vbar, dbar=dbar, model=model) %>%
  # merge average genetic variance and dispersal into a single column
  mutate(parameterization=paste0("V=", vbar, " d=", dbar)) %>%
  # create regions
  mutate(region=case_when(
    (patch<=round(max(patch)/3))   ~ "polar", # top third of patches are "polar"
    (patch>=round(2*max(patch)/3)) ~ "tropical", # bottom third are "tropical"
    TRUE                           ~ "temperate")) # the rest are "temperate"


# --------------------------- generate output ----------------------------------

if (outfile!="") { # if data file to save to was not specified as empty (""):
  write_csv(dat, path=outfile) # save data to specified file
} # otherwise, create a plot:
  # replace function below with any function from "plotting_functions.R"
  t = min(after_cc$time)
  plot_timeseries(dat %>% filter(time %in% c(tstart,0,t/2,t,max(dat$time))))
  if(min(after_cc$time)-max(during_cc$time)!=during_step) {
    print("Not Continuous")
  }else print("Continuous")
  print(Sys.time()-start)