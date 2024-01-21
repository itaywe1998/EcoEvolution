# Copyright (C) 2021 György Barabás ----
# This program comes with ABSOLUTELY NO WARRANTY. This is free software, and
# you are welcome to redistribute it under certain conditions. for details,
# see the GNU General Public License Agreement (in the file COPYING.txt).

# To run, either execute within R or enter the following at the command prompt:
# Rscript ecoevo.R [vbar] [dbar] [model] [replicate] [outfile]---- 
setwd("~/EcoEvolution/spatial_ecoevo_v1.0.0/ecoevo_costumT")
suppressPackageStartupMessages({
  suppressWarnings({
    rm(list = ls())
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
    sourceCpp("./rhs_eval.cpp") # compile external C functions
  })
})

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
  model <- "Tdep" # 2 trophic levels & temperature-dependent competition
  id <-"LargeAdaptTimeMildAlive"
  seed <- 3690
  vbar <- 2e-3 # average genetic variance in Celsius squared 
  dbar <- 1e-5 # average dispersal (1e-7 <=> 1 meter per year)
  # more precisely, in units of pole to equator distance , which is ~100,000 km (1e7 meter)
}
S <- 4 # fifty species per trophic level
replicate <- 1 # replicate number = 1
set.seed(NULL) #unsets random seed, it is set again before the important stuff 
run_indicator <- sample(1:20000,1)
file <- paste("v",toString(format(vbar, scientific = TRUE)),"_d",toString(dbar),"id",toString(id),toString(run_indicator),sep ="")
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

kappa <- 0.1 # intrinsic mortality parameter
venv <- vbar # environmental variance
vmat <- matrix(rep(v, L), S, L) # genetic variances at each patch
s <- v + venv # species' total phenotypic variances
eta <- 1 # competition width (centigrade; only for Tdep and Tdep_trophic)
eps <- c(rep(0, SR), rep(0.3, SC)) # feeding efficiency of consumers
nmin <- 1e-5 # below this threshold density, genetic variances are reduced
aw <- 0.1 # (negative) slope of trait-dependence of tolerance width
bw <- 4 # intercept of trait-dependence of tolerance width

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

# Temperatures----
old_profile <- TRUE
if (old_profile){
  wksp_name <- "LargeAdaptTimeMild"
  kozai_wksp <- paste("~/EcoEvolution/Kozai_parameters/",wksp_name, sep="")
  tmp.env <- new.env() # create a temporary environment
  load(kozai_wksp, envir=tmp.env) # load workspace into temporary environment
  T_kozai <- tmp.env$Tvec 
  rm(tmp.env) 
}else{
  T_kozai <- kozai()
}

lT <- length(T_kozai[,1])
save.image(file = workspace)
Tmin <- unname(T_kozai[1,2])
Tmax <- unname(T_kozai[1,3])
# initial conditions----
ninit <- matrix(0, S, L) # reserve memory for initial densities
muinit <- matrix(seq(Tmin, Tmin, l=SR), SR, L) # initial trait means
# Edit ! all initial species start with same location controlled de-facto by muninit 
# initial temperatures
Tempinit <- Temp(seq(from=0, to=1, l=L), Tmax, Tmin)
for (i in 1:SR) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
# initial traits and densities for consumers
if (model %in% c("trophic", "Tdep_trophic")) {
  muinit <- rbind(muinit, matrix(seq(Tmin, Tmin, l=SC), SC, L))
  for (i in (SR+1):S) ninit[i,] <- exp(-(muinit[i,1]-Tempinit)^2/(2*2^2))
}
ic <- c(ninit, muinit) # merge initial conditions into a vector

# coerce parameters into a list----
pars <- list(SR=SR, SC=SC, S=S, L=L, rho=rho, kappa=kappa, a=a, eta=eta,
             eps=eps, W=W, venv=venv, vmat=vmat, s=s, nmin=nmin, aw=aw, bw=bw,
             Th=Th, arate=arate,d=d, mig=mig, model=model,T_kozai=T_kozai, lT=lT)


# --------------------------- integrate ODEs -----------------------------------
#consider changing rtol and atol
at <-1e-8
rt <-1e-8
maxsteps <- 7000
tE <-tail(T_kozai, n=1)[1]
step <- unname(T_kozai[2,1])-unname(T_kozai[1,1])
#step <-1e5
fail_time <- 0
original_tE <- tE
tryCatch({results <-ode(y=ic, times=seq(0, tE, by=step), func=eqs, parms=pars,
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
             # if needed in another place will move to a function
             results <-ode(y=ic, times=seq(0, tE, by=step), func=eqs, parms=pars,
                             method = "bdf",atol  = at, rtol = rt, maxsteps = maxsteps) 
           }
           diagnostics(results)
           results <- results %>% # put during-climate-change solution into tidy tibble:
             organize_data(times=seq(from=0, to=tE, by=step), pars = pars) #%>%
           
           temperature<-rep(seq(Tmin, Tmax, l=L),each = S)
           limit <- nrow(results)/(S*L)
           for (i in 2:limit){
             Tlow <- unname(T_kozai[i,2])
             Thigh <- unname(T_kozai[i,3])
             temperature<-c(temperature, rep(seq(Tlow, Thigh, l=L),each = S))
           }
           dat <-results%>%mutate(Tenv=temperature)
         }) 
# --------------------------- generate output ----------------------------------
print(original_tE-max(dat$time))
temp <-(dat %>% filter(time %in% c(max(dat$time))))
print(mean(temp$n))
# if data file to save to was not specified as empty (""):
suppressWarnings(write_csv(dat, path=outfile)) # save data to specified file
# plot_timeseries(dat %>% filter(time %in% c(0,step)))
plot_timeseries(dat %>% filter(time %in% seq(from=tE-step,to=tE,by=step)))
toSave <- FALSE
if (toSave){
  plt <- plot_timeseries(dat %>% filter(time %in% seq(from=0,to=tE,by=100*step)))
  #plot_timeseries(dat %>% filter(time %in% seq(from=9e8,to=max(dat$time),by=2*step)))
  ggsave(filename =  paste("plots/v",toString(format(vbar, scientific = TRUE)),"_d",
                           toString(dbar),"id",toString(id),".png",sep =""), plot = plt,
         dpi=300, height = 7, width = 10, units = "in")
  
  
}
print("Final Runtime")
print(Sys.time()-start)

