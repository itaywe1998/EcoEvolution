# Recreating the TPO .Even though here there is no test particle, the full octopule should be close
# to the aproxiamtions results.
# There was an attempt on Naoz 2016, Figure 5 [SMBH perturbuter] but due to
# numeric representations complexities (larger than 1e308 and smaller than 1e-308, intermittently)
# I chose to go with a solar system scaled configuration, figure 18.

suppressPackageStartupMessages({
  suppressWarnings({
    rm(list = ls())
    start <- Sys.time()
    require(deSolve) # solving ordinary differential equations (ODEs)
    library(ggplot2)
    library(patchwork)
    library(pracma)
  })
})
# --- Functions ----
present_profile <- function(name){
  folder <-"/home/itay/EcoEvolution/Kozai_parameters/"
  file <- paste(folder,name,sep = "")
  dat <- load(file)
  p1+p2+p3+plot_layout(ncol=1)
}

kozai_osc <- function(t, state, params){
  # global parameters
  L1 <- unlist(params[1])
  L2 <- unlist(params[2])
  Gtot <- unlist(params[3])
  C3_noG <- unlist(params[4])
  C2_coeff <-unlist(params[5])
  m3 <- unlist(params[6])
  
  # state variables
  e1 <- state[1]
  e2 <- state[2]
  G1 <- state[3]
  G2 <- state[4]
  cosi1 <- state[5]
  cosi2 <- state[6]
  omega1 <- state[7]
  omega2 <-state[8]
  
  if(cosi1>1){
    disp(t/31556926) #for some reason cosi1 is increasing up to >1, which is undefined under acos(), therfore all nans
   disp("Y")
     # I would try different itot seperation, since it worked well with previous examples
  }
  
  # useful parameters per iteration
  C3 <-  C3_noG / G2^5  #V
  C2 <- (e2 ^2 -1) * C3 / C2_coeff #V
  i1 <- acos(cosi1)
  i2 <- acos(cosi2)
  itot <- i1 + i2
  
  cit <- cos(itot)
  co1 <- cos(omega1)
  co2 <- cos(omega2)
  c2o1 <- cos(2*omega1)
  s2o1 <- sin(2*omega1)
  e12 <- e1^2
  e22 <-e2^2
  so1 <- sin(omega1)
  so2 <- sin(omega2)
  sit <- sin(itot)
  sit2 <- sin(itot^2)
  
  B <- 2 +5*e12 - 7 * e12 * c2o1 #V
  A <- 4 + 3*e12 -2.5 * B * sit^2 #V
  #A <- 4 + 3*e12 -2.5 * B * sit2
  cosphi <- -co1 * co2 - cit*so1*so2 #V
  
  # Rates
  
  omega1_dt <-6* C2*(1/G1 * (4*cit^2 + (5*c2o1-1)*(1-e12-cit^2)+cit/G2 * (2+e12* (3-5*c2o1)))) -
    C3*e2*(e1*(1/G2 + cit/G1)* (so1 * so2 * (10 * (3*cit^2 - 1)*(1-e12)+A) - 5*B * cit * cosphi)
          -(1-e12)/(e1 * G1) * (so1*so2 *10 * cit * sit^2 * (1-3*e12) + cosphi*(3*A-10*cit^2 +2))) #V
  #omega1_dt <-6* C2*(1/G1 * (4*cit^2 + (5*c2o1-1)*(1-e12-cit^2)+cit/G2 * (2+e12* (3-5*c2o1)))) -
  #  C3*e2*(e1*(1/G2 + cit/G1)* (so1 * so2 * (10 * (3*cit^2 - 1)*(1-e12)+A) - 5*B * cit * cosphi)
   #        -(1-e12)/(e1 * G1) * (so1*so2 *10 * cit * sit2 * (1-3*e12) + cosphi*(3*A-10*cit^2 +2)))
  
  omega2_dt <- 3*C2*(2*cit/G1 * (2+e12*(3-5*c2o1)) + 1/G2 * (4+6*e12+(5*cit^2-3)*(2+e12*(3-5*c2o1)))) +
    C3*e1*(so1*so2*((4*e22+1)/(e2*G2) * 10 * cit * sit^2 *(1-e12)  -
                      e2*(1/G1 + cit/G2)*(A+10*(3*cit^2-1)*(1-e12))) +
             cosphi*(5*B*cit*e2*(1/G1 +cit/G2) + (4*e22+1)/(e2*G2) * A)  
    ) #V
  
  e1_dt <- C2 * (1-e12)/G1 * (30*e1*sit^2 * s2o1) + C3*e2*(1-e12)/G1 * (35*cosphi*sit^2 *e12 * s2o1
                                                                        -10*cit*sit^2 * co1 * so2 * (1-e12) - A*(so1*co2-cit*co1*so2)) #V
  e2_dt <- -C3*e1*(1-e12)/G2 * (10*cit*sit^2 * (1-e12)*so1*co2+ A*(co1*so2-cit*so1*co2)) #V
  G1_dt <- -C2*30*e12*s2o1*sit^2 +C3*e1*e2*(-35*e12*sit^2 *s2o1 * cosphi +A*(so1*co2-cit*co1*so2) + 10*cit*sit^2 * (1-e12)*co1*so2) #V
  G2_dt <- C3 * e1*e2*(A*(co1*so2-cit*so1*co2) +10*cit*sit^2 * (1-e12)*so1*co2) #V
  ### Not a part of the state vector but yet useful intermediate equations
  H1_dt <- G1/Gtot * G1_dt - G2/Gtot * G2_dt #V
  H2_dt <- -H1_dt #V
  ####
  cosi1_dt <- H1_dt/G1 - G1_dt/G1 * cosi1 #V
  cosi2_dt <- H2_dt/G2 - G2_dt/G2 * cosi2 #V
  
  
  list(c(e1_dt,e2_dt,G1_dt,G2_dt,cosi1_dt,cosi2_dt,omega1_dt,omega2_dt))
}

lyr_to_AU <- function(l){
  return(l*6.3241e4)
}

# Convert mili-arcseconds to distance in the same units as the radius is given in
mas_to_dist <- function(angle, r){
  p<-2*pi*r
  as <- (p/360)/(60*60)
  mas <- 1e-3 * as
  return(angle*mas)
}

to_radians<-function(deg){
  return(deg*pi/180)
}

to_deg<-function(rad){
  return(rad*180/pi)
}



#* A little "Controls & sensitivities" Guide
#*   1.if desired to decrease the gap between the pole and equator - increase obliquity (eps)
#*   even possible to revert the relation above eps=45 degrees, making the pole hotter than the equator
#*   2. the Jupiter-like perturbuter's mass (m3) affects the order of behavior for the eccentricity (e1)
#*    certain values will generate a periodic (ordered) pattern, others will create a rising frequency
#*   3. i1 (inclination for planet) heavily impacts - increasing makes the eccentricity much more extreme and hence the temperature difference over time
#*   4. lowering i2 as well depresses to a well ordered and low difference eccentricity, so it is the combination of both inclinations who rules  
#*   5. Additionally, i1-i2 matters, the larger it is the less moderate the system (e1) is
#*   6. As planets' mass tend to be lower compared to the other masses, it is more vulnerable to extreme changes
#*   the bigger it is the more ordered behavior occurs. In small masses - rising frequency
#*   7. Omegas don't change much (as expected), just the initial point of the same behaviour
#*   8. initial e1 changes the stability of the e1 solution, leading to various resulting 
#*   trends, not very linearly responding.
#*   9. Star mass does not change much of its own, but in order to keep sense,
#*   the mass-luminosity relations obliges to increase much of Ls with minor changes to Ms
#*   10. In a similar manner - a1 also controls the average heat, similarly to Ls
#*   and so both parameters are recommended to stay the same.
#*   11. Albedo is a rather straight forward tuning to average heat 
#*   12. the more distant (a2) the perturbater is the more ordered behavior 
#*   extremely close Jupiter-like breaks the solution (as expected)
#*
#*
#*
#*
#*
#*   
kozai <-function(){
  # -----Global Constants-----
  # working in mks
  AU <- 1.495978707e11 # m
  pc <- 3.08567758e16 # m
  yr <- 31556926 # s
  Ls <- 3.846e26 # W
  Ms <- 1.989e30 # kg
  Mj <- 1.89813e27 # kg
  Me <- 5.972e24 #kg
  GC <- 6.6743e-11
  sigma <- 5.670373e-8 # Stephan-Boltzmann constant [W m^-2 K^-4]
  albedo <- 0.49 # Earth's Albedo
  eps <- to_radians(40) # Earth's Obliquity
  T0 <- -273.15 # Kelvin to Celsius conversion
  
  #---System data ------
  m1 <- 1.4 * Ms # main
  m2 <- 0.3 * Ms # companion
  m3 <- 0.01 * Ms # small perturber
  # 1 is for inner binary
  a1 <- 5 * AU
  e1 <- 	0.5 # THE FULL EQUATION SET IS NOT DEFINED FOR INITIAL ECCENTICITIES ZERO!!!!
  x <-6.75
  i1 <- to_radians(x)
  omega1 <- to_radians(120)
  # 2 is for outer binary
  a2 <- 50 * AU
  e2 <- 0.01 # THE FULL EQUATION SET IS NOT DEFINED FOR INITIAL ECCENTICITIES ZERO!!!!
  i2 <- to_radians(70-x)
  omega2 <- to_radians(0) # NOT GIVEN , will have to play with until stable or reasonable results occur
  
  rel <- a1/a2
  #-------Kozai Inputs-----------
  # constants
  L1 <-m1*m2/(m1+m2) * (GC * (m1+m2) * a1)^0.5 # V
  L2 <-m3*(m1+m2)/(m1+m2+m3) * (GC * (m1+m2+m3) * a2)^0.5 # V
  C3_noG<-  -15/16 * GC^2 /4 *  (m1+m2)^9 / (m1+m2+m3)^4 * (m1-m2)/(m1*m2)^5 * L1^6 / L2^3  * m3^9 #V
  C2_coeff <- 15/4 * (m1-m2)/(m1+m2) * (a1/a2)
  # initial values for state vector
  G1 <-L1 * (1-e1^2)^0.5 # V
  G2 <-L2 * (1-e2^2)^0.5 # V
  H1 <- G1 * cos(i1) # V
  H2 <- G2 * cos(i2)  # V
  # constant
  Gtot <- H1 + H2 # V
  
  ic <- c(e1, e2, G1, G2, cos(i1), cos(i2), omega1, omega2)
  pars <- list(L1 = L1 , L2 = L2, Gtot = Gtot, C3_noG = C3_noG, C2_coeff = C2_coeff, m3 = m3)
  at <- 1e-10
  rt <- 1e-10
  tE <- 1e6 * yr
  step <- tE/1000
  
  workspace_name <- "KozaiSanity_NaozFig4Blue"
  workspace <- paste("~/EcoEvolution/Kozai_parameters/",workspace_name, sep="")
  
  #---- Differential Equation -------
  results <-ode(y=ic, times=seq(0, tE, by=step), func=kozai_osc, parms=pars,
                method="bdf", atol  = at, rtol = rt, maxsteps = 5000)
  diagnostics(results)
  
  disp_results <- results
  disp_results[,1] <- results[,1]/yr
  disp_results[,6:7] <-to_deg(acos(results[,6:7]))
  disp_results[,8:9] <-results[,8:9] %% 360
  times <- disp_results[,1]
  ecc_vec <- disp_results[,2] # e1
  itot_vec <-disp_results[,6]+disp_results[,7] # itot
  disp_results <-as.data.frame(disp_results)
  colnames(disp_results) <- c("Time(years)", "e1", "e2", "G1", "G2", "i1", "i2", "omega1", "omega2")
  
  p1 <-ggplot(data = as.data.frame(cbind(times,itot_vec)))+aes(x= times, y= itot_vec) +geom_line()
  p2<-ggplot(data = as.data.frame(cbind(times,ecc_vec)))+aes(x= times, y= ecc_vec) +geom_line()
  
  p1+p2+plot_layout(ncol=1)
  
  save.image(file = workspace)
  
}
