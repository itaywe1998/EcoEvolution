suppressPackageStartupMessages({
  suppressWarnings({
    rm(list = ls())
    start <- Sys.time()
    require(deSolve) # solving ordinary differential equations (ODEs)
    library(ggplot2)
    library(patchwork)
  })
})
# working in mks
AU <- 1.495978707e11 # m
yr <- 31556926 # s
Ls <- 3.846e26 # W
Ms <- 1.989e30 # kg
Mj <- 1.89813e27 # kg
Me <- 5.972e24 #kg
GC <- 6.6743e-11
sigma <- 5.670373e-8 # Stephan-Boltzmann constant [W m^-2 K^-4]
albedo <- 0.29 # Earth's Albedo
T0 <- -273.15 # Kelvin to Celsius conversion
# --- Functions ----
kozai_osc <- function(t, state, params){
  # global parameters
  L1 <- unlist(params[1])
  L2 <- unlist(params[2])
  Gtot <- unlist(params[3])
  C3_noG <- unlist(params[4])
  C2_coeff <-unlist(params[5])
  
  # state variables
  e1 <- state[1]
  e2 <- state[2]
  G1 <- state[3]
  G2 <- state[4]
  cosi1 <- state[5]
  cosi2 <- state[6]
  omega1 <- state[7]
  omega2 <-state[8]
  
  # useful parameters per iteration
  C3 <- C3_noG / G2^5
  C2 <- (e2 ^2 -1) * C3 / C2_coeff
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
  
  B <- 2 +5*e12 - 7 * e12 * c2o1
  A <- 4 + 3*e12 -2.5 * B * sit^2
  cosphi <- -co1 * co2 - cit*so1*so2
  
  # Rates
  
  omega1_dt <-6* C2*(1/G1 * (4*cit^2 + (5*c2o1-1)*(1-e12-cit^2)+cit/G2 * (2+e12* (3-5*co1)))) -
              C3*e2*(e1*(1/G2 + cit/G1)* (so1 * so2 * (10 * (3*cit^2 - 1)*(1-e12)+A) - 5*B * cit * cosphi)
                       -(1-e12)/(e1 * G1) * (so1*so2 *10 * cit * sit^2 * (1-3*e12) + cosphi*(3*A-10*cit^2 +2)))
  
  omega2_dt <- 3*C2*(2*cit/G1 * (2+e12*(3-5*c2o1)) + 1/G2 * (4+6*e12+(5*cit^2-3)*(2+e12*(3-5*c2o1)))) +
               C3*e1*(so1*so2*((4*e22+1)/(e2*G2) * 10 * cit * sit^2 *(1-e12)  -
                      e2*(1/G1 + cit/G2)*(A+10*(3*cit^2-1)*(1-e12))) +
                      cosphi*(5*B*cit*e2*(1/G1 +cit/G2) + (4*e22+1)/(e2*G2) * A)  
                        )
  
  e1_dt <- C2 * (1-e12)/G1 * (30*e1*sit^2 * s2o1) + C3*e2*(1-e12)/G1 * (35*cosphi*sit^2 *e12 * s2o1
            -10*cit*sit^2 * co1 * so2 * (1-e12) - A*(so1*co2-cit*co1*so2))
  e2_dt <- -C3*e1*(1-e12)/G2 * (10*cit*sit^2 * (1-e12)*so1*co2+ A*(co1*so2-cit*so1*co2))
  G1_dt <- -C2*30*e12*s2o1*sit^2 +C3*e1*e2*(-35*e12*sit^2 *s2o1 * cosphi +A*(so1*co2-cit*co1*so2) + 10*cit*sit^2 * (1-e12)*co1*so2)
  G2_dt <- C3 * e1*e2*(A*(co1*so2-cit*so1*co2) +10*cit*sit^2 * (1-e12)*so1*co2)
  ### Not a part of the state vector but yet useful intermediate equations
  H1_dt <- G1/Gtot * G1_dt - G2/Gtot * G2_dt
  H2_dt <- -H1_dt
  ####
  cosi1_dt <- H1_dt/G1 - G1_dt/G1 * cosi1
  cosi2_dt <- H2_dt/G2 - G2_dt/G2 * cosi2
  
  
  list(c(e1_dt,e2_dt,G1_dt,G2_dt,cosi1_dt,cosi2_dt,omega1_dt,omega2_dt))
}

lyr_to_AU <- function(l){
  return(l*6.3241e4)
}

to_radians<-function(deg){
  return(deg*pi/180)
}
to_deg<-function(rad){
  return(rad*180/pi)
}

# Convert mili-arcseconds to distance in the same units as the radius is given in
mas_to_dist <- function(angle, r){
  p<-2*pi*r
  as <- (p/360)/(60*60)
  mas <- 1e-3 * as
  return(angle*mas)
}

#---HD 202206 data ------
m1 <- 1 * Ms # Star
m2 <- 1 * Me # solid-planet
m3 <- 1 * Mj# gas-planet
a1_mas <-0.1 #mas
d_system <-150 # lyr from sun
# 1 is for inner binary
 #a1 <- mas_to_dist(a1_mas, lyr_to_AU(d_system)) * AU
a1 <- 1 * AU
e1 <- 	0.01
i1 <- to_radians(32)
omega1 <- to_radians(40)
# 2 is for outer planet (HD 202206-c)
a2 <- 3.9 * AU
e2 <- 0.3
i2 <- to_radians(15)
omega2 <- to_radians(0) # NOT GIVEN , will have to play with until stable or reasonable results occur

rel <- a1/a2
#-------Kozai Inputs-----------
# constants
L1 <-m1*m2/(m1+m2) * (GC * (m1+m2) * a1)^0.5
L2 <-m3*(m1+m2)/(m1+m2+m3) * (GC * (m1+m2+m3) * a2)^0.5
C3_noG <-  -15/16 * GC^2 /4 * (m1+m2)^9 / (m1+m2+m3)^4 * (m1-m2)/(m1*m2)^5 * L1^6 / L2^3  * m3^9
C2_coeff <- 15/4 * (m1-m2)/(m1+m2) * (a1/a2)
# initial values for state vector
G1 <-L1 * (1-e1^2)^0.5
G2 <-L2 * (1-e2^2)^0.5
H1 <- G1 * cos(i1)
H2 <- G2 * cos(i2) 
# constant
Gtot <- H1 + H2

ic <- c(e1, e2, G1, G2, cos(i1), cos(i2), omega1, omega2)
pars <- list(L1 = L1 , L2 = L2, Gtot = Gtot, C3_noG = C3_noG , C2_coeff = C2_coeff)
at <- 1e-10
rt <- 1e-10
tE <- 1e8 * yr
step <- tE/500

#---- Differential Equation -------
results <-ode(y=ic, times=seq(0, tE, by=step), func=kozai_osc, parms=pars,
                method="bdf", atol  = at, rtol = rt, maxsteps = 5000)
diagnostics(results)

disp_results <- results
disp_results[,1] <- results[,1]/yr
disp_results[,6:7] <-to_deg(acos(results[,6:7]))
disp_results[,8:9] <-results[,8:9] %% 360
times <- disp_results[,1]
ecc_vec <- disp_results[,2]
disp_results <-as.data.frame(disp_results)
colnames(disp_results) <- c("Time(years)", "e1", "e2", "G1", "G2", "i1", "i2", "omega1", "omega2")

p1<-ggplot(data = as.data.frame(cbind(times,ecc_vec)))+aes(x= times, y= ecc_vec) +
  geom_line()




#----- Eccentricity to Temperature ----------
d <- a1
Lum <- 1 * Ls # originaly 1.084
e <- ecc_vec # This should come from kozai results
# source : https://www.sciencedirect.com/science/article/pii/S1631071310000052
# the relation of average energy received over entire orbit to eccentricity
E <- (Lum*(1-albedo)/(16*sigma*pi*d^2))*(1-e^2)^(-0.5) # energy per area per time
T <- (E)^(1/4) + T0

 p2<-ggplot(data = as.data.frame(cbind(times,T)))+aes(x= times, y= T) +
   geom_line()
 
 p1+p2+plot_layout(ncol=1)