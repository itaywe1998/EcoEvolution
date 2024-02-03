setwd("~/EcoEvolution/spatial_ecoevo_v1.0.0/ecoevo_original")

suppressPackageStartupMessages({
  suppressWarnings({
    library(gridExtra)
    require(ggpmisc) # adding statistics to plots
    require(Rcpp) # importing C functions
    library(ggplot2)
  })
})

smoothstep<-function(x){
  return(x*x*x*(10.0+x*(-15.0+6.0*x)));
}
periodic <- FALSE
cycles <- 0
t <- seq(0, tE, by=tE/200)
if (periodic){
  T_mid <- (Tmax-Tmin)*0.5+Tmin+((Cmin-Cmax)*0.5+Cmax)*sin((t/tE)*cycles*pi)
  T_polar <- Tmin+Cmax*sin((t/tE)*cycles*pi)
  T_equator <- Tmax+Cmin*sin((t/tE)*cycles*pi)
}else{
  T_mid <- (Tmax-Tmin)*0.5+Tmin+((Cmin-Cmax)*0.5+Cmax)*smoothstep((t/tE))
  T_polar <- Tmin+Cmax*smoothstep((t/tE))
  T_equator <- Tmax+Cmin*smoothstep((t/tE))
}

df <- data.frame(t =t, T_polar=T_polar, T_equator=T_equator,T_mid=T_mid)
colors <- c("Polar" = "blue", "Equator" = "red", "Mid" = "green")
s <- 0.5
p <- ggplot(df, aes(x = t)) +
   geom_line(aes(y = T_polar, color = "Polar"), size = s) +
   geom_line(aes(y =  T_equator, color = "Equator"), size = s) +
   geom_line(aes(y = T_mid, color = "Mid"), size = s) +
   labs(x = "Time (years)",
        y = "Temprature (\u00B0C)",
        color = "") +
   scale_color_manual(values = colors)+
  scale_y_continuous(breaks = seq(-40, 30, by = 10))
print(p)

#v3e-05_d1e-07idBiggerCompetition_Seed7234

dest_fold <- "/home/itay/EcoEvolution/spatial_ecoevo_v1.0.0/Ecoevo Paper Examples/Plots/"
runname <- "v3e-05_d1e-07idBiggerCompetition_Seed7234_Paper"
ggsave(filename =  paste(dest_fold, "WeatherProfile_",runname,".png",sep ="")
       , plot = p, width = 6, height = 4, units = "in")
