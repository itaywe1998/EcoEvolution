heavy <- FALSE
if(heavy){
  #  LargeAdaptTimeMildExtendedSpanMemHeavy results
  x <-c(50,200,500,1000,2000,4000) # tE in kyrs
  y <- c(1.08,1.31,1.92,3.05,6.62,13.27) #R part runtime
  z <- c(1.69,28.7,1.07 * 60,3.82*60, 7.71 * 60, 21.5 * 60) #cpp part runtime
}else{
  #  LargeAdaptTimeMildExtendedSpan results
  x <-c(50,200,500,1000,2000,4000, 5000) # tE in kyrs
  y <- c(0.72,0.86,1.49,2.6,5.09,10.52, 14.89) #R part runtime
  z <- c(1.56,22.36,48.27,2.44 * 60, 4.01 * 60, 9.52*60, 11.85 * 60) #cpp part runtime
}



df <- data.frame(x, y) 
df2 <- data.frame(x,z) 

##### 
#fit1 <- lm(y~x, data=df)
# fit2 <- lm(y~poly(x,2,raw=TRUE), data=df)
# fit3 <- lm(y~poly(x,3,raw=TRUE), data=df)
# fit4 <- lm(y~poly(x,4,raw=TRUE), data=df)
# fit5 <- lm(y~poly(x,5,raw=TRUE), data=df)
#plot(df$x, df$y, pch=19, xlab='tE (kyrs)', ylab='R runtime (sec)')

#####
fit1 <- lm(z~x, data=df2)
fit2 <- lm(z~poly(x,2,raw=TRUE), data=df2)
fit3 <- lm(z~poly(x,3,raw=TRUE), data=df2)
# fit4 <- lm(z~poly(x,4,raw=TRUE), data=df2)
# fit5 <- lm(z~poly(x,5,raw=TRUE), data=df2)

plot(df2$x, df2$z, pch=19, xlab='tE (kyrs)', ylab='Cpp runtime (sec)')

#define x-axis values
x_axis <- x

#add curve of each model to plot
lines(x_axis, predict(fit1, data.frame(x=x_axis)), col='green')
lines(x_axis, predict(fit2, data.frame(x=x_axis)), col='red')
# lines(x_axis, predict(fit3, data.frame(x=x_axis)), col='purple')
# lines(x_axis, predict(fit4, data.frame(x=x_axis)), col='blue')
# lines(x_axis, predict(fit5, data.frame(x=x_axis)), col='orange')