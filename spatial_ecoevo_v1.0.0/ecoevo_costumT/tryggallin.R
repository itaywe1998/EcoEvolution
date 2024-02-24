library(ggplot2)
library(ggallin)

x <- c(-1,-10,-100,-1000,-1e5, 1, 10, 100,1000,1e5)
y <- c(1,2,3, 1,2,3 , 4, 5 , 6 , 7)

df = data.frame(x = x, y = y)

My_Plot = ggplot(
  df, 
  aes(x=x, y=y)) + 
  geom_point() + 
  scale_x_continuous(trans = pseudolog10_trans)

My_Plot