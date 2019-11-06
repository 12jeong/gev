rm(list=ls()); gc()
setwd("~/GitHub/gev")
source("./lib/pack.R")
source("./lib/sgevlibrary.R")
load("./kma_data/Pr_46.RData")
load("./paper_co/result_spatial_real.RData")

lambdaset

ns <- length(unique(Pr_46$stnlds))
fit = result[[21]]
x = unique(Pr_46$lat)
y = unique(Pr_46$long)
z = drop(fit[1]+Z%*%tail(fit,p))

library(fields)
Tpsout = fastTps(x=data.frame(x,y),Y=z,theta=1)
x2 = seq(min(x),max(x),length.out=50)
y2 = seq(min(y),max(y),length.out=50)
data = expand.grid(x2,y2)
fitted = matrix(predict(Tpsout,data),50,50)


plot_ly(x=x2,y=y2,z=fitted, type="surface") %>%
  layout(scene = list(
    xaxis = list(title = "latitude"),
    yaxis = list(title = "longitude"),
    zaxis = list(title = "")))
