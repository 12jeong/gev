# tvec 추정이 어떻게 되는지 확인하는 코드

rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")
library(distrEx)
library(RobExtremes) 

S_num = 1
eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/result_scenario',S_num, '.RData',"')"))))
eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/AIC_scenario',S_num, '.RData',"')"))))

lam.grid2[lam.min.vec,]

min.ind =  unlist(lapply(AIC_list,function(x) which.min(x)))
lam.min = as.numeric(names(which.max(table(min.ind))))
lam.min.vec = lam.grid2[lam.min,]
lam.min.vec

i=1; l=132
tvec=result[[i]][[l]]
loc.vec.reg = tvec[3+(1:p)]
sc.vec.reg = tvec[(3+p)+(1:p)]
sh.vec.reg = tvec[(3+2*p)+(1:p)]
loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
par(mfrow=c(2,3))
plot(x1,loc); points(x1,loc.vec, col="blue")
plot(x1,sc); points(x1,sc.vec, col="blue")
plot(x1,sh); points(x1,sh.vec, col="blue")

plot3d(x1,x2,sc)

l=251
tvec=result[[i]][[l]]
loc.vec.reg = tvec[3+(1:p)]
sc.vec.reg = tvec[(3+p)+(1:p)]
sh.vec.reg = tvec[(3+2*p)+(1:p)]
loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
plot(x1,loc); points(x1,loc.vec, col="blue")
plot(x1,sc); points(x1,sc.vec, col="blue")
plot(x1,sh); points(x1,sh.vec, col="blue")

