rm( list = ls()); gc()
setwd("C:/Users/uos/Dropbox/Extreme value/gev")

library("Deriv")
library("evd")
source("fgevlibrary.R")


n = 1000
true_theta=c(100,40,0.1)
true_beta = c(1)

result <- list()
for ( s in c(1:100) ) {
set.seed(s)
eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
Z = as.matrix(rnorm(n,0,1))
Y = Z%*%true_beta + eps 

result1 <- tryCatch(GEVnewtonRaphson_reg(x=Y, z=Z,  theta0=true_theta, expr=expr_reg, step_beta = 0.1), error=function(e){},  warning = function(e){}  )

est_beta <- lm( Y ~ Z-1 )$coefficient
Y_tilda <- Y - Z %*% est_beta
est_theta <- GEVnewtonRaphson(x=Y_tilda,theta0=true_theta,expr=expr_mle,maxiter=1)
result2 <- c(est_theta$root,est_beta)

est_theta <- fgev(Y_tilda)
result3 <- c(est_theta$estimate,est_beta)

result[[s]] <- rbind(result1$root,result2,result3)
}


save(n,true_theta,true_beta,result,file="simulation_Error.Rdata")
# load("simulation_Error.Rdata")
# result


getwd()

