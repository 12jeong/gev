rm( list = ls()); gc()
setwd("C:/Users/uos/Dropbox/Extreme value/gev")

library("Deriv")
library("evd")
library("quantreg")
source("fgevlibrary.R")
source("testlibrary.R")

n = 1000
p=1
true_theta=c(100,40,0.1)
true_beta = c(1)

result <- list()
for ( s in c(1:50) ) {
set.seed(s)
eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
Z = matrix(rnorm(p*n),n,p)
Y = Z%*%true_beta + eps 

##### theta, beta together
result_newton <- tryCatch(GEV_regfull(x=Y, z=Z, theta0=true_theta, beta0=rep(0,dim(Z)[2]), expr=expr_reg, alpha=1, maxiter=100),
                          error=function(e){}, warning=function(e){})$root

##### beta OLS, theta newton
est_beta1 <- lm( Y ~ Z )$coefficient
Y_tilda1 <- Y - Z %*% est_beta1[-1]
# est_beta1 <-  lm( Y-mean(eps) ~ Z )$coefficient
# Y_tilda1 <-   Y - cbind(rep(1,n),Z) %*% est_beta1
est_theta1 <- GEVnewtonRaphson(x=Y_tilda1,theta0=true_theta,expr=expr_mle,maxiter=1,step_theta=1)
result_lm <- c(est_theta1$root,est_beta1[-1])

##### beta QR, theta newton
est_beta2 <- rq( Y ~ Z )$coefficient
Y_tilda2 <- Y - Z %*% est_beta2[-1]
# est_beta2 <- rq( Y-median(eps) ~ Z )$coefficient
# Y_tilda2 <- Y - cbind(rep(1,n),Z) %*% est_beta2
est_theta2 <- GEVnewtonRaphson(x=Y_tilda2,theta0=true_theta,expr=expr_mle,maxiter=1,step_theta=1)
result_rq <- c(est_theta2$root,est_beta2[-1])

#####
result_lm_second <- GEV_regfull(x=Y,z=Z,theta0=est_theta1$root,beta0=est_beta1[-1],alpha=1,maxiter=1)$root
result_rq_second <- GEV_regfull(x=Y,z=Z,theta0=est_theta2$root,beta0=est_beta2[-1],alpha=1,maxiter=1)$root

result[[s]] <- rbind(result_newton,result_lm,result_lm_second,result_rq,result_rq_second)
}

result


# save(result,file="simulation1101.Rdata")

# load("C:/Users/uos/Dropbox/Extreme value/gev/simulation_p1.Rdata")
# load("C:/Users/uos/Dropbox/Extreme value/gev/simulation_p2.Rdata")



