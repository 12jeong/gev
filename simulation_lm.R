rm( list = ls()); gc()
setwd("C:/Users/uos/Dropbox/Extreme value/gev")
library("Deriv")
library("evd")
source("C:/Users/uos/Dropbox/Extreme value/library/fgev.R")


set.seed(1)
n = 10000
true_theta=c(100,40,0.1)
# true_beta = 1
Z = as.matrix(rnorm(n,0,1))
eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
Y = Z + eps # true_beta %*% Z + eps


result1 <- GEVnewtonRaphson_reg(x=Y, z=Z, theta0=true_theta, expr=expr_reg, step_theta=0.5, step_beta=0.1, maxiter=1000)

est_beta <- lm( Y ~ Z-1 )$coefficient
Y_tilda <- Y - Z %*% est_beta
est_theta <- fgev(Y_tilda)
result2 <- c(est_theta$estimate,est_beta)

result1$root
result2


GEVnewtonRaphson_reg <- function (x, z, theta0, expr, step_theta=1, step_beta=0.1, maxiter=1000, tol = 1e-6)
{
  old_theta <- theta0 
  old_beta <- rep(0,ncol(z))
  niter <- 0
  alp <- seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  for (i in 1:maxiter) {
    niter = niter + 1
    
    # theta update
    x_d = x- z %*% old_beta 
    mu=old_theta[1]; sigma=old_theta[2]; k=old_theta[3]
    fit_mle <- GEVnewtonRaphson(x = x_d, theta0 = c(mu,sigma,k), step_theta=step_theta)
    new_theta = fit_mle$root
    
    # beta update
    mu=c(new_theta[1]+z%*%old_beta); sigma=new_theta[2]; k=new_theta[3]
    jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
    hmat_beta = eval(Hmat)[1,1] *t(z)%*%z
    new_beta <- old_beta - step_beta *solve(hmat_beta)%*%jvec_beta
    
    # break rule
    if ( abs(max(jvec_beta)) < tol ) {   
      break
    }    
    
    old_theta = new_theta
    old_beta = new_beta
    
    if (niter == maxiter) {
      warning("Maximum number of iterations 'maxiter' was reached.")
    }
  }
  return(list(initial = theta0, root = c(old_theta,old_beta), step = niter))
}

