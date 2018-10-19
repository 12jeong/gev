rm(list=ls())
setwd("C:/Users/uos/Dropbox/Extreme value/gev")

library("Deriv")
library("evd")
source("fgevlibrary.R")
source("testlibrary.R")

n = 1000
p=1
true_theta=c(100,40,0.1)
true_beta = c(1)

set.seed(42)
eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
Z = matrix(rnorm(p*n),n,p)
Y = Z%*%true_beta + eps 

est_beta <- lm( Y ~ Z )$coefficient
Y_tilda <- Y - Z %*% est_beta[-1]
est_theta <- GEVnewtonRaphson(x=Y_tilda,theta0=true_theta,expr=expr_mle,maxiter=1,step_theta=1)
result_lm <- c(est_theta$root,est_beta[-1])

GEVnewtonregiter1 <- function(xdata, zdata, theta0, beta0, expr=expr_reg){
  
  # xdata=Y;zdata=Z;theta0=est_theta$root;beta0=est_beta[-1];expr=expr_reg
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  
  x <- xdata
  z <- as.matrix(zdata)
  old_solution <- c(theta0,beta0)
  
  mu=c(old_solution[1]+ z%*%old_solution[-c(1,2,3)]); sigma=old_solution[2]; k=old_solution[3]
  grad = apply(cbind(eval(Jaco),eval(Jaco)[,1]*z),2,sum)
  h1=matrix(apply(eval(Hmat),2,sum),3,3)
  h2=rbind(apply(h1[1,1]*z,2,sum), apply(h1[1,2]*z,2,sum), apply(h1[1,3]*z,2,sum))
  h3=t(h2)
  h4=h1[1,1]*t(z)%*%z
  hmat=rbind(cbind(h1,h2),cbind(h3,h4)) 
  # grad = GEVgradient(j=Jaco,x=x,z=z,mu=mu,sigma=sigma,k=k)
  # hmat = GEVhessian(h=Hmat,x=x,z=z,mu=mu,sigma=sigma,k=k)
  
  new_solution <- old_solution - solve(hmat)%*%grad
  mu=c(new_solution[1]+ z%*%new_solution[-c(1,2,3)]); sigma=new_solution[2]; k=new_solution[3]
  new_grad <- GEVgradient(j=Jaco,x=x,z=z,mu=mu,sigma=sigma,k=k)
  
  return(list(root = c(new_solution), gradient=new_grad))
  
}

# GEVnewtonRaphson(x=Y_tilda,theta0=true_theta,expr=expr_mle,maxiter=1,step_theta=1)
GEVnewtonregiter1(xdata=Y, zdata=Z, theta0=est_theta$root, beta0=est_beta[-1])

# old_solution <- c(96.70437 ,39.44810, 0.1106132, 1.5780940)
# mu=c(old_solution[1]+ Z%*%old_solution[-c(1,2,3)]); sigma=old_solution[2]; k=old_solution[3]
# GEVgradient(Jaco,z=Z)
# apply(cbind(eval(Jaco),eval(Jaco)[,1]*Z),2,sum)
