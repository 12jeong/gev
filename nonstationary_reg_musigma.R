rm( list = ls()); gc()
setwd("C:/Users/uos/Dropbox/Extreme value/gev")
library("Deriv")
library("evd")
source("fgevlibrary.R")


#######
n = 1000
p = 3
true.vec = c(120,40,-0.1)
true.beta1 = c(50,0,0)
true.beta2 = c(10,0,0)
theta0 <- true.vec
set.seed(141)
z1 = matrix(rnorm(p*n),n,p)
z2 = matrix(rnorm(p*n),n,p) 
x1= rgev(n,loc=true.vec[1] ,scale=true.vec[2],shape=true.vec[3]) 
x2 = x1 + z1 %*% true.beta1
x3= exp(z2%*%true.beta2)*x1+ z1%*%true.beta1


hist(x1)
hist(x2 - z1 %*% c(0,0,0))
hist(x3-z1%*%c(1,0,0)/exp(z2%*%c(1,0,0)))
hist((x3 - z1 %*% c(50,0,0)) / exp(z2 %*% c(10,0,0)))


#####################################################
#### Beta newton method without GEVnewtonRaphson ####

xdata = x2
theta0=true.vec
expr=expr_reg
step=c(1,0.1,0.1)
maxiter=1000
tol=1e-6

GEVnewtonRaphson_reg3 <- function(xdata, z1, z2, theta0, expr=expr_reg, step=c(1,0.1,0.1), maxiter=1000, tol = 1e-6) {

  old_theta <- theta0 #### 여기서 sigma log? 
  old_beta1 <- rep(0,ncol(z1))
  old_beta2 <- rep(0,ncol(z2))
  niter <- 0
  # alp <- seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  
  for (i in 1:maxiter) {
    #cat("iter::", i, '\n')
    niter = niter + 1
    
    # theta update
    x = (xdata- z1 %*% old_beta1) / exp(z2%*%old_beta2)
    mu=old_theta[1]; sigma=old_theta[2] ; k=old_theta[3]
    jvec_theta = apply(eval(Jaco) * matrix( c(rep(1,n), sigma, rep(1,n)), n, 3) ,2,sum)
    hmat_theta = matrix(apply(eval(Hmat) * cbind(rep(1,n), sigma, rep(1,n), sigma, sigma^2 , sigma , rep(1,n), sigma, rep(1,n)),2,sum),3,3) 
    new_theta <- old_theta - step[1] *solve(hmat_theta)%*%jvec_theta
    
    #cat("grad_theta::", jvec_theta,"\n")
    #cat("theta:", new_theta, '\n')
    #v = lossfun(x = x - c(z%*%old_beta), new_theta[1], new_theta[2] ,new_theta[3])
    #cat("theta_update:::", v, '\n')
    
    # beta1 update
    x = xdata
    mu=c(new_theta[1]+z1%*%old_beta1); sigma=c(exp(log(new_theta[2])+z2%*%old_beta2)); k=new_theta[3]
    jvec_beta1 = apply(eval(Jaco)[,1]*z1,2,sum)
    hmat_beta1 = eval(Hmat)[1,1] *t(z1)%*%z1
    new_beta1 <- old_beta1 - step[2]*solve(hmat_beta1)%*%jvec_beta1
    
    #cat("grad_beta::", jvec_beta,"\n")
    #cat("beta:", new_beta, '\n')
    #v = lossfun(x = x - c(z%*%new_beta), new_theta[1], new_theta[2] ,new_theta[3])
    #cat("beta_update:::", v, '\n')
    
    # beta2 update
    mu=c(new_theta[1]+z1%*%new_beta1); sigma=c(exp(log(new_theta[2])+z2%*%old_beta2)); k=new_theta[3]
    jvec_beta2 = apply(eval(Jaco)[,2]*z2*sigma,2,sum)
    hmat_beta2 = eval(Jaco)[,2]*t(z2)%*%z2*sigma + eval(Hmat)[2,2]*t(z2)%*%z2*sigma^2
    new_beta2 <- old_beta2 - step[3]*solve(hmat_beta2)%*%jvec_beta2  
    

    old_theta1 = new_theta1
    old_beta1 = new_beta1
    old_beta2 = new_beta2
    
    if (niter == maxiter) {
      warning("Maximum number of iterations 'maxiter' was reached.")
    }
    #v = lossfun(x = x - c(z%*%new_beta), new_theta[1], new_theta[2] ,new_theta[3])
    #cat("V:::", v, '\n')
  }
  return(list(initial = theta0, root = c(old_theta,old_beta1,old_beta2), step = niter))
}

GEVnewtonRaphson_reg3(xdata=x2, z1, z2, theta0=true.vec)
