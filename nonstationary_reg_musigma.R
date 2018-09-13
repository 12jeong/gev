rm( list = ls()); gc()
setwd("C:/Users/uos/Dropbox/Extreme value/gev")
library("Deriv")
library("evd")
source("C:/Users/uos/Dropbox/Extreme value/library/fgev.R")

# derivatives for regression 
logl=expression(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k))
Jaco2=Deriv(logl,c("mu","sigma","k"),combine="cbind")
Hmat2=Deriv(logl,c("mu","sigma","k"),n=c(hessian=2),combine="cbind")
expr_reg <- list (Jaco = Jaco2, Hmat = Hmat2)


#######
n = 1000
p = 3
true.vec = c(120,40,-0.1)
true.beta = c(50,0,0)
theta0 <- true.vec
set.seed(141)
z = matrix(rnorm(p*n),n,p)
x1= rgev(n,loc=true.vec[1] ,scale=true.vec[2],shape=true.vec[3]) 
x2= exp(z%*%true.beta)*x1+ z%*%true.beta

alpha = 1
maxiter = 1000
tol = 1e-4
expr = expr_reg

#######
GEVnewtonRaphson_reg3 <- function(xdata, z, theta0, expr, alpha=1, maxiter = 1000, tol = 1e-05) {
  old_theta <- theta0 
  old_beta <- rep(0,ncol(z))
  niter <- 0
  alp <- seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  for (i in 1:maxiter) {
    #cat("iter::", i, '\n')
    niter = niter + 1
    
    # theta update
    x = ( xdata- z %*% old_beta ) / exp(z%*%old_beta)
    mu=old_theta[1]; sigma=old_theta[2]; k=old_theta[3]
    jvec_theta = apply(eval(Jaco)+ matrix( c(rep(0,n), exp(sigma+z%*%old_beta), rep(0,n)), n, 3) ,2,sum)
    hmat_theta = matrix(apply(eval(Hmat),2,sum),3,3)
    new_theta <- old_theta - alpha*solve(hmat_theta)%*%jvec_theta
    del_theta = new_theta - old_theta
    
    #cat("grad_theta::", jvec_theta,"\n")
    #cat("theta:", new_theta, '\n')
    #v = lossfun(x = x - c(z%*%old_beta), new_theta[1], new_theta[2] ,new_theta[3])
    #cat("theta_update:::", v, '\n')
    
    # beta1 update
    mu=c(new_theta[1]+z%*%old_beta); sigma=new_theta[2]; k=new_theta[3]
    x = xdata 
    jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
    hmat_beta = eval(Hmat)[1,1] *t(z)%*%z
    new_beta <- old_beta - 0.1*solve(hmat_beta)%*%jvec_beta
    del_beta = new_beta - old_beta
    
    #cat("grad_beta::", jvec_beta,"\n")
    #cat("beta:", new_beta, '\n')
    #v = lossfun(x = x - c(z%*%new_beta), new_theta[1], new_theta[2] ,new_theta[3])
    #cat("beta_update:::", v, '\n')
    
    # beta1 update
    mu=c(new_theta[1]+z%*%old_beta); sigma=new_theta[2]; k=new_theta[3]
    x = xdata 
    jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
    hmat_beta = eval(Hmat)[1,1] *t(z)%*%z
    new_beta <- old_beta - 0.1*solve(hmat_beta)%*%jvec_beta
    del_beta = new_beta - old_beta
    
    # break rule
    delta=c(del_theta,del_beta)
    if ( norm(as.matrix(delta)) < tol ) {    ### if (max(abs(del_theta)) < tol )
      break
    }    
    old_theta = new_theta
    old_beta = new_beta
    
    if (niter == maxiter) {
      warning("Maximum number of iterations 'maxiter' was reached.")
    }
    #v = lossfun(x = x - c(z%*%new_beta), new_theta[1], new_theta[2] ,new_theta[3])
    #cat("V:::", v, '\n')
  }
  return(list(initial = theta0, root = c(new_theta,new_beta), step = niter, del = delta))
}

GEVnewtonRaphson_reg3(x2, z, theta0, expr, alpha=1, maxiter=1000)
