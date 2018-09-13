rm( list = ls()); gc()
setwd("C:/Users/uos/Dropbox/Extreme value/gev")
library("Deriv")
library("evd")
source("C:/Users/uos/Dropbox/Extreme value/library/fgev.R")

# derivatives for mles
Jaco1=Deriv(expression(f_density_gev(mu,sigma,k,x)),c("mu","sigma","k"))
Hmat1=Deriv(expression(f_density_gev(mu,sigma,k,x)),c("mu","sigma","k"),n=c(hessian=2))
expr_mle <- list (Jaco = Jaco1, Hmat = Hmat1)

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
x2= x1+ z%*%true.beta
x = x2 

alpha = 1
maxiter = 1000
tol = 1e-4
expr = expr_reg


#####################################################
#### Beta gradient method using GEVnewtonRaphson ####
GEVnewtonRaphson_reg1 <- function (x, z, theta0, expr, alpha=1, maxiter=1000, tol = 1e-4)
{
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
      x_d = x- z %*% old_beta 
      mu=old_theta[1]; sigma=old_theta[2]; k=old_theta[3]
      fit_mle <- GEVnewtonRaphson(x = x_d, theta0 = c(mu,sigma,k), 
                                  alpha=1, expr = expr_mle, maxiter = 10, tol = 1e-05)
      new_theta = fit_mle$root
      #cat("theta:", new_theta, '\n')
      del_theta = new_theta - old_theta
      #v = lossfun(x = x - c(z%*%old_beta), new_theta[1], new_theta[2] ,new_theta[3])
      #cat("theta_update:::", v, '\n')
      
    # beta update
      mu=c(new_theta[1]+z%*%old_beta); sigma=new_theta[2]; k=new_theta[3]
      jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
      #jvec_beta =colSums(z*c(grad_mu_mat(x , mu+z%*%old_beta , sigma, k)))
      
      # modify jvec
      #cat("grad_beta::", jvec_beta,"\n")
      new_beta <- old_beta - 0.1*jvec_beta

      #cat("beta:", new_beta, '\n')
      #v = lossfun(x = x - c(z%*%new_beta), new_theta[1], new_theta[2] ,new_theta[3])
      #cat("beta_update:::", v, '\n')
     del_beta <- new_beta - old_beta
    
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
    # v = lossfun(x = x - c(z%*%new_beta), new_theta[1], new_theta[2] ,new_theta[3])
    #cat("V:::", v, '\n')
  }
  return(list(initial = theta0, root = c(new_theta,new_beta), step = niter, del = delta))
}

GEVnewtonRaphson_reg1(x, z, theta0, expr, alpha=1, maxiter=1000)

###################################################
#### Beta newton method using GEVnewtonRaphson ####

GEVnewtonRaphson_reg2 <- function (x, z, theta0, expr, alpha=1, maxiter=1000, tol = 1e-4)
{
  old_theta <- theta0 
  old_beta <- rep(0,ncol(z))
  niter <- 0
  alp <- seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  for (i in 1:maxiter) {
    cat("iter::", i, '\n')
    niter = niter + 1
    # theta update
    x_d = x- z %*% old_beta 
    mu=old_theta[1]; sigma=old_theta[2]; k=old_theta[3]
    fit_mle <- GEVnewtonRaphson(x = x_d, theta0 = c(mu,sigma,k), 
                                alpha=1, expr = expr_mle, maxiter = 10, tol = 1e-05)
    new_theta = fit_mle$root
    cat("theta:", new_theta, '\n')
    del_theta = new_theta - old_theta
    v = lossfun(x = x - c(z%*%old_beta), new_theta[1], new_theta[2] ,new_theta[3])
    cat("theta_update:::", v, '\n')
    
    # beta update
    mu=c(new_theta[1]+z%*%old_beta); sigma=new_theta[2]; k=new_theta[3]
    jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
    #jvec_beta =colSums(z*c(grad_mu_mat(x , mu+z%*%old_beta , sigma, k)))
    hmat_beta = eval(Hmat)[1,1] *t(z)%*%z
    new_beta <- old_beta - 0.1*solve(hmat_beta)%*%jvec_beta
    
    cat("grad_beta::", jvec_beta,"\n")
    cat("beta:", new_beta, '\n')
    v = lossfun(x = x - c(z%*%new_beta), new_theta[1], new_theta[2] ,new_theta[3])
    cat("beta_update:::", v, '\n')
    del_beta <- new_beta - old_beta
    
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

GEVnewtonRaphson_reg2(x, z, theta0, expr, alpha=1, maxiter=1000)


#####################################################
#### Beta newton method without GEVnewtonRaphson ####

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
      x = xdata- z %*% old_beta
      mu=old_theta[1]; sigma=old_theta[2]; k=old_theta[3]
      jvec_theta = apply(eval(Jaco),2,sum)
      hmat_theta = matrix(apply(eval(Hmat),2,sum),3,3)
      new_theta <- old_theta - alpha*solve(hmat_theta)%*%jvec_theta
      del_theta = new_theta - old_theta

      #cat("grad_theta::", jvec_theta,"\n")
      #cat("theta:", new_theta, '\n')
      #v = lossfun(x = x - c(z%*%old_beta), new_theta[1], new_theta[2] ,new_theta[3])
      #cat("theta_update:::", v, '\n')
      
      # beta update
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

GEVnewtonRaphson_reg3(x, z, theta0, expr, alpha=1, maxiter=1000)


#####################################################
#### Beta newton method without GEVnewtonRaphson ####
#### added dealing with Error term               ####

GEVnewtonRaphson_reg4 <- function(xdata, z, theta0, expr, alpha=1, maxiter = 1000, tol = 1e-05) {
  
  if ( !all(1+theta0[3]*(x-theta0[1])/theta0[2]>0) ){
    warning("NaN produced in log by 'theta0'")
  }
  
  old_theta <- theta0 
  old_beta <- rep(0,ncol(z))
  niter <- 0
  alp <- seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  for (i in 1:maxiter) {
    niter = niter + 1
    
    # theta update
    x = xdata- z %*% old_beta
    mu=old_theta[1]; sigma=old_theta[2]; k=old_theta[3]
    jvec_theta = apply(eval(Jaco),2,sum)
    hmat_theta = matrix(apply(eval(Hmat),2,sum),3,3)
    new_theta <- old_theta - alpha* tryCatch(solve(hmat_theta),
                                             error=function(e){solve(hmat_theta+diag(1e-05,nrow(hmat_theta)))
                                             }) %*% jvec_theta
    del_theta = new_theta - old_theta
    
    # beta update
    mu=c(new_theta[1]+z%*%old_beta); sigma=new_theta[2]; k=new_theta[3]
    x = xdata 
    jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
    hmat_beta = eval(Hmat)[1,1] *t(z)%*%z
    new_beta <- old_beta - alpha* tryCatch(solve(hmat_beta),
                                             error=function(e){solve(hmat_beta+diag(1e-05,nrow(hmat_beta)))
                                             }) %*% jvec_beta
    del_beta = new_beta - old_beta
  
    # break rule
    delta=c(del_theta,del_beta)
    if ( norm(as.matrix(delta)) < tol ) { 
      break
    }    
    old_theta = new_theta
    old_beta = new_beta
    
    if (niter == maxiter) {
      warning("Maximum number of iterations 'maxiter' was reached.")
    }
  }
  return(list(initial = theta0, root = c(new_theta,new_beta), step = niter, del = delta))
}

GEVnewtonRaphson_reg4(x, z, theta0, expr, alpha=1, maxiter=1000)



