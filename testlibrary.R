# rm( list = ls()); gc()
# setwd("C:/Users/uos/Dropbox/Extreme value/gev")
# 
# library("Deriv")
# library("evd")
# source("fgevlibrary.R")

## check gradient and loss in 1st order approximation
# GEVnewtonRaphson_reg_test2(x=Y, z=Z, theta0=true_theta, expr=expr_reg, step_theta=1, step_beta=0.5, maxiter=1000)

# x=Y; z=Z; theta0=true_theta; expr=expr_reg ; step_tehta=1; step_beta=0.5 ;maxiter=1000; tol=1e-6

GEVnewtonRaphson_test <- function (x, theta0, step_theta=1, expr = expr_mle, maxiter = 5000, tol = 1e-06) {
  
  if ( !all(1+theta0[3]*(x-theta0[1])/theta0[2]>0) ){
    warning("NaN produced in log by 'theta0'")
  }
  
  old_theta <- theta0
  niter <- 0
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  for (i in 1:maxiter) {
    niter = niter + 1
    mu=old_theta[1]; sigma=old_theta[2]; k=old_theta[3]
    jvec = eval(Jaco)
    cat("theta gradient:::",jvec,'\n')
    if ( abs(max(jvec)) < tol ) {    
      break
    }
    hmat = matrix(eval(Hmat),3,3)
    new_theta <- old_theta - step_theta * tryCatch(solve(hmat),
                                                   error= function(e) {solve(hmat+diag(tol,nrow(hmat)))}) %*%jvec
    
    old_theta = new_theta
  }
  return(list(initial = theta0, root = c(old_theta), step = niter, grad=jvec)) 
}


GEVnewtonRaphson_reg_test <- function (x, z, theta0, expr, step_theta=1, step_beta=1, maxiter=10, tol = 1e-6)
{
  old_theta <- theta0 
  old_beta <- rep(0,ncol(z))
  niter <- 0
  alp <- seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  for (i in 1:maxiter) {
    niter = niter + 1
    cat("niter:::", niter,  '\n')
    
    # theta update
    x_d = x- z %*% old_beta 
    mu=old_theta[1]; sigma=old_theta[2]; k=old_theta[3]
    fit_mle <- GEVnewtonRaphson_test(x = x_d, theta0 = c(mu,sigma,k), step_theta=step_theta, maxiter=5000)
    new_theta = fit_mle$root
    cat("new_theta:::",new_theta,'\n')
    # cat("theta gradient:::",fit_mle$grad,'\n')
    # v1 = lossfun(x = x - c(z%*%old_beta), new_theta[1], new_theta[2] ,new_theta[3])
    # cat("theta_update:::",v1,'\n"')
    cat("- - - - - - - - - - - - - - - - - - - - - - - -",'\n')
    # beta update
    for (j in 1:5000){
    mu=c(new_theta[1]+z%*%old_beta); sigma=new_theta[2]; k=new_theta[3]
    jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
    cat("beta gradient:::",jvec_beta,'\n')
      if ( abs(max(jvec_beta)) < tol ) {   
        break
      }  
    hmat_beta = eval(Hmat)[1,1] *t(z)%*%z
    new_beta <- old_beta - step_beta *solve(hmat_beta)%*%jvec_beta
    old_beta <- new_beta
    }
    
    old_theta <- new_theta
    
    # v2 = lossfun(x = x - c(z%*%new_beta), new_theta[1], new_theta[2] ,new_theta[3])
    # cat("beta_update:::",v2,'\n',"==========================",'\n')
    cat("===========================================",'\n')
  }
  return(list(step = niter, initial = theta0, root = c(old_theta,old_beta), grad=c(fit_mle$grad,jvec_beta)))
}



###########################################################################################################
###########################################################################################################

# n = 1000
# p=2
# true_theta=c(100,40,0.1)
# true_beta = c(1,-1)
# set.seed(1)
# eps = rgev(n,loc=true_theta[1], scale=true_theta[2], shape=true_theta[3])
# Z = matrix(rnorm(p*n),n,p)
# Y = Z%*%true_beta + eps 
# 
# GEVnewtonRaphson_reg_test(x=Y, z=Z, theta0=true_theta, expr=expr_reg, step_theta=1, step_beta=0.1, maxiter=2)

