# rm( list = ls()); gc()
# setwd("C:/Users/uos/Dropbox/Extreme value/gev")
# 
# library("Deriv")
# library("evd")
# source("fgevlibrary.R")

# x=Y; z=Z; theta0=true_theta; expr=expr_reg ; step_tehta=1; step_beta=0.5 ;maxiter=1000; tol=1e-6

GEVnewtonRaphson_test <- function (x, theta0, step_theta=1, expr = expr_mle, maxiter = 5000, tol = 1e-06) {
  if ( !all(1+theta0[3]*(x-theta0[1])/theta0[2]>0) ){
    warning("NaN produced in log by 'theta0'")
  }
  old_theta <- theta0
  niter <- 0
  alp <- 1-seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  for (i in 1:maxiter) {
    niter = niter + 1
    mu=old_theta[1]; s=old_theta[2]; k=old_theta[3]
    jvec = eval(Jaco)
    cat("theta gradient:::",jvec,'\n')
    if ( abs(max(jvec)) < tol ) {    
      break
    }
    hmat = matrix(eval(Hmat),3,3)
    new_theta <- old_theta - step_theta * tryCatch(solve(hmat),
                                                   error= function(e) {solve(hmat+diag(tol,nrow(hmat)))}) %*%jvec

    if (any(gev_positive(x=x, mu=new_theta[1],s=new_theta[2],k=new_theta[3])<0)) {
      del_theta <- new_theta-old_theta
      for (i in 1:length(alp)){
        alpha=alp[i]
        if(all(gev_positive(x=x, mu=old_theta[1]+alpha*del_theta[1],s=old_theta[2]+alpha*del_theta[2],k=old_theta[3]+alpha*del_theta[3])>0)){
          break
        }
      }
      new_theta = old_theta - alpha*solve(hmat)%*%jvec
    }
  
    old_theta = new_theta
  }
  return(list(initial = theta0, root = c(old_theta), step = niter, grad=jvec)) 
}


GEVnewtonRaphson_reg_test <- function (x, z, theta0, expr, step_theta=1, step_beta=1, maxiter=10, tol = 1e-6)
{
  z <- as.matrix(z)
  old_theta <- theta0 
  old_beta <- rep(0,ncol(z))
  niter <- 0
  alp <- 1-seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  for (i in 1:maxiter) {
    niter = niter + 1
    cat("niter:::", niter,  '\n')
    
    # theta update
    x_d = x- z %*% old_beta 
    mu=old_theta[1]; s=old_theta[2]; k=old_theta[3]
    fit_mle <- GEVnewtonRaphson_test(x = x_d, theta0 = c(mu,s,k), step_theta=step_theta, maxiter=5000)
    new_theta = fit_mle$root
    cat("new_theta:::",new_theta,'\n')
    # cat("theta gradient:::",fit_mle$grad,'\n')
    # v1 = lossfun(x = x - c(z%*%old_beta), new_theta[1], new_theta[2] ,new_theta[3])
    # cat("theta_update:::",v1,'\n"')
    cat("- - - - - - - - - - - - - - - - - - - - - - - -",'\n')
    # beta update
    for (j in 1:100){
    mu=c(new_theta[1]+z%*%old_beta); s=new_theta[2]; k=new_theta[3]
    jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
    cat("beta gradient:::",jvec_beta,'\n')
      if ( abs(max(jvec_beta)) < tol ) {   
        break
      }  
    hmat_beta = sum(eval(Hmat)[,1]) *t(z)%*%z
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



GEV_regfull_test <- function (x, z, theta0, beta0, expr=expr_reg, alpha=1, maxiter = 1000, tol = 1e-05) {
  old_theta <- c(theta0,beta0)
  niter <- 0
  alp <- 1-seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  
  for (i in 1:maxiter) {
    niter = niter + 1
    cat("niter:::", niter,  '\n')
    mu=c(old_theta[1]+z%*%old_theta[-c(1,2,3)]); s=old_theta[2]; k=old_theta[3]
    
    grad=apply(cbind(eval(Jaco),eval(Jaco)[,1]*z),2,sum)
    cat("grad:::",grad,'\n')
    cat("===============================================================",'\n')    
    
    if (max(abs(grad))<1e-07){
      break
    }
    
    hess = GEVhessian(x,z,mu,s,k)
    
    new_theta = old_theta - alpha*solve(hess)%*%grad
    
    if (any(gev_positive(x=x, mu=new_theta[1]+z%*%new_theta[-c(1,2,3)],s=new_theta[2],k=new_theta[3])<0)) {
      del_theta <- new_theta-old_theta
      for (i in 1:length(alp)){
        alpha=alp[i]
        if(all(gev_positive(x=x, mu=old_theta[1]+alpha*del_theta[1]+z%*%old_theta[-c(1,2,3)]+alpha*del_theta[-c(1,2,3)],
                            s=old_theta[2]+alpha*del_theta[2],k=old_theta[3]+alpha*del_theta[3])>0)){
          break
        }
      }
      new_theta = old_theta - alpha*solve(hess)%*%grad
    }
    cat("new_theta:::",new_theta,'\n')
    
    
    v = lossfun(x-z%*%new_theta[-c(1,2,3)],mu=new_theta[1],s=new_theta[2],k=new_theta[3])
    cat("loss_update:::", v, '\n')
    
    old_theta = new_theta
  }
  return(list(initial = theta0, root = c(old_theta), step = niter))
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





GEV_regfull <- function (x, z, theta0, beta0, expr=expr_reg, alpha=1, maxiter = 1000, tol = 1e-05) {
  old_theta <- c(theta0,beta0)
  niter <- 0
  alp <- 1-seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  
  for (i in 1:maxiter) {
    niter = niter + 1
    mu=c(old_theta[1]+z%*%old_theta[-c(1,2,3)]); s=old_theta[2]; k=old_theta[3]
    
    grad=apply(cbind(eval(Jaco),eval(Jaco)[,1]*z),2,sum)

    if (max(abs(grad))<1e-07){
      break
    }
    
    hess = GEVhessian(x,z,mu,s,k)
    
    new_theta = old_theta - alpha*solve(hess)%*%grad
    
    if (any(gev_positive(x=x, mu=new_theta[1]+z%*%new_theta[-c(1,2,3)],s=new_theta[2],k=new_theta[3])<0)) {
      del_theta <- new_theta-old_theta
      for (i in 1:length(alp)){
        alpha=alp[length(alp)-i]
        if(all(gev_positive(x=x, mu=old_theta[1]+alpha*del_theta[1]+z%*%old_theta[-c(1,2,3)]+alpha*del_theta[-c(1,2,3)],s=old_theta[2]+alpha*del_theta[2],k=old_theta[3]+alpha*del_theta[3])>0)){
          break
        }
      }
      new_theta = old_theta - alpha*solve(hess)%*%grad
    }
    
    
    # # loss가 초반에 너무 뛸경우 방지 -수정 필요
    # oldloss <- lossfun(x-z%*%old_theta[-c(1,2,3)],mu=old_theta[1],s=old_theta[2],k=old_theta[3]) 
    # newloss <- lossfun(x-z%*%new_theta[-c(1,2,3)],mu=new_theta[1],s=new_theta[2],k=new_theta[3]) 
    # if (newloss > oldloss || is.na(newloss)){
    #   new_theta = old_theta - 0.05*solve(hess)%*%grad
    # }
    # 
    # new_theta = old_theta - alpha*grad

    old_theta = new_theta
  }
  return(list(initial = theta0, root = c(old_theta), step = niter))
}
