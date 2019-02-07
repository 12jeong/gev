gev_positive = function(x,mu,s,k){
  1+k*(x-mu)/s
}

f_density_gev=function(mu,sigma,k,x){
  sum(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k))
}

lossfun = function(x,mu,s,k){
  v = log(s)+(1+1/k)*log(1+k*(x-mu)/s)+(1+k*(x-mu)/s)^(-1/k)  
  sum(v)
}

# derivatives for mles
Jaco1=Deriv(expression(f_density_gev(mu,s,k,x)),c("mu","s","k"))
Hmat1=Deriv(expression(f_density_gev(mu,s,k,x)),c("mu","s","k"),n=c(hessian=2))
expr_mle <- list (Jaco = Jaco1, Hmat = Hmat1)

# derivatives for regression 
logl=expression(log(s)+(1+1/k)*log(1+k*(x-mu)/s)+(1+k*(x-mu)/s)^(-1/k))
Jaco2=Deriv(logl,c("mu","s","k"),combine="cbind")
Hmat2=Deriv(logl,c("mu","s","k"),n=c(hessian=2),combine="cbind")
expr_reg <- list (Jaco = Jaco2, Hmat = Hmat2)

# gradient and hessian for full update in newtonraphson + reg
GEVgradient <- function(x,z,mu,s,k){
  apply(cbind(eval(expr_reg$Jaco),eval(expr_reg$Jaco)[,1]*z),2,sum)
}

GEVhessian <- function(x,z,mu,s,k){
  evalh = eval(expr_reg$Hmat)
  h1=matrix(apply(evalh,2,sum),3,3)
  h2=rbind(apply(evalh[,1]*z,2,sum), apply(evalh[,2]*z,2,sum), apply(evalh[,3]*z,2,sum))
  h3=t(h2)
  # h4=apply((eval(expr_reg$Hmat)[,1]*z%*%t(z)),2,sum)
  h4=matrix(0,dim(z)[2],dim(z)[2])
  for (i in 1:dim(z)[1]){
    h4row=evalh[,1][i] * z[i,]%*%t(z[i,])
    h4=h4+h4row
  }
  rbind(cbind(h1,h2),cbind(h3,h4)) 
}

# Finding MLE for stationary
GEVnewtonRaphson <- function (x, theta0, step_theta=1, expr = expr_mle, maxiter = 5000, tol = 1e-06) {
  
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
    if ( abs(max(jvec)) < tol ) {    
      break
    }
    hmat = matrix(eval(Hmat),3,3)
    
    # new_theta <- old_theta - step_theta * solve(hmat) %*% jvec
    ## debuging for computationally singular
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
    
    
    # ## debuging for log1p NaNs produced
    # if ( new_theta[2] < 0 || any(1+new_theta[3]*(x-new_theta[1])/new_theta[2]<0) ) {
    #   del_theta <- new_theta - old_theta
    #   step_theta <- min(alp[alp>(tol-old_theta[2])/del_theta[2] &&
    #                           (old_theta[3]+alp*del_theta[3])*(min(x)-(old_theta[1]-alp*del_theta[1])) > -(old_theta[2]+alp*del_theta[2])*tol])
    #   new_theta <- old_theta - step_theta * solve(hmat) %*% jvec
  }
  return(list(initial = theta0, root = c(old_theta), step = niter, grad=jvec))
}



# Finding MLE for nonstationary
GEVnewtonRaphson_reg <- function (x, z, theta0, expr=expr_reg, step_theta=1, step_beta=1, maxiter=10, tol = 1e-6)
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
    
    # theta update
    x_d = x- z %*% old_beta 
    mu=old_theta[1]; s=old_theta[2]; k=old_theta[3]
    fit_mle <- GEVnewtonRaphson(x = x_d, theta0 = c(mu,s,k), step_theta=step_theta, maxiter=5000)
    new_theta = fit_mle$root
    
    # beta update
    for (j in 1:100){
      mu=c(new_theta[1]+z%*%old_beta); s=new_theta[2]; k=new_theta[3]
      jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
      if ( abs(max(jvec_beta)) < tol ) {   
        break
      }  
      hmat_beta = sum(eval(Hmat)[,1]) *t(z)%*%z
      new_beta <- old_beta - step_beta *solve(hmat_beta)%*%jvec_beta
      old_beta <- new_beta
    }
    old_theta <- new_theta
  }
  return(list(step = niter, initial = theta0, root = c(old_theta,old_beta), grad=c(fit_mle$grad,jvec_beta)))
}



# Finding MLE for nonstationary - fullhessian method
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
        alpha=alp[i]
        if(all(gev_positive(x=x, mu=old_theta[1]+alpha*del_theta[1]+z%*%old_theta[-c(1,2,3)]+alpha*del_theta[-c(1,2,3)],s=old_theta[2]+alpha*del_theta[2],k=old_theta[3]+alpha*del_theta[3])>0)){
          break
        }
      }
      new_theta = old_theta - alpha*solve(hess)%*%grad
    }
    
    
    
    
    old_theta = new_theta
  
  }
  return(list(initial = theta0, root = c(old_theta), step = niter))
}




