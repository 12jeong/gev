gev_positive = function(x,mu,sigma,k){
  1+k*(x-mu)/sigma
}

f_density_gev=function(mu,sigma,k,x){
  log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k)
}


grad_mu_mat <- function(x,mu, sigma, k)
{
  y=(x-mu)/sigma
  z=1+k*y
  phi=z^(-1/k)
  psi=k+1-phi
  # partial derivariate
  grad_mu_mat =-psi/(sigma*z)
  return(grad_mu_mat)
}


lossfun = function(x,mu,sigma,k){
  v = log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k)  
  sum(v)
}

# derivatives for mles
Jaco1=Deriv(expression(f_density_gev(mu,sigma,k,x)),c("mu","sigma","k"))
Hmat1=Deriv(expression(f_density_gev(mu,sigma,k,x)),c("mu","sigma","k"),n=c(hessian=2))
expr_mle <- list (Jaco = Jaco1, Hmat = Hmat1)

# derivatives for regression 
logl=expression(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k))

Jaco2=Deriv(logl,c("mu","sigma","k"),combine="cbind")
Hmat2=Deriv(logl,c("mu","sigma","k"),n=c(hessian=2),combine="cbind")
expr_reg <- list (Jaco = Jaco2, Hmat = Hmat2)


# Finding MLE for stationary
GEVnewtonRaphson <- function (x, theta0, step_theta=1, expr = expr_mle, maxiter = 5000, tol = 1e-06) {
  
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
    if ( abs(max(jvec)) < tol ) {    
      break
    }
    hmat = matrix(eval(Hmat),3,3)
    
    # new_theta <- old_theta - step_theta * solve(hmat) %*% jvec
    ## debuging for computationally singular
    new_theta <- old_theta - step_theta * tryCatch(solve(hmat),
                                                   error= function(e) {solve(hmat+diag(tol,nrow(hmat)))}) %*%jvec
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
    fit_mle <- GEVnewtonRaphson(x = x_d, theta0 = c(mu,sigma,k), step_theta=step_theta, maxiter=5000)
    new_theta = fit_mle$root
    
    # beta update
    for (j in 1:5000){
      mu=c(new_theta[1]+z%*%old_beta); sigma=new_theta[2]; k=new_theta[3]
      jvec_beta = apply(eval(Jaco)[,1]*z,2,sum)
      if ( abs(max(jvec_beta)) < tol ) {   
        break
      }  
      hmat_beta = eval(Hmat)[1,1] *t(z)%*%z
      new_beta <- old_beta - step_beta *solve(hmat_beta)%*%jvec_beta
      old_beta <- new_beta
    }
    old_theta <- new_theta
  }
  return(list(step = niter, initial = theta0, root = c(old_theta,old_beta), grad=c(fit_mle$grad,jvec_beta)))
}

