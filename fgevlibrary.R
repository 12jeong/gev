f_density_gev=function(mu,sigma,k,x){
  sum(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k))
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



GEVnewtonRaphson <- function (x, theta0, alpha=1, expr = expr_mle, maxiter = 1000, tol = 1e-05) {
  
  if ( !all(1+theta0[3]*(x-theta0[1])/theta0[2]>0) ){
    warning("NaN produced in log by 'theta0'")
  }
  
  old_theta <- theta0
  niter <- 0
  alp <- seq(from=0,to=100,by=1)/100
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  for (i in 1:maxiter) {
    niter = niter + 1
    mu=old_theta[1]; sigma=old_theta[2]; k=old_theta[3]
    jvec = eval(Jaco)
    hmat = matrix(eval(Hmat),3,3)
    new_theta <- old_theta - alpha*solve(hmat)%*%jvec
    del_theta = new_theta - old_theta
    if ( new_theta[2] < 0 || !all(1+new_theta[3]*(x-new_theta[1])/new_theta[2]>0) ) {
      alpha = min(alp[alp>(tol-old_theta[2])/del_theta[2] && (old_theta[3]+alp*del_theta[3])*
                        (min(x)-(old_theta[1]-alp*del_theta[1])) > -(old_theta[2]+alp*del_theta[2])*tol])
      new_theta <- old_theta - alpha*solve(hmat)%*%jvec
      del_theta = new_theta - old_theta
    }
    if ( norm(as.matrix(del_theta)) < tol ) {    ### if (max(abs(del_theta)) < tol )
      break
    }
    old_theta = new_theta
    if (niter == maxiter) {
      warning("Maximum number of iterations 'maxiter' was reached.")
    }
  }
  return(list(initial = theta0, root = c(new_theta), step = niter, delta = del_theta))
}
