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

# Finding MLE for stationary
GEVnewtonRaphson_step1 <- function (x, theta0, step_theta=1, expr = expr_mle, tol = 1e-06) {
  
  if ( !all(1+theta0[3]*(x-theta0[1])/theta0[2]>0) ){
    warning("NaN produced in log by 'theta0'")
  }
  old_theta <- theta0
  # niter <- 0
  alp <- 1-seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  
    # niter = niter + 1
    mu=old_theta[1]; s=old_theta[2]; k=old_theta[3]
    jvec = eval(Jaco)
    # if ( abs(max(jvec)) < tol ) {    
    #   break
    # }
    hmat = matrix(eval(Hmat),3,3)
    
    # new_theta <- old_theta - step_theta * solve(hmat) %*% jvec
    ## debuging for computationally singular
    new_theta <- old_theta - step_theta * tryCatch(solve(hmat),
                                                   error= function(e) {solve(hmat+diag(tol,nrow(hmat)))}) %*%jvec
    
    # if (any(gev_positive(x=x, mu=new_theta[1],s=new_theta[2],k=new_theta[3])<0)) {
    #   del_theta <- new_theta-old_theta
    #   for (i in 1:length(alp)){
    #     alpha=alp[i]
    #     if(all(gev_positive(x=x, mu=old_theta[1]+alpha*del_theta[1],s=old_theta[2]+alpha*del_theta[2],k=old_theta[3]+alpha*del_theta[3])>0)){
    #       break
    #     }
    #   }
    #   new_theta = old_theta - alpha*solve(hmat)%*%jvec
    # }
    

    oldloss <- lossfun(x,mu=old_theta[1],s=old_theta[2],k=old_theta[3])
    newloss <- tryCatch(lossfun(x,mu=new_theta[1],s=new_theta[2],k=new_theta[3]),warning=function(e){})
    
    # loss가 초반에 너무 뛸경우 방지 -수정 필요
    # if (newloss > oldloss || is.null(newloss)){
    #   alpha = 0.05
    #   new_theta = old_theta - alpha*solve(hmat)%*%jvec
    #   newloss <- tryCatch(lossfun(x,mu=new_theta[1],s=new_theta[2],k=new_theta[3]),warning=function(e){})
    # }
    
  
  return(list(initial = theta0, root = c(new_theta), oldloss = oldloss, newloss = newloss , hess=hmat))
}



# Finding MLE for nonstationary - fullhessian method
GEV_regfull_step1 <- function (x, z, theta0, beta0, expr=expr_reg, alpha=1, tol = 1e-05) {
  old_theta <- c(theta0,beta0)
  # niter <- 0
  alp <- 1-seq(from=0,to=100,by=1)/100
  
  Jaco <- expr$Jaco
  Hmat <- expr$Hmat
  
    mu=c(old_theta[1]+z%*%old_theta[-c(1,2,3)]); s=old_theta[2]; k=old_theta[3]
    
    grad=apply(cbind(eval(Jaco),eval(Jaco)[,1]*z),2,sum)
    hess = GEVhessian(x,z,mu,s,k)
    new_theta = old_theta - alpha*solve(hess)%*%grad
    
    # if (any(gev_positive(x=x, mu=new_theta[1]+z%*%new_theta[-c(1,2,3)],s=new_theta[2],k=new_theta[3])<0)) {
    #   del_theta <- new_theta-old_theta
    #   for (i in 1:length(alp)){
    #     alpha=alp[i]
    #     if(all(gev_positive(x=x, mu=old_theta[1]+alpha*del_theta[1]+z%*%old_theta[-c(1,2,3)]+alpha*del_theta[-c(1,2,3)],s=old_theta[2]+alpha*del_theta[2],k=old_theta[3]+alpha*del_theta[3])>0)){
    #       break
    #     }
    #   }
    #   new_theta = old_theta - alpha*solve(hess)%*%grad
    # }
    
    # loss가 초반에 너무 뛸경우 방지 -수정 필요
    oldloss <- lossfun(x-z%*%old_theta[-c(1,2,3)],mu=old_theta[1],s=old_theta[2],k=old_theta[3])
    newloss <- tryCatch(lossfun(x-z%*%new_theta[-c(1,2,3)],mu=new_theta[1],s=new_theta[2],k=new_theta[3]),warning=function(e){})
    # if (newloss > oldloss || is.na(newloss)){
    #   new_theta = old_theta - 0.05*solve(hess)%*%grad
    #   newloss <- tryCatch(lossfun(x,mu=new_theta[1],s=new_theta[2],k=new_theta[3]),warning=function(e){})
    # }

    
  
  return(list(initial = theta0, root = c(new_theta), oldloss = oldloss, newloss = newloss ,hess=hess))
}




