gevreg_m = function(xlist, zlist, lambda = 0, lambda2=0, Om=NULL, mat=NULL,
                    method = c('linear', 'B-spline'))
{
  p = ncol(zlist[[1]])
  tvec = rep(0, ns*3 + p) # ns*mu,sigma,kappa(3), basis(p)
  ns = length(xlist)
  Om = Om
  mat = mat
  
  l2gev_m = function (tvec, lambda, xlist, zlist, Om, lambda2, mat)
  {
    ns = length(xlist)
    v1 = 0
    for ( i in 1: length(xlist))
    {
      x = xlist[[i]]
      z = zlist[[i]]
      loc.vec.reg = drop(z%*%tail(tvec, p))
      loc.vec = tvec[3*(i-1)+1] + loc.vec.reg
      sc.vec = tvec[3*(i-1)+2]
      sh.vec = tvec[3*(i-1)+3]
      v1 = v1 - sum(lgev(x, loc = loc.vec, 
                         scale = sc.vec, shape = sh.vec))   # loss : +(-loglikelihood)
    }
    # loglikelihood
    
    # regularization
    if (method == 'B-spline') {
      v2 = lambda*t(tail(tvec,p))%*%Om%*%tail(tvec,p)
      v3 = lambda2 * t(tvec[(0:(ns-1))*3+1] ) %*% mat %*% tvec[(0:(ns-1))*3+1] 
      v = v1 + v2 + v3
    }
    
    if (method == 'linear') {
      v2 = lambda*sum(tail(tvec,p)^2)
      v = v1 + v2
    }
    
    return(v)
  }
  
  lgev = function (x, loc = 0, scale = 1, shape = 0) 
  {
    if (min(scale) <= 0) 
      return( - 1e+6)
    if (length(shape) != 1) 
      stop("invalid shape")
    x <- (x - loc)/scale
    if (shape == 0) 
      d <- log(1/scale) - x - exp(-x)
    else {
      nn <- length(x)
      xx <- 1 + shape * x
      xxpos <- xx[xx > 0 | is.na(xx)]
      scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
      d <- numeric(nn)
      d[xx > 0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) - 
        (1/shape + 1) * log(xxpos)
      d[xx <= 0 & !is.na(xx)] <- -(1e+6)
    }
    return(d)
  }
  
  
  for ( i in 1:ns)
  {
    x = xlist[[i]]
    start <- list()
    start$scale <- sqrt(6 * var(x))/pi
    start$loc <- mean(x) - 0.58 * start$scale
    tvec[3*(i-1)+1] = start$loc
    tvec[3*(i-1)+2] = start$scale
  }
  
  
  return( optim(tvec, l2gev_m, lambda = lambda, lambda2=lambda2, Om=Om, mat=mat,
                method = "BFGS", 
                xlist = xlist,
                zlist = zlist)$par) 
}




lossfun = function (x, loc = 0, scale = 1, shape = 0) 
{
  if (min(scale) <= 0) 
    return( - 1e+6)
  if (length(shape) != 1) 
    stop("invalid shape")
  x <- (x - loc)/scale
  if (shape == 0) 
    d <- log(1/scale) - x - exp(-x)
  else {
    nn <- length(x)
    xx <- 1 + shape * x
    xxpos <- xx[xx > 0 | is.na(xx)]
    scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
    d <- numeric(nn)
    d[xx > 0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) - 
      (1/shape + 1) * log(xxpos)
    d[xx <= 0 & !is.na(xx)] <- -(1e+6)
  }
  return(d)
}
