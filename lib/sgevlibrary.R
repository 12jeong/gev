gevreg_m = function(xlist, zlist, lambda = 0, Om=NULL,
                    method = c('linear', 'B-spline'))
{
  p = ncol(zlist[[1]])
  tvec = rep(0, 1 + ns*2 + p) # m_0, ns*(sigma,kappa), basis(p)
  ns = length(xlist)
  Om = Om
  
  l2gev_m = function (tvec, lambda, xlist, zlist, Om)
  {
    ns = length(xlist)
    v1 = 0
    for ( i in 1: length(xlist))
    {
      x = xlist[[i]]
      z = zlist[[i]]
      loc.vec.reg = drop(z%*%tail(tvec, p))
      loc.vec = tvec[1] + loc.vec.reg
      sc.vec = tvec[2*i]
      sh.vec = tvec[2*i+1]
      v1 = v1 - sum(lgev(x, loc = loc.vec, 
                          scale = sc.vec, shape = sh.vec))   # loss : +(-loglikelihood)
    }
    # loglikelihood
    # regularization
    if (method == 'B-spline') {
      v2 = lambda*t(tail(tvec,p))%*%Om%*%tail(tvec,p)
      v = v1 + v2
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
  
  tmp_loc = 0
  for ( i in 1:ns)
  {
    x = xlist[[i]]
    fit = fgev(x)
    tmp_loc = tmp_loc + fit$estimate[1]
    tvec[2*i] = fit$estimate[2]
    tvec[2*i+1] = fit$estimate[3]
  }
  tmp_loc = tmp_loc/ns
  tvec[1] = tmp_loc
  
  fit = optim(tvec, l2gev_m, lambda = lambda,method = "BFGS", 
              xlist = xlist, zlist = zlist, Om = Om)
  return(fit) 
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
