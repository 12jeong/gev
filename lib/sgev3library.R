gevreg_3m = function(xlist, zlist, lambda = c(0,0,0), Om=NULL)
{
  lambda = c(lambda)
  p = ncol(zlist[[1]])
  ns = length(xlist)
  tvec = rep(0,3*(p+1)) 
  Om = Om
  
  l2gev_m = function (tvec, lambda, xlist, zlist, Om)
  {
    ns = length(xlist)
    v1 = 0
    loc.vec.reg = tvec[3+(1:p)]
    sc.vec.reg = tvec[(3+p)+(1:p)]
    sh.vec.reg = tvec[(3+2*p)+(1:p)]
    for ( i in 1: length(xlist))
    {
      x = xlist[[i]]
      z = zlist[[i]]
      loc.vec = tvec[1] + drop(z%*%loc.vec.reg)
      sc.vec = tvec[2] + drop(z%*%sc.vec.reg)
      sh.vec = tvec[3] + drop(z%*%sh.vec.reg)
      v1 = v1 - sum(lgev(x, loc = loc.vec, 
                         scale = exp(sc.vec), shape = sh.vec))   # loss : +(-loglikelihood)
    }
    # regularization
      v2 = lambda[1]*t(loc.vec.reg)%*%Om%*%loc.vec.reg
      v3 = lambda[2]*t(sc.vec.reg)%*%Om%*%sc.vec.reg
      v4 = lambda[3]*t(sh.vec.reg)%*%Om%*%sh.vec.reg
      v = v1 + v2 + v3+ v4


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
  
  tmp_loc = 0; tmp_sc = 0; tmp_sh = 0
  for ( i in 1:ns)
  {
    x = xlist[[i]]
    fit = fgev(x)
    tmp_loc = tmp_loc + fit$estimate[1]
    tmp_sc = tmp_sc + fit$estimate[2]
    tmp_sh = tmp_sh + fit$estimate[3]
  }
  tmp_loc = tmp_loc/ns
  tmp_sc = tmp_sc/ns
  tmp_sh = tmp_sh/ns
  tvec[1] = tmp_loc
  tvec[2] = log(tmp_sc)
  tvec[3] = tmp_sh
  
  fit = optim(tvec, l2gev_m, lambda = lambda,method = "BFGS", 
              xlist = xlist, zlist = zlist, Om = Om, control = list(maxit=10000))
  return(fit) 
}


