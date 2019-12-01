# 시뮬레이션 셋팅 만들기
rm(list=ls()); gc()
setwd("~/GitHub/gev")
source("./lib/pack.R")

ns = 20; nt = 40
library(metaSEM)
library(genlasso)


## setting2. unimodal
xyrange = c(-10,10)
nBS = 3
n = 30
x1 = seq(xyrange[1],xyrange[2],length.out=n)
x2 = seq(xyrange[1],xyrange[2],length.out=n)
mu = c(0,0)
sig = matrix(c(30,0,0,30),nrow=2)
fx = outer(x1, x2, function(x1,x2) {
  dmvnorm(cbind(x1,x2), mean=mu, sigma=sig)
})

set.seed(2)
s_ind = sort(sample(n^2,ns))
zval = rep(0, ns)
k = 1
for ( i in s_ind)
{
  zval[k] =fx[i]
  k = k + 1
}
s_col = s_ind%/%n + 1 ; s_col[s_ind%/%n == n] = n
s_row = s_ind%%n + 1 ; s_row[s_row == 0] = n

lat = x1[s_col]
long = x2[s_row]
plot(lat,long)

mu_tmp = zval*6000

par_mu = list()
set.seed(3)
for (s in 1:ns){
  m0 =  runif(1,90,100)
  par_mu[[s]] = mu_tmp[s] + seq(m0,m0+30,length.out = nt)
}

set.seed(4)
par_scale = rtruncnorm(ns, a=35, b=45, mean=40, sd=2)
par_shape = runif(ns,0.1,0.25)

set.seed(5)
xlist = list()
for ( s in 1:ns){
  xlist[[s]] = rgev(nt,loc=par_mu[[s]],scale=par_scale[s],shape=par_shape[s])
}

par(mfrow=c(2,2))
plot(par_mu[[1]])
plot(long,lapply(par_mu,function(x) x[1]))
plot(long,lapply(par_mu,function(x) x[2]))
plot(long,lapply(par_mu,function(x) x[40]))


mat_func = function(n) {
  m = matrix(0,n-2,n)   ## dmatrix : (n-2)*n 
  for (i in 1:(n-2)) {
    m[i,i] = 1
    m[i,i+1] = -2
    m[i,i+2] = 1
  }
  return(m)
}

dmatrix = mat_func(nt)

x_bsobj <- create.bspline.basis(range(long),
                                breaks=quantile(long,prob = seq(0, 1, length = 3)))
y_bsobj <- create.bspline.basis(range(lat),
                                breaks=quantile(lat,prob = seq(0, 1, length = 3)))

zlist <- list()
Dlist <- list()
for (i in 1:ns){
  xbs <- eval.basis(long[i],x_bsobj)
  ybs <- eval.basis(lat[i],y_bsobj)
  tensorbs <- do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs))  # row-wise kronecker product
  ztmp = list()
  for (ii in 1:nt){
    ztmp[[ii]] = tensorbs
  }
  zlist[[i]] = bdiagMat(ztmp) 
  # zlist[[i]] <- matrix(rep(tensorbs, nt),byrow=T,nrow=nt)
  Dlist[[i]] <- dmatrix %*% zlist[[i]]
}
dim(zlist[[1]])
newZ = do.call("rbind",zlist)
# View(head(newZ))
dim(newZ)



Y= unlist(xlist)
# lam_t = 1
# newD = matrix(rep(t(dmatrix),56),byrow=T,ncol=46) %*% newZ
newD = do.call("rbind",Dlist)
dim(newD)
# ?genlasso

out = genlasso(y=Y,X=newZ,D=newD)
summary(out$lambda)

est_beta = coef(out, lambda=0)$beta
Yhat = newZ%*%est_beta
summary(Yhat)

par(mfrow=c(3,3))
for ( s in 1:ns){
  plot(xlist[[s]])
  lines(Yhat[(1:nt)+nt*(s-1)],col="blue")
  # predict = newZ[(1:nt)+nt*(s-1),]%*%est_beta
  # lines(predict,col="blue")
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
  # if (shape > 1)
  #   d <- -(1e+6)
  else {
    nn <- length(x)
    xx <- 1 + shape * x
    xxpos <- xx[xx > 0 | is.na(xx)]
    scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
    d <- numeric(nn)
    d[xx > 0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) - 
      (1/shape + 1) * log(xxpos)
    d[xx <= 0 & !is.na(xx)] <- -(1e+6)
    # d[shape>1] <- -(1e+6)
  }
  return(d)
}

l2gev_t = function (x, tvec_t)
{  
  tvec = tvec_t
  # loglikelihood
  v = - sum(lgev(x, loc = tvec[1], scale = tvec[2], shape = tvec[3]))
  
  return(v)
}


ctr_list = list()
ctr_list$maxit = 1

Y2 = Y-Yhat
tout = list()
for (s in 1:ns){
  y = xlist[[s]]
  tvec_t = fgev(y)$estimate
  tout[[s]] = optim(tvec_t, l2gev_t, x=y, control=ctr_list, method="BFGS")$par  
}
unlist(tout)

ss = split.data.frame(Pr_46,Pr_46$stnlds)  # stnlds로 dataframe 쪼개서 list에 분배
xlist = lapply(ss,"[[","pr")               # 강수량(pr) 변수로만 이루어진 list 생성

Ymat =matrix(Yhat,nrow=ns,byrow=T)
mu_list = list()
for ( s in 1:ns){
  mu_list[[s]] = tout[[s]][1]+Ymat[s,]
}

ll = 0 
for (s in 1:ns){
  y = xlist[[s]]
  mu = mu_list[[s]]
  sc = tout[[s]][2]
  sh = tout[[s]][3]
  ll = ll - sum(lgev(x=y,loc=mu,scale=sc,shape=sh))
}

df1 = length(which(abs(newD%*%est_beta) > 1e-04)) + ns
df2 = 0
df3 = 2*ns
DF = df1+df2+df3

AIC = 2*ll + 2*DF
AIC

