rm(list=ls()); gc()
setwd("~/GitHub/gev")
source("./lib/pack.R")
load("~/GITHUB/gev/kma_data/Pr_46.RData")
# install.packages("genlasso")
# install.packages("metaSEM")
library(metaSEM)
library(genlasso)

set.seed(119)
Pr_46 = Pr_46 %>% filter(stnlds %in% sample(unique(Pr_46$stnlds),27))
# Pr_46 = Pr_46 %>% filter(obsyear < 2010 )
ns = length(unique(Pr_46$stnlds)) ;ns
nt = length(unique(Pr_46$obsyear)) ;nt

mat_func = function(n) {
  m = matrix(0,n-2,n)   ## dmatrix : (n-2)*n 
  for (i in 1:(n-2)) {
    m[i,i] = 1
    m[i,i+1] = -2
    m[i,i+2] = 1
  }
  return(m)
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

dmatrix = mat_func(nt)

# Y= Pr_46$pr[Pr_46$stnlds==90]
# out = genlasso(y=Y,X=diag(1,nt),D=dmatrix)
# Y2 <-  Y - coef(out, lambda=100)$beta
# plot(Y); points(coef(out, lambda=100)$beta,col="blue")
# tvec_t = fgev(Y)$estimate[-1]

x_bsobj <- create.bspline.basis(range(Pr_46$long),
                                breaks=quantile(Pr_46$long,prob = seq(0, 1, length = 3)))
y_bsobj <- create.bspline.basis(range(Pr_46$lat),
                                breaks=quantile(Pr_46$lat,prob = seq(0, 1, length = 3)))

zlist <- list()
Dlist <- list()
for (i in 1:ns){
  xbs <- eval.basis(unique(Pr_46$long)[i],x_bsobj)
  ybs <- eval.basis(unique(Pr_46$lat)[i],y_bsobj)
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


Y= Pr_46$pr
# lam_t = 1
# newD = matrix(rep(t(dmatrix),56),byrow=T,ncol=46) %*% newZ
newD = do.call("rbind",Dlist)
dim(newD)
# ?genlasso

out = genlasso(y=Y,X=newZ,D=newD,maxsteps=5500)
summary(out$lambda)

est_beta = coef(out, lambda=5500)$beta
Yhat = newZ%*%est_beta
summary(Yhat )

par(mfrow=c(3,3))
for ( s in 1:ns){
  plot(Y[(1:nt)+nt*(s-1)], main=paste("s:",s))
  lines(Yhat[(1:nt)+nt*(s-1)],col="blue")
  # predict = newZ[(1:nt)+nt*(s-1),]%*%est_beta
  # lines(predict,col="blue")
}

ctr_list = list()
ctr_list$maxit = 1

Y2 = Y-Yhat
tout = list()
for (s in 1:ns){
  y = Y2[(1:nt)+nt*(s-1)]
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

for (s in 1:ns){
  plot(xlist[[s]])
  lines(mu_list[[s]],col="red")
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
ll; AIC
