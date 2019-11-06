rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")
load("./sgev3_train1003.RData")

result_df = as.data.frame(expand.grid(lam_set,lam_set,lam_set))

df.tps = c()
for (i in 1:length(lam_set)){
  lam_s = lam_set[i]
  Hatmat= Z%*%solve(t(Z)%*%Z +lam_s*Om+diag(1e-08,nrow(Om)))%*%t(Z)
  df.tps[i] =  sum(diag(Hatmat))
}
DF = rowSums(as.matrix(expand.grid( df.tps, df.tps, df.tps)))

nll = c()
for (i in 1:length(result_list)){
  tvec=result_list[[i]]
  
  loc.vec.reg = tvec[3+(1:p)]
  sc.vec.reg = tvec[(3+p)+(1:p)]
  sh.vec.reg = tvec[(3+2*p)+(1:p)]
  loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
  sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
  sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
  
  ll = 0
  for (s in 1:length(xlist)){
  ll = ll - sum(dgev(x=xlist[[s]],loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s], log=TRUE))
  }
  nll[i] = ll
}

AIC = 2*nll + 2*DF
BIC = 2*nll + log(length(train_46))*DF
which.min(AIC)
which.min(BIC)
result_df[541,]

par(mfrow=c(1,1))
plot(AIC,type="l")
plot(BIC,type="l")

tvec = gevreg_3m(xlist=xlist, zlist=zlist, lambda = c(0,0,0), Om=Om)$par
loc.vec.reg = tvec[3+(1:p)]
sc.vec.reg = tvec[(3+p)+(1:p)]
sh.vec.reg = tvec[(3+2*p)+(1:p)]
loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)

point.est = unlist(lapply(xlist,function(x) fgev(x)$estimate))
point.loc = point.est[3*(1:ns)-2]
point.sc =  point.est[3*(1:ns)-1]
point.sh =  point.est[3*(1:ns)]

par(mfrow=c(1,3))
plot(point.loc,loc.vec)
plot(point.sc,sc.vec)
plot(point.sh,sh.vec)

plot(unique(train_46$lat),point.loc); points(unique(train_46$lat),loc.vec,col="red")
plot(unique(train_46$lat),point.sc); points(unique(train_46$lat),sc.vec,col="red")
plot(unique(train_46$lat),point.sh); points(unique(train_46$lat),sh.vec,col="red")

