rm(list=ls()); gc()
setwd("~/GitHub/gev")
source("./lib/pack.R")
load("C:/Users/UOS/Downloads/trend_filtering_2/result_lamt3.RData")

ns = length(xlist)
x_bsobj <- create.bspline.basis(range(Pr_46$long),breaks=quantile(Pr_46$long,prob = seq(0, 1, length = 3)))
y_bsobj <- create.bspline.basis(range(Pr_46$lat),breaks=quantile(Pr_46$lat,prob = seq(0, 1, length = 3)))
zlist <- list()
for (i in 1:length(xlist)){
  xbs <- eval.basis(ss[[i]]$long[1],x_bsobj)
  ybs <- eval.basis(ss[[i]]$lat[1],y_bsobj)
  tensorbs <- do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs))  # row-wise kronecker product
  zlist[[i]] <- tensorbs 
}
Z = do.call("rbind",zlist)

# 2-D splines penlaty matrix
Fmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om <- Fmat+Gmat+2*Hmat

Mu = do.call("rbind",old_mu_list)
beta_list = list()
for ( i in 1:nt){
  beta_list[[i]] = lm(Mu[,i]~Z-1)$coefficients
}

mu_hat_list = list()
for ( s in 1:ns){
  mu_vec = c()
  for ( i in 1:nt){
  mu_vec[i] = (Z%*%beta_list[[i]])[s]
  }
  mu_hat_list[[s]] = mu_vec
}

par(mfrow=c(3,3))
for ( s in 1:ns){
  plot(xlist[[s]])
  lines(old_mu_list[[s]],col="green") # trend filtering만 한 것
  lines(mu_hat_list[[s]],col="blue") # Z%*%beta_t 로 추정한것
}

par(mfrow=c(4,4))
for ( i in 1:nt){
  plot(unique(Pr_46$lat),Z%*%beta_list[[i]])
  points(unique(Pr_46$lat),Mu[,i],col="gray")
}

par(mfrow=c(4,4))
for ( i in 1:nt){
  plot(unique(Pr_46$lat),lapply(xlist,function(x)x[i]),col="gray",main=paste0("t:",i)) #data_t
  points(unique(Pr_46$lat),Z%*%beta_list[[i]],col="red") # Z%*%beta_t
}



library(fields)
par(mfrow=c(4,4))
for ( i in 1:nt){
  Tpsfit = Tps(x=data.frame(unique(Pr_46$lat),unique(Pr_46$long)),Y=Mu[,i],m=3)
  fit1 = predict(Tpsfit,lambda=12)
  plot(unique(Pr_46$lat),lapply(xlist,function(x)x[i]),col="gray",main=paste0("t:",i)) #data_t
  # points(unique(Pr_46$lat),fit1,col="blue") # Z%*%beta_t
  # plot(unique(Pr_46$lat),Z%*%beta_list[[i]],col="gray") # Z%*%beta_t
  points(unique(Pr_46$lat),fit1,col="blue") # Z%*%beta_t
}
