rm( list = ls()); gc()
setwd("C:\\Users\\UOS\\Documents\\GITHUB\\gev")
source("sgevlibrary.R")
load("kma_data\\Pr_46.Rdata")
library(fda)
library(dplyr)

ylat = unique(Pr_46$lat)
xlong = unique(Pr_46$long)
plot(xlong,ylat)

ns = length(unique(Pr_46$stnlds)) # 56개

# 1973~2012 for train, 2013~2018 for test
train = Pr_46 %>% filter(obsyear < 2013)
test = Pr_46 %>% filter(obsyear >= 2013)

# design matrix for v3 (m0 statinary?)
matfunc = function(ns){
  matt = matrix(0,nrow=ns*(ns-1),ncol=ns)
  k = 1
  for (i in (1: (ns-1))){
    for (j in ((i+1) :ns)){
      matt[k,i] = 1
      matt[k,j] = -1
      k = k+1
    }  
  }
  mat = t(matt) %*% matt
  return(mat)
}
mat = matfunc(ns)
dim(mat)

ss = split.data.frame(train,train$stnlds)
xlist = lapply(ss,"[[","pr")
x_bsobj = create.bspline.basis(range(train$long),breaks=quantile(train$long,prob = seq(0, 1, length = 3)))
y_bsobj = create.bspline.basis(range(train$lat),breaks=quantile(train$lat,prob = seq(0, 1, length = 3)))

zlist = list()
for (i in 1:ns){
  xbs = eval.basis(ss[[i]]$long,x_bsobj)
  ybs = eval.basis(ss[[i]]$lat,y_bsobj)
  tensorbs = do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs)) 
  zlist[[i]] = tensorbs 
}

dim(zlist[[1]]) # (frist)stnlds 2D-splines tensor, nbasis = df x df

# Omega matrix
Fmat = kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat = kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat = kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om = Fmat+2*Gmat+Hmat

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

######### k-folds cross validation ##########
## 순서대로 
# ksize = 10 # number of K
# gsize = 40/ksize
# folds = rep(c(1:ksize),each=gsize)
# train$k = rep(folds,times=ns)

## 랜덤하게 
set.seed(517)
ksize = 5 # number of K
folds = sample(cut(1:40,breaks=ksize,labels=FALSE),40)
train$k = rep(folds,times=ns)

lambdaset = seq(0,1,length=11)
result_train = list()

s_time <- Sys.time()
for (i in c(1:5)){
  tmp = list()
  trainx = lapply(1:ns, function(k) xlist[[k]][folds!=i] )
  trainz = lapply(1:ns, function(k) zlist[[k]][1:(40-sum(folds==i)),])
  for ( j in c(1:length(lambdaset))){
    tmp[[j]] <- gevreg_m(xlist=trainx,zlist=trainz, 
                         lambda = lambdaset[j], lambda2=1, Om=Om, mat=mat, method="B-spline")
      cat("i=",i," j=",j,"\n" )
  }
  result_train[[i]] <- tmp
}
e_time <- Sys.time()
e_time-s_time

# tval <- gevreg_m(xlist=trainx,zlist=trainz,
#                 lambda = lambdaset[1], lambda2=1, Om=Om, mat=mat, method="B-spline")
# 
# tmp = list()
# for ( j in 1:11){
#   tmp[[j]] <- tval
# }
# for ( i in 1:5){
#   result_train[[i]] <- tmp
# }



##########################################

lossfun = function(x,mu,s,k){
  v = log(s)+(1+1/k)*log(1+k*(x-mu)/s)+(1+k*(x-mu)/s)^(-1/k)  
  sum(v)
}

i=1 # validation set
l=1 # lambda 
s=1 # location

# load("kfolds=random2018_knots=4.RData")
result_validation = list()
for (i in 1:5){
  lset = c()
  for( l in c(1:length(lambdaset))){
    v = c()
    est_z = tail(result_train[[i]][[l]],dim(zlist[[1]])[2]) 
    for( s in 1:ns){
      est_s = result_train[[i]][[l]][(1:3)+(3*s-3)]
      v[s] = lossfun(xlist[[s]][c(1:8)*i],mu=c(est_s[1]+zlist[[s]][1,]%*%est_z),s=est_s[2],k=est_s[3])
    }
    if (any(is.na(v))) cat("NaN produced in:","validation set",i,"lambda",lambdaset[l],"location",which(is.na(v)),"\n")
    lset[l] <- sum(v)
  }
  result_validation[[i]] <- lset
}

result_validation

# save(result_train,result_validation,file="kfolds=random527_knots=4.RData")

# result_validation = list()
# for (i in 1:5){
#   lset = c()
#   for( l in c(1:length(lambdaset))){
#     v = c()
#     for( s in 1:ns){
#       v[s] = 48.01299
#     }
#     lset[l] <- sum(v)
#   }
#   result_validation[[i]] <- lset
# }

# lambda1 min 찾기 
validation.mat <- do.call('rbind',result_validation)
par(mfrow=c(1,1))
plot(lambdaset,colMeans(validation.mat,na.rm=TRUE))
lambda.min <- lambdaset[which.min(colMeans(validation.mat,na.rm=TRUE))]


# 최종 적합 - train set
test_est <- gevreg_m(xlist,zlist,lambda = lambda.min , lambda2=1, Om=Om, mat=mat, method="B-spline")

# test likelihood
test_ss = split.data.frame(test,test$stnlds)
test_xlist = lapply(test_ss,"[[","pr")

v=c()
est_z = tail(test_est,dim(zlist[[1]])) 
for (s in 1:ns){
  est_s = test_est[(1:3)+(3*s-3)]
  v[s] = lossfun(test_xlist[[s]],mu=c(est_s[1]+zlist[[s]][1,]%*%est_z),s=est_s[2],k=est_s[3])
}
if (any(is.na(v))) cat("NaN produced in:","location",which(is.na(v)),"\n")
result_test <- sum(v)
