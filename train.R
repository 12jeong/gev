rm( list = ls()); gc()
setwd("C:\\Users\\UOS\\Documents\\GITHUB\\gev")
source("sgevlibrary.R")
load("kma_data\\Pr_46.Rdata")
library(fda)
library(dplyr)

ylat = unique(Pr_46$lat)
xlong = unique(Pr_46$long)

plot(xlong,ylat)
summary(Pr_46$lat)
summary(Pr_46$long)
ns = length(unique(Pr_46$stnlds)) # 56개

set.seed(3)
kcluster = kmeans(data.frame(ylat,xlong), 5)
plot(xlong,ylat,col=kcluster$cluster)
kcluster$size

train = Pr_46 %>% filter(obsyear < 2013)
train$k = rep(kcluster$cluster,each=nrow(train)/ns)
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
dim(matfunc(ns))

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

dim(zlist[[1]]) # (frist)stnlds 2D-splines tensor, nbasis = 5 x 5

# Omega matrix
Fmat = kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat = kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat = kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om = Fmat+2*Gmat+Hmat

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3






######### k-folds cross validation ##########
lambdaset = seq(0,1,length=11)
result_train = list()

for (i in c(1:5)){
  tmp = list()
  for ( j in c(1:length(lambdaset))){
    ns = sum(kcluster$size[-i])
    mat = matfunc(ns)
    tmp[[j]] <- gevreg_m(xlist[which(kcluster$cluster!=i)], zlist[which(kcluster$cluster!=i)], 
                         lambda = lambdaset[[j]], lambda2=1, Om=Om, mat=mat, method="B-spline")
  }
  result_train[[i]] <- tmp
}

result_train = result

result_validation = list()
for (i in c(1:5)){
  tmp = list()
  for ( j in c(1:kcluster$size[i])){
    tmp[[j]] <- fgev(xlist[[which(kcluster$cluster==i)[j]]])$estimate
  }
  result_validation[[i]] <- tmp
}


loclambda = list()
for (whichlambda in c(1:length(lambdaset))){
  loclust = list()
  for(clust in c(1:5)) {
    loc = c()
    for ( i in c(1:kcluster$size[clust])){
      j = which(kcluster$cluster==clust)[i]
      loc[i] =  (mean(result_train[[clust]][[whichlambda]][(3*c(1:sum(kcluster$size[-clust])))-2]) + 
              zlist[[j]]%*%tail(result_train[[clust]][[whichlambda]],25))[1]
    }
    loclust[[clust]] = loc
  }
  loclambda[[whichlambda]] = loclust
}
loclambda



mleclust = list()
for(clust in c(1:5)) {
  mleclust[[clust]] = unlist(lapply(1:kcluster$size[clust], function(i) {result_validation[[clust]][[i]][1]}))
}

mseloc = c() 
for ( whichlambda in c(1:11)){
  mseloc[whichlambda] = mean(unlist(lapply(1:5,function(i) mean((loclambda[[whichlambda]][[i]]-mleclust[[i]])^2))))
}
plot(lambdaset,mseloc)

# save.image("train.RData")

s.lambda = lambdaset[which.min(mseloc)]
ns = length(unique(train$stnlds)) 
train_mle = gevreg_m(xlist, zlist, lambda = s.lambda , 
                     lambda2=1, Om=Om, mat=matfunc(ns) , method="B-spline")

test_ss = split.data.frame(test,test$stnlds)
test_xlist = lapply(test_ss,"[[","pr")

test_mle = list()
for ( i in c(1:8)) {
  test_mle[[i]] = fgev(test_xlist[[i]])$estimate
}
test_xlist
