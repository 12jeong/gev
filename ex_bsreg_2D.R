rm(list=ls())
library(fda)
library(evd)
library(Deriv)
setwd("~/GITHUB/gev")
source("./sgevlibrary.R")
load("./kma_data/Pr_46.RData")

ss <- split.data.frame(Pr_46,Pr_46$stnlds)  # stnlds로 dataframe 쪼개서 list에 분배
xlist <- lapply(ss,"[[","pr")               # 강수량(pr) 변수로만 이루어진 list 생성
ns <- length(unique(Pr_46$stnlds))

# create.bspline.basis : nbasis = norder + length(breaks) -2
# ex1 : order=4, quantile knots, nbasis=7 
# ex2 : order=4, (0,0.5,1) knots, nbasis=5

# x : long (경도) , y : lat (위도) 
x_bsobj <- create.bspline.basis(range(Pr_46$long),breaks=quantile(Pr_46$long,prob = seq(0, 1, length = 3)))
y_bsobj <- create.bspline.basis(range(Pr_46$lat),breaks=quantile(Pr_46$lat,prob = seq(0, 1, length = 3)))

zlist <- list()
for (i in 1:ns){
  xbs <- eval.basis(ss[[i]]$long,x_bsobj)
  ybs <- eval.basis(ss[[i]]$lat,y_bsobj)
  tensorbs <- do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs))  # row-wise kronecker product
  zlist[[i]] <- tensorbs 
}
dim(zlist[[1]]) # (frist)stnlds 2D-splines tensor, nbasis = 5 x 5

# 2-D splines penlaty matrix
Fmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om <- Fmat+Gmat+2*Hmat

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

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





## example 
result1 <- round(gevreg_m(xlist, zlist, lambda = 0, lambda2=1, Om=Om, mat=mat, method="B-spline"),3)
result2 <- round(gevreg_m(xlist, zlist, lambda = 0.2, lambda2=1, Om=Om, mat=mat, method="B-spline"),3)
result3 <- round(gevreg_m(xlist, zlist, lambda = 0.4, lambda2=1, Om=Om, mat=mat, method="B-spline"),3)
result4 <- round(gevreg_m(xlist, zlist, lambda = 0.6, lambda2=1, Om=Om, mat=mat, method="B-spline"),3)
result5 <- round(gevreg_m(xlist, zlist, lambda = 0.8, lambda2=1, Om=Om, mat=mat, method="B-spline"),3)
result6 <- round(gevreg_m(xlist, zlist, lambda = 1, lambda2=1, Om=Om, mat=mat, method="B-spline"),3)

p <- ncol(zlist[[1]])
x <- unique(Pr_46$long)
y <- unique(Pr_46$lat)
xbss <- eval.basis(x,x_bsobj)
ybss <- eval.basis(y,y_bsobj)
tensorbss <- do.call('cbind', lapply(1:ncol(xbss), function(i) xbss[, i] * ybss)) 
z1 <- tensorbss %*% tail(result1,p)
z2 <- tensorbss %*% tail(result2,p)
z3 <- tensorbss %*% tail(result3,p)
z4 <- tensorbss %*% tail(result4,p)
z5 <- tensorbss %*% tail(result5,p)
z6 <- tensorbss %*% tail(result6,p)

# bsmat <- data.frame(x,y,z1,z2,z3)
save(result1,result2,result3,result4,result5,result6,file="ex_bsreg_2D.RData")
# load("ex_bsreg_2D.RData")
library(scatterplot3d)
par(mfrow=c(2,3))
scatterplot3d(x,y,z1,scale.y=0.6,angle=10,main=0)
scatterplot3d(x,y,z2,scale.y=0.6,angle=10,main=0.2)
scatterplot3d(x,y,z3,scale.y=0.6,angle=10,main=0.4)
scatterplot3d(x,y,z4,scale.y=0.6,angle=10,main=0.6)
scatterplot3d(x,y,z5,scale.y=0.6,angle=10,main=0.8)
scatterplot3d(x,y,z6,scale.y=0.6,angle=10,main=1)

library(rgl)
plot3d(x,y,z3)
plot3d(x,y,z4)
