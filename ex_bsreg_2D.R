rm(list=ls())
library(fda)
library(evd)
library(Deriv)
source("sgevlibrary.R")
setwd("C:\\Users\\UOS\\Documents\\GITHUB\\gev")
load("kma_data\\Pr_46.RData")

ss <- split.data.frame(Pr_46,Pr_46$stnlds)
xlist <- lapply(ss,"[[","pr")
ns <- length(unique(Pr_46$stnlds))

# create.bspline.basis : nbasis = norder + length(breaks) -2
# order=4, quantile knots, nbasis=7, 
x_bsobj <- create.bspline.basis(range(Pr_46$long),breaks=quantile(Pr_46$long,prob = seq(0, 1, length = 5)))
y_bsobj <- create.bspline.basis(range(Pr_46$lat),breaks=quantile(Pr_46$lat,prob = seq(0, 1, length = 5)))

zlist <- list()
for (i in 1:ns){
  xbs <- eval.basis(ss[[i]]$long,x_bsobj)
  ybs <- eval.basis(ss[[i]]$lat,y_bsobj)
  tensorbs <- do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs)) 
  zlist[[i]] <- tensorbs 
}

dim(zlist[[1]]) # (frist)stnlds 2D-splines tensor, nbasis = 7 x 7

# n = nrow(zlist[[1]])
# p = ncol(zlist[[1]])
# true_theta=c(100,30,0.1)
# true_beta = rep(1,p) # setting?
# xlist = list()
# for (i in 1:ns){
#   z = zlist[[i]]
#   eps = rgev(n,loc=true_theta[1], 
#              scale=true_theta[2], 
#              shape=true_theta[3])
#   xlist[[i]] = z%*%true_beta + eps 
# }

# 2-D splines penlaty matrix
Fmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om <- Fmat+2*Gmat+Hmat

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

# result1 <- round(gevreg_m(xlist, zlist, lambda = 0, Om=Om, method="B-spline"),3)
# result2 <- round(gevreg_m(xlist, zlist, lambda = 0.01, Om=Om, method="B-spline"),3)
result3 <- round(gevreg_m(xlist, zlist, lambda = 1, lambda2=1, Om=Om, mat=mat, method="B-spline"),3)
result4 <- round(gevreg_m(xlist, zlist, lambda = 0.1, lambda2=1, Om=Om, mat=mat, method="B-spline"),3)

# result1[1:(ns*3)]
# result2[1:(ns*3)]
# result3[1:(ns*3)]

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
scatterplot3d(x,y,z4,scale.y=0.6,angle=10)

plot3d(x,y,z4)
# bsmat <- data.frame(x,y,z1,z2,z3)

# x,y,z plot 그리기 
# library(plotly)
# library(rgl)
# plot3d(x,y,z1,col="#999999")
# points3d(x,y,z2,col="#E69F00")
# points3d(x,y,z3,col="#56B4E9")

library(scatterplot3d)
par(mfrow=c(1,3))
scatterplot3d(x,y,z1,scale.y=0.6,angle=10)
scatterplot3d(x,y,z2,scale.y=0.6,angle=10)
scatterplot3d(x,y,z3,scale.y=0.6,angle=10)


save(result1,result2,result3,result4,file="ex_bsreg_2D.RData")
# load("ex_bsreg_2D.RData")
