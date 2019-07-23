rm(list=ls())

########### 이변량정규분포3차원그래프 ##############
if(!require(mvtnorm)) install.packages('mvtnorm'); require(mvtnorm)
if(!require(plotly)) install.packages('plotly'); require(plotly)
if(!require(evd)) install.packages('evd'); require(evd)

# data generation
set.seed(2019)
n=100
x1 = seq(-5,5,length.out=n)
x2 = seq(-5,5,length.out=n)
mu1 = c(0,-2)
mu2 = c(0,2)
sig = matrix(c(2,1,1,2),nrow=2)
fx = outer(x1, x2, function(x1,x2) { 
  0.4*dmvnorm(cbind(x1,x2),mean=mu1, sigma=sig) +
  0.6*dmvnorm(cbind(x1,x2),mean=mu2, sigma=sig)
  })
plot_ly(x = x1, y = x2, z = fx) %>% add_surface()

# matrix(1:n^2,nrow=n)
ns = 50
set.seed(21)
s_ind = sample(n^2,ns)
s_row = trunc(s_ind/100)+1 # 1, 100, 101, 10000
s_row[s_ind%%100==0] = s_ind[s_ind%%100==0]/100
s_col = (s_ind/100-s_row+1)*100 
zval = c()
for ( i in 1:ns) { zval[i] =fx[s_row[i],s_col[i]] ;cat(i,fx[s_row[i],s_col[i]],"\n") }
df_mu = data.frame(x1=x1[s_col],x2=x2[s_row],z=100+zval*100)
plot3d(x=df_mu$x1,y=df_mu$x2,z=df_mu$z)

xlist = list()
for (i in 1:ns){
xlist[[i]] = rgev(n,loc=df_mu$z[i],scale=40,shape=0.1)
}


# train
setwd("C:\\Users\\UOS\\Documents\\GITHUB\\gev")
source("sgevlibrary.R")
if(!require(fda)) install.packages('fda'); library(fda)

# zlist - 2-D basis matrix by each location
x_bsobj = create.bspline.basis(range(x1),breaks=quantile(x1,prob = seq(0, 1, length = 5)))
y_bsobj = create.bspline.basis(range(x2),breaks=quantile(x2,prob = seq(0, 1, length = 5)))
zlist = list()
for (i in 1:ns){
  xbs = eval.basis(df_mu$x1[i],x_bsobj)
  ybs = eval.basis(df_mu$x2[i],y_bsobj)
  tensorbs = do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs)) 
  zlist[[i]] = tensorbs 
}
# dim(zlist[[1]]) # (frist)stnlds 2D-splines tensor, nbasis = df x df

# Omega matrix
Fmat = kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat = kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat = kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om = Fmat+Gmat+2*Hmat

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

# design matrix for v3 (mu_0 stationary)
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

# result1 = gevreg_m(xlist,zlist,
#                    lambda = 1, lambda2=1, Om=Om, mat=mat, method="B-spline")
# result2 = gevreg_m(xlist,zlist,
#                    lambda = 0, lambda2=1, Om=Om, mat=mat, method="B-spline")
# result3 = gevreg_m(xlist,zlist,
#                    lambda = 0.5, lambda2=1, Om=Om, mat=mat, method="B-spline")

# save(result1,result2,result3,file="result0705.RData")
load("result0705.RData")

p <- ncol(zlist[[1]])
xbss <- eval.basis(x1[s_row],x_bsobj)
ybss <- eval.basis(x2[s_col],y_bsobj)
tensorbss <- do.call('cbind', lapply(1:ncol(xbss), function(i) xbss[, i] * ybss)) 
z1 <- 100+tensorbss %*% tail(result1,p)
z2 <- 100+tensorbss %*% tail(result2,p)
z3 <- 100+tensorbss %*% tail(result3,p)

library(scatterplot3d)
par(mfrow=c(2,2))
scatterplot3d(df_mu$x1,df_mu$x2,df_mu$z,scale.y=1,angle=40,main="True")
scatterplot3d(df_mu$x1,df_mu$x2,z2,scale.y=1,angle=40,main="Est,0")
scatterplot3d(df_mu$x1,df_mu$x2,z3,scale.y=1,angle=40,main="Est,0.5")
scatterplot3d(df_mu$x1,df_mu$x2,z1,scale.y=1,angle=40,main="Est,1")
scatterplot3d(df_mu$x2,df_mu$x1,df_mu$z,scale.y=1,angle=40,main="True")
scatterplot3d(df_mu$x2,df_mu$x1,z2,scale.y=1,angle=40,main="Est,0")
scatterplot3d(df_mu$x2,df_mu$x1,z3,scale.y=1,angle=40,main="Est,0.5")
scatterplot3d(df_mu$x2,df_mu$x1,z1,scale.y=1,angle=40,main="Est,1")


library(rgl)
plot3d(df_mu$x1,df_mu$x2,df_mu$z)
plot3d(df_mu$x1,df_mu$x2,z2,main="Est,0")
plot3d(df_mu$x1,df_mu$x2,z3,main="Est,0.5")
plot3d(df_mu$x1,df_mu$x2,z1,main="Est,1")
