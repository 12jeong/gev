rm(list=ls())
setwd("~/GitHub/gev")
source("./lib/pack.R")
source("./lib/sgevlibrary.R")
# data generation
set.seed(21)
nobs = 100
ns= 50
nBS = 5
# x1 = seq(-5,5,length=100)
# x2 = seq(-5,5,length=100)
# fx = outer(x1, x2, function(x1,x2) {
#   zval = 100+(2*x1-3*x2)/4
# })
# # plot_ly(x = x1, y = x2, z = fx) %>% add_surface()

x1 = runif(ns,-5,5)
x2 = runif(ns,-5,5)
zval = 100+(2*x1-3*x2)/4 # 평면
plot3d(x1,x2,zval)

# matrix(1:n^2,nrow=n)

df_mu = data.frame(x1=x1,x2=x2,z=zval)
plot3d(x=df_mu$x1,y=df_mu$x2,z=df_mu$z)

xlist = list()
for (i in 1:ns){
  xlist[[i]] = rgev(nobs,loc=df_mu$z[i],scale=40,shape=0.1)
}


# train
# zlist - 2-D basis matrix by each location
# (fda package)
x_bsobj = create.bspline.basis(c(-5,5),norder=4, 
                               breaks=quantile(x1,prob = seq(0, 1, length = nBS)))
y_bsobj = create.bspline.basis(c(-5,5),norder=4, 
                               breaks=quantile(x2,prob = seq(0, 1, length = nBS)))
zlist = list()
for (i in 1:ns){
  xbs = eval.basis(df_mu$x1[i],x_bsobj)
  ybs = eval.basis(df_mu$x2[i],y_bsobj)
  tensorbs = do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs)) 
  zlist[[i]] = tensorbs 
}
sum(zlist[[3]])
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
#                    lambda = 0.3, lambda2=1, Om=Om, mat=mat, method="B-spline")


# save(result1,result2,result3,file="result0724.RDa")
load("result0724.RDa")

p <- ncol(zlist[[1]])
xbss <- eval.basis(x1,x_bsobj)
ybss <- eval.basis(x2,y_bsobj)
tensorbss <- do.call('cbind', lapply(1:ncol(xbss), function(i) xbss[, i] * ybss)) 
z1 <- result1[1]+tensorbss %*% tail(result1,p)
z2 <- result2[1]+tensorbss %*% tail(result2,p)
z3 <- result3[1]+tensorbss %*% tail(result3,p)

#
par(mfrow = c(3,1))
plot(df_mu$z, ylim = c(90,150), main="red,0")
points(z2, col = 'red')       # 0
norm(z2 - df_mu$z,"2")

plot(df_mu$z, ylim = c(90,150), main="red,0.3")
points(z3, col = 'red')       # 0.3
norm(z3 - df_mu$z,"2")

plot(df_mu$z, ylim = c(90,150),main="red,1")
points(z1, col = 'red')       # 1
norm(z1 - df_mu$z,"2")


plot3d(df_mu$x1,df_mu$x2,df_mu$z)
plot3d(df_mu$x1,df_mu$x2,z2,main="Est,0")
plot3d(df_mu$x1,df_mu$x2,z1,main="Est,1")
plot3d(df_mu$x1,df_mu$x2,z3,main="Est,0.3")

# 서버에서 lambda1 0,...,2로 돌려보기


par(mfrow=c(2,2))
scatterplot3d(df_mu$x1,df_mu$x2,df_mu$z,scale.y=0.4,angle=20,main="True")
scatterplot3d(df_mu$x1,df_mu$x2,z2,scale.y=0.4,angle=20,main="Est,0")
scatterplot3d(df_mu$x1,df_mu$x2,z3,scale.y=0.4,angle=20,main="Est,0.3")
scatterplot3d(df_mu$x1,df_mu$x2,z1,scale.y=0.4,angle=20,main="Est,1")
