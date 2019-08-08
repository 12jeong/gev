rm(list=ls())
setwd("~/GitHub/gev")
source("./lib/pack.R")
source("./lib/sgevlibrary.R")

# data generation
set.seed(1)
nobs = 100
ns= 50
######################## mu setting ##########################
### setting1 : 평면 ###
xyrange = c(-10,10)
nBS = 3
x1 = seq(xyrange[1],xyrange[2],length=100)
x2 = seq(xyrange[1],xyrange[2],length=100)
fx = outer(x1, x2, function(x1,x2) {
  par_mu = 100+(-2*x1-3*x2)
})
plot_ly(x = x1, y = x2, z = fx) %>% add_surface()
range(fx)
x1 = runif(ns,xyrange[1],xyrange[2])
x2 = runif(ns,xyrange[1],xyrange[2])
par_mu = 100+(-2*x1-3*x2) # 평면
range(par_mu)
df_mu = data.frame(x1=x1,x2=x2,par_mu=par_mu)
plot3d(x=df_mu$x1,y=df_mu$x2,z=par_mu)

### setting2 : 일봉분포 ###
# xyrange = c(-10,10)
# nBS = 3
# n = 30
# x1 = seq(xyrange[1],xyrange[2],length.out=n)
# x2 = seq(xyrange[1],xyrange[2],length.out=n)
# mu = c(0,0)
# sig = matrix(c(30,5,5,30),nrow=2)
# fx = outer(x1, x2, function(x1,x2) {
#   dmvnorm(cbind(x1,x2), mean=mu, sigma=sig)
# })
# plot_ly(x = x1, y = x2, z = 90+fx*4000) %>% add_surface()
# 
# set.seed(3)
# s_ind = sort(sample(n^2,ns))
# zval = rep(0, ns)
# k = 1
# for ( i in s_ind)
# {
#   zval[k] =fx[i]
#   k = k + 1
# }
# s_col = s_ind%/%n + 1 ; s_col[s_ind%/%n == n] = n
# s_row = s_ind%%n + 1 ; s_row[s_row == 0] = n
# # par_mu = 85+zval*6000 # setting2_1 : 86~115
# par_mu = 90+zval*4000 # setting2_2 : 91~110
# df_mu = data.frame(x1=x1[s_col],x2=x2[s_row],par_mu=par_mu)
# range(par_mu)
# plot3d(x=df_mu$x1,y=df_mu$x2,z=par_mu)

### setting3 : 이봉분포 ###
# set.seed(1)
# nBS = 5
# xyrange = c(-10,10)
# n=100
# x1 = seq(xyrange[1],xyrange[2],length=n)
# x2 = seq(xyrange[1],xyrange[2],length=n)
# mu1 = c(5,0)
# mu2 = c(-5,0)
# sig = matrix(c(10,5,5,10),nrow=2)
# fx = outer(x1, x2, function(x1,x2) {
#   0.4*dmvnorm(cbind(x1,x2),mean=mu1, sigma=sig*1) +
#   0.6*dmvnorm(cbind(x1,x2),mean=mu2, sigma=sig*2)
# })
# plot_ly(x = x1, y = x2, z = 95+2000*fx) %>% add_surface()
# 
# set.seed(3)
# s_ind = sort(sample(n^2,ns))
# zval = rep(0, ns)
# k = 1
# for ( i in s_ind)
# {
#   zval[k] =fx[i]
#   k = k + 1
# }
# s_col = s_ind%/%n + 1 ; s_col[s_ind%/%n == n] = n
# s_row = s_ind%%n + 1 ; s_row[s_row == 0] = n
# # par_mu = 85+zval*5000 # setting3_1 : 85~113, nBS=5
# # par_mu = 90+zval*4000 # setting3_2 : 90~112, nBS=5
# par_mu = 95+zval*2000 # setting3_3 : 95~107, nBS=5
# range(par_mu)
# df_mu = data.frame(x1=x1[s_col],x2=x2[s_row],par_mu=par_mu)
# plot3d(x=df_mu$x1,y=df_mu$x2,z=par_mu)

######################################################################
set.seed(2)
par_scale = rtruncnorm(ns, a=35, b=45, mean=40, sd=2)
par_shape = runif(ns,0.1,0.25)

xlist = list()
for (i in 1:ns){
  xlist[[i]] = rgev(nobs,loc=par_mu[i],scale=par_scale[i],shape=par_shape[i])
}
sum(unlist(lapply(1:length(xlist),function(x) sum(is.na(xlist[[x]])))))

# zlist - 2-D basis matrix by each location
# (fda package)
x_bsobj = create.bspline.basis(xyrange,norder=4, 
                               breaks=quantile(x1,prob = seq(0, 1, length = nBS)))
y_bsobj = create.bspline.basis(xyrange,norder=4, 
                               breaks=quantile(x2,prob = seq(0, 1, length = nBS)))

zlist = list()
for (i in 1:ns){
  xbs = eval.basis(df_mu$x1[i],x_bsobj)
  ybs = eval.basis(df_mu$x2[i],y_bsobj)
  tensorbs = do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs)) 
  zlist[[i]] = tensorbs 
}
xbss = eval.basis(df_mu$x1,x_bsobj)
ybss = eval.basis(df_mu$x2,y_bsobj)
Z = do.call('cbind', lapply(1:ncol(xbss), function(i) xbss[, i] * ybss)) 

# sum(zlist[[3]])
# dim(zlist[[1]]) # (frist)stnlds 2D-splines tensor, nbasis = df x df

# Omega matrix
Fmat = kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat = kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat = kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om = Fmat+Gmat+2*Hmat

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

# to save setting values
# save(list=ls(),file="./numerical_setting.RDa")

# to etimate model
# lambdaset = c(seq(0,2,length=21),5,10,100)
lambda = 10

fit = gevreg_m(xlist,zlist,
               lambda = lambda, Om= Om, 
               method="B-spline")
result = fit$par
# scale and shape;
a = matrix(result[2:(1+ns*2)], ns, 2, byrow = T) 
plot(a[,1], par_scale)
plot(a[,2], par_shape)

p = length(drop(zlist[[1]]))
b = c() # Z %*% beta 
for (i in 1:length(zlist))  b[i] = sum( drop(zlist[[i]])*tail(result,p) )
est_mus = b + result[1]  # ms = m0 + Z%*%beta, length(est_mus)=ns
plot(est_mus, par_mu)

# BIC and AIC
like = 0
for (i in 1:ns)
{
  x = xlist[[i]]
  m = est_mus[i]
  s = a[i,1]
  k = a[i,2]
  like = like - sum(dgev(x, m, s, k, log = T)) #  loss + (-log likelihood)
}
like
Hatmat= Z%*%solve(t(Z)%*%Z +lambda*Om+diag(1e-08,nrow(Om)))%*%t(Z) # 가능한가?
DF = sum(diag(Hatmat))+2*ns

AIC = 2*like + 2*DF ; AIC
BIC = 2*like + log(ns*nobs)*DF ; BIC

# save.file for each lambda
# eval(parse(text = paste0('result', lambda_idx, ' = result')))
# 
# eval(parse(text = paste0('save(result', lambda_idx,
#                          paste0(",file =","'", paste0('./result_train',lambda_idx, '.RDa',"')")))))


plot3d(df_mu$x1,df_mu$x2,est_mus)

plot(df_mu$par_mu)
points(est_mus, col = 'red')
