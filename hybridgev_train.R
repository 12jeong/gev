rm(list=ls())
setwd("~/GitHub/gev")
source("./lib/pack.R")
source("./lib/sgevlibrary.R")

# data generation
set.seed(1)
nobs = 100
ns = 50
m0 = 80+ (sin(seq(0,2*pi,length = 100))*20 + seq(1,nobs)) / 3
plot(m0); range(m0)
### setting1 : 평면 ###
# xyrange = c(-10,10)
# nBS = 3
# x1 = seq(xyrange[1],xyrange[2],length=100)
# x2 = seq(xyrange[1],xyrange[2],length=100)
# fx = outer(x1, x2, function(x1,x2) {
#   par_mu = (-2*x1-3*x2)/5
# })
# plot_ly(x = x1, y = x2, z = fx) %>% add_surface()
# range(fx)
# x1 = runif(ns,xyrange[1],xyrange[2])
# x2 = runif(ns,xyrange[1],xyrange[2])
# par_mu = (-2*x1-3*x2)/5 # 평면
# range(par_mu)
# df_mu = data.frame(x1=x1,x2=x2,par_mu=par_mu)
# 
# set.seed(2)
# par_scale = rtruncnorm(ns, a=35, b=45, mean=40, sd=2)
# par_shape = runif(ns,0.1,0.25)
# 
# xlist = list()
# for (i in 1:ns){
#   xlist[[i]] = rgev(nobs,loc= m0+par_mu[i],scale=par_scale[i],shape=par_shape[i])
# }
# sum(unlist(lapply(1:length(xlist),function(x) sum(is.na(xlist[[x]])))))
# 
# ns_sample = sample(1:ns,4)
# par(mfrow=c(2,2))
# plot(m0+par_mu[ns_sample[1]]);plot(m0+par_mu[ns_sample[2]])
# plot(m0+par_mu[ns_sample[3]]);plot(m0+par_mu[ns_sample[4]])
# par(mfrow=c(1,1))
# 
# par(mfrow=c(2,2))
# plot(xlist[[ns_sample[1]]]);plot(xlist[[ns_sample[2]]])
# plot(xlist[[ns_sample[3]]]);plot(xlist[[ns_sample[4]]])
# par(mfrow=c(1,1))


## setting2 : 일봉분포 ###
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
# par_mu = -10 + zval*4000
# range(par_mu)
# df_mu = data.frame(x1=x1[s_col],x2=x2[s_row],par_mu=par_mu)
# 
# set.seed(2)
# par_scale = rtruncnorm(ns, a=35, b=45, mean=40, sd=2)
# par_shape = runif(ns,0.1,0.25)
# 
# xlist = list()
# for (i in 1:ns){
#   xlist[[i]] = rgev(nobs,loc= m0+par_mu[i],scale=par_scale[i],shape=par_shape[i])
# }
# sum(unlist(lapply(1:length(xlist),function(x) sum(is.na(xlist[[x]])))))
# 
# ns_sample = sample(1:ns,4)
# par(mfrow=c(2,2))
# plot(m0+par_mu[ns_sample[1]]);plot(m0+par_mu[ns_sample[2]])
# plot(m0+par_mu[ns_sample[3]]);plot(m0+par_mu[ns_sample[4]])
# par(mfrow=c(1,1))
# 
# par(mfrow=c(2,2))
# plot(xlist[[ns_sample[1]]]);plot(xlist[[ns_sample[2]]])
# plot(xlist[[ns_sample[3]]]);plot(xlist[[ns_sample[4]]])
# par(mfrow=c(1,1))


### setting3 : 이봉분포 ###
set.seed(1)
nBS = 5
xyrange = c(-10,10)
n=100
x1 = seq(xyrange[1],xyrange[2],length=n)
x2 = seq(xyrange[1],xyrange[2],length=n)
mu1 = c(5,0)
mu2 = c(-5,0)
sig = matrix(c(10,5,5,10),nrow=2)
fx = outer(x1, x2, function(x1,x2) {
  0.4*dmvnorm(cbind(x1,x2),mean=mu1, sigma=sig*1) +
    0.6*dmvnorm(cbind(x1,x2),mean=mu2, sigma=sig*2)
})
# matrix(1:n^2,nrow=n)
set.seed(3)
s_ind = sort(sample(n^2,ns))
zval = rep(0, ns)
k = 1
for ( i in s_ind)
{
  zval[k] =fx[i]
  k = k + 1
}
s_col = s_ind%/%n + 1 ; s_col[s_ind%/%n == n] = n
s_row = s_ind%%n + 1 ; s_row[s_row == 0] = n
par_mu = -6+zval*2000
range(par_mu)
df_mu = data.frame(x1=x1[s_col],x2=x2[s_row],par_mu=par_mu)

set.seed(2)
par_scale = rtruncnorm(ns, a=35, b=45, mean=40, sd=2)
par_shape = runif(ns,0.1,0.25)

xlist = list()
for (i in 1:ns){
  xlist[[i]] = rgev(nobs,loc= m0+par_mu[i],scale=par_scale[i],shape=par_shape[i])
}
sum(unlist(lapply(1:length(xlist),function(x) sum(is.na(xlist[[x]])))))

ns_sample = sample(1:ns,4)
par(mfrow=c(2,2))
plot(m0+par_mu[ns_sample[1]]);plot(m0+par_mu[ns_sample[2]])
plot(m0+par_mu[ns_sample[3]]);plot(m0+par_mu[ns_sample[4]])
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(xlist[[ns_sample[1]]]);plot(xlist[[ns_sample[2]]])
plot(xlist[[ns_sample[3]]]);plot(xlist[[ns_sample[4]]])
par(mfrow=c(1,1))

