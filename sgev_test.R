rm(list=ls())
setwd("~/GitHub/gev")
source("./lib/pack.R")
source("./lib/sgevlibrary.R")

# load data
load("./Rexport/RData_setting3_3/numerical_setting.RDa")
file_count = length(lambdaset)
for (i in 1:file_count){
  eval(parse(text = paste0('load("./Rexport/RData_setting3_3/result_train', i, '.RDa")')))
}

# combine result by list
eval(parse(text = paste0('result_list = list(', paste0("result",1:file_count,collapse=",") , ')')))

# for Z %*% beta
p <- ncol(zlist[[1]])
xbss = eval.basis(df_mu$x1,x_bsobj)
ybss = eval.basis(df_mu$x2,y_bsobj)
Z = do.call('cbind', lapply(1:ncol(xbss), function(i) xbss[, i] * ybss)) 
est_mus = lapply(1:length(result_list), function(i) result_list[[i]][1] + Z%*% tail(result_list[[i]],p))

# scale and shape;
a = lapply(1:length(result_list) ,function(i) matrix(result_list[[i]][2:(1+ns*2)], ns, 2, byrow = T) )
plot(a[[10]][,1], par_scale)
plot(a[[10]][,2], par_shape)

# sum of squares of error : mu_s
mu_sse = lapply(1:length(result_list), function(i) norm(est_mus[[i]]-df_mu$par_mu, "2") )
plot(unlist(mu_sse),type="l")
lambdaset[which.min(unlist(mu_sse))]

# plot 
# plot3d(x=df_mu$x1,y=df_mu$x2,z=df_mu$par_mu)        # 3d plot for true parameter
# plot3d(x=df_mu$x1,y=df_mu$x2,z=est_mus[[4]])       # 3d plot for mus

plot(df_mu$par_mu)
points(est_mus[[1]], col = 'red')
plot(df_mu$par_mu)
points(est_mus[[12]], col = 'red')

# plot(df_mu$par_mu, ylim = c(90,150), main="red,0")
# points(est_mus[[1]], col = 'red')  
# plot(par_scale,  main="red,0" , ylim=c(30,50))
# points(result_list[[1]][3*(0:(ns-1))+2], col = 'red')  
# plot(par_shape,  main="red,0", ylim=c(0,0.4))
# points(result_list[[1]][3*(0:(ns-1))+3], col = 'red')  


# train likelihood
like_vec = c()
for (l in 1:length(result_list)){
  like = 0
  for (i in 1:ns)
  {
    x = xlist[[i]]
    m = est_mus[[l]][i]
    s = a[[l]][i,1]
    k = a[[l]][i,2]
    like = like - sum(dgev(x, m, s, k, log = T)) #  loss + (-log likelihood)
  }
  like_vec[l]=like
}
# lambdaset[which.min(like_vec)]
# plot(like_vec,main="negative loglikelihood")

# for AIC, BIC
AICvec = c()
for (l in 1:length(result_list)){
  # Hatmat= Z%*%solve(t(Z)%*%Z +lambdaset[l]*Om)%*%t(Z)
  Hatmat= Z%*%solve(t(Z)%*%Z +lambdaset[l]*Om+diag(1e-08,nrow(Om)))%*%t(Z)
  DF = sum(diag(Hatmat))+2*ns
  AICvec[l] = 2*like_vec[l] + 2*DF 
}
BICvec = c()  
for (l in 1:length(result_list)){
  # Hatmat= Z%*%solve(t(Z)%*%Z +lambdaset[l]*Om)%*%t(Z)
  Hatmat= Z%*%solve(t(Z)%*%Z +lambdaset[l]*Om+diag(1e-08,nrow(Om)))%*%t(Z)
  DF = sum(diag(Hatmat))+2*ns
  BICvec[l] = 2*like_vec[l] + log(ns*nobs)*DF 
}
plot(AICvec,type="l")
plot(BICvec,type="l")
lambdaset[which.min(AICvec)]
lambdaset[which.min(BICvec)]


# rgev for test
testlist = list()
set.seed(4)
for (s in 1:ns){
  testlist[[s]] = rgev(10000,loc=par_mu[s],scale=par_scale[s],shape=par_shape[s])
}

loss=c()
for (l in 1:length(result_list)){
  like = 0
  for (i in 1:ns)
  {
    x = testlist[[i]]
    m = est_mus[[l]][i]
    s = a[[l]][i,1]
    k = a[[l]][i,2]
    like = like - sum(dgev(x, m, s, k, log = T)) #  loss + (-log likelihood)
  }
  loss[l]=like
}
plot(loss,type="l")
lambdaset[which.min(loss)]

