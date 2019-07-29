rm(list=ls())
setwd("~/GitHub/gev")
source("./lib/pack.R")
source("./lib/sgevlibrary.R")

# load data
load("./numerical_setting.RDa")
# file_count = length(lambdaset)
file_count=1
for (i in 1:file_count){
  eval(parse(text = paste0('load("./result_train', i, '.RDa")')))
}
# combine result by list
eval(parse(text = paste0('result_list = list(', paste0("result",1,collapse=",") , ')')))

# for Z %*% beta
p <- ncol(zlist[[1]])
xbss <- eval.basis(x1,x_bsobj)
ybss <- eval.basis(x2,y_bsobj)
tensorbss <- do.call('cbind', lapply(1:ncol(xbss), function(i) xbss[, i] * ybss)) 
result_mus = lapply(1:length(result_list), function(i) result_list[[i]][1] + tensorbss%*% tail(result_list[[i]],p))


# sum of squares of error : mu_s
mu_sse = lapply(1:length(result_list), function(i) norm(result_mus[[i]]-df_mu$z, "2") )
# plot(c(seq(0,2,length=21),3,4,5),unlist(mu_sse),type="l")
# lambdaset[which.min(unlist(mu_sse))]

# plot 
plot3d(x=df_mu$x1,y=df_mu$x2,z=df_mu$z)         # 3d plot for true parameter
plot3d(x=df_mu$x1,y=df_mu$x2,z=result_mus[[1]]) # 3d plot for mus

plot(df_mu$z, ylim = c(90,150), main="red,0")
points(result_mus[[1]], col = 'red')  

plot(par_scale,  main="red,0" , ylim=c(30,50))
points(result_list[[1]][3*(0:(ns-1))+2], col = 'red')  

plot(par_shape,  main="red,0", ylim=c(0,0.4))
points(result_list[[1]][3*(0:(ns-1))+3], col = 'red')  

# plot(df_mu$z, ylim = c(90,150), main="red,0.1")
# points(result_mus[[2]], col = 'red')
# plot(df_mu$z, ylim = c(90,150), main="red,0.5")
# points(result_mus[[6]], col = 'red')  
# plot(df_mu$z, ylim = c(90,150), main="red,1")
# points(result_mus[[11]], col = 'red')  
# plot(df_mu$z, ylim = c(90,150), main="red,2")
# points(result_mus[[21]], col = 'red')  
# plot(df_mu$z, ylim = c(90,150), main="red,5")
# points(result_mus[[22]], col = 'red')  
# plot(df_mu$z, ylim = c(90,150), main="red,10")
# points(result_mus[[23]], col = 'red')  
# plot(df_mu$z, ylim = c(90,150), main="red,100")
# points(result_mus[[24]], col = 'red')  

# v=c()
# loss=c()
# for (i in 1:length(result_list)){
#   for (s in 1:ns){
#     est_s = result_list[[i]][(1:3)+(3*s-3)]
#     v[s] =  sum(lossfun(x=xlist[[s]], loc=result_mus[[i]][s], scale=est_s[2], shape=est_s[3])) # s번째 지역 loss
#   }
#   loss[i] = sum(v)
# }
# plot(c(1:24),-loss,type="l")
# lambdaset[which.min(-loss)]

# rgev for test
testlist = list()
for (s in 1:ns){
  testlist[[s]] = rgev(10000,loc=df_mu$z[s],scale=par_scale[s],shape=par_shape[s])
}
# loss
v=c()
loss=c()
for (i in 1:length(result_list)){
  for (s in 1:ns){
    est_s = result_list[[i]][(1:3)+(3*s-3)]
    v[s] =  sum(lossfun(x=testlist[[s]], loc=result_mus[[i]][s], scale=est_s[2], shape=est_s[3])) # s번째 지역 loss
  }
  loss[i] = sum(v)
}
loss
# plot(c(1:24),-loss ,type="l")
lambdaset[which.min(-loss)]