rm( list = ls()); gc()
setwd("~/GITHUB/gev/Rexport/Spatial_simulation")
library(evd)
library(distrEx)
library(RobExtremes) 

# load data
load("./result3/numerical_setting3.RDa") 
lambdaset

file_count=100
for (i in 1:file_count){
  eval(parse(text = paste0('load("./result3/result3_seed', i, '.RDa")')))
}

# combine result by list
eval(parse(text = paste0('result_list = list(', paste0("result",1:file_count,collapse=",") , ')')))

# loc for Z %*% beta
p = length(drop(zlist[[1]]))
xbss = eval.basis(df_mu$x1,x_bsobj)
ybss = eval.basis(df_mu$x2,y_bsobj)
Z = do.call('cbind', lapply(1:ncol(xbss), function(i) xbss[, i] * ybss)) 

# sum of squares of error : mu_s
est_mus = lapply(1:file_count, function(i) lapply(1:length(lambdaset),
                                                  function(j) result_list[[i]][[j]][1] + drop(Z%*% tail(result_list[[i]][[j]],p))))

# scale and shape;
a = lapply(1:length(result_list) ,function(i) lapply(1:length(lambdaset),
                                                     function(j) matrix(result_list[[i]][[j]][2:(1+ns*2)], ns, 2, byrow = T) ))
# plot(a[[10]][,1], par_scale)
# plot(a[[10]][,2], par_shape)


SSElist_mu = list()
for ( i in 1:file_count){
  SSEvec_mu=c(); 
  for (j in 1:length(lambdaset)) {
    SSEvec_mu[j]= base::norm(est_mus[[i]][[j]] - df_mu$par_mu,"2")
  }
  SSElist_mu[[i]]= SSEvec_mu
}
SSEmat_mu = do.call('rbind',lapply(1:length(lambdaset), 
                                   function(j) unlist(lapply(1:file_count,function(i) SSElist_mu[[i]][j]))))


SSElist_all = list()
for ( i in 1:file_count){
  SSEvec_all=c(); 
  for (j in 1:length(lambdaset)) {
    SSEvec_all[j]= base::norm(c(est_mus[[i]][[j]],a[[i]][[j]][,1],a[[i]][[j]][,2]) - 
                                c(df_mu$par_mu,par_scale,par_shape),"2")
  }
  SSElist_all[[i]]= SSEvec_all
}
SSEmat_all = do.call('rbind',lapply(1:length(lambdaset), 
                                    function(j) unlist(lapply(1:file_count,function(i) SSElist_all[[i]][j]))))


# train likelihood
like_list = list()
for (l in 1:length(result_list)){
  like_vec = c()
  set.seed(l)
  xlist = list()
  for (i in 1:ns){
    xlist[[i]] = rgev(nobs,loc=par_mu[i],scale=par_scale[i],shape=par_shape[i])
  }
  tt = sum(unlist(lapply(1:length(xlist),function(x) sum(is.na(xlist[[x]])))))==0
  # cat("seed",l,"/",tt,"\n")
  for (j in 1:length(lambdaset)){
    like = 0
    for (i in 1:ns)
    {
      x = xlist[[i]]
      m = est_mus[[l]][[j]][i]
      s = a[[l]][[j]][i,1]
      k = a[[l]][[j]][i,2]
      like = like - sum(dgev(x, m, s, k, log = T)) #  loss + (-log likelihood)
    }
    like_vec[j]=like
  }
  like_list[[l]] = like_vec
}

# for AIC, BIC
AIC_list = list() ; BIC_list = list() 
for (l in 1:length(result_list)){
  AICvec = c(); BICvec = c()
  for (j in 1:length(lambdaset)){
    Hatmat= Z%*%solve(t(Z)%*%Z +lambdaset[j]*Om+diag(1e-08,nrow(Om)))%*%t(Z)
    DF = sum(diag(Hatmat))+2*ns
    AICvec[j] = 2*like_list[[l]][j] + 2*DF 
    BICvec[j] = 2*like_list[[l]][j] + log(ns*nobs)*DF
  }
  AIC_list[[l]] = AICvec
  BIC_list[[l]] = BICvec
}
AICmat = do.call('rbind',lapply(1:length(lambdaset), function(j) unlist(lapply(1:file_count,function(i) AIC_list[[i]][j]))))
BICmat = do.call('rbind',lapply(1:length(lambdaset), function(j) unlist(lapply(1:file_count,function(i) BIC_list[[i]][j]))))

h_dist_list = list()
for (i in 1:file_count){
  h_dist_vec = c()
  for (l in 1:length(lambdaset)){
    h_dist = 0 
    for (s in 1:ns){
      x = GEV(loc=par_mu[s],scale=par_scale[s],shape=par_shape[s]) # RobExtremes
      y = GEV(loc=est_mus[[i]][[l]][s],scale=a[[i]][[l]][s,1],shape=a[[i]][[l]][s,2])
      h_dist = h_dist + HellingerDist(x,y,smooth) # distrEx
    }
    h_dist_vec[l] = h_dist
  }
  h_dist_list[[i]]=h_dist_vec
}
distmat = do.call('rbind',lapply(1:length(lambdaset), function(j) unlist(lapply(1:file_count,function(i) h_dist_list[[i]][j]))))



save(result_list,lambdaset,Z, est_mus,a,
     SSEmat_mu,SSEmat_all,AICmat,BICmat,distmat,like_list,
     par_scale,par_mu,par_shape,df_mu,p
     ,file="result3_analysis.RData")


# rm( list = ls()); gc()
# setwd("/home/hyj0827")
# source("./lib/sgevlibrary.R")
# load("RData/result3_analysis.RData")

