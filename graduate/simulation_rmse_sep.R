# mu, sigma, kappa에 대해 따로 rmse 보는 코드

rm(list=ls())
setwd("~/GITHUB/gev")
# setwd("C:/Users/UOS/Downloads")
source("./lib/sgev3library.R")
source("./lib/pack.R")
# source("C:/Users/UOS/Documents/GITHUB/gev/lib/sgev3library.R")
# source("C:/Users/UOS/Documents/GITHUB/gev/lib/pack.R")

library(distrEx)
library(RobExtremes) 

S_num = 1

eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/AIC_scenario',S_num, '.RData',"')"))))
# eval(parse(text = paste0("load(file =","'", paste0('./AIC_scenario',S_num, '.RData',"')"))))
# lam.grid2 = expand.grid(lam_set1,lam_set2,lam_set3)

# true paramter
loc = mu_set[,S_map[S_num,]$Var1]
sc =  sc_set[,S_map[S_num,]$Var2]
sh =  sh_set[,S_map[S_num,]$Var3]

rmse_list = list()
# hdist_list = list()
for (i in 1:length(result)){
  rmse_vec_mu = c()
  rmse_vec_sc = c()
  rmse_vec_sh = c()
  # hdist_vec = c()
  for (l in 1:length(result[[1]])){
    tvec=result[[i]][[l]]
    loc.vec.reg = tvec[3+(1:p)]
    sc.vec.reg = tvec[(3+p)+(1:p)]
    sh.vec.reg = tvec[(3+2*p)+(1:p)]
    loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
    sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
    sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
    # h_dist = 0 
    # for (s in 1:ns){
    #   x = GEV(loc=loc[s],scale=sc[s],shape=sh[s]) # RobExtremes
    #   y = GEV(loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s])
    #   h_dist = h_dist + HellingerDist(x,y,smooth) # distrEx
    # }
    rmse_vec_mu[l] = base::norm(loc-loc.vec, "2")/sqrt(ns)
    rmse_vec_sc[l] = base::norm(sc-sc.vec, "2")/sqrt(ns)
    rmse_vec_sh[l] = base::norm(sh-sh.vec, "2")/sqrt(ns)
    # hdist_vec[l] = h_dist
  }
  rmse_list[[i]] = list(rmse_vec_mu,rmse_vec_sc,rmse_vec_sh)
  # hdist_list[[i]] = hdist_vec
}

rmse_mu_mat = do.call("rbind",lapply(rmse_list,function(x) x[[1]]))
rmse_sc_mat = do.call("rbind",lapply(rmse_list,function(x) x[[2]]))
rmse_sh_mat = do.call("rbind",lapply(rmse_list,function(x) x[[3]]))

min.ind =  unlist(lapply(AIC_list,function(x) which.min(x)))
lam.min = as.numeric(names(which.max(table(min.ind))))
# lam.grid2[as.numeric(names(table(min.ind))),]
lam.min.vec = lam.grid2[lam.min,]

a1=c(); a2=c(); a3=c()
for(i in 1:nrow(rmse_mu_mat)) {
  a1[i] = rmse_mu_mat[i,min.ind[i]]
  a2[i] = rmse_sc_mat[i,min.ind[i]]
  a3[i] = rmse_sh_mat[i,min.ind[i]]
}

loc.map = c(0,max(lam_set1),lam.min.vec[1])
sc.map = c(0,max(lam_set2),lam.min.vec[2])
sh.map = c(0,max(lam_set3),lam.min.vec[3])

loc.df = data.frame(c(rmse_mu_mat[,which(lam.grid2[,1]== loc.map[1])]),
                    c(rmse_mu_mat[,which(lam.grid2[,1]== loc.map[2])]),
                    c(rmse_mu_mat[,which(lam.grid2[,1]== loc.map[3])]),
                    a1)
sc.df = data.frame(c(rmse_sc_mat[,which(lam.grid2[,2]== sc.map[1])]),
                   c(rmse_sc_mat[,which(lam.grid2[,2]== sc.map[2])]),
                   c(rmse_sc_mat[,which(lam.grid2[,2]== sc.map[3])]),
                   a2)
sh.df = data.frame(c(rmse_sh_mat[,which(lam.grid2[,3]== sh.map[1])]),
                   c(rmse_sh_mat[,which(lam.grid2[,3]== sh.map[2])]),
                   c(rmse_sh_mat[,which(lam.grid2[,3]== sh.map[3])]),
                   a3)
par(mfrow=c(2,3))
boxplot(loc.df); boxplot(sc.df); boxplot(sh.df)

################# for test ########################
# sampling X coordinates 
ns = 10
xyrange = c(-10,10)
set.seed(201)
x1.test = runif(ns,xyrange[1],xyrange[2])
x2.test = runif(ns,xyrange[1],xyrange[2])

# surface base setting
mean_vec = c(0,0) 
sig_mat = matrix(c(30,0,0,30),nrow=2)
# set.seed(202)
set_uni = dmvnorm(cbind(x1.test,x2.test), mean=mean_vec, sigma=sig_mat)
mean_vec1 = c(5,0); mean_vec2 = c(-5,0) 
sig_mat = matrix(c(10,0,0,10),nrow=2)
# set.seed(203)
set_bi = 0.4*dmvnorm(cbind(x1.test,x2.test),mean=mean_vec1, sigma=sig_mat*1) +0.6*dmvnorm(cbind(x1.test,x2.test),mean=mean_vec2, sigma=sig_mat*2)

# location setting
mu1.test = 100 + (-3*x1.test + 3*x2.test) 
mu2.test = 90 + set_uni*4000
mu3.test = 90 + set_bi*3000 
mu_set.test = data.frame(plane=mu1.test,unimodal=mu2.test,bimodal=mu3.test)

# scale setting
sc1.test = 40 + (-3*x1.test + 3*x2.test)*0.5
sc2.test = 30 + set_uni*4000 
sc3.test = 30 + set_bi*3000 
sc_set.test = data.frame(plane=sc1.test,unimodal=sc2.test,bimodal=sc3.test)

# shape setting
sh1.test = 0.1+(-3*x1.test + 3*x2.test)*0.005
sh2.test = set_uni*50 
sh3.test = set_bi*50 
sh_set.test = data.frame(plane=sh1.test,unimodal=sh2.test,bimodal=sh3.test)

# generate test data
loc.test = mu_set.test[,S_map[S_num,]$Var1]
sc.test =  sc_set.test[,S_map[S_num,]$Var2]
sh.test =  sh_set.test[,S_map[S_num,]$Var3]

# TPS basis for test data
nBS=3
x_bsobj = create.bspline.basis(xyrange,norder=4, 
                               breaks=quantile(x1,prob = seq(0, 1, length = nBS)))
y_bsobj = create.bspline.basis(xyrange,norder=4, 
                               breaks=quantile(x2,prob = seq(0, 1, length = nBS)))
zlist.test = list()
for (i in 1:ns){
  xbs = eval.basis(x1.test[i],x_bsobj)
  ybs = eval.basis(x2.test[i],y_bsobj)
  tensorbs = do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs)) 
  zlist.test[[i]] = tensorbs 
}
Z.test = do.call('rbind', zlist.test) 
p = ncol(Z.test)


rmse_list = list()
# hdist_list = list()
for (i in 1:length(result)){
  rmse_vec_mu = c()
  rmse_vec_sc = c()
  rmse_vec_sh = c()
  # hdist_vec = c()
  for (l in 1:length(result[[1]])){
    tvec=result[[i]][[l]]
    loc.vec.reg = tvec[3+(1:p)]
    sc.vec.reg = tvec[(3+p)+(1:p)]
    sh.vec.reg = tvec[(3+2*p)+(1:p)]
    loc.vec = tvec[1] + drop(Z.test%*%loc.vec.reg)
    sc.vec = exp(tvec[2] + drop(Z.test%*%sc.vec.reg))
    sh.vec = tvec[3] + drop(Z.test%*%sh.vec.reg)
    # h_dist = 0 
    # for (s in 1:ns){
    #   x = GEV(loc=loc[s],scale=sc[s],shape=sh[s]) # RobExtremes
    #   y = GEV(loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s])
    #   h_dist = h_dist + HellingerDist(x,y,smooth) # distrEx
    # }
    rmse_vec_mu[l] = base::norm(loc.test-loc.vec, "2")/sqrt(ns)
    rmse_vec_sc[l] = base::norm(sc.test-sc.vec, "2")/sqrt(ns)
    rmse_vec_sh[l] = base::norm(sh.test-sh.vec, "2")/sqrt(ns)
    # hdist_vec[l] = h_dist
  }
  rmse_list[[i]] = list(rmse_vec_mu,rmse_vec_sc,rmse_vec_sh)
  # hdist_list[[i]] = hdist_vec
}

rmse_mu_mat = do.call("rbind",lapply(rmse_list,function(x) x[[1]]))
rmse_sc_mat = do.call("rbind",lapply(rmse_list,function(x) x[[2]]))
rmse_sh_mat = do.call("rbind",lapply(rmse_list,function(x) x[[3]]))

min.ind =  unlist(lapply(AIC_list,function(x) which.min(x)))
lam.min = as.numeric(names(which.max(table(min.ind))))
# lam.grid2[as.numeric(names(table(min.ind))),]
lam.min.vec = lam.grid2[lam.min,]
a1=c(); a2=c(); a3=c()
for(i in 1:nrow(rmse_mu_mat)) {
  a1[i] = rmse_mu_mat[i,min.ind[i]]
  a2[i] = rmse_sc_mat[i,min.ind[i]]
  a3[i] = rmse_sh_mat[i,min.ind[i]]
}

loc.map = c(0,max(lam_set1),lam.min.vec[1])
sc.map = c(0,max(lam_set2),lam.min.vec[2])
sh.map = c(0,max(lam_set3),lam.min.vec[3])

loc.df = data.frame(c(rmse_mu_mat[,which(lam.grid2[,1]== loc.map[1])]),
                    c(rmse_mu_mat[,which(lam.grid2[,1]== loc.map[2])]),
                    c(rmse_mu_mat[,which(lam.grid2[,1]== loc.map[3])]),
                    a1)
sc.df = data.frame(c(rmse_sc_mat[,which(lam.grid2[,2]== sc.map[1])]),
                   c(rmse_sc_mat[,which(lam.grid2[,2]== sc.map[2])]),
                   c(rmse_sc_mat[,which(lam.grid2[,2]== sc.map[3])]),
                   a2)
sh.df = data.frame(c(rmse_sh_mat[,which(lam.grid2[,3]== sh.map[1])]),
                   c(rmse_sh_mat[,which(lam.grid2[,3]== sh.map[2])]),
                   c(rmse_sh_mat[,which(lam.grid2[,3]== sh.map[3])]),
                   a3)
# par(mfrow=c(2,3))
boxplot(loc.df); boxplot(sc.df,outline = F); boxplot(sh.df,outline = F)


apply(loc.df,2,mean)
apply(sc.df,2,mean)
apply(sh.df,2,mean)
