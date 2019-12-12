# 100개 seed에 대해 216 개의 모형에 대한 loss를 계산하는 것이 오래걸린다.
# 필요한 set에 대한 loss만 계산해보자.

rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")

# setwd("C:/Users/UOS/Downloads")
# source("C:/Users/UOS/Documents/GITHUB/gev/lib/sgev3library.R")
# source("C:/Users/UOS/Documents/GITHUB/gev/lib/pack.R")
library(distrEx)
library(RobExtremes) 

S_num =8

eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/AIC_scenario',S_num, '.RData',"')"))))
# eval(parse(text = paste0("load(file =","'", paste0('./AIC_scenario',S_num, '.RData',"')"))))

min.ind =  unlist(lapply(AIC_list,function(x) which.min(x)))
lam.min.ind = as.numeric(names(which.max(table(min.ind))))
lam.min.star = lam.grid2[lam.min.ind,] # 가장 빈번하게 선택된 lam 조합

loc.map = c(0,max(lam_set1),lam.min.star[1])
sc.map = c(0,max(lam_set2),lam.min.star[2])
sh.map = c(0,max(lam_set3),lam.min.star[3])

match.map = expand.grid(loc.map,sc.map,sh.map)
match.vec = c()
for (i in 1:nrow(match.map)){
  vec_tmp = match.map[i,]
  match.vec[i] = which(colSums(apply(lam.grid2,1,function(x) x== vec_tmp))==3)
}
match.vec # match.map에 해당하는 model number


############################### test set ###############################
# sampling X coordinates 
ns = 10
xyrange = c(-10,10)
set.seed(201)
x1.test = runif(ns,xyrange[1],xyrange[2])
x2.test = runif(ns,xyrange[1],xyrange[2])

# surface base setting
mean_vec = c(0,0) 
sig_mat = matrix(c(30,0,0,30),nrow=2)
set_uni = dmvnorm(cbind(x1.test,x2.test), mean=mean_vec, sigma=sig_mat)
mean_vec1 = c(5,0); mean_vec2 = c(-5,0) 
sig_mat = matrix(c(10,0,0,10),nrow=2)
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

set.seed(202)
xlist.test = list()
for (s in 1:ns){
  xlist.test[[s]] = rgev(100,loc=loc.test[s], scale=sc.test[s], shape=sh.test[s])
}

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

########################################################################

start_time = Sys.time()
rmse_list = list()
hdist_list = list()
for (i in 1:length(result)){
  rmse_vec = c()
  hdist_vec = c()
  for (j in 1:length(match.vec)){
    l = match.vec[j]
    tvec=result[[i]][[l]]
    loc.vec.reg = tvec[3+(1:p)]
    sc.vec.reg = tvec[(3+p)+(1:p)]
    sh.vec.reg = tvec[(3+2*p)+(1:p)]
    loc.vec = tvec[1] + drop(Z.test%*%loc.vec.reg)
    sc.vec = exp(tvec[2] + drop(Z.test%*%sc.vec.reg))
    sh.vec = tvec[3] + drop(Z.test%*%sh.vec.reg)
    
    h_dist = 0 
    z_hat_vec = c()
    for (s in 1:ns){
      x = GEV(loc=loc.test[s],scale=sc.test[s],shape=sh.test[s]) # RobExtremes
      y = GEV(loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s])
      h_dist = h_dist + HellingerDist(x,y,smooth) # distrEx
      
      Fn = ecdf(xlist.test[[s]])
      z_hat = qgev(p=Fn(xlist.test[[s]])*n/(n+1), loc=loc.vec[s], sc=sc.vec[s], sh=sh.vec[s])    
      z_hat_vec[s] = (base::norm(z_hat - xlist.test[[s]],"2"))^2
    }
    hdist_vec[j] = h_dist
    rmse_vec[j] = sqrt(sum(z_hat_vec))/sqrt(n*ns)
  }
  rmse_list[[i]] = rmse_vec
  hdist_list[[i]] = hdist_vec
  cat(i," ")
}

# 각 seed의 lam_AIC에 대한 loss도 계산
min.ind 
rmse_vec = c()
hdist_vec = c()
for (i in 1:length(result)){
  l = min.ind[i]
  tvec=result[[i]][[l]]
  loc.vec.reg = tvec[3+(1:p)]
  sc.vec.reg = tvec[(3+p)+(1:p)]
  sh.vec.reg = tvec[(3+2*p)+(1:p)]
  loc.vec = tvec[1] + drop(Z.test%*%loc.vec.reg)
  sc.vec = exp(tvec[2] + drop(Z.test%*%sc.vec.reg))
  sh.vec = tvec[3] + drop(Z.test%*%sh.vec.reg)
  
  h_dist = 0 
  z_hat_vec = c()
  for (s in 1:ns){
    x = GEV(loc=loc.test[s],scale=sc.test[s],shape=sh.test[s]) # RobExtremes
    y = GEV(loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s])
    h_dist = h_dist + HellingerDist(x,y,smooth) # distrEx
    
    Fn = ecdf(xlist.test[[s]])
    z_hat = qgev(p=Fn(xlist.test[[s]])*n/(n+1), loc=loc.vec[s], sc=sc.vec[s], sh=sh.vec[s])    
    z_hat_vec[s] = (base::norm(z_hat - xlist.test[[s]],"2"))^2
  }
  hdist_vec[i] = h_dist
  rmse_vec[i] = sqrt(sum(z_hat_vec))/sqrt(n*ns)
  
  cat(i," ")
}
end_time = Sys.time()

rmse_mat = cbind(do.call('rbind',rmse_list),rmse_vec)
hdist_mat = cbind(do.call('rbind',hdist_list),hdist_vec)


end_time -start_time
eval(parse(text = paste0("save.image(file =","'", paste0('./Rexport/RData_sgev3_simulation/testloss_',S_num, '.RData',"')"))))
# eval(parse(text = paste0("save.image(file =","'", paste0('./testloss_',S_num, '.RData',"')"))))
boxplot(rmse_mat, outline=F)
boxplot(hdist_mat, outline=F)

par(mfrow=c(2,1))
boxplot(rmse_mat[,c(13:18,22:28)], outline=F)
boxplot(hdist_mat[,c(13:18,22:28)], outline=F)

which.min(apply(rmse_mat,2,mean))
which.min(apply(rmse_mat,2,median))
which.min(apply(hdist_mat,2,mean))
which.min(apply(hdist_mat,2,median))
