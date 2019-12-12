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

S_num = 8

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

start_time = Sys.time()
rmse_list = list()
hdist_list = list()
for (i in 1:length(result)){
  set.seed(i)
  xlist = list()
  for (s in 1:ns){
    xlist[[s]] = rgev(n,loc=loc[s], scale=sc[s], shape=sh[s])
  }
  rmse_vec = c()
  hdist_vec = c()
  for (j in 1:length(match.vec)){
    l = match.vec[j]
    tvec=result[[i]][[l]]
    loc.vec.reg = tvec[3+(1:p)]
    sc.vec.reg = tvec[(3+p)+(1:p)]
    sh.vec.reg = tvec[(3+2*p)+(1:p)]
    loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
    sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
    sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
    
    h_dist = 0 
    z_hat_vec = c()
    for (s in 1:ns){
      x = GEV(loc=loc[s],scale=sc[s],shape=sh[s]) # RobExtremes
      y = GEV(loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s])
      h_dist = h_dist + HellingerDist(x,y,smooth) # distrEx
      
      Fn = ecdf(xlist[[s]])
      z_hat = qgev(p=Fn(xlist[[s]])*n/(n+1), loc=loc.vec[s], sc=sc.vec[s], sh=sh.vec[s])    
      z_hat_vec[s] = (base::norm(z_hat - xlist[[s]],"2"))^2
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
  set.seed(i)
  l = min.ind[i]
  xlist = list()
  for (s in 1:ns){
    xlist[[s]] = rgev(n,loc=loc[s], scale=sc[s], shape=sh[s])
  }
  tvec=result[[i]][[l]]
  loc.vec.reg = tvec[3+(1:p)]
  sc.vec.reg = tvec[(3+p)+(1:p)]
  sh.vec.reg = tvec[(3+2*p)+(1:p)]
  loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
  sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
  sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
    
  h_dist = 0 
  z_hat_vec = c()
  for (s in 1:ns){
    x = GEV(loc=loc[s],scale=sc[s],shape=sh[s]) # RobExtremes
    y = GEV(loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s])
    h_dist = h_dist + HellingerDist(x,y,smooth) # distrEx
    
    Fn = ecdf(xlist[[s]])
    z_hat = qgev(p=Fn(xlist[[s]])*n/(n+1), loc=loc.vec[s], sc=sc.vec[s], sh=sh.vec[s])    
    z_hat_vec[s] = (base::norm(z_hat - xlist[[s]],"2"))^2
  }
  hdist_vec[i] = h_dist
  rmse_vec[i] = sqrt(sum(z_hat_vec))/sqrt(n*ns)

  cat(i," ")
}
end_time = Sys.time()

rmse_mat = cbind(do.call('rbind',rmse_list),rmse_vec)
hdist_mat = cbind(do.call('rbind',hdist_list),hdist_vec)

end_time -start_time
eval(parse(text = paste0("save.image(file =","'", paste0('./Rexport/RData_sgev3_simulation/trainloss_',S_num, '.RData',"')"))))
# eval(parse(text = paste0("save.image(file =","'", paste0('./trainloss_',S_num, '.RData',"')"))))
boxplot(rmse_mat, outline=F)
boxplot(hdist_mat, outline=F)
