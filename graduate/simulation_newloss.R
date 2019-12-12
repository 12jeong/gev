rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")
library(distrEx)
library(RobExtremes) 


S_num=27
eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/result_scenario',S_num, '.RData',"')"))))


# scenario Mapping
mapa = c("Plane","Unimodal","Bimodal")
S_map = expand.grid(1:3, 1:3, 1:3)
S_map2 = expand.grid(mu=mapa, sc=mapa, sh=mapa)

# sampling X coordinates 
ns = 10
xyrange = c(-10,10)
set.seed(201)
x1.test = runif(ns,xyrange[1],xyrange[2])
x2.test = runif(ns,xyrange[1],xyrange[2])

# surface base setting
mean_vec = c(0,0) 
sig_mat = matrix(c(30,0,0,30),nrow=2)
set.seed(202)
set_uni = dmvnorm(cbind(x1.test,x2.test), mean=mean_vec, sigma=sig_mat)
mean_vec1 = c(5,0); mean_vec2 = c(-5,0) 
sig_mat = matrix(c(10,0,0,10),nrow=2)
set.seed(203)
set_bi = 0.4*dmvnorm(cbind(x1.test,x2.test),mean=mean_vec1, sigma=sig_mat*1) +0.6*dmvnorm(cbind(x1.test,x2.test),mean=mean_vec2, sigma=sig_mat*2)

# location setting
mu1.test = 100 + (-2*x1.test + 3*x2.test) 
mu2.test = 90 + set_uni*4000
mu3.test = 90 + set_bi*3000 
mu_set.test = data.frame(plane=mu1.test,unimodal=mu2.test,bimodal=mu3.test)

# scale setting
sc1.test = 40 + (-2*x1.test + 3*x2.test)*0.3 
sc2.test = 30 + set_uni*4000 
sc3.test = 30 + set_bi*3000 
sc_set.test = data.frame(plane=sc1.test,unimodal=sc2.test,bimodal=sc3.test)

# shape setting
sh1.test = 0.1+(-2*x1.test + 3*x2.test)*0.004
sh2.test = set_uni*50 
sh3.test = set_bi*50 
sh_set.test = data.frame(plane=sh1.test,unimodal=sh2.test,bimodal=sh3.test)

# generate test data
loc.test = mu_set.test[,S_map[S_num,]$Var1]
sc.test =  sc_set.test[,S_map[S_num,]$Var2]
sh.test =  sh_set.test[,S_map[S_num,]$Var3]

xlist.test = list()
for (s in 1:ns){
  xlist.test[[s]] = rgev(10000,loc=loc.test[s], scale=sc.test[s], shape=sh.test[s])
}
fit.test = list()
for ( s in 1:ns){
  x = xlist.test[[s]]
  fit.test[[s]] = fgev(x)$estimate
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

# compute new test loss
hdist.test_list = list()
for (i in 1:length(result)){
  hdist_vec = c()
  for (l in 1:length(result[[1]])){
  tvec=result[[i]][[l]]
  loc.vec.reg = tvec[3+(1:p)]
  sc.vec.reg = tvec[(3+p)+(1:p)]
  sh.vec.reg = tvec[(3+2*p)+(1:p)]
  loc.vec.hat = tvec[1] + drop(Z.test%*%loc.vec.reg)
  sc.vec.hat = exp(tvec[2] + drop(Z.test%*%sc.vec.reg))
  sh.vec.hat = tvec[3] + drop(Z.test%*%sh.vec.reg)
  h_dist = 0 
    for (s in 1:ns){
      x = GEV(loc=loc.test[s],scale=sc.test[s],shape=sh.test[s]) # RobExtremes
      y = GEV(loc=loc.vec.hat[s],scale=sc.vec.hat[s],shape=sh.vec.hat[s])
      h_dist = h_dist + HellingerDist(x,y,smooth) # distrEx
    }
  hdist_vec[l] = h_dist
  }
  hdist.test_list[[i]] = hdist_vec
  cat(i,"/")
}

eval(parse(text = paste0("save.image(file =","'", paste0('./Rexport/RData_sgev3_simulation/newloss_scenario',S_num, '.RData',"')"))))

