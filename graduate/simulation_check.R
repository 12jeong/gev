# To verify that the initial value (fgev) is working

rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")

# scenario Mapping
mapa = c("Plane","Unimodal","Bimodal")
S_map = expand.grid(1:3, 1:3, 1:3)
S_map2 = expand.grid(mu=mapa, sc=mapa, sh=mapa)


# sampling X coordinates 
ns = 30
xyrange = c(-10,10)
set.seed(101)
x1 = runif(ns,xyrange[1],xyrange[2])
x2 = runif(ns,xyrange[1],xyrange[2])

# surface base setting
mean_vec = c(0,0) 
sig_mat = matrix(c(30,0,0,30),nrow=2)
set.seed(102)
set_uni = dmvnorm(cbind(x1,x2), mean=mean_vec, sigma=sig_mat)
mean_vec1 = c(5,0); mean_vec2 = c(-5,0) 
sig_mat = matrix(c(10,0,0,10),nrow=2)
set.seed(103)
set_bi = 0.4*dmvnorm(cbind(x1,x2),mean=mean_vec1, sigma=sig_mat*1) +0.6*dmvnorm(cbind(x1,x2),mean=mean_vec2, sigma=sig_mat*2)

# location setting
mu1 = 100 + (-2*x1 + 3*x2) 
mu2 = 90 + set_uni*4000
mu3 = 90 + set_bi*3000 
mu_set = data.frame(plane=mu1,unimodal=mu2,bimodal=mu3)

# scale setting
sc1 = 40 + (-2*x1 + 3*x2)*0.3 
sc2 = 30 + set_uni*4000 
sc3 = 30 + set_bi*3000 
sc_set = data.frame(plane=sc1,unimodal=sc2,bimodal=sc3)

# shape setting
# set.seed(104)
# sh1 = runif(ns,0,0.3)
sh1 = 0.1+(-2*x1 + 3*x2)*0.004
range(sh1)
plot(x1,sh1)

sh2 = set_uni*50 
sh3 = set_bi*50 
sh_set = data.frame(plane=sh1,unimodal=sh2,bimodal=sh3)



for (S_num in 1:9){
  n = 100
  loc = mu_set[,S_map[S_num,]$Var1]
  sc =  sc_set[,S_map[S_num,]$Var2]
  sh =  sh_set[,S_map[S_num,]$Var3]
  
  for ( seed in 1:100){
    set.seed(seed)
    cat(seed,"/")
    xlist = list()
    for (s in 1:ns){
      xlist[[s]] = rgev(n,loc=loc[s], scale=sc[s], shape=sh[s])
    }
    
    for ( i in 1:ns)
    {
      x = xlist[[i]]
      fit = fgev(x)
    }
  }
}
