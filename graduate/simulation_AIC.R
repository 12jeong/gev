rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")

# setwd("C:/Users/UOS/Downloads")
# source("C:/Users/UOS/Documents/GITHUB/gev/lib/sgev3library.R")
# source("C:/Users/UOS/Documents/GITHUB/gev/lib/pack.R")

S_num = 18

eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/result_scenario',S_num, '.RData',"')"))))
# eval(parse(text = paste0("load(file =","'", paste0('./result_scenario',S_num, '.RData',"')"))))

p=length(zlist[[1]])

df.1 = c()
df.2 = c()
df.3 = c()
for (i in 1:length(lam_set1)){
  lam_s1 = lam_set1[i] ; lam_s2 = lam_set2[i] ; lam_s3 = lam_set3[i]
  Hatmat1 = Z%*%solve(t(Z)%*%Z +lam_s1*Om+diag(1e-08,nrow(Om)))%*%t(Z)
  Hatmat2 = Z%*%solve(t(Z)%*%Z +lam_s2*Om+diag(1e-08,nrow(Om)))%*%t(Z)
  Hatmat3 = Z%*%solve(t(Z)%*%Z +lam_s3*Om+diag(1e-08,nrow(Om)))%*%t(Z)
  df.1[i] =  sum(diag(Hatmat1))
  df.2[i] =  sum(diag(Hatmat2))
  df.3[i] =  sum(diag(Hatmat3))
}
DF = rowSums(as.matrix(expand.grid( df.1, df.2, df.3)))


nll_list = list()
for (i in 1:length(result)){
  set.seed(i)
  xlist = list()
  for (s in 1:ns){
    xlist[[s]] = rgev(n,loc=loc[s], scale=sc[s], shape=sh[s])
  }
  nll = c()
  for (l in 1:length(result[[1]])){
    tvec=result[[i]][[l]]
    loc.vec.reg = tvec[3+(1:p)]
    sc.vec.reg = tvec[(3+p)+(1:p)]
    sh.vec.reg = tvec[(3+2*p)+(1:p)]
    loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
    sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
    sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
    
    ll = 0
    for (s in 1:length(xlist)){
      ll = ll - sum(dgev(x=xlist[[s]],loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s], log=TRUE))
    }
    nll[l] = ll
  }
  nll_list[[i]] = nll
}

AIC_list = lapply(nll_list,function(x) 2*x+2*DF)
lam.grid2 = expand.grid(lam_set1,lam_set2,lam_set3)
lam.min.table = table(unlist(lapply(AIC_list,function(x) which.min(x))))
lam.min.vec = as.numeric(rownames(lam.min.table))
lam.grid2[lam.min.vec,] # 선택된 lambda 범위 

eval(parse(text = paste0("save.image(file =","'", paste0('./Rexport/RData_sgev3_simulation/AIC_scenario',S_num, '.RData',"')"))))
# eval(parse(text = paste0("save.image(file =","'", paste0('./AIC_scenario',S_num, '.RData',"')"))))
