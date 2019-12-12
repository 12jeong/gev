rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")

# setwd("C:/Users/UOS/Downloads")
# source("C:/Users/UOS/Documents/GITHUB/gev/lib/sgev3library.R")
# source("C:/Users/UOS/Documents/GITHUB/gev/lib/pack.R")
library(distrEx)
library(RobExtremes) 

S_num = 1

eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/AIC_scenario',S_num, '.RData',"')"))))
# eval(parse(text = paste0("load(file =","'", paste0('./AIC_scenario',S_num, '.RData',"')"))))

# true paramter
loc = mu_set[,S_map[S_num,]$Var1]
sc =  sc_set[,S_map[S_num,]$Var2]
sh =  sh_set[,S_map[S_num,]$Var3]
Sys.time()
rmse_list = list()
hdist_list = list()
for (i in 1:length(result)){
  rmse_vec = c()
  hdist_vec = c()
  for (l in 1:length(result[[1]])){
    tvec=result[[i]][[l]]
    loc.vec.reg = tvec[3+(1:p)]
    sc.vec.reg = tvec[(3+p)+(1:p)]
    sh.vec.reg = tvec[(3+2*p)+(1:p)]
    loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
    sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
    sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
    
    h_dist = 0 
    for (s in 1:ns){
      x = GEV(loc=loc[s],scale=sc[s],shape=sh[s]) # RobExtremes
      y = GEV(loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s])
      h_dist = h_dist + HellingerDist(x,y,smooth) # distrEx
    }

    rmse_vec[l] = base::norm(c(loc,sc,sh)-c(loc.vec,sc.vec,sh.vec), "2")/sqrt(3*ns)
    hdist_vec[l] = h_dist
  }
  rmse_list[[i]] = rmse_vec
  hdist_list[[i]] = hdist_vec
  cat(i," ")
}

rmse_mat = do.call('rbind',rmse_list)
hdist_mat = do.call('rbind',hdist_list)


eval(parse(text = paste0("save.image(file =","'", paste0('./Rexport/RData_sgev3_simulation/testloss_scenario',S_num, '.RData',"')"))))
# eval(parse(text = paste0("save.image(file =","'", paste0('./testloss_scenario',S_num, '.RData',"')"))))