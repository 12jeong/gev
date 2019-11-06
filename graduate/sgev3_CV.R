rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")
# load("./sgev3_cvtrain.RData")
load("./sgev3_cvtrain1001.RData")
library(distrEx)
library(RobExtremes)
library(statip)

point.est = unlist(lapply(xlist,function(x) fgev(x)$estimate))
point.loc = point.est[3*(1:ns)-2]
point.sc =  point.est[3*(1:ns)-1]
point.sh =  point.est[3*(1:ns)]

kk=1 # validation folds index
ll=1 # lambda index

error_list = list()
for (kk in 1:k){
  testZ = do.call("rbind",zlist[folds==kk])
  
  error_vec = c()
  for (ll in 1:length(result_list[[kk]])){
    tvec = result_list[[kk]][[ll]]
    loc.vec.reg = tvec[3+(1:p)]
    sc.vec.reg = tvec[(3+p)+(1:p)]
    sh.vec.reg = tvec[(3+2*p)+(1:p)]
    loc.vec = tvec[1] + drop(testZ%*%loc.vec.reg)
    sc.vec = exp(tvec[2] + drop(testZ%*%sc.vec.reg))
    sh.vec = tvec[3] + drop(testZ%*%sh.vec.reg)
    
    theta.hat = c(loc.vec,sc.vec,sh.vec)
    theta = c(point.loc[folds==kk],point.sc[folds==kk],point.sh[folds==kk])
    
    # tl = 4 # test location index
    error_tmp = 0
    for (tl in 1:sum(folds==kk)){
      x = c(unlist(xlist[folds==kk][tl]))
      # validation error in kk th folds
      # error_tmp = error_tmp + sum(dgev(x=x, loc=loc.vec[tl], scale=sc.vec[tl], shape=sh.vec[tl], log=T))
      # y = rgev(10000,loc=loc.vec[tl],scale=sc.vec[tl],shape=sh.vec[tl])
      # error_tmp = error_tmp + hellinger(x,y) 
      e1 = GEV(loc=loc.vec[tl],scale=sc.vec[tl],shape=sh.vec[tl])
      e2 = GEV(loc=point.loc[folds==kk][tl],scale=point.sc[folds==kk][tl],shape=point.sh[folds==kk][tl])
      error_tmp = error_tmp + HellingerDist(e1,e2,"smooth")
      # error_tmp = error_tmp + HellingerDist(x,e1,"smooth")
      # SSE = sum((theta-theta.hat)^2) # validation error in kk th folds
    }
    error_vec[ll] = error_tmp
  }
  error_list[[kk]] = error_vec
}

# save.image("result_sgev3_CV.RData")

error_mat = do.call("rbind",error_list)

apply(error_mat,1,which.min)
table.grid[582,]
table.grid[52,]
table.grid[482,]
table.grid[782,]
table.grid[83,]

which.min(colMeans(error_mat))
table.grid[652,]
plot(colMeans(error_mat),type="l")
abline(v=1:10*100,lty=2)
abline(v=1:100*10,lty=2,col="blue")

lam_set[2]
lam_set[9]
lam_set[6]


# x=rnorm(100)
# y=Norm(0,1)
# HellingerDist(x,y,"smooth")
# HellingerDist(x,y)
