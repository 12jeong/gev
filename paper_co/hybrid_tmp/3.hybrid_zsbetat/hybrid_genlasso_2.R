rm(list=ls()); gc()
setwd("~/GitHub/gev/")
source("./lib/pack.R")

# if(!require(lubridate)){install.packages("lubridate")}; require(lubridate)
if(!require(data.table)){install.packages("data.table")}; require(data.table)
if(!require(genlasso)) install.packages('genlasso'); require(genlasso) # genlasso 
if(!require(dplyr)){install.packages("dplyr")}; require(dplyr)
# source("miss/onestep_library.R")
load("./hybrid_genlasso/hybrid_genlasso_ns56.RData")

stn_data <- fread("./kma_data/stnlds.csv", stringsAsFactors = T)
data_stnlds <- stn_data %>% filter(duplicated(지점) == F) %>% select(지점, 지점명)
stnlds_name <-left_join(data.frame(stnlds=unique(Pr_46$stnlds)),data_stnlds,by=c("stnlds"="지점"))

p = ncol(newZ)/nt
ss = split.data.frame(Pr_46,Pr_46$stnlds)  # stnlds로 dataframe 쪼개서 list에 분배
xlist = lapply(ss,"[[","pr")               # 강수량(pr) 변수로만 이루어진 list 생성

# c(0,1000,2000,2500,3000,3500,3510,3518,3519,3520)
# for(lam_t in c(0,1000,2000,2100,2200,2300,2400,2500)){
lam_t = 3600
  est_beta = coef(out, lambda=lam_t)$beta
predictlist = lapply(1:ns,function(i) drop(zlist[[i]]%*%est_beta))

de_xlist =  lapply(1:ns,function(i) xlist[[i]] - predictlist[[i]] )
est_tmp = lapply(de_xlist, function(x) fgev(x)$estimate)
est_mulist = lapply(1:ns, function(i) est_tmp[[i]][1] + predictlist[[i]])

# par(mfrow=c(2,2))
# for(s in 1:ns){
#   plot(xlist[[s]], main =paste0("s:",s))
#   lines(predictlist[[s]],col="blue")
#   lines(est_mulist[[s]],col="red")
# }

set.seed(1025)
par(mfrow=c(2,2))
for(s in sample(1:ns,4)){
  plot(x=unique(Pr_46$obsyear),y=xlist[[s]], 
       main = stnlds_name[s,2], xlab="", ylab="연최대강수량")
  # lines(predictlist[[s]],col="blue")
  lines(x=unique(Pr_46$obsyear),y=est_mulist[[s]],col="red")
}

par(mfrow=c(1,2))
plot(unique(Pr_46$lat),unlist(lapply(predictlist,function(x) x[1])),xlab="경도",ylab="z_s%*%beta_t",main="1973")
plot(unique(Pr_46$lat),unlist(lapply(predictlist,function(x) x[45])),xlab="경도",ylab="z_s%*%beta_t",main="2017")
# plot3d(unique(Pr_46$lat),unique(Pr_46$long),unlist(lapply(predictlist,function(x) x[1])))

ll=0
for (s in 1:ns){
  # ll=ll-sum(dgev(x=xlist[[s]],loc=est_mulist[[s]], scale=est_tmp[[s]][2], shape=est_tmp[[s]][3],log=T))
  ll=ll-sum(dgev(x=xlist[[s]],loc=est_mulist[[s]], scale=sd(xlist[[s]]), shape=est_tmp[[s]][3],log=T))
}
ll

df1 = length(which(abs(newD%*%est_beta) > 1e-04)) +ns
AIC = 2*ll + 2*(df1+2*ns)
BIC = 2*ll + log(ns*nt)*(df1*2*ns)

cat (lam_t , ":", df1, ll,AIC,BIC,"\n")
}

lgev = function (x, loc = 0, scale = 1, shape = 0) 
{
  if (min(scale) <= 0) 
    return( - 1e+6)
  if (length(shape) != 1) 
    stop("invalid shape")
  x <- (x - loc)/scale
  if (shape == 0) 
    d <- log(1/scale) - x - exp(-x)
  # if (shape > 1)
  #   d <- -(1e+6)
  else {
    nn <- length(x)
    xx <- 1 + shape * x
    xxpos <- xx[xx > 0 | is.na(xx)]
    scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
    d <- numeric(nn)
    d[xx > 0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) - 
      (1/shape + 1) * log(xxpos)
    d[xx <= 0 & !is.na(xx)] <- -(1e+6)
    # d[shape>1] <- -(1e+6)
  }
  return(d)
}

l2gev_t = function (x, tvec_t)
{  
  tvec = tvec_t
  # loglikelihood
  v = - sum(lgev(x, loc = tvec[1], scale = tvec[2], shape = tvec[3]))
  
  return(v)
}

x=xlist[[s]]
theta1 = fgev(x)$estimate; theta1
-sum(dgev(x,loc=est_mulist[[s]],scale=theta1[2],shape=theta1[3],log=T))
theta2 = fgev(de_xlist[[s]])$estimate; theta2
-sum(dgev(x,loc=est_mulist[[s]],scale=theta2[2],shape=theta2[3],log=T))

x=de_xlist[[s]]
theta0=c(mean(x),sd(x),-0.1)
theta0[3] = fgev(x)$estimate[3]
theta0[2] = fgev(x)$estimate[2]
tvec_t = theta0
fgev(x)$estimate
GEVnewtonRaphson_step1(x=x,theta0=theta0)$root
ctr_list = list()
ctr_list$maxit = 1
optim(tvec_t, l2gev_t, x=x, control=ctr_list, method="BFGS")$par  
