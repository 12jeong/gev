rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")
load("./kma_data/Pr_46.RData")

set.seed(2019)
train_stnlds = sample(unique(Pr_46$stnlds),45)
train_46 = Pr_46 %>% filter(stnlds %in% train_stnlds)
test_46 = Pr_46 %>% filter(!stnlds %in% train_stnlds)
ss <- split.data.frame(train_46,train_46$stnlds)  # stnlds로 dataframe 쪼개서 list에 분배
xlist <- lapply(ss,"[[","pr")               # 강수량(pr) 변수로만 이루어진 list 생성
ns <- length(unique(train_46$stnlds))

# x : long (경도) , y : lat (위도) 
x_bsobj <- create.bspline.basis(range(train_46$long),breaks=quantile(train_46$long,prob = seq(0, 1, length = 3)))
y_bsobj <- create.bspline.basis(range(train_46$lat),breaks=quantile(train_46$lat,prob = seq(0, 1, length = 3)))

zlist <- list()
for (i in 1:ns){
  xbs <- eval.basis(unique(ss[[i]]$long),x_bsobj)
  ybs <- eval.basis(unique(ss[[i]]$lat),y_bsobj)
  tensorbs <- do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs))  # row-wise kronecker product
  zlist[[i]] <- tensorbs 
}
dim(zlist[[1]]) # (frist)stnlds 2D-splines tensor, nbasis = 5 x 5

# 2-D splines penlaty matrix
Fmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om <- Fmat+Gmat+2*Hmat

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

## example 
p=length(zlist[[1]])
Z = do.call("rbind",zlist)

lam_set = c(0,0.002,0.005, 0.02, 0.1, 0.3, 0.5,1, 5,100)
lam_len = length(lam_set)
table.grid = expand.grid(1:lam_len,1:lam_len,1:lam_len)

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

# 최종 선택 lambda로 추정
lambda = c(0.002,5,0.3)
result = gevreg_3m(xlist=xlist, zlist=zlist, lambda = lambda, Om=Om)
tvec = result$par

loc.vec.reg = tvec[3+(1:p)]
sc.vec.reg = tvec[(3+p)+(1:p)]
sh.vec.reg = tvec[(3+2*p)+(1:p)]
loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)

point.est = unlist(lapply(xlist,function(x) fgev(x)$estimate))
point.loc = point.est[3*(1:ns)-2]
point.sc =  point.est[3*(1:ns)-1]
point.sh =  point.est[3*(1:ns)]

par(mfrow=c(1,3))
plot(unique(train_46$lat),point.loc); points(unique(train_46$lat),loc.vec,col="red")
plot(unique(train_46$lat),point.sc); points(unique(train_46$lat),sc.vec,col="red")
plot(unique(train_46$lat),point.sh); points(unique(train_46$lat),sh.vec,col="red")
