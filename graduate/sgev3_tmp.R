rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")
load("./kma_data/Pr_46.RData")

# set.seed(518)
# train_stnlds = sample(unique(Pr_46$stnlds),length(unique(Pr_46$stnlds))*0.8)
# train_46 = Pr_46 %>% filter(stnlds %in% train_stnlds)
# test_46 = Pr_46 %>% filter(!stnlds %in% train_stnlds)
# ss <- split.data.frame(train_46,train_46$stnlds)
ss <- split.data.frame(Pr_46,Pr_46$stnlds)  # stnlds로 dataframe 쪼개서 list에 분배
xlist <- lapply(ss,"[[","pr")               # 강수량(pr) 변수로만 이루어진 list 생성
ns <- length(unique(Pr_46$stnlds))
# ns <- length(unique(train_46$stnlds))

# x : long (경도) , y : lat (위도) 
x_bsobj <- create.bspline.basis(range(Pr_46$long),breaks=quantile(Pr_46$long,prob = seq(0, 1, length = 3)))
y_bsobj <- create.bspline.basis(range(Pr_46$lat),breaks=quantile(Pr_46$lat,prob = seq(0, 1, length = 3)))

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

lam_set = c(0,0.002,0.005, 0.02, 0.1, 0.3, 0.5)
lam_len = length(lam_set)
table.grid = expand.grid(1:lam_len,1:lam_len,1:lam_len)

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3


##### train set - AIC #######

s_time = Sys.time()
result_list = list()
for (i in 1:lam_len^3){
  lambda = c(lam_set[table.grid[i,1]],lam_set[table.grid[i,2]],lam_set[table.grid[i,3]])
  result = gevreg_3m(xlist=xlist, zlist=zlist, lambda = lambda, Om=Om)
  result_list[[i]] = result$par
  cat(i,"/")
}

e_time = Sys.time()
run_time = e_time - s_time

# save.image("~/GitHub/gev/sgev3_train1204.RData")
save.image("~/GitHub/gev/sgev3_real1209.RData")

######## for cross validation ######### 

# set.seed(930)
# ksize = 5 
# folds = sample(cut(1:45,breaks=ksize,labels=FALSE),45)
# 
# s_time = Sys.time()
# result_list = list()
# for (k in 1:ksize){
#   trainx = xlist[folds!=k]
#   trainz = zlist[folds!=k]
#   result_tmp = list()
#   for (i in 1:lam_len^3){
#     lambda = c(lam_set[table.grid[i,1]],lam_set[table.grid[i,2]],lam_set[table.grid[i,3]])
#     result = gevreg_3m(xlist=trainx, zlist=trainz, lambda = lambda, Om=Om)
#     result_tmp[[i]] = result$par
#     cat(i,"/")
#   }
#   result_list[[k]] = result_tmp
#   cat("\n")
# }
# 
# e_time = Sys.time()
# run_time = e_time - s_time

# save.image("~/GitHub/gev/sgev3_cvtrain1001.RData")

