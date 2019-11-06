rm(list=ls()); gc()
source("./lib/pack.R")
load("./Pr_46.RData")

if(!require(genlasso)) install.packages('genlasso'); require(genlasso) # genlasso 
if(!require(metaSEM)) install.packages('metaSEM'); require(metaSEM)    # bdiagMat
# if(!require(dplyr)) install.packages('dplyr'); require(dplyr)

# set.seed(119)
# Pr_46 = Pr_46 %>% filter(stnlds %in% sample(unique(Pr_46$stnlds),27))
# Pr_46 = Pr_46 %>% filter(obsyear < 2010 )
ns = length(unique(Pr_46$stnlds)) ;ns
nt = length(unique(Pr_46$obsyear)) ;nt

mat_func = function(n) {
  m = matrix(0,n-2,n)   ## dmatrix : (n-2)*n 
  for (i in 1:(n-2)) {
    m[i,i] = 1
    m[i,i+1] = -2
    m[i,i+2] = 1
  }
  return(m)
}

dmatrix = mat_func(nt)

x_bsobj <- create.bspline.basis(range(Pr_46$long),
                                breaks=quantile(Pr_46$long,prob = seq(0, 1, length = 3)))
y_bsobj <- create.bspline.basis(range(Pr_46$lat),
                                breaks=quantile(Pr_46$lat,prob = seq(0, 1, length = 3)))
zlist <- list()
Dlist <- list()
for (i in 1:ns){
  xbs <- eval.basis(unique(Pr_46$long)[i],x_bsobj)
  ybs <- eval.basis(unique(Pr_46$lat)[i],y_bsobj)
  tensorbs <- do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs))  # row-wise kronecker product
  ztmp = list()
  for (ii in 1:nt){
    ztmp[[ii]] = tensorbs
  }
  zlist[[i]] = bdiagMat(ztmp) 
  Dlist[[i]] <- dmatrix %*% zlist[[i]]
}
dim(zlist[[1]])
newZ = do.call("rbind",zlist)
dim(newZ)

Y= Pr_46$pr
newD = do.call("rbind",Dlist)
dim(newD)

out = genlasso(y=Y,X=newZ,D=newD,maxsteps=5500)
summary(out$lambda)


