rm( list = ls()); gc()
setwd("C:\\Users\\UOS\\Documents\\GITHUB\\gev")
source("sgevlibrary.R")
if(!require(fda)) install.packages('fda'); library(fda)

# zlist - 2-D basis matrix by each location
x_bsobj = create.bspline.basis(range(x1),breaks=quantile(x1,prob = seq(0, 1, length = 3)))
y_bsobj = create.bspline.basis(range(x2),breaks=quantile(x2,prob = seq(0, 1, length = 3)))
zlist = list()
for (i in 1:ns){
  xbs = eval.basis(df$row[i],x_bsobj)
  ybs = eval.basis(df$col[i],y_bsobj)
  tensorbs = do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs)) 
  zlist[[i]] = tensorbs 
}
# dim(zlist[[1]]) # (frist)stnlds 2D-splines tensor, nbasis = df x df

# Omega matrix
Fmat = kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat = kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat = kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om = Fmat+Gmat+2*Hmat

# optim control
optim_controlList = list()
optim_controlList$maxit = 1e+3

# design matrix for v3 (mu_0 stationary)
matfunc = function(ns){
  matt = matrix(0,nrow=ns*(ns-1),ncol=ns)
  k = 1
  for (i in (1: (ns-1))){
    for (j in ((i+1) :ns)){
      matt[k,i] = 1
      matt[k,j] = -1
      k = k+1
    }  
  }
  mat = t(matt) %*% matt
  return(mat)
}
mat = matfunc(ns)

result1 = gevreg_m(xlist,zlist,
          lambda = 1, lambda2=1000, Om=Om, mat=mat, method="B-spline")


p <- ncol(zlist[[1]])
xbss <- eval.basis(x1[s_row],x_bsobj)
ybss <- eval.basis(x2[s_col],y_bsobj)
tensorbss <- do.call('cbind', lapply(1:ncol(xbss), function(i) xbss[, i] * ybss)) 
z1 <- 100+tensorbss %*% tail(result1,p)

library(rgl)
plot3d(df$row,df$col,z1)
plot3d(df$row,df$col,df$z)

######### k-folds cross validation ##########
lambdaset = seq(0,1,length=11)

## 순서대로 
ksize = 5 # number of K
gsize = length(xlist[[1]])/ksize
folds = rep(c(1:ksize),each=gsize)
# train$k = rep(folds,times=ns) ; train$k

## 랜덤하게 
# set.seed(seed)
# ksize = 5 # number of K
# folds = sample(cut(1:40,breaks=ksize,labels=FALSE),40)


# eval(parse(text = paste0('result_par', file_idx, ' = list()')))
result_par = list()
for (i in 1:ksize){
  trainx = lapply(1:ns, function(k) xlist[[k]][folds!=i] )
  trainz = lapply(1:ns, function(k) zlist[[k]][1:(40-sum(folds==i)),])
  tmp <- gevreg_m(xlist=trainx,zlist=trainz,
                  lambda = lambdaset[lambda_idx], lambda2=1, Om=Om, mat=mat, method="B-spline")
  result_par[[i]] = tmp
  # eval(parse(text = paste0('result_par', file_idx, '[[i]] = tmp')))
}

# calculate validation loss
result_loss = c()
for ( i in 1:ksize){
  v = c()
  est_z = tail(result_par[[i]], dim(zlist[[1]])[2])
  for( s in 1:ns){
    est_s = result_par[[i]][(1:3)+(3*s-3)]
    gsize = length(xlist[[s]][folds==i])
    v[s] = sum(lossfun(x=xlist[[s]][folds==i], loc=c(est_s[1]+zlist[[s]][1,]%*%est_z), scale=est_s[2], shape=est_s[3]))
  }
  result_loss[i] = sum(v)
}

# save.file
eval(parse(text = paste0('result_par', file_idx, ' = result_par')))
eval(parse(text = paste0('result_loss', file_idx, ' = result_loss')))

# lambda 값 하나당 5folds 실행
eval(parse(text = paste0('save(result_par', file_idx, ',', 'result_loss', file_idx ,
                         ',', paste0("file =","'", paste0('result_train2',file_idx, '.RData',"')")))))
