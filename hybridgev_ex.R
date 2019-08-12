rm(list=ls()); gc()
setwd("~/GITHUB/gev")
source("./lib/pack.R")

# ns = 2, T = 46, p = 25

load("./kma_data/Pr_46.RData")
ss = split.data.frame(Pr_46,Pr_46$stnlds)  # stnlds로 dataframe 쪼개서 list에 분배
xlist = lapply(ss,"[[","pr")               # 강수량(pr) 변수로만 이루어진 list 생성
# xlist = list(xlist[[1]],xlist[[2]])

x_bsobj <- create.bspline.basis(range(Pr_46$long),breaks=quantile(Pr_46$long,prob = seq(0, 1, length = 3)))
y_bsobj <- create.bspline.basis(range(Pr_46$lat),breaks=quantile(Pr_46$lat,prob = seq(0, 1, length = 3)))
zlist <- list()
for (i in 1:length(xlist)){
  xbs <- eval.basis(ss[[i]]$long[1],x_bsobj)
  ybs <- eval.basis(ss[[i]]$lat[1],y_bsobj)
  tensorbs <- do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs))  # row-wise kronecker product
  zlist[[i]] <- tensorbs 
}
dim(zlist[[1]])

# 2-D splines penlaty matrix
Fmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat <- kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om <- Fmat+Gmat+2*Hmat


lgev = function (x, loc = 0, scale = 1, shape = 0) 
{
  if (min(scale) <= 0) 
    return( - 1e+6)
  if (length(shape) != 1) 
    stop("invalid shape")
  x <- (x - loc)/scale
  if (shape == 0) 
    d <- log(1/scale) - x - exp(-x)
  else {
    nn <- length(x)
    xx <- 1 + shape * x
    xxpos <- xx[xx > 0 | is.na(xx)]
    scale <- rep(scale, length.out = nn)[xx > 0 | is.na(xx)]
    d <- numeric(nn)
    d[xx > 0 | is.na(xx)] <- log(1/scale) - xxpos^(-1/shape) - 
      (1/shape + 1) * log(xxpos)
    d[xx <= 0 & !is.na(xx)] <- -(1e+6)
  }
  return(d)
}

#z_function
func_z <- function(dmatrix,mu,u,lam,rho){
  z = ifelse(abs(dmatrix %*% mu + u) > (lam/rho) ,
             (dmatrix %*% mu) + u - sign(u + (dmatrix %*% mu)) * (lam /rho), 0)
  return(z)
}


#u_function
func_u <- function(dmatrix,mu,z,u) {
  u <- u + (dmatrix %*% mu) - z
  return(u)
}

#dmatrix
mat_func = function(n) {
  m = matrix(0,n-2,n)   ## dmatrix : (n-2)*n 
  for (i in 1:(n-2)) {
    m[i,i] = 1
    m[i,i+1] = -2
    m[i,i+2] = 1
  }
  return(m)
}

# initial value
p = ncol(zlist[[1]]); p
ns = length(xlist) ; ns
nt = length(xlist[[1]]) ; nt
tvec = rep(0, nt+ 2*ns+ p) ; length(tvec)
tmp_loc = 0
for ( i in 1:ns)
{
  x = xlist[[i]]
  fit = fgev(x)
  tmp_loc = tmp_loc + fit$estimate[1]
  tvec[nt+2*i-1] = fit$estimate[2]
  tvec[nt+2*i] = fit$estimate[3]
}
tvec[1:nt] = tmp_loc/ns

z_init = rep(0,nt-2) 
u_init = rep(1,nt-2)

dmatrix = mat_func(nt)
lam_s = 0.5
lam_t = 5
rho = 0.5

hygev_m = function(tvec, lamb_t, lamb_s, xlist, zlist, Om, dmatrix, rho, z_init, u_init)  
{
  p = ncol(zlist[[1]])
  v1 = 0
  for ( i in 1: ns){
    x = xlist[[i]]
    z = zlist[[i]]
    loc.vec.reg = drop(z%*%tail(tvec,p)) # z_s %*% beta : scalar
    loc.vec = tvec[1:nt] 
    sc.vec = tvec[nt+2*i-1]
    sh.vec = tvec[nt+2*i]
    v1 = v1 - sum(lgev(x, loc = loc.vec + loc.vec.reg, 
                       scale = sc.vec, shape = sh.vec))     # 모든 ns에 대하여 -log likelihood
    }
  # tps penalty term
  v2 = lam_s *t(tail(tvec,p))%*%Om%*%tail(tvec,p)
  # lagrangian term
  v3 = (rho/2)*sum(((dmatrix %*% loc.vec) - z_init + u_init)^2)
  v = v1 + v2 + v3
  return(v)
}

# parameter used in optim()
ctr_list = list()
ctr_list$maxit = 20
ctr_list$reltol = 1e-6

# resolution in ADMM
epri = 1e-4
edul = 1e-4
A = dmatrix
B = -diag(nt-2)

start_time = Sys.time()
# ADMM
for (iter in 1:10000){
  if (iter %% 100 == 0){
    cat("*****iter:::", iter, "\n")
  }
  #mu, sigma, kappa optim
  old_tvec = optim(tvec, hygev_m, lamb_t=lamb_t, lamb_s=lamb_s,
                   xlist = xlist, zlist = zlist, Om = Om, dmatrix = dmatrix, 
                   rho=rho, z_init = z_init, u_init = u_init, method = "BFGS",
                   control = ctr_list)$par
  
  old_mu = old_tvec[1:nt]
  
  # z update in ADMM
  tmp_z = func_z(dmatrix = dmatrix, mu = old_mu, u = u_init, lam = lam_t, rho = rho)
  
  AA = rho*t(A)%*%B
  sk = norm(drop(AA%*%(tmp_z -z_init)),"2")
  rk = norm(drop(A%*%old_mu + B%*%tmp_z),"2")
  
  if (rk<epri & sk <edul) break
  
  # u update in ADMM
  tmp_u = func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init)
  
  u_init = tmp_u ; z_init = tmp_z
  tvec = old_tvec 
}
end_time = Sys.time()

start_time - end_time

tvec[1:nt]
tvec[nt+(2*(1:ns)-1)]
tvec[nt+(2*(1:ns))]
tail(tvec,p)
plot(tvec[1:nt] + rep(zlist[[1]]%*% tail(tvec,p),nt))
plot(tvec[1:nt] + rep(zlist[[56]]%*% tail(tvec,p),nt))


p <- ncol(zlist[[1]])
x <- unique(Pr_46$long)
y <- unique(Pr_46$lat)
xbss <- eval.basis(x,x_bsobj)
ybss <- eval.basis(y,y_bsobj)
tensorbss <- do.call('cbind', lapply(1:ncol(xbss), function(i) xbss[, i] * ybss)) 
z <- tensorbss %*% tail(tvec,p)
library(rgl)
plot3d(x,y,z)
