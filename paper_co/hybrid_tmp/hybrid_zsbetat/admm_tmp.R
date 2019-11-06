
rm(list=ls()); gc()
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")
load("./kma_data/Pr_46.RData")

library(dplyr)
set.seed(1)
stnlds_tmp = sample(unique(Pr_46$stnlds),43)
Pr_46 = Pr_46 %>% filter(stnlds %in% stnlds_tmp)

ss = split.data.frame(Pr_46,Pr_46$stnlds)  # stnlds로 dataframe 쪼개서 list에 분배
xlist = lapply(ss,"[[","pr")               # 강수량(pr) 변수로만 이루어진 list 생성

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
Z=do.call("rbind",zlist)

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
tvec = rep(NA,p*ns+2*ns); length(tvec)
tmp_loc = c()
for ( i in 1:ns)
{
  x = xlist[[i]]
  fit = fgev(x)
  tmp_loc[i] = fit$estimate[1]
  tvec[p*ns+2*i-1] = fit$estimate[2]
  tvec[p*ns+2*i] = fit$estimate[3]
}
tmp_beta = lm(tmp_loc~Z-1)$coefficients
tvec[1:(p*ns)] = rep(tmp_beta,ns)

z_init = rep(0,nt-2) 
u_init = rep(1,nt-2)
z_init_list = list() ; u_init_list = list()
for (i in 1:ns){
  z_init_list[[i]] = z_init
  u_init_list[[i]] = u_init
}
dmatrix = mat_func(nt)
rho = 0.5

hygev_m = function(tvec, lam_s, xlist, zlist, Om, dmatrix, rho, z_init_list, u_init_list)  
{
  nS = length(xlist)
  nt = length(xlist[[1]])
  p = ncol(zlist[[1]])
  v1 = 0
  v2 = 0
  for ( i in 1: nS){
    x = xlist[[i]]; z = zlist[[i]]
    z_init = z_init_list[[i]]; u_init = u_init_list[[i]]
    loc.vec = c()
    for ( t in 1:nt){
      beta.t = tvec[(1:p)+p*(t-1)]
      loc.vec[t] = z%*%beta.t
    }
    sc.vec = tvec[p*ns+2*i-1]
    sh.vec = tvec[p*ns+2*i]
    
    # 모든 ns에 대하여 -log likelihood
    v1 = v1 - sum(lgev(x, loc = loc.vec, 
                       scale = sc.vec, shape = sh.vec))   
    # lagrangian term
    v2 = v2 + (rho/2)*sum(((dmatrix %*% loc.vec) - z_init + u_init)^2)
  }
  # tps penalty term
  vtps = 0 
  for (t in 1:nt){
    beta.t =tvec[(1:p)+p*(t-1)]
    vtps = vtps + t(beta.t)%*%Om%*%beta.t
  }
  v3 = lam_s * vtps
  
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

lam_s = 0
lam_t = 0.5

start_time = Sys.time()
# ADMM
maxiiter= 600
for (iiter in 1:maxiiter){
  # mu, sigma, kappa optim
  tvec = optim(tvec, hygev_m, lam_s=lam_s,
               xlist = xlist, zlist = zlist, Om = Om, dmatrix = dmatrix, 
               rho=rho, z_init_list = z_init_list, u_init_list = u_init_list, method = "BFGS",
               control = ctr_list)$par
  
  for (i in 1:ns) {
    xlist[[i]]; z = zlist[[i]]
    z_init = z_init_list[[i]]
    u_init = u_init_list[[i]]
    
    old_mu = c()
    for ( t in 1:nt){
      beta.t = tvec[(1:p)+p*(t-1)]
      old_mu[t] = z%*%beta.t
    }
    
    # z update in ADMM
    tmp_z = drop(func_z(dmatrix = dmatrix, mu = old_mu, u = u_init, lam = lam_t, rho = rho))
    
    AA = rho*t(A)%*%B
    sk = base::norm(drop(AA%*%(tmp_z -z_init)),"2")   # dual
    rk = base::norm(drop(A%*%old_mu + B%*%tmp_z),"2") # primal
    
    # u update in ADMM
    tmp_u = drop(func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init))
    u_init_list[[i]] = tmp_u ; z_init_list[[i]] = tmp_z
    
  }
  # if (rk<epri & sk <edul) break
  # if (iiter %% 100 == 0){
  mid_time = Sys.time()
  cat("*****iter:::", iiter, "/log:", mid_time-start_time  ,"\n")
  # }
}
end_time = Sys.time()

run_time = end_time - start_time

par (mfrow=c(2,5))
for (s in 1:ns){
  z=zlist[[s]]
  for ( t in 1:nt){
    beta.t = tvec[(1:p)+p*(t-1)]
    old_mu[t] = z%*%beta.t
  }
  plot(old_mu,type="l")
}

