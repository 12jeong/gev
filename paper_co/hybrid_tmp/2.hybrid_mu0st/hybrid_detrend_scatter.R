rm(list=ls()); gc()
setwd("~/GITHUB/gev")
source("./lib/pack.R")

# ns = 2, T = 46, p = 25
load("./kma_data/Pr_46.RData")
ss = split.data.frame(Pr_46,Pr_46$stnlds)  # stnlds로 dataframe 쪼개서 list에 분배
xlist = lapply(ss,"[[","pr")               # 강수량(pr) 변수로만 이루어진 list 생성
# xlist = list(xlist[[1]],xlist[[2]])

x_bsobj <- create.bspline.basis(range(Pr_46$long),
                                breaks=quantile(Pr_46$long,prob = seq(0, 1, length = 3)))
y_bsobj <- create.bspline.basis(range(Pr_46$lat),
                                breaks=quantile(Pr_46$lat,prob = seq(0, 1, length = 3)))
zlist <- list()
for (i in 1:length(xlist)){
  xbs <- eval.basis(ss[[i]]$long[1],x_bsobj)
  ybs <- eval.basis(ss[[i]]$lat[1],y_bsobj)
  tensorbs <- do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs))  # row-wise kronecker product
  zlist[[i]] <- tensorbs 
}
dim(zlist[[1]])
Z = do.call('rbind',zlist)

# 2-D splines penlaty matrix
Fmat <- kronecker(eval.penalty(x_bsobj,Lfdobj=2),eval.penalty(y_bsobj,Lfdobj=0))
Gmat <- kronecker(eval.penalty(x_bsobj,Lfdobj=0),eval.penalty(y_bsobj,Lfdobj=2))
Hmat <- kronecker(eval.penalty(x_bsobj,Lfdobj=1),eval.penalty(y_bsobj,Lfdobj=1))
Om <- Fmat+Gmat+2*Hmat

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


lgev = function (x, loc = 0, scale = 1, shape = 0) 
{
  if (min(scale) <= 0) 
    return( - 1e+6)
  if (length(shape) != 1) 
    stop("invalid shape")
  x <- (x - loc)/scale
  if (shape == 0) 
    d <- log(1/scale) - x - exp(-x)
  if (shape > 1)
    d <- -(1e+6)
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

l2gev_t = function (x, tvec_t, dmatrix, rho, z_init, u_init)
{  
  tvec = tvec_t
  nt = length(x)
  loc.vec = tvec[1:nt]
  # loglikelihood
  v1 = - sum(lgev(x, loc = loc.vec, scale = tvec[nt+1], shape = tvec[nt+2]))
  
  # lagrangian term
  v2 = (rho/2)*sum(((dmatrix %*% loc.vec) - z_init + u_init)^2)
  v = v1 + v2
  return(v)
}

l2gev_s = function (xlist, zlist, tvec_s, lambda, Om)
{
  tvec = tvec_s
  nS = length(xlist)
  p = length(zlist[[1]])
  v1 = 0
  for ( i in 1: length(xlist))
  {
    x = xlist[[i]]
    z = zlist[[i]]
    loc.vec.reg = drop(z%*%tail(tvec, p))
    sc.vec = tvec[2*i-1]
    sh.vec = tvec[2*i]
    v1 = v1 - sum(lgev(x, loc = loc.vec.reg, 
                       scale = sc.vec, shape = sh.vec))   # loss : +(-loglikelihood)
  }
  # regularization
  v2 = lambda*t(tail(tvec,p))%*%Om%*%tail(tvec,p)
  v = v1 + v2
  
  return(v)
}

# initial value
p = ncol(zlist[[1]]); p
nS = length(xlist) ; nS
nt = length(xlist[[1]]) ; nt
tvec0 = rep(0, nt+ 2*nS+ p) ; length(tvec0)
tmp_loc = 0
for ( i in 1:nS)
{
  x = xlist[[i]]
  fit = fgev(x)
  tmp_loc = tmp_loc + fit$estimate[1]
  tvec0[nt+2*i-1] = fit$estimate[2]
  tvec0[nt+2*i] = fit$estimate[3]
}
tvec0[1:nt] = tmp_loc/nS

z_init = rep(0,nt-2) 
u_init = rep(1,nt-2)

dmatrix = mat_func(nt)
lam_s = 0
lam_t = 0
rho = 0.5  

# parameter used in optim()
ctr_list = list()
ctr_list$maxit = 20
ctr_list$reltol = 1e-6

# resolution in ADMM
epri = 1e-3
edul = 1e-3
A = dmatrix
B = -diag(nt-2)

s_time = Sys.time()
# ADMM
tvec_s = tvec0[-c(1:nt)]
t_tvec_list = lapply(1:length(xlist),function(i) c(tvec0[1:nt],tvec0[nt+2*i-1],tvec0[nt+2*i]))
old_mu_list = lapply(1:length(xlist), function(i) t_tvec_list[[i]][1:nt])

tlist = list(); slist=list()
for (iterr in 1:10){
  # spatial part
  xlist_s = lapply(1:nS, function(i) xlist[[i]]-old_mu_list[[i]])
  tvec_s = optim(tvec_s, l2gev_s, lambda=lam_s, xlist=xlist_s, zlist=zlist, Om=Om, method="BFGS")$par
  for(s in 1:nS){
    t_tvec_list[[s]][c(nt+1,nt+2)] = tvec_s[c(2*s-1,2*s)] 
  }
  xlist_t = lapply(1:nS, function(i) xlist[[i]]-drop(zlist[[i]]%*%tail(tvec_s,p)))
  for(s in 1:nS){
    x = xlist_t[[s]]
    tvec_t = t_tvec_list[[s]][1:(nt+2)]
    for (iter in 1:10000){
      #mu, sigma, kappa optim
      tvec_t =  optim(tvec_t, l2gev_t, x=x, dmatrix = dmatrix,
                        rho=rho, z_init = z_init, u_init = u_init, method = "BFGS",
                        control = ctr_list)$par
      # tvec[1:(nt+2*nS)] = tvec_top
      old_mu = tvec_t[1:nt]
      
      # z update in ADMM
      tmp_z = func_z(dmatrix = dmatrix, mu = old_mu, u = u_init, lam = lam_t, rho = rho)
      
      AA = rho*t(A)%*%B
      sk = base::norm(drop(AA%*%(tmp_z -z_init)),"2")
      rk = base::norm(drop(A%*%old_mu + B%*%tmp_z),"2")
      
      if (rk<epri & sk <edul) break
      
      # u update in ADMM
      tmp_u = func_u(dmatrix = dmatrix , mu = old_mu, z = tmp_z, u = u_init)
      
      u_init = tmp_u ; z_init = tmp_z
    }
    t_tvec_list[[s]][1:(nt+2)] = tvec_t
    tvec_s[c(2*s-1,2*s)] = tvec_t[c(nt+1,nt+2)]
    old_mu_list[[s]] = old_mu
    cat("*****iter:::", iterr, "/", s,"/", iter, "sk::", sk, "  rk::", rk, "\n")
  }
  # 수렴기준잡기
  slist[[iterr]] = tvec_s
  tlist[[iterr]] = t_tvec_list
}
e_time = Sys.time()
e_time - s_time

# save.image("detrend_0916_tmp.RData")

# matrix(unlist(slist),ncol=10)
# matrix(unlist(lapply(1:10,function(i) tlist[[i]][[1]])),ncol=10) # 첫번째 지역?
