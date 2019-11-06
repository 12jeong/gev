rm(list=ls())
setwd("~/GitHub/gev")
source("./lib/pack.R")
source("./lib/sgevlibrary.R")

# lambda setting
lam_t = 0
lam_s = seq(0,1,length=11)

# data generation
set.seed(1)
nobs = 100
ns = 50

## 1.stationay
m0 = rep(100,nobs)
## 2. increasing
# m0 = seq(1,10,length.out = nobs)
## 3.sin
# y = sin(seq(0,2*pi,length = nobs))*20; m0 = y +100

## setting2. unimodal
xyrange = c(-10,10)
nBS = 3
n = 30
x1 = seq(xyrange[1],xyrange[2],length.out=n)
x2 = seq(xyrange[1],xyrange[2],length.out=n)
mu = c(0,0)
sig = matrix(c(30,0,0,30),nrow=2)
fx = outer(x1, x2, function(x1,x2) {
  dmvnorm(cbind(x1,x2), mean=mu, sigma=sig)
})

set.seed(2)
s_ind = sort(sample(n^2,ns))
zval = rep(0, ns)
k = 1
for ( i in s_ind)
{
  zval[k] =fx[i]
  k = k + 1
}
s_col = s_ind%/%n + 1 ; s_col[s_ind%/%n == n] = n
s_row = s_ind%%n + 1 ; s_row[s_row == 0] = n
par_mu = zval*4000 
df_mu = data.frame(x1=x1[s_col],x2=x2[s_row],par_mu=par_mu)

set.seed(2)
par_scale = rtruncnorm(ns, a=35, b=45, mean=40, sd=2)
par_shape = runif(ns,0.1,0.25)

xlist = list()
for (i in 1:ns){
  xlist[[i]] = rgev(nobs,loc= m0+par_mu[i],scale=par_scale[i],shape=par_shape[i])
}
sum(unlist(lapply(1:length(xlist),function(x) sum(is.na(xlist[[x]])))))

# ns_sample = sample(1:ns,4)
# par(mfrow=c(2,2))
# plot(m0+par_mu[ns_sample[1]]);plot(m0+par_mu[ns_sample[2]])
# plot(m0+par_mu[ns_sample[3]]);plot(m0+par_mu[ns_sample[4]])
# par(mfrow=c(1,1))

############################################################
# zlist - 2-D basis matrix by each location
# (fda package)
x_bsobj = create.bspline.basis(xyrange,norder=4, 
                               breaks=quantile(x1,prob = seq(0, 1, length = nBS)))
y_bsobj = create.bspline.basis(xyrange,norder=4, 
                               breaks=quantile(x2,prob = seq(0, 1, length = nBS)))
zlist = list()
for (i in 1:ns){
  xbs = eval.basis(df_mu$x1[i],x_bsobj)
  ybs = eval.basis(df_mu$x2[i],y_bsobj)
  tensorbs = do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs)) 
  zlist[[i]] = tensorbs 
}

# Omega matrix
Fmat = kronecker(bsplinepen(x_bsobj,Lfdobj=2),bsplinepen(y_bsobj,Lfdobj=0))
Gmat = kronecker(bsplinepen(x_bsobj,Lfdobj=0),bsplinepen(y_bsobj,Lfdobj=2))
Hmat = kronecker(bsplinepen(x_bsobj,Lfdobj=1),bsplinepen(y_bsobj,Lfdobj=1))
Om = Fmat+Gmat+2*Hmat

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
rho = 0.5

hygev_m = function(tvec, lam_s, xlist, zlist, Om, dmatrix, rho, z_init, u_init)  
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
  old_tvec = optim(tvec, hygev_m, lam_s=lam_s,
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

