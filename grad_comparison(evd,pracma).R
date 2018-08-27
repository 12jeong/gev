#detach("package:fExtremes")
#search()
library(evd)
library(pracma)
library(Deriv)

mu=0
sigma=1
k=0.3
theta=c(mu,sigma,k)

set.seed(1882)
x=rgev(n=1000,loc=0,scale=1,shape=0.3)

## exact form
# substitution
y=(x-mu)/sigma
z=1+k*y
phi=z^(-1/k)
psi=k+1-phi
# partial derivariate
pd_ym=-1/sigma
pd_ys=-y/sigma
pd_yk=0
pd_zm=-k/sigma
pd_zs=-k*y/sigma
pd_zk=y
pd_phm= phi/(sigma*z)
pd_phs= y*phi/(sigma*z)
pd_phk= -phi*(log(phi)+y/z)/k
pd_psm= -pd_phm
pd_pss= -pd_phs
pd_psk= 1-pd_phk

grad_mu =sum(psi/(sigma*z))
grad_sigma= -sum(1-y*psi/z)/sigma
grad_k= -sum(y*psi/z+(1-phi)*log(phi))/k

?grad
gev_loglh <- function(theta){
  mu=theta[1]; sigma=theta[2]; k=theta[3]
  return(-sum(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k)))
}

f=function(mu,sigma,k){-sum(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k))}

# Gradient Comparison
grad(f=gev_loglh,theta)  # numerical (central differential) - pracma packages
c(grad_mu,grad_sigma,grad_k)  # exact form by hand
#Deriv(expression(f(mu,sigma,k)))  #  symbolic differentiation - Deriv packages
eval(Deriv(expression(f(mu,sigma,k))))

# Hessian
h11=sum((psi*k-phi)/z^2)/sigma^2 # mumu
h22=sum(1-2*psi*y/z-y^2*(phi-psi*k)/z^2)/sigma^2 # sigmasigma
h33=sum(2*(1-phi-phi*log(phi))*y/z  + 2*(1-phi)*log(phi)-phi*log(phi)^2 + y^2 *(psi*k-phi)/z^2) /k^2 #kk
h12=sum((psi*k*y-psi*z-phi*y)/(sigma*z)^2) # musigma
h13=sum((k*(z-psi*y)+phi*(z*log(phi)+y))/(sigma*k*z^2) ) # muk
h23=sum(y*(y*(phi-k*psi)+z*(k+phi*log(phi)))/(k*sigma*z^2)) # sigmak

# Comparison
hessian(f=gev_loglh,theta) # numerical (central)
matrix(c(h11,h12,h13,0,h22,h23,0,0,h33),3,3) # exact
matrix(eval(Deriv(expression(f(mu,sigma,k)),n=c(hessian=2))),3,3) # symbolic

