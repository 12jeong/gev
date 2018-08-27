rm(list=ls())

library("Deriv")
library("evd")

n = 100
p = 2
true.vec = c(120,40,-0.1)
true.beta = c(5,0)
mu=true.vec[1]
sigma=true.vec[2]
k=true.vec[3]

set.seed(141)
z = matrix(rnorm(p*n),n,p)
x=rgev(n,loc=true.vec[1] ,scale=true.vec[2],shape=true.vec[3]) + z%*%true.beta

logl=expression(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k))

####### Stationary Sum #######
logll=expression(sum(log(sigma)+(1+1/k)*log(1+k*(x-mu)/sigma)+(1+k*(x-mu)/sigma)^(-1/k)))

Jaco=Deriv(logll,c("mu","sigma","k")) ; eval(Jaco)
Hmat=Deriv(logll,c("mu","sigma","k"),n=c(hessian=2)) ;eval(Hmat)

######### Stationary #########
Jaco1=Deriv(logl,c("mu","sigma","k"),combine="cbind") ; Jaco1
Hmat1=Deriv(logl,c("mu","sigma","k"),n=c(hessian=2),combine="cbind") ; Hmat1

head(eval(Jaco1))
head(eval(Hmat1))

jvec1=apply(eval(Jaco1),2,sum) ;jvec1
hmat1=matrix(apply(eval(Hmat1),2,sum),3,3) ;hmat1

######## NonStationary ########
jvec2=apply(cbind(eval(Jaco1),eval(Jaco1)[,1]*z),2,sum) ;jvec2
h1=hmat1
h2=rbind(apply(hmat1[1,1]*z,2,sum), apply(hmat1[1,2]*z,2,sum), apply(hmat1[1,3]*z,2,sum))
h3=t(h2)
h4=hmat1[1,1]*t(z)%*%z
hmat2=rbind(cbind(h1,h2),cbind(h3,h4)) ;hmat2

