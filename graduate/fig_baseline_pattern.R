# To verify that the initial value (fgev) is working

rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")

# surface base setting
mean_vec = c(0,0) 
sig_mat = matrix(c(30,0,0,30),nrow=2)
set_uni = dmvnorm(cbind(x1,x2), mean=mean_vec, sigma=sig_mat)
mean_vec1 = c(5,0); mean_vec2 = c(-5,0) 
sig_mat = matrix(c(10,0,0,10),nrow=2)
set_bi = 0.4*dmvnorm(cbind(x1,x2),mean=mean_vec1, sigma=sig_mat*1) +0.6*dmvnorm(cbind(x1,x2),mean=mean_vec2, sigma=sig_mat*2)

xyrange = c(-10,10)
# x = seq(xyrange[1],xyrange[2],length=20)
# y = seq(xyrange[1],xyrange[2],length=20)
x = seq(xyrange[1],xyrange[2],length=1000)
y = seq(xyrange[1],xyrange[2],length=1000)


par(mfrow=c(1,1))
fxy =  outer(x, y, function(x,y) -3*x+3*y) 
range(100+fxy)
range(40+0.5*fxy)
range(0.1+0.005*fxy)
pmat1 <- persp(x,y, fxy, phi = 10, theta = 30,
               xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
               xlab="x1",ylab="x2",zlab="",main="Plane")

fxy =  outer(x, y, function(x,y) dmvnorm(cbind(x,y), mean=mean_vec, sigma=sig_mat) )
pmat2 <- persp(x,y, fxy, phi = 10, theta = 30,
               xlim=c(-10,10), ylim=c(-10,10),ticktype="detailed",
               xlab="x1",ylab="x2",zlab="",main="Unimodal")
range(90+4000*fxy)
range(30+4000*fxy)
range(50*fxy)
fxy =  outer(x, y, function(x,y) 0.4*dmvnorm(cbind(x,y),mean=mean_vec1, sigma=sig_mat*1) +0.6*dmvnorm(cbind(x,y),mean=mean_vec2, sigma=sig_mat*2) )
pmat3 <- persp(x,y, fxy, phi = 10, theta = 30,
               xlim=c(-10,10), ylim=c(-10,10),ticktype="detailed",
               xlab="x1",ylab="x2",zlab="",main="Bimodal")
range(90+3000*fxy)
range(30+3000*fxy)
range(50*fxy)
