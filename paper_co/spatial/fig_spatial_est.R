rm(list=ls()); gc()
setwd("~/GitHub/gev")
source("./lib/pack.R")
source("./lib/sgevlibrary.R")

par(mfrow=c(1,3))

set.seed(1)
xyrange = c(-10,10)
n = 100
ns = 50
x = seq(xyrange[1],xyrange[2],length=20)
y = seq(xyrange[1],xyrange[2],length=20)


### setting1 : 평면 ###
load("./Rexport/Spatial_simulation/result1/numerical_setting1.RDa")
load("./Rexport/Spatial_simulation/result1/result1_seed1.RDa")
lambdaset
fxy =  outer(x, y, function(x,y) {par_mu = 100+(-2*x-3*y)})

fit = result1[[12]]
Z = do.call('rbind',zlist)
p = length(zlist[[1]])
# plot3d(df_mu$x1,df_mu$x2,fit[1]+drop(Z%*%tail(fit,p)),
#        xlim=c(-12,12),ylim=c(-12,12),ticktype="detailed",
#        xlab="x1",ylab="x2",zlab="")
pmat1 <- persp(x,y, fxy, phi = 10, theta = 200,
               xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
               xlab="x1",ylab="x2",zlab="",border="gray")
mypoints1 <- trans3d(df_mu$x1, df_mu$x2, fit[1]+drop(Z%*%tail(fit,p)), pmat=pmat1)
points(mypoints1, pch=16, col=1)


### setting2 : 일봉분포 ###
load("./Rexport/Spatial_simulation/result2/numerical_setting2.RDa")
load("./Rexport/Spatial_simulation/result2/result2_seed1.RDa")
lambdaset
mu = c(0,0)
sig = matrix(c(30,0,0,30),nrow=2)
fxy = outer(x, y, function(x,y) {
  dmvnorm(cbind(x,y), mean=mu, sigma=sig)
})
fit = result1[[3]]
Z = do.call('rbind',zlist)
p = length(zlist[[1]])
# plot3d(df_mu$x1,df_mu$x2,fit[1]+drop(Z%*%tail(fit,p)),
#        xlim=c(-12,12),ylim=c(-12,12),ticktype="detailed",
#        xlab="x1",ylab="x2",zlab="")
pmat2 <- persp(x,y, 90+4000*fxy, phi = 10, theta = 200,
               xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
               xlab="x1",ylab="x2",zlab="",border="gray")
mypoints2 <- trans3d(df_mu$x1, df_mu$x2, fit[1]+drop(Z%*%tail(fit,p)), pmat=pmat2)
points(mypoints2, pch=16, col=1)


### setting3 : 이봉분포 ###
load("./Rexport/Spatial_simulation/result3/numerical_setting3.RDa")
load("./Rexport/Spatial_simulation/result3/result3_seed1.RDa")
mu1 = c(5,0)
mu2 = c(-5,0)
sig = matrix(c(10,0,0,10),nrow=2)
fxy = outer(x, y, function(x,y) {
  0.4*dmvnorm(cbind(x,y),mean=mu1, sigma=sig*1) +
    0.6*dmvnorm(cbind(x,y),mean=mu2, sigma=sig*2)
})
lambdaset
fit = result1[[3]]
Z = do.call('rbind',zlist)
p = length(zlist[[1]])
# plot3d(df_mu$x1,df_mu$x2,fit[1]+drop(Z%*%tail(fit,p)),
#        xlim=c(-12,12),ylim=c(-12,12),ticktype="detailed",
#        xlab="x1",ylab="x2",zlab="")
pmat3 <- persp(x,y, 90+fxy*3000, phi = 10, theta = 20,  
               xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
               xlab="x1",ylab="x2",zlab="",border="gray")
mypoints3 <- trans3d(df_mu$x2, df_mu$x1, fit[1]+drop(Z%*%tail(fit,p)), pmat=pmat3)
points(mypoints3, pch=16, col=1)


# plot(df_mu$x1,fit[1]+drop(Z%*%tail(fit,p)))
# plot(df_mu$x2,fit[1]+drop(Z%*%tail(fit,p)))
