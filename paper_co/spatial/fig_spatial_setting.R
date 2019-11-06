rm(list=ls())
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
fxy =  outer(x, y, function(x,y) {par_mu = 100+(-2*x-3*y)})
range(fxy)
set.seed(4)
x1 = runif(ns,xyrange[1],xyrange[2])
x2 = runif(ns,xyrange[1],xyrange[2])
par_mu = 100+(-2*x1-3*x2) # 평면
df_mu = data.frame(x1=x1,x2=x2,par_mu=par_mu)

## persp 
# pmat1 <- persp(x,y, fxy, phi = 10, theta = 120, zlim=c(20,max(fxy)), 
#                xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
#                xlab="x1",ylab="x2",zlab="")
# mypoints1 <- trans3d(df_mu$x1, df_mu$x2, rep(20,ns), pmat=pmat1)
# points(mypoints1, pch=16, col=1)

pmat2 <- persp(x,y, fxy, phi = 10, theta = 200,
               xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
               xlab="x1",ylab="x2",zlab="")
mypoints2 <- trans3d(df_mu$x1, df_mu$x2, df_mu$par_mu, pmat=pmat2)
points(mypoints2, pch=16, col=1)
# plot_ly
# plot_ly(x=x,y=y,z=fxy, type="surface") %>%
#   add_trace(data=df_mu, x=x1, y=x2, z=rep(40,ns), mode="markers",
#             type="scatter3d",
#             marker = list(size = 5, color = "black", symbol = 104, opacity=0.6)) %>%
#   layout(scene = list(
#     xaxis = list(title = "x1"),
#     yaxis = list(title = "x2"),
#     zaxis = list(title = "z" )))



### setting2 : 일봉분포 ###
mu = c(0,0)
sig = matrix(c(30,0,0,30),nrow=2)
fxy = outer(x, y, function(x,y) {
  dmvnorm(cbind(x,y), mean=mu, sigma=sig)
})
x1 = seq(xyrange[1],xyrange[2],length=n)
x2 = seq(xyrange[1],xyrange[2],length=n)
fx = outer(x1, x2, function(x1,x2) {
  dmvnorm(cbind(x1,x2), mean=mu, sigma=sig)
})
set.seed(3)
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
par_mu = 90+zval*4000 # setting2_2 : 91~110
df_mu = data.frame(x1=x1[s_col],x2=x2[s_row],par_mu=par_mu)
## persp 
# range(90+4000*fxy)
# pmat3 <- persp(x,y, 90+4000*fxy, phi = 10, theta = 120, zlim=c(80,max(90+4000*fxy)),
#                xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
#                xlab="x1",ylab="x2",zlab="")
# mypoints3 <- trans3d(df_mu$x2, df_mu$x1, rep(80,ns), pmat=pmat3)
# points(mypoints3, pch=16, col=1)

pmat4 <- persp(x,y, 90+4000*fxy, phi = 10, theta = 20,
               xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
               xlab="x1",ylab="x2",zlab="",border="gray")
mypoints4 <- trans3d(df_mu$x2, df_mu$x1, df_mu$par_mu, pmat=pmat4)
points(mypoints4, pch=16, col=1)

### setting3 : 이봉분포 ###
mu1 = c(5,0)
mu2 = c(-5,0)
sig = matrix(c(10,0,0,10),nrow=2)
fxy = outer(x, y, function(x,y) {
  0.4*dmvnorm(cbind(x,y),mean=mu1, sigma=sig*1) +
    0.6*dmvnorm(cbind(x,y),mean=mu2, sigma=sig*2)
})
x1 = seq(xyrange[1],xyrange[2],length=n)
x2 = seq(xyrange[1],xyrange[2],length=n)
fx = outer(x1, x2, function(x1,x2) {
  0.4*dmvnorm(cbind(x1,x2),mean=mu1, sigma=sig*1) +
    0.6*dmvnorm(cbind(x1,x2),mean=mu2, sigma=sig*2)
})
set.seed(3)
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
par_mu = 90+zval*3000 
df_mu = data.frame(x1=x1[s_col],x2=x2[s_row],par_mu=par_mu)
## persp 
range(90+fx*3000)
# pmat5 <- persp(x,y, 95+fxy*2000, phi = 5, theta = 20, zlim=c(90.5,max(95+fxy*2000)),
#                xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
#                xlab="x1",ylab="x2",zlab="")
# mypoints5 <- trans3d(df_mu$x2, df_mu$x1, rep(90.5,ns), pmat=pmat5)
# points(mypoints5, pch=16, col=1)

pmat6 <- persp(x,y, 90+fxy*3000, phi = 10, theta = 20,  
               xlim=c(-12,12), ylim=c(-12,12),ticktype="detailed",
               xlab="x1",ylab="x2",zlab="",border="gray")
mypoints6 <- trans3d(df_mu$x2, df_mu$x1, df_mu$par_mu, pmat=pmat6)
points(mypoints6, pch=16, col=1)


