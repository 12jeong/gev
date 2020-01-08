rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")
library(xtable)
# colour.plot = c("royalblue1","forestgreen","orange") 
rainbow.n = rainbow(8, s = 1, v = 1, start = 0, end = max(1,8 - 1)/8, alpha = 0.7)
# colour.plot = rainbow(8, s = 1, v = 1, start = 0, end = max(1,8 - 1)/8, alpha = 0.7)[c(1,2,3)]
colour.plot  = c("grey95","grey70","grey30")

par(mfrow=c(1,1))
S_num = 26
eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/trainloss_',S_num, '.RData',"')"))))
# eval(parse(text = paste0("load(file =","'", paste0('C:/Users/UOS/Downloads/trainloss_',S_num, '.RData',"')"))))

xtable(data.frame(apply(hdist_mat,2,mean)[c(1,14,27,28)],apply(hdist_mat,2,median)[c(1,14,27,28)]))
xtable(data.frame(apply(rmse_mat,2,mean)[c(1,14,27,28)],apply(rmse_mat,2,median)[c(1,14,27,28)]))



### train loss
# dev.off()

# png(filename=paste0("./graduate/fig/fig_contrast/HD_train",S_num,".png"),
#     width=800,height=500)
par(cex=1.5)
plot.new()
l <- legend(0,0,legend=c(expression(lambda[mu]^"min"),expression(lambda[mu]^"max"),expression(lambda[paste(mu,"*")])),plot=F,
           col=colour.plot,pch=rep(20,3),cex=0.9,bty="n")
w <- grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc')
par(omd=c(0, 1-w, 0, 1))
boxplots.triple = boxplot(hdist_mat,outline=F,ylab="HD",xaxt="n",col=c(rep(colour.plot,9),"white"),xlim=c(1,28))
rect(xleft = 27.6, xright = 28.4, 
     ybottom = boxplots.triple$stats[4,28],
     ytop = boxplots.triple$stats[2,28],
     density=15, angle=45)
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:3)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c(expression(lambda[sigma]^"min"),expression(lambda[sigma]^"max"),expression(lambda[sigma]^"*")))
axis(1, at=c(1:9)*3-1, labels=rep(c(expression(lambda[kappa]^"min"),expression(lambda[kappa]^"max"),expression(lambda[kappa]^"*")),3))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
       c(expression(lambda[mu^"min"]),expression(lambda[mu^"max"]),expression(lambda[paste(mu,"*")]),expression(lambda[AIC])),
       border=rep("black",4),fill=c("grey95","grey70","grey30","white"),
       col=c(NA,NA,NA,"black"),
       cex=0.9,angle=c(NA,NA,NA,45),density=c(NA,NA,NA,20))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
       c("","","",""),
       col=c(colour.plot,"black"),cex=0.9,angle=c(0,0,0,45),density=c(0,0,0,20))
# dev.off()

# boxplot(rmse_mat,outline=F,ylab="RMSE",xaxt="n",col=c(rep(colour.plot,9),rainbow.n[6]),xlim=c(1,28))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# abline(v=c(1:3)*9+0.5,lwd=1)
# axis(3, at=c(5,14,23), labels=c(expression(lambda[sigma]^"min"),expression(lambda[sigma]^"max"),expression(lambda[sigma]^"*")))
# axis(1, at=c(1:9)*3-1, labels=rep(c(expression(lambda[kappa]^"min"),expression(lambda[kappa]^"max"),expression(lambda[kappa]^"*")),3))
# legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
#        c(expression(lambda[mu^"min"]),expression(lambda[mu^"max"]),expression(lambda[paste(mu,"*")]),expression(lambda[AIC])),
#        col=c(colour.plot,rainbow.n[6]),pch=rep(20,3),cex=0.9)


### test loss 
eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/testloss_',S_num, '.RData',"')"))))

# eval(parse(text = paste0("load(file =","'", paste0('C:/Users/UOS/Downloads/testloss_',S_num, '.RData',"')"))))
xtable(data.frame(apply(hdist_mat,2,mean)[c(1,14,27,28)],apply(hdist_mat,2,median)[c(1,14,27,28)]))
xtable(data.frame(apply(rmse_mat,2,mean)[c(1,14,27,28)],apply(rmse_mat,2,median)[c(1,14,27,28)]))



apply(hdist_mat,2,mean)[c(1,14,27,28)]
apply(rmse_mat,2,mean)[c(1,14,27,28)]

apply(hdist_mat,2,median)[c(1,14,27,28)]
apply(rmse_mat,2,median)[c(1,14,27,28)]

# png(filename=paste0("./graduate/fig/fig_contrast/HD_test",S_num,".png"),width=800,height=500)
# par(cex=1.5)
boxplots.triple = boxplot(log(hdist_mat),outline=F,ylab="HD (log scale)",xaxt="n",col=c(rep(colour.plot,9),"white"),xlim=c(1,28))
rect(xleft = 27.6, xright = 28.4, 
     ybottom = boxplots.triple$stats[4,28],
     ytop = boxplots.triple$stats[2,28],
     density=15, angle=45)
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:3)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c(expression(lambda[sigma]^"min"),expression(lambda[sigma]^"max"),expression(lambda[sigma]^"*")))
axis(1, at=c(1:9)*3-1, labels=rep(c(expression(lambda[kappa]^"min"),expression(lambda[kappa]^"max"),expression(lambda[kappa]^"*")),3))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
       c(expression(lambda[mu^"min"]),expression(lambda[mu^"max"]),expression(lambda[paste(mu,"*")]),expression(lambda[AIC])),
       border=rep("black",4),fill=c("grey95","grey70","grey30","white"),
       col=c(NA,NA,NA,"black"),
       cex=0.9,angle=c(NA,NA,NA,45),density=c(NA,NA,NA,20))
legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
       c("","","",""),
       col=c(colour.plot,"black"),cex=0.9,angle=c(0,0,0,45),density=c(0,0,0,20))
# dev.off()
apply(hdist_mat,2,mean)

# boxplot(log(rmse_mat),outline=F,ylab="RMSE (log scale)",xaxt="n",col=c(rep(colour.plot,9),rainbow.n[6]),xlim=c(1,28))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# abline(v=c(1:3)*9+0.5,lwd=1)
# axis(3, at=c(5,14,23), labels=c(expression(lambda[sigma]^"min"),expression(lambda[sigma]^"max"),expression(lambda[sigma]^"*")))
# axis(1, at=c(1:9)*3-1, labels=rep(c(expression(lambda[kappa]^"min"),expression(lambda[kappa]^"max"),expression(lambda[kappa]^"*")),3))
# legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA,
#        c(expression(lambda[mu^"min"]),expression(lambda[mu^"max"]),expression(lambda[paste(mu,"*")]),expression(lambda[AIC])),
#        col=c(colour.plot,rainbow.n[6]),pch=rep(20,3),cex=0.9)
# 
# boxplot(hdist_mat[,c(13:18,22:28)],outline=F,ylab="HDllinger Distance",xaxt="n",col=rep(colour.plot,3))
# abline(v=c(1:4)*3+0.5,lty=2,col="gray")
# boxplot(rmse_mat[,c(13:18,22:28)],outline=F,ylab="RMSE",xaxt="n",col=rep(colour.plot,3))
# abline(v=c(1:4)*3+0.5,lty=2,col="gray")

apply(hdist_mat,2,mean)[c(1,14,27,28)]
apply(rmse_mat,2,mean)[c(1,14,27,28)]

apply(hdist_mat,2,median)[c(1,14,27,28)]
apply(rmse_mat,2,median)[c(1,14,27,28)]
