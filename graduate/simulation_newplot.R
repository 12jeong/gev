rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")

S_num = 27

eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/AIC_scenario',S_num, '.RData',"')"))))
eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/newloss_scenario',S_num, '.RData',"')"))))

min.ind =  unlist(lapply(AIC_list,function(x) which.min(x)))
lam.min = as.numeric(names(which.max(table(min.ind))))
lam.min.vec = lam.grid2[lam.min,]

loc.map = c(0,max(lam_set1),lam.min.vec[1])
sc.map = c(0,max(lam_set2),lam.min.vec[2])
sh.map = c(0,max(lam_set3),lam.min.vec[3])

hdist.test_mat = do.call("rbind",hdist.test_list)



par(mfrow=c(2,1))
############### Hellinger Distance
#### 1
match.map = expand.grid(loc=loc.map,sh=sh.map,sc=sc.map)
match.vec = c()
for (i in 1:nrow(match.map)){
  match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
                                     function(x) x[c(1,3,2)]== match.map[i,]))==3)
}
df.hdist.test = data.frame(hdist.test_mat[,match.vec])
which.min(apply(df.hdist.test,2,mean))
t.test(df.hdist.test[,14],df.hdist.test[,27])
apply(df.hdist.test,2,sd)

boxplot(df.hdist.test,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_2(min)","lam_2(max)","lam_2(*)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_3(min)","lam_3(max)","lam_3(*)"),3))
legend("topright",legend=c("lam_1(min)","lam_1(max)","lam_1(*)"),col=c(2:4),pch=rep(20,3))
#### 2
match.map = expand.grid(loc=loc.map,sc=sc.map,sh=sh.map)
match.vec = c()
for (i in 1:nrow(match.map)){
  match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
                                     function(x) x[c(1,2,3)]== match.map[i,]))==3)
}
df.hdist.test = data.frame(hdist.test_mat[,match.vec])
boxplot(df.hdist.test,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_3(min)","lam_3(max)","lam_3(*)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_2(min)","lam_2(max)","lam_2(*)"),3))
legend("topright",legend=c("lam_1(min)","lam_1(max)","lam_1(*)"),col=c(2:4),pch=rep(20,3))
legend("topright",legend=c("lam_2(min)","lam_2(max)","lam_2(*)"),col=c(2:4),pch=rep(20,3))
#### 3
match.map = expand.grid(sc=sc.map,loc=loc.map,sh=sh.map)
match.vec = c()
for (i in 1:nrow(match.map)){
  match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
                                     function(x) x[c(2,1,3)]== match.map[i,]))==3)
}
df.hdist.test = data.frame(hdist.test_mat[,match.vec])
boxplot(df.hdist.test,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_3(min)","lam_3(max)","lam_3(*)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_1(min)","lam_1(max)","lam_1(*)"),3))
legend("topright",legend=c("lam_2(min)","lam_2(max)","lam_2(*)"),col=c(2:4),pch=rep(20,3))
#### 4
match.map = expand.grid(sc=sc.map,sh=sh.map,loc=loc.map)
match.vec = c()
for (i in 1:nrow(match.map)){
  match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
                                     function(x) x[c(2,3,1)]== match.map[i,]))==3)
}
df.hdist.test = data.frame(hdist.test_mat[,match.vec])
boxplot(df.hdist.test,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_1(min)","lam_1(max)","lam_1(*)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_3(min)","lam_3(max)","lam_3(*)"),3))
legend("topright",legend=c("lam_2(min)","lam_2(max)","lam_2(*)"),col=c(2:4),pch=rep(20,3))

#### 5
match.map = expand.grid(sh=sh.map,sc=sc.map,loc=loc.map)
match.vec = c()
for (i in 1:nrow(match.map)){
  match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
                                     function(x) x[c(3,2,1)]== match.map[i,]))==3)
}
df.hdist.test = data.frame(hdist.test_mat[,match.vec])
boxplot(df.hdist.test,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_1(min)","lam_1(max)","lam_1(*)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_2(min)","lam_2(max)","lam_2(*)"),3))
legend("topright",legend=c("lam_3(min)","lam_3(max)","lam_3(*)"),col=c(2:4),pch=rep(20,3))
#### 6
match.map = expand.grid(sh=sh.map,loc=loc.map,sc=sc.map)
match.vec = c()
for (i in 1:nrow(match.map)){
  match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
                                     function(x) x[c(3,1,2)]== match.map[i,]))==3)
}
df.hdist.test = data.frame(hdist.test_mat[,match.vec])
boxplot(df.hdist.test,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_2(min)","lam_2(max)","lam_2(*)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_1(min)","lam_1(max)","lam_1(*)"),3))
legend("topright",legend=c("lam_3(min)","lam_3(max)","lam_3(*)"),col=c(2:4),pch=rep(20,3))