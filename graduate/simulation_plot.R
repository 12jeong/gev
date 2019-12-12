rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")

S_num = 2

eval(parse(text = paste0("load(file =","'", paste0('./Rexport/RData_sgev3_simulation/testloss_scenario',S_num, '.RData',"')"))))
lam.grid2 = expand.grid(lam_set1,lam_set2,lam_set3)
rmse_mean = apply(rmse_mat,2,mean)

min.ind =  unlist(lapply(AIC_list,function(x) which.min(x)))
lam.min = as.numeric(names(which.max(table(min.ind))))
# lam.grid2[as.numeric(names(table(min.ind))),]
lam.min.vec = lam.grid2[lam.min,]
lam.min.vec

loc.map = c(0,max(lam_set1),lam.min.vec[1])
sc.map = c(0,max(lam_set2),lam.min.vec[2])
sh.map = c(0,max(lam_set3),lam.min.vec[3])


# ######################### lambda_* ######################
# par(mfrow=c(2,1))
# ################### RMSE
# #### 1
# match.map = expand.grid(loc=loc.map,sh=sh.map,sc=sc.map)
# match.vec = c()
# for (i in 1:nrow(match.map)){
#   match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
#                                      function(x) x[c(1,3,2)]== match.map[i,]))==3)
# }
# df.rmse = data.frame(rmse_mat[,match.vec])
# boxplot(df.rmse,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# abline(v=c(1:2)*9+0.5,lwd=1)
# axis(3, at=c(5,14,23), labels=c("lam_2(min)","lam_2(max)","lam_2(*)"))
# axis(1, at=c(1:9)*3-1, labels=rep(c("lam_3(min)","lam_3(max)","lam_3(*)"),3))
# legend("topright",legend=c("lam_1(min)","lam_1(max)","lam_1(*)"),col=c(2:4),pch=rep(20,3))
# #### 2
# match.map = expand.grid(loc=loc.map,sc=sc.map,sh=sh.map)
# match.vec = c()
# for (i in 1:nrow(match.map)){
#   match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
#                                      function(x) x[c(1,2,3)]== match.map[i,]))==3)
# }
# df.rmse = data.frame(rmse_mat[,match.vec])
# boxplot(df.rmse,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# abline(v=c(1:2)*9+0.5,lwd=1)
# axis(3, at=c(5,14,23), labels=c("lam_3(min)","lam_3(max)","lam_3(*)"))
# axis(1, at=c(1:9)*3-1, labels=rep(c("lam_2(min)","lam_2(max)","lam_2(*)"),3))
# legend("topright",legend=c("lam_1(min)","lam_1(max)","lam_1(*)"),col=c(2:4),pch=rep(20,3))
# #### 3
# match.map = expand.grid(sc=sc.map,loc=loc.map,sh=sh.map)
# match.vec = c()
# for (i in 1:nrow(match.map)){
#   match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
#                                      function(x) x[c(2,1,3)]== match.map[i,]))==3)
# }
# df.rmse = data.frame(rmse_mat[,match.vec])
# boxplot(df.rmse,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# abline(v=c(1:2)*9+0.5,lwd=1)
# axis(3, at=c(5,14,23), labels=c("lam_3(min)","lam_3(max)","lam_3(*)"))
# axis(1, at=c(1:9)*3-1, labels=rep(c("lam_1(min)","lam_1(max)","lam_1(*)"),3))
# legend("topright",legend=c("lam_2(min)","lam_2(max)","lam_2(*)"),col=c(2:4),pch=rep(20,3))
# #### 4
# match.map = expand.grid(sc=sc.map,sh=sh.map,loc=loc.map)
# match.vec = c()
# for (i in 1:nrow(match.map)){
#   match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
#                                      function(x) x[c(2,3,1)]== match.map[i,]))==3)
# }
# df.rmse = data.frame(rmse_mat[,match.vec])
# boxplot(df.rmse,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# abline(v=c(1:2)*9+0.5,lwd=1)
# axis(3, at=c(5,14,23), labels=c("lam_1(min)","lam_1(max)","lam_1(*)"))
# axis(1, at=c(1:9)*3-1, labels=rep(c("lam_3(min)","lam_3(max)","lam_3(*)"),3))
# legend("topright",legend=c("lam_2(min)","lam_2(max)","lam_2(*)"),col=c(2:4),pch=rep(20,3))
# 
# #### 5
# match.map = expand.grid(sh=sh.map,sc=sc.map,loc=loc.map)
# match.vec = c()
# for (i in 1:nrow(match.map)){
#   match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
#                                      function(x) x[c(3,2,1)]== match.map[i,]))==3)
# }
# df.rmse = data.frame(rmse_mat[,match.vec])
# boxplot(df.rmse,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# abline(v=c(1:2)*9+0.5,lwd=1)
# axis(3, at=c(5,14,23), labels=c("lam_1(min)","lam_1(max)","lam_1(*)"))
# axis(1, at=c(1:9)*3-1, labels=rep(c("lam_2(min)","lam_2(max)","lam_2(*)"),3))
# legend("topright",legend=c("lam_3(min)","lam_3(max)","lam_3(*)"),col=c(2:4),pch=rep(20,3))
# #### 6
# match.map = expand.grid(sh=sh.map,loc=loc.map,sc=sc.map)
# match.vec = c()
# for (i in 1:nrow(match.map)){
#   match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
#                                      function(x) x[c(3,1,2)]== match.map[i,]))==3)
# }
# df.rmse = data.frame(rmse_mat[,match.vec])
# boxplot(df.rmse,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# abline(v=c(1:2)*9+0.5,lwd=1)
# axis(3, at=c(5,14,23), labels=c("lam_2(min)","lam_2(max)","lam_2(*)"))
# axis(1, at=c(1:9)*3-1, labels=rep(c("lam_1(min)","lam_1(max)","lam_1(*)"),3))
# legend("topright",legend=c("lam_3(min)","lam_3(max)","lam_3(*)"),col=c(2:4),pch=rep(20,3))

############### Hellinger Distance
#### 1
match.map = expand.grid(loc=loc.map,sh=sh.map,sc=sc.map)
match.vec = c()
for (i in 1:nrow(match.map)){
  match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
                                     function(x) x[c(1,3,2)]== match.map[i,]))==3)
}
df.hdist = data.frame(hdist_mat[,match.vec])
boxplot(df.hdist,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
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
df.hdist = data.frame(hdist_mat[,match.vec])
boxplot(df.hdist,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_3(min)","lam_3(max)","lam_3(*)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_2(min)","lam_2(max)","lam_2(*)"),3))
legend("topright",legend=c("lam_1(min)","lam_1(max)","lam_1(*)"),col=c(2:4),pch=rep(20,3))

#### 3
match.map = expand.grid(sc=sc.map,loc=loc.map,sh=sh.map)
match.vec = c()
for (i in 1:nrow(match.map)){
  match.vec[i] = which(colSums(apply(expand.grid(lam_set1,lam_set2,lam_set3) ,1,
                                     function(x) x[c(2,1,3)]== match.map[i,]))==3)
}
df.hdist = data.frame(hdist_mat[,match.vec])
boxplot(df.hdist,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
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
df.hdist = data.frame(hdist_mat[,match.vec])
boxplot(df.hdist,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
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
df.hdist = data.frame(hdist_mat[,match.vec])
boxplot(df.hdist,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
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
df.hdist = data.frame(hdist_mat[,match.vec])
boxplot(df.hdist,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_2(min)","lam_2(max)","lam_2(*)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_1(min)","lam_1(max)","lam_1(*)"),3))
legend("topright",legend=c("lam_3(min)","lam_3(max)","lam_3(*)"),col=c(2:4),pch=rep(20,3))


# three row
# df.rmse = data.frame(rmse_mat[,match.vec])
# par(mfrow=c(3,1))
# boxplot(df.rmse[,c(1:9)],outline=F,ylab="RMSE",xlab="Model Number",names=rep("",9))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# boxplot(df.rmse[,c(1:9)+9],outline=F,ylab="RMSE",xlab="Model Number",names=rep("",9))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")
# boxplot(df.rmse[,c(1:9)+18],outline=F,ylab="RMSE",xlab="Model Number",names=rep("",9))
# abline(v=c(1:8)*3+0.5,lty=2,col="gray")



###################### lambda_AIC 추가

par(mfrow=c(1,2))
# a=c()
# for(i in 1:nrow(rmse_mat)) a[i] = rmse_mat[i,min.ind[i]]
# boxplot(df.rmse[,ncol(df.rmse)],a)
# b=c()
# for(i in 1:nrow(hdist_mat)) b[i] = hdist_mat[i,min.ind[i]]
# boxplot(df.hdist[,ncol(df.hdist)],b)

min.ind =  unlist(lapply(AIC_list,function(x) which.min(x)))
lam.min.mat = lam.grid2[min.ind,]

loc.map = c(0,max(lam_set1),NA)
sc.map = c(0,max(lam_set2),NA)
sh.map = c(0,max(lam_set3),NA)
# map_tmp = expand.grid(loc=loc.map, sc=sc.map, sh=sh.map)
map_tmp = expand.grid(loc=loc.map, sc=sc.map, sh=sh.map)

box_rmse_list = list()
for (i in 1:nrow(map_tmp)){
  rmse_tmp = c()
  vec_tmp = rep(NA,3)
    # i=27
    if (sum(is.na(map_tmp[i,]))==3){
      for(j in 1:nrow(rmse_mat)) rmse_tmp[j] = rmse_mat[j,min.ind[j]]
      box_rmse_list[[i]] = rmse_tmp
    }
    if (sum(is.na(map_tmp[i,]))==1){
      tmp_ind = which(!is.na(map_tmp[i,]))
      tmp_ind2 = which(is.na(map_tmp[i,]))
      vec_tmp[tmp_ind] = unlist(map_tmp[i,tmp_ind])
      for (k in 1:nrow(rmse_mat)){
        vec_tmp[tmp_ind2] = unlist(lam.grid2[min.ind,tmp_ind2][k])
        match.ind = which(colSums(apply(lam.grid2,1,function(x) x== vec_tmp)) ==3)
        rmse_tmp[k] = rmse_mat[k,match.ind]
      }
      box_rmse_list[[i]] = rmse_tmp
    }
    if (sum(is.na(map_tmp[i,]))==2){
      tmp_ind = which(!is.na(map_tmp[i,]))
      tmp_ind2 = which(is.na(map_tmp[i,]))
      vec_tmp[tmp_ind] = unlist(map_tmp[i,tmp_ind])
      for (k in 1:nrow(rmse_mat)){
        vec_tmp[tmp_ind2[1]] = unlist(lam.grid2[min.ind,tmp_ind2[1]][k])
        vec_tmp[tmp_ind2[2]] = unlist(lam.grid2[min.ind,tmp_ind2[2]][k])
        match.ind = which(colSums(apply(lam.grid2,1,function(x) x== vec_tmp)) ==3)
        rmse_tmp[k] = rmse_mat[k,match.ind]
      }
      box_rmse_list[[i]] = rmse_tmp
    }
  if (sum(is.na(map_tmp[i,]))==0){
    vec_tmp = map_tmp[i,]
    match.vec = which(colSums(apply(lam.grid2,1,function(x) x== vec_tmp))==3)
    box_rmse_list[[i]] = rmse_mat[,match.vec]
  }
}
box_rmse_tmp = do.call("cbind",box_rmse_list)
############ RMSE
par(mfrow=c(2,1))
#### 1
boxplot(box_rmse_tmp,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_3(min)","lam_3(max)","lam_3(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_2(min)","lam_2(max)","lam_2(aic)"),3))
legend("topright",legend=c("lam_1(min)","lam_1(max)","lam_1(aic)"),col=c(2:4),pch=rep(20,3))
# par(mfrow=c(1,3))
# boxplot(box_rmse_df[,c(1:9)],outline=F,ylab="RMSE",xlab="",main="lam3(min)",xaxt="n",ylim=c(range(box_rmse_df)))
# axis(1, at=c(2,5,8), labels=c("lam_2(min)","lam_2(max)","lam_2(*)"))
# abline(v=c(1:2)*3+0.5,lty=2,col="gray")
# boxplot(box_rmse_df[,c(1:9)+9],outline=F,ylab="RMSE",xlab="",main="lam3(max)",xaxt="n",ylim=c(range(box_rmse_df)))
# axis(1, at=c(2,5,8), labels=c("lam_2(min)","lam_2(max)","lam_2(*)"))
# abline(v=c(1:2)*3+0.5,lty=2,col="gray")
# boxplot(box_rmse_df[,c(1:9)+18],outline=F,ylab="RMSE",xlab="",main="lam3(*)",xaxt="n",ylim=c(range(box_rmse_df)))
# axis(1, at=c(2,5,8), labels=c("lam_2(min)","lam_2(max)","lam_2(*)"))
# abline(v=c(1:2)*3+0.5,lty=2,col="gray")

#### 2
box_rmse_df = data.frame(box_rmse_tmp[,c(1:3,10:12,19:21)],box_rmse_tmp[,c(1:3,10:12,19:21)+3],box_rmse_tmp[,c(1:3,10:12,19:21)+6])
boxplot(box_rmse_df,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_2(min)","lam_2(max)","lam_2(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_3(min)","lam_3(max)","lam_3(aic)"),3))
legend("topright",legend=c("lam_1(min)","lam_1(max)","lam_1(aic)"),col=c(2:4),pch=rep(20,3))

#### 3
box_rmse_df = data.frame(box_rmse_tmp[,c(1:9)*3-2],box_rmse_tmp[,c(1:9)*3-1],box_rmse_tmp[,c(1:9)*3])
boxplot(box_rmse_df,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_1(min)","lam_1(max)","lam_1(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_3(min)","lam_3(max)","lam_3(aic)"),3))
legend("topright",legend=c("lam_2(min)","lam_2(max)","lam_2(aic)"),col=c(2:4),pch=rep(20,3))

#### 4
box_rmse_df = data.frame(box_rmse_tmp[,c(1,4,7,2,5,8,3,6,9)],box_rmse_tmp[,c(1,4,7,2,5,8,3,6,9)+9],box_rmse_tmp[,c(1,4,7,2,5,8,3,6,9)+18])
boxplot(box_rmse_df,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_3(min)","lam_3(max)","lam_3(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_1(min)","lam_1(max)","lam_1(aic)"),3))
legend("topright",legend=c("lam_2(min)","lam_2(max)","lam_2(aic)"),col=c(2:4),pch=rep(20,3))

#### 5
box_rmse_df = data.frame(box_rmse_tmp[,c(1,10,19,2,11,20,3,12,21)],
                         box_rmse_tmp[,c(1,10,19,2,11,20,3,12,21)+3],
                         box_rmse_tmp[,c(1,10,19,2,11,20,3,12,21)+6])
boxplot(box_rmse_df,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_2(min)","lam_2(max)","lam_2(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_1(min)","lam_1(max)","lam_1(aic)"),3))
legend("topright",legend=c("lam_3(min)","lam_3(max)","lam_3(aic)"),col=c(2:4),pch=rep(20,3))

#### 6
box_rmse_df = data.frame(box_rmse_tmp[,c(1,10,19,4,13,22,7,16,25)],
                         box_rmse_tmp[,c(1,10,19,4,13,22,7,16,25)+1],
                         box_rmse_tmp[,c(1,10,19,4,13,22,7,16,25)+2])
boxplot(box_rmse_df,outline=F,ylab="RMSE",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_1(min)","lam_1(max)","lam_1(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_2(min)","lam_2(max)","lam_2(aic)"),3))
legend("topright",legend=c("lam_3(min)","lam_3(max)","lam_3(aic)"),col=c(2:4),pch=rep(20,3))




# lambda_AIC 그림 - Hellinger Distance
box_hdist_list = list()
for (i in 1:nrow(map_tmp)){
  hdist_tmp = c()
  vec_tmp = rep(NA,3)
  # i=27
  if (sum(is.na(map_tmp[i,]))==3){
    for(j in 1:nrow(hdist_mat)) hdist_tmp[j] = hdist_mat[j,min.ind[j]]
    box_hdist_list[[i]] = hdist_tmp
  }
  if (sum(is.na(map_tmp[i,]))==1){
    tmp_ind = which(!is.na(map_tmp[i,]))
    tmp_ind2 = which(is.na(map_tmp[i,]))
    vec_tmp[tmp_ind] = unlist(map_tmp[i,tmp_ind])
    for (k in 1:nrow(hdist_mat)){
      vec_tmp[tmp_ind2] = unlist(lam.grid2[min.ind,tmp_ind2][k])
      match.ind = which(colSums(apply(lam.grid2,1,function(x) x== vec_tmp)) ==3)
      hdist_tmp[k] = hdist_mat[k,match.ind]
    }
    box_hdist_list[[i]] = hdist_tmp
  }
  if (sum(is.na(map_tmp[i,]))==2){
    tmp_ind = which(!is.na(map_tmp[i,]))
    tmp_ind2 = which(is.na(map_tmp[i,]))
    vec_tmp[tmp_ind] = unlist(map_tmp[i,tmp_ind])
    for (k in 1:nrow(hdist_mat)){
      vec_tmp[tmp_ind2[1]] = unlist(lam.grid2[min.ind,tmp_ind2[1]][k])
      vec_tmp[tmp_ind2[2]] = unlist(lam.grid2[min.ind,tmp_ind2[2]][k])
      match.ind = which(colSums(apply(lam.grid2,1,function(x) x== vec_tmp)) ==3)
      hdist_tmp[k] = hdist_mat[k,match.ind]
    }
    box_hdist_list[[i]] = hdist_tmp
  }
  if (sum(is.na(map_tmp[i,]))==0){
    vec_tmp = map_tmp[i,]
    match.vec = which(colSums(apply(lam.grid2,1,function(x) x== vec_tmp))==3)
    box_hdist_list[[i]] = hdist_mat[,match.vec]
  }
}
match.map
box_hdist_tmp = do.call("cbind",box_hdist_list)
############ hdist
par(mfrow=c(2,1))
#### 1
boxplot(box_hdist_tmp,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_3(min)","lam_3(max)","lam_3(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_2(min)","lam_2(max)","lam_2(aic)"),3))
legend("topright",legend=c("lam_1(min)","lam_1(max)","lam_1(aic)"),col=c(2:4),pch=rep(20,3))

#### 2
box_hdist_df = data.frame(box_hdist_tmp[,c(1:3,10:12,19:21)],box_hdist_tmp[,c(1:3,10:12,19:21)+3],box_hdist_tmp[,c(1:3,10:12,19:21)+6])
boxplot(box_hdist_df,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_2(min)","lam_2(max)","lam_2(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_3(min)","lam_3(max)","lam_3(aic)"),3))
legend("topright",legend=c("lam_1(min)","lam_1(max)","lam_1(aic)"),col=c(2:4),pch=rep(20,3))

#### 3
box_hdist_df = data.frame(box_hdist_tmp[,c(1:9)*3-2],box_hdist_tmp[,c(1:9)*3-1],box_hdist_tmp[,c(1:9)*3])
boxplot(box_hdist_df,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_1(min)","lam_1(max)","lam_1(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_3(min)","lam_3(max)","lam_3(aic)"),3))
legend("topright",legend=c("lam_2(min)","lam_2(max)","lam_2(aic)"),col=c(2:4),pch=rep(20,3))

#### 4
box_hdist_df = data.frame(box_hdist_tmp[,c(1,4,7,2,5,8,3,6,9)],box_hdist_tmp[,c(1,4,7,2,5,8,3,6,9)+9],box_hdist_tmp[,c(1,4,7,2,5,8,3,6,9)+18])
boxplot(box_hdist_df,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_3(min)","lam_3(max)","lam_3(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_1(min)","lam_1(max)","lam_1(aic)"),3))
legend("topright",legend=c("lam_2(min)","lam_2(max)","lam_2(aic)"),col=c(2:4),pch=rep(20,3))

#### 5
box_hdist_df = data.frame(box_hdist_tmp[,c(1,10,19,2,11,20,3,12,21)],
                          box_hdist_tmp[,c(1,10,19,2,11,20,3,12,21)+3],
                          box_hdist_tmp[,c(1,10,19,2,11,20,3,12,21)+6])
boxplot(box_hdist_df,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_2(min)","lam_2(max)","lam_2(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_1(min)","lam_1(max)","lam_1(aic)"),3))
legend("topright",legend=c("lam_3(min)","lam_3(max)","lam_3(aic)"),col=c(2:4),pch=rep(20,3))

#### 6
box_hdist_df = data.frame(box_hdist_tmp[,c(1,10,19,4,13,22,7,16,25)],
                          box_hdist_tmp[,c(1,10,19,4,13,22,7,16,25)+1],
                          box_hdist_tmp[,c(1,10,19,4,13,22,7,16,25)+2])
boxplot(box_hdist_df,outline=F,ylab="Hellinger Distance",xaxt="n",col=rep(c(2:4),3))
abline(v=c(1:8)*3+0.5,lty=2,col="gray")
abline(v=c(1:2)*9+0.5,lwd=1)
axis(3, at=c(5,14,23), labels=c("lam_1(min)","lam_1(max)","lam_1(aic)"))
axis(1, at=c(1:9)*3-1, labels=rep(c("lam_2(min)","lam_2(max)","lam_2(aic)"),3))
legend("topright",legend=c("lam_3(min)","lam_3(max)","lam_3(aic)"),col=c(2:4),pch=rep(20,3))


# * point * #
S_map2[S_num,]
lam.min.vec


### fixed
rmse_mean = apply(rmse_mat,2,mean)
col_sc =c("red","orange","yellow","green","blue","purple","pink")
# kappa를 고정했을 때 sigma에 따른 mu의 변화 - 뚜렷함
par(mfrow=c(2,3))
for (j in 1:6){
plot(rmse_mean[c(1:6)],type="n",ylim=c(1,4),col="gray",main=paste("kappa=",lam_set2[j]),xlab="mu")
  for (i in 1:6){
    lines(rmse_mean[6^2*(j-1)+c(1:6)+6*(i-1)],type="l",col=col_sc[i])
  }
}
# sigma를 고정했을 때 kappa에 따른 mu의 변화 - 뚜렷하지 않음
par(mfrow=c(2,3))
for (j in 1:6){
  plot(rmse_mean[c(1:6)],type="n",ylim=c(1,4),col="gray",main=paste("sigma=",lam_set3[j]),xlab="mu")
  for (i in 1:6){
    lines(rmse_mean[6^2*(i-1)+c(1:6)+6*(j-1)],type="l",col=col_sc[i])
  }
}




# mu를 고정했을 때 kappa에 따른 sigma의 변화 - 뚜렷하지않음
par(mfrow=c(2,3))
for (j in 1:6){
  plot(rmse_mean[c(1:6)],type="n",ylim=c(1,4),col="gray",main=paste("kappa=",lam_set3[j]),xlab="sigma")
  for (i in 1:6){
    lines(rmse_mean[6^2*(i-1)+c(1:6)*6-5 +(j-1)],type="l",col=col_sc[i])
  }
}
# kappa를 고정했을 때 mu에 따른 sigma의 변화 - 뚜렷함
par(mfrow=c(2,3))
for (j in 1:6){
  plot(rmse_mean[c(1:6)],type="n",ylim=c(1,4),col="gray",main=paste("mu=",lam_set1[j]),xlab="sigma")
  for (i in 1:6){
    lines(rmse_mean[6^2*(j-1)+c(1:6)*6-5 +(i-1) ],type="l",col=col_sc[i])
  }
}



# mu를 고정했을 때 sigma에 따른 kappa의 변화- sigma=100일때 뚜렷함
par(mfrow=c(2,3))
for (j in 1:6){
  plot(rmse_mean[c(1:6)],type="n",ylim=c(1,4),col="gray",main=paste("sigma=",lam_set2[j]),xlab="kappa")
  for (i in 1:6){
    lines(rmse_mean[6*(i-1)+ c(1:6)*6^2-35 +(j-1) ],type="l",col=col_sc[i])
  }
}
# sigma를 고정했을 때 mu에 따른 kappa의 변화 - sigma = 0 일때 뚜렷함
par(mfrow=c(2,3))
for (j in 1:6){
  plot(rmse_mean[c(1:6)],type="n",ylim=c(1,4),col="gray",main=paste("mu=",lam_set1[j]),xlab="kappa")
  for (i in 1:6){
    lines(rmse_mean[6*(j-1)+ c(1:6)*6^2-35 +(i-1)],type="l",col=col_sc[i])
  }
}



