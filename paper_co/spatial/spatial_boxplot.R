rm( list = ls()); gc()
setwd("~/GITHUB/gev")

# To load the file and make RMSE / AIC / HD mat
Ns=50
load(file="./Rexport/Spatial_simulation/result1_analysis.RData")
SSEmat1= SSEmat_all
RMSEmat1 = SSEmat1/sqrt(3*Ns)
AICmat1 = AICmat
# BICmat1 = BICmat
distmat1 = distmat
load(file="./Rexport/Spatial_simulation/result2_analysis.RData")
SSEmat2 = SSEmat_all
RMSEmat2 =  SSEmat2/sqrt(3*Ns)
AICmat2 = AICmat[-2,]
# BICmat2 = BICmat
distmat2 = distmat
load(file="./Rexport/Spatial_simulation/result3_analysis.RData")
SSEmat3 = SSEmat_all
RMSEmat3 =  SSEmat3/sqrt(3*Ns)
AICmat3 = AICmat[-2,]
# BICmat3 = BICmat
distmat3 = distmat

# To select best lambda
par(mfrow=c(1,3))
barplot(table(apply(AICmat1,2,which.min)))
barplot(table(apply(AICmat2,2,which.min)))
barplot(table(apply(AICmat3,2,which.min)))
table(apply(AICmat1,2,which.min))/100
table(apply(AICmat2,2,which.min))/100
table(apply(AICmat3,2,which.min))/100
lam.min1 = as.numeric(names(which.max(table(apply(AICmat1,2,which.min)))))
lam.min2 = as.numeric(names(which.max(table(apply(AICmat2,2,which.min)))))
lam.min3 = as.numeric(names(which.max(table(apply(AICmat3,2,which.min)))))

# Data frame for RMSE boxplot
df_RMSE1 = data.frame(RMSEmat1[1,],RMSEmat1[12,],RMSEmat1[lam.min1,])
df_RMSE2 = data.frame(RMSEmat2[1,],RMSEmat2[10,],RMSEmat2[lam.min2,])
df_RMSE3 = data.frame(RMSEmat3[1,],RMSEmat3[10,],RMSEmat3[lam.min3,])

par(mfrow=c(1,3))
op <- par(cex = 0.75)
boxplot(df_RMSE1,main="Hyperplane",ylim=c(3.5,5.8),
        ylab="RMSE",
        names=c(expression(paste(lambda,'=',0)),
                expression(paste(lambda,'=',100)),
                bquote(paste(lambda,"*"))))
boxplot(df_RMSE2,main="Unimodal",
        ylab="RMSE",
        names=c(expression(paste(lambda,'=',0)),
                expression(paste(lambda,'=',100)),
                bquote(paste(lambda,"*"))))
boxplot(df_RMSE3,main="Bimodal",
        ylab="RMSE",
        names=c(expression(paste(lambda,'=',0)),
                expression(paste(lambda,'=',100)),
                bquote(paste(lambda,"*"))))

# Best Lambda Index according to each example
id1 = apply(AICmat1,2,which.min)
id2 = apply(AICmat2,2,which.min)
id3 = apply(AICmat3,2,which.min)

# Data frame for AIC boxplot
a1=c(); a2=c(); a3=c()
for(i in 1:ncol(AICmat1)) a1[i] = AICmat1[id1[i],i]
for(i in 1:ncol(AICmat2)) a2[i] = AICmat2[id2[i],i]
for(i in 1:ncol(AICmat3)) a3[i] = AICmat3[id3[i],i]
df_AIC1 = data.frame(AICmat1[1,],AICmat1[12,],AICmat1[lam.min1,],a1)
df_AIC2 = data.frame(AICmat2[1,],AICmat2[9,],AICmat2[lam.min2,],a2)
df_AIC3 = data.frame(AICmat3[1,],AICmat3[9,],AICmat3[lam.min3,],a3)
boxplot(df_AIC1,main="Hyperplane",
        ylab="AIC",
        names=c(expression(paste(lambda,'=',0)),
                expression(paste(lambda,'=',100)),
                bquote(lambda[AIC]),
                bquote(paste(lambda,"*"))))
boxplot(df_AIC2,main="Unimodal",
        ylab="AIC",
        names=c(expression(paste(lambda,'=',0)),
                expression(paste(lambda,'=',100)),
                bquote(lambda[AIC]),
                bquote(paste(lambda,"*"))))
boxplot(df_AIC3,main="Bimodal",
        ylab="AIC",
        names=c(expression(paste(lambda,'=',0)),
                expression(paste(lambda,'=',100)),
                bquote(lambda[AIC]),
                bquote(paste(lambda,"*"))))

# Dataframe for Hellinger distance 
d1=c(); d2=c(); d3=c()
for(i in 1:ncol(distmat1)) d1[i] = distmat1[id1[i],i]
for(i in 1:ncol(distmat2)) d2[i] = distmat2[id2[i],i]
for(i in 1:ncol(distmat3)) d3[i] = distmat3[id3[i],i]
df_dist1 = data.frame(distmat1[1,],distmat1[12,],distmat1[lam.min1,],d1)
df_dist2 = data.frame(distmat2[1,],distmat2[10,],distmat2[lam.min2,],d2)
df_dist3 = data.frame(distmat3[1,],distmat3[10,],distmat3[lam.min3,],d3)
boxplot(df_dist1, main = "Hyperplane", outline=FALSE,
        ylab="Hellinger distance",
        names=c(expression(paste(lambda,'=',0)),
                expression(paste(lambda,'=',100)),
                bquote(lambda[AIC]),
                bquote(paste(lambda,"*"))))
boxplot(df_dist2, main = "Unimodal", 
        ylab="Hellinger distance",
        names=c(expression(paste(lambda,'=',0)),
                expression(paste(lambda,'=',100)),
                bquote(lambda[AIC]),
                bquote(paste(lambda,"*"))))
boxplot(df_dist3, main = "Bimodal", 
        ylab="Hellinger distance",
        names=c(expression(paste(lambda,'=',0)),
                expression(paste(lambda,'=',100)),
                bquote(lambda[AIC]),
                bquote(paste(lambda,"*"))))

