rm(list=ls())
setwd("~/GITHUB/gev")
source("./lib/sgev3library.R")
source("./lib/pack.R")
# load("./graduate/sgev3_train1204.RData")
load("./graduate/sgev3_real1209.RData")

n = length(xlist[[1]])
p = 25
lam.grid = as.data.frame(expand.grid(lam_set,lam_set,lam_set))

df.tps = c()
for (i in 1:length(lam_set)){
  lam_s = lam_set[i]
  Hatmat= Z%*%solve(t(Z)%*%Z +lam_s*Om+diag(1e-08,nrow(Om)))%*%t(Z)
  df.tps[i] =  sum(diag(Hatmat))
}
DF = rowSums(as.matrix(expand.grid( df.tps, df.tps, df.tps)))

nll = c()
for (i in 1:length(result_list)){
  tvec=result_list[[i]]
  
  loc.vec.reg = tvec[3+(1:p)]
  sc.vec.reg = tvec[(3+p)+(1:p)]
  sh.vec.reg = tvec[(3+2*p)+(1:p)]
  loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
  sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
  sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
  
  ll = 0
  for (s in 1:length(xlist)){
    ll = ll - sum(dgev(x=xlist[[s]],loc=loc.vec[s],scale=sc.vec[s],shape=sh.vec[s], log=TRUE))
  }
  nll[i] = ll
}

AIC = 2*nll + 2*DF
# BIC = 2*nll + log(length(train_46))*DF
# par(mfrow=c(1,1))
# plot(AIC,type="l")
# plot(BIC,type="l")

min.ind = which.min(AIC)
lam.grid[min.ind,]

### test set
# unique(test_46$stnlds)

# unique(test_46[test_46$long > max(train_46$long),]$stnlds)
# test_46_tmp = test_46
# test_46 = test_46_tmp %>% filter(!stnlds %in% c(130,277))
# unique(test_46[test_46$long > max(train_46$long),]$stnlds)

# ss <- split.data.frame(test_46,test_46$stnlds)  # stnlds로 dataframe 쪼개서 list에 분배
# xlist.test <- lapply(ss,"[[","pr")    
# ns.test <- length(unique(test_46$stnlds))
# 
# point.est = unlist(lapply(xlist.test,function(x) fgev(x)$estimate))
# point.loc = point.est[3*(1:ns.test)-2]
# point.sc =  point.est[3*(1:ns.test)-1]
# point.sh =  point.est[3*(1:ns.test)]
# 
# zlist <- list()
# for (i in 1:ns.test){
#   xbs <- eval.basis(unique(ss[[i]]$long),x_bsobj)
#   ybs <- eval.basis(unique(ss[[i]]$lat),y_bsobj)
#   tensorbs <- do.call('cbind', lapply(1:ncol(xbs), function(i) xbs[, i] * ybs))  # row-wise kronecker product
#   zlist[[i]] <- tensorbs 
# }
# Z.test = do.call('rbind',zlist)
# 
# p=25
# tvec=result_list[[min.ind]]
# loc.vec.reg = tvec[3+(1:p)]
# sc.vec.reg = tvec[(3+p)+(1:p)]
# sh.vec.reg = tvec[(3+2*p)+(1:p)]
# loc.vec = tvec[1] + drop(Z%*%loc.vec.reg)
# sc.vec = exp(tvec[2] + drop(Z%*%sc.vec.reg))
# sh.vec = tvec[3] + drop(Z%*%sh.vec.reg)
# par(mfrow=c(4,3))
# for ( s in 1:ns){    
#   Fn = ecdf(xlist[[s]])
#   z_hat = qgev(p=Fn(xlist[[s]])*n/(n+1), loc=loc.vec[s], sc=sc.vec[s], sh=sh.vec[s])    
#   plot(z_hat,xlist[[s]],xlab="Theoretical quantiles",ylab="Sample quantiles",
#        col="blue",main=unique(train_46$stnlds)[s])
#   abline(a=0,b=1,col="black")
# }
# 
# lam.min.star = lam.grid[min.ind,]
# 
#     tvec=result_list[[min.ind]]
#     loc.vec.reg = tvec[3+(1:p)]
#     sc.vec.reg = tvec[(3+p)+(1:p)]
#     sh.vec.reg = tvec[(3+2*p)+(1:p)]
#     loc.vec = tvec[1] + drop(Z.test%*%loc.vec.reg)
#     sc.vec = exp(tvec[2] + drop(Z.test%*%sc.vec.reg))
#     sh.vec = tvec[3] + drop(Z.test%*%sh.vec.reg)
# 
# par(mfrow=c(4,3))
#   for ( s in 1:ns.test){    
# 
#   Fn = ecdf(xlist.test[[s]])
#   z_hat = qgev(p=Fn(xlist.test[[s]])*n/(n+1), loc=loc.vec[s], sc=sc.vec[s], sh=sh.vec[s])    
#   
#   plot(z_hat,xlist.test[[s]],xlab="Theoretical quantiles",ylab="Sample quantiles",
#        col="blue",main=unique(test_46$stnlds)[s])
#   abline(a=0,b=1,col="black")
#   }
# 
# par(mfrow=c(1,1))
# plot(train_46$long,train_46$lat)
# points(test_46$long,test_46$lat,col="red")
# text(test_46$long,test_46$lat, labels= test_46$stnlds , cex= 0.7)


plot.long = seq(x_bsobj$rangeval[1],x_bsobj$rangeval[2],length=100)
plot.lat = seq(y_bsobj$rangeval[1],y_bsobj$rangeval[2],length=100)
gg = expand.grid(plot.long,plot.lat)
xbs = eval.basis(gg[,1],x_bsobj)
ybs = eval.basis(gg[,2],y_bsobj)
z <- do.call('cbind', lapply(1:ncol(xbs), function(k) xbs[,k] * ybs))  # row-wise kronecker product
tvec=result_list[[min.ind]]
loc.vec.reg = tvec[3+(1:p)]
sc.vec.reg = tvec[(3+p)+(1:p)]
sh.vec.reg = tvec[(3+2*p)+(1:p)]
gg$loc = tvec[1] + drop(z %*%loc.vec.reg)
gg$sc = exp(tvec[2] + drop(z %*%sc.vec.reg))
gg$sh = tvec[3] + drop(z %*%sh.vec.reg)

# plot_ly(x=gg$Var1,y=gg$Var2,z=gg$loc, type="heatmap")
# plot_ly(x=gg$Var1,y=gg$Var2,z=gg$sc, type="heatmap")
# plot_ly(x=gg$Var1,y=gg$Var2,z=gg$sh, type="heatmap")

library(sf)
library(sp)
korea_shp <- st_read("./shapefile/SouthKoreaboundary96.shp") # sf packages
# korea_shp$CTP_KOR_NM <- iconv(korea_shp$CTP_KOR_NM, from = "CP949", to = "UTF-8", sub = NA, mark = TRUE, toRaw = FALSE)
# korea_shp %>% select(id) %>% plot()
st_crs(korea_shp)  # proj4stiring
korea_shp_df <-  as(korea_shp, 'Spatial')  # sf dataframe을 shapefile로 변환
class(korea_shp_df)

# Removing data outside country map boundary
spdf <- SpatialPointsDataFrame(coords = gg[, c("Var1", "Var2")], data = gg,
                               proj4string = CRS( "+proj=longlat +datum=WGS84 +no_defs"))
whatever <- spdf[!is.na(over(spdf, as(korea_shp_df, "SpatialPolygons"))), ]
whatever <- as.data.frame(whatever)

# library(RColorBrewer)
# display.brewer.all()

g1 <- ggplot() +
  geom_polygon(data = korea_shp_df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.5, fill = "white") +
  geom_tile(data = whatever, aes(x = Var1, y = Var2, fill = loc), alpha = 0.8) +
  labs(title = "", x = "경도", y = "위도") +
  scale_fill_distiller(type = "div", palette =  "Spectral") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 15)) +
  theme(legend.title=element_blank()) +
  coord_map(xlim=c(125,130),ylim=c(33.7,39))


g2 <- ggplot() +
  geom_polygon(data = korea_shp_df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.5, fill = "white") +
  geom_tile(data = whatever, aes(x = Var1, y = Var2, fill = sc), alpha = 0.8) +
  labs(title = "", x = "경도", y = "위도") +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 15)) +
  theme(legend.title=element_blank()) +
  coord_map(xlim=c(125,130),ylim=c(33.7,39))



g3 <- ggplot() +
  geom_polygon(data = korea_shp_df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.5, fill = "white") +
  geom_tile(data = whatever, aes(x = Var1, y = Var2, fill = sh), alpha = 0.8) +
  labs(title = "", x = "경도", y = "위도") +
  scale_fill_distiller(type = "div", palette = "Spectral" ) +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 15)) +
  theme(legend.title=element_blank()) +
  coord_map(xlim=c(125,130),ylim=c(33.7,39))

# library(gridExtra)
# grid.arrange(g1, g2, g3, ncol=3)

g1
g2
g3

# Return level
# install.packages("mev")
library(mev)
head(gg)
r10 = c(); r20=c(); r50=c(); r100=c(); r5=c()
for ( i in 1:nrow(gg)) {
  r5[i] = gev.retlev(par=gg[i,3:5],1/2)
  # r10[i] = gev.retlev(par=gg[i,3:5],1/10)
  # r20[i] = gev.retlev(par=gg[i,3:5],1/20)
  # r50[i] = gev.retlev(par=gg[i,3:5],1/50)
  # r100[i] = gev.retlev(par=gg[i,3:5],1/100)
}
gg$return10 = unlist(r10)
gg$return20 = unlist(r20)
gg$return50 = unlist(r50)
gg$return100 = unlist(r100)
gg$return1 = unlist(r5)

spdf <- SpatialPointsDataFrame(coords = gg[, c("Var1", "Var2")], data = gg,
                               proj4string = CRS( "+proj=longlat +datum=WGS84 +no_defs"))
whatever <- spdf[!is.na(over(spdf, as(korea_shp_df, "SpatialPolygons"))), ]
whatever <- as.data.frame(whatever)

r1 <- ggplot() +
  geom_polygon(data = korea_shp_df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.5, fill = "white") +
  geom_tile(data = whatever, aes(x = Var1, y = Var2, fill = return10), alpha = 0.8) +
  labs(title = "", x = "경도", y = "위도") +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 15)) +
  theme(legend.title=element_blank()) +
  coord_map(xlim=c(125,130),ylim=c(33.7,39))


r2 <- ggplot() +
  geom_polygon(data = korea_shp_df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.5, fill = "white") +
  geom_tile(data = whatever, aes(x = Var1, y = Var2, fill = return20), alpha = 0.8) +
  labs(title = "", x = "경도", y = "위도") +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 15)) +
  theme(legend.title=element_blank()) +
  coord_map(xlim=c(125,130),ylim=c(33.7,39))


r3 <- ggplot() +
  geom_polygon(data = korea_shp_df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.5, fill = "white") +
  geom_tile(data = whatever, aes(x = Var1, y = Var2, fill = return50), alpha = 0.8) +
  labs(title = "", x = "경도", y = "위도") +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 15)) +
  theme(legend.title=element_blank()) +
  coord_map(xlim=c(125,130),ylim=c(33.7,39))


r4 <- ggplot() +
  geom_polygon(data = korea_shp_df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.5, fill = "white") +
  geom_tile(data = whatever, aes(x = Var1, y = Var2, fill = return100), alpha = 0.8) +
  labs(title = "", x = "경도", y = "위도") +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 15)) +
  theme(legend.title=element_blank()) +
  coord_map(xlim=c(125,130),ylim=c(33.7,39))



r5 <- ggplot() +
  geom_polygon(data = korea_shp_df, aes(x = long, y = lat, group = group),
               colour = "black", size = 0.5, fill = "white") +
  geom_tile(data = whatever, aes(x = Var1, y = Var2, fill = return1), alpha = 0.8) +
  labs(title = "", x = "경도", y = "위도") +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 15)) +
  theme(legend.title=element_blank()) +
  coord_map(xlim=c(125,130),ylim=c(33.7,39))
