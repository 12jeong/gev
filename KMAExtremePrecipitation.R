rm( list = ls()); gc()

if(!require(dplyr)){ install.packages('dplyr')}; require(dplyr)
if(!require(data.table)){ install.packages('data.table')}; require(data.table)
# if(!require(chron)){ install.packages('chron')}; require(chron)
# if(!require(lattice)){ install.packages('lattice')}; require(lattice)
# if(!require(RColorBrewer)){ install.packages('RColorBrewer')}; require(RColorBrewer)
if(!require(maps)){ install.packages('maps')}; require(maps) 
# if(!require(sos)){ install.packages('sos')}; require(sos)
# ???tps
# if(!require(ggformula)){ install.packages('ggformula')}; require(ggformula)
# if(!require(splines)){ install.packages('splines')}; require(splines)
# detach(package:splines)
# if(!require(ggplot2)){ install.packages('ggplot2')}; require(ggplot2)
# if(!require(reshape2)){ install.packages('reshpae2')}; require(reshape2)

setwd("C:\\Users\\UOS\\Dropbox\\Extreme value\\gev")

# load("C:\\Users\\UOS\\서울시립대학교\\전종준 - lab_work\\Lab_process\\KEY\\kma_data\\kma_data.Rdata")
# 
# head(kma_data)    # sum_rn : 일강수량 (mm)
# 
# table(kma_data$stnlds) # 102개 지점
# 
# table(substr(kma_data$time,1,4))
# 
# kma_data$stnlds <- as.Date(kma_data$stnlds)
# 
# summary(kma_data)
# 
# ### 지점별 연최대강수량을 뽑자
# Pr_data <- kma_data %>% group_by(stnlds,obsyear=substr(kma_data$time,1,4),lat,long) %>% summarise(pr=max(sum_rn,na.rm=TRUE))
# Pr_data <- Pr_data[!is.infinite(Pr_data$pr),]
# 
# ### 지점별로 몇년도부터 데이터가 있는지 알고싶다
# table(Pr_data$stnlds) # 102개 지점
# plot(Pr_data$obsyear,Pr_data$stnlds)
# 
# ### 1973년부터 2018년 : 46개년
# table(Pr_data$stnlds)
# stn1 <- names(table(Pr_data$stnlds))[c(unname(which(table(Pr_data$stnlds)>45)))] # 46개년 이상인 지점번호
# Pr_48 <- Pr_data %>% filter(stnlds %in% stn1 && obsyear > 1972)
# table(Pr_48$stnlds) # 60개 지점
# Pr_48 %>% group_by(stnlds,lat,long) %>% summarise(minyear=min(as.numeric(obsyear)), maxyear=max(as.numeric(obsyear))) %>% View
# Pr_48 <- as.data.frame(Pr_48)

# save(Pr_48, file="Pr_48.RData")


##################################################################################
##################################################################################
load("Pr_48.RData")

### 지점별 연최대강수량 평균
# Pr_avg <- Pr_48 %>% group_by(stnlds,lat,long) %>% summarise(avgpr = mean(pr))
# 
# if(!require(akima)){ install.packages('akima')}; require(akima)
# fld <- with(Pr_avg,interp(x=long,y=lat,z=avgpr))
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette = colorRampPalette(c("white", "blue")),
#                xlab = "Longitude",
#                ylab = "Latitude",
#                main = "Maximum Precipitation",
#                key.title = title(main = "Rain (mm)", cex.main = 1))
# 
# maps:::map(database='world',region = "South Korea",lwd=2,add=T)
# points(Pr_48$long,Pr_48$lat,pch=4)
# 
# ?filled.contour

### 년도별로 봐볼까
# if(!require(akima)){ install.packages('akima')}; require(akima)
# year_contour = function(vyear){
# pdt <- Pr_48 %>% filter(obsyear==vyear)
# fld <- with(pdt,interp(x=long,y=lat,z=pr))
# 
# filled.contour(x = fld$x,
#                y = fld$y,
#                z = fld$z,
#                color.palette = colorRampPalette(c("white", "blue")),
#                xlab = "Longitude",
#                ylab = "Latitude",
#                main =vyear,
#                key.title = title(main = "Rain (mm)", cex.main = 1))
# 
# maps:::map(database='world',region = "South Korea",lwd=2,add=T)
# points(Pr_48$long,Pr_48$lat,pch=4)
# }

# year_contour(1980) # 1980년 7월 21일 집중호우
# year_contour(1984) # 1984년 8월 31일 서울 대홍수
# year_contour(1986) # 1986년 8월 27일 태풍 베라
# year_contour(1990) # 1990 9월 한강 대홍수
# year_contour(1996) # 대홍수
# year_contour(2000) # 2000년 8월 23일 태풍 프라피룬
# year_contour(2002) # 태풍 루사 (강릉지방 8월 30일 극값)
# year_contour(2011) # 수도권 폭우
# year_contour(2014) # 동남권 폭우
# year_contour(2017) # 중부권, 동남권 폭우
# year_contour(1984) # 1984년 8월 31일 서울 대홍수
# year_contour(2018) # 국지성폭우

# Pr_avg <- Pr_avg[with(Pr_avg,order(lat,long)),] 
# image(Pr_avg$long,Pr_avg$lat,Pr_avg$avgpr, col=rev(brewer.pal(10,"RdBu")))


##################################################################################
##################################################################################

### library - fields
if(!require(fields)){ install.packages('fields')}; require(fields)

load("Pr_48.RData")

##################################
############ 이거 아님 ###########
fit <- Tps(x=cbind(Pr_48$long,Pr_48$lat),Y=Pr_48$pr)
summary(fit)

set.panel(2,2)
plot(fit)

# matrix(as.vector(fit$fitted.values),nrow=46)[1:3,1:3]
# matrix(as.vector(fit$fitted.values),nrow=46)[1:3,58:60]
set.panel()
surface(fit)
map(database='world',region = "South Korea",lwd=2,add=T)

# hatz <- c()
# for(i in 1:60){
#   hatz[i] <- fit$fitted.values[i*46-45]
# }

# ?surface
# 
# cbind(Pr_avg$avgpr,hatz)
# summary(Pr_avg$avgpr)
# summary(hatz)

cbind(Pr_48,fit$fitted.values)
Tps_fit <- matrix(as.vector(fit$fitted.values),nrow=46)
colnames(Tps_fit) <- paste0("L",unique(Pr_48$stnlds))
dim(Tps_fit)
Tps_fit[1:3,1:3]
Tps_fit[44:46,1:3]


fgev(Pr_48$pr - fit$fitted.values)
localno=108
fgev( Pr_48[Pr_48$stnlds==localno,]$pr - fit$fitted.values[which(unique(Pr_48$stnlds)==localno)*46] )
##################################
##################################

Tps_fitted <- matrix(ncol=60,nrow=46)
for (i in 1:46){
Pr_oneyear <- Pr_48 %>% filter(obsyear==1972+i) 
fit <- Tps(x=cbind(Pr_oneyear$long,Pr_oneyear$lat),Y=Pr_oneyear$pr)  # 연도마다 60개 지점 mu 추정
Tps_fitted[i,] =as.vector(fit$fitted.values)
}
colnames(Tps_fitted) <- paste0("L",unique(Pr_48$stnlds))
dim(Tps_fitted)
Tps_fitted[1:3,1:3]
Tps_fitted[44:46,1:3]

localno=108
fgev( Pr_48[Pr_48$stnlds==localno,]$pr - Tps_fitted[ ,which(unique(Pr_48$stnlds)==localno)] )

fgev( Pr_48$pr - as.vector(Tps_fitted))


##################################################################################
### library - mgcv
if(!require(mgcv)){ install.packages('mgcv')}; require(mgcv) 

# model <- gam(pr~ s(long, bs = 'tp', k = 4, fx = TRUE)
#                   + s(lat, bs = 'tp', k = 4, fx = TRUE)
#                   + ti(long, lat, bs = 'tp', k = c(4, 4), d = c(1, 1), fx = TRUE),
#                     data = Pr_48)

model <- gam(pr~ s(long)+ s(lat), data = Pr_48)


vis.gam(model, plot.type="contour",color="topo")
maps:::map(database='world',region = "South Korea",lwd=2,add=T)

?vis.gam
surface(model)

?gam





##################################################################################
##################################################################################
# 울릉도 제주도 빼봅시다

unique( Pr_48[Pr_48$long==max(Pr_48$long),]$stnlds) #울릉도 지점코드 = 115
unique( Pr_48[Pr_48$lat==min(Pr_48$lat),]$stnlds)   #제주도 지점코드 = 184, 188, 189

Pr2 <- Pr_48 %>% filter(Pr_48$long <= 130, Pr_48$lat >= 34)

Pr_avg2 <- Pr2 %>% group_by(stnlds,lat,long) %>% summarise(avgpr = mean(pr))

# if(!require(akima)){ install.packages('akima')}; require(akima)
fld2 <- with(Pr_avg2,interp(x=long,y=lat,z=avgpr))

filled.contour(x = fld2$x,
               y = fld2$y,
               z = fld2$z,
               color.palette = colorRampPalette(c("white", "blue")),
               xlab = "Longitude",
               ylab = "Latitude",
               main = "Maximum Precipitation",
               key.title = title(main = "Rain (mm)", cex.main = 1))

maps:::map(database='world',region = "South Korea",lwd=2,add=T)
points(Pr2$long,Pr2$lat,pch=4)

##################################################################################
fit2 <- Tps(x=cbind(Pr2$long,Pr2$lat),Y=Pr2$pr)
summary(fit2)

set.panel(2,2)
plot(fit2)

set.panel()
surface(fit2)
map(database='world',add=T,lwd=2)

localno=90
fgev(  Pr2[Pr2$stnlds==localno,]$pr - fit2$fitted.values[which(unique(Pr2$stnlds)==localno):which(unique(Pr2$stnlds)==localno)*46])

