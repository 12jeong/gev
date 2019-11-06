rm( list = ls()); gc()
setwd("C:\\Users\\UOS\\Documents\\GITHUB\\gev\\kma_data")
load("kma_data.Rdata")

if(!require(dplyr)){ install.packages('dplyr')}; require(dplyr)
if(!require(data.table)){ install.packages('data.table')}; require(data.table)
if(!require(maps)){ install.packages('maps')}; require(maps) 

head(kma_data)    # sum_rn : 일강수량 (mm)
table(kma_data$stnlds) 
length(table(kma_data$stnlds)) # 102개 지점
table(substr(kma_data$time,1,4))
length(table(substr(kma_data$time,1,4))) # 115개 년도


### 지점별 연최대강수량을 뽑자
Pr_data <- kma_data %>% group_by(stnlds,obsyear=substr(kma_data$time,1,4),lat,long) %>% summarise(pr=max(sum_rn,na.rm=TRUE))
Pr_data <- Pr_data[!is.infinite(Pr_data$pr),]

### 지점별로 몇년도부터 데이터가 있는지 알고싶다
table(Pr_data$stnlds) # 102개 지점
plot(Pr_data$obsyear,Pr_data$stnlds)

sum(table(Pr_data$stnlds) > 49) # 50개년 이상인것
which(table(Pr_data$stnlds) == 50) 
range(Pr_data[Pr_data$stnlds==133,]$obsyear)

stn1 <- names(table(Pr_data$stnlds))[c(unname(which(table(Pr_data$stnlds)>=50)))] # 50개년 이상인 지점번호
Pr_50 <- Pr_data %>% filter(stnlds %in% stn1 && obsyear >= 1969)
table(Pr_50$stnlds) # 24개 지점
Pr_50 %>% group_by(stnlds,lat,long) %>% summarise(minyear=min(as.numeric(obsyear)), maxyear=max(as.numeric(obsyear))) %>% View
Pr_50 <- as.data.frame(Pr_50)

# 제주도 2개 지점, 울릉도 1개 지점 제외
unique(Pr_50[Pr_50$lat < 34,]$stnlds)
unique(Pr_50[Pr_50$long > 130,]$stnlds)
Pr_50 <- Pr_50 %>% filter(!stnlds %in% c(184,189,115))

length(table(Pr_50$stnlds)) # 21개 지점

plot(Pr_50$lat,Pr_50$long)
# save(Pr_50, file="Pr_50.RData")


