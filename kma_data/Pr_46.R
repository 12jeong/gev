rm( list = ls()); gc()
setwd("C:\\Users\\UOS\\Dropbox\\Extreme value\\kma_data")
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

range(Pr_data[Pr_data$stnlds==285,]$obsyear)
### 1973년부터 2018년 : 46개년
stn1 <- names(table(Pr_data$stnlds))[c(unname(which(table(Pr_data$stnlds)>=46)))] # 46개년 이상인 지점번호
Pr_46 <- Pr_data %>% filter(stnlds %in% stn1 && obsyear >= 1973)
table(Pr_46$stnlds) # 60개 지점
Pr_46 %>% group_by(stnlds,lat,long) %>% summarise(minyear=min(as.numeric(obsyear)), maxyear=max(as.numeric(obsyear))) %>% View
Pr_46 <- as.data.frame(Pr_46)

# save(Pr_46, file="Pr_46.RData")

