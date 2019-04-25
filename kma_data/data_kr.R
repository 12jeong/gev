rm(list=ls())
setwd("C:\\Users\\UOS\\Dropbox\\Extreme value\\kma_data")

if(!require(lubridate)){install.packages("lubridate")}; require(lubridate)
if(!require(dplyr)){install.packages("dplyr")}; require(dplyr)
if(!require(data.table)){install.packages("data.table")}; require(data.table)

file_1 <- read.csv("1900-1909.csv")
file_2 <- read.csv("1910-1919.csv")
file_3 <- read.csv("1920-1929.csv")
file_4 <- read.csv("1930-1939.csv")
file_5 <- read.csv("1940-1949.csv")
file_6 <- read.csv("1950-1959.csv")
file_7 <- read.csv("1960-1969.csv")
file_8 <- read.csv("1970-1979.csv")
file_9 <- read.csv("1980-1989.csv")
file_10 <- read.csv("1990-1999.csv")
file_11 <- read.csv("2000-2009.csv")
file_12 <- read.csv("2010-2018.csv")
data_kr_1 <- rbind(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_9,file_10,file_11,file_12)
colnames(data_kr_1) <- c("stnlds","time","min_ta","max_ta","sum_rn","avg_ws","avg_rhm","sum_gsr")
colnames(file_8) <- c("stnlds","time","min_ta","max_ta","sum_rn","avg_ws","avg_rhm","sum_gsr")
data_kr<- rbind(data_kr_1,file_8)
head(data_kr)

# data_gsr <- data_kr %>% filter(is.na(sum_gsr)==FALSE)
# unique(year(data_gsr$time))


# 지점파일과 merge

stnlds <- fread("stnlds.csv", stringsAsFactors = T)
data_stnlds <- stnlds %>% filter(duplicated(지점) == F) %>% select(지점, 위도, 경도)
colnames(data_stnlds) <- c("stnlds", "lat", "long")

kma_data <- merge(x = data_kr, y = data_stnlds, by = "stnlds", all.x = T)

save(kma_data, file = "kma_data.RData")

