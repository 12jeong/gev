rm(list=ls())
setwd("~/GITHUB/gev")
library(ggplot2)
library(ggmap)
library(raster)
library(rgeos)
library(maptools)
library(maps)
library(knitr)
library(rgdal)
load("./kma_data/Pr_46.RData")
# library(sf)
# korea_shp <- st_read("./shapefile/SouthKoreaboundary96.shp") 
# st_crs(korea_shp) 
korea_shp <- shapefile("./shapefile/TL_SCCO_CTPRVN.shp") 
korea <- spTransform(korea_shp, CRS("+proj=longlat +datum=WGS84 +no_defs"))
korea_map <- fortify(korea)
# korea <- map_data(map = 'world', region = 'South Korea')

df = data.frame(lat=unique(Pr_46$lat),long=unique(Pr_46$long),stnlds=unique(Pr_46$stnlds))
require("ggrepel")
rainbow.n = rainbow(8, s = 1, v = 1, start = 0, end = max(1,8 - 1)/8, alpha = 0.7)
set.seed(100)
ggplot(data=df, aes(x=long,y=lat)) +
  geom_polygon(data =korea_map, mapping = aes(x = long,y = lat,group = group), 
               fill = 'white', color="gray40") +
  geom_point(color="blue", size=3) +
  geom_point(colour="white", size=1.5) +
  geom_label_repel(aes(label=stnlds),  segment.color="blue",
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.3, "lines")) +
  theme_classic(base_size = 16) +
  # geom_text(data=df, aes(x=long+0.02, y=lat+0.02, label=stnlds), 
  #           colour="blue", hjust=0, vjust=0, size=4) + 
  labs(title = "", x = "경도", y = "위도") +
  coord_map(xlim=c(125,130),ylim=c(33.7,39)) +
  theme_bw() +
  theme(plot.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2)) 
