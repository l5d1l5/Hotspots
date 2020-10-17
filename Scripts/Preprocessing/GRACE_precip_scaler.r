## this script plots the cumulative area percentage encapsulated
## by various area-weighted percentile thresholds in absolute TWSt/MAP trends
## to select the max clip threshold for indicator plotting/development

library(tidyverse); library(magrittr); library(here)
library(raster); library(sf); library(spatstat)
library(fasterize); library(ggplot2); 

# Import grid areas, TWSt, and MAP
twst <- raster(here::here("Data", "Rodell_GRACE_TWSt_raw.tif"))*10 # cm to mm
MAP <- raster(here::here("ProducedData", "MeanAnnualPrecip_1985_2014.tif"))
ga <-  raster(here::here("ProducedData", "WGS_gridArea_halfdegree.tif"))
Mask <- raster(here::here("ProducedData", "GlobalMask.tif"))

twst[Mask == 0] <- NA
twst[1,] <- NA;twst[,1] <- NA;
MAP[Mask == 0] <- NA
MAP[1,] <- NA;MAP[,1] <- NA;
TWSt_MAP <- twst/(MAP)
TWSt_MAP[MAP == 0] <- NA

df.stat <- cbind(TWSt_MAP[], ga[]) %>% as.data.frame() %>% set_colnames(c("twst", "ga"))
df.stat <- df.stat[complete.cases(df.stat),]
# calculate inter 10th and 90th percentiles instead of stdev as skewed
x = 0.4
lim = Hmisc::wtd.quantile(x = df.stat$twst, weights = df.stat$ga, probs = 0.5+x) - 
  Hmisc::wtd.quantile(x = df.stat$twst, weights = df.stat$ga, probs = 0.5-x) 

TWSt_MAP_scaler <- raster(ga)
TWSt_MAP_scaler[] <- TWSt_MAP[]/lim
TWSt_MAP_scaler[TWSt_MAP_scaler < -1] <- -1
TWSt_MAP_scaler[TWSt_MAP_scaler > 1] <- 1
plot(TWSt_MAP_scaler)

plot(TWSt_MAP_scaler, col = plt, zlim = c(-1, 1))

writeRaster(TWSt_MAP_scaler, here::here("ProducedData", "TWS_MAP_sclaer_p10p90_1985_2014map.tif"),
            format = "GTiff", overwrite = T)

# below if elbow method we want to use
################################ 
TWSt.precip <- TWSt/MAP
TWSt.precip[MAP == 0] <- NA
TWSt.precip <- abs(TWSt.precip)

# load earthquake interference regions to mask
EQ_int <- sf::read_sf(here::here("Data", "EQ_int", "EQ_int.shp"))
EQ_int.r <- fasterize(sf = EQ_int, raster = ga, fun = "last")
TWSt.precip[EQ_int.r == 1] <- NA

# couple rasters for analysis
m.df <- cbind(as.data.frame(TWSt.precip), as.data.frame(ga)) %>% 
  set_colnames(c("ratio", "area"))
m.df <- m.df[complete.cases(m.df$ratio),]

# create results dataframe
r.df <- data.frame(x = c(seq(1,100, by = 1)), 
                   y = c(rep(NA, 100)),
                   z = c(rep(NA, 100))) %>% 
  set_colnames(c("Ptile", "PctArea", "Value"))

for (i in 1:nrow(r.df)) {
  tv <- r.df$Ptile[i]
  
  thrsh <- spatstat::weighted.quantile(x = m.df$ratio, w = m.df$area,
                                       probs = tv/100, na.rm = T)
  
  ar <- m.df %>% filter(ratio <= thrsh) %>% pull(area) %>% sum()
  p.ar <- ar/sum(m.df$area)
  r.df$PctArea[i] <- p.ar
  r.df$Value[i] <- thrsh
  print(i/100)
  
}

# confirm 1:1 percentile to percent area results
# plot(r.df$PctArea ~ r.df$Ptile)

# Create straight line between 1st and 99th percentile to find "elbow"
for(i in 1:nrow(r.df)){
  r.df$line[i] <- r.df$Value[1] + (r.df$Value[99] - r.df$Value[1])*((i-1)/(99-1))
}
r.df$line[100] <- NA

# Find nearest distance to line for all points
X1 = r.df$Ptile[1]
X2 = r.df$Ptile[99]
Y1 = r.df$Value[1]
Y2 = r.df$Value[99]

for(i in 1:nrow(r.df)){
  num = (X2-X1)*(Y1-r.df$Value[i]) - (X1-r.df$Ptile[i])*(Y2-Y1)
  den = sqrt( (X2-X1)^2 + (Y2-Y1)^2 )
  r.df$dist[i] <- abs(num)/den
}
r.df$dist[100] <- NA

max(r.df$dist, na.rm = T)

mv = r.df[which.max(r.df$dist),]
ggplot(data = r.df, aes(x = Ptile))+
  geom_point( aes(y = Value), size = 3) +
  geom_line( aes(y = line), size = 2, col = "red") +
  geom_point( aes(x = mv$Ptile, y = mv$Value), size = 3, col = "red") +
  coord_cartesian(ylim = c(0,0.8), xlim = c(1,99))

## now that elbow value is established, scale TWSt.precip raster to these limits
TWSt.pre.scale <- (TWSt/MAP)/mv$Value
TWSt.pre.scale[MAP == 0] <- NA
TWSt.pre.scale[EQ_int.r == 1] <- NA
TWSt.pre.scale[TWSt.pre.scale < -1] <- -1
TWSt.pre.scale[TWSt.pre.scale > 1] <- 1
plot(TWSt.pre.scale)

writeRaster(TWSt.pre.scale, here::here("ProducedData", "TWSt_precip_scale.tif"),
            format = "GTiff", overwrite = T)
