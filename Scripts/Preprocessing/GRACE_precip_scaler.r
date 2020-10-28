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

# calculate IDR (p90 - p10) instead of stdev as skewed
x = 0.4
lim = Hmisc::wtd.quantile(x = df.stat$twst, weights = df.stat$ga, probs = 0.5+x) - 
  Hmisc::wtd.quantile(x = df.stat$twst, weights = df.stat$ga, probs = 0.5-x) 

TWSt_MAP_scaler <- raster(ga)
TWSt_MAP_scaler[] <- TWSt_MAP[]/lim
TWSt_MAP_scaler[TWSt_MAP_scaler < -1] <- -1
TWSt_MAP_scaler[TWSt_MAP_scaler > 1] <- 1
plot(TWSt_MAP_scaler)

plot(TWSt_MAP_scaler, col = plt, zlim = c(-1, 1))

writeRaster(TWSt_MAP_scaler, here::here("ProducedData", "TWS_MAP_scaler_p10p90_1985_2014map.tif"),
            format = "GTiff", overwrite = T)
