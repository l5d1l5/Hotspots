## this script derives mean agricultural intensification for 
## wheat, rice, and Rice from Mueller et al. 2012 data on EarthStat

library(tidyverse); library(magrittr); library(raster); library(here)

## Maize
Maize_gap <- raster(here::here("Data", "YieldGaps", "maize_yieldgap_geotiff",
                              "maize_yieldgap.tif"))
Maize_pot <- raster(here::here("Data", "YieldGaps", "maize_yieldgap_geotiff",
                              "maize_yieldpotential.tif"))
Maize_frac <- Maize_gap/Maize_pot

## Rice
Rice_gap <- raster(here::here("Data", "YieldGaps", "rice_yieldgap_geotiff",
                               "rice_yieldgap.tif"))
Rice_pot <- raster(here::here("Data", "YieldGaps", "rice_yieldgap_geotiff",
                               "rice_yieldpotential.tif"))
Rice_frac <- Rice_gap/Rice_pot

## Wheat
Wheat_gap <- raster(here::here("Data", "YieldGaps", "wheat_yieldgap_geotiff",
                              "wheat_yieldgap1.tif"))
Wheat_pot <- raster(here::here("Data", "YieldGaps", "wheat_yieldgap_geotiff",
                              "wheat_yieldpotential.tif"))
Wheat_frac <- Wheat_gap/Wheat_pot

## calculate average yield gap fraction for three crops
Grains_stack <- stack(Maize_frac, Rice_frac, Wheat_frac)
Grains_mean <- calc(Grains_stack, mean, na.rm = T)
# Grains_mean[is.na(Grains_mean)] <- 0

## import grid area raster at half degree for resampling
ga <- raster(here::here("ProducedData", "WGS_gridArea_halfdegree.tif"))

Grains_mean_0d5 <- raster::resample(x = Grains_mean, y = ga, method = "bilinear")

writeRaster(Grains_mean_0d5, here::here("ProducedData", "Grains_YieldGap.tif"),
            format = "GTiff", overwrite = T)