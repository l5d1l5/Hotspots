## this script creates a mask for all analysis that:
## (1) masks to global land extend
## (2) removes Antarctica (all land south of 60*S)
## (3) Removes earthquake interference signals in Indonesia and Japan

library(tidyverse); library(magrittr); library(here)
library(raster); library(sf); library(ncdf4); library(fasterize)
source(here::here("Scripts", "gen_funs.r"))

Land_mask <- nc_open(here::here("Data", "LAND_MASK.CRI.nc"))
lon <- ncvar_get(Land_mask, "lon")
lat <- ncvar_get(Land_mask, "lat")
dname <- "land_mask"
Land_mask_get <- ncvar_get(Land_mask, dname)
Land_RAS <- raster(t(Land_mask_get[,]), xmn=min(lon), xmx=max(lon), 
                   ymn=min(lat), ymx=max(lat), 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ 
                           towgs84=0,0,0"))
Land_RAS <- flip(Land_RAS, direction = "y") # data is flipped in y direction
Mask <- raster(Land_RAS)
Mask[,1:360] <- Land_RAS[,361:720]
Mask[,361:720] <- Land_RAS[,1:360]
extent(Mask) <- c(-180, 180, -90, 90)

Mask_buffer <- RasterGridBuffer(Mask)

plot(Mask)
plot(Mask_buffer)

## Now import EQ interference
EQ_int <- sf::read_sf(here::here("Data", "EQ_int", "EQ_int.shp"))
EQ_int.r <- fasterize(sf = EQ_int, raster = Mask_buffer, fun = "last")

EQ_buffer <- RasterGridBuffer(EQ_int.r)

FinalMask <- Mask_buffer
FinalMask[EQ_buffer == 1] <- 0
FinalMask[(300:360),] <- 0
plot(FinalMask)
FinalMask[1, ] <- 0
FinalMask[, 1] <- 0

writeRaster(FinalMask, here::here("ProducedData", "GlobalMask"),
            format = "GTiff", overwrite = T)