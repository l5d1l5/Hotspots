library(tidyverse); library(magrittr); library(raster)

Rodell_raw <- readr::read_csv(here::here("Data", "41586_2018_123_MOESM1_ESM.csv"),
                              col_names = F) %>%
  as.matrix()

TWSt.ras <- raster(Rodell_raw)
extent(TWSt.ras) <- c(-180, 180, -90, 90)
crs(TWSt.ras) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"
TWSt.ras <- flip(TWSt.ras, direction = "y")

writeRaster(TWSt.ras, here::here("Data", "Rodell_GRACE_TWSt_raw.tif"),
            format = "GTiff", overwrite = T)
