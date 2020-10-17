library(tidyverse)
library(magrittr)
library(ncdf4)
library(here)
library(raster)

### - Now calculate the total withdrawals in the year 2010
### - for the non-agricultural economic sectors: (energy, livestock, mining, manufacturing)
# initalize summing rasters
Elec_2010 <- raster(ext = extent(-180, 180, -90, 90),
                   res=  c(0.5, 0.5),
                   crs = "+proj=longlat")
Liv_2010 <- raster(Elec_2010)
Mfg_2010 <- raster(Elec_2010)
Min_2010 <- raster(Elec_2010)

Elec_2010[] <- 0
Liv_2010[] <- 0
Mfg_2010[] <- 0
Min_2010[] <- 0

setwd("C:/Users/xande/Desktop/Database/Huang_WaterUse")
for (i in 1:12) {
  Elec_m <- raster::stack("./Withdrawal_3d_ncdf/withd_elec_3d.nc")[[468+i]]
  Liv_m <- raster::stack("./Withdrawal_3d_ncdf/withd_liv_3d.nc")[[468+i]]
  Mfg_m <- raster::stack("./Withdrawal_3d_ncdf/withd_mfg_3d.nc")[[468+i]]
  Min_m <- raster::stack("./Withdrawal_3d_ncdf/withd_min_3d.nc")[[468+i]]
  
  crs(Elec_m) <- "+proj=longlat"
  crs(Mfg_m) <- "+proj=longlat"
  crs(Elec_m) <- "+proj=longlat"
  crs(Min_m) <- "+proj=longlat"
  
  Elec_m[is.na(Elec_m)] <- 0
  Liv_m[is.na(Liv_m)] <- 0
  Mfg_m[is.na(Mfg_m)] <- 0
  Min_m[is.na(Min_m)] <- 0
  
  Elec_2010 <- Elec_2010 + Elec_m
  Liv_2010 <- Liv_2010 + Liv_m
  Mfg_2010 <- Mfg_2010 + Mfg_m
  Min_2010 <- Min_2010 + Min_m
  
  print(468+i)
}


# Calculate total water withdrawls per 0.5 grid cell for 2010
TotalWithd_2010 <- raster(Elec_m)
TotalWithd_2010 <- Elec_2010 + Liv_2010 + Mfg_2010 + Min_2010

writeRaster(TotalWithd_2010,
            "./Withdrawal_3d_ncdf/EconSectors_2010Withdrawals_sum.tif",
            format = "GTiff",
            overwrite = T)
