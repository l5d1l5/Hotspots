library(tidyverse); library(magrittr); library(ncdf4); library(here)
library(raster)

setwd("C:/Users/xande/Desktop/Database/WaterUse_Cons")
a <- ncdf4::nc_open("./2d netcdf/withd_min.nc")

# read lat and lon
lon <- ncdf4::ncvar_get(a, "lon")
lat <- ncdf4::ncvar_get(a, "lat")

# assign variable to parameter to extract
a$var$
dname <- "withd_min"

# extract precipitation data
Data <- ncdf4::ncvar_get(a, dname)
dim(Data)

# create empty raster stack to populate
r = stack()
for (j in 1:480) {
  Template.matrix <- matrix(nrow = 360, ncol = 720)
  Template.matrix[] <- NA
  
  for (i in 1:dim(Data)[1]) {
    x = (lon[i] + 180.25)/0.5
    y = (lat[i] + 90.25)/0.5
    d = Data[i, j]
    Template.matrix[y, x] <- d
  }
  
  Template.ras <- raster(Template.matrix)
  Template.ras <- flip(Template.ras, direction = 'y')
  extent(Template.ras) <- c(-180, 180, -90, 90)
  crs(Template.ras) <- CRS("+init=epsg:4326")
  r = stack(r, Template.ras)
  print(j)
}
names(r) <- c(seq(1,480, by = 1))
crs(r) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

writeRaster(r, "./Withdrawal_3d_ncdf/withd_min_3d.nc", 
            overwrite=TRUE, 
            format="CDF", 
            varname= dname, 
            varunit="mm/mo", 
            longname="mm/mo of water consumption", 
            xname="lon", 
            yname="lat")


### - Now calculate the total withdrawals in the year 2010, first by sector

# initalize summing rasters
Dom_2010 <- raster(ext = extent(-180, 180, -90, 90),
                   res=  c(0.5, 0.5),
                   crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
Elec_2010 <- raster(Dom_2010)
Irr_h08_2010 <- raster(Dom_2010)
Irr_lpjml_2010 <- raster(Dom_2010)
Irr_pcr_2010 <- raster(Dom_2010)
Irr_wgap_2010 <- raster(Dom_2010)
Liv_2010 <- raster(Dom_2010)
Mfg_2010 <- raster(Dom_2010)
Min_2010 <- raster(Dom_2010)

Dom_2010[] <- 0
Elec_2010[] <- 0
Irr_h08_2010[] <- 0
Irr_lpjml_2010[] <- 0
Irr_pcr_2010[] <- 0
Irr_wgap_2010[] <- 0
Liv_2010[] <- 0
Mfg_2010[] <- 0
Min_2010[] <- 0

for (i in 1:12) {
  Dom_m <- raster::stack("./Withdrawal_3d_ncdf/withd_dom_3d.nc")[[468+i]]
  Elec_m <- raster::stack("./Withdrawal_3d_ncdf/withd_elec_3d.nc")[[468+i]]
  Irr_h08_m <- raster::stack("./Withdrawal_3d_ncdf/withd_irr_h08_3d.nc")[[468+i]]
  Irr_lpj_m <- raster::stack("./Withdrawal_3d_ncdf/withd_irr_lpjml_3d.nc")[[468+i]]
  Irr_pcr_m <- raster::stack("./Withdrawal_3d_ncdf/withd_irr_pcrglobwb_3d.nc")[[468+i]]
  Irr_wgap_m <- raster::stack("./Withdrawal_3d_ncdf/withd_irr_watergap_3d.nc")[[468+i]]
  Liv_m <- raster::stack("./Withdrawal_3d_ncdf/withd_liv_3d.nc")[[468+i]]
  Mfg_m <- raster::stack("./Withdrawal_3d_ncdf/withd_mfg_3d.nc")[[468+i]]
  Min_m <- raster::stack("./Withdrawal_3d_ncdf/withd_min_3d.nc")[[468+i]]
  
  Dom_m[is.na(Dom_m)] <- 0
  Elec_m[is.na(Elec_m)] <- 0
  Irr_h08_m[is.na(Irr_h08_m)] <- 0
  Irr_lpj_m[is.na(Irr_lpj_m)] <- 0
  Irr_pcr_m[is.na(Irr_pcr_m)] <- 0
  Irr_wgap_m[is.na(Irr_wgap_m)] <- 0
  Liv_m[is.na(Liv_m)] <- 0
  Mfg_m[is.na(Mfg_m)] <- 0
  Min_m[is.na(Min_m)] <- 0
  
  Dom_2010 <- Dom_2010 + Dom_m
  Elec_2010 <- Elec_2010 + Elec_m
  Irr_h08_2010 <- Irr_h08_2010 + Irr_h08_m
  Irr_lpjml_2010 <- Irr_lpjml_2010 + Irr_lpj_m
  Irr_pcr_2010 <- Irr_pcr_2010 + Irr_pcr_m
  Irr_wgap_2010 <- Irr_wgap_2010 + Irr_wgap_m
  Liv_2010 <- Liv_2010 + Liv_m
  Mfg_2010 <- Mfg_2010 + Mfg_m
  Min_2010 <- Min_2010 + Min_m
  
  print(468+i)
}

# Calculate mean irrigation withdrawls from 4 models
Irr_ensm <- stack(Irr_h08_2010, Irr_lpjml_2010, Irr_pcr_2010, Irr_wgap_2010)
Irr_avg <- calc(Irr_ensm, fun = mean, na.rm = T)

# Calculate total water withdrawls per 0.5 grid cell for 2010
TotalWithd_2010 <- raster(Irr_avg)
TotalWithd_2010 <- Dom_2010 + Elec_2010 + Irr_avg + Liv_2010 + Mfg_2010 + Min_2010

writeRaster(TotalWithd_2010,
            "./Withdrawal_3d_ncdf/TotalWithdrawls_2010.tif",
            format = "GTiff",
            overwrite = T)
