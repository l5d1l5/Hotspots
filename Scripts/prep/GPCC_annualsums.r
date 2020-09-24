## this script produces a raster containing the long-term (30-year) mean annual
## precipitation at 0.5 degrees, using 1972-2001 GPCC monthly data

library(tidyverse); library(magrittr); library(raster); library(ncdf4)
library(here)

# import GPCC net cdf
GPCC.ncdf <- nc_open(here::here("Data", "precip.mon.total.v2018.nc"))
lon = ncvar_get(GPCC.ncdf, "lon")
lat = ncvar_get(GPCC.ncdf, "lat")

v.name = "precip"

# see metadata
ncatt_get(GPCC.ncdf, v.name, "long_name")
ncatt_get(GPCC.ncdf, v.name, "units")
fv = ncatt_get(GPCC.ncdf, v.name, "_FillValue")

# extract precip data
GPCC.precip <- ncvar_get(GPCC.ncdf, v.name)
dim(GPCC.precip)

# run loop to calculate annual precipitation for entire record
for(i in 1:(1512/12)){
  j = 1 + ((i-1)*12)
  
  # extract each month in calendar year
  m1 <- raster(t(GPCC.precip[,,j]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m2 <- raster(t(GPCC.precip[,,j+1]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m3 <- raster(t(GPCC.precip[,,j+2]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m4 <- raster(t(GPCC.precip[,,j+3]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m5 <- raster(t(GPCC.precip[,,j+4]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m6 <- raster(t(GPCC.precip[,,j+5]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m7 <- raster(t(GPCC.precip[,,j+6]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m8 <- raster(t(GPCC.precip[,,j+7]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m9 <- raster(t(GPCC.precip[,,j+8]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m10 <- raster(t(GPCC.precip[,,j+9]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m11 <- raster(t(GPCC.precip[,,j+10]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  m12 <- raster(t(GPCC.precip[,,j+11]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  
  # sum month sums to year
  yr <- m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12 
  
  # CS is 0-360, not -180-180, thus needs to be re-arranged
  Arrange <- raster(yr)
  Arrange[,1:360] <- yr[,361:720]
  Arrange[,361:720] <- yr[,1:360]
  
  # ensure correct extent
  extent(Arrange) <- c(-180, 180, -90, 90)
  
  # write each year as .tif file
  writeRaster(Arrange, paste(here::here("Data", "Precip_AnnualSums"), "/yr_", 1890+i, ".tif", sep=""), 
              format = "GTiff",overwrite = T)
}

# call list of years and select only years 1972 through 2001 (30 years pre GRACE)
List_years = list.files(path = here::here("Data", "Precip_AnnualSums"), 
                                          pattern = "*.tif")
List_years = List_years[(1985+1-1891):(2014+1-1891)]

Sum_ras <- raster(ext = extent(-180, 180, -90, 90), res = c(0.5, 0.5), 
                  crs = CRS("+proj=longlat"))
Sum_ras[] <- 0

for (i in 1:length(List_years)) {
  tmp.Ras = raster(paste(here::here("Data", "Precip_AnnualSums"), "/", List_years[i], sep = ""))
  tmp.Ras[is.na(tmp.Ras)] <- 0
  Sum_ras <- Sum_ras + tmp.Ras
  print(1984+i)
}

Sum_ras <- Sum_ras/length(List_years) # get average annual rainfall

writeRaster(Sum_ras, here::here("ProducedData", "MeanAnnualPrecip_1985_2014.tif"),
                                format = "GTiff", overwrite = T)
            