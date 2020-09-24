## this script harmonizes dimension analysis data

library(raster); library(gdalUtils); library(ncdf4); library(sf)
library(fasterize); library(here)

# load each dataset to see resolution
######pop
Pop <- raster(here::here("Data", "Dimensions",
    "gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif"))
# Pop is already in 0.5d resolution

######cal
Cal <- raster(here::here("Data", "Dimensions", "Glbkcal.tif"))
x = 0.5/res(Cal)
Cal_0d5 <- raster::aggregate(x = Cal, fact = x, fun = sum, expand = F, na.rm = T, 
                             filename = here::here("Data", "Dimensions", "kcal_0d5.tif"))

###### GDP
GDP <- nc_open(here::here("Data", "Dimensions", "GDP_PPP_1990_2015_5arcmin_v2.nc"))
lon = ncvar_get(GDP, "longitude")
lat = ncvar_get(GDP, "latitude")

v.name = "GDP_PPP"

# see metadata
ncatt_get(GDP, v.name, "long_name")
ncatt_get(GDP, v.name, "units")
fv = ncatt_get(GDP, v.name, "_FillValue")

# extract precip data
GDP_years <- ncvar_get(GDP, v.name)
dim(GDP_years)

# 2015 is year # 26 (last year) in dataset
GDP_2015 <- raster(t(GDP_years[,,26]), xmn=min(lon), xmx=max(lon), 
                   ymn=min(lat), ymx=max(lat), 
                   crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
extent(GDP_2015) <- c(-180, 180, -90, 90)
x = 0.5/res(GDP_2015)
GDP_2015_0d5 <- raster::aggregate(x = GDP_2015, fact = x, fun = sum, expand = F, na.rm = T, 
                             filename = here::here("Data", "Dimensions", 
                                                   "GDP_2015_0d5.tif"))

###### Biodiversity (1) amphibians, (2) mammals, (3) ecoregions, (4) VSI, (5) Eflows
ga <- raster(here::here("ProducedData", "WGS_gridArea_halfdegree.tif"))

###(1) amphibians
amph <- raster(here::here("Data", "Dimensions", "all_amphibians.tif"))
extent(amph) = c(-180, 180, -90, 90)
x = 0.5/res(amph)
# use maximum TD value per grid cell
# use maximum TD value per grid cell... because high res, use gdalwarp instead
gdalwarp(srcfile = here::here("Data", "Dimensions", "all_amphibians.tif"),
         dstfile = here::here("Data", "Dimensions", "amph_0d5.tif"),
         te = c(-180, -90, 180, 90),
         tr = c(0.5, 0.5),
         r = "max",
         output_Raster = TRUE,
         overwrite = TRUE,
         verbose = TRUE)

###(2) mammals
mamm <- raster(here::here("Data", "Dimensions", "all_mammals.tif"))
extent(mamm) = c(-180, 180, -90, 90)
x = 0.5/res(mamm)
# use maximum TD value per grid cell... because high res, use gdalwarp instead
gdalwarp(srcfile = here::here("Data", "Dimensions", "all_mammals.tif"),
         dstfile = here::here("Data", "Dimensions", "mamm_0d5.tif"),
         te = c(-180, -90, 180, 90),
         tr = c(0.5, 0.5),
         r = "max",
         output_Raster = TRUE,
         overwrite = TRUE,
         verbose = TRUE)

###(3) Global200
g200.fw <- sf::read_sf(here::here("Data", "Dimensions", "G200", "g200_fw.shp"))
g200.tr <- sf::read_sf(here::here("Data", "Dimensions", "G200", "g200_terr.shp"))

g200.fw.r <- fasterize(sf = g200.fw, raster = ga, fun = "last")
g200.fw.r[is.na(g200.fw.r)] <- 0
g200.tr.r <- fasterize(sf = g200.tr, raster = ga, fun = "last")
g200.tr.r[is.na(g200.tr.r)] <- 0

g200.s <- raster::stack(g200.fw.r, g200.tr.r)
g200.r <- raster::calc(g200.s, fun = max, na.rm = T)
writeRaster(g200.r, filename = here::here("Data", "Dimensions", "g200_fw_tr_0d5.tif"),
            format = "GTiff", overwrite = T)

###(4) VSI
VSI <- raster(here::here("Data", "Dimensions", "SensAETW.tif"))
gdalwarp(srcfile = here::here("Data", "Dimensions", "SensAETW.tif"),
         dstfile = here::here("Data", "Dimensions", "VSI_aetw_0d5.tif"),
         te = c(-180, -90, 180, 90),
         tr = c(0.5, 0.5),
         r = "average",
         output_Raster = TRUE,
         overwrite = TRUE,
         verbose = TRUE)                                     
plot(raster(here::here("Data", "Dimensions", "VSI_aetw_0d5.tif")))

###(5) EFNs
EFN.dg <- raster(here::here("Data", "Dimensions", "headdrop2limit_hydr06.map"))
x = 0.5/res(EFN.dg)
EFN.dg_0d5 <- raster::aggregate(x = EFN.dg, fact = x, fun = mean, expand = F, na.rm = T, 
                                  filename = here::here("Data", "Dimensions", 
                                                        "EFN_0d5.tif"),
                                overwrite = T)
