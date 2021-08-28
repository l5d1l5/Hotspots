################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - Preprocess HydroBASINS discretization schemes for use. 
################################################################################

# Load 'here' library for easy path management.  
library(here)

# Import all setup and user-defined functions in R/setup and R/udfs folders
invisible(sapply(paste0(here("R/setup"), "/", list.files(here("R/setup"))), source)) 
invisible(sapply(paste0(here("R/udfs"), "/", list.files(here("R/udfs"))), source))

# NOTE: While only the process for Level 3 is shown below, it is repeated for Levels 4 and 5 by only changing the imported shapefiles and write locations accordingly, and no other changes are made to the code.

# Import all world regions
af <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/hybas_af_lev03_v1c.shp')
ar <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/hybas_ar_lev03_v1c.shp')
as <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/hybas_as_lev03_v1c.shp')
au <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/hybas_au_lev03_v1c.shp')
eu <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/hybas_eu_lev03_v1c.shp')
gr <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/hybas_gr_lev03_v1c.shp')
na <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/hybas_na_lev03_v1c.shp')
sa <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/hybas_sa_lev03_v1c.shp')
si <- sf::read_sf('D:/!! Geodatabase/Basins/hybas/hybas_si_lev03_v1c.shp')

# Merge into single file
hybas_list <- list(af, ar, as, au, eu, gr, na, sa, si)
hybas_l_set   <- do.call(rbind, hybas_list) 

# Write combined shapefile
write_sf(obj = hybas_l_set, dsn = here('Data/'), layer = 'hybas_l3.shp', 
         driver = 'ESRI Shapefile')

# Rasterize, using finer resolution then aggregating to avoid clipping
hybas_ras <- fasterize(sf = hybas_l_set, raster = WGS84_areaRaster(0.05),  
                          field = 'PFAF_ID')

# Aggregate using modal value, ignoring NAs
hybas_rasag <- aggregate(hybas_ras, fact = 10,
                          expand = FALSE, fun = modal, na.rm = T)

writeRaster(hybas_rasag, here('Data/hybas_l3_pfafid.tif'),
            format = 'GTiff', overwrite = T)

# Convert raster to shapefile, for plotting
hybas_0d5 <- rasterToPolygons(hybas_rasag, dissolve = T)
hybas_0d5 <- sf::st_as_sf(hybas_rasag)
write_sf(obj = hybas_0d5, dsn = here('Data/'), layer = 'hybas_l3_0d5.shp', 
         driver = 'ESRI Shapefile')

# Clean basins for analysis, by masking regions with partial input data coverage, removing Greenland

# GRUN has limited data coverage in high latitudes, so use as extent limiting dataset
GRUN <- raster(here("Data", "GRUN", "MeanGRUN.tif"))

# Create binary raster indicating where GRUN exists
GRUNexist <- raster(GRUN)
GRUNexist[] <- 0
GRUNexist[GRUN >= 0] <- 1

# Calculate GRUN coverage per basin
GRUN_coverage <- FeatureAreaAverage(FT.id = hybas_rasag,
                                    RawDS = GRUNexist,
                                    AreaDS = WGS84_areaRaster(0.5),
                                    operation = "mean",
                                    varnam = "GRUNcoverage")

# Exclude basins with less than 67% area GRUN coverage
GRUN_clip <- GRUN_coverage
GRUN_clip[GRUN_clip < 0.67] <- NA
GRUN_clip[GRUN_clip > 0] <- 1
plot(GRUN_clip)

# Manually remove other problematic basins

# All of Greenland
manualclip <- raster(hybas_rasag)
manualclip[] <- 1
manualclip[hybas_rasag >= 900] <- NA # Greenland PFAF ids

# Combine clipping extents
ClipExtent <- GRUN_clip * manualclip

# Clip the hydroBASIN ID raster
hybas_rasag_cleaned <- hybas_rasag*ClipExtent


# Write cleaned raster
writeRaster(hybas_rasag_cleaned, here::here("Data", "HyBas3_cleaned.tif"), 
            format = "GTiff", overwrite = T)