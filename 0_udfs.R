################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - All user defined functions used in the project.
################################################################################

# In no particular order ...

# Decompress 2-dimensional netcdf to 3-dimensional netcdf ----
Decompress2dto3dNCDf <- function(Raw2dNCDF, varname, var.unit, CRS.set, WriteLocation) {
  # @Raw2dNCDF: The raw 2-dimensional netcdf 
  # @varname: Variable name of netcdf to extract
  # @var.unit: Variable units
  # @CRS.set: Coordinate reference system to use
  # @WriteLocation: Path to write 3-dimensional netcdf.
  
  # load netcdf and attributes
  TempNCDF <- Raw2dNCDF
  lon <- ncdf4::ncvar_get(TempNCDF, "lon")
  lat <- ncdf4::ncvar_get(TempNCDF, "lat")
  vname <- varname
  message("Data attributes found")
  
  # extract data
  DataExtr <- ncvar_get(TempNCDF, vname)
  message("Data extracted")
  
  # create empty raster stack to populate
  target_stack = stack()
  
  # loop through each month, converting 2d ncdf to 3d
  for (j in 1:dim(DataExtr)[2]) { # repeat for 480 months in dataset
    Template.matrix <- matrix(nrow = 360, ncol = 720)
    Template.matrix[] <- NA
    
    for (i in 1:dim(DataExtr)[1]) { # loop through each entry in 2d ncdf
      x = (lon[i] + 180.25)/0.5
      y = (lat[i] + 90.25)/0.5
      d = DataExtr[i, j]
      Template.matrix[y, x] <- d
    }
    
    Filled.ras <- raster(Template.matrix)
    Filled.ras <- flip(Filled.ras, direction = 'y')
    extent(Filled.ras) <- c(-180, 180, -90, 90)
    crs(Filled.ras) <- CRS("+init=epsg:4326")
    target_stack = stack(target_stack, Filled.ras)
    
    if (j %% 20 == 0) { print(j)}
    
  }
  message("Data decompressed")
  
  names(target_stack) <- sprintf(paste0(varname,"_%d", sep = ""), seq(1:dim(DataExtr)[2]))
  crs(target_stack) <- CRS.set
  
  writeRaster(target_stack, WriteLocation, 
              overwrite=TRUE, 
              format="CDF", 
              varname= vname, 
              varunit= var.unit, 
              xname="lon", 
              yname="lat")
  
  message("3d netcdf stack written")
}

# Area-weighted feature average ----
FeatureAreaAverage <- function(FT.id, RawDS, AreaDS, operation, varnam) {
  # @FT.id: Input raster with ID values unique to each feature
  # @RawDS: Raw raster dataset to calculate area-weighted averages of
  # @AreaDS: Raster with cell values representing the cell area
  # @operation: Specifies summary function (average or sum)
  # @varnam: Variable name to apply to output raster
  
  temp <- raster::stack(FT.id, RawDS, AreaDS) %>% as.data.frame() %>% 
    set_colnames(c("Feature", "RawData", "Area")) 
  
  if (operation == "mean") {
    Sumtab <- temp %>%
      group_by(Feature) %>%
      summarise(
        Stat1 = weighted.mean(x = RawData, w = Area, na.rm = T)
      )  
  }
  
  if (operation == "sum") {
    Sumtab <- temp %>%
      group_by(Feature) %>%
      summarise(
        Stat1 = sum(x = RawData, w = Area, na.rm = T)
      )  
  }
  
  # Format for raster reclassification
  Sumtab$uplim <- Sumtab$Feature+0.5
  Sumtab$lowlim <- Sumtab$Feature-0.5
  rclmtx <- matrix(c(Sumtab$lowlim, Sumtab$uplim, Sumtab$Stat1), ncol = 3)
  
  # Reclassify with summary statistic
  ProdRas <- reclassify(FT.id, rclmtx)
  ProdRas[is.na(FT.id)] <- NA
  names(ProdRas) <- as.character(varnam)
  plot(ProdRas)
  
  # Return raster
  return(ProdRas)
  
}

# Raster percentile reclassification function ----
RasterAreaPercentiles <- function(RasterToClassify, WeightRaster, MaskRaster, clipToExtent, CRS.set, ext.set){
  # @RasterToClassify: Input raster
  # @WeightRaster: Typically representing cell area
  # @MaskRaster: Binary raster indicating cells to mask
  # @clipToExtent: If == "clip", use MaskRaster
  # @CRS.set: Coordinate reference system
  # @ext.set: Optional, if wanting to specify the extent of the output raster
  
  PercentileRaster <- raster(RasterToClassify) # Initialize percentile raster
  
  crs(RasterToClassify) <- CRS.set
  extent(RasterToClassify) <- ext.set
  message(paste0("CRS and extent set to: ", crs(RasterToClassify), "  &  ",
                 extent(RasterToClassify), sep = ""))
  
  RasterToClassify[MaskRaster != 1] <- NA
  
  m.df <- raster::stack(RasterToClassify, WeightRaster, MaskRaster) %>% 
    as.data.frame() %>% 
    set_colnames(c("Input", "Weight", "Mask"))
  m.df <- m.df[complete.cases(m.df$Input),]
  
  if(clipToExtent == "clip"){
    m.df <- m.df %>% dplyr::filter(Mask == 1)
  }
  
  pb <- txtProgressBar(min = 0, max = 99, style = 3)
  
  for(i in 0:99){
    j = 1 - (i*0.01)
    k = 0.99 - (i*0.01)
    
    # Upper bound
    ub = as.numeric(unname(spatstat::weighted.quantile(m.df$Input, 
                                                       m.df$Weight, 
                                                       j, 
                                                       na.rm = TRUE)))
    
    # Lower bound
    lb = as.numeric(unname(spatstat::weighted.quantile(m.df$Input, 
                                                       m.df$Weight, 
                                                       k, 
                                                       na.rm = TRUE)))
    
    PercentileRaster[RasterToClassify <= ub & RasterToClassify > lb] <- j
    setTxtProgressBar(pb, i)
  }
  
  PercentileRaster[is.na(PercentileRaster)] <- 0
  PercentileRaster[is.na(MaskRaster) | MaskRaster != 1] <- NA # mask classified by mask raster
  plot(PercentileRaster)
  
  return(PercentileRaster)
}

# Raster grid buffer function ----
RasterGridBuffer <- function(RastertoBuffer){
  # @RastertoBuffer: Rater to buffer by 1 grid cell, using custom logic shown below
  # Created to enlarge the JPL land mask 
  
  RastertoBuffer[is.na(RastertoBuffer[])] <- 0
  RastertoBuffer[RastertoBuffer[] != 0] <- 1
  
  Input.as.matrix <- as.matrix(RastertoBuffer)
  
  # Buffer JPL land mask by 1 grid cell to make inclusive of other data products 
  Raster.buffer <- matrix(nrow = nrow(Input.as.matrix),
                          ncol = ncol(Input.as.matrix))
  
  for (i in 2:(nrow(as.matrix(Input.as.matrix))-1)) {
    for (j in 2:(ncol(as.matrix(Input.as.matrix))-1)) {
      surrounding_sum = 0
      surrounding_sum <- sum(Input.as.matrix[i+1,j],  
                             Input.as.matrix[i+1,j+1], 
                             Input.as.matrix[i,j+1],
                             Input.as.matrix[i-1,j+1],
                             Input.as.matrix[i-1,j],   
                             Input.as.matrix[i-1,j-1],
                             Input.as.matrix[i,j-1],  
                             Input.as.matrix[i+1,j-1], na.rm = T)
      
      
      if(surrounding_sum != 0 & Input.as.matrix[i,j] == 0){
        Raster.buffer[i,j] <- 1
      } else {
        Raster.buffer[i,j] <- Input.as.matrix[i,j]
      }
    }
  }
  Raster.buffer <- raster(Raster.buffer)
  crs(Raster.buffer) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ 
                           towgs84=0,0,0"
  extent(Raster.buffer) <- extent(RastertoBuffer)
  return(Raster.buffer)
}

# Clip and reproject raster for tmap plotting ----
tmap_clipproj <- function(InputRaster){
  # @InputRaster: Raster to convert
  
  clip.r <- raster::crop(InputRaster, extent(-179, 180, -60, 88))
  rpj.r = projectRaster(clip.r, crs = crs("+proj=robin"), method = 'ngb')
  return(rpj.r)
}

# WGS84 cell area ---- 
WGS84_areaRaster <- function(ResT) {
  # @ResT: Desired resolution
  
  library(tidyverse)
  require(raster)
  
  pi.g <- 3.14159265358979  
  Fl <- 0.00335281066474 # Flattening
  SMA <- 6378137.0 # Semi major axis
  e <- sqrt((2*Fl) - (Fl^2)) # Eccentricity  
  RES <- ResT
  
  # error check entry
  if ((90/RES) %% 1 > 1e-6) { stop("'ResT' must a factor of 90") }
  if (90/RES > (90/(1/24))) { stop("'ResT' is too fine (will require more memory than can be allocated") }
  
  # initialize dataframe with geodetic latitudes
  df_a <- data.frame(LowLAT = seq(-90, 90-RES, by = RES), 
                     UppLAT = seq(-90+RES, 90, by = RES))
  
  # Convert geodetic latitudes degrees to radians
  df_a$LowLATrad <- df_a$LowLAT * pi.g / 180
  df_a$UppLATrad <- df_a$UppLAT * pi.g / 180
  
  # Calculate q1 and q2
  df_a$q1 <- (1-e*e)* ((sin(df_a$LowLATrad)/(1-e*e*sin(df_a$LowLATrad)^2)) - ((1/(2*e))*log((1-e*sin(df_a$LowLATrad))/(1+e*sin(df_a$LowLATrad)))))
  
  df_a$q2 <- (1-e*e)* ((sin(df_a$UppLATrad)/(1-e*e*sin(df_a$UppLATrad)^2)) - ((1/(2*e))*log((1-e*sin(df_a$UppLATrad))/(1+e*sin(df_a$UppLATrad)))))
  
  # calculate q constant
  q_const <- (1-e*e)* ((sin(pi.g/2)/(1-e*e*sin(pi.g/2)^2)) - ((1/(2*e))*log((1-e*sin(pi.g/2))/(1 + e*sin(pi.g/2)))))
  
  # Calculate authaltic latitudes
  df_a$phi1 <- asin(df_a$q1 / q_const)
  df_a$phi2 <- asin(df_a$q2 / q_const)
  
  # Calculate authaltic radius
  R.adius <- sqrt(SMA*SMA*q_const/2)
  
  # Calculate cell size in radians
  CS <- (RES) * pi.g/180
  
  # Calculate cell area in m2
  df_a$area_m2 <- R.adius*R.adius*CS*(sin(df_a$phi2)-sin(df_a$phi1))
  
  # Convert to raster, and replicate column throughout global longitude domain
  WGS84area_km2 <- matrix(df_a$area_m2/1e6, nrow = 180/RES, ncol = 360/RES, 
                          byrow = FALSE, dimnames = NULL) %>% raster()
  extent(WGS84area_km2) <- c(-180, 180, -90, 90) # set extent of raster
  crs(WGS84area_km2) <- crs("+proj=longlat") # set CRS of raster
  
  WGS84area_km2 <- raster::flip(WGS84area_km2, direction = "y")
  
  message(paste0("Calculated global surface area at: ", RES, 
                 "deg. is ", sum(WGS84area_km2[]), " km2.", sep = ""))
  
  return(WGS84area_km2)
}

# Head/tail breaks classification scheme ----
ht_breaks2.0 <- function(x, tsh){
  # @x: vector to classify
  # @tsh: 'head' proportion threshold to end recursive function
  
  allheads <- c()
  
  # Breaks identification
  ht_inner <- function(x, mu){
    n <- length(x)
    mu <- c(mu, mean(x))
    
    # Clip to just 'head' values
    h <- x[x >= mean(x)]
    
    # Threshold check
    headfrac <- length(h)/n
    allheads <- c(allheads, headfrac)
    message(mean(allheads))
    
    # Recursive function 
    if(length(h) > 1 && mean(allheads) <= tsh){
      ht_inner(h, mu)
    } else mu
  }
  
  ht_inner(x, NULL)
}

# Calculate mode of vector ----
getmode <- function(v) {
  # @v: Vector to calculate mode of
  
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Calculate 2nd most frequent value; see above ----
getmode2 <- function(v) {
  uniqv <- unique(v)
  modev <- uniqv[which.max(tabulate(match(v, uniqv)))]
  
  dropmode <- v[v != modev]
  
  if (length(dropmode) > 0) {
    uniqv2 <- unique(dropmode)
    uniqv2[which.max(tabulate(match(dropmode, uniqv2)))]
  } else {
    return(NA)
  }
  
}

### Jenks natural breaks ----
BinMaker <- function(VectorSet, Threshold) {
  # @VectorSet: Vector to identify natural breaks of
  # @Threshold: Goodness of fit variance threshold to stop iteration.
  
  # can't have less than 2 bins
  n_bin = 2
  
  repeat{
    CIobj = classInt::classIntervals(var = VectorSet,
                                     n = n_bin,
                                     style = "jenks")
    
    GoF = classInt::jenks.tests(CIobj)
    
    if(GoF[2] >= Threshold){
      break
    }
    
    n_bin = n_bin + 1
    
  }
  
  return(GoF)
}

# Functions to calculate social-ecological activity statistics for various FW stress and TWS trend combinations ---- 
StressDry <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(tws < -3 & fwstrs > 0.1) %>% pull(col.in) %>% 
    sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

UnstressDry <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(tws < -3 & fwstrs <= 0.1) %>% pull(col.in) %>% 
    sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

StressWet <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(tws > 3 & fwstrs > 0.1) %>% pull(col.in) %>% 
    sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

UnstressWet <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(tws > 3 & fwstrs <= 0.1) %>% pull(col.in) %>% 
    sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}


# Functions to calculate social-ecological activity statistics for various combined freshwater stress indicator and adaptive capacity combinations ---- 

Modstress_badAC <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(cind > ptls$ind[67] & cind <= ptls$ind[80] & ac > ptls$ac[67]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

Modstress_goodAC <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(cind > ptls$ind[67] & cind <= ptls$ind[80] & ac <= ptls$ac[67]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

Highstress_badAC <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(cind > ptls$ind[80] & ac > ptls$ac[67]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

Highstress_goodAC <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(cind > ptls$ind[80] & ac <= ptls$ac[67]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}

# Functions to calculate social-ecological activity statistics for vulnerability classes ----

ModerateVul <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(overallprod >= htb_o[1] & overallprod < htb_o[2]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}
HighVul <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(overallprod >= htb_o[2] & overallprod < htb_o[3]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}
VeryHighVul <- function(df.in, col.in, norm.by) {
  glob.sum <- df.in %>% pull(col.in) %>% sum(na.rm = T)
  
  tot <- df.in %>% filter(overallprod >= htb_o[3]) %>% 
    pull(col.in) %>% sum(na.rm = T)/norm.by
  
  pct <- 100 * (tot*norm.by/glob.sum)
  
  message(paste0("Total is: ", round(tot, 2), "by", norm.by, sep = ""))
  message(paste0("Pct is: ", round(pct, 1), sep = ""))
  
}