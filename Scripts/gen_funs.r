require(spatstat)

################## 1 - Raster percentile function
RasterAreaPercentiles <- function(RasterToClassify, WeightRaster, MaskRaster, clipToExtent){
  PercentileRaster <- raster(RasterToClassify) # initalize a percentile raster
  # PercentileRaster[] <- 0 # set inital values to 0 
  
  m.df <- cbind(RasterToClassify[], WeightRaster[], MaskRaster[]) %>% 
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
    
    # upper bound
    ub = as.numeric(unname(spatstat::weighted.quantile(m.df$Input, m.df$Weight, j, na.rm = TRUE)))
    lb = as.numeric(unname(spatstat::weighted.quantile(m.df$Input, m.df$Weight, k, na.rm = TRUE)))
    PercentileRaster[RasterToClassify <= ub & RasterToClassify > lb] <- j
    setTxtProgressBar(pb, i)
  }
  PercentileRaster[is.na(PercentileRaster)] <- 0
  PercentileRaster[MaskRaster != 1] <- NA # mask classified by mask raster
  return(PercentileRaster)
}

################## 2 - Raster grid buffer function
RasterGridBuffer <- function(RastertoBuffer){
  RastertoBuffer[is.na(RastertoBuffer[])] <- 0
  RastertoBuffer[RastertoBuffer[] != 0] <- 1
  
  Input.as.matrix <- as.matrix(RastertoBuffer)
  
  # Buffer JPL land mask by 1 grid cell to make inclusive of other data products 
  Raster.buffer <- matrix(nrow = nrow(Input.as.matrix),
                        ncol = ncol(Input.as.matrix))
  
  for (i in 2:(nrow(as.matrix(Input.as.matrix))-1)) {
    for (j in 2:(ncol(as.matrix(Input.as.matrix))-1)) {
      surrounding_sum = 0
      surrounding_sum <- sum(Input.as.matrix[i+1,j],  Input.as.matrix[i+1,j+1], Input.as.matrix[i,j+1],
                             Input.as.matrix[i-1,j+1],Input.as.matrix[i-1,j],   Input.as.matrix[i-1,j-1],
                             Input.as.matrix[i,j-1],  Input.as.matrix[i+1,j-1], na.rm = T)
      
      
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

################## 3 - tmap raster clip and reproject function
tmap_clipproj <- function(InputRaster){
  clip.r <- raster::crop(InputRaster, extent(-179, 180, -60, 88))
  rpj.r = projectRaster(clip.r, crs = crs("+proj=robin"))
  return(rpj.r)
}
 
################## 4 - Hotspot filter to (1) fill holes, (2) remove islands
hotspot_id_smoother <- function(RawRaster, WeightRaster, MaskRaster, Method, Multiplier){
  
  library(Weighted.Desc.Stat); library(Hmisc)
  
  Int.ras <- raster(WeightRaster)
  ID.ras <- raster(WeightRaster)
  ID.ras[] <- 0
  
  df.temp <- cbind(RawRaster[], WeightRaster[]) %>% as.data.frame() %>%
    set_colnames(c("Value", "Area"))
  df.temp <- df.temp[complete.cases(df.temp),]
  
  if (Method == "IQR") {
    calc.IQR <- spatstat::weighted.quantile(x = df.temp$Value, w = df.temp$Area,
                                            probs = 0.75, na.rm = T) -
      spatstat::weighted.quantile(x = df.temp$Value, w = df.temp$Area,
                                  probs = 0.25, na.rm = T)
    calc.IQR < abs(calc.IQR)
    
    Int.ras <- RawRaster/(calc.IQR*Multiplier)
    ID.ras[abs(Int.ras[]) >= 1] <- 1
  }
  
  if (Method == "IDR") {
    calc.IDR <- spatstat::weighted.quantile(x = df.temp$Value, w = df.temp$Area,
                                            probs = 0.9, na.rm = T) -
      spatstat::weighted.quantile(x = df.temp$Value, w = df.temp$Area,
                                  probs = 0.1, na.rm = T)
    calc.IDR < abs(calc.IDR)
    
    Int.ras <- RawRaster/(calc.IDR*Multiplier)
    ID.ras[abs(Int.ras[]) >= 1] <- 1
  }
  
  if (Method == "SD") {
    calc.sd <- Weighted.Desc.Stat::w.sd(x = df.temp$Value, mu = df.temp$Area)
    
    Int.ras <- RawRaster/(calc.sd*Multiplier)
    ID.ras[abs(Int.ras[]) >= 1] <- 1
  }
  
  Matrix.1 <- as.matrix(ID.ras)
  
  Smooth.1 <- matrix(nrow = nrow(Matrix.1),
                           ncol = ncol(Matrix.1))
  
  for (i in 2:(nrow(as.matrix(Matrix.1))-1)) {
    for (j in 2:(ncol(as.matrix(Matrix.1))-1)) {
      
      # i = 150
      # j = 150
      sur_sum <- sum(Matrix.1[i+1,j],   Matrix.1[i+1,j+1], Matrix.1[i,j+1],
                     Matrix.1[i-1,j+1], Matrix.1[i-1,j],   Matrix.1[i-1,j-1],
                     Matrix.1[i,j-1],   Matrix.1[i+1,j-1], na.rm = T)
      
      if(is.na(Matrix.1[i,j])){
        Smooth.1[i,j] <- NA
      } else {
        if(Matrix.1[i,j] == 1 & sur_sum < 2){ # if hotspot but not touching 2 other hotspot cells
          Smooth.1[i,j] <- 0              # smooth to make not a hotshot 
        }
        
        if(Matrix.1[i,j] == 1 & sur_sum >= 2){# if hotspot and touching min. 2 other hotspot cells 
          Smooth.1[i,j] <- 1              # keep as hotspot     
        }
        
        if(Matrix.1[i,j] == 0 & sur_sum > 3){# if not a hotspot but touching min. 3 other hotspot cells 
          Smooth.1[i,j] <- 1            # make a hotspot (i.e. fill gaps)  
        }
        
        if(Matrix.1[i,j] == 0 & sur_sum <= 3){# if not a hotspot and not touching others
          Smooth.1[i,j] <- 0              # keep as not a hotspot
        }
      }
    }
  }
  
  Smooth.1 <- raster(Smooth.1)
  crs(Smooth.1) <- "+proj=longlat"
  extent(Smooth.1) <- c(-180, 180, -90, 90)
  Smooth.1[1,] <- NA
  Smooth.1[,1] <- NA
  
  Matrix.2 <- as.matrix(Smooth.1)
  
  Smooth.2 <- matrix(nrow = nrow(Matrix.2),
                             ncol = ncol(Matrix.2))
  
  for (i in 2:(nrow(as.matrix(Matrix.2))-1)) {
    for (j in 2:(ncol(as.matrix(Matrix.2))-1)) {
      
      # i = 150
      # j = 150
      sur_sum <- sum(Matrix.2[i+1,j],   Matrix.2[i+1,j+1], Matrix.2[i,j+1],
                     Matrix.2[i-1,j+1], Matrix.2[i-1,j],   Matrix.2[i-1,j-1],
                     Matrix.2[i,j-1],   Matrix.2[i+1,j-1], na.rm = T)
      
      
      if(is.na(Matrix.2[i,j] | is.na(sur_sum))){
        Smooth.2[i,j] <- NA
      } else {
        if(Matrix.2[i,j] == 1 & sur_sum < 1) { # remove any "island" cells 
          Smooth.2[i,j] <- 0
        } else {
          Smooth.2[i,j] <- Matrix.2[i,j]
        }
      }
    }
  }
  
  
  Smooth.2 <- raster(Smooth.2)
  crs(Smooth.2) <- "+proj=longlat"
  extent(Smooth.2) <- c(-180, 180, -90, 90)
  Smooth.2[1,] <- NA
  Smooth.2[,1] <- NA
  
  Smooth.2[MaskRaster != 1] <- NA
  
  return(Smooth.2)
}