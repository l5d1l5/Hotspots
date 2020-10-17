library(shiny)
library(leaflet)
library(tidyverse)
library(magrittr)
library(sf)
library(proj4)
library(raster)
library(readxl)
library(leaflet.esri)
library(wesanderson)
library(RColorBrewer)
library(googlesheets4)
library(rsconnect)
library(htmltools)
library(rgdal)
library(rgeos)
library(shinyjs)
library(shinyWidgets)
library(spatstat)
library(tmap)

## Comment out line below before publishing to shinyapps.io
# setwd("C:/Users/xande/Desktop/Scripts/twst_ses_project/Scripts/Shiny")
# 
# source functions
##### 

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

################## 5 - Hotspot legend function
addLegendCustom1 <- function(map, colors, labels, sizes, shapes, borders, location, opacity = 1.0){
  
  make_shapes <- function(colors, sizes, borders, shapes) {
    shapes <- gsub("square", "0%", shapes)
    paste0(colors, "; width:", sizes, "px; height:", sizes, "px; border:3px solid ", borders, "; border-radius:", shapes)
  }
  make_labels <- function(sizes, labels) {
    paste0("<div style='display: inline-block;height: ", 
           sizes, "px;margin-top: 0px;line-height: ", 
           sizes, "px;'>", labels, "</div>")
  }
  
  legend_colors <- make_shapes(colors, sizes, borders, shapes)
  legend_labels <- make_labels(sizes, labels)
  
  return(addLegend(map, colors = legend_colors, labels = legend_labels, 
                   position = location, opacity = opacity))
}

################## 6 - WGS area grid
WGS84_areaRaster <- function(ResT) {
  library(tidyverse)
  require(raster)
  
  pi.g <- 3.14159265358979  
  Fl <- 0.0033528106811823171 # flattening
  SMA <- 6378137.0 # semi major axis
  e <- sqrt((2*Fl) - (Fl^2)) # eccentricity  
  RES <- ResT
  
  # error check entry
  if ((90/RES) %% 1 > 1e-6) { stop("'ResT' must a factor of 90") }
  if (90/RES > (90/(1/24))) { stop("'ResT' is too fine (will require more memory than can be allocated") }
  
  # initialize dataframe with geodetic latitudes
  df.a <- data.frame(LowLAT = seq(-90, 90-RES, by = RES), 
                     UppLAT = seq(-90+RES, 90, by = RES))
  
  # Convert geodetic latitudes degrees to radians
  df.a$LowLATrad <- df.a$LowLAT * pi.g / 180
  df.a$UppLATrad <- df.a$UppLAT * pi.g / 180
  
  # Calculate q1 and q2
  df.a$q1 <- (1-e*e)* ((sin(df.a$LowLATrad)/(1-e*e*sin(df.a$LowLATrad)^2)) - ((1/(2*e))*log((1-e*sin(df.a$LowLATrad))/(1+e*sin(df.a$LowLATrad)))))
  df.a$q2 <- (1-e*e)* ((sin(df.a$UppLATrad)/(1-e*e*sin(df.a$UppLATrad)^2)) - ((1/(2*e))*log((1-e*sin(df.a$UppLATrad))/(1+e*sin(df.a$UppLATrad)))))
  
  # calculate q constant
  q.const <- (1-e*e)* ((sin(pi.g/2)/(1-e*e*sin(pi.g/2)^2)) - ((1/(2*e))*log((1-e*sin(pi.g/2))/(1+e*sin(pi.g/2)))))
  
  # Calculate authaltic latitudes
  df.a$phi1 <- asin(df.a$q1 / q.const)
  df.a$phi2 <- asin(df.a$q2 / q.const)
  
  # Calculate authaltic radius
  R.adius <- sqrt(SMA*SMA*q.const/2)
  
  # Calculate cell size in radians
  CS <- (RES) * pi.g/180
  
  # Calculate cell area in m2
  df.a$area_m2 <- R.adius*R.adius*CS*(sin(df.a$phi2)-sin(df.a$phi1))
  
  # Convert to raster, and replicate column throughout global longitude domain
  WGS84area_km2 <- matrix(df.a$area_m2/1e6, nrow = 180/RES, ncol = 360/RES, 
                          byrow = FALSE, dimnames = NULL) %>% raster()
  extent(WGS84area_km2) <- c(-180, 180, -90, 90) # set extent of raster
  crs(WGS84area_km2) <- crs("+proj=longlat") # set CRS of raster
  WGS84area_km2 <- raster::flip(WGS84area_km2, direction = "y") # flip raster (inconsequential, but technically correct)
  message(paste0("Calculated global surface area at: ", RES, "deg. is ", sum(WGS84area_km2[]), " km2.", sep = ""))
  return(WGS84area_km2)
}

#####

# Import individual dimensions
Pop <- raster("./Pop_indicator.tif")
Cal <- raster("./Cal_indicator.tif")
GDP <- raster("./Econ_indicator.tif")
Eco <- raster("./Eco_sens_indicator.tif")
ga <- WGS84_areaRaster(0.5)
Mask <- raster("./GlobalMask.tif")
Trends <- raster("./Rodell_GRACE_TWSt_raw.tif")

data(World)
World$id <- c(seq(1, nrow(World), by = 1))
World <- sf::st_as_sf(World)
World84 <- st_transform(World, 4326)

OverallComposite <- raster("./DefaultComposite.tif")
Hotspots_all.sf <- sf::read_sf("./DefaultHotspots.shp")

# initialize temp directory for session and remove existing files
temp_dir <- tmpDir()

## end temp directory initialization

df <- data.frame(Dimension = c("Population", "Crop calories", "GDP (PPP)", 
                               "Water sensitive environ."), 
                 `In Drying Hotspots` = c(803.50, 879.54, 10.43, 1.23),
                 `In Wetting Hotspots` = c(54.92, 215.08, 0.82, 0.57),
                 Units = c("million people", "trillion calories",
                           "trillion USD (2011)", "million sq. km"))
names(df) <- c("Dimension", "In Drying Hotspots", "In Wetting Hotspots",
               "Units/Magnitude")

############ Start making shiny app below ... 
plt = pals::ocean.balance(20)[3:18] %>% rev()

pal <- colorNumeric(palette = plt, domain = c(-1,1), na.color = "transparent")
palr <- colorNumeric(palette = rev(plt), domain = c(-1,1), na.color = "transparent")

ui <- bootstrapPage(
  shinyjs::useShinyjs(),
  # Formatting 
  tags$head(
    # tags$link(href = "https://fonts.googleapis.com/css?family=Neuton", rel = "stylesheet"),
    tags$style(type = "text/css", "html, body {width:100%;height:100%; color: black;}"),
    tags$style(type = 'text/css', '.well {background-color: #00244a;}'),
    tags$style(HTML('#panel {background-color: grey;}')),
    tags$head(tags$style(
      "#sumCheck {
      color: red; font-size: 12px; font-weight: bold;
      }",
      ".leaflet .legend {
                 line-height: 13px;
                 font-size: 13px;
                 }"))),
  
  # Leaflet map
  
  leafletOutput("mymap", width = "100%", height = "100%"),
  
  # input toolbar to adjust map
  absolutePanel(
    top = 50, left = 10, draggable = FALSE, width = "15%", style = "z-index:500; min-width: 200px;",
    id = "panel", style='background-color: rgba(255, 255, 255,  0.95); padding: 5px;',
    actionButton(inputId = "Recalc", label = "Update map", class = "btn-primary",
                 style = "padding-bottom: 5px;", width = "100%"),
    p(), p(),
    progressBar(id = "pb4", value = 0, display_pct = T),
    hr(),
    sliderInput(inputId = "Pop.w", 
                label = "Population Weight", 
                min = 0, max = 1, step = 0.01, 
                value = 0.22),
    sliderInput(inputId = "Ag.w", 
                label = "Agricultural Weight", 
                min = 0, max = 1, step = 0.01, 
                value = 0.22),
    sliderInput(inputId = "Econ.w", 
                label = "Economic Weight", 
                min = 0, max = 1, step = 0.01, 
                value = 0.11),
    sliderInput(inputId = "Env.w", 
                label = "Environmental Weight", 
                min = 0, max = 1, step = 0.01, 
                value = 0.45),
    textOutput('sumCheck'), 
    p(), p(),
    selectInput("Method", "Statistic to use:",
                c("Standard deviation" = "SD",
                  "IQR (p75-p25)" = "IQR",
                  "IDR (p90-p10)" = "IDR")),
    sliderInput(inputId = "Mult", 
                label = "Statistic multiplier", 
                min = 1, max = 3, step = 0.05, 
                value = 2),
    p(),p(),
    sliderInput(inputId = "opac", 
                label = "Raster opacity", 
                min = 0, max = 1, step = 0.05, 
                value = 0.75)),
  
  # stats re dimension distributions against new map
  absolutePanel(
    top = 10, right = 10, draggable = FALSE, width = "400px", style = "z-index:500; min-width: 200px;",
    id = "stats", style='background-color: rgba(255, 255, 255,  0.95); border-color: rgba(0,0,0,1); border-width: 10px; padding: 5px;',
    h6("(stats for current map)", align = "center"),
    tableOutput('statsTable'),
    cursor = "inherit"),
  
  absolutePanel(
    top = 10, left = 10, draggable = FALSE,
    id = "downld", 
    downloadButton(outputId = "DOWNLOADtile", label = "Download indicator raster"),
    downloadButton(outputId = "DOWNLOADshape", label = "Download hotspot shapefile"),
    
    conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                     tags$div("Loading...",id="loadmessage"),
                     style = "z-index:600;",
                     tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               bottom: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #FFD300;
               z-index: 105;
             }
          "))))
)
  
  # absolutePanel(
  #   top = "55%", right = 10, draggable = FALSE,
  #   id = "downld", 
  #   downloadButton(outputId = "DOWNLOADshape", label = "Download hotspot shapefile")
  # )
  

# 
# lat <- c(rep(60, 5), rep(30, 5), rep(0, 5), rep(-30, 5), rep(-60, 5))
# lon <- c(rep(c(-150, -60, 0, 60, 150), 5))
# label = rep("Lodaing...", 25)
# label.df <- data.frame(lat, lon, label)
# label.shp <- st_as_sf(label.df, coords = c("lon", "lat"), crs = 4326)

# options = leafletOptions(zoomControl = FALSE)
# Create server
server <- function(input, output, session) {
  
  # Create default map
  output$mymap <- renderLeaflet({
    leaflet(options = leafletOptions(zoomControl = FALSE, worldCopyJump = FALSE,
                                     minZoom = 2)) %>%
      addProviderTiles(providers$CartoDB) %>%
      addRasterImage(OverallComposite, colors = plt, opacity = 0.75, project = TRUE) %>% 
      addPolylines(data=World84,
                   layerId = World84$id,
                   color = "black", 
                   weight = 0.5,
                   opacity = 0.8) %>% 
      addLegend(position = "bottomright", pal = palr, values = values(OverallComposite),
                opacity = 1, title = "Hotspot <br> index",
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))) %>%
      addPolylines(data=Hotspots_all.sf,
                   layerId = Hotspots_all.sf$layer,
                   color = "black", 
                   weight = 3,
                   opacity = 1) %>% 
      addLegendCustom1("white", "Original hotspots", 15, "square", "black",
                       location = "bottomright") %>%
      setView(-30, 30, zoom = 3) %>%
      setMaxBounds(lng1 = -210, lat1 = -80,
                   lng2 = 210, lat2 = 90) 

      } )
  
  # Check to see if weightings sum to 1
  output$sumCheck <- renderText({
    sum.v <- input$Pop.w + input$Ag.w + input$Econ.w + input$Env.w
    validate(need(sum.v < 1.01 & sum.v > 0.99, "The weights need to sum to 1"))
    })
  
  observe({
    if((input$Pop.w + input$Ag.w + input$Econ.w + input$Env.w) < 0.99 |
       (input$Pop.w + input$Ag.w + input$Econ.w + input$Env.w) > 1.01 ) {
      disable("Recalc")
    }
    else{
      enable("Recalc")
    }
  })
  
  output$statsTable <- renderTable(df)
   
  v <- reactiveValues(counter = 0)
  v$counter = 0
  
  # Recalculate raster, hotspots, and stats on click
  observeEvent(input$Recalc, {
    Sys.sleep(1.0)  
    v$counter <- v$counter + 1L
    maxi = 12
    for(i in 0:maxi) {
      if (i == 0) {
        disable("downld")
        OverallComposite.new <- NA
        
        
        
      }
      if (i == 1) {
        OverallComposite.new <- ((input$Pop.w * Pop) + (input$Ag.w * Cal) + 
                                   (input$Econ.w * GDP) + (input$Env.w * Eco))
      }
      if (i == 2) {
        Hotspots_all.new <- hotspot_id_smoother(RawRaster = OverallComposite.new,
                                                WeightRaster = ga,
                                                MaskRaster = Mask,
                                                Method = input$Method,
                                                Multiplier = input$Mult)
      }
      if (i == 3) {
        Hotspots_all.new[Hotspots_all.new == 0] <- NA
      }
      if (i == 4) {
        Hotspots_all.sf.new <- rasterToPolygons(Hotspots_all.new, dissolve = TRUE)
      }
      if (i == 5) {
        Hotspots_all.sf.new <- sf::st_as_sf(Hotspots_all.sf.new)
        Hotspots_all.sf.new <- st_buffer(Hotspots_all.sf.new, dist = 0)
      }
      if (i == 6) {
        leafletProxy("mymap") %>% 
          clearTiles %>% clearImages() %>% clearShapes() %>%  clearMarkers() %>% clearControls() %>%
          addProviderTiles(providers$CartoDB) %>%
          addRasterImage(OverallComposite.new, colors = plt, opacity = input$opac, project = TRUE) %>% 
          addPolylines(data=World84,
                       layerId = World84$id,
                       color = "black", 
                       weight = 0.5,
                       opacity = 0.8) %>% 
          addLegend(position = "bottomright", pal = palr, values = values(OverallComposite.new),
                    opacity = 1, title = "Hotspot <br> index",
                    labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE))) %>%
          addPolylines(data=Hotspots_all.sf,
                       layerId = Hotspots_all.sf$layer,
                       color = "#00FFFF", 
                       weight = 1,
                       opacity = 1) %>% 
          addPolylines(data=Hotspots_all.sf.new,
                       layerId = Hotspots_all.sf.new$layer,
                       color = "black", 
                       weight = 3,
                       opacity = 1) %>% 
          addLegendCustom1("white", "Hotspots based on <br> current weightings", 15, "square", "black",
                           location = "bottomright") %>% 
        addLegendCustom1("white", "Original hotspots", 15, "square", "#00FFFF",
                         location = "bottomright")
    }
      if (i == 7) {
        temp_dir <- tmpDir()
        # write raster to temp directory
        writeRaster(x = OverallComposite.new, 
                    filename = paste0(temp_dir, "/CustomHotspotIndicator.tif", sep = ""), 
                    format = "GTiff", overwrite = TRUE)
        
        # write shapefile to temp directory
        poly_hotspot <- sf::st_as_sf(Hotspots_all.sf.new)
        
        sf::write_sf(poly_hotspot, paste0(temp_dir, "/CustomHotspotShapefile.shp", sep = ""),
                     driver = "ESRI Shapefile", overwrite = T)
                
        # writeOGR(obj = poly_hotspot,
        #          dsn = paste0(temp_dir, "/CustomHotspotShapefile.shp", sep = ""),
        #          layer = "CustomHotspotShapefile",
        #          driver = "ESRI Shapefile", 
        #          overwrite_layer = TRUE)
      } 
      
      if (i == 8) {
        
        PopCount <- raster("./gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif")
        PopCount[is.na(Hotspots_all.new)] <- NA
        PopAll <- round(sum(PopCount[], na.rm = T)/1e6, 4)
        
        PopCount <- raster("./gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif")
        PopCount[is.na(Hotspots_all.new) | Trends >= 0] <- NA
        PopDry <- round(sum(PopCount[], na.rm = T)/1e6, 4)
        
        PopCount <- raster("./gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif")
        PopCount[is.na(Hotspots_all.new) | Trends < 0] <- NA
        PopWet <- round(sum(PopCount[], na.rm = T)/1e6, 4)
        }
      if (i == 9) {
        CalCount <- raster("./kcal_0d5.tif")
        CalCount[is.na(Hotspots_all.new)] <- NA
        CalAll <- round(sum(CalCount[], na.rm = T)/1e12, 4)
        
        CalCount <- raster("./kcal_0d5.tif")
        CalCount[is.na(Hotspots_all.new) | Trends >= 0] <- NA
        CalDry <- round(sum(CalCount[], na.rm = T)/1e12, 4)
        
        CalCount <- raster("./kcal_0d5.tif")
        CalCount[is.na(Hotspots_all.new) | Trends < 0] <- NA
        CalWet <- round(sum(CalCount[], na.rm = T)/1e12, 4)
        
      }
      if (i == 10) {
        GDPCount <- raster("./GDP_2015_0d5.tif")
        GDPCount[is.na(Hotspots_all.new)] <- NA
        GDPAll <- round(sum(GDPCount[], na.rm = T)/1e12, 4)
        
        GDPCount <- raster("./GDP_2015_0d5.tif")
        GDPCount[is.na(Hotspots_all.new) | Trends >= 0] <- NA
        GDPDry <- round(sum(GDPCount[], na.rm = T)/1e12, 4)
        
        GDPCount <- raster("./GDP_2015_0d5.tif")
        GDPCount[is.na(Hotspots_all.new) | Trends < 0] <- NA
        GDPWet <- round(sum(GDPCount[], na.rm = T)/1e12, 4)
      }
      if (i == 11) {
        
        EcoSens <- raster("./EcologicalSensitivity_Indicator.tif")
        ga <- WGS84_areaRaster(0.5)
        ga[is.na(Hotspots_all.new) | EcoSens < 0.9] <- NA
        EcoAll <- round(sum(ga[], na.rm = T)/1e6, 4)
        
        EcoSens <- raster("./EcologicalSensitivity_Indicator.tif")
        ga <- WGS84_areaRaster(0.5)
        ga[is.na(Hotspots_all.new) |  EcoSens < 0.9 | Trends >= 0] <- NA
        EcoDry <- round(sum(ga[], na.rm = T)/1e6, 4)
        
        EcoSens <- raster("./EcologicalSensitivity_Indicator.tif")
        ga <- WGS84_areaRaster(0.5)
        ga[is.na(Hotspots_all.new) |  EcoSens < 0.9 | Trends < 0] <- NA
        EcoWet <- round(sum(ga[], na.rm = T)/1e6, 4)
        }
      if (i == 12) {
        df <- data.frame(Dimension = c("Population", "Crop calories", "GDP (PPP)", 
                                  "Water sensitive environ."), 
                         `In Drying Hotspots` = c(PopDry, CalDry, GDPDry, EcoDry),
                         `In Wetting Hotspots` = c(PopWet, CalWet, GDPWet, EcoWet),
                         Units = c("million people", "trillion crop calories",
                                   "trillion USD (2011)", "million sq. km"))
        names(df) <- c("Dimension", "In Drying Hotspots", "In Wetting Hotspots",
                       "Units")
        output$statsTable <- renderTable(df)
        enable("downld")
      }
      updateProgressBar(session = session, id = "pb4", value = (i/maxi)*100)
      Sys.sleep(0.1)  
      } 
    },
    ignoreNULL = TRUE
    )
  
  output$DOWNLOADtile <- downloadHandler(
    
    filename = function() {
      if (v$counter == 0) {
        paste("OriginalIndicator_", Sys.Date(), ".tif", sep="") 
      } else if (v$counter!= 0 ) {
        paste("CustomIndicator_", Sys.Date(), ".tif", sep="") 
      } },
      
    content = function(file) {
      if (v$counter== 0) {
        raster.file <- "./DefaultComposite.tif"
        file.copy(raster.file, file)
      } else if (v$counter!= 0) {
        
        raster.file <- list.files(temp_dir,
                                  "CustomHotspotIndicator.*",
                                  full.names = TRUE)
        
        file.copy(raster.file, file)
      } 
      }
    )
  
  output$DOWNLOADshape <- downloadHandler(
    filename = function() {
      paste0("HotspotShapefile.zip")
      },
    
    content = function(file) {
      
      if (v$counter== 0) {
        zip(zipfile='HotspotShapefile.zip', files = "./DefaultHotspots.*")  
        file.copy("./HotspotShapefile.zip", file)
      } else if (v$counter!= 0) {
        zip(zipfile='CustomHotspotShapefile.zip', files = list.files(temp_dir,
                                                               "CustomHotspotShapefile.*",
                                                               full.names = TRUE))
        file.copy("./CustomHotspotShapefile.zip", file, overwrite = TRUE)
      }
    }
  )
    
}
    
# Run app
shinyApp(ui, server)
