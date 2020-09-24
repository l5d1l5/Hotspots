## this script maps the combined hotspots with all combinations of weightings
## and selects the hotpots based on the elbow method 

library(tidyverse); library(magrittr); library(ncdf4); library(here)
library(raster); library(sf); library(tmap); library(tmaptools); 
library(fasterize); library(RColorBrewer); library(pals);
library(ggplot2); library(scales); library(spatstat); library(ahp)
source(here::here("Scripts", "gen_funs.r"))

# Import each of the dimensions
Pop <- raster(here::here("ProducedData", "IndividualDims", "Pop_indicator.tif"))
Cal <- raster(here::here("ProducedData", "IndividualDims", "Cal_indicator.tif"))
GDP <- raster(here::here("ProducedData", "IndividualDims", "GDP_indicator.tif"))
Eco <- raster(here::here("ProducedData", "IndividualDims", "Eco_sens_indicator.tif"))
ga <- raster(here::here("ProducedData", "WGS_gridArea_halfdegree.tif"))
Mask <- raster(here::here("ProducedData", "GlobalMask.tif"))

# Establish AHP preferences, where 1 = Ecology, 2 = population, 3 = ag, 4 = econ.
BDC.pw <- rbind(c(1, 2, 2),
                c(1, 3, 2),
                c(1, 4, 4),
                c(2, 3, 1),
                c(2, 4, 2),
                c(3, 4, 2)) %>% as.data.frame()

BdayCake <- ahp::AhpMatrix(BDC.pw)
BdayCake.weights <- ahp::PrioritiesFromPairwiseMatrixEigenvalues(mat = BdayCake)
BdayCake.weights

Eco.w <- BdayCake.weights$priority[1]
Pop.w <- BdayCake.weights$priority[2]
Ag.w <- BdayCake.weights$priority[3]
Econ.w <- BdayCake.weights$priority[4]

OverallComposite <- Eco.w*Eco + Pop.w*Pop + Ag.w*Cal + Econ.w*GDP 

Hotspots_all <- hotspot_id_smoother(RawRaster = OverallComposite,
                                    WeightRaster = ga,
                                    MaskRaster = Mask,
                                    Method = "SD",
                                    Multiplier = 2)
Hotspots_pop <- hotspot_id_smoother(RawRaster = Pop,
                                    WeightRaster = ga,
                                    MaskRaster = Mask,
                                    Method = "SD",
                                    Multiplier = 2)
Hotspots_ag <- hotspot_id_smoother(RawRaster = Cal,
                                   WeightRaster = ga,
                                   MaskRaster = Mask,
                                   Method = "SD",
                                   Multiplier = 2)
Hotspots_Econ <- hotspot_id_smoother(RawRaster = GDP,
                                     WeightRaster = ga,
                                     MaskRaster = Mask,
                                     Method = "SD",
                                     Multiplier = 2)
Hotspots_Ecol <- hotspot_id_smoother(RawRaster = Eco,
                                     WeightRaster = ga,
                                     MaskRaster = Mask,
                                     Method = "SD",
                                     Multiplier = 2)

Hotspots_nonSubstitutable <- Hotspots_pop + Hotspots_ag + Hotspots_Econ + Hotspots_Ecol
Hotspots_nonSubstitutable[Hotspots_nonSubstitutable >= 1] <- 1
plot(Hotspots_nonSubstitutable)

################# Now create maps... 
# Convert all hotspots to shapefiles for plotting
Hotspots_all[Hotspots_all == 0] <- NA
Hotspots_pop[Hotspots_pop == 0] <- NA
Hotspots_ag[Hotspots_ag == 0] <- NA
Hotspots_Econ[Hotspots_Econ == 0] <- NA
Hotspots_Ecol[Hotspots_Ecol == 0] <- NA
Hotspots_nonSubstitutable[Hotspots_nonSubstitutable == 0] <- NA

Hotspots_all.sf <- rasterToPolygons(Hotspots_all, dissolve = T)
Hotspots_all_ns.sf <- rasterToPolygons(Hotspots_nonSubstitutable, dissolve = T)
Hotspots_pop.sf <- rasterToPolygons(Hotspots_pop, dissolve = T)
Hotspots_ag.sf <- rasterToPolygons(Hotspots_ag, dissolve = T)
Hotspots_Econ.sf <- rasterToPolygons(Hotspots_Econ, dissolve = T)
Hotspots_Ecol.sf <- rasterToPolygons(Hotspots_Ecol, dissolve = T)

## Create tmap elements
data(World)
lakes <- sf::read_sf(here::here("Data", "NE_lakes", "ne_10m_lakes.shp"))
lakes$scalerank %<>% as.numeric()
lakes.l <- subset(lakes, scalerank < 2)
plt = pals::ocean.balance(20)[3:18] %>% rev()

Pop[Mask != 1] <- NA
Cal[Mask != 1] <- NA
GDP[Mask != 1] <- NA
Eco[Mask != 1] <- NA

# crop and reproject all for plotting
Pop.rpj <- tmap_clipproj(Pop)
Cal.rpj <- tmap_clipproj(Cal)
GDP.rpj <- tmap_clipproj(GDP)
Eco.rpj <- tmap_clipproj(Eco)
Overall.rpj <- tmap_clipproj(OverallComposite)

##### - POPULATION DIMENSION
plot_pop <- tm_shape(Pop.rpj, projection = "+proj=robin") +
  tm_raster(style = "cont", palette = plt, midpoint = 0, breaks = c(-1, 1)) +
  tm_shape(World) + tm_borders(lwd = 0.7) +
  tm_shape(lakes.l) + tm_polygons(col = "white", border.col = "grey50", lwd = 0.8) +
  tm_shape(Hotspots_pop.sf) + tm_borders(lwd = 1.8, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
plot_pop

##### - AGRICULTURAL DIMENSION
plot_cal <- tm_shape(Cal.rpj, projection = "+proj=robin") +
  tm_raster(style = "cont", palette = plt, midpoint = 0, breaks = c(-1, 1)) +
  tm_shape(World) + tm_borders(lwd = 0.7) +
  tm_shape(lakes.l) + tm_polygons(col = "white", border.col = "grey50", lwd = 0.8) +
  tm_shape(Hotspots_ag.sf) + tm_borders(lwd = 1.8, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
plot_cal

##### - ECONOMIC DIMENSION
plot_gdp <- tm_shape(GDP.rpj, projection = "+proj=robin") +
  tm_raster(style = "cont", palette = plt, midpoint = 0, breaks = c(-1, 1)) +
  tm_shape(World) + tm_borders(lwd = 0.7) +
  tm_shape(lakes.l) + tm_polygons(col = "white", border.col = "grey50", lwd = 0.8) +
  tm_shape(Hotspots_Econ.sf) + tm_borders(lwd = 1.8, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
plot_gdp

##### - ECOLOGICAL DIMENSION
plot_ecol <- tm_shape(Eco.rpj, projection = "+proj=robin") +
  tm_raster(style = "cont", palette = plt, midpoint = 0, breaks = c(-1, 1)) +
  tm_shape(World) + tm_borders(lwd = 0.7) +
  tm_shape(lakes.l) + tm_polygons(col = "white", border.col = "grey50", lwd = 0.8) +
  tm_shape(Hotspots_Ecol.sf) + tm_borders(lwd = 1.8, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
plot_ecol

######################################################.
######### - OVERALL HOTSPOTS - #######################.
######################################################.
plot_all <- tm_shape(Overall.rpj, projection = "+proj=robin") +
  tm_raster(style = "cont", palette = plt, midpoint = 0, breaks = c(-1, 1)) +
  tm_shape(World) + tm_borders(lwd = 0.7) +
  tm_shape(lakes.l) + tm_polygons(col = "white", border.col = "grey50", lwd = 0.8) +
  tm_shape(Hotspots_all_ns.sf) + tm_borders(lwd = 1, col = "black", lty = "dotted") +
  tm_shape(Hotspots_all.sf) + tm_borders(lwd = 1.8, col = "black") +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
plot_all

tmap::tmap_save(plot_pop, here::here("Figures", "Hotspot_maps", "Pop_hotspots.svg"), 
                units = "in", dpi = 500) 
tmap_save(plot_cal, here::here("Figures", "Hotspot_maps", "Cal_hotspots.svg"), 
          units = "in", dpi = 500) 
tmap_save(plot_gdp, here::here("Figures", "Hotspot_maps", "GDP_hotspots.svg"), 
          units = "in", dpi = 500) 
tmap_save(plot_ecol, here::here("Figures", "Hotspot_maps", "Ecol_hotspots.svg"), 
          units = "in", dpi = 500) 
tmap_save(plot_all, here::here("Figures", "Hotspot_maps", "Composite_hotspots_sdgs.svg"), 
          units = "in", dpi = 500) 

##################################################
####### Calculate statistics in hotspots #########
##################################################
Hotspots_all
Pop.r <- raster(here::here("Data", "Dimensions", 
                           "gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif"))
Pop.r[is.na(Hotspots_all)] <- NA
sum(Pop.r[], na.rm = T)/1e9

Cal.r <- raster(here::here("Data", "Dimensions", "kcal_0d5.tif"))
Cal.r[is.na(Hotspots_all)] <- NA
sum(Cal.r[], na.rm = T)/1e12

GDP.r <- raster(here::here("Data", "Dimensions", "GDP_2015_0d5.tif"))
GDP.r[is.na(Hotspots_all)] <- NA
sum(GDP.r[], na.rm = T)/1e12

Eco.r <- raster(here::here("ProducedData", "EcologicalSensitivity_Indicator.tif"))
ga <- raster(here::here("ProducedData", "WGS_gridArea_halfdegree.tif"))
ga[is.na(Hotspots_all) | Eco.r < 0.9] <- NA
ga[Eco.r < 0.9 | Mask != 1] <- NA
sum(ga[], na.rm = T)/1e6

ga <- raster(here::here("ProducedData", "WGS_gridArea_halfdegree.tif"))
ga[is.na(Hotspots_all) | Mask != 1] <- NA
ga[Mask != 1] <- NA
sum(ga[], na.rm = T)/1e6

twst_i <- raster(here::here("ProducedData", "TWS_MAP_sclaer_p10p90_1985_2014map.tif"))
ga <- raster(here::here("ProducedData", "WGS_gridArea_halfdegree.tif"))
ga[is.na(Hotspots_all) | twst_i >= 0] <- NA 
a = sum(ga[], na.rm = T)

ga <- raster(here::here("ProducedData", "WGS_gridArea_halfdegree.tif"))
ga[is.na(Hotspots_all)] <- NA
b = sum(ga[], na.rm = T)
a/b


###################################################
####### Calculate countries per dimension #########
###################################################

Admin0 <- sf::read_sf(here::here("Data", "Admin", "g2008_1.shp"))
Admin0.r <- fasterize(sf = Admin0, raster = ga, field = "ADM0_CODE")

count.df <- cbind(Admin0.r[], Hotspots_pop[], Hotspots_ag[],
                 Hotspots_Econ[], Hotspots_Ecol[]) %>%
  as.data.frame() %>%
  set_colnames(c("ADM", "POP", "AG", "ECON", "ECOL"))

count.df %>% filter(POP == 1) %>% pull(ADM) %>% unique() %>% length()
count.df %>% filter(AG == 1) %>% pull(ADM) %>% unique() %>% length()
count.df %>% filter(ECON == 1) %>% pull(ADM) %>% unique() %>% length()
count.df %>% filter(ECOL == 1) %>% pull(ADM) %>% unique() %>% length()

#############################################################
####### Calculate drying area percent per dimension #########
#############################################################

twst_i <- raster(here::here("ProducedData", "TWS_MAP_sclaer_p10p90_1985_2014map.tif"))
trend.df <- cbind(twst_i[], Hotspots_pop[], Hotspots_ag[],
                  Hotspots_Econ[], Hotspots_Ecol[], ga[]) %>%
  as.data.frame() %>%
  set_colnames(c("trend", "POP", "AG", "ECON", "ECOL", "AREA"))

trend.df %>% filter(ECOL == 1 & trend < 0) %>% pull(AREA) %>% sum(na.rm = T) /
  trend.df %>% filter(ECOL == 1) %>% pull(AREA) %>% sum(na.rm = T)