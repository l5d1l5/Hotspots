## This plot takes population, calories, and GDP data
## reclassifies based on density percentiles
## multiplies by the TWSt/MAP developed scalar
## and creates an exposure indicator for each dimension

library(tidyverse); library(magrittr); library(here)
library(raster); library(sf); library(spatstat); library(ncdf4)
library(fasterize); library(ggplot2); 

# load general files
source(here::here("Scripts", "gen_funs.r"))
ga <- WGS84_areaRaster(0.5)
twst_i <- raster(here::here("ProducedData", "TWS_MAP_scaler_p10p90_1985_2014map.tif")) 
Mask <- raster(here::here("ProducedData", "GlobalMask.tif"))

####################1. population dimension
Pop.r <- raster(here::here("Data", "Dimensions", 
                           "gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif"))
Pop.r[Pop.r < 1] <- NA
Pop.dens <- Pop.r/ga
pop.pctl <- RasterAreaPercentiles(RasterToClassify = Pop.dens,
                                  WeightRaster = ga, 
                                  MaskRaster =  Mask, 
                                  clipToExtent =  "clip")
# plot global distribution
tmp <- cbind(pop.pctl[], ga[], Mask[], Pop.r[]) %>% as.data.frame() %>% set_colnames(c("a", "b", "c", "d"))
sm <- tmp[complete.cases(tmp$d),]
sm <- sm %>% filter(c == 1) %>% group_by(a) %>% summarize(area = sum(b, na.rm = T)) %>% as.data.frame()
sm$sum <- (cumsum(sm$area)/sum(sm$area,na.rm = T)) %>% round(digits = 2)
plot(sm$sum  ~ sm$a); sm
pop.pctl[is.na(pop.pctl[]) & Mask[] == 1] <- 0
plot(pop.pctl)

# create population hotspot indicator
pop_ind <- pop.pctl*twst_i
plot(pop_ind)

####################2. agricultural dimension
Cal.r <- raster(here::here("Data", "Dimensions", "kcal_0d5.tif"))
Cal.r[Cal.r[] < 1] <- NA
Cal.dens <- Cal.r/ga
Cal.pctl <- RasterAreaPercentiles(RasterToClassify = Cal.dens,
                                  WeightRaster = ga, 
                                  MaskRaster =  Mask, 
                                  clipToExtent =  "clip")

# plot global distribution
tmp <- cbind(Cal.pctl[], ga[], Mask[], Cal.r[]) %>% as.data.frame() %>% set_colnames(c("a", "b", "c", "d"))
sm <- tmp[complete.cases(tmp$d),]
sm <- sm %>% filter(c == 1) %>% group_by(a) %>% summarize(area = sum(b, na.rm = T)) %>% as.data.frame()
sm$sum <- (cumsum(sm$area)/sum(sm$area,na.rm = T)) %>% round(digits = 2)
plot(sm$sum  ~ sm$a); sm
Cal.pctl[is.na(Cal.pctl[]) & Mask[] == 1] <- 0
plot(Cal.pctl)

# create agricultural indicator
Cal_ind <- Cal.pctl*twst_i
plot(Cal_ind)

#3. economic dimension
EconSect_withd <- raster(here::here("Data", "Dimensions", "EconSectors_2010Withdrawals_sum.tif"))
EconSect_withd[EconSect_withd == 0] <- NA
Econ.pctl <- RasterAreaPercentiles(RasterToClassify = EconSect_withd,
                                   WeightRaster = ga, 
                                   MaskRaster =  Mask, 
                                   clipToExtent =  "clip")

# plot global distribution
tmp <- cbind(Econ.pctl[], ga[], Mask[], EconSect_withd[]) %>% as.data.frame() %>% 
  set_colnames(c("a", "b", "c", "d"))
sm <- tmp[complete.cases(tmp$d),]
sm <- sm %>% filter(c == 1) %>% group_by(a) %>% summarize(area = sum(b, na.rm = T)) %>% as.data.frame()
sm$sum <- (cumsum(sm$area)/sum(sm$area,na.rm = T)) %>% round(digits = 2)
plot(sm$sum  ~ sm$a); sm
Econ.pctl[is.na(Econ.pctl[]) & Mask[] == 1] <- 0
plot(Econ.pctl)

# create economic indicator
Econ_ind <- Econ.pctl*twst_i
plot(Econ_ind)

#4. Ecological dimension
## Sensitivity (vegetation to soil moisture and Eflos to GW head decline)
VSI_sm <- raster(here::here("Data", "Dimensions", "VSI_aetw_0d5.tif"))
crs(VSI_sm) <- crs("+proj=longlat")
VSI.pctl <- RasterAreaPercentiles(RasterToClassify = VSI_sm,
                                  WeightRaster = ga, 
                                  MaskRaster =  Mask, 
                                  clipToExtent =  "clip")
VSI.pctl[1,] <- NA; VSI.pctl[,1] <- NA

# plot global distribution
tmp <- cbind(VSI.pctl[], ga[], Mask[], VSI_sm[]) %>% as.data.frame() %>% set_colnames(c("a", "b", "c", "d"))
sm <- tmp[complete.cases(tmp$d),]
sm <- sm %>% filter(c == 1) %>% group_by(a) %>% summarize(area = sum(b, na.rm = T)) %>% as.data.frame()
sm$sum <- (cumsum(sm$area)/sum(sm$area,na.rm = T)) %>% round(digits = 2)
plot(sm$sum  ~ sm$a); sm

VSI.pctl[is.na(VSI.pctl[]) & Mask[] == 1] <- 0
plot(VSI.pctl)

## Sensitivity (GW head decline and eflows)
EFs <- raster(here::here("Data", "Dimensions", "EFN_0d5.tif"))
crs(EFs) <- crs("+proj=longlat")
extent(EFs) <- extent(-180, 180, -90, 90)
EFN.pctl <- RasterAreaPercentiles(RasterToClassify = EFs,
                                  WeightRaster = ga, 
                                  MaskRaster =  Mask, 
                                  clipToExtent =  "clip")
EFN.pctl <- 1 - EFN.pctl # invert as lower values = more sensitive
EFN.pctl[1,] <- NA; EFN.pctl[,1] <- NA
EFN.pctl[is.na(EFs)] <- 0
EFN.pctl[is.na(EFN.pctl[]) & Mask[] == 1] <- 0
EFN.pctl[Mask[] != 1] <- NA

# plot global distribution
tmp <- cbind(EFN.pctl[], ga[], Mask[], EFs[]) %>% as.data.frame() %>% set_colnames(c("a", "b", "c", "d"))
sm <- tmp[complete.cases(tmp$d),]
sm <- sm %>% filter(c == 1) %>% group_by(a) %>% summarize(area = sum(b, na.rm = T)) %>% as.data.frame()
sm$sum <- (cumsum(sm$area)/sum(sm$area,na.rm = T)) %>% round(digits = 4)
plot(sm$sum  ~ sm$a); sm
plot(EFN.pctl)

# Combine ecological indicators and calculate average per grid cell
Sens.stack <- stack(EFN.pctl, VSI.pctl)
Sens.ind <- calc(Sens.stack, fun = mean, na.rm = T)
plot(Sens.ind)

# Sens.ind[Sens.ind < 0.0025] <- NA
Sens.ind[1,] <- NA; Sens.ind[,1] <- NA
Sens.ind[is.nan(Sens.ind) | Sens.ind == 0] <- NA
Sens.ind.pctl <- RasterAreaPercentiles(RasterToClassify = Sens.ind,
                                       WeightRaster = ga, 
                                       MaskRaster =  Mask, 
                                       clipToExtent =  "clip")

# plot global distribution
tmp <- cbind(Sens.ind.pctl[], ga[], Mask[], Sens.ind[]) %>% as.data.frame() %>% set_colnames(c("a", "b", "c", "d"))
sm <- tmp[complete.cases(tmp$d),]
sm <- sm %>% filter(c == 1) %>% group_by(a) %>% summarize(area = sum(b, na.rm = T)) %>% as.data.frame()
sm$sum <- (cumsum(sm$area)/sum(sm$area,na.rm = T)) %>% round(digits = 2)
plot(sm$sum  ~ sm$a); sm
Sens.ind.pctl[is.na(Sens.ind.pctl[]) & Mask[] == 1] <- 0
Sens.ind.pctl[1,] <- NA; Sens.ind.pctl[,1] <- NA
plot(Sens.ind.pctl)

writeRaster(Sens.ind.pctl, here::here("ProducedData", "EcologicalSensitivity_Indicator"),
            format = "GTiff", overwrite = T)

# now make overall sensitivity and priority indicator
Eco_ind <- Sens.ind.pctl*twst_i
plot(Eco_ind)

# write all rasters
write.stack = stack(pop_ind, Cal_ind, Econ_ind, Eco_ind)
names(write.stack) <- c("Pop", "Cal", "Econ", "Eco_sens")
for (i in 1:length(write.stack)) {
  writeRaster(write.stack[[i]], 
              filename = paste(here::here("ProducedData", "IndividualDims"),
                               "/", names(write.stack[[i]]), "_indicator.tif", sep = ""), 
              by.layer = TRUE, format = "GTiff", overwrite = T)
}

######################################################.
######### - make plots for SI - #######################.
######################################################.
twst_i <- tmap_clipproj(twst_i)
pop.pctl <- tmap_clipproj(pop.pctl)
Cal.pctl <- tmap_clipproj(Cal.pctl)
Econ.pctl <- tmap_clipproj(Econ.pctl)
Sens.ind.pctl <- tmap_clipproj(Sens.ind.pctl)

GDP.r <- tmap_clipproj(GDP.r)
Wat_withd <- tmap_clipproj(Wat_withd)
EconWithd <- tmap_clipproj(EconSect_withd)
VSI_sm <- tmap_clipproj(VSI_sm)
VSI.pctl <- tmap_clipproj(VSI.pctl)
EFs <- tmap_clipproj(EFs)
EFN.pctl <- tmap_clipproj(EFN.pctl)
Mean.p <- tmap_clipproj(Sens.ind)

plt = pals::ocean.balance(20)[3:18] %>% rev()
plt1 = pals::jet(20)

library(tmap); library(tmaptools)
data(World)
lakes <- sf::read_sf(here::here("Data", "NE_lakes", "ne_10m_lakes.shp"))
lakes$scalerank %<>% as.numeric()
lakes.l <- subset(lakes, scalerank < 2)

plot <- tm_shape(log10(EconWithd), projection = "+proj=robin") +
  # tm_raster(style = "cont", palette = plt, midpoint = 0, breaks = c(-1, 1)) +
  tm_raster(style = "cont", palette = plt1, breaks = c(0, 2)) +
  tm_shape(World) + tm_borders(lwd = 0.7) +
  tm_shape(lakes.l) + tm_polygons(col = "white", border.col = "grey50", lwd = 0.8) +
  tm_layout(legend.show = F, earth.boundary = c(-179, -60, 179, 88),
            earth.boundary.color = "white", space.color = "white",
            legend.frame = F, frame = F,
            outer.margins = c(-0, -0.09, -0, -0.04)) # B, L, T, R
plot

tmap::tmap_save(plot, here::here("Figures", "SI_figures", "EconWithd_log10.svg"), 
                units = "in", dpi = 500) 

## plot sm lines (throughaway code for SI figures)
ggplot(data = sm, aes(x = a, y = sum, color = a)) +
  geom_line(size = 2) +
  scale_color_gradientn(colours = plt1) +
  theme_bw() +
  theme(legend.position = "none") +
  coord_cartesian(expand = 0) +
  xlab("Percentile assignment") +
  ylab("Cumulative area (%)")

ggsave(plot = last_plot(), here::here("Figures", "SI_figures", "Ecol_curveline.eps"),
       dpi = 500, width = 6, height = 4, units = "in")
