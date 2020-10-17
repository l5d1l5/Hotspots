library(tidyverse); library(magrittr); library(here)
library(raster); library(sf); library(spatstat); library(ncdf4); library(scales)
library(fasterize); library(ggplot2); library(cartography)
library(Weighted.Desc.Stat); library(Hmisc)
source(here::here("Scripts", "gen_funs.r"))

# Import data
ga <- WGS84_areaRaster(0.5)
twst <- raster(here::here("Data", "Rodell_GRACE_TWSt_raw.tif"))*10 # cm to mm
MAP <- raster(here::here("ProducedData", "MeanAnnualPrecip_1985_2014.tif"))
Mask <- raster(here::here("ProducedData", "GlobalMask.tif"))
Pop.r <- raster(here::here("Data", "Dimensions", 
                           "gpw_v4_population_count_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_min.tif"))
Cal.r <- raster(here::here("Data", "Dimensions", "kcal_0d5.tif"))
GDP.r <- raster(here::here("Data", "Dimensions", "GDP_2015_0d5.tif"))
EcoSens <- raster(here::here("ProducedData", "EcologicalSensitivity_Indicator.tif"))

#1 Reclassify TWSt.rt raster into bins
twst[Mask == 0] <- NA
twst[1,] <- NA;twst[,1] <- NA;
MAP[Mask == 0] <- NA
MAP[1,] <- NA;MAP[,1] <- NA;
TWSt_MAP <- twst/(MAP)
TWSt_MAP[MAP == 0] <- NA

# Determine standard deviation of trends
df.stat <- cbind(TWSt_MAP[], ga[]) %>% as.data.frame() %>% set_colnames(c("twst", "ga"))
df.stat <- df.stat[complete.cases(df.stat),]

glbmd <- Hmisc::wtd.quantile(x = df.stat$twst, weights = df.stat$ga, probs = 0.5)
# calculate 1.5*IQR instead of stdev as skewed
x = 0.25
lim = Hmisc::wtd.quantile(x = df.stat$twst, weights = df.stat$ga, probs = 0.5+x) - 
  Hmisc::wtd.quantile(x = df.stat$twst, weights = df.stat$ga, probs = 0.5-x) 
lim = 1.5*lim

l = -0.08
h = 0.08
b = (h-l)/(70-2)
low = c(-Inf, seq(l, h, by = b))
high = c(seq(l, h, by = b), Inf)
ReCLASS = c(seq(l-0.5*b, h+0.5*b, by = b))
rcl.mx <- data.frame(low, high, ReCLASS)
TWSt.bin <- raster::reclassify(TWSt_MAP, rcl.mx)
# plot(TWSt.bin)

master.df <- cbind(TWSt_MAP[], TWSt.bin[], Pop.r[], Cal.r[], GDP.r[], 
                   EcoSens[], ga[]) %>% as.data.frame() %>%
  set_colnames(c("trend", "Bin", "Pop", "Cal", "GDP", "EcoSens", "ga"))

plot.df <- master.df %>% group_by(Bin) %>%
  summarise(Pop.c = sum(Pop, na.rm = T),
            Cal.c = sum(Cal, na.rm = T),
            GDP.c = sum(GDP, na.rm = T)) %>% as.data.frame()
plot.df$Ext <- ifelse(plot.df$Bin <= -(lim-0.25*b), -1, ifelse(plot.df$Bin >= (lim-0.25*b), 1, 0))
# write.csv(x = plot.df, file = here::here("Figures", "Figure2_plots", "AbsoluteStats_3.csv"))

theme1 = theme_minimal()+
  theme(legend.title = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

library(pals)
plt = pals::ocean.balance(20) %>% rev()
cols = c("-1" = plt[5], "0" = "grey70", "1" = plt[15])

ggplot(data = plot.df, aes(x = Bin, y = Pop.c/1e6, fill = as.factor(Ext))) +
  geom_bar(stat = "identity", lwd = 0, width =b*1.05) +
  scale_fill_manual(values = cols)+
  geom_vline(xintercept = -lim, lwd = 1, lty = "dashed", col = "grey40") +
  geom_vline(xintercept = lim, lwd = 1, lty = "dashed", col = "grey40") +
  geom_vline(xintercept = 0, lwd = 1, lty = "dotted", col = "grey40") +
  geom_vline(xintercept = glbmd, lwd = 1, col = "grey40") +
  coord_cartesian(xlim = c(l - b, h + b), ylim = c(0,1100), expand = c(0,0))+
  theme1 + theme(axis.ticks.x = element_line(size = 1))
ggsave(plot = last_plot(), here::here("Figures", "Histograms", "Pop_relative.eps"),
       dpi = 500, width = 6, height = 4, units = "in")

ggplot(data = plot.df, aes(x = Bin, y = Cal.c/1e14, fill = as.factor(Ext))) +
  geom_bar(stat = "identity", lwd = 0, width =b*1.05) +
  scale_fill_manual(values = cols) +
  geom_vline(xintercept = -lim, lwd = 1, lty = "dashed", col = "grey40") +
  geom_vline(xintercept = lim, lwd = 1, lty = "dashed", col = "grey40") +
  geom_vline(xintercept = 0, lwd = 1, lty = "dotted", col = "grey40") +
  geom_vline(xintercept = glbmd, lwd = 1, col = "grey40") +
  coord_cartesian(xlim = c(l - b, h + b), ylim = c(0,13), expand = c(0,0))+
  theme1 + theme(axis.ticks.x = element_line(size = 1))
ggsave(plot = last_plot(), here::here("Figures", "Histograms", "Cal_relative.eps"),
       dpi = 500, width = 6, height = 4, units = "in")

ggplot(data = plot.df, aes(x = Bin, y = GDP.c/1e12, fill = as.factor(Ext))) +
  geom_bar(stat = "identity", lwd = 0, width =b*1.05) +
  scale_fill_manual(values = cols) +
  geom_vline(xintercept = -lim, lwd = 1, lty = "dashed", col = "grey40") +
  geom_vline(xintercept = lim, lwd = 1, lty = "dashed", col = "grey40") +
  geom_vline(xintercept = 0, lwd = 1, lty = "dotted", col = "grey40") +
  geom_vline(xintercept = glbmd, lwd = 1, col = "grey40") +
  coord_cartesian(xlim = c(l - b, h + b), ylim = c(0,18), expand = c(0,0))+
  theme1 + theme(axis.ticks.x = element_line(size = 1))
ggsave(plot = last_plot(), here::here("Figures", "Histograms", "GDP_relative.eps"),
       dpi = 500, width = 6, height = 4, units = "in")

## Now calculate area distribution of top 50% of water sensitivit ecological areas
plot.df <- master.df %>% filter(EcoSens > 0.5) %>% group_by(Bin) %>%
  summarise(Area = sum(ga, na.rm = T)) %>% as.data.frame()

plot.df$Ext <- ifelse(plot.df$Bin <= -(lim-0.25*b), -1, ifelse(plot.df$Bin >= (lim-0.25*b), 1, 0))
# write.csv(x = plot.df, file = here::here("Figures", "Figure2_plots", "AbsoluteStats_eco.csv"))

ggplot(data = plot.df, aes(x = Bin, y = Area/1e6, fill = as.factor(Ext))) +
  geom_bar(stat = "identity", lwd = 0, width =b*1.05) +
  scale_fill_manual(values = cols) +
  geom_vline(xintercept = -lim, lwd = 1, lty = "dashed", col = "grey40") +
  geom_vline(xintercept = lim, lwd = 1, lty = "dashed", col = "grey40") +
  geom_vline(xintercept = 0, lwd = 1, lty = "dotted", col = "grey40") +
  geom_vline(xintercept = glbmd, lwd = 1, col = "grey40") +
  coord_cartesian(xlim = c(l - b, h + b), ylim = c(0,10), expand = c(0,0))+
  theme1 + theme(axis.ticks.x = element_line(size = 1))
ggsave(plot = last_plot(), here::here("Figures", "Histograms", "Eco_relative.eps"),
       dpi = 500, width = 6, height = 4, units = "in")


### Calculate dimension counts in each range
All.df <- cbind(Pop.r[], Cal.r[], GDP.r[], EcoSens[], ga[], TWSt_MAP[]) %>%
  set_colnames(c("pop", "cal", "gdp", "eco", "area", "trend")) %>% 
  as.data.frame()

All.df %>% filter(trend < -lim) %>% pull(pop) %>% sum(na.rm = T)/1e6
All.df %>% filter(trend > lim) %>% pull(pop) %>% sum(na.rm = T)/1e6

All.df %>% filter(trend < -lim) %>% pull(cal) %>% sum(na.rm = T)/1e12
All.df %>% filter(trend > lim) %>% pull(cal) %>% sum(na.rm = T)/1e12

All.df %>% filter(trend < -lim) %>% pull(gdp) %>% sum(na.rm = T)/1e12
All.df %>% filter(trend > lim) %>% pull(gdp) %>% sum(na.rm = T)/1e12

All.df %>% filter(trend < -lim & eco >= 0.5) %>% pull(area) %>% sum(na.rm = T)/1e6
All.df %>% filter(trend > lim & eco >= 0.5) %>% pull(area) %>% sum(na.rm = T)/1e6
