################################################################################
# @Manuscript - "Hotspots of social and ecological impacts from freshwater stress and storage loss" (Huggins et al.) 
# @Description - Load all necessary libraries for the project, including custom user defined functions, and set custom ggplot theme customizations. 
################################################################################

# Load libraries
library(tidyverse)
library(magrittr) 
library(raster) 
library(sf) 
library(spatstat) 
library(ncdf4)
library(fasterize) 
library(ggplot2) 
library(tmap) 
library(tmaptools) 
library(gdalUtils)
library(scico)
library(ggrepel)
library(rgdal)
library(RColorBrewer)
library(BAMMtools)

# Load user defined functions
source(here::here("0_udfs.R"))

# Set customizations for ggplot theme_minimal()
theme1 = theme_minimal()+
  theme(legend.title = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())