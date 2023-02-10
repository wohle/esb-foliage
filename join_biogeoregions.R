######################################################################################################|

# This code assigns EEA biogeographical regions to plots of the ICP Forests foliage survey
# EEA biogeographical regions: https://www.eea.europa.eu/data-and-maps/data/biogeographical-regions-europe-3
# The foliage survey is included in the ICP Forests database: https://www.icp-forests.org/data/fm_start.php

######################################################################################################|
## Set working directory ####

setwd("~/Documents/esb-foliage")

######################################################################################################|
## Load packages ####

library(dplyr)
library(rgdal)
library(sf)
library(sp)

######################################################################################################
## Load data ####

plf_f0_raw <- read.table("~/Documents/esb-foliage/01_Data/440_fo_20230206132239/fo_plf.csv",
                         header = T, sep = ";", 
                         colClasses = list(latitude = "character", longitude = "character"))
plf_f1_raw <- read.table("~/Documents/esb-foliage/01_Data/440_f1_20230206144750/f1_plf.csv",
                         header = T, sep = ";", 
                         colClasses = list(latitude = "character", longitude = "character"))
ecozones <- st_read("~/Documents/esb-foliage/01_Data/BiogeoRegions2016_shapefile/BiogeoRegions2016.shp")


######################################################################################################|
## Prepare datasets ####

plf_f0 <- plf_f0_raw %>% mutate(site = paste(code_country, partner_code,
                                             code_plot, sep = "_"))
plf_f1 <- plf_f1_raw %>% mutate(site = paste(code_country, partner_code,
                                             code_plot, sep = "_"))

coords_f0 <- plf_f0 %>% select(site, longitude, latitude)
coords_f1 <- plf_f1 %>% select(site, longitude, latitude)

coords <- rbind(coords_f0, coords_f1)
coords <- coords %>% distinct() %>% rename(lon_char = longitude,
                                           lat_char = latitude)

remove(coords_f0, coords_f1)

# convert ICP Forests coordinates (D°M’S“) to decimal degrees
coords <- coords %>% mutate(longitude = as.numeric(substr(lon_char, 1, 3)) +
                              as.numeric(substr(lon_char, 4, 5))/60 +
                              as.numeric(substr(lon_char, 6, 7))/3600,
                            latitude = as.numeric(substr(lat_char, 1, 3)) +
                              as.numeric(substr(lat_char, 4, 5))/60 +
                              as.numeric(substr(lat_char, 6, 7))/3600) %>%
  select(site, longitude, latitude)

# convert dataframe to spatial object (crs: 4326 is WGS84)
coords_sf <- st_as_sf(coords, coords = c("longitude", "latitude"), 
                      crs = 4326) 

# transform ecozone coordinates to decimal degrees 
ecozones_dd <- ecozones %>% st_transform(4326)


######################################################################################################|
## Overlay spatial objects ####

dat <- st_join(coords_sf, ecozones_dd[c("short_name", "name", "code")], 
               join = st_intersects, left = T)

######################################################################################################|
## Solve issue with coastal forest plots ####

# Data issue: ICP Forests coordinates from the ICP Forests database are round to
# minutes, which is why coordinates of plots are not exact
# round coordinates of plots close to the coast sometimes point to the sea where
# the EEA biozones are not defined
# -> find closest biozone polygon for coastal forest plots

# plots for which problem exists
dat_coastal <- dat %>% filter(is.na(short_name))

# create sf object for coastal plots
sites_coastal_sf <- coords_sf %>% filter(site %in% dat_coastal$site)

# filter out missing values from data
dat <- dat %>% filter(!is.na(short_name))

remove(dat_coastal)

# find nearest polygon of ecozones to coastal sites
# requires some computing power!
nearest <- ecozones_dd[st_nearest_feature(sites_coastal_sf, ecozones_dd),]  %>% 
  dplyr::mutate(site = sites_coastal_sf$site)
nearest

plots_coastal_biogeo <- nearest %>% select(site, short_name, name, code)
plots_coastal_biogeo$geometry <- NULL

######################################################################################################|
## Join coastal and non-coastal plots ####

dat <- dat %>% select(site, short_name, name, code)
dat$geometry <- NULL

dat <- rbind(dat, plots_coastal_biogeo)

#write.csv(dat, "Foliage_plots_biogeozones.csv", row.names = F)
