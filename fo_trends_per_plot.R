######################################################################################################|
# Author: Lena Wohlgemuth

# This code analyzes ICP Forests foliage data for temporal trends individually for each Level II plot
# Data can be found in the ICP Forests database: https://www.icp-forests.org/data/fm_start.php 
# and will be provided by Pasi Rautio (Luke) upon request

######################################################################################################|
## Set working directory ####

setwd("~/Documents/esb-foliage")

######################################################################################################|
## Load packages ####

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(rkt)

######################################################################################################|
## Load functions ####

source("~/Documents/esb-foliage/03_Scripts/lm_individual_fun.R")

######################################################################################################|
## Load data ####

# Final data for analysis contains one dataset provided by Pasi Rautio and one dataset from the ICP Forests DB
# after the survey year 2009, which were merged by merge_datasets.R

dat <- read.csv2("~/Documents/esb-foliage/01_Data/20230425_joined_data_trend_analysis.csv",
                 header = T, sep = ",")

coords <- read.csv2("~/Documents/esb-foliage/01_Data/foliage_levelII_plot_coords.csv",
                    header = T, sep = ",") %>%
  rename(key_plot = site) %>% select(key_plot, latitude, longitude)
coords <- coords[!duplicated(coords[c("key_plot")]), ] #delete duplicate entries in coords

######################################################################################################|
## Add tree species information to dat ####

# Load dictionary
d_tree_spec <- read.csv2("~/Documents/esb-foliage/01_Data/d_tree_spec.csv",
                         header = T, sep = ";") %>% 
  rename(code_tree_species = code, tree_species_name = description) %>%
  select(code_tree_species, tree_species_name, grp_tree_species)

# Join dat with dictionary
dat <- left_join(dat, d_tree_spec, by = "code_tree_species")
remove(d_tree_spec)

######################################################################################################|
## Calculate ratios and contents from concentrations ####

# Filter out relevant missing values from data set
dat <- dat %>% filter(!is.na(code_tree_species))

# Make sure that all broadleaf age classes are set to 0
dat$code_leaves_type[dat$grp_tree_species == "broadleaves"] <- 0 

# Aggregate data for Quercus petraea (code: 48) and Quercus robur (code: 51)
dat$code_tree_species <- 
  ifelse(dat$code_tree_species == 51, 48, dat$code_tree_species)
# tree species code 98: Quercus petraea_or_robur
dat$code_tree_species <- 
  ifelse(dat$code_tree_species == 98, 48, dat$code_tree_species)

# Filter for tree species of interest
# for now these species are Fagus sylvatica (20), Quercus petraea & robur (48), 
# Picea abies (118), Pinus sylvestris (134), Abies alba (100)
dat <- dat %>% filter(code_tree_species %in% c(20, 48, 100, 118, 134))

# Set values to numeric, initialize columns, calculate ratios
dat <- dat %>% mutate(across(mass_100_leaves:al, as.numeric)) %>%
  mutate(key_plot = paste(code_country, partner_code, code_plot, sep = "_")) %>%
  mutate(cn = c/n, np = n/p, cp = c/p, .after = al) %>%
  mutate(n_cont = NA, c_cont = NA, p_cont = NA, k_cont = NA, s_cont = NA,
         mg_cont = NA, ca_cont = NA, .after = cp)

# Calculate nutrient contents (all contents in mg, except for C, which is g/100)
dat <- dat %>% 
  mutate(n_cont = case_when(grp_tree_species == "broadleaves" ~ n*mass_100_leaves, 
                            grp_tree_species == "conifers" ~ n*mass_1000_needles),
         c_cont = case_when(grp_tree_species == "broadleaves" ~ c*mass_100_leaves, 
                            grp_tree_species == "conifers" ~ c*mass_1000_needles),
         p_cont = case_when(grp_tree_species == "broadleaves" ~ p*mass_100_leaves, 
                            grp_tree_species == "conifers" ~ p*mass_1000_needles),
         k_cont = case_when(grp_tree_species == "broadleaves" ~ k*mass_100_leaves, 
                            grp_tree_species == "conifers" ~ k*mass_1000_needles),
         s_cont = case_when(grp_tree_species == "broadleaves" ~ s*mass_100_leaves, 
                            grp_tree_species == "conifers" ~ s*mass_1000_needles),
         mg_cont = case_when(grp_tree_species == "broadleaves" ~ mg*mass_100_leaves,
                             grp_tree_species == "conifers" ~ mg*mass_1000_needles),
         ca_cont = case_when(grp_tree_species == "broadleaves" ~ ca*mass_100_leaves,
                             grp_tree_species == "conifers" ~ ca*mass_1000_needles))


######################################################################################################|
## Calculate average value for each forest plot ####

# Explanation: Foliage values per plot and year are either derived from individual trees 
# or from a pooled sample of multiple trees (usually 5) on one plot
# -> calculate average value per plot, year, parameter, tree species, foliage age class

# Move all measurement values (param_value) to one column
dat <- dat %>% 
  pivot_longer(cols = names(dat)[which(colnames(dat) == "mass_100_leaves"):which(colnames(dat) == "ca_cont")],
               names_to = "param", values_to = "param_value")
dat$param_value <- as.numeric(dat$param_value)

# Filter out negative values and NA
dat <- dat %>% filter(!is.na(param_value), param_value > 0)

# Key for grouping before calculating average values per group
dat <- dat %>% 
  mutate(key_group_avg = paste(code_country, partner_code, code_plot, survey_year, param,
                               code_tree_species, code_leaves_type, sep = "_"))

#test_missing_leaftype <- dat %>% filter(grepl("_NA$", key_group_avg))

# Calculate average foliage value per plot, year, parameter, tree species, leaf type
dat <- dat %>% group_by(key_group_avg) %>%
  summarize(param_value_avg = mean(param_value, na.rm = T),
            param = first(param), survey_year = first(survey_year), code_tree_species = 
              first(code_tree_species), code_leaves_type = first(code_leaves_type), key_plot =
              first(key_plot)) %>% ungroup() %>% select(key_plot, survey_year, code_tree_species,
                                                        code_leaves_type, param, param_value_avg)

######################################################################################################|
## Perform an outlier correction ####

length_uncorrected <- length(dat$param_value_avg)

dat <- dat %>%
  group_by(param, code_tree_species, code_leaves_type) %>% 
  mutate(z_score = (param_value_avg - mean(param_value_avg, na.rm = T))/
           sd(param_value_avg, na.rm = T)) %>%
  ungroup() %>%
  filter(z_score > -3 & z_score < 3) #remove outliers, for which cutoff z-score is +/-3

length_outlier_removed <- length(dat$param_value_avg)

percent_outlier <- 1 - (length_outlier_removed/length_uncorrected)

remove(length_uncorrected, length_outlier_removed)

# ######################################################################################################|
# ## Insertion: what is the avg. time series length? ####
# 
# # Identify plots with the time series >= 6 years
# n_series_datapoints <- dat %>% 
#   filter(!is.na(param_value_avg)) %>%
#   filter(param != "al", param != "b") %>%
#   group_by(param, code_tree_species, code_leaves_type, key_plot) %>%
#   summarize(n_years_collected = n_distinct(survey_year), 
#             year_min = min(survey_year), year_max = max(survey_year)) %>%
#   filter(n_years_collected >= 6) #random: filter time series >= 6 years
#   
# avg_timeseries_length <- mean(n_series_datapoints$n_years_collected)
# sd_timeseries_length <- sd(n_series_datapoints$n_years_collected)

######################################################################################################|
## Plot-wise LM for different tree species and needle age classes ####

# Prepare parameters of interest
params_oi_broadleaves <- 
  data.frame(params = c("mass_100_leaves", "n", "n_cont", "c", "c_cont", "p", 
                        "p_cont", "s", "s_cont", "k", "k_cont", "cn", "np", "cp",
                        "mg", "mg_cont", "ca", "ca_cont"))
params_oi_conifers <- 
  data.frame(params = c("mass_1000_needles", "n", "n_cont", "c", "c_cont", "p", 
                        "p_cont", "s", "s_cont", "k", "k_cont", "cn", "np", "cp",
                        "mg", "mg_cont", "ca", "ca_cont"))


########## Fagus sylvatica

beech_plots_output <- lm_indiv_calc(20, 0, dat, params_oi_broadleaves) %>%
  mutate(species_name = "Fagus sylvatica", code_leaves_type = 0)

# add coordinates
beech_plots_output <- left_join(beech_plots_output, coords, by = "key_plot")

# save to rdata
#saveRDS(beech_plots_output, file = "~/Documents/esb-foliage/01_Data/beech_y0_lm.rds")


########## Quercus petraea & robur

oak_plots_output <- lm_indiv_calc(48, 0, dat, params_oi_broadleaves) %>%
  mutate(species_name = "Quercus petraea & robur", code_leaves_type = 0)

# add coordinates
oak_plots_output <- left_join(oak_plots_output, coords, by = "key_plot")


########## Abies alba

#### current

fir_y0_plots_output <- lm_indiv_calc(100, 0, dat, params_oi_conifers) %>%
  mutate(species_name = "Abies alba", code_leaves_type = 0)

# add coordinates
fir_y0_plots_output <- left_join(fir_y0_plots_output, coords, by = "key_plot")

#### c + 1

fir_y1_plots_output <- lm_indiv_calc(100, 1, dat, params_oi_conifers) %>%
  mutate(species_name = "Abies alba", code_leaves_type = 1)

# add coordinates
fir_y1_plots_output <- left_join(fir_y1_plots_output, coords, by = "key_plot")


########## Picea abies

#### current

spruce_y0_plots_output <- lm_indiv_calc(118, 0, dat, params_oi_conifers) %>%
  mutate(species_name = "Picea abies", code_leaves_type = 0)

# add coordinates
spruce_y0_plots_output <- left_join(spruce_y0_plots_output, coords, by = "key_plot")

#### c + 1

spruce_y1_plots_output <- lm_indiv_calc(118, 1, dat, params_oi_conifers) %>%
  mutate(species_name = "Picea abies", code_leaves_type = 1)

# add coordinates
spruce_y1_plots_output <- left_join(spruce_y1_plots_output, coords, by = "key_plot")


########## Pinus sylvestris

#### current

pine_y0_plots_output <- lm_indiv_calc(134, 0, dat, params_oi_conifers) %>%
  mutate(species_name = "Pinus sylvestris", code_leaves_type = 0)

# add coordinates
pine_y0_plots_output <- left_join(pine_y0_plots_output, coords, by = "key_plot")

#### c + 1

pine_y1_plots_output <- lm_indiv_calc(134, 1, dat, params_oi_conifers) %>%
  mutate(species_name = "Pinus sylvestris", code_leaves_type = 1)

# add coordinates
pine_y1_plots_output <- left_join(pine_y1_plots_output, coords, by = "key_plot")


######################################################################################################|
## Dataset containing all tree species ####

dat_lm_plots_trends <- 
  rbind(beech_plots_output, oak_plots_output, fir_y0_plots_output,
        fir_y1_plots_output, spruce_y0_plots_output, spruce_y1_plots_output,
        pine_y0_plots_output, pine_y1_plots_output)

# save to rdata
#saveRDS(dat_lm_plots_trends, file = "~/Documents/esb-foliage/01_Data/20230814_trends_species_lm.rds")





