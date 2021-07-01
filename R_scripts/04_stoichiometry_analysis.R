## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script is a piece of the 
## overarching analysis and looks mostly at relating 
## periphyton Carbon, Nitrogen, and Phosphorus stoichiometries
## with wastewater treatment infrastructure
## and a sampling point's temporal position to the tourist 
## season. Questions about this script should be directed to 
## Michael F. Meyer (michael.f.meyer@wsu.edu). 


# 1. Load packages and data -----------------------------------------------

library(tidyverse)
library(vegan)
library(ggrepel)
library(viridis)
library(vegan)
library(car)
library(ggpubr)

stoichiometry_orig <- read.csv(file = "../cleaned_data/stoichiometry.csv")



# 2. Convert mass to moles ------------------------------------------------

stoich <- stoichiometry_orig %>%
  mutate(carbon_moles = (carbon_mg*0.001)/12.0107,
         carbon = carbon_moles/(total_dry_mass_carbon_mg*0.001),
         nitrogen_moles = (nitrogen_mg*0.001)/14.0067,
         nitrogen = nitrogen_moles/(total_dry_mass_nitrogen_mg*0.001),
         phosphorus_mg = phosphorus_ug * 0.001,
         phosphorus_moles = (phosphorus_mg*0.001)/39.973762,
         phosphorus = phosphorus_moles/(total_dry_mass_phosphorus_mg*0.001)) %>%
  select(site, month, stt, tourist_season, carbon, nitrogen, phosphorus)


# 3. Create stoichiometric plots and analyses -----------------------------

carbon_nitrogen <- ggplot(stoich, aes(tourist_season, carbon/nitrogen)) +
  geom_boxplot(width = 0.2, outlier.alpha = 0) +
  geom_jitter(size = 5) +
  facet_wrap(~stt) +
  theme_bw()

carbon_phosphorus <- ggplot(stoich, aes(tourist_season, carbon/phosphorus)) +
  geom_boxplot(width = 0.2, outlier.alpha = 0) +
  geom_jitter(size = 5) +
  facet_wrap(~stt) +
  theme_bw()

nitrogen_phosphorus <- ggplot(stoich, aes(tourist_season, nitrogen/phosphorus)) +
  geom_boxplot(width = 0.2, outlier.alpha = 0) +
  geom_jitter(size = 5) +
  facet_wrap(~stt) +
  theme_bw()
