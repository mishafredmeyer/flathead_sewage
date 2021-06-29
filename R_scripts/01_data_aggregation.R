## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script takes data from a cleaned,
## disaggregated format prepares them for an analysis-ready format. Questions about 
##this script should be directed to 
## Michael F. Meyer (michael.f.meyer@wsu.edu). 
## The script can be divided into the following sections:
## 1. Load Packages
## 2. Ash Free Dry Mass
## 3. Nutrients
## 4. Pharmaceuticals and Personal Care Products (PPCPs)
## 5. Stoichiometry
## 6. State Park Attendance
## 7. Fatty Acids
## 8. Temporally Scaled Inverse Distance Weighted Population

# 1. Load Packages --------------------------------------------------------

library(tidyverse)
library(lubridate)
library(stringi)
library(janitor)

# 2. Ash Free Dry Mass ----------------------------------------------------

afdm_orig <- read.csv("../cleaned_disaggregated_data/afdm.csv", header = TRUE)

afdm_clean <- afdm_orig 

head(afdm_orig)

# The data appear to be properly structured for analysis. 

write.csv(x = afdm_clean, 
          file = "../cleaned_data/afdm.csv", 
          row.names = FALSE)


# 3. Nutrients ------------------------------------------------------------

nutrients_orig <- read.csv("../cleaned_disaggregated_data/nutrients.csv", header = TRUE)

nutrients_cleaned <- nutrients_orig %>%
  select(-collection_date) %>%
  pivot_longer(cols = c(nh3_n:total_p), names_to = "nutrient_type", values_to = "concentration") %>%
  mutate(tourist_season = ifelse(month %in% c("june", "july", "august"), "In Season", 
                                 "Out of Season"),
         stt = ifelse(site %in% c("HO", "FLBS", "WF", "WB"), "Centralized", "Decentralized")) %>%
  pivot_wider(names_from = "nutrient_type", values_from = "concentration")
  
write.csv(x = nutrients_cleaned, 
          file = "../cleaned_data/nutrients.csv", 
          row.names = FALSE)

# 4. Pharmaceuticals and Personal Care Products (PPCPs) -------------------

ppcp_orig <- read.csv(file = "../cleaned_disaggregated_data/ppcp.csv", header= TRUE)

head(ppcp_orig)

ppcp_cleaned <- ppcp_orig

write.csv(ppcp_cleaned, "../cleaned_data/ppcp.csv", row.names = FALSE)

# 5. Stoichiometry --------------------------------------------------------

stoich_orig <- read.csv(file = "../cleaned_disaggregated_data/stoichiometry.csv",
                        header = TRUE)

stoich_cleaned <- stoich_orig %>%
  select(site, month, total_dry_mass_carbon_mg:phosphorus_ug, collection_data_formatted)

write.csv(x = stoich_cleaned, 
          file = "../cleaned_data/stoichiometry.csv", 
          row.names = FALSE)


# 6. Fatty Acids ----------------------------------------------------------

fatty_acids_orig <- read.csv(file = "../cleaned_disaggregated_data/fatty_acids.csv", 
                             header = TRUE)

head(fatty_acids_orig)

write.csv(x = fatty_acids_orig, 
          file = "../cleaned_data/fatty_acids.csv", 
          row.names = FALSE)

# 7. Temporally Scaled Inverse Distance Weighted Population ---------------

tsidw_pop <- read.csv(file = "../cleaned_disaggregated_data/temporally_scaled_inverse_distance_weighted_population_metrics.csv",
                      header = TRUE)

head(tsidw_pop)

write.csv(x = tsidw_pop,
          file = "../cleaned_data/temporally_scaled_inverse_distance_weighted_population_metrics.csv",
          row.names = FALSE)

