## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script takes data from a "raw" format
## and prepare them in a disaggregated format. Questions about 
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
library(vegan)
library(lubridate)
library(stringi)
library(janitor)
library(readxl)

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
  select(-collection_date)
  
write.csv(x = nutrients_cleaned, 
          file = "../cleaned_data/nutrients.csv", 
          row.names = FALSE)

# 4. Pharmaceuticals and Personal Care Products (PPCPs) -------------------

ppcp_orig <- read.csv(file = "../cleaned_disaggregated_data/ppcp.csv", header= TRUE)

ppcp_cleaned <- ppcp_orig %>%
  filter(!(site %in% c("FI2", "DU2", "HO1"))) %>%
  mutate(time = as.numeric(time),
         sampling_event = ifelse(month == "may", "may_1", NA),
         sampling_event = ifelse(month == "june" & between(time, 1,6), "may_1", sampling_event),
         sampling_event = ifelse(month == "june" & between(time, 7, 26), "june_1", sampling_event),
         sampling_event = ifelse(month == "june" & between(time, 27, 30), "june_2", sampling_event),
         sampling_event = ifelse(month == "july" & between(time, 1,11), "june_2", sampling_event),
         sampling_event = ifelse(month == "july" & between(time, 12,24), "july_1", sampling_event),
         sampling_event = ifelse(month == "july" & between(time, 25,31), "july_2", sampling_event),
         sampling_event = ifelse(month == "august" & between(time, 1,8), "july_2", sampling_event),
         sampling_event = ifelse(month == "august" & between(time, 9,21), "august_1", sampling_event),
         sampling_event = ifelse(month == "august" & between(time, 22,31), "august_2", sampling_event),
         sampling_event = ifelse(month == "september" & between(time, 1,5), "august_2", sampling_event),
         sampling_event = ifelse(month == "september" & between(time, 6,18), "september_1", sampling_event),
         sampling_event = ifelse(month == "september" & between(time, 19,30), "september_2", sampling_event),
         sampling_event = ifelse(month == "october", "september_2", sampling_event)) %>%
  separate(col = sampling_event, into = c("month_rep", "sampling_event"), sep = "_", remove = TRUE) %>%
  mutate(sampling_event = as.factor(sampling_event)) %>%
  filter(month != "xxxx") %>%
  select(site, month, sampling_event, ppcp, concentration) %>%
  replace_na(replace = list(concentration = 0))

head(ppcp_cleaned)

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

write.csv(x = tsidw_pop,
          file = "../cleaned_data/temporally_scaled_inverse_distance_weighted_population_metrics.csv",
          row.names = FALSE)

