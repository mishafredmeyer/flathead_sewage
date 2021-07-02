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
library(lubridate)
library(stringi)
library(janitor)
library(readxl)
library(sf)
library(OpenStreetMap)

# 2. Ash Free Dry Mass ----------------------------------------------------

afdm <- read.csv("../raw_data/AFDM_final_mfm_20180810.csv", header = TRUE)

afdm_clean <- afdm %>%
  clean_names() %>%
  mutate(site = trimws(site),
         site = ifelse(site == "Beardance", "BD", site),
         site = ifelse(site == "Blue Bay", "BB", site),
         site = ifelse(site == "Boettcher", "BO", site),
         site = ifelse(site == "Dayton", "DA", site),
         site = ifelse(site == "Finley Point", "FI", site),
         site = ifelse(site == "FLBS", "FLBS", site),
         site = ifelse(site == "Lakeside", "LK", site),
         site = ifelse(site == "Sacajawea", "SJ", site),
         site = ifelse(site == "Salish Park", "SA", site),
         site = ifelse(site == "Wayfarers", "WF", site),
         site = ifelse(site == "West Shore", "WS", site),
         site = ifelse(site == "Woods Bay", "WB", site),
         site = ifelse(site == "Yellow Bay", "YB", site),
         month = tolower(month)) %>%
  mutate(afdm = dry_weight - after_ignition,
         stt = ifelse(test = site %in% c("WF", "YB", "FLBS", "HO"), 
                      yes = "centralized", no = "decentralized"),
         tourist_season = ifelse(test = month %in% c("june", "july", "august"), 
                                 yes = "In Season", no = "Out of Season"))

write.csv(x = afdm_clean, 
          file = "../cleaned_disaggregated_data/afdm.csv", 
          row.names = FALSE)


# 3. Nutrients ------------------------------------------------------------

nutrients <- read.csv("../raw_data/nutrients.csv", header = TRUE)

nutrients_cleaned <- nutrients[-c(1:2) , ] %>%
  clean_names() %>%
  select(-contains("X"), -analyte) %>% 
  pivot_longer(cols = c(nh3_n:total_p), 
               names_to = "nutrient_type", 
               values_to = "concentration") %>%
  mutate(concentration = ifelse(test = grepl(pattern = "<", x = concentration), 
                                yes = as.numeric(trimws(gsub(pattern = "<", replacement = "", x = concentration)))/2,
                                no = concentration)) %>%
  pivot_wider(names_from = "nutrient_type", values_from = "concentration") %>%
  mutate(site = trimws(sample_id)) %>%
  separate(sample_id, c("site", "rep"), sep = "-") %>%
  mutate(site = ifelse(site == "Beardance", "BD", site),
         site = ifelse(site == "Blue Bay", "BB", site),
         site = ifelse(site == "Boettcher", "BO", site),
         site = ifelse(site == "Dayton", "DA", site),
         site = ifelse(site == "Finley Point", "FI", site),
         site = ifelse(site == "FLBS", "FLBS", site),
         site = ifelse(site == "Lakeside", "LK", site),
         site = ifelse(site == "Sacajawea", "SJ", site),
         site = ifelse(site == "Salish Park", "SA", site),
         site = ifelse(site == "Wayfarers", "WF", site),
         site = ifelse(site == "West Shore", "WS", site),
         site = ifelse(site == "Woods Bay", "WB", site),
         site = ifelse(site == "Yellow Bay", "YB", site),
         site = ifelse(site == "Holt", "HO", site),
         site = ifelse(site == "Ducharme", "DU", site),
         month = month(mdy(collection_date), label = TRUE, abbr = FALSE),
         month = tolower(as.character(month)),
         collection_data_formatted = paste0(year(mdy(collection_date)), 0,
                                            month(mdy(collection_date)),
                                            day(mdy(collection_date))))

write.csv(x = nutrients_cleaned, 
          file = "../cleaned_disaggregated_data/nutrients.csv", 
          row.names = FALSE)

# 4. Pharmaceuticals and Personal Care Products (PPCPs) -------------------

ppcp_outputs <- list.files(pattern=".CSV",
                          path = "../raw_data/hplc_raw_output/",
                          recursive = TRUE)
ppcp_outputs <- paste0("../raw_data/hplc_raw_output/", ppcp_outputs)

dat_full <- data.frame()

for(i in 1:length(ppcp_outputs)){
  if(ppcp_outputs[i] != "../raw_data/hplc_raw_output/Flathead_PPCP_20180904/STD_200UG_LINEA_20180830_ROW4.CSV") {
    temp <- read.csv(paste(ppcp_outputs[i]), header = FALSE, fileEncoding= "UTF-16LE", sep = "\t")
  }
  else {
    temp <- read.csv(paste(ppcp_outputs[i]), header = FALSE)
  }
  temp$FILE <- paste(ppcp_outputs[i])
  dat_full <- rbind(dat_full, temp)
}
colnames(dat_full) <- c("PEAK_TIME", "AREA", "FILE")

dat_full$FILE <- gsub(pattern = "../raw_data/hplc_raw_output/", replacement = "", dat_full$FILE)
dat_full$FILE <- gsub(pattern = "Flathead_PPCP_20180920/", replacement = "", dat_full$FILE)
dat_full$FILE <- gsub(pattern = "Flathead_PPCP_20180904/", replacement = "", dat_full$FILE)
dat_full$FILE <- gsub(pattern = "Flathead_PPCP_20181001/", replacement = "", dat_full$FILE)
dat_full$FILE <- gsub(pattern = "Flathead_PPCP_20181008/", replacement = "", dat_full$FILE)

samples <- dat_full %>%
  filter(!grepl("STD", FILE) & !grepl("BLANK", FILE)) %>%
  mutate(FILE = gsub("\\.CSV", "", FILE)) %>%
  separate(FILE, into = c("SITE", "MONTH", "TIME",  "LINE", "DATE_HPLC", "ROW"), sep = "_") %>%
  mutate(PPCP = ifelse(ROW == "ROW2" & PEAK_TIME < 10, "acetaminophen", NA), 
         PPCP = ifelse(ROW == "ROW2" & PEAK_TIME > 10, "diphenhydramine", PPCP),
         PPCP = ifelse(ROW == "ROW3" & PEAK_TIME < 10, "cotinine", PPCP),
         PPCP = ifelse(ROW == "ROW3" & PEAK_TIME > 10, "carbamezapine", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME < 5, "paraxanthine", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME > 5 & PEAK_TIME < 10, "caffeine", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME > 10 & PEAK_TIME < 20, "sulfamexthazole", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME > 20 & PEAK_TIME < 30, "warfarin", PPCP))
samples$PPCP <- as.factor(samples$PPCP)
samples$SITE <- as.factor(samples$SITE)
samples$MONTH <- as.factor(samples$MONTH)
#samples$TIME <- as.factor(samples$TIME)
samples$LINE <- as.factor(samples$LINE)
samples$DATE_HPLC <- as.factor(samples$DATE_HPLC)
samples$ROW <- as.factor(samples$ROW)

blanks <- dat_full %>%
  filter(grepl("BLANK", FILE)) %>%
  mutate(FILE = gsub("\\.CSV", "", FILE),
         FILE = gsub("\\.txt", "", FILE)) %>%
  separate(FILE, into = c("SAMPLE", "LINE", "DATE_HPLC", "ROW"), sep = "_") %>%
  mutate(PPCP = ifelse(ROW == "ROW2" & PEAK_TIME < 10, "acetaminophen", NA), 
         PPCP = ifelse(ROW == "ROW2" & PEAK_TIME > 10, "diphenhydramine", PPCP),
         PPCP = ifelse(ROW == "ROW3" & PEAK_TIME < 10, "cotinine", PPCP),
         PPCP = ifelse(ROW == "ROW3" & PEAK_TIME > 10, "carbamezapine", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME < 5, "paraxanthine", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME > 5 & PEAK_TIME < 10, "caffeine", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME > 10 & PEAK_TIME < 20, "sulfamexthazole", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME > 20 & PEAK_TIME < 30, "warfarin", PPCP))

blanks$PPCP <- as.factor(blanks$PPCP)
blanks$SAMPLE <- as.factor(blanks$SAMPLE)
blanks$LINE <- as.factor(blanks$LINE)
blanks$DATE_HPLC <- as.factor(blanks$DATE_HPLC)
blanks$ROW <- as.factor(blanks$ROW)

blanks.join <- select(blanks, PPCP, LINE, AREA)
samples <- left_join(samples, blanks.join, by = c("PPCP", "LINE")) %>%
  mutate(AREA.y = ifelse(is.na(AREA.y), 0, AREA.y),
         AREA = AREA.x - AREA.y)

standards <- dat_full %>%
  filter(grepl("STD", FILE)) %>%
  mutate(FILE = gsub("\\.CSV", "", FILE)) %>%
  separate(FILE, into = c("SAMPLE", "CONCENTRATION", "LINE", "DATE_HPLC", "ROW"), sep = "_") %>%
  mutate(PPCP = ifelse(ROW == "ROW2" & PEAK_TIME < 10, "acetaminophen", NA), 
         PPCP = ifelse(ROW == "ROW2" & PEAK_TIME > 10, "diphenhydramine", PPCP),
         PPCP = ifelse(ROW == "ROW3" & PEAK_TIME < 10, "cotinine", PPCP),
         PPCP = ifelse(ROW == "ROW3" & PEAK_TIME > 10, "carbamezapine", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME < 5, "paraxanthine", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME > 5 & PEAK_TIME < 10, "caffeine", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME > 10 & PEAK_TIME < 20, "sulfamexthazole", PPCP),
         PPCP = ifelse(ROW == "ROW4" & PEAK_TIME > 20 & PEAK_TIME < 30, "warfarin", PPCP)) %>%
  filter(!is.na(PPCP)) %>%
  mutate(CONCENTRATION_FIX = ifelse(CONCENTRATION == "20NG", "2NG", CONCENTRATION),
         CONCENTRATION_FIX = ifelse(CONCENTRATION == "200NG", "20NG", CONCENTRATION_FIX),
         CONCENTRATION_FIX = ifelse(CONCENTRATION == "2NG", "200NG", CONCENTRATION_FIX)) %>%
  separate(CONCENTRATION_FIX, into = c("CONCENTRATION_FIX", "UNITS"), sep = -2, remove = TRUE) %>%
  mutate(CONCENTRATION_FIX = ifelse(CONCENTRATION_FIX == 2 & UNITS == "UG", 2000, CONCENTRATION_FIX),
         CONCENTRATION_FIX = ifelse(CONCENTRATION_FIX == 20 & UNITS == "UG", 20000, CONCENTRATION_FIX),
         CONCENTRATION_FIX = ifelse(CONCENTRATION_FIX == 200 & UNITS == "UG", 200000, CONCENTRATION_FIX),
         CONCENTRATION_FIX = ifelse(CONCENTRATION_FIX == 2 & UNITS == "MG", 2000000, CONCENTRATION_FIX)) %>%
  filter(CONCENTRATION_FIX != 200000)

standards$PPCP <- as.factor(standards$PPCP)
standards$SAMPLE <- as.factor(standards$SAMPLE)
standards$LINE <- as.factor(standards$LINE)
standards$DATE_HPLC <- as.factor(standards$DATE_HPLC)
standards$ROW <- as.factor(standards$ROW)
standards$UNITS <- as.factor(standards$UNITS)
standards$CONCENTRATION_FIX <- as.numeric(standards$CONCENTRATION_FIX)

standards %>%
  filter(CONCENTRATION_FIX != 200000) %>%
  filter(!(PPCP == "warfarin" & log10(AREA) <= 5)) %>%
  ggplot(aes(log10(CONCENTRATION_FIX), log10(AREA))) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ PPCP)

aceta.lm <- lm(log10(CONCENTRATION_FIX) ~ log10(AREA) , data = standards[standards$PPCP == "acetaminophen",])
summary(aceta.lm)

caff.lm <- lm(log10(CONCENTRATION_FIX) ~ log10(AREA), data = standards[standards$PPCP == "caffeine",])
summary(caff.lm)

carb.lm <- lm(log10(CONCENTRATION_FIX) ~ log10(AREA), data = standards[standards$PPCP == "carbamezapine",])
summary(carb.lm)

cot.lm <- lm(log10(CONCENTRATION_FIX) ~ log10(AREA), data = standards[standards$PPCP == "cotinine",])
summary(cot.lm)

diph.lm <- lm(log10(CONCENTRATION_FIX) ~ log10(AREA), data = standards[standards$PPCP == "diphenhydramine",])
summary(diph.lm)

para.lm <- lm(log10(CONCENTRATION_FIX) ~ log10(AREA), data = standards[standards$PPCP == "paraxanthine",])
summary(para.lm)

sulfa.lm <- lm(log10(CONCENTRATION_FIX) ~ log10(AREA), data = standards[standards$PPCP == "sulfamexthazole",])
summary(sulfa.lm)

warf.lm <- lm(log10(CONCENTRATION_FIX) ~ log10(AREA), data = standards[standards$PPCP == "warfarin" &
                                                                         log10(standards$AREA) >= 5,])
summary(warf.lm)


aceta_data <- filter(samples, PPCP == "acetaminophen")
new <- data.frame(AREA = aceta_data$AREA)
aceta_data$CONCENTRATION <- predict.lm(aceta.lm, newdata = new)
aceta_data$CONCENTRATION <- 10^(aceta_data$CONCENTRATION)

caff_data <- filter(samples, PPCP == "caffeine")
new <- data.frame(AREA = caff_data$AREA)
caff_data$CONCENTRATION <- predict.lm(caff.lm, newdata = new)
caff_data$CONCENTRATION <- 10^(caff_data$CONCENTRATION)

carb_data <- filter(samples, PPCP == "carbamezapine")
new <- data.frame(AREA = carb_data$AREA)
carb_data$CONCENTRATION <- predict.lm(carb.lm, newdata = new)
carb_data$CONCENTRATION <- 10^(carb_data$CONCENTRATION)

cot_data <- filter(samples, PPCP == "cotinine")
new <- data.frame(AREA = cot_data$AREA)
cot_data$CONCENTRATION <- predict.lm(cot.lm, newdata = new)
cot_data$CONCENTRATION <- 10^(cot_data$CONCENTRATION)

diph_data <- filter(samples, PPCP == "diphenhdyramine")
new <- data.frame(AREA = diph_data$AREA)
diph_data$CONCENTRATION <- predict.lm(diph.lm, newdata = new)
diph_data$CONCENTRATION <- 10^(diph_data$CONCENTRATION)

para_data <- filter(samples, PPCP == "paraxanthine")
new <- data.frame(AREA = para_data$AREA)
para_data$CONCENTRATION <- predict.lm(para.lm, newdata = new)
para_data$CONCENTRATION <- 10^(para_data$CONCENTRATION)

sulfa_data <- filter(samples, PPCP == "sulfamexathazole")
new <- data.frame(AREA = sulfa_data$AREA)
sulfa_data$CONCENTRATION <- predict.lm(sulfa.lm, newdata = new)
sulfa_data$CONCENTRATION <- 10^(sulfa_data$CONCENTRATION)

warf_data <- filter(samples, PPCP == "warfarin")
new <- data.frame(AREA = warf_data$AREA)
warf_data$CONCENTRATION <- predict.lm(warf.lm, newdata = new)
warf_data$CONCENTRATION <- 10^(warf_data$CONCENTRATION)

samples_full <- rbind(aceta_data, caff_data, carb_data, cot_data, diph_data, para_data, sulfa_data, warf_data)
summary(samples_full)

samples_formatted <- samples_full %>%
  mutate(MONTH = as.character(MONTH),
         MONTH = trimws(tolower(MONTH)),
         MONTH = ifelse(MONTH == "juny", "june", MONTH),
         MONTH = ifelse(MONTH == "sept", "september", MONTH),
         #MONTH = ifelse(MONTH == "october", "september", MONTH),
         SITE = as.character(SITE),
         SITE = ifelse(SITE == "SA", "SL", SITE)) %>%
  clean_names(case = "snake") %>%
  rename("area_observed" = "area_x",
         "area_blank" = "area_y",
         "area_corrected" = "area") %>%
  filter(!(site %in% c("FI2", "DU2", "HO1"))) %>%
  filter(month != "xxxx") %>%
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
  mutate(peri_sampling = ifelse(sampling_event %in% c("june_1", "june_2"), "peri_1", NA),
         peri_sampling = ifelse(sampling_event %in% c("july_1", "july_2"), "peri_2", peri_sampling),
         peri_sampling = ifelse(sampling_event %in% c("august_1", "august_2"), "peri_3", peri_sampling),
         peri_sampling = ifelse(sampling_event %in% c("september_1", "september_2"), "peri_4", peri_sampling)) %>%
  replace_na(list(concentration = 0)) %>%
  #select(-month) %>%
  #separate(col = "sampling_event", into = c("month", "sampling_time"), remove = TRUE) %>%
  mutate(stt = ifelse(site %in% c("WF", "HO", "FLBS", "YB"), "Centralized", "Decentralized")) %>%
  select(site, sampling_event, ppcp, concentration, peri_sampling, month, time) %>%
  mutate(stt = ifelse(test = site %in% c("WF", "YB", "FLBS", "HO"), 
                      yes = "centralized", no = "decentralized"),
         tourist_season = ifelse(peri_sampling %in% c("peri_4"), "Out of Season", "In Season"))
  

write.csv(samples_formatted, "../cleaned_disaggregated_data/ppcp.csv", row.names = FALSE)

# 5. Stoichiometry -------------------------------------------------------

stoich <- read.csv("../raw_data/biofilm_stoichiometry_20171219.csv", 
                   header = TRUE, stringsAsFactors = FALSE)

stoich_cleaned <- stoich[-c(1:3),] %>%
  rename("SITE_full" = "X",
         "date_collected" = "X.1",
         "total_dry_mass_carbon_mg" = "Particulate.Carbon",
         "carbon_mg" = "X.3",
         "total_dry_mass_nitrogen_mg" = "Particulate.Nitrogen",
         "nitrogen_mg" = "X.4",
         "total_dry_mass_phosphorus_mg" = "Particulate.Phosphorus",
         "phosphorus_ug" = "X.5") %>%
  select(-contains("X")) %>%
  mutate(MONTH = tolower(month(mdy(date_collected), label = TRUE, abbr = FALSE)),
         SITE = unlist(str_extract_all(SITE_full,  "(?<=\\().+?(?=\\))"))) %>%
  pivot_longer(cols = c(total_dry_mass_carbon_mg:phosphorus_ug), 
               names_to = "nutrient_type", 
               values_to = "mass") %>%
  mutate(mass = ifelse(test = grepl(pattern = "<", x = mass), 
                                yes = as.numeric(trimws(gsub(pattern = "<", replacement = "", x = mass)))/2,
                                no = mass)) %>%
  pivot_wider(names_from = "nutrient_type", values_from = "mass") %>%
  mutate(across(c("total_dry_mass_carbon_mg", "carbon_mg",
                      "total_dry_mass_nitrogen_mg", "nitrogen_mg",
                      "total_dry_mass_phosphorus_mg", "phosphorus_ug"), 
            as.numeric)) %>%
  select(-SITE_full) %>%
  mutate(collection_data_formatted = paste0(year(mdy(date_collected)), 0,
                                            month(mdy(date_collected)))) %>%
  clean_names(case = "snake") %>%
  mutate(stt = ifelse(test = site %in% c("WF", "YB", "FLBS", "HO"), 
                      yes = "centralized", no = "decentralized"),
         tourist_season = ifelse(test = month %in% c("june", "july", "august"), 
                                 yes = "In Season", no = "Out of Season"))

write.csv(stoich_cleaned, "../cleaned_disaggregated_data/stoichiometry.csv", row.names = FALSE)

# 6. State Park Attendance ------------------------------------------------

state_park_attendance <- read.csv("../raw_data/state_park_attendance.csv", header = TRUE)

state_park <- state_park_attendance %>%
  rename("FI" = "finley_point",
         "WF" = "wayfarers",
         "WS" = "west_shore",
         "YB" = "yellow_bay",
         "MONTH" = "Month") %>%
  select(-big_arm) %>%
  gather(SITE, visitors, FI:YB) %>%
  mutate(MONTH = toupper(MONTH))

write.csv(x = state_park, 
          file = "../cleaned_disaggregated_data/state_park_attendance.csv", 
          row.names = FALSE)


# 7. Fatty Acids ----------------------------------------------------------

july_fatty_acids <- read_excel(path = "../raw_data/Flathead_FA_orig_2018April.xlsx", 
                               sheet = "July2017", trim_ws = TRUE, col_names = TRUE) %>%
  row_to_names(row_number = 1) %>%
  clean_names() %>%
  separate(col = na, into = c("site", "month", "id", "sample_number"), sep = "_") %>%
  filter(month == "07") %>%
  select(site, month, c12_0:c28_0) %>%
  mutate(month = "july", 
         tourist_season = ifelse(month %in% c("june", "july", "august"), "In Season", 
                                 "Out of Season"),
         stt = ifelse(site %in% c("HO", "FLBS", "WF", "YB"), "Centralized", "Decentralized"),
         c22_1w9c = ifelse(c22_1w9c == "-----", 0, c22_1w9c))

august_fatty_acids <- read_excel(path = "../raw_data/Flathead_FA_orig_2018April.xlsx", 
                               sheet = "August2017", trim_ws = TRUE, col_names = TRUE) %>%
  row_to_names(row_number = 1) %>%
  clean_names() %>%
  separate(col = na, into = c("site", "month", "id", "sample_number"), sep = "_") %>%
  filter(month == "08") %>%
  select(site, month, c12_0:c28_0) %>%
  mutate(month = "august", 
         tourist_season = ifelse(month %in% c("june", "july", "august"), "In Season", 
                                 "Out of Season"),
         stt = ifelse(site %in% c("HO", "FLBS", "WF", "YB"), "Centralized", "Decentralized"),
         c22_1w9c = ifelse(c22_1w9c == "-----", 0, c22_1w9c))

september_fatty_acids <- read_excel(path = "../raw_data/Flathead_FA_orig_2018April.xlsx", 
                                 sheet = "Sept2017", trim_ws = TRUE, col_names = TRUE) %>%
  row_to_names(row_number = 1) %>%
  clean_names() %>%
  separate(col = na, into = c("site", "month", "id", "sample_number"), sep = "_") %>%
  filter(month == "09") %>%
  select(site, month, c12_0:c28_0) %>%
  mutate(month = "september", 
         tourist_season = ifelse(month %in% c("june", "july", "august"), "In Season", 
                                 "Out of Season"),
         stt = ifelse(site %in% c("HO", "FLBS", "WF", "YB"), "Centralized", "Decentralized"),
         c22_1w9c = ifelse(c22_1w9c == "-----", 0, c22_1w9c))

fatty_acids <- rbind(july_fatty_acids, august_fatty_acids, september_fatty_acids)

write.csv(x = fatty_acids, 
          file = "../cleaned_disaggregated_data/fatty_acids.csv", 
          row.names = FALSE)


# 8. Temporally Scaled Inverse Distance Weighted Population ---------------

# Load the shapefiles and metadata containing information of the development
# polygons and the lake shorelines.

flathead_shapefile <- sf::st_read(dsn = "../raw_data/Flathead sampling locations.kml")

sampling_loc <- c("LK", "WS", "DA", "SJ", "SL", "BO", "DU", 
                  "FI", "BB", "FLBS", "YB", "BD", "WB", "WF", "HO")

# Extract and format development area data 
loc_areas <- flathead_shapefile %>%
  filter(grepl("shapefile", Name)) %>%
  st_area() %>%
  enframe(name = NULL) %>%
  rename(development_area_m2 = value) %>%
  mutate(development_area_m2 = as.numeric(development_area_m2),
         development_area_km2 = development_area_m2 / 1000000) %>%
  cbind(flathead_shapefile %>% as_tibble() %>%
          filter(grepl("shapefile", Name)) %>%
          mutate(Name = gsub(pattern = "_shapefile", 
                             replacement = "", x = Name),
                 Name = gsub(pattern = "_", 
                             replacement = " ", x = Name)) %>%
          dplyr::select(Name) %>%
          as.vector()) %>%
  rename("development" = Name)

loc_areas_centroids <- flathead_shapefile %>% 
  filter(grepl("shapefile", Name)) %>%
  st_centroid() %>%
  as_tibble() %>%
  mutate(Name = gsub(pattern = "_shapefile", 
                     replacement = "", x = Name),
         Name = gsub(pattern = "_", 
                     replacement = " ", x = Name)) %>%
  dplyr::select(development = Name, centroid_loc = geometry) %>%
  inner_join(. , y = loc_areas)

# Extract shoreline length info and join with area
loc_shoreline_area_length <- flathead_shapefile %>%
  filter(grepl("shoreline", Name)) %>%
  st_length() %>%
  as_tibble() %>%
  rename(development_shoreline_length_m = value) %>%
  mutate(development_shoreline_length_m = as.numeric(development_shoreline_length_m),
         development_shoreline_length_km = development_shoreline_length_m / 1000) %>%
  cbind(flathead_shapefile %>% as_tibble() %>%
          filter(grepl("shoreline", Name)) %>%
          mutate(Name = gsub(pattern = "_shoreline", 
                             replacement = "", x = Name),
                 Name = gsub(pattern = "_", 
                             replacement = " ", x = Name)) %>%
          dplyr::select(Name) %>%
          as.vector()) %>%
  rename(development = Name) %>%
  full_join(x = ., y = loc_areas_centroids, by = c("development")) %>%
  dplyr::select(development, centroid_loc, development_shoreline_length_km, development_area_km2) %>%
  mutate(neighbor = gsub(pattern = " ", replacement = "_", development),
         neighbor = tolower(neighbor))

# Convert metadata into a spatial object 
site_loc <- flathead_shapefile %>%
  filter(Name %in% sampling_loc) 


# Join all data together and format into a dataframe
locs_centroids <- st_distance(x = site_loc,
                              y = loc_areas_centroids$centroid_loc) %>%
  as_tibble() %>%
  rename("Bigfork" = V1,
         "Woods Bay" = V2,
         "Bear Dance" = V3,
         "FLBS" = V4,
         "Yellow Bay" = V5,
         "Blue Bay" = V6,
         "Finley" = V7,
         "Polson" = V8,
         "Dayton" = V9,
         "West Shore" = V10,
         "Lakeside" = V11) %>%
  clean_names() %>%
  cbind(., site_loc$Name) %>%
  as_tibble() %>%
  gather(key = neighbor, value = distance, bigfork:lakeside) %>%
  rename("sampling_site" = `site_loc$Name`) %>%
  mutate(distance_km = distance / 1000,
         population = ifelse(test = neighbor == "bigfork",
                             yes = 4797, no = NA),
         population = ifelse(test = neighbor == "woods_bay",
                             yes = 688, no = population),
         population = ifelse(test = neighbor == "bear_dance",
                             yes = 212, no = population),
         population = ifelse(test = neighbor == "flbs",
                             yes = 90, no = population),
         population = ifelse(test = neighbor == "yellow_bay",
                             yes = 90, no = population),
         population = ifelse(test = neighbor == "blue_bay",
                             yes = 110, no = population),
         population = ifelse(test = neighbor == "finley",
                             yes = 453, no = population),
         population = ifelse(test = neighbor == "polson",
                             yes = 4702, no = population),
         population = ifelse(test = neighbor == "dayton",
                             yes = 112, no = population),
         population = ifelse(test = neighbor == "west_shore",
                             yes = 66, no = population),
         population = ifelse(test = neighbor == "lakeside",
                             yes = 2459, no = population)) %>%
  filter(!(neighbor %in% c("bear_dance", "flbs", "blue_bay", "finley", "west_shore", "yellow_bay"))) %>%
  left_join(x = ., y = loc_shoreline_area_length) %>%
  mutate(distance_weighted_population = ((population * development_shoreline_length_km) /
                                           development_area_km2) / distance_km) %>%
  group_by(sampling_site) %>%
  summarize(distance_weighted_population = (sum(as.numeric(distance_weighted_population)))) %>%
  arrange(distance_weighted_population)

write.csv(x = locs_centroids,
          file = "../cleaned_disaggregated_data/static_distance_weighted_population_metrics.csv",
          row.names = FALSE)

### Calculate temporal scaling factors

### Start with the park admittance data

park_data_orig <- read.csv("../raw_data/state_park_attendance.csv")

## We will standardize the scaling metric by May attendance, because this is likely 
## before the tourism season started 

park_scalars <- park_data_orig %>%
  pivot_longer(cols = big_arm:yellow_bay, names_to = "state_park", 
               values_to = "admitted_persons") %>%
  pivot_wider(names_from = Month, values_from = admitted_persons) %>%
  mutate(june_scalar = June/May,
         july_scalar = July/May,
         august_scalar = August/May,
         september_scalar = September/May,
         september_scalar = ifelse(state_park == "finley_point", 
                                   (September+August)/2/May, 
                                   september_scalar)) %>%
  dplyr::select(state_park, june_scalar:september_scalar)

flathead_parks_locs <- flathead_shapefile %>%
  filter(grepl(pattern = "State_Park", x = Name)) %>%
  mutate(state_park = gsub(pattern = "_State_Park", replacement = "", x = Name),
         state_park = tolower(state_park)) %>%
  full_join(., y = park_scalars)

flathead_park_sampling_locs <- st_distance(x = site_loc,
                                           y = flathead_parks_locs) %>%
  as_tibble() %>%
  rename("west_shore" = V1,
         "big_arm" = V2,
         "finley_point" = V3,
         "wayfarers" = V4,
         "yellow_bay" = V5) %>%
  clean_names() %>%
  cbind(., site_loc$Name) %>%
  as_tibble() %>%
  rename(sampling_site = `site_loc$Name`) %>%
  pivot_longer(cols = west_shore:yellow_bay, names_to = "state_park", values_to = "distance") %>%
  inner_join(. , y = park_scalars) %>%
  pivot_longer(cols = june_scalar:september_scalar, names_to = "month_scalar", values_to = "scalar_value") %>%
  mutate(weighted_scalar = scalar_value/(distance/1000)) %>%
  group_by(sampling_site, month_scalar) %>%
  summarize(inverse_distance_weighted_scalar = log10(sum(as.numeric(weighted_scalar))+1)) %>%
  arrange(((inverse_distance_weighted_scalar)))

## Combined space and time component

flathead_park_space_time <- inner_join(x = locs_centroids,
                                       y = flathead_park_sampling_locs) %>%
  mutate(scaled_idw_population = ((distance_weighted_population) * inverse_distance_weighted_scalar)) %>%
  arrange(((scaled_idw_population))) %>%
  separate(col = "month_scalar", into = c("month", "scalar"), remove = FALSE) %>%
  mutate(stt = ifelse(test = sampling_site %in% c("WF", "YB", "FLBS", "HO"), 
                      yes = "centralized", no = "decentralized"),
         tourist_season = ifelse(test = month %in% c("june", "july", "august"), 
                                 yes = "In Season", no = "Out of Season")) %>%
  select(-scalar)
  


write.csv(x = flathead_park_space_time,
          file = "../cleaned_disaggregated_data/temporally_scaled_inverse_distance_weighted_population_metrics.csv",
          row.names = FALSE)


## Make a map

combioned_shp <- full_join(x = flathead_shapefile, 
                           y = flathead_park_space_time, 
                           by = c("Name" = "sampling_site"))


base_map <- openmap(upperLeft = c(48.987618, -116.054273),
                    lowerRight = c(45.029625, -104.077145),
                    type = "esri") %>%
  openproj()


base_map_zoom <- openmap(upperLeft = c(48.119469, -114.398000),
                         lowerRight = c(47.663195, -113.947497),
                         type = "bing", zoom = 11) %>%
  openproj()

inset_map <- autoplot.OpenStreetMap(base_map) +
  geom_point(data = site_loc, 
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y)),
             stat = "sf_coordinates", 
             alpha = 1,  color = "grey20", 
             size = 3, shape = 21) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.background = element_rect(fill = "snow1")) +
  #annotate(geom = "text", label = "Irkutsk Oblast",
  #         x = 105, y = 54, color = "white", size = 5) +
  #annotate(geom = "text", label = "Republic of\nBuryatiya",
  #         x = 109, y = 52.35, color = "white", size = 5) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

zoom_map <- autoplot.OpenStreetMap(base_map_zoom) +
  geom_point(data = combioned_shp %>%
               filter(Name %in% c("BB", "BD", "FLBS", "YB", "SL",
                                  "SJ", "BO", "DU", "WS", "LK", "DA",
                                  "WF", "FI", "WB")) %>%
               mutate(month_scalar = factor(month_scalar, 
                                            levels = c("june_scalar", "july_scalar", "august_scalar", "september_scalar"), 
                                            labels = c("June", "July", "August", "September"))), 
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y),
                 size = log10(scaled_idw_population),
                 fill = log10(scaled_idw_population)),
             stat = "sf_coordinates", 
             alpha = 0.8,  color = "grey70", shape = 21,
             stroke = 2.5) + 
  scale_fill_viridis(option = "plasma", name = "log10(IDW Pop)") +
  scale_size_continuous(range = c(1, 12), guide = "none") +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~month_scalar) +
  theme(legend.key.height = unit(1, "in"),
        legend.key.width = unit(0.65, "in"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        panel.background = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 20))

sites_ts <- ggplot(data = flathead_park_space_time %>%
                     mutate(month_scalar = factor(month_scalar, 
                                                  levels = c("june_scalar", "july_scalar", "august_scalar", "september_scalar"), 
                                                  labels = c("June", "July", "August", "September"))), 
                   aes(month_scalar, (scaled_idw_population)))+
  geom_point(size = 5) +
  facet_wrap(~sampling_site)

## Alternative method for creating an aggregate scalar across parks

park_scalars <- park_data_orig %>%
  pivot_longer(cols = big_arm:yellow_bay, names_to = "state_park", 
               values_to = "admitted_persons") %>%
  pivot_wider(names_from = Month, values_from = admitted_persons) %>%
  mutate(may_scalar = May/May, 
         june_scalar = June/May,
         july_scalar = July/May,
         august_scalar = August/May,
         september_scalar = September/May,
         september_scalar = ifelse(state_park == "finley_point", 
                                   (September+August)/2/May, 
                                   september_scalar)) %>%
  pivot_longer(cols = may_scalar:september_scalar, 
               names_to = "month", 
               values_to = "scalar_value") %>%
  group_by(month) %>%
  summarize(average_scalar = mean(scalar_value)) %>%
  pivot_wider(names_from = "month", values_from = "average_scalar")

park_scalars %>%
  pivot_longer(cols = c(august_scalar:september_scalar), names_to = "month", values_to = "scalar") %>%
  mutate(month = factor(x = month, 
                        levels = c("june_scalar", "july_scalar", "august_scalar", "september_scalar"),
                        labels = c("June", "July", "August", "September"))) %>%
  filter(!is.na(month)) %>%
  ggplot() +
  geom_point(aes(x = month, y = scalar), size = 5) +
  geom_hline(aes(yintercept = 1.0), linetype = "dashed", size = 2) +
  ggtitle("Mean Aggregate Scalar") +
  ylab("Mean Aggregate Scalar") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 24))

locs_centroids_scaled <- locs_centroids %>%
  mutate(may_idw_pop = distance_weighted_population * park_scalars$may_scalar,
         june_idw_pop = distance_weighted_population * park_scalars$june_scalar, 
         july_idw_pop = distance_weighted_population * park_scalars$july_scalar,
         august_idw_pop = distance_weighted_population * park_scalars$august_scalar,
         september_idw_pop = distance_weighted_population * park_scalars$september_scalar) %>%
  select(-distance_weighted_population) %>%
  pivot_longer(cols = c(may_idw_pop:september_idw_pop), names_to = "month", values_to = "idw_pop") %>%
  mutate(month = gsub("_idw_pop", "", month),
         month = tolower(month)) %>%
  mutate(stt = ifelse(test = sampling_site %in% c("WF", "YB", "FLBS", "HO"), 
                      yes = "centralized", no = "decentralized"),
         tourist_season = ifelse(test = month %in% c("june", "july", "august"), 
                                 yes = "In Season", no = "Out of Season"))

write.csv(x = locs_centroids_scaled, 
          file = "../cleaned_disaggregated_data/temporally_scaled_inverse_distance_weighted_population_metrics.csv", 
          row.names = FALSE)

## Make a map

combioned_shp <- full_join(x = flathead_shapefile, 
                           y = locs_centroids_scaled, 
                           by = c("Name" = "sampling_site"))


base_map <- openmap(upperLeft = c(48.987618, -116.054273),
                    lowerRight = c(45.029625, -104.077145),
                    type = "esri") %>%
  openproj()


base_map_zoom <- openmap(upperLeft = c(48.119469, -114.398000),
                         lowerRight = c(47.663195, -113.947497),
                         type = "bing", zoom = 11) %>%
  openproj()

inset_map <- autoplot.OpenStreetMap(base_map) +
  geom_point(data = site_loc, 
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y)),
             stat = "sf_coordinates", 
             alpha = 1,  color = "grey20", 
             size = 3, shape = 21) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.background = element_rect(fill = "snow1")) +
  #annotate(geom = "text", label = "Irkutsk Oblast",
  #         x = 105, y = 54, color = "white", size = 5) +
  #annotate(geom = "text", label = "Republic of\nBuryatiya",
  #         x = 109, y = 52.35, color = "white", size = 5) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

zoom_map <- autoplot.OpenStreetMap(base_map_zoom) +
  geom_point(data = combioned_shp %>%
               filter(Name %in% c("BB", "BD", "FLBS", "YB", "SL",
                                  "SJ", "BO", "DU", "WS", "LK", "DA",
                                  "WF", "FI", "WB")) %>%
               mutate(month = factor(month,levels = c("May", "June", "July", "August", "September"), 
                                     labels = c("May", "June", "July", "August", "September"))), 
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y),
                 size = log10(idw_pop),
                 fill = log10(idw_pop)),
             stat = "sf_coordinates", 
             alpha = 0.8,  color = "grey70", shape = 21,
             stroke = 2.5) + 
  scale_fill_viridis(option = "plasma", name = "log10(IDW Pop)") +
  scale_size_continuous(range = c(1, 12), guide = "none") +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~month) +
  theme(legend.key.height = unit(1, "in"),
        legend.key.width = unit(0.65, "in"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        panel.background = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 20))


#### Boxplot to complemenet other analyses

locs_centroids_scaled <- locs_centroids_scaled %>%
  filter(month != "May") %>%
  mutate(tourist_season = ifelse(month %in% c("June", "July", "August"), "In Season", "Out of Season"),
         stt = ifelse(sampling_site %in% c("WF", "HO", "FLBS", "YB"), "Centralized", "Decentralized"))

ggplot(locs_centroids_scaled, aes(tourist_season, log10(idw_pop))) +
  geom_boxplot(alpha = 0.33, outlier.alpha = 0, width = 0.33) +
  geom_jitter() +
  xlab("") +
  ylab("log10(IDW Population)") +
  facet_wrap(~stt) +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 18))


##### Maps for Scalars

park_scalars <- park_data_orig %>%
  pivot_longer(cols = big_arm:yellow_bay, names_to = "state_park", 
               values_to = "admitted_persons") %>%
  pivot_wider(names_from = Month, values_from = admitted_persons) %>%
  mutate(june_scalar = June/May,
         july_scalar = July/May,
         august_scalar = August/May,
         september_scalar = September/May,
         september_scalar = ifelse(state_park == "finley_point", 
                                   (September+August)/2/May, 
                                   september_scalar)) %>%
  dplyr::select(state_park, june_scalar:september_scalar)
  

flathead_parks_locs <- flathead_shapefile %>%
  filter(grepl(pattern = "State_Park", x = Name)) %>%
  mutate(state_park = gsub(pattern = "_State_Park", replacement = "", x = Name),
         state_park = tolower(state_park)) %>%
  cbind(., st_coordinates(flathead_shapefile %>%
                            filter(grepl(pattern = "State_Park", x = Name)))) %>%
  st_drop_geometry() %>%
  full_join(., y = park_scalars) %>%
  pivot_longer(cols = c(june_scalar:september_scalar), names_to = "month", values_to = "scalar") %>%
  mutate(month = factor(x = month, 
                        levels = c("june_scalar", "july_scalar", "august_scalar", "september_scalar"), 
                        labels = c("June", "July", "August", "September"))) %>%
  

zoom_map <- autoplot.OpenStreetMap(base_map_zoom) +
  geom_point(data = flathead_parks_locs, 
             aes(x = X,
                 y = Y,
                 size = (scalar),
                 fill = (scalar)),
             alpha = 0.8,  color = "grey70", shape = 21,
             stroke = 2.5) + 
  scale_fill_viridis(option = "plasma", name = "IDW Pop Scalar") +
  scale_size_continuous(range = c(1, 12), guide = "none") +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~month) +
  theme(legend.key.height = unit(1, "in"),
        legend.key.width = unit(0.65, "in"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        panel.background = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 20))


# 9. Periphyton -----------------------------------------------------------

periphyton <- read.csv(file = "../raw_data/flathead_periphyton_abundance.csv")

head(periphyton)

periphyton_clean <- periphyton %>%
  rename("chlorophyta_filamentts" = "chlorophyta_filamnets")

write.csv(x = periphyton_clean, 
          file = "../cleaned_disaggregated_data/periphyton_abundance.csv", 
          row.names = FALSE)
