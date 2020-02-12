library(tidyverse)
library(vegan)
library(lubridate)
library(stringi)

afdm <- read.csv("../cleaned_data/AFDM.csv", header = TRUE)
stoich <- read.csv("../cleaned_data/stoichiometry.csv", header = TRUE)
ppcp <- read.csv("../cleaned_data/ppcp.csv", header = TRUE)
utilities <- read.csv("../cleaned_data/utilities.csv", header = TRUE)
nutrients <- read.csv("../cleaned_data/nutrients.csv", header = TRUE)
fatty_acids <- read.csv("../cleaned_data/fatty_acids.csv")

afdm <- afdm %>%
  rename("MONTH" = "month") %>%
  mutate(SITE = trimws(site)) %>%
  mutate(SITE = ifelse(site == "Beardance", "BD", SITE),
         SITE = ifelse(site == "Blue Bay", "BB", SITE),
         SITE = ifelse(site == "Boettcher", "BO", SITE),
         SITE = ifelse(site == "Dayton", "DA", SITE),
         SITE = ifelse(site == "Finley Point", "FI", SITE),
         SITE = ifelse(site == "FLBS", "FLBS", SITE),
         SITE = ifelse(site == "Lakeside", "LK", SITE),
         SITE = ifelse(site == "Sacajawea", "SJ", SITE),
         SITE = ifelse(site == "Salish Park", "SA", SITE),
         SITE = ifelse(site == "Wayfarers", "WF", SITE),
         SITE = ifelse(site == "West Shore", "WS", SITE),
         SITE = ifelse(site == "Woods Bay", "WB", SITE),
         SITE = ifelse(site == "Yellow Bay", "YB", SITE)) %>%
  mutate(AFDM = dry_weight - after_ignition) %>%
  group_by(SITE, MONTH) %>%
  summarise(mean_afdm = mean(AFDM)) %>%
  mutate(z_scored_afdm = scale(mean_afdm))

fatty_acids <- fatty_acids %>%
  rename("SITE" = "LOC") %>%
  mutate(MONTH = ifelse(MONTH == 9, "SEPTEMBER", MONTH),
         MONTH = ifelse(MONTH == 8, "AUGUST", MONTH),
         MONTH = ifelse(MONTH == 7, "JULY", MONTH))

fatty_acids_reduced <- fatty_acids %>%
  select(-C19.0) %>%
  gather(FATTYACID, CONC, C12.0:C28.0) %>%
  dplyr::mutate(TYPE = ifelse(FATTYACID %in% SAFA, "SAFA", "OTHER"),
                TYPE = ifelse(FATTYACID %in% MUFA, "MUFA", TYPE),
                TYPE = ifelse(FATTYACID %in% SCUFA_LUFA, "SCUFA_LUFA", TYPE),
                TYPE = ifelse(FATTYACID %in% SCUFA_HUFA, "SCUFA_HUFA", TYPE),
                TYPE = ifelse(FATTYACID %in% LCUFA_HUFA, "LCUFA_HUFA", TYPE),
                TYPE = ifelse(FATTYACID %in% LCUFA_LUFA, "LCUFA_LUFA", TYPE)) %>%
  dplyr::group_by(SITE, MONTH, TYPE) %>%
  dplyr::summarize(TOTAL_FA_TYPE = sum(CONC)) %>%
  ungroup() %>%
  group_by(SITE, MONTH) %>%
  mutate(TOTAL_FA_OVERALL = sum(TOTAL_FA_TYPE)) %>%
  ungroup() %>%
  mutate(PROP_FA = TOTAL_FA_TYPE/TOTAL_FA_OVERALL) %>%
  dplyr::select(-TOTAL_FA_OVERALL, -TOTAL_FA_TYPE) %>%
  spread(TYPE, PROP_FA)

utilities <- utilities %>%
  mutate(SITE = ifelse(Site == "Big_Fork", "WF", NA),
         SITE = ifelse(Site == "Lakeside", "LK", SITE),
         SITE = ifelse(Site == "Woods_Bay", "WB", SITE),
         MONTH = ifelse(Time == "Aug", "AUGUST", Time),
         MONTH = ifelse(Time == "Jun", "JUNE", MONTH),
         MONTH = ifelse(Time == "Jul", "JULY", MONTH),
         MONTH = ifelse(Time == "Sep", "SEPTEMBER", MONTH),
         MONTH = ifelse(Time == "May", "MAY", MONTH)) %>%
  group_by(SITE) %>%
  mutate(z_Scored_usage = scale(Cum_usage))

ppcp <- ppcp %>%
  group_by(SITE, MONTH) %>%
  summarize(total_ppcp = sum(CONCENTRATION, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(SITE) %>%
  mutate(z_scored_conc = scale(total_ppcp))

nutrients <- nutrients %>%
  mutate(SITE = trimws(Sample_ID)) %>%
  separate(SITE, c("SITE", "REP")) %>%
  mutate(SITE = ifelse(SITE == "Beardance", "BD", SITE),
         SITE = ifelse(SITE == "Blue Bay", "BB", SITE),
         SITE = ifelse(SITE == "Boettcher", "BO", SITE),
         SITE = ifelse(SITE == "Dayton", "DA", SITE),
         SITE = ifelse(SITE == "Finley Point", "FI", SITE),
         SITE = ifelse(SITE == "FLBS", "FLBS", SITE),
         SITE = ifelse(SITE == "Lakeside", "LK", SITE),
         SITE = ifelse(SITE == "Sacajawea", "SJ", SITE),
         SITE = ifelse(SITE == "Salish Park", "SA", SITE),
         SITE = ifelse(SITE == "Wayfarers", "WF", SITE),
         SITE = ifelse(SITE == "West Shore", "WS", SITE),
         SITE = ifelse(SITE == "Woods", "WB", SITE),
         SITE = ifelse(SITE == "Yellow Bay", "YB", SITE),
         MONTH = month(mdy(Collection.Date), label = TRUE, abbr = FALSE),
         MONTH = toupper(MONTH)) %>%
  group_by(SITE, MONTH) %>%
  summarize(mean_nh3 = mean(NH3),
            mean_NO3 = mean(NO3),
            mean_TN = mean(TN),
            mean_SRP = mean(SRP),
            mean_TP = mean(TP))

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
  mutate(MONTH = toupper(month(mdy(date_collected), label = TRUE, abbr = FALSE)),
         SITE = unlist(str_extract_all(SITE_full,  "(?<=\\().+?(?=\\))"))) %>%
  mutate_at(.vars = c("total_dry_mass_carbon_mg", "carbon_mg",
                      "total_dry_mass_nitrogen_mg", "nitrogen_mg",
                      "total_dry_mass_phosphorus_mg", "phosphorus_ug"), 
            as.numeric)

rda_data <- inner_join(fatty_acids, utilities) %>%
  inner_join(., ppcp) %>%
  inner_join(., afdm) %>%
  inner_join(., nutrients) %>%
  inner_join(., fatty_acids_reduced) %>%
  inner_join(., stoich_cleaned) %>%
  mutate(carbon_mol_per_mg_tissue = log10(((carbon_mg/1000)/12.01)/total_dry_mass_carbon_mg),
         nitrogen_mol_per_mg_tissue = log10(((nitrogen_mg/1000)/14)/total_dry_mass_nitrogen_mg),
         phosphorus_mol_per_mg_tissue = log10(((phosphorus_ug/1000000)/30.97)/total_dry_mass_phosphorus_mg)) %>%
  select(SITE, MONTH, total_ppcp, z_scored_conc, mean_afdm, z_scored_afdm, Cum_usage,
                mean_nh3, mean_NO3, mean_TN, mean_SRP, mean_TP, LCUFA_HUFA:SCUFA_LUFA,
         carbon_mol_per_mg_tissue, nitrogen_mol_per_mg_tissue, phosphorus_mol_per_mg_tissue, z_Scored_usage) %>%
  drop_na()

rda_results <- dbrda(rda_data[ , c(12, 13:15)] ~  z_Scored_usage + mean_NO3 + z_scored_conc + mean_SRP, 
                     data = rda_data, distance = "euclidean")
rda_results
plot(rda_results)
anova(rda_results, permutations = 5000)
anova(rda_results, by="axis")
anova(rda_results, by="terms", permu=5000)

combined_data <- #inner_join(fatty_acids, utilities) %>%
  inner_join(fatty_acids, ppcp) %>%
  inner_join(., afdm) %>%
  inner_join(., nutrients) %>%
  inner_join(., fatty_acids_reduced) %>%
  inner_join(., stoich_cleaned) %>%
  mutate(carbon_mol_per_mg_tissue = log10(((carbon_mg/1000)/12.01)/total_dry_mass_carbon_mg),
         nitrogen_mol_per_mg_tissue = log10(((nitrogen_mg/1000)/14)/total_dry_mass_nitrogen_mg),
         phosphorus_mol_per_mg_tissue = log10(((phosphorus_ug/1000000)/30.97)/total_dry_mass_phosphorus_mg)) %>%
  select(SITE, MONTH, total_ppcp, z_scored_conc, mean_afdm, z_scored_afdm, 
          mean_NO3, mean_SRP, LCUFA_HUFA:SCUFA_LUFA,
         carbon_mol_per_mg_tissue, nitrogen_mol_per_mg_tissue, phosphorus_mol_per_mg_tissue) %>%
  drop_na()

combined_data <- inner_join(fatty_acids, utilities) %>%
  inner_join(nutrients, ppcp) %>%
  inner_join(., afdm) %>%
  inner_join(., nutrients) %>%
  inner_join(., fatty_acids_reduced) %>%
  inner_join(., stoich_cleaned) %>%
  mutate(carbon_mol_per_mg_tissue = log10(((carbon_mg/1000)/12.01)/total_dry_mass_carbon_mg),
        nitrogen_mol_per_mg_tissue = log10(((nitrogen_mg/1000)/14)/total_dry_mass_nitrogen_mg),
        phosphorus_mol_per_mg_tissue = log10(((phosphorus_ug/1000000)/30.97)/total_dry_mass_phosphorus_mg)) %>%
  select(SITE, MONTH, total_ppcp, z_scored_conc, mean_afdm, z_scored_afdm,
      LCUFA_HUFA:SCUFA_LUFA, mean_SRP, mean_NO3,
        carbon_mol_per_mg_tissue, nitrogen_mol_per_mg_tissue, phosphorus_mol_per_mg_tissue) %>%
  drop_na()
ggplot(combined_data, aes(x = log10(total_ppcp),
                          y = SCUFA_LUFA,
                          color = MONTH)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(rda_data %>% filter(z_scored_conc < 3), aes(x = log10(Cum_usage),
                          y = log10(total_ppcp))) +
  geom_point(aes(color = SITE, shape = MONTH)) +
  geom_smooth(method = "lm")

kendall_data <- combined_data %>% 
  filter(MONTH == "AUGUST") %>%
  mutate(log10_ppcp = log10(total_ppcp),
         FA_response = (nitrogen_mol_per_mg_tissue/phosphorus_mol_per_mg_tissue))

kendall_results <- Kendall(kendall_data$FA_response, 
                           kendall_data$log10_ppcp)
summary(kendall_results)

summary(glm(FA_response ~ log10_ppcp, data = kendall_data))
