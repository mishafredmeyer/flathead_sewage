## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script is a piece of the 
## overarching analysis and looks mostly at relating 
## sewage indicators with wastewater treatment infrastructure
## and a sampling point's temporal position to the tourist 
## season. Questions about this script should be directed to 
## Michael F. Meyer (michael.f.meyer@wsu.edu). 

ggplot(afdm_edit, #%>%
         #filter(MONTH != "JUNE",
        #        SITE != "YB"), 
       aes(tourist_season, mean_afdm)) +
  geom_boxplot(alpha = 0.33, outlier.alpha = 0) +
  geom_jitter() +
  facet_wrap(~stt) +
  theme_bw()

afdm_lm <- lm(mean_afdm ~ stt*tourist_season, 
              data = afdm_edit)

Anova(afdm_lm, type = "II")

fatty_acids <- fatty_acids %>%
  rename("SITE" = "LOC") %>%
  mutate(MONTH = ifelse(MONTH == 9, "SEPTEMBER", MONTH),
         MONTH = ifelse(MONTH == 8, "AUGUST", MONTH),
         MONTH = ifelse(MONTH == 7, "JULY", MONTH))

SAFA <- c("C12.0", "C14.0", "C15.0", "C16.0", "C17.0", "C18.0", "C20.0", "C22.0", "iso.C15.0", "C24.0", "C26.0", "C28.0")
MUFA <- c( "C15.1", "C14.1n5", "C15.1w7", "C16.1w5", "C16.1w6", "C16.1w7", "C16.1w7c",  "C16.1w8", "C16.1w9", "C17.1n7", "C18.1w7", "C18.1w7c", 
           "C18.1w9", "C18.1w9c", "C20.1w7", "C20.1w9", "C22.1w7", "C22.1w9", "C22.1w9c")
LUFA <- c("C16.2", "C16.2w4",  "C16.2w6",  "C16.2w7",  "C16.3w3",  "C16.3w4",  "C16.3w6",  "C18.2w6",  "C18.2w6t", "C18.2w6c",
          "C18.3w3",  "C18.3w6", "C20.2w6",  "C20.3w3",  "C20.3w6",  "C22.2w6",  "C22.3w3")
HUFA <- c("C16.4w1", "C16.4w3", "C18.4w3", "C18.4w4", "C18.5w3", "C20.4w2", "C20.4w3", "C20.4w6", "C20.5w3", "C22.4w3", 
          "C22.4w6", "C22.5w3", "C22.5w6", "C22.6w3")
SCUFA_LUFA <- c("C16.2", "C16.2w4",  "C16.2w6",  "C16.2w7",  "C16.3w3",  "C16.3w4",  "C16.3w6",  "C18.2w6", "C18.2w6c", 
                "C18.2w6t", "C18.3w3",  "C18.3w6")
LCUFA_LUFA <- c("C20.2w6",  "C20.3w3",  "C20.3w6",  "C22.2w6",  "C22.3w3")
SCUFA_HUFA <- c("C16.4w1", "C16.4w3", "C18.4w3", "C18.4w4", "C18.5w3")
LCUFA_HUFA <- c( "C20.4w2", "C20.4w3", "C20.4w6", "C20.5w3", "C22.4w3", 
                 "C22.4w6", "C22.5w3", "C22.5w6", "C22.6w3")
SCUFA <- c("C16.2w4",  "C16.2w6",  "C16.2w7",  "C16.3w3",  "C16.3w4",  "C16.3w6",  "C18.2w6",  "C18.2w6t", "C18.2w6t",
           "C18.3w3",  "C18.3w6", "C16.4w1", "C16.4w3", "C18.4w3", "C18.4w4", "C18.5w3")
LCUFA <- c( "C20.4w2", "C20.4w3", "C20.4w6", "C20.5w3", "C22.4w3", 
            "C22.4w6", "C22.5w3", "C22.5w6", "C22.6w3", "C20.2w6",  "C20.3w3",  "C20.3w6",  "C22.2w6",  "C22.3w3")

C18PUFA <- c("C18.2w6",  "C18.2w6t", "C18.3w3",  "C18.3w6", "C18.4w3", "C18.4w4", "C18.5w3")
C20PUFA <- c("C20.4w3", "C20.4w6", "C20.5w3", "C20.2w6",  "C20.3w3",  "C20.3w6")

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


nutrients_edit <- nutrients %>%
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
            mean_TP = mean(TP)) %>%
  pivot_longer(cols = c(mean_nh3:mean_TP), names_to = "nutrient_type", values_to = "concentration") %>%
  mutate(stt = ifelse(SITE %in% c("WF", "YB", "WB"), "Centralized", "Decentralized"),
         tourist_season = ifelse(MONTH %in% c("JUNE", "JULY", "AUGUST"), "In Season", "Out of Season"))


ggplot(nutrients_edit, 
       aes(tourist_season, log10(concentration), fill = stt)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(aes(color = stt)) +
  scale_y_log10() +
  facet_wrap(~nutrient_type) +
  theme_bw()

phos_lm <- lm(log10(concentration) ~ stt*tourist_season, 
              data = nutrients_edit %>%
                filter(nutrient_type == "mean_NO3"))

Anova(phos_lm, type = "II")

