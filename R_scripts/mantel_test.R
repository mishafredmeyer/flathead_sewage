## This script is to perform a partial mantel test
## There may be some upstream data cleaning necessary

library(tidyverse)
library(vegan)
library(sf)

## Step 1: Get the distance matrix for the spatial data

flathead_shapefile <- sf::st_read(dsn = "../raw_data/Flathead sampling locations.kml")

sampling_loc <- c("LK", "WS", "DA", "SJ", "SL", "BO", "DU", 
                  "FI", "BB", "FLBS", "YB", "BD", "WB", "WF", "HO")

# Convert metadata into a spatial object 
site_loc <- flathead_shapefile %>%
  filter(Name %in% sampling_loc) %>%
  st_coordinates() %>%
  data.frame()  %>%
  cbind(., flathead_shapefile %>%
          filter(Name %in% sampling_loc) %>%
          st_drop_geometry() %>%
          select(Name)) %>%
  as_tibble() %>%
  select(-Z) %>%
  rename("sampling_site" = "Name",
         "Long" = "X",
         "Lat" = "Y")

site_dist <- vegdist(x = site_loc[, 1:2], method = "euclidean", diag = TRUE, upper = TRUE)

## Step 2: Get the distance matrix for the temporal, IDW pop, and PPCP data

ppcp <- read.csv("../cleaned_data/ppcp.csv")

ppcp_formatted <- ppcp %>%
  filter(!(SITE %in% c("FI2", "DU2", "HO1"))) %>%
  mutate(TIME = as.numeric(TIME),
         sampling_event = ifelse(MONTH == "MAY", "MAY_1", NA),
         sampling_event = ifelse(MONTH == "JUNE" & between(TIME, 1,6), "MAY_1", sampling_event),
         sampling_event = ifelse(MONTH == "JUNE" & between(TIME, 7, 26), "JUNE_1", sampling_event),
         sampling_event = ifelse(MONTH == "JUNE" & between(TIME, 27, 30), "JUNE_2", sampling_event),
         sampling_event = ifelse(MONTH == "JULY" & between(TIME, 1,11), "JUNE_2", sampling_event),
         sampling_event = ifelse(MONTH == "JULY" & between(TIME, 12,24), "JULY_1", sampling_event),
         sampling_event = ifelse(MONTH == "JULY" & between(TIME, 25,31), "JULY_2", sampling_event),
         sampling_event = ifelse(MONTH == "AUGUST" & between(TIME, 1,8), "JULY_2", sampling_event),
         sampling_event = ifelse(MONTH == "AUGUST" & between(TIME, 9,21), "AUGUST_1", sampling_event),
         sampling_event = ifelse(MONTH == "AUGUST" & between(TIME, 22,31), "AUGUST_2", sampling_event),
         sampling_event = ifelse(MONTH == "SEPTEMBER" & between(TIME, 1,5), "AUGUST_2", sampling_event),
         sampling_event = ifelse(MONTH == "SEPTEMBER" & between(TIME, 6,18), "SEPTEMBER_1", sampling_event),
         sampling_event = ifelse(MONTH == "SEPTEMBER" & between(TIME, 19,30), "SEPTEMBER_2", sampling_event),
         sampling_event = ifelse(MONTH == "OCTOBER", "SEPTEMBER_2", sampling_event)) %>%
  mutate(peri_sampling = ifelse(sampling_event %in% c("JUNE_1", "JUNE_2"), "Peri_1", NA),
         peri_sampling = ifelse(sampling_event %in% c("JULY_1", "JULY_2"), "Peri_2", peri_sampling),
         peri_sampling = ifelse(sampling_event %in% c("AUGUST_1", "AUGUST_2"), "Peri_3", peri_sampling),
         peri_sampling = ifelse(sampling_event %in% c("SEPTEMBER_1", "SEPTEMBER_2"), "Peri_4", peri_sampling)) %>%
  group_by(SITE, MONTH, PPCP) %>%
  summarize(mean_conc = mean(CONCENTRATION, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(SITE, MONTH) %>%
  summarize(total_conc = sum(mean_conc, na.rm = TRUE),
            mean_conc = mean(mean_conc, na.rm = TRUE)) %>%
  rename("sampling_site" = "SITE",
         "month" = "MONTH") %>%
  mutate(month = str_to_sentence(month))

idw_pop <- read.csv("../cleaned_data/temporal_inverse_distance_weighted_population_metrics.csv") %>%
  select(sampling_site, month_scalar, scaled_idw_population) %>%
  mutate(month = ifelse(grepl(pattern = "june", month_scalar), "June", NA),
         month = ifelse(grepl(pattern = "july", month_scalar), "July", month),
         month = ifelse(grepl(pattern = "august", month_scalar), "August", month),
         month = ifelse(grepl(pattern = "september", month_scalar), "September", month)) %>%
  select(-month_scalar)

month_ppcp_idw <- inner_join(idw_pop, ppcp_formatted) %>%
  mutate(log_idw_pop = log10(scaled_idw_population),
         log_tot_ppcp = log10(total_conc+1),
         log_mean_ppcp = log10(mean_conc+1))

## Step 3: Load in and format the periphyton counts and fatty acid data

periphyton_abundance <- read.csv("../cleaned_data/periphyton.csv")

periphyton_prop <- periphyton_abundance %>%
  select(-volume_ul, -volume_rep, -chlorophyta_filamnets, -date, -notes) %>%
  pivot_longer(cols = c(diatom, chlorophyta, crysophyta, cryptophyta, cyanobacteria),
               names_to = "taxon",
               values_to = "abundance") %>%
  group_by(site, month, rep) %>%
  mutate(total_site_month_rep = sum(abundance)) %>%
  ungroup() %>%
  mutate(prop_site_month_rep = abundance/total_site_month_rep) %>%
  select(-total_site_month_rep, -abundance, -rep) %>%
  group_by(site, month, taxon) %>%
  summarize(average_prop = mean(prop_site_month_rep, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(site, month) %>%
  mutate(average_prop_sum = sum(average_prop),
         average_prop_norm = average_prop/average_prop_sum) %>%
  select(-average_prop, -average_prop_sum) %>%
  pivot_wider(names_from = taxon, values_from = average_prop_norm) %>%
  rename("sampling_site" = "site")


fatty_acid_prop <- read.csv("../raw_data/Flathead_FA_cleaned_mfm_20180720.csv", header = TRUE) %>%
  gather(FATTYACID, CONC, C12.0:C28.0) %>%
  dplyr::mutate(LOC = as.character(LOC),
                MONTH = ifelse(MONTH == 7, "July", MONTH), 
                MONTH = ifelse(MONTH == 8, "August", MONTH),
                MONTH = ifelse(MONTH == 9, "September", MONTH))%>%
  dplyr::group_by(LOC, MONTH) %>%
  dplyr::mutate(TOTAL.FA = sum(CONC),
                PROP.FA = CONC/TOTAL.FA) %>%
  dplyr::select(-CONC, -TOTAL.FA) %>%
  spread(FATTYACID, PROP.FA) %>%
  as.data.frame() %>%
  rename(month = MONTH,
         sampling_site = LOC)

## Step 4: Join the three dataframes

## Step 4.1: Try with periphyton abundance 

peri_idwpop_ppcp_month <- inner_join(x = periphyton_prop,
                                     y = month_ppcp_idw)

peri_idwpop_ppcp_month_locs <- inner_join(x = peri_idwpop_ppcp_month,
                                          y = site_loc) %>%
  mutate(month = as.numeric(factor(month, 
                                   levels = c("July", "August", "September"))),
         stt = ifelse(sampling_site %in% c("WF", "FLBS", "YB"), "centralized", "decentralized"),
         stt = as.numeric(factor(stt))) %>%
  arrange(month)

partial_mantel_results <- mantel.partial(xdis = vegdist(peri_idwpop_ppcp_month_locs[,c(3,6,7)], method = "bray", diag = TRUE, upper = TRUE),
                                 ydis = vegdist(peri_idwpop_ppcp_month_locs[,14:15], method = "euclidean", diag = TRUE, upper = TRUE),
                                 zdis = vegdist(peri_idwpop_ppcp_month_locs[,c(2,13,16)], method = "mahalanobis", diag = TRUE, upper = TRUE), 
                                 method = "kendall", permutations = 999)
partial_mantel_results

mantel_results <- mantel(xdis = vegdist(peri_idwpop_ppcp_month_locs[,c(3,6,7)], method = "bray", diag = TRUE, upper = TRUE),
                                 ydis = vegdist(peri_idwpop_ppcp_month_locs[,c(2,13,16)], method = "mahalanobis", diag = TRUE, upper = TRUE), 
                                 method = "kendall", permutations = 1999)

mantel_results

mantel_results <- mantel(xdis = vegdist(peri_idwpop_ppcp_month_locs[,c(3,6,7)], method = "bray", diag = TRUE, upper = TRUE),
                         ydis = vegdist(peri_idwpop_ppcp_month_locs[,14:15], method = "euclidean", diag = TRUE, upper = TRUE), 
                         method = "kendall", permutations = 1999)

mantel_results

rda <- capscale(peri_idwpop_ppcp_month_locs[,c(3,6,7)] ~ log_mean_ppcp + log_idw_pop + month + stt, 
             distance = "bray", data = peri_idwpop_ppcp_month_locs)
sppscores(rda) <- peri_idwpop_ppcp_month_locs[,c(3,6,7)]
ordiplot(rda)
scores(rda)

## Step 4.2: Try with periphyton fatty acids

fa_idwpop_ppcp_month <- inner_join(x = fatty_acid_prop,
                                     y = month_ppcp_idw)

fa_idwpop_ppcp_month_locs <- inner_join(x = fa_idwpop_ppcp_month,
                                          y = site_loc) %>%
  mutate(month = as.numeric(factor(month, 
                                   levels = c("July", "August", "September"))),
         stt = ifelse(sampling_site %in% c("WF"), "centralized", "decentralized"),
         stt = as.numeric(factor(stt))) %>%
  arrange(month)

partial_mantel_results <- mantel.partial(xdis = vegdist(fa_idwpop_ppcp_month_locs[,c(8,16,18,20,25)], method = "bray", diag = TRUE, upper = TRUE),
                                         ydis = vegdist(fa_idwpop_ppcp_month_locs[,39:40], method = "euclidean", diag = TRUE, upper = TRUE),
                                         zdis = vegdist(fa_idwpop_ppcp_month_locs[,c(2,36)], method = "euclidean", diag = TRUE, upper = TRUE), 
                                         method = "pearson", permutations = 999)
partial_mantel_results

mantel_results <- mantel(xdis = vegdist(fa_idwpop_ppcp_month_locs[,c(8,16,18,25)], method = "bray", diag = TRUE, upper = TRUE),
                         ydis = vegdist(fa_idwpop_ppcp_month_locs[,c(2,36,38,41)], method = "euclidean", diag = TRUE, upper = TRUE), 
                         method = "kendall", permutations = 999)
mantel_results

mantel_results <- mantel(xdis = vegdist(fa_idwpop_ppcp_month_locs[,c(8,16,18,20,25)], method = "bray", diag = TRUE, upper = TRUE),
                         ydis = vegdist(fa_idwpop_ppcp_month_locs[,39:40], method = "euclidean", diag = TRUE, upper = TRUE), 
                         method = "kendall", permutations = 999)
mantel_results

## dbRDA version

rda <- dbrda(fa_idwpop_ppcp_month_locs[,c(8,16,18,20,25)] ~ log_mean_ppcp + month + stt, 
             distance = "bray", data = fa_idwpop_ppcp_month_locs)
sppscores(rda) <- fa_idwpop_ppcp_month_locs[,c(8, 16,18,20,25)]
ordiplot(rda)
scores(rda)

## STI models

stimodels(fa_idwpop_ppcp_month_locs[fa_idwpop_ppcp_month_locs$sampling_site %in% c("BD", "BO", "FI", "FLBS", "SJ", "SL", "WB", "YB"),c(8,16,18,20,25)], 
          S = unique(fa_idwpop_ppcp_month_locs[fa_idwpop_ppcp_month_locs$sampling_site %in% c("BD", "BO", "FI", "FLBS", "SJ", "SL", "WB", "YB"),c(39:40)]), 
          Ti = 3, nperm = 199, model = "5")

stimodels(peri_idwpop_ppcp_month_locs[peri_idwpop_ppcp_month_locs$sampling_site %in% c("BO", "FI", "FLBS", "SJ", "SL", "WB", "YB"),c(3,6,7)], 
          S = unique(peri_idwpop_ppcp_month_locs[peri_idwpop_ppcp_month_locs$sampling_site %in% c("BO", "FI", "FLBS", "SJ", "SL", "WB", "YB"),c(14:15)]), 
          Ti = 3, nperm = 199, model = "5")

test <- dbmem(
  xyORdist = unique(peri_idwpop_ppcp_month_locs[peri_idwpop_ppcp_month_locs$sampling_site %in% c("BO", "FI", "FLBS", "SJ", "SL", "WB", "YB"),c(14:15)]),
  thresh = NULL,
  #MEM.autocor = c("positive", "non-null", "all", "negative"),
  store.listw = TRUE,
  silent = TRUE)
