## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script is a piece of the 
## overarching analysis and looks mostly at relating 
## sewage indicators with wastewater treatment infrastructure
## and a sampling point's temporal position to the tourist 
## season. Questions about this script should be directed to 
## Michael F. Meyer (michael.f.meyer@wsu.edu). 


# 1. Load Packages and Data -----------------------------------------------

library(tidyverse)
library(viridis)
library(car)
library(ggpubr)

afdm <- read.csv(file = "../cleaned_data/afdm.csv",
                 header = TRUE)
  
fatty_acids <- read.csv(file = "../cleaned_data/fatty_acids.csv",
                        header = TRUE)  

nutrients <- read.csv(file = "../cleaned_data/nutrients.csv",
                      header = TRUE)

ppcp <- read.csv(file = "../cleaned_data/ppcp.csv",
                 header = TRUE)

tsidw_pop <- read.csv(file = "../cleaned_data/averaged_temporally_scaled_inverse_distance_weighted_population_metrics.csv",
                      header = TRUE)


# 2. Nutrient Analysis ----------------------------------------------------

nutrient_labels <- c("Ammonia as Nitrogen", "Nitrate-Nitrite", "Total Nitrogen", "SRP",  "Total Phosphorus")
names(nutrient_labels) <- c("nh3_n", "no3_no2", "total_n", "srp", "total_p")

nutrient_plot <- nutrients %>%
  pivot_longer(cols = c(nh3_n:total_p), names_to = "nutrient_type", values_to = "concentration") %>%
  filter(!site %in% c("DU", "HO")) %>%
  ggplot(aes(tourist_season, concentration, fill = stt, color = stt)) +
  geom_boxplot(alpha = 0.7, width = 0.33, outlier.alpha = 0) +
  geom_jitter() +
  scale_y_log10() +
  ylab(expression(paste("Concentration (\u03BCg/L- N or P)"))) +
  scale_fill_manual(values = plasma(30)[c(5, 19)], name = "STT") +
  scale_color_manual(values = plasma(30)[c(5, 19)], name = "STT") +
  facet_wrap(~nutrient_type, scales = "free", labeller = labeller(nutrient_type = nutrient_labels)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, units = "in"),
        legend.key.width = unit(0.75, units = "in"),
        strip.text = element_text(size = 20))
  
ggsave(filename = "nutrient_boxplots.png", plot = nutrient_plot, 
       device = "png", path = "../figures_tables", 
       width = 15, height = 7, units = "in")


# Ammonium Model
nh3_lm <- lm(nh3_n ~ stt*tourist_season, 
              data = nutrients %>% filter(!site %in% c("DU", "HO")))

Anova(nh3_lm, type = "II")

# Nitrate/Nitrite Model
no3_no2_lm <- lm(no3_no2 ~ stt*tourist_season, 
             data = nutrients %>% filter(!site %in% c("DU", "HO")))

Anova(no3_no2_lm, type = "II")

# SRP
srp_lm <- lm(srp ~ stt*tourist_season, 
             data = nutrients %>% filter(!site %in% c("DU", "HO")))

Anova(srp_lm, type = "II")

# Total Nitrogen
tn_lm <- lm(total_n ~ stt*tourist_season, 
             data = nutrients %>% filter(!site %in% c("DU", "HO")))

Anova(tn_lm, type = "II")

# Total Phosphorus
tp_lm <- lm(total_p ~ stt*tourist_season, 
            data = nutrients %>% filter(!site %in% c("DU", "HO")))

Anova(tp_lm, type = "II")



# 2. Ash Free Dry Mass ----------------------------------------------------

stt_labels <- c("Centralized", "Decentralized")
names(stt_labels) <- c("centralized", "decentralized")

afdm_plot <- ggplot(afdm, aes(tourist_season, afdm)) +
  geom_boxplot(alpha = 0.33, outlier.alpha = 0, width = 0.2) +
  geom_jitter(width = 0.25) +
  facet_wrap(~stt, labeller = labeller(stt = stt_labels)) +
  ylab("Ash Free Dry Mass (mg)") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, units = "in"),
        legend.key.width = unit(0.75, units = "in"),
        strip.text = element_text(size = 20))

ggsave(filename = "afdm_boxplots.png", plot = afdm_plot, 
       device = "png", path = "../figures_tables", 
       width = 8, height = 6, units = "in")

afdm_lm <- lm(afdm ~ stt*tourist_season, 
              data = afdm)

Anova(afdm_lm, type = "II")


# 3. Branched and Odd Chain Fatty Acids -----------------------------------

branched_odd_chain_fatty_acids_plot <- fatty_acids %>%
  select(-c19_0) %>%
  pivot_longer(cols = c(c12_0:c28_0), names_to = "fatty_acid", values_to = "concentration") %>%
  group_by(tourist_season, stt, site, month) %>%
  mutate(total_fatty_acid = sum(concentration, na.rm = TRUE),
                prop_fatty_acid = concentration/total_fatty_acid) %>%
  select(-concentration, -total_fatty_acid) %>%
  pivot_wider(names_from = "fatty_acid", values_from = "prop_fatty_acid") %>%
  select(site, month, stt, tourist_season, contains("15"), contains("17")) %>%
  ungroup() %>%
  group_by(stt, tourist_season, site, month) %>%
  summarize(total_branched_odd = rowSums(across(.cols = iso_c15_0:c17_0))) %>%
  ggplot(aes(x = tourist_season, y = total_branched_odd)) +
  geom_boxplot(alpha = 0.33, width = 0.2, outlier.alpha = 0) +
  xlab("")+
  ylab("Odd-chain + Branched Fatty Acids") +
  geom_jitter(size = 2, width = 0.25) +
  facet_wrap(~stt) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"))


ggsave(filename = "branched_odd_chain_fatty_acid_boxplots.png", plot = branched_odd_chain_fatty_acids_plot, 
       device = "png", path = "../figures_tables", 
       width = 8, height = 6, units = "in")


branched_odd_chain_fatty_acids <- fatty_acids %>%
  select(-c19_0) %>%
  pivot_longer(cols = c(c12_0:c28_0), names_to = "fatty_acid", values_to = "concentration") %>%
  group_by(tourist_season, stt, site, month) %>%
  mutate(total_fatty_acid = sum(concentration, na.rm = TRUE),
         prop_fatty_acid = concentration/total_fatty_acid) %>%
  select(-concentration, -total_fatty_acid) %>%
  pivot_wider(names_from = "fatty_acid", values_from = "prop_fatty_acid") %>%
  select(site, month, stt, tourist_season, contains("15"), contains("17")) %>%
  ungroup() %>%
  group_by(stt, tourist_season, site, month) %>%
  summarize(total_branched_odd = rowSums(across(.cols = iso_c15_0:c17_0)))


branched_odd_chain_fatty_acids_lm <- lm(total_branched_odd ~ stt * tourist_season,
                                        data = branched_odd_chain_fatty_acids)

Anova(branched_odd_chain_fatty_acids_lm, type = "II")


Anova(lm(total_branched_odd ~  tourist_season,
         data = branched_odd_chain_fatty_acids %>%
           filter(stt == "Decentralized")))

Anova(lm(total_branched_odd ~  tourist_season,
         data = branched_odd_chain_fatty_acids %>%
           filter(stt == "Centralized")))


# 4. PPCPs ----------------------------------------------------------------

ppcp_reduced <- ppcp %>%
  #filter(!site %in% c("HO", "DU")) %>%
  #filter(month != "may") %>%
  # group_by(tourist_season, stt, site, month, sampling_event) %>%
  # summarize(mean_conc = sum(concentration, na.rm = TRUE)) %>%
  group_by(tourist_season, stt, site, peri_sampling, sampling_event, ppcp) %>%
  summarize(mean_conc = mean(concentration, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(tourist_season, stt, site, peri_sampling, sampling_event) %>%
  summarize(total_concentration = sum(mean_conc, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(tourist_season, stt, site, peri_sampling) %>%
  summarize(mean_conc = mean(total_concentration, na.rm = TRUE)) 

  
ppcp_plot <- ggplot(ppcp_reduced, 
       aes(tourist_season, (log10(mean_conc)))) +
  geom_boxplot(alpha = 0.33, outlier.alpha = 0, width = 0.2) +
  geom_jitter() +
  facet_wrap(~stt, labeller = labeller(stt = stt_labels)) +
  ylab("log10([Total PPCP])") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, units = "in"),
        legend.key.width = unit(0.75, units = "in"),
        strip.text = element_text(size = 20))

ggsave(filename = "ppcp_boxplots.png", plot = ppcp_plot, 
       device = "png", path = "../figures_tables", 
       width = 8, height = 6, units = "in")

## ANOVA with stt and tourist season

Anova(lm(log10(mean_conc) ~ stt * tourist_season, 
         data = ppcp_reduced),
      type = "II")

## ANOVA just with Decentralized treatment

Anova(lm(log10(mean_conc) ~ tourist_season, 
         data = ppcp_reduced %>%
           filter(stt == "decentralized")),
      type = "II")

## ANOVA just with Centralized treatment

Anova(lm(log10(mean_conc) ~ tourist_season, 
         data = ppcp_reduced %>%
           filter(stt == "centralized")),
      type = "II")


# 5. TSIDW ----------------------------------------------------------------

tsidw_plot <- ggplot(tsidw_pop, 
                    aes(tourist_season, (log10(atsidw_pop)))) +
  geom_boxplot(alpha = 0.33, outlier.alpha = 0, width = 0.2) +
  geom_jitter() +
  facet_wrap(~stt, labeller = labeller(stt = stt_labels)) +
  ylab("log10(TSIDW Population)") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, units = "in"),
        legend.key.width = unit(0.75, units = "in"),
        strip.text = element_text(size = 20))

ggsave(filename = "tsidw_population_boxplots.png", plot = tsidw_plot, 
       device = "png", path = "../figures_tables", 
       width = 8, height = 6, units = "in")

## ANOVA with stt and tourist season

Anova(lm(log10(atsidw_pop) ~ stt * tourist_season, 
         data = locs_centroids_scaled),
      type = "II")

## ANOVA just with Decentralized treatment

Anova(lm(log10(atsidw_pop) ~ tourist_season, 
         data = tsidw_pop %>%
           filter(stt == "decentralized")),
      type = "II")

## ANOVA just with Centralized treatment

Anova(lm(log10(atsidw_pop) ~ tourist_season, 
         data = tsidw_pop %>%
           filter(stt == "centralized")),
      type = "II")


# 6. Combine Figures ------------------------------------------------------

arranged_plots <- ggarrange(plotlist = list(tsidw_plot, ppcp_plot,
                                         branched_odd_chain_fatty_acids_plot, afdm_plot), 
                            ncol = 2, nrow = 2, labels = "AUTO",
                            font.label = list(size = 24))

ggsave(filename = "combined_boxplots.png", plot = arranged_plots, 
       device = "png", path = "../figures_tables", 
       width = 16, height = 12, units = "in")
