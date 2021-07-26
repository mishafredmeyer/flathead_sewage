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


stt_labels <- c("Centralized", "Decentralized")
names(stt_labels) <- c("centralized", "decentralized")

carbon_nitrogen <- ggplot(stoich, aes(tourist_season, carbon/nitrogen)) +
  geom_boxplot(alpha = 0.33, outlier.alpha = 0) +
  geom_jitter() +
  geom_hline(aes(yintercept = 119/17), size = 2, linetype = "dotted") + 
  facet_wrap(~stt, labeller = labeller(stt = stt_labels)) +
  ylab("Carbon:Nitrogen") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, units = "in"),
        legend.key.width = unit(0.75, units = "in"),
        strip.text = element_text(size = 20))

ggsave(filename = "carbon_nitrogen_boxplots.png", plot = carbon_nitrogen, 
       device = "png", path = "../figures_tables", 
       width = 8, height = 6, units = "in")

carbon_phosphorus <- ggplot(stoich, aes(tourist_season, carbon/phosphorus)) +
  geom_boxplot(alpha = 0.33, outlier.alpha = 0) +
  geom_jitter() +
  geom_hline(aes(yintercept = 119), size = 2, linetype = "dotted") + 
  facet_wrap(~stt, labeller = labeller(stt = stt_labels)) +
  ylab("Carbon:Phosphorus") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, units = "in"),
        legend.key.width = unit(0.75, units = "in"),
        strip.text = element_text(size = 20))

ggsave(filename = "carbon_phosphorus_boxplots.png", plot = carbon_phosphorus, 
       device = "png", path = "../figures_tables", 
       width = 8, height = 6, units = "in")

nitrogen_phosphorus <- ggplot(stoich, aes(tourist_season, nitrogen/phosphorus)) +
  geom_boxplot(alpha = 0.33, outlier.alpha = 0) +
  geom_jitter() +
  geom_hline(aes(yintercept = 17), size = 2, linetype = "dotted") + 
  facet_wrap(~stt, labeller = labeller(stt = stt_labels)) +
  ylab("Nitrogen:Phosphorus") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, units = "in"),
        legend.key.width = unit(0.75, units = "in"),
        strip.text = element_text(size = 20))

ggsave(filename = "nitrogen_phosphorus_boxplots.png", plot = nitrogen_phosphorus, 
       device = "png", path = "../figures_tables", 
       width = 8, height = 6, units = "in")

combined_plots <- ggarrange(plotlist = list(carbon_nitrogen, carbon_phosphorus, nitrogen_phosphorus), 
                            ncol = 1, labels = "AUTO",
                            font.label = list(size = 24))

ggsave(filename = "combined_stoich_boxplots.png", plot = combined_plots, 
       device = "png", path = "../figures_tables", 
       width = 8, height = 16, units = "in")


Anova(lm(carbon/nitrogen ~ tourist_season*stt, data = stoich), type = "II")

Anova(lm(carbon/phosphorus ~ tourist_season*stt, data = stoich), type = "II")

Anova(lm(nitrogen/phosphorus ~ tourist_season*stt, data = stoich), type = "II")
