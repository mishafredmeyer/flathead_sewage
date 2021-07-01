## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script is a piece of the 
## overarching analysis and looks mostly at relating 
## periphyton fatty acid profiles with wastewater treatment infrastructure
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

fatty_acids_orig <- read.csv(file = "../cleaned_data/fatty_acids.csv",
                            header = TRUE)

## Define the fatty acid categories in this analysis

safa <- c("c12_0", "c14_0", "c15_0", "c16_0", "c17_0", "c18_0", "c20_0", "c22_0", 
          "iso_c15_0", "c24_0", "c26_0", "c28_0")

mufa <- c( "c15_1", "c14_1n5", "c15_1w7", "c16_1w5", "c16_1w6", "c16_1w7", 
           "c16_1w7c",  "c16_1w8", "c16_1w9", "c17_1n7", "c18_1w7", "c18_1w7c", 
           "c18_1w9", "c18_1w9c", "c20_1w7", "c20_1w9", "c22_1w7", "c22_1w9", "c22_1w9c")

scufa <- c("c16_2", "c16_2w4",  "c16_2w6",  "c16_2w7",  "c16_3w3",  "c16_3w4",  "c16_3w6",  
           "c18_2w6",  "c18_2w6t", "c18_2w6c",
           "c18_3w3",  "c18_3w6", "c16_4w1", "c16_4w3", "c18_4w3", "c18_4w4", "c18_5w3")

lcufa <- c( "c20_4w2", "c20_4w3", "c20_4w6", "c20_5w3", "c22_4w3", 
            "c22_4w6", "c22_5w3", "c22_5w6", "c22_6w3", "c20_2w6",  "c20_3w3",  "c20_3w6", 
            "c22_2w6",  "c22_3w3")

# 2. Convert data into proportions ----------------------------------------

fatty_acid_prop <- fatty_acids_orig %>%
  select(-c19_0) %>%
  pivot_longer(cols = c(c12_0:c28_0), 
               names_to = "fatty_acid", 
               values_to = "concentration") %>%
  group_by(site, month) %>%
  mutate(total_fatty_acids = sum(concentration, na.rm = TRUE),
         proportion = concentration/total_fatty_acids) %>%
  dplyr::select(-concentration, -total_fatty_acids, -type) %>%
  pivot_wider(names_from = "fatty_acid", values_from = "proportion")

mean <- as.vector(sapply(fatty_acid_prop[,5:33], mean))
var <- as.vector(sapply(fatty_acid_prop[,5:33], var))
mean.var <- data.frame(cbind(mean[1:29], var[1:29]))
colnames(mean.var)[colnames(mean.var) == "X1"] <- "Mean"
colnames(mean.var)[colnames(mean.var) == "X2"] <- "Variance"
mean.var <- dplyr::mutate(mean.var, Var.Mean.RATIO = Variance/Mean)
row.names(mean.var) <- colnames(fatty_acid_prop[,5:33])
mean.var[order(-mean.var$Var.Mean.RATIO),]

summary(fatty_acid_prop)

# 3. Univariate Analysis --------------------------------------------------

stt_labels <- c("Centralized", "Decentralized")
names(stt_labels) <- c("centralized", "decentralized")

fatty_acid_type_boxplot <- fatty_acid_prop %>%
  pivot_longer(cols = c(c12_0:c28_0), 
               names_to = "fatty_acid", 
               values_to = "proportion") %>%
  mutate(type = ifelse(fatty_acid %in% safa, "SAFA", "other"),
         type = ifelse(fatty_acid %in% mufa, "MUFA", type),
         type = ifelse(fatty_acid %in% c(scufa, lcufa), "PUFA", type),
         type = factor(type, levels = c("SAFA", "MUFA", "PUFA"))) %>%
  group_by(tourist_season, stt, month, site, type) %>%
  summarize(proportion_sum = sum(proportion)) %>%
  ggplot(aes(type, proportion_sum, fill = tourist_season)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
  geom_jitter(aes(color = tourist_season), shape = 21) +
  scale_color_manual(values = viridis(20)[c(5, 14)], 
                     name = "Tourist Season") +
  scale_fill_manual(values = viridis(20)[c(5, 14)], 
                    name = "Tourist Season") +
  ylab("Proportion") +
  xlab("") +
  facet_wrap(~stt) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"))

ggsave(filename = "fatty_acid_type_boxplot.png", plot = fatty_acid_type_boxplot, 
       device = "png", path = "../figures_tables", 
       width = 12, height = 7, units = "in")

fatty_acid_type_data <- fatty_acid_prop %>%
  pivot_longer(cols = c(c12_0:c28_0), 
               names_to = "fatty_acid", 
               values_to = "proportion") %>%
  mutate(type = ifelse(fatty_acid %in% safa, "SAFA", "other"),
         type = ifelse(fatty_acid %in% mufa, "MUFA", type),
         type = ifelse(fatty_acid %in% c(scufa, lcufa), "PUFA", type),
         type = factor(type, levels = c("SAFA", "MUFA", "PUFA"))) %>%
  group_by(tourist_season, stt, month, site, type) %>%
  summarize(proportion_sum = sum(proportion)) %>%
  pivot_wider(names_from = "type", values_from = "proportion_sum")

Anova(lm(PUFA ~ tourist_season*stt, 
         data = fatty_acid_type_data), 
      type = "II")

efa_boxplot <- fatty_acid_prop %>%
  select(site:stt, c18_2w6c, c18_3w3, c18_4w3, 
         c20_5w3, c22_6w3, c20_4w6) %>%
  pivot_longer(cols = c(c18_2w6c:c20_4w6), 
               names_to = "EFA", 
               values_to = "proportion") %>%
  mutate(EFA = ifelse(EFA == "c18_3w3", paste0("18:3", '\u03C9', "3"), EFA),
         EFA = ifelse(EFA == "c18_4w3", paste0("18:4", '\u03C9', "3"), EFA),
         EFA = ifelse(EFA == "c20_5w3", paste0("20:5", '\u03C9', "3"), EFA),
         EFA = ifelse(EFA == "c22_6w3", paste0("22:6", '\u03C9', "3"), EFA),
         EFA = ifelse(EFA == "c18_2w6c", paste0("18:2", '\u03C9', "6"), EFA),
         EFA = ifelse(EFA == "c20_4w6", paste0("20:4", '\u03C9', "6"), EFA)) %>%
  ggplot(aes(tourist_season, proportion, fill = EFA)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
  #geom_jitter(aes(color = EFA), shape = 21) +
  # scale_color_manual(values = viridis(6), 
  #                    name = "EFA") +
  scale_fill_manual(values = viridis(6),
                    name = "Tourist Season") +
  ylab("Proportion") +
  xlab("") +
  facet_wrap(~stt) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"))

ggsave(filename = "efa_boxplot.png", plot = efa_boxplot, 
       device = "png", path = "../figures_tables", 
       width = 12, height = 7, units = "in")


# 4. Multivariate Analysis ------------------------------------------------

efa_nmds <- metaMDS(comm = (fatty_acid_prop[, c(20, 22, 23, 25, 27, 30)]), 
                           distance = "bray", autotransform = TRUE, k = 2)

efa_nmds

adonis2(formula = fatty_acid_prop[, c(20, 22, 23, 25, 27, 30)] ~ stt+tourist_season, 
        data = fatty_acid_prop, method = "bray", by = "margin", 
        sqrt.dist = TRUE)

summary(simper(fatty_acid_prop[, c(20, 22, 23, 25, 27, 30)], 
               group = fatty_acid_prop$tourist_season,
               permutations = 999))

# Pull scores from NMDS and add site data
data_scores <- as.data.frame(scores(x = efa_nmds))
data_scores$stt <- fatty_acid_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(stt)
data_scores$tourist_season <- fatty_acid_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(tourist_season)

data_scores$stt <- fatty_acid_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(stt)

# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = efa_nmds, display = "species"))
species_scores$species <- rownames(species_scores)

# Plot the NMDS
efa_nmds_plot <- ggplot() +
  geom_point(data = data_scores,
             aes(x = NMDS1, y = NMDS2, 
                 fill = stt, shape = tourist_season), 
             size = 10, stroke = 2, color = "grey60", alpha = 0.8) +
  scale_shape_manual(values = c(21, 23), name = "Tourist Season") +
  scale_fill_manual(values = plasma(30)[c(5, 19)], 
                    name = "STT", labels = c("Centralized", "Decentralized")) +
  geom_text_repel(data = species_scores %>%
                    mutate(species = ifelse(species == "c18_3w3", paste0("18:3", '\u03C9', "3"), species),
                           species = ifelse(species == "c18_4w3", paste0("18:4", '\u03C9', "3"), species),
                           species = ifelse(species == "c20_5w3", paste0("20:5", '\u03C9', "3"), species),
                           species = ifelse(species == "c22_6w3", paste0("22:6", '\u03C9', "3"), species),
                           species = ifelse(species == "c18_2w6c", paste0("18:2", '\u03C9', "6"), species),
                           species = ifelse(species == "c20_4w6", paste0("20:4", '\u03C9', "6"), species)), 
                  aes(x = NMDS1, y = NMDS2, label = species), size = 9) + 
  annotate("label", x = -0.25, y = -0.3, size = 10,
           label = paste("Stress: ", round(efa_nmds$stress, digits = 3))) +
  guides(fill = guide_legend(override.aes = list(size=10,
                                                 color = plasma(30)[c(5, 19)])),
         shape = guide_legend(override.aes = list(size=10))) + 
  theme_minimal() +
  theme(legend.position = "right",
        legend.key=element_blank(),
        text = element_text(size = 24))

efa_nmds_plot

ggsave(filename = "efa_nmds.png", plot = efa_nmds_plot, 
       device = "png", path = "../figures_tables", 
       width = 16, height = 12, units = "in")

nmds_stt_boxplot <- data_scores %>%
  select(NMDS1, NMDS2, stt, tourist_season) %>%
  pivot_longer(cols = c(NMDS1, NMDS2), names_to = "NMDS", values_to = "scores") %>%
  ggplot(aes(NMDS, scores, fill = stt)) +
  geom_boxplot(alpha = 0.33, width = 0.2) +
  geom_jitter(aes(color = stt), size = 2, width = 0.25) +
  scale_color_manual(values = plasma(30)[c(5, 19)], name = "STT", labels = c("Centralized", "Decentralized")) +
  scale_fill_manual(values = plasma(30)[c(5, 19)], name = "STT", labels = c("Centralized", "Decentralized")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"),
        legend.title = element_text(size = 20))

nmds_tourist_season_boxplot <- data_scores %>%
  select(NMDS1, NMDS2, stt, tourist_season) %>%
  pivot_longer(cols = c(NMDS1, NMDS2), names_to = "NMDS", values_to = "scores") %>%
  ggplot(aes(NMDS, scores, fill = tourist_season)) +
  geom_boxplot(alpha = 0.33, width = 0.2) +
  geom_jitter(aes(color = tourist_season), size = 2, width = 0.25) +
  scale_color_manual(values = viridis(30)[c(5, 19)], name = "Tourist Season") +
  scale_fill_manual(values = viridis(30)[c(5, 19)], name = "Tourist Season") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"),
        legend.title = element_text(size = 20))

Anova(lm(NMDS2 ~ stt*tourist_season, 
         data = data_scores), type = "II")

Anova(lm(NMDS1 ~ stt*tourist_season, 
         data = data_scores), type = "II")

ggsave(filename = "efa_nmds_stt_boxplot.png", plot = nmds_stt_boxplot, 
       device = "png", path = "../figures_tables", 
       width = 10, height = 8, units = "in")

ggsave(filename = "efa_nmds_tourist_season_boxplot.png", plot = nmds_tourist_season_boxplot, 
       device = "png", path = "../figures_tables", 
       width = 10, height = 8, units = "in")

arranged_plots <- ggarrange(plotlist = list(nmds_stt_boxplot, nmds_tourist_season_boxplot),
                            ncol = 1, labels = "AUTO")

ggsave(filename = "combined_efa_nmds_boxplot.png", plot = arranged_plots, 
       device = "png", path = "../figures_tables", 
       width = 10, height = 12, units = "in")
