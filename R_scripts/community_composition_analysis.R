library(tidyverse)
library(vegan)
library(ggrepel)
library(viridis)

periphyton_orig <- read.csv(file = "../cleaned_data/periphyton.csv",
                            header = TRUE)

space_time_population <- read.csv(file = "../cleaned_data/temporal_inverse_distance_weighted_population_metrics.csv",
                                  header = TRUE) %>%
  mutate(month = ifelse(grepl(pattern = "june", x = month_scalar), yes = "June", no = NA),
         month = ifelse(grepl(pattern = "july", x = month_scalar), yes = "July", no = month),
         month = ifelse(grepl(pattern = "august", x = month_scalar), yes = "August", no = month),
         month = ifelse(grepl(pattern = "september", x = month_scalar), yes = "September", no = month)) %>%
  select(site = sampling_site,
         distance_weighted_population,
         scaled_idw_population,
         month) 

samples_full <- read.csv("../cleaned_data/ppcp.csv") %>%
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
  mutate(month = ifelse(sampling_event %in% c("JUNE_1", "JUNE_2"), "June", NA),
         month = ifelse(sampling_event %in% c("JULY_1", "JULY_2"), "July", month),
         month = ifelse(sampling_event %in% c("AUGUST_1", "AUGUST_2"), "August", month),
         month = ifelse(sampling_event %in% c("SEPTEMBER_1", "SEPTEMBER_2"), "September", month)) %>%
  mutate(stt = ifelse(SITE %in% c("WF", "HO", "FLBS", "YB"), "Centralized", "Decentralized")) %>%
  group_by(stt, SITE, month, sampling_event, PPCP) %>%
  summarize(mean_conc = mean(CONCENTRATION, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(stt, SITE, month, sampling_event) %>%
  summarize(total_conc = sum(mean_conc, na.rm = TRUE)) %>%
  group_by(stt, SITE, month) %>%
  summarize(mean_conc = mean(total_conc, na.rm = TRUE),
            var_conc = var(total_conc, na.rm = TRUE)) %>%
  mutate(tourist_season = ifelse(month %in% c("September"), "Out of Season", "In Season")) %>%
  rename("site" = "SITE")


summary(periphyton_orig)

periphyton_prop <- periphyton_orig %>%
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
  inner_join(., y = space_time_population) %>%
  inner_join(., y = samples_full) %>%
  ungroup() %>%
  mutate(month = factor(month, levels = c("July", "August", "September")))

summary(periphyton_prop)

periphyton_prop %>%
  pivot_longer(cols = chlorophyta:diatom, names_to = "taxon" ,values_to = "proportion") %>%
  group_by(stt, month, taxon) %>%
  summarize(mean_proportion = mean(proportion)) %>%
  ungroup() %>%
  group_by(stt, month) %>%
  mutate(sum_mean_prop = sum(mean_proportion),
         proportion = mean_proportion/sum_mean_prop) %>%
  ggplot(aes(x = month, y = proportion, fill = taxon)) +
  geom_bar(stat = "identity") +
  ylab("Mean Proportion") +
  xlab("") +
  scale_fill_manual(values = viridis(40)[c(34,28,19,16,5)], name = NULL) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~stt) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"),
        panel.grid = element_blank())


periphyton_prop %>%
  #filter(site %in% c("WF", "BD", "WB", "FLBS", "YB")) %>%
ggplot(aes(log10(scaled_idw_population), chlorophyta/diatom)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  facet_grid(~month)

periphyton_prop %>%
  pivot_longer(cols = c(chlorophyta:diatom), names_to = "taxon", values_to = "proportion") %>%
  filter(!taxon  %in% c("crysophyta", "cryptophyta")) %>%
  mutate(tourist_season = ifelse(month %in% c("July", "August"), "In Season", "Out of Season")) %>%
  #filter(site %in% c("WF", "BD", "WB", "FLBS", "YB")) %>%
  ggplot(aes(month, proportion, fill = stt, color = stt)) +
  geom_boxplot(alpha = 0.33, width = 0.33, outlier.alpha = 0) +
  geom_jitter() +
  ylab("Proportion") +
  xlab("") +
  geom_jitter(size = 2) +
  scale_color_manual(values = plasma(30)[c(5, 19)], name = NULL) +
  scale_fill_manual(values = plasma(30)[c(5, 19)], name = NULL) +
  facet_grid(~taxon) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"))


car::Anova(lm(NMDS2 ~ stt*tourist_season, 
           data = data_scores), type = "II")


summary(lm(chlorophyta/(diatom) ~ log10(scaled_idw_population), data = periphyton_prop %>% filter(stt == "Decentralized")))

library(lme4)
library(merTools)
library(lmerTest)

chloro_mixed <- lmer(chlorophyta/diatom ~ log10(scaled_idw_population)*stt + (1|site), data = periphyton_prop)

plotREsim(REsim(chloro_mixed))

ranova(chloro_mixed)
anova(chloro_mixed, type = 2)

periphyton_nmds <- metaMDS(comm = (periphyton_prop[, c(3, 6, 7)]), 
                           distance = "bray", autotransform = TRUE, k = 3)

periphyton_nmds

perm <- how(nperm = 999)
setBlocks(perm) <- with(periphyton_prop, site)
adonis2(formula = periphyton_prop[, c(3, 6, 7)] ~ stt+tourist_season, 
        data = periphyton_prop, method = "bray", by = "margin", 
        sqrt.dist = TRUE)

# Pull scores from NMDS and add site data
data_scores <- as.data.frame(scores(x = periphyton_nmds))
data_scores$stt <- periphyton_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(stt)
data_scores$tourist_season <- periphyton_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(tourist_season)
data_scores$scaled_idw_population <- periphyton_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(scaled_idw_population)

data_scores$mean_conc <- periphyton_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(mean_conc)

data_scores$stt <- periphyton_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(stt)

# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = periphyton_nmds, display = "species"))
species_scores$species <- rownames(species_scores)

# Plot the NMDS
periphyton_nmds_plot <- ggplot() +
  geom_point(data = data_scores,
             aes(x = NMDS2, y = NMDS3, 
                 color = stt, shape = tourist_season)) +
  scale_size_continuous(range = c(5, 10), guide = FALSE) +
  geom_text_repel(data = species_scores, 
                  aes(x = NMDS2, y = NMDS3, label = species), 
                  size = 9, parse = TRUE) + 
  annotate("label", x = -0.25, y = 0.125, size = 10,
           label = paste("Stress: ", round(periphyton_nmds$stress, digits = 3))) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.key=element_blank(),
        text = element_text(size = 24))

periphyton_nmds_plot

ggsave(filename = "../figures_tables/nmds_periphtyon.png", plot = periphyton_nmds_plot, device = "png", width = 10, height = 8, units = "in")


### Overall number of cells

periphyton_abundance <- periphyton_orig %>%
  mutate(volume_sampled = volume_ul * volume_rep) %>%
  dplyr::select(-volume_ul, -volume_rep, -chlorophyta_filamnets, -date, -notes) %>%
  pivot_longer(cols = c(diatom, chlorophyta, crysophyta, cryptophyta, cyanobacteria),
               names_to = "taxon",
               values_to = "abundance") %>%
  mutate(abundance_per_ul = abundance/volume_sampled) %>%
  dplyr::select(-abundance, -rep) %>%
  group_by(site, month, taxon) %>%
  summarize(average_prop = mean(abundance_per_ul)) %>%
  pivot_wider(names_from = taxon, values_from = average_prop) %>%
  inner_join(., y = space_time_population) %>%
  inner_join(., y = samples_full) %>%
  ungroup() %>%
  mutate(month = factor(month, levels = c("July", "August", "September")))


summary(periphyton_abundance)

ggplot(periphyton_abundance, aes(log10(mean_conc), cyanobacteria/chlorophyta, color = month)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE)

model <- lmer(diatom/(chlorophyta+cyanobacteria) ~ month*log10(scaled_idw_population) + (1|site), data = periphyton_prop)

summary(model)
anova(model)

model <- lm(diatom/(chlorophyta+cyanobacteria) ~ 0 + month*log10(scaled_idw_population), data = periphyton_prop)

summary(model)


periphyton_nmds <- metaMDS(comm = periphyton_abundance[, c(3, 6, 7)], distance = "bray")

periphyton_nmds

adonis2(formula = (periphyton_abundance[, c(3, 6, 7)]) ~ stt+month+log10(mean_conc), 
        data = periphyton_abundance, method = "bray", by = "margin")

# Pull scores from NMDS and add site data
data_scores <- as.data.frame(scores(x = periphyton_nmds))
data_scores$site <- periphyton_abundance %>%
  #filter(!site %in% c("YB","WB", "BB")) %>%
  pull(site)
data_scores$month <- periphyton_abundance %>%
  #filter(!site %in% c("YB","WB", "BB")) %>%
  pull(month)
data_scores$scaled_idw_population <- periphyton_abundance %>%
  #filter(!site %in% c("YB","WB", "BB")) %>%
  pull(scaled_idw_population)
data_scores$mean_conc <- periphyton_abundance %>%
  #filter(!site %in% c("YB","WB", "BB")) %>%
  pull(mean_conc)
data_scores$stt <- periphyton_abundance %>%
  #filter(!site %in% c("YB","WB", "BB")) %>%
  pull(stt)


# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = periphyton_nmds, display = "species"))
species_scores$species <- rownames(species_scores)

# Plot the NMDS
periphyton_nmds_plot <- ggplot() +
  geom_point(data = data_scores,
            aes(x = NMDS1, y = NMDS2, 
                color = stt, size = mean_conc)) +
  scale_size_continuous(range = c(3, 10), guide = FALSE) +
  geom_text_repel(data = species_scores, 
                  aes(x = NMDS1, y = NMDS2, label = species), 
                  size = 9, parse = TRUE) + 
  annotate("label", x = -0.5, y = -0.25, size = 10,
           label = paste("Stress: ", round(periphyton_nmds$stress, digits = 3))) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.key=element_blank(),
        text = element_text(size = 24))

periphyton_nmds_plot

ggsave(filename = "../figures_tables/nmds_periphtyon_abundance.png", plot = periphyton_nmds_plot, device = "png", width = 14, height = 12, units = "in")


