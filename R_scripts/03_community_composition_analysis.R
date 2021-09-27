## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script is a piece of the 
## overarching analysis and looks mostly at relating 
## periphyton community compositions with wastewater treatment infrastructure
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
library(ape)
library(janitor)
library(MVN)

periphyton_orig <- read.csv(file = "../cleaned_data/periphyton_abundance.csv",
                            header = TRUE)


# 2. Convert data into proportions ----------------------------------------

periphyton_prop <- periphyton_orig %>%
  select(-volume_ul, -volume_rep, -chlorophyta_filamentts, -date, -notes) %>%
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
  ungroup() %>%
  mutate(month = factor(month, levels = c("July", "August", "September"))) %>%
  mutate(stt = ifelse(test = site %in% c("WF", "YB", "FLBS", "HO"), 
                      yes = "centralized", no = "decentralized"),
         tourist_season = ifelse(test = month %in% c("June", "July", "August"), 
                                 yes = "In Season", no = "Out of Season"))

summary(periphyton_prop)

peri_prop_table <- periphyton_prop %>%
  select(tourist_season, stt, chlorophyta, cyanobacteria, diatom) %>%
  pivot_longer(cols = c(chlorophyta, cyanobacteria, diatom),
               names_to = "taxon", values_to = "proportion") %>%
  group_by(tourist_season, stt, taxon) %>%
  summarize(mean_prop = mean(proportion),
            sd_prop = sd(proportion),
            cv_prop = sd_prop/mean_prop) %>%
  pivot_wider(names_from = "stt", 
              values_from = c("mean_prop", "sd_prop", "cv_prop")) %>%
  select(taxon, tourist_season, mean_prop_centralized:cv_prop_decentralized) %>%
  arrange(taxon)

write.csv(x = peri_prop_table, 
          file = "../figures_tables/periphtyon_summary_stats.csv", 
          row.names = FALSE)

# 3. Univariate Analysis --------------------------------------------------

stt_labels <- c("Centralized", "Decentralized")
names(stt_labels) <- c("centralized", "decentralized")

stacked_bar <- periphyton_prop %>%
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
  facet_wrap(~stt, labeller = labeller(stt = stt_labels)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"),
        panel.grid = element_blank())

ggsave(filename = "periphyton_stacked_barchart.png", plot = stacked_bar, 
       device = "png", path = "../figures_tables", 
       width = 12, height = 7, units = "in")

peri_boxplot <- periphyton_prop %>%
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
  scale_color_manual(values = plasma(30)[c(5, 19)], name = "STT", labels = c("Centralized", "Decentralized")) +
  scale_fill_manual(values = plasma(30)[c(5, 19)], name = "STT", labels = c("Centralized", "Decentralized")) +
  facet_grid(~taxon) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"))

ggsave(filename = "periphyton_boxplot.png", plot = peri_boxplot, 
       device = "png", path = "../figures_tables", 
       width = 16, height = 7, units = "in")


# 4. Multivariate Analysis ------------------------------------------------

periphyton_nmds <- metaMDS(comm = (periphyton_prop[, c(3, 6, 7)]), 
                           distance = "bray", autotransform = TRUE, k = 2)

periphyton_nmds

permanova_results <- adonis2(formula = (periphyton_prop[, c(3, 6, 7)]) ~ stt+tourist_season, 
                             data = periphyton_prop, 
                             method = "bray", by = "margin", 
                             sqrt.dist = TRUE, 
                             permutations = 4999)

write.csv(x = permanova_results, file = "../figures_tables/periphyton_permanova.csv")

# Pull scores from NMDS and add site data
data_scores <- as.data.frame(scores(x = periphyton_nmds))
data_scores$stt <- periphyton_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(stt)
data_scores$tourist_season <- periphyton_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(tourist_season)

data_scores$stt <- periphyton_prop %>%
  #filter(!site %in% c("WB", "BB", "YB")) %>%
  pull(stt)

# Pull species scores from NMDS
species_scores <- as.data.frame(scores(x = periphyton_nmds, display = "species"))
species_scores$species <- rownames(species_scores)

# Plot the NMDS
periphyton_nmds_plot <- ggplot() +
  geom_point(data = data_scores,
             aes(x = NMDS1, y = NMDS2, 
                 fill = stt, shape = tourist_season), 
             size = 10, stroke = 2, color = "grey60", alpha = 0.8) +
  scale_shape_manual(values = c(21, 23), name = "Tourism Season") +
  scale_fill_manual(values = plasma(30)[c(5, 19)], 
                    name = "STT", labels = c("Centralized", "Decentralized")) +
  geom_text_repel(data = species_scores, 
                  aes(x = NMDS1, y = NMDS2, label = species), 
                  size = 9, parse = TRUE) + 
  annotate("label", x = -0.25, y = -0.2, size = 10,
           label = paste("Stress: ", round(periphyton_nmds$stress, digits = 3))) +
  guides(fill = guide_legend(override.aes = list(size=10,
                                                 color = plasma(30)[c(5, 19)])),
         shape = guide_legend(override.aes = list(size=10))) + 
  theme_minimal() +
  theme(legend.position = "right",
        legend.key=element_blank(),
        text = element_text(size = 24))

periphyton_nmds_plot

ggsave(filename = "periphyton_nmds.png", plot = periphyton_nmds_plot, 
       device = "png", path = "../figures_tables", 
       width = 16, height = 12, units = "in")

## Now try with Principal Coordinate Analysis
## This method is meant to confirm the patterns
## observed with NMDS but with a completely 
## different approach (i.e., numerical vs eigencalculations)

## This approach is NOT included in the manuscript but is included
## in the R code so to serve as an archive or how a similar but 
## different technique got us qualitatively similar results. 

periphyton_pcoa <- pcoa(vegdist(x = periphyton_prop[, c(3, 6, 7)], method = "bray"))

diag_dir <- diag(c(1,1))

periphyton_pcoa$vectors[,c(1,2)] <- periphyton_pcoa$vectors[,c(1,2)] %*% diag_dir

## This step extracts the eigenvalues 

eigenvalues_periphyton <- periphyton_pcoa$vectors[,c(1,2)] %>%
  data.frame() %>%
  clean_names() %>%
  cbind(., periphyton_prop[, c(8,9)])

## This step scales the eigenvectors by the covariance matrix
n <- nrow(periphyton_prop)
standardize_points <- scale(eigenvalues_periphyton[, c(1:2)], center = TRUE, scale = TRUE)
covariance_matrix <- cov(periphyton_prop[,c(3,6,7)], standardize_points)
species_scores <- covariance_matrix %*% diag((periphyton_pcoa$values$Eigenvalues[c(1:2)]/(n-1))^(-0.5))
colnames(species_scores) <- colnames(periphyton_pcoa$vectors[,c(1:2)])

species_scores_formatted <- species_scores %>%
  data.frame() %>%
  rename("axis_1" = "Axis.1",
         "axis_2" = "Axis.2") %>%
  mutate(axis_2 = axis_2 * max(abs(eigenvalues_periphyton$axis_2)),
         axis_1 = axis_1 * max(abs(eigenvalues_periphyton$axis_1))) %>%
  rownames_to_column(var = "taxon")

pcoa_biplot <- ggplot() +
  geom_point(data = eigenvalues_periphyton,
             aes(x = axis_1, y = axis_2, 
                 fill = stt, shape = tourist_season), 
             size = 10, stroke = 2, color = "grey60", alpha = 0.8) +
  scale_shape_manual(values = c(21, 23), name = "Tourist Season") +
  scale_fill_manual(values = plasma(30)[c(5, 19)], 
                    name = "STT", labels = c("Centralized", "Decentralized")) +
  geom_text_repel(data =  species_scores_formatted, 
                  aes(x = axis_1, y = axis_2, label = taxon), 
                  size = 9, parse = TRUE) + 
  xlab(paste("PCo 1 (", round(periphyton_pcoa$values[1,2]*100, 2), "% of Variance)", sep = "")) +
  ylab(paste("PCo 2 (", round(periphyton_pcoa$values[2,2]*100, 2), "% of Variance)", sep = "")) +
  #annotate("label", x = -0.25, y = -0.3, size = 10,
  #         label = paste("Stress: ", round(periphyton_nmds$stress, digits = 3))) +
  guides(fill = guide_legend(override.aes = list(size=10,
                                                 color = plasma(30)[c(5, 19)])),
         shape = guide_legend(override.aes = list(size=10))) + 
  theme_minimal() +
  theme(legend.position = "right",
        legend.key=element_blank(),
        text = element_text(size = 24))

pcoa_biplot

## To test for significant differences, we perform a MANOVA with the 
## 

manova_model <- manova(as.matrix(periphyton_prop[, c(3, 6, 7)]) ~ periphyton_prop$tourist_season *  periphyton_prop$stt)

summary(manova_model)

mvn(data = as.matrix(periphyton_prop[, c(3, 6, 7)]), 
    mvnTest = "mardia", multivariatePlot = "qq")

Anova(lm(axis_2 ~ stt*tourist_season, 
         data = eigenvalues_periphyton), type = "II")

Anova(lm(axis_1 ~ stt*tourist_season, 
         data = eigenvalues_periphyton), type = "II")

pcoa_stt_boxplot <- eigenvalues_periphyton %>%
  select(axis_1, axis_2, stt, tourist_season) %>%
  pivot_longer(cols = c(axis_1, axis_2), names_to = "PCo", values_to = "scores") %>%
  mutate(PCo = gsub(pattern = "axis_", replacement = "PCo ", x = PCo)) %>%
  ggplot(aes(PCo, scores, fill = stt)) +
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

pcoa_tourist_season_boxplot <- eigenvalues_periphyton %>%
  select(axis_1, axis_2, stt, tourist_season) %>%
  pivot_longer(cols = c(axis_1, axis_2), names_to = "PCo", values_to = "scores") %>%
  mutate(PCo = gsub(pattern = "axis_", replacement = "PCo ", x = PCo)) %>%
  ggplot(aes(PCo, scores, fill = tourist_season)) +
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


ggsave(filename = "pcoa_stt_boxplot.png", plot = pcoa_stt_boxplot, 
       device = "png", path = "../figures_tables", 
       width = 10, height = 8, units = "in")

ggsave(filename = "pcoa_tourist_season_boxplot.png", plot = pcoa_tourist_season_boxplot, 
       device = "png", path = "../figures_tables", 
       width = 10, height = 8, units = "in")

arranged_plots <- ggarrange(pcoa_biplot, 
                            ggarrange(plotlist = list(pcoa_stt_boxplot, pcoa_tourist_season_boxplot),
                            ncol = 1, labels = c("B", "C"), font.label = list(size = 24)),
                            ncol = 2, widths = c(1.65,1), labels = c("A"), font.label = list(size = 24))

ggsave(filename = "combined_pcoa_boxplot.png", plot = arranged_plots, 
       device = "png", path = "../figures_tables", 
       width = 20, height = 10, units = "in")
