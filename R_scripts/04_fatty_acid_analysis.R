## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script is a piece of the 
## overarching analysis and looks mostly at relating 
## periphyton fatty acid profiles with wastewater treatment infrastructure
## and a sampling point's temporal position to the tourist 
## season. Questions about this script should be directed to 
## Michael F. Meyer (michael.f.meyer@wsu.edu). 


# 0. Define a function for permutational analysis -------------------------


permute_data_analytics <- function(data, metric, full_model, metric_plot_title, transform_response){
  for(i in 1:5000){
    
    # First permute the response variable. The variable is supplied by the user.
    permuted_data <- data %>%
      ungroup() %>%
      select(paste(metric)) %>%
      rename("permuted_metric" = paste(metric)) %>%
      sample_frac(size = 1) %>%
      as.vector() %>%
      cbind(., data.frame(data))
    
    # If the user has specified to log-transform the variable, this step will 
    # actually perform the log-transform. If not, then  
    if(transform_response == "log10"){
      permuted_model <- Anova(lm(log10(permuted_metric) ~ stt * tourist_season,
                                 data = permuted_data), type = "II")
    }else if(transform_response == "asin_sqrt"){
      permuted_model <- Anova(lm(asin(sqrt(permuted_metric)) ~ stt * tourist_season,
                                 data = permuted_data), type = "II")
    }else if(transform_response == "none"){
      permuted_model <- Anova(lm(permuted_metric ~ stt * tourist_season,
                                 data = permuted_data), type = "II")
    }else{
      cat("You entered in the wrong transformation")
      stop()
    }
    
    # If this iteration is the first, then function creates two repos for the 
    # p-value and r-squared values. Note that this step requires the broom package
    # be installed. 
    if(i == 1){
      tidy_repo <- tidy(permuted_model) %>%
        filter(term != "Residuals")
    } else {
      tidy_repo <- rbind(tidy_repo, tidy(permuted_model)) %>%
        filter(term != "Residuals")
    }
  }
  
  # This step removes all intercept coefficients from the repo. 
  tidy_full_model <- tidy(full_model) %>%
    filter(term != "Residuals")
  
  # Add facet labels for plot
  tidy_repo_formatted <- tidy_repo %>%
    mutate(facet_label = ifelse(test = term == "stt", 
                                yes = paste("STT \n(", nrow(filter(tidy_repo, 
                                                                   term == "stt" & 
                                                                     statistic >= tidy_full_model[tidy_full_model$term == "stt", "statistic"]$statistic))/5000*100,
                                            "% of models had higher F values)", sep = ""),
                                no = "something wrong"), 
           facet_label = ifelse(test = term == "tourist_season", 
                                yes = paste("Tourist Season \n(", nrow(filter(tidy_repo, 
                                                                              term == "tourist_season" & 
                                                                                statistic >= tidy_full_model[tidy_full_model$term == "tourist_season", "statistic"]$statistic))/5000*100,
                                            "% of models had higher F values)", sep = ""),
                                no = facet_label),
           facet_label = ifelse(test = term == "stt:tourist_season", 
                                yes = paste("Tourist Season:STT \n(", nrow(filter(tidy_repo, 
                                                                                  term == "stt:tourist_season" & 
                                                                                    statistic >= tidy_full_model[tidy_full_model$term == "stt:tourist_season", "statistic"]$statistic))/5000*100,
                                            "% of models had higher F values)", sep = ""),
                                no = facet_label))
  
  tidy_full_model_formatted <- tidy_full_model %>%
    mutate(facet_label = ifelse(test = term == "stt", 
                                yes = paste("STT \n(", nrow(filter(tidy_repo, 
                                                                   term == "stt" & 
                                                                     statistic >= tidy_full_model[tidy_full_model$term == "stt", "statistic"]$statistic))/5000*100,
                                            "% of models had higher F values)", sep = ""),
                                no = "something wrong"), 
           facet_label = ifelse(test = term == "tourist_season", 
                                yes = paste("Tourist Season \n(", nrow(filter(tidy_repo, 
                                                                              term == "tourist_season" & 
                                                                                statistic >= tidy_full_model[tidy_full_model$term == "tourist_season", "statistic"]$statistic))/5000*100,
                                            "% of models had higher F values)", sep = ""),
                                no = facet_label),
           facet_label = ifelse(test = term == "stt:tourist_season", 
                                yes = paste("Tourist Season:STT \n(", nrow(filter(tidy_repo, 
                                                                                  term == "stt:tourist_season" & 
                                                                                    statistic >= tidy_full_model[tidy_full_model$term == "stt:tourist_season", "statistic"]$statistic))/5000*100,
                                            "% of models had higher F values)", sep = ""),
                                no = facet_label))
  
  # This step plots the p-value histogram figure. 
  permuted_plot <- ggplot() +
    geom_histogram(data = tidy_repo_formatted,
                   aes(statistic), bins = 40, color = "white") +
    facet_grid(~facet_label) +
    geom_vline(data = tidy_full_model_formatted,
               aes(xintercept = statistic), linetype = "dashed", size = 2, color = viridis(8)[4]) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    xlab("F-value") +
    ylab("Frequency") +
    ggtitle(paste(metric_plot_title)) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 18),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20))
  
  # The two resulting figures are returned as a list. 
  return(list(permuted_plot))
}


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
library(broom)
library(emmeans)
library(rstatix)


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
  dplyr::select(-concentration, -total_fatty_acids) %>%
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
                     name = "Tourism Season") +
  scale_fill_manual(values = viridis(20)[c(5, 14)], 
                    name = "Tourism Season") +
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

pufa_lm <- Anova(lm(asin(sqrt(PUFA)) ~ stt*tourist_season, 
                    data = fatty_acid_type_data), 
                 type = "II")

pufa_permutation <- permute_data_analytics(data = fatty_acid_type_data, 
                                         metric = "PUFA", 
                                         full_model = pufa_lm, 
                                         metric_plot_title = "PUFA ~ STT * Tourist Season", 
                                         transform_response = "asin_sqrt")

stt_bf_pufa <- fatty_acid_type_data %>% 
  group_by(stt) %>%
  emmeans_test(asin(sqrt(PUFA)) ~ tourist_season, p.adjust.method = "bonferroni") 

ts_bf_pufa <- fatty_acid_type_data %>% 
  group_by(tourist_season) %>%
  emmeans_test(asin(sqrt(PUFA)) ~ stt, p.adjust.method = "bonferroni") 

safa_lm <- Anova(lm(asin(sqrt(SAFA)) ~ stt*tourist_season, 
                    data = fatty_acid_type_data), 
                 type = "II")

safa_permutation <- permute_data_analytics(data = fatty_acid_type_data, 
                                           metric = "SAFA", 
                                           full_model = safa_lm, 
                                           metric_plot_title = "SAFA ~ STT * Tourist Season", 
                                           transform_response = "asin_sqrt")

stt_bf_safa <- fatty_acid_type_data %>% 
  group_by(stt) %>%
  emmeans_test(asin(sqrt(SAFA)) ~ tourist_season, p.adjust.method = "bonferroni") 

ts_bf_safa <- fatty_acid_type_data %>% 
  group_by(tourist_season) %>%
  emmeans_test(asin(sqrt(SAFA)) ~ stt, p.adjust.method = "bonferroni") 

mufa_lm <- Anova(lm(asin(sqrt(MUFA)) ~ stt*tourist_season, 
                    data = fatty_acid_type_data), 
                 type = "II")

mufa_permutation <- permute_data_analytics(data = fatty_acid_type_data, 
                                           metric = "MUFA", 
                                           full_model = mufa_lm, 
                                           metric_plot_title = "MUFA ~ STT * Tourist Season", 
                                           transform_response = "asin_sqrt")

stt_bf_mufa <- fatty_acid_type_data %>% 
  group_by(stt) %>%
  emmeans_test(asin(sqrt(MUFA)) ~ tourist_season, p.adjust.method = "bonferroni") 

ts_bf_mufa <- fatty_acid_type_data %>% 
  group_by(tourist_season) %>%
  emmeans_test(asin(sqrt(MUFA)) ~ stt, p.adjust.method = "bonferroni") 

## Save Figures 

permutation_plots <- c(pufa_permutation, safa_permutation, mufa_permutation)

permutation_plot_names <- c("pufa_permutation", "safa_permutation", "mufa_permutation")

walk2(.x = permutation_plots,
      .y = permutation_plot_names,
      .f = ~ ggsave(filename = paste(.y, "_histogram.png", sep = ""), 
                    plot = .x, 
                    device = "png", path = "../figures_tables", 
                    width = 16, height = 8, units = "in"))

anova_table <- list(clean_names(data.frame(pufa_lm)), 
                    clean_names(data.frame(safa_lm)), 
                    clean_names(data.frame(mufa_lm)))

anova_table_names <- c("pufa", "safa", "mufa")

walk2(.x = anova_table,
      .y = anova_table_names,
      .f = ~ write.csv(file = paste("../figures_tables/", .y, "_anova_table.csv", sep = ""), 
                       x = .x, row.names = TRUE))

stt_bf <- rbind(stt_bf_pufa, stt_bf_safa, stt_bf_mufa) %>%
  rename("variable" = ".y.")

write.csv(x = stt_bf, 
          file = "../figures_tables/stt_fa_stoich_results.csv", 
          row.names = FALSE)

ts_bf <- rbind(ts_bf_pufa, ts_bf_safa, ts_bf_mufa) %>%
  rename("variable" = ".y.")

write.csv(x = ts_bf, 
          file = "../figures_tables/ts_bf_fa_results.csv", 
          row.names = FALSE)


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
  ggplot(aes(EFA, proportion, fill = tourist_season)) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) +
  geom_jitter(aes(color = tourist_season), width = 0.2) +
  scale_color_manual(values = viridis(20)[c(5, 14)], 
                     name = "Tourism Season") +
  scale_fill_manual(values = viridis(20)[c(5, 14)], 
                    name = "Tourism Season") +
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
       width = 16, height = 7, units = "in")


# 4. Multivariate Analysis ------------------------------------------------

## First try NMDS with fatty acid data

efa_nmds <- metaMDS(comm = (fatty_acid_prop[, c(20, 22, 23, 25, 27, 30)]), 
                           distance = "bray", autotransform = TRUE, k = 2)

efa_nmds

permanova_results <- adonis2(formula = fatty_acid_prop[, c(20, 22, 23, 25, 27, 30)] ~ stt+tourist_season, 
                             data = fatty_acid_prop, method = "bray", by = "margin", permutations = 4999,
                             sqrt.dist = TRUE)

write.csv(x = permanova_results, file = "../figures_tables/fatty_acid_permanova.csv")

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
  scale_shape_manual(values = c(21, 23), name = "Tourism Season") +
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

## Now try again with principal coorindate analysis

fatty_acid_pcoa <- pcoa(vegdist(x = fatty_acid_prop[, c(20, 22, 23, 25, 27, 30)], 
                                method = "bray"))

diag_dir <- diag(c(1,1))

fatty_acid_pcoa$vectors[,c(1,2)] <- fatty_acid_pcoa$vectors[,c(1,2)] %*% diag_dir

eignenvalues_periphyton <- fatty_acid_pcoa$vectors[,c(1,2)] %>%
  data.frame() %>%
  clean_names() %>%
  cbind(., fatty_acid_prop[, c(3,4)])

n <- nrow(fatty_acid_prop)
standardize_points <- scale(eignenvalues_periphyton[, c(1:2)], center = TRUE, scale = TRUE)
covariance_matrix <- cov(fatty_acid_prop[, c(20, 22, 23, 25, 27, 30)], standardize_points)
species_scores <- covariance_matrix %*% diag((fatty_acid_pcoa$values$Eigenvalues[c(1:2)]/(n-1))^(-0.5))
colnames(species_scores) <- colnames(fatty_acid_pcoa$vectors[,c(1:2)])

species_scores_formatted <- species_scores %>%
  data.frame() %>%
  rename("axis_1" = "Axis.1",
         "axis_2" = "Axis.2") %>%
  # mutate(axis_2 = axis_2 * max(abs(eignenvalues_periphyton$axis_2)),
  #        axis_1 = axis_1 * max(abs(eignenvalues_periphyton$axis_1))) %>%
  rownames_to_column(var = "species")

pcoa_biplot <- ggplot() +
  geom_point(data = eignenvalues_periphyton,
             aes(x = axis_1, y = axis_2, 
                 fill = stt, shape = tourist_season), 
             size = 10, stroke = 2, color = "grey60", alpha = 0.8) +
  scale_shape_manual(values = c(21, 23), name = "Tourist Season") +
  scale_fill_manual(values = plasma(30)[c(5, 19)], 
                    name = "STT", labels = c("Centralized", "Decentralized")) +
  geom_text_repel(data =  species_scores_formatted %>%
                    mutate(species = ifelse(species == "c18_3w3", paste0("18:3", '\u03C9', "3"), species),
                           species = ifelse(species == "c18_4w3", paste0("18:4", '\u03C9', "3"), species),
                           species = ifelse(species == "c20_5w3", paste0("20:5", '\u03C9', "3"), species),
                           species = ifelse(species == "c22_6w3", paste0("22:6", '\u03C9', "3"), species),
                           species = ifelse(species == "c18_2w6c", paste0("18:2", '\u03C9', "6"), species),
                           species = ifelse(species == "c20_4w6", paste0("20:4", '\u03C9', "6"), species)), 
                  aes(x = axis_1, y = axis_2, label = species), 
                  size = 9) + 
  xlab(paste("PCo 1 (", round(fatty_acid_pcoa$values[1,2]*100, 2), "% of Variance)", sep = "")) +
  ylab(paste("PCo 2 (", round(fatty_acid_pcoa$values[2,2]*100, 2), "% of Variance)", sep = "")) +
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

manova_model <- manova(fatty_acid_pcoa$vectors[, c(1:2)] ~ eignenvalues_periphyton$stt * eignenvalues_periphyton$tourist_season)

summary(manova_model)

library(MVN)

mvn(data = manova_model$residuals, 
    mvnTest = "hz", multivariatePlot = "qq")

Anova(lm(axis_2 ~ stt*tourist_season, 
         data = eignenvalues_periphyton), type = "II")

Anova(lm(axis_1 ~ stt*tourist_season, 
         data = eignenvalues_periphyton), type = "II")

pcoa_stt_boxplot <- eignenvalues_periphyton %>%
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

pcoa_tourist_season_boxplot <- eignenvalues_periphyton %>%
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

Anova(lm(axis_2 ~ stt*tourist_season, 
         data = eignenvalues_periphyton), type = "II")

Anova(lm(axis_1 ~ stt*tourist_season, 
         data = eignenvalues_periphyton), type = "II")

ggsave(filename = "efa_pcoa_stt_boxplot.png", plot = pcoa_stt_boxplot, 
       device = "png", path = "../figures_tables", 
       width = 10, height = 8, units = "in")

ggsave(filename = "efa_pcoa_tourist_season_boxplot.png", plot = pcoa_tourist_season_boxplot, 
       device = "png", path = "../figures_tables", 
       width = 10, height = 8, units = "in")

arranged_plots <- ggarrange(plotlist = list(pcoa_stt_boxplot, pcoa_tourist_season_boxplot),
                            ncol = 1, labels = "AUTO")

ggsave(filename = "combined_efa_pcoa_boxplot.png", plot = arranged_plots, 
       device = "png", path = "../figures_tables", 
       width = 10, height = 12, units = "in")


arranged_plots <- ggarrange(pcoa_biplot, 
                            ggarrange(plotlist = list(pcoa_stt_boxplot, pcoa_tourist_season_boxplot),
                                      ncol = 1, labels = c("B", "C"), font.label = list(size = 24)),
                            ncol = 2, widths = c(2,1), labels = c("A"), font.label = list(size = 24))

ggsave(filename = "combined_efa_pcoa_boxplot.png", plot = arranged_plots, 
       device = "png", path = "../figures_tables", 
       width = 20, height = 10, units = "in")
