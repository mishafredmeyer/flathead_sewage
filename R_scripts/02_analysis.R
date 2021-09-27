## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script is a piece of the 
## overarching analysis and looks mostly at relating 
## sewage indicators with wastewater treatment infrastructure
## and a sampling point's temporal position to the tourist 
## season. Questions about this script should be directed to 
## Michael F. Meyer (michael.f.meyer@wsu.edu). 

# 0. Define a function for permutation of sewage indicators ---------------

# Here, we will define a function to perform a permutational analysis of the 
# of the data. In general, this function permutes the response variable 5,000 times, 
# creates an ANOVA model with the permuted data, and extracts the F value
# from the model. Once the interation is finished, the 5,000 F values are plotted 
# on a histogram and compared to the model generated with non-permuted data. 

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


nutrients_formatted <- nutrients %>%
  rowwise() %>%
  mutate(srp = ifelse(srp == 0.4, runif(n = 1, min = 0, max = 0.4), srp),
         nh3_n = ifelse(nh3_n == 1.5, runif(n = 1, min = 0, max = 1.5), nh3_n),
         total_p = ifelse(total_p == 0.75, runif(n = 1, min = 0, max = 0.75), total_p),
         no3_no2 = ifelse(no3_no2 == 0.75, runif(n = 1, min = 0, max = 0.75), no3_no2)) %>%
  ungroup() %>%
  pivot_longer(cols = c(nh3_n:total_p), names_to = "nutrient_type", values_to = "concentration") %>%
  filter(!site %in% c("DU")) %>%
  group_by(tourist_season, stt, site, nutrient_type) %>%
  summarize(mean_concentration = mean(concentration, na.rm = TRUE))

nutrient_plot <- nutrients_formatted %>%
  ungroup() %>%
  mutate(nutrient_type = factor(nutrient_type, 
                                levels = c("nh3_n", "no3_no2", "total_n",
                                           "srp", "total_p"))) %>%
  ggplot(aes(tourist_season, (mean_concentration), fill = stt, color = stt)) +
  geom_boxplot(alpha = 0.33, width = 0.25, outlier.alpha = 0) +
  geom_jitter(width = 0.3) +
  ylab(expression(paste("Concentration (\u03BCg/L- N or P)"))) +
  scale_fill_manual(values = plasma(30)[c(5, 19)], name = "STT") +
  scale_color_manual(values = plasma(30)[c(5, 19)], name = "STT") +
  scale_y_log10() +
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
       width = 15, height = 8, units = "in")


# Ammonium Model
nh3_lm <- Anova(lm(log10(nh3_n) ~ stt * tourist_season, 
              data = nutrients_formatted_wide), 
              type = "II")

nh3_permutation <- permute_data_analytics(data = nutrients_formatted_wide, 
                                               metric = "nh3_n", 
                                               full_model = nh3_lm, 
                                               metric_plot_title = "Ammonia ~ STT * Tourist Season", 
                                               transform_response = "log10")


# Nitrate/Nitrite Model
no3_no2_lm <- Anova(lm(log10(no3_no2) ~ stt*tourist_season, 
                       data = nutrients_formatted_wide),
                    type = "II")

no3_no2_permutation <- permute_data_analytics(data = nutrients_formatted_wide, 
                                               metric = "no3_no2", 
                                               full_model = no3_no2_lm, 
                                               metric_plot_title = "Nitrate/Nitrite ~ STT * Tourist Season", 
                                               transform_response = "log10")

# SRP
srp_lm <- Anova(lm(log10(srp) ~ stt * tourist_season,
                   data = nutrients_formatted_wide), 
             type = "II")

srp_permutation <- permute_data_analytics(data = nutrients_formatted_wide, 
                                              metric = "srp", 
                                              full_model = srp_lm, 
                                              metric_plot_title = "SRP ~ STT * Tourist Season", 
                                              transform_response = "log10")

# Total Nitrogen
tn_lm <- Anova(lm(log10(total_n) ~ stt*tourist_season,
                  data = nutrients_formatted_wide), 
             type = "II")

tn_permutation <- permute_data_analytics(data = nutrients_formatted_wide, 
                                          metric = "total_n", 
                                          full_model = tn_lm, 
                                          metric_plot_title = "Total Nitrogen ~ STT * Tourist Season", 
                                          transform_response = "log10")

# Total Phosphorus
tp_lm <- Anova(lm(log10(total_p) ~ stt*tourist_season, 
                  data = nutrients_formatted_wide), 
               type = "II")

tp_permutation <- permute_data_analytics(data = nutrients_formatted_wide, 
                                         metric = "total_p", 
                                         full_model = tp_lm, 
                                         metric_plot_title = "Total Phosphorus ~ STT * Tourist Season", 
                                         transform_response = "log10")

# 2. Ash Free Dry Mass ----------------------------------------------------

stt_labels <- c("Centralized", "Decentralized")
names(stt_labels) <- c("centralized", "decentralized")

afdm_mean <- afdm %>%
  group_by(tourist_season, stt, site) %>%
  summarize(mean_afdm = mean(afdm, na.rm = TRUE))

afdm_plot <- ggplot(afdm_mean, aes(tourist_season, mean_afdm)) +
  geom_boxplot(alpha = 0.33, width = 0.25, outlier.alpha = 0) +
  geom_jitter(width = 0.3) +
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

afdm_lm <- Anova(lm(mean_afdm ~ stt*tourist_season, 
                    data = afdm_mean), 
                 type = "II")

afdm_permutation <- permute_data_analytics(data = afdm_mean, 
                                           metric = "mean_afdm", 
                                           full_model = afdm_lm, 
                                           metric_plot_title = "AFDM ~ STT * Tourist Season", 
                                           transform_response = "none")


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
  ungroup() %>%
  group_by(tourist_season, stt, site) %>%
  summarize(mean_total_branched_odd = mean(total_branched_odd, na.rm = TRUE)) %>%
  ggplot(aes(x = tourist_season, y = asin(sqrt(mean_total_branched_odd)))) +
  geom_boxplot(alpha = 0.33, width = 0.25, outlier.alpha = 0) +
  geom_jitter(width = 0.3) +
  xlab("")+
  ylab("Branched- & Odd-chain Fatty Acids") +
  facet_wrap(~stt) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.key.height = unit(0.75, "in"),
        legend.key.width = unit(0.33, "in"))

ggsave(filename = "branched_odd_chain_fatty_acid_boxplots.png", 
       plot = branched_odd_chain_fatty_acids_plot, 
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
  summarize(total_branched_odd = rowSums(across(.cols = iso_c15_0:c17_0))) %>%
  ungroup() %>%
  group_by(tourist_season, stt, site) %>%
  summarize(mean_total_branched_odd = mean(total_branched_odd))


branched_odd_chain_fatty_acids_lm <- Anova(lm(asin(sqrt(mean_total_branched_odd)) ~ stt * tourist_season,
                                              data = branched_odd_chain_fatty_acids), 
                                           type = "II")

branched_odd_fa_permutation <- permute_data_analytics(data = branched_odd_chain_fatty_acids,
                                                       metric = "mean_total_branched_odd", 
                                                       full_model = branched_odd_chain_fatty_acids_lm, 
                                                       metric_plot_title = "Branched- & Odd Chain Fatty Acids ~ STT * Tourist Season", 
                                                       transform_response = "asin_sqrt")



# 4. PPCPs ----------------------------------------------------------------

ppcp_reduced <- ppcp %>%
  filter(!site %in% c("HO", "DU")) %>%
  # filter(month != "may") %>%
  # group_by(tourist_season, stt, site, month, sampling_event) %>%
  # summarize(mean_conc = sum(concentration, na.rm = TRUE)) %>%
  group_by(tourist_season, stt, site, peri_sampling, sampling_event, ppcp) %>%
  summarize(mean_conc = mean(concentration, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(tourist_season, stt, site, peri_sampling, sampling_event) %>%
  summarize(total_concentration = sum(mean_conc, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(tourist_season, stt, site) %>%
  summarize(mean_conc = mean(total_concentration, na.rm = TRUE)) 

ppcp_plot <- ggplot(ppcp_reduced , 
       aes(tourist_season, (mean_conc))) +
  geom_boxplot(alpha = 0.33, width = 0.25, outlier.alpha = 0) +
  geom_jitter(width = 0.3) +
  facet_wrap(~stt, labeller = labeller(stt = stt_labels)) +
  ylab("Total PPCP (ng/L)") +
  scale_y_log10() +
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


ppcp_lm <- Anova(lm(log10(mean_conc) ~ stt * tourist_season,
                    data = ppcp_reduced), 
                 type = "II")

ppcp_permutation <- permute_data_analytics(data = ppcp_reduced,
                                           metric = "mean_conc", 
                                           full_model = ppcp_lm, 
                                           metric_plot_title = "PPCP ~ STT * Tourist Season", 
                                           transform_response = "log10")

# 5. TSIDW ----------------------------------------------------------------

tsidw_formatted <- tsidw_pop %>%
  group_by(tourist_season, stt, site) %>%
  summarize(mean_tsidw_pop = mean(atsidw_pop))

tsidw_plot <- ggplot(tsidw_formatted, 
                    aes(tourist_season, ((mean_tsidw_pop)))) +
  geom_boxplot(alpha = 0.33, width = 0.25, outlier.alpha = 0) +
  geom_jitter(width = 0.3) +
  scale_y_log10() +
  facet_wrap(~stt, labeller = labeller(stt = stt_labels)) +
  ylab("TSIDW Population (persons)") + 
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


# 6. Combine Figures ------------------------------------------------------

arranged_plots <- ggarrange(plotlist = list(tsidw_plot, ppcp_plot,
                                         branched_odd_chain_fatty_acids_plot, afdm_plot), 
                            ncol = 2, nrow = 2, labels = "AUTO",
                            font.label = list(size = 24))

ggsave(filename = "combined_boxplots.png", plot = arranged_plots, 
       device = "png", path = "../figures_tables", 
       width = 16, height = 12, units = "in")

permutation_plots <- c(nh3_permutation, no3_no2_permutation, srp_permutation, 
                       tn_permutation, tp_permutation, afdm_permutation,
                       branched_odd_fa_permutation, ppcp_permutation)

permutation_plot_names <- c("nh3_permutation", "no3_no2_permutation", "srp_permutation", 
                       "tn_permutation", "tp_permutation", "afdm_permutation",
                       "branched_odd_fa_permutation", "ppcp_permutation")

walk2(.x = permutation_plots,
      .y = permutation_plot_names,
      .f = ~ ggsave(filename = paste(.y, "_histogram.png", sep = ""), 
                   plot = .x, 
                   device = "png", path = "../figures_tables", 
                   width = 16, height = 8, units = "in"))

anova_table <- list(clean_names(data.frame(nh3_lm)), 
                    clean_names(data.frame(no3_no2_lm)), 
                    clean_names(data.frame(srp_lm)), 
                    clean_names(data.frame(tn_lm)), 
                    clean_names(data.frame(tp_lm)), 
                    clean_names(data.frame(afdm_lm)),
                    clean_names(data.frame(branched_odd_chain_fatty_acids_lm)), 
                    clean_names(data.frame(ppcp_lm)))

anova_table_names <- c("nh3", "no3_no2", "srp", 
                            "tn", "tp", "afdm",
                            "branched_odd_chain_fatty_acids", 
                            "ppcp")

walk2(.x = anova_table,
      .y = anova_table_names,
      .f = ~ write.csv(file = paste("../figures_tables/", .y, "_anova_table.csv", sep = ""), 
                    x = .x, row.names = TRUE))


# 7. Post revision analysis -----------------------------------------------

## Following review by the committee, we decided to add a multi-panel figure
## of sewage indicators regressed against continuous TSIDW metrics.

## Nutrients

nutrients_table <- nutrients_formatted %>%
  ungroup() %>%
  group_by(tourist_season, stt, nutrient_type) %>%
  summarize(mean_conc = mean(mean_concentration),
            sd_conc = sd(mean_concentration),
            cv_conc = sd_conc/mean_conc) %>%
  pivot_wider(names_from = "stt", 
              values_from = c("mean_conc", "sd_conc", "cv_conc")) %>%
  arrange(nutrient_type)

write.csv(x = nutrients_table, 
          file = "../figures_tables/nutrients_summary_stats.csv", 
          row.names = FALSE)

nutrients_tsidw <- inner_join(tsidw_formatted, nutrients_formatted) %>%
  rename("analyte" = "nutrient_type",
         "mean_conc" = "mean_concentration") %>%
  select(tourist_season, stt, site, analyte, mean_tsidw_pop, mean_conc)

## AFDM

afdm_table <- afdm_mean %>%
  ungroup() %>%
  group_by(tourist_season, stt) %>%
  summarize(mean_conc = mean(mean_afdm, na.rm = TRUE),
            sd_conc = sd(mean_afdm, na.rm = TRUE),
            cv_conc = sd_conc/mean_conc) %>%
  pivot_wider(names_from = "stt", 
              values_from = c("mean_conc", "sd_conc", "cv_conc"))

write.csv(x = afdm_table, 
          file = "../figures_tables/afdm_summary_stats.csv", 
          row.names = FALSE)

afdm_tsidw <- inner_join(afdm_mean, tsidw_formatted) %>%
  mutate(analyte = "afdm") %>%
  rename("mean_conc" = "mean_afdm") %>%
  select(tourist_season, stt, site, analyte, mean_tsidw_pop, mean_conc)

## Branched and Odd-Chain Fatty Acids

branched_odd_chain_fatty_acids_table <- fatty_acids %>%
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
  ungroup() %>%
  group_by(tourist_season, stt, site) %>%
  summarize(mean_total_branched_odd = mean(total_branched_odd)) %>%
  ungroup() %>%
  group_by(tourist_season, stt) %>%
  summarize(mean_total_branched_odd_new = mean(mean_total_branched_odd),
            sd_total_branched_odd_new = sd(mean_total_branched_odd)) %>%
  pivot_wider(names_from = "stt", values_from = c("mean_total_branched_odd_new",
                                                  "sd_total_branched_odd_new"))

branched_odd_tsidw <- inner_join(branched_odd_chain_fatty_acids %>%
                                   mutate(stt = tolower(stt)), tsidw_formatted) %>%
  mutate(analyte = "branched_odd_chains") %>%
  rename("mean_conc" = "mean_total_branched_odd") %>%
  select(tourist_season, stt, site, analyte, mean_tsidw_pop, mean_conc)

branched_odd_table <- branched_odd_chain_fatty_acids %>%
  mutate(stt = tolower(stt)) %>%
  ungroup() %>%
  group_by(tourist_season, stt) %>%
  summarize(mean_conc = mean(mean_total_branched_odd, na.rm = TRUE),
            sd_conc = sd(mean_total_branched_odd, na.rm = TRUE),
            cv_conc = sd_conc/mean_conc) %>%
  pivot_wider(names_from = "stt", 
              values_from = c("mean_conc", "sd_conc", "cv_conc"))

write.csv(x = branched_odd_table, 
          file = "../figures_tables/branched_odd_fa_summary_stats.csv", 
          row.names = FALSE)

## PPCPs

ppcp_reduced_formatted <- ppcp_reduced %>%
  group_by(tourist_season, stt) %>%
  summarize(mean_ppcp_conc = mean(mean_conc),
            sd_conc = sd(mean_conc),
            cv_conc = sd_conc/mean_ppcp_conc) %>%
  rename("mean_conc" = "mean_ppcp_conc") %>%
  pivot_wider(names_from = "stt", 
              values_from = c("mean_conc", "sd_conc", "cv_conc"))

write.csv(x = ppcp_reduced_formatted, 
          file = "../figures_tables/ppcp_summary_stats.csv", 
          row.names = FALSE)

ppcp_tsidw <- left_join(ppcp_reduced, tsidw_formatted) %>%
  mutate(analyte = "ppcp") %>%
  select(tourist_season, stt, site, analyte, mean_tsidw_pop, mean_conc)

## Combine all the data frames

combined_tsidw_df <- rbind(nutrients_tsidw, afdm_tsidw, 
                           branched_odd_tsidw, ppcp_tsidw)

combined_tsidw_lm_log_mod <- combined_tsidw_df %>%
  group_by(analyte) %>%
  filter(analyte != "branched_odd_chains") %>%
  nest() %>%
  mutate(model = map(data, function(df) lm(log10(mean_conc) ~ log10(mean_tsidw_pop), 
                                           data = df)))

combined_tsidw_lm_asinsqrt_mod <- combined_tsidw_df %>%
  group_by(analyte) %>%
  filter(analyte == "branched_odd_chains") %>%
  nest() %>%
  mutate(model = map(data, function(df) lm(asin(sqrt(mean_conc)) ~ log10(mean_tsidw_pop), 
                                           data = df)))

combined_tsidw_mod <- rbind(combined_tsidw_lm_log_mod, 
                            combined_tsidw_lm_asinsqrt_mod) %>%
  mutate(tidy_model = map(.x = model,
                          .f = ~ glance(.x))) %>%
  select(-data, -model) %>%
    unnest(cols = c(tidy_model))

plots_combined <- map2(.x = unique(combined_tsidw_df$analyte),
                       .y = c("Ammonia (mg/L)", "Nitrate/Nitrite (mg/L)",
                             "SRP (mg/L)", "Total Nitrogen (mg/L)", 
                             "Total Phosphorus (mg/L)","AFDM (mg)",
                             "Branched + Odd Chain Fatty Acids",
                             "Total PPCP (ng/L)"),
     .f = ~ if(.x != "branched_odd_chain"){
     combined_tsidw_df %>%
      filter(analyte == .x) %>%
      ggplot(aes(x = (mean_tsidw_pop),
                 y = (mean_conc))) +
      geom_point(size = 6) +
      geom_smooth(method = "lm", se = TRUE) +
      scale_y_log10()+
      scale_x_log10() +
      ylab(paste(.y)) +
      xlab("TSIDW Population (people)") +
      ggtitle(paste(.y, " vs. TSIDW Population")) +
      geom_richtext(data = combined_tsidw_mod %>%
                      filter(analyte == .x), 
                    mapping = aes(label = paste0("p-value: ",
                                                 round(p.value, 3),
                                                 "<br>R<sup>2</sup>: ",
                                                 round(r.squared, 3),
                                                 "<br> N = ", nobs), 
                                  x = 6000, 
                                  y = max(combined_tsidw_df %>%
                                        ungroup() %>%
                                        filter(analyte == .x) %>%
                                        select(mean_conc))*0.7), size = 7) +
      theme_minimal() +
      theme(strip.text = element_markdown(),
            axis.text = element_text(size = 24),
            axis.title = element_text(size = 28),
            plot.title = element_text(size = 34))
     } else if (.x == "branched_odd_chain") {
       combined_tsidw_df %>%
         filter(analyte == .x) %>%
         ggplot(aes(x = (mean_tsidw_pop),
                    y = asin(sqrt(mean_conc)))) +
         geom_point(size = 6) +
         geom_smooth(method = "lm", se = TRUE) +
         #scale_y_log10()+
         scale_x_log10() +
         ylab("Branched + Odd Chain Fatty Acids (\u00b5g/mg)") +
         xlab("TSIDW Population (people)") +
         ggtitle("Branched + Odd Chain FA vs. TSIDW Population") +
         geom_richtext(data = combined_tsidw_mod %>%
                         filter(analyte == .x), 
                       mapping = aes(label = paste0("p-value: ",
                                                    round(p.value, 3),
                                                    "<br>R<sup>2</sup>: ",
                                                    round(r.squared, 3),
                                                    "<br> N = ", nobs), 
                                     x = 6000, 
                                     y = max(combined_tsidw_df %>%
                                               ungroup() %>%
                                               filter(analyte == .x) %>%
                                               select(mean_conc))*0.7), size = 7) +
         theme_minimal() +
         theme(strip.text = element_markdown(),
               axis.text = element_text(size = 24),
               axis.title = element_text(size = 28),
               plot.title = element_text(size = 34))
     }
       )

ggarrange(plotlist = plots_combined, 
          ncol = 2, nrow = 4,
          labels = "AUTO", 
          font.label = list(size = 36)) %>%
  ggexport(filename = "../figures_tables/combined_plot_continuous.png",
           height = 4000, width = 3200, res = 120)
