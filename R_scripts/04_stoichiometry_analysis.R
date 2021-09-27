## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script is a piece of the 
## overarching analysis and looks mostly at relating 
## periphyton Carbon, Nitrogen, and Phosphorus stoichiometries
## with wastewater treatment infrastructure
## and a sampling point's temporal position to the tourist 
## season. Questions about this script should be directed to 
## Michael F. Meyer (michael.f.meyer@wsu.edu). 


# 0. Create function to run permutational analysis ------------------------

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
library(janitor)

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

stoich_table <- stoich %>%
  mutate(carbon_nitrogen = carbon/nitrogen,
         carbon_phosphorus = carbon/phosphorus,
         nitrogen_phosphorus = nitrogen/phosphorus) %>%
  select(tourist_season, stt, carbon_nitrogen:nitrogen_phosphorus) %>%
  pivot_longer(cols = c(carbon_nitrogen:nitrogen_phosphorus), 
               names_to = "elements", values_to = "ratios") %>%
  group_by(tourist_season, stt, elements) %>%
  summarize(mean_ratio = mean(ratios),
            sd_ratio = sd(ratios),
            cv_ratio = sd_ratio/mean_ratio) %>%
  pivot_wider(names_from = "stt", 
              values_from = c("mean_ratio", "sd_ratio", "cv_ratio")) %>%
  select(elements, tourist_season, mean_ratio_centralized:cv_ratio_decentralized) %>%
  arrange(elements)
  
write.csv(x = stoich_table, 
          file = "../figures_tables/stoich_summary_stats.csv", 
          row.names = FALSE)
  
# 3. Create stoichiometric plots and analyses -----------------------------


stt_labels <- c("Centralized", "Decentralized")
names(stt_labels) <- c("centralized", "decentralized")

carbon_nitrogen <- ggplot(data = stoich, aes(tourist_season, carbon/nitrogen)) +
  geom_boxplot(stoich %>%
                 filter(!(stt == "centralized" & tourist_season == "Out of Season")), 
               mapping = aes(tourist_season, carbon/nitrogen),
               alpha = 0.33, outlier.alpha = 0) +
  geom_jitter(stoich, 
              mapping = aes(tourist_season, carbon/nitrogen), size = 3) +
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

carbon_phosphorus <- ggplot(data = stoich, aes(tourist_season, carbon/phosphorus)) +
  geom_boxplot(stoich %>%
                 filte`r(!(stt == "centralized" & tourist_season == "Out of Season")), 
               mapping = aes(tourist_season, carbon/phosphorus),
               alpha = 0.33, outlier.alpha = 0) +
  geom_jitter(stoich, 
              mapping = aes(tourist_season, carbon/phosphorus), size = 3) +
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

nitrogen_phosphorus <- ggplot(data = stoich, aes(tourist_season, nitrogen/phosphorus)) +
  geom_boxplot(stoich %>%
                 filter(!(stt == "centralized" & tourist_season == "Out of Season")), 
               mapping = aes(tourist_season, nitrogen/phosphorus),
               alpha = 0.33, outlier.alpha = 0) +
  geom_jitter(stoich, 
              mapping = aes(tourist_season, nitrogen/phosphorus), size = 3) +
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


cn_lm <- Anova(lm(carbon/nitrogen ~ stt * tourist_season, 
                  data = stoich), 
               type = "II")

cn_permutation <- permute_data_analytics(data = stoich %>%
                                           mutate(carbon_nitrogen = carbon/nitrogen), 
                                         metric = "carbon_nitrogen", 
                                         full_model = cn_lm, 
                                         metric_plot_title = "Carbon:Nitrogen ~ STT * Tourist Season", 
                                         transform_response = "none")
                                              
stt_bf_cn <- stoich %>% 
  group_by(stt) %>%
  emmeans_test(carbon/nitrogen ~ tourist_season, p.adjust.method = "bonferroni") 

ts_bf_cn <- stoich %>% 
  group_by(tourist_season) %>%
  emmeans_test(carbon/nitrogen ~ stt, p.adjust.method = "bonferroni") 


cp_lm <- Anova(lm(carbon/phosphorus ~ stt * tourist_season,
                  data = stoich), 
               type = "II")

cp_permutation <- permute_data_analytics(data = stoich %>%
                                           mutate(carbon_phosphorus = carbon/phosphorus), 
                                         metric = "carbon_phosphorus", 
                                         full_model = cp_lm, 
                                         metric_plot_title = "Carbon:Phosphorus ~ STT * Tourist Season", 
                                         transform_response = "none")

stt_bf_cp <- stoich %>% 
  group_by(stt) %>%
  emmeans_test(carbon/phosphorus ~ tourist_season, p.adjust.method = "bonferroni") 

ts_bf_cp <- stoich %>% 
  group_by(tourist_season) %>%
  emmeans_test(carbon/phosphorus ~ stt, p.adjust.method = "bonferroni") 

np_lm <- Anova(lm(nitrogen/phosphorus ~ stt * tourist_season, 
                  data = stoich), 
               type = "II")

np_permutation <- permute_data_analytics(data = stoich %>%
                                           mutate(nitrogen_phosphorus = nitrogen/phosphorus), 
                                         metric = "nitrogen_phosphorus", 
                                         full_model = np_lm, 
                                         metric_plot_title = "Nitrogen:Phosphorus ~ STT * Tourist Season", 
                                         transform_response = "none")

stt_bf_np <- stoich %>% 
  group_by(stt) %>%
  emmeans_test(nitrogen/phosphorus ~ tourist_season, p.adjust.method = "bonferroni") 

ts_bf_np <- stoich %>% 
  group_by(tourist_season) %>%
  emmeans_test(nitrogen/phosphorus ~ stt, p.adjust.method = "bonferroni") 


# 4. Save Figures ------------------------------------------------------

permutation_plots <- c(cn_permutation, cp_permutation, np_permutation)

permutation_plot_names <- c("cn_permutation", "cp_permutation", "np_permutation")

walk2(.x = permutation_plots,
      .y = permutation_plot_names,
      .f = ~ ggsave(filename = paste(.y, "_histogram.png", sep = ""), 
                    plot = .x, 
                    device = "png", path = "../figures_tables", 
                    width = 16, height = 8, units = "in"))

anova_table <- list(clean_names(data.frame(cn_lm)), 
                    clean_names(data.frame(cp_lm)), 
                    clean_names(data.frame(np_lm)))

anova_table_names <- c("cn", "cp", "np")

walk2(.x = anova_table,
      .y = anova_table_names,
      .f = ~ write.csv(file = paste("../figures_tables/", .y, "_anova_table.csv", sep = ""), 
                       x = .x, row.names = TRUE))

stt_bf <- rbind(stt_bf_cn, stt_bf_cp, stt_bf_np) %>%
  rename("variable" = ".y.")

write.csv(x = stt_bf, 
          file = "../figures_tables/stt_bf_stoich_results.csv", 
          row.names = FALSE)

ts_bf <- rbind(ts_bf_cn, ts_bf_cp, ts_bf_np) %>%
  rename("variable" = ".y.")

write.csv(x = ts_bf, 
          file = "../figures_tables/ts_bf_stoich_results.csv", 
          row.names = FALSE)
