str(mcmcOut)

mcmc_df <- data.frame(mcmcOut)


test_plot <- mcmc_df %>%
  gather(parameter_est, value, sigma.e:sigma.site.month) %>%
  mutate(parameter_est = as.factor(parameter_est),
         sigma_label = ifelse(parameter_est == "sigma.month", "temporal", NA),
         sigma_label = ifelse(parameter_est == "sigma.site", "spatial", sigma_label),
         sigma_label = ifelse(parameter_est == "sigma.site.month", "spatio-temporal", sigma_label),
         sigma_label = ifelse(parameter_est == "sigma.e", "residual", sigma_label)) %>%
  ggplot(aes(sigma_label, value)) +
  geom_boxplot() +
  ylab("Sigma")

ggsave("../figures_tables/test_plot.png", test_plot)


test_df <- mcmc_df %>%
  gather(parameter_est, value, sigma.site:sigma.site.month) %>%
  mutate(parameter_est = as.factor(parameter_est))

aov_results <- aov(value ~ parameter_est, data = test_df)
anova(aov_results)
