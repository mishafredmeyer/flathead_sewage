library(tidyverse)
library(lubridate)
library(data.table)
library(ggpubr)

in_situ_temperature <- read_csv("../cleaned_data/flathead_temperature.csv") %>%
  mutate(date = ymd(date),
         site = ifelse(site == "SL", "SA", site))

landsat_temperature <- read_csv("../cleaned_data/Water_Temp_Extraction_L7_L8.csv") %>%
  filter(date > as.Date('2017-05-01') & date < as.Date('2017-10-01'),
         satellite == "LANDSAT_8") %>%
  mutate(site = case_when(id == 0 ~ "WB",
                          id == 1 ~ "BD", 
                          id == 2 ~ "FLBS",
                          id == 3 ~ "YB",
                          id == 4 ~ "BB",
                          id == 5 ~ "FI",
                          id == 6 ~ "DU",
                          id == 7 ~ "BO",
                          id == 8 ~ "SA",
                          id == 9 ~ "SJ",
                          id == 10 ~ "DA",
                          id == 11 ~ "WS",
                          id == 12 ~ "LK",
                          id == 13 ~ "HO",
                          id == 14 ~ "WF"),
         month = month(ymd(date))) 

# Convert both data.frames to data.tables
setDT(in_situ_temperature)
setDT(landsat_temperature)

# Duplicate the time column in df2 so we can verify the match later
landsat_temperature[, matched_date := date]

# Perform the rolling join
# Syntax: df2[df1, on = .(join_columns), roll = "nearest"]
result <- landsat_temperature[in_situ_temperature, on = .(site, date), roll = "nearest"]

prior_correction_plot <- ggplot(result) +
  geom_point(aes(y = surface_temp, x = mean)) +
  geom_smooth(aes(y = surface_temp, x = mean), method = "lm") +
  geom_abline(slope = 1, yintercept = 0) +
  ylab("In Situ Surface Temperature (\u00B0C)") +
  xlab("Remotely Sensed Surface Temperature (\u00B0C)") +
  ggtitle("Temperature Data Prior to Correcrtion") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 16))

ggsave("../figures_tables/temperature_bias_plot.png", plot = prior_correction_plot,
       width = 6, height = 6, units = "in")

temp_correction <- lm(surface_temp ~ mean, data = result)

result$predicted_temp <- predict(object = temp_correction, newdata = result)

post_correction <- ggplot(result) +
  geom_point(aes(y = surface_temp, x = predicted_temp)) +
  geom_smooth(aes(y = surface_temp, x = predicted_temp), method = "lm") +
  geom_abline(slope = 1, yintercept = 0) +
  ylab("In Situ Surface Temperature (\u00B0C)") +
  xlab("Remotely Sensed Surface Temperature (\u00B0C)") +
  ggtitle("Temperature Data After Correcrtion") +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 16))

ggsave("../figures_tables/temperature_corrected_plot.png", plot = post_correction,
       width = 6, height = 6, units = "in")

ggarrange(plotlist = c(prior_correction_plot, post_correction + theme(axis.title.y = element_blank())), 
          ncol = 2, nrow = 1, align = "hv",
          labels = "AUTO") %>%
  ggexport(filename = "../figures_tables/combined_plot_temp_correction.png", height = 600, width = 1250, res = 120)

landsat_temperature$predicted_temp <- predict(object = temp_correction, newdata = landsat_temperature)

target_dates <- data.table(
  date = as.Date(c("2017-06-28", 
                   "2017-06-29",
                   "2017-07-24", 
                   "2017-07-25",
                   "2017-07-26",
                   "2017-08-23",
                   "2017-08-24",
                   "2017-09-21"))
)

setDT(target_dates)

landsat_temperature %>%
  mutate(day = day(date)) %>%
  filter(month %in% c(5:9)) %>%
  # group_by(site, month) %>%
  # summarize(mean_temp = mean(predicted_temp)) %>%
  ggplot() +
  geom_point(aes(x = date, y = predicted_temp)) +
  facet_wrap(~site) + 
  ylab("Remotely Sensed Surface Temperature (\u00B0C)") +
  xlab("Date") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

ggsave("../figures_tables/site_level_temperature.png", height = 8, width = 14, units = "in")


write_csv(x = landsat_temperature, "../cleaned_data/temperature.csv")

landsat_temperature %>%
  mutate(day = day(date)) %>%
  filter(month %in% c(6:9),
         day >= 20,
         !(site %in% c("DU", "HO"))) %>%
  mutate(season = ifelse(month %in% c(6:8), "In Season", "Out of Season"),
         stt = ifelse(site %in% c("FLBS", "YB", "WF"), "Centralized", "Decentralized")) %>%
  ggplot() +
  geom_boxplot(aes(x = season, y = predicted_temp)) +
  facet_wrap(~stt)

in_situ_temperature %>%
  mutate(day = day(date)) %>%
  filter(!(site %in% c("DU", "HO"))) %>%
  mutate(stt = ifelse(site %in% c("FLBS", "YB", "WF"), "Centralized", "Decentralized")) %>%
  rowwise() %>%
  mutate(row_mean = mean(c(surface_temp, mid_temp, bottom_temp))) %>%
  ungroup() %>%
  ggplot() +
  geom_boxplot(aes(x = stt, y = row_mean)) 
