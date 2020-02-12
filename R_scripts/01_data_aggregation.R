library(dplyr)
library(tidyr)
library(ggplot2)
library(xlsx)

file_names <- list.files(path = "../raw_data/flathead_coop_data/")

for( i in 1:length(file_names)){
  temp <- read.xlsx(paste("../raw_data/flathead_coop_data/", file_names[i], sep = ""), sheetIndex = 1)
  temp <- temp %>%
    select(May2017_Usage, Jun2017_Usage, Jul2017_Usage, Aug2017_Usage, Sep2017_Usage) %>%
    mutate(Pixel = paste(file_names[i]))
  if(i != 1){
    df_whole <- rbind(df_whole, temp)
  } else {
    df_whole <- temp
  }
}

utilities <- df_whole %>%
  mutate(Pixel = gsub(".xlsx", "", Pixel),
         Site = ifelse(Pixel %in% c("262007", "262018"), "Lakeside", NA),
         Site = ifelse(Pixel %in% c("272026", "272025", "272036", "272035"), "Big_Fork", Site),
         Site = ifelse(Pixel %in% c("261919", "261920", "261918"), "Woods_Bay", Site)) %>%
  gather(Time, Usage, May2017_Usage:Sep2017_Usage) %>%
  mutate(Time = gsub("2017_Usage", "", Time)) %>%
  filter(Usage != "MT") %>%
  filter(!is.na(Usage)) %>%
  group_by(Site, Time) %>%
  mutate(Usage = as.numeric(Usage)) %>%
  summarize(Cum_usage = sum(Usage))

utilities$Time <- factor(utilities$Time, levels = c("May", "Jun", "Jul", "Aug", "Sep"))
ggplot(utilities %>% filter (Time != "May"), aes(Time, log10(Cum_usage))) +
  geom_point(size = 3) +
  facet_wrap(~ Site)

write.csv(x = utilities, "../cleaned_data/utilities.csv", row.names = FALSE)
