library(tidyverse)
library(lubridate)
library(stringr)
library(janitor)
library(sf)
library(spdplyr)
library(OpenStreetMap)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(viridis)

# Load the shapefiles and metadata containing information of the development
# polygons and the lake shorelines.

flathead_shapefile <- sf::st_read(dsn = "../raw_data/Flathead sampling locations.kml")

sampling_loc <- c("LK", "WS", "DA", "SJ", "SL", "BO", "DU", 
                  "FI", "BB", "FLBS", "YB", "BD", "WB", "WF", "HO")

# Extract and format development area data 
loc_areas <- flathead_shapefile %>%
  filter(grepl("shapefile", Name)) %>%
  st_area() %>%
  enframe(name = NULL) %>%
  rename(development_area_m2 = value) %>%
  mutate(development_area_m2 = as.numeric(development_area_m2),
         development_area_km2 = development_area_m2 / 1000000) %>%
  cbind(flathead_shapefile %>% as_tibble() %>%
          filter(grepl("shapefile", Name)) %>%
          mutate(Name = gsub(pattern = "_shapefile", 
                             replacement = "", x = Name),
                 Name = gsub(pattern = "_", 
                             replacement = " ", x = Name)) %>%
          dplyr::select(Name) %>%
          as.vector()) %>%
  rename("development" = Name)

loc_areas_centroids <- flathead_shapefile %>% 
  filter(grepl("shapefile", Name)) %>%
  st_centroid() %>%
  as_tibble() %>%
  mutate(Name = gsub(pattern = "_shapefile", 
                     replacement = "", x = Name),
         Name = gsub(pattern = "_", 
                     replacement = " ", x = Name)) %>%
  dplyr::select(development = Name, centroid_loc = geometry) %>%
  inner_join(. , y = loc_areas)

# Extract shoreline length info and join with area
loc_shoreline_area_length <- flathead_shapefile %>%
  filter(grepl("shoreline", Name)) %>%
  st_length() %>%
  as_tibble() %>%
  rename(development_shoreline_length_m = value) %>%
  mutate(development_shoreline_length_m = as.numeric(development_shoreline_length_m),
         development_shoreline_length_km = development_shoreline_length_m / 1000) %>%
  cbind(flathead_shapefile %>% as_tibble() %>%
          filter(grepl("shoreline", Name)) %>%
          mutate(Name = gsub(pattern = "_shoreline", 
                             replacement = "", x = Name),
                 Name = gsub(pattern = "_", 
                             replacement = " ", x = Name)) %>%
          dplyr::select(Name) %>%
          as.vector()) %>%
  rename(development = Name) %>%
  full_join(x = ., y = loc_areas_centroids, by = c("development")) %>%
  dplyr::select(development, centroid_loc, development_shoreline_length_km, development_area_km2) %>%
  mutate(neighbor = gsub(pattern = " ", replacement = "_", development),
         neighbor = tolower(neighbor))

# Convert metadata into a spatial object 
site_loc <- flathead_shapefile %>%
  filter(Name %in% sampling_loc) 


# Join all data together and format into a dataframe
locs_centroids <- st_distance(x = site_loc,
                              y = loc_areas_centroids$centroid_loc) %>%
  as_tibble() %>%
  rename("Bigfork" = V1,
         "Woods Bay" = V2,
         "Bear Dance" = V3,
         "FLBS" = V4,
         "Yellow Bay" = V5,
         "Blue Bay" = V6,
         "Finley" = V7,
         "Polson" = V8,
         "Dayton" = V9,
         "West Shore" = V10,
         "Lakeside" = V11) %>%
  clean_names() %>%
  cbind(., site_loc$Name) %>%
  as_tibble() %>%
  gather(key = neighbor, value = distance, bigfork:lakeside) %>%
  rename("sampling_site" = `site_loc$Name`) %>%
  mutate(distance_km = distance / 1000,
         population = ifelse(test = neighbor == "bigfork",
                             yes = 4797, no = NA),
         population = ifelse(test = neighbor == "woods_bay",
                             yes = 688, no = population),
         population = ifelse(test = neighbor == "bear_dance",
                             yes = 212, no = population),
         population = ifelse(test = neighbor == "flbs",
                             yes = 90, no = population),
         population = ifelse(test = neighbor == "yellow_bay",
                             yes = 90, no = population),
         population = ifelse(test = neighbor == "blue_bay",
                             yes = 110, no = population),
         population = ifelse(test = neighbor == "finley",
                             yes = 453, no = population),
         population = ifelse(test = neighbor == "polson",
                             yes = 4702, no = population),
         population = ifelse(test = neighbor == "dayton",
                             yes = 112, no = population),
         population = ifelse(test = neighbor == "west_shore",
                             yes = 66, no = population),
         population = ifelse(test = neighbor == "lakeside",
                             yes = 2459, no = population)) %>%
  filter(!(neighbor %in% c("bear_dance", "flbs", "blue_bay", "finley", "west_shore", "yellow_bay"))) %>%
  left_join(x = ., y = loc_shoreline_area_length) %>%
  mutate(distance_weighted_population = ((population * development_shoreline_length_km) /
                                           development_area_km2) / distance_km) %>%
  group_by(sampling_site) %>%
  summarize(distance_weighted_population = (sum(as.numeric(distance_weighted_population)))) %>%
  arrange(distance_weighted_population)

write.csv(x = locs_centroids,
          file = "../cleaned_data/static_distance_weighted_population_metrics.csv",
          row.names = FALSE)

### Calculate temporal scaling factors

### Start with the park admittance data

park_data_orig <- read.csv("../raw_data/state_park_attendance.csv")

## We will standardize the scaling metric by May attendance, because this is likely 
## before the tourism season started 

park_scalars <- park_data_orig %>%
  pivot_longer(cols = big_arm:yellow_bay, names_to = "state_park", 
               values_to = "admitted_persons") %>%
  pivot_wider(names_from = Month, values_from = admitted_persons) %>%
  mutate(june_scalar = June/May,
         july_scalar = July/May,
         august_scalar = August/May,
         september_scalar = September/May,
         september_scalar = ifelse(state_park == "finley_point", 
                                   (September+August)/2/May, 
                                   september_scalar)) %>%
  dplyr::select(state_park, june_scalar:september_scalar)

flathead_parks_locs <- flathead_shapefile %>%
  filter(grepl(pattern = "State_Park", x = Name)) %>%
  mutate(state_park = gsub(pattern = "_State_Park", replacement = "", x = Name),
         state_park = tolower(state_park)) %>%
  full_join(., y = park_scalars)

flathead_park_sampling_locs <- st_distance(x = site_loc,
                                           y = flathead_parks_locs) %>%
  as_tibble() %>%
  rename("west_shore" = V1,
         "big_arm" = V2,
         "finley_point" = V3,
         "wayfarers" = V4,
         "yellow_bay" = V5) %>%
  clean_names() %>%
  cbind(., site_loc$Name) %>%
  as_tibble() %>%
  rename(sampling_site = `site_loc$Name`) %>%
  pivot_longer(cols = west_shore:yellow_bay, names_to = "state_park", values_to = "distance") %>%
  inner_join(. , y = park_scalars) %>%
  pivot_longer(cols = june_scalar:september_scalar, names_to = "month_scalar", values_to = "scalar_value") %>%
  mutate(weighted_scalar = scalar_value/(distance/1000)) %>%
  group_by(sampling_site, month_scalar) %>%
  summarize(inverse_distance_weighted_scalar = log10(sum(as.numeric(weighted_scalar))+1)) %>%
  arrange(((inverse_distance_weighted_scalar)))

## Combined space and time component

flathead_park_space_time <- inner_join(x = locs_centroids,
                                       y = flathead_park_sampling_locs) %>%
  mutate(scaled_idw_population = ((distance_weighted_population) * inverse_distance_weighted_scalar)) %>%
  arrange(((scaled_idw_population)))


write.csv(x = flathead_park_space_time,
          file = "../cleaned_data/temporal_inverse_distance_weighted_population_metrics.csv",
          row.names = FALSE)


## Make a map

combioned_shp <- full_join(x = flathead_shapefile, 
                           y = flathead_park_space_time, 
                           by = c("Name" = "sampling_site"))


base_map <- openmap(upperLeft = c(48.987618, -116.054273),
                    lowerRight = c(45.029625, -104.077145),
                    type = "esri") %>%
  openproj()


base_map_zoom <- openmap(upperLeft = c(48.119469, -114.398000),
                         lowerRight = c(47.663195, -113.947497),
                         type = "bing", zoom = 11) %>%
  openproj()

inset_map <- autoplot.OpenStreetMap(base_map) +
  geom_point(data = site_loc, 
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y)),
             stat = "sf_coordinates", 
             alpha = 1,  color = "grey20", 
             size = 3, shape = 21) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.background = element_rect(fill = "snow1")) +
  #annotate(geom = "text", label = "Irkutsk Oblast",
  #         x = 105, y = 54, color = "white", size = 5) +
  #annotate(geom = "text", label = "Republic of\nBuryatiya",
  #         x = 109, y = 52.35, color = "white", size = 5) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

zoom_map <- autoplot.OpenStreetMap(base_map_zoom) +
  geom_point(data = combioned_shp %>%
               filter(Name %in% c("BB", "BD", "FLBS", "YB", "SL",
                                  "SJ", "BO", "DU", "WS", "LK", "DA",
                                  "WF", "FI", "WB")) %>%
               mutate(month_scalar = factor(month_scalar, 
                                            levels = c("june_scalar", "july_scalar", "august_scalar", "september_scalar"), 
                                            labels = c("June", "July", "August", "September"))), 
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y),
                 size = log10(scaled_idw_population),
                 fill = log10(scaled_idw_population)),
             stat = "sf_coordinates", 
             alpha = 0.8,  color = "grey70", shape = 21,
             stroke = 2.5) + 
  scale_fill_viridis(option = "plasma", name = "log10(IDW Pop)") +
  scale_size_continuous(range = c(1, 12), guide = "none") +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~month_scalar) +
  theme(legend.key.height = unit(1, "in"),
        legend.key.width = unit(0.65, "in"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        panel.background = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 20))

sites_ts <- ggplot(data = flathead_park_space_time %>%
                     mutate(month_scalar = factor(month_scalar, 
                                                  levels = c("june_scalar", "july_scalar", "august_scalar", "september_scalar"), 
                                                  labels = c("June", "July", "August", "September"))), 
                   aes(month_scalar, (scaled_idw_population)))+
  geom_point(size = 5) +
  facet_wrap(~sampling_site)

## Alternative method for creating an aggregate scalar across parks

park_scalars <- park_data_orig %>%
  pivot_longer(cols = big_arm:yellow_bay, names_to = "state_park", 
               values_to = "admitted_persons") %>%
  pivot_wider(names_from = Month, values_from = admitted_persons) %>%
  mutate(may_scalar = May/May, 
         june_scalar = June/May,
         july_scalar = July/May,
         august_scalar = August/May,
         september_scalar = September/May,
         september_scalar = ifelse(state_park == "finley_point", 
                                   (September+August)/2/May, 
                                   september_scalar)) %>%
  pivot_longer(cols = may_scalar:september_scalar, 
               names_to = "month", 
               values_to = "scalar_value") %>%
  group_by(month) %>%
  summarize(average_scalar = mean(scalar_value)) %>%
  pivot_wider(names_from = "month", values_from = "average_scalar")

park_scalars %>%
  pivot_longer(cols = c(august_scalar:september_scalar), names_to = "month", values_to = "scalar") %>%
  mutate(month = factor(x = month, 
                        levels = c("june_scalar", "july_scalar", "august_scalar", "september_scalar"),
                        labels = c("June", "July", "August", "September"))) %>%
  filter(!is.na(month)) %>%
  ggplot() +
  geom_point(aes(x = month, y = scalar), size = 5) +
  geom_hline(aes(yintercept = 1.0), linetype = "dashed", size = 2) +
  ggtitle("Mean Aggregate Scalar") +
  ylab("Mean Aggregate Scalar") +
  xlab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size = 24))

locs_centroids_scaled <- locs_centroids %>%
  mutate(may_idw_pop = distance_weighted_population * park_scalars$may_scalar,
         june_idw_pop = distance_weighted_population * park_scalars$june_scalar, 
         july_idw_pop = distance_weighted_population * park_scalars$july_scalar,
         august_idw_pop = distance_weighted_population * park_scalars$august_scalar,
         september_idw_pop = distance_weighted_population * park_scalars$september_scalar) %>%
  select(-distance_weighted_population) %>%
  pivot_longer(cols = c(may_idw_pop:september_idw_pop), names_to = "month", values_to = "idw_pop") %>%
  mutate(month = gsub("_idw_pop", "", month),
         month = str_to_sentence(month))

## Make a map

combioned_shp <- full_join(x = flathead_shapefile, 
                           y = locs_centroids_scaled, 
                           by = c("Name" = "sampling_site"))


base_map <- openmap(upperLeft = c(48.987618, -116.054273),
                    lowerRight = c(45.029625, -104.077145),
                    type = "esri") %>%
  openproj()


base_map_zoom <- openmap(upperLeft = c(48.119469, -114.398000),
                         lowerRight = c(47.663195, -113.947497),
                         type = "bing", zoom = 11) %>%
  openproj()

inset_map <- autoplot.OpenStreetMap(base_map) +
  geom_point(data = site_loc, 
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y)),
             stat = "sf_coordinates", 
             alpha = 1,  color = "grey20", 
             size = 3, shape = 21) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.background = element_rect(fill = "snow1")) +
  #annotate(geom = "text", label = "Irkutsk Oblast",
  #         x = 105, y = 54, color = "white", size = 5) +
  #annotate(geom = "text", label = "Republic of\nBuryatiya",
  #         x = 109, y = 52.35, color = "white", size = 5) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

zoom_map <- autoplot.OpenStreetMap(base_map_zoom) +
  geom_point(data = combioned_shp %>%
               filter(Name %in% c("BB", "BD", "FLBS", "YB", "SL",
                                  "SJ", "BO", "DU", "WS", "LK", "DA",
                                  "WF", "FI", "WB")) %>%
               mutate(month = factor(month,levels = c("May", "June", "July", "August", "September"), 
                                           labels = c("May", "June", "July", "August", "September"))), 
             aes(geometry = geometry,
                 x = after_stat(x),
                 y = after_stat(y),
                 size = log10(idw_pop),
                 fill = log10(idw_pop)),
             stat = "sf_coordinates", 
             alpha = 0.8,  color = "grey70", shape = 21,
             stroke = 2.5) + 
  scale_fill_viridis(option = "plasma", name = "log10(IDW Pop)") +
  scale_size_continuous(range = c(1, 12), guide = "none") +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~month) +
  theme(legend.key.height = unit(1, "in"),
        legend.key.width = unit(0.65, "in"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        panel.background = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 20))


#### Boxplot to complemenet other analyses

locs_centroids_scaled <- locs_centroids_scaled %>%
  filter(month != "May") %>%
  mutate(tourist_season = ifelse(month %in% c("June", "July", "August"), "In Season", "Out of Season"),
         stt = ifelse(sampling_site %in% c("WF", "HO", "FLBS", "YB"), "Centralized", "Decentralized"))

ggplot(locs_centroids_scaled, aes(tourist_season, log10(idw_pop))) +
  geom_boxplot(alpha = 0.33, outlier.alpha = 0, width = 0.33) +
  geom_jitter() +
  xlab("") +
  ylab("log10(IDW Population)") +
  facet_wrap(~stt) +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        strip.text = element_text(size = 18))


##### Maps for Scalars

park_scalars <- park_data_orig %>%
  pivot_longer(cols = big_arm:yellow_bay, names_to = "state_park", 
               values_to = "admitted_persons") %>%
  pivot_wider(names_from = Month, values_from = admitted_persons) %>%
  mutate(june_scalar = June/May,
         july_scalar = July/May,
         august_scalar = August/May,
         september_scalar = September/May,
         september_scalar = ifelse(state_park == "finley_point", 
                                   (September+August)/2/May, 
                                   september_scalar)) %>%
  dplyr::select(state_park, june_scalar:september_scalar)

flathead_parks_locs <- flathead_shapefile %>%
  filter(grepl(pattern = "State_Park", x = Name)) %>%
  mutate(state_park = gsub(pattern = "_State_Park", replacement = "", x = Name),
         state_park = tolower(state_park)) %>%
  cbind(., st_coordinates(flathead_shapefile %>%
                            filter(grepl(pattern = "State_Park", x = Name)))) %>%
  st_drop_geometry() %>%
  full_join(., y = park_scalars) %>%
  pivot_longer(cols = c(june_scalar:september_scalar), names_to = "month", values_to = "scalar") %>%
  mutate(month = factor(x = month, 
                        levels = c("june_scalar", "july_scalar", "august_scalar", "september_scalar"), 
                        labels = c("June", "July", "August", "September")))

zoom_map <- autoplot.OpenStreetMap(base_map_zoom) +
  geom_point(data = flathead_parks_locs, 
             aes(x = X,
                 y = Y,
                 size = (scalar),
                 fill = (scalar)),
             alpha = 0.8,  color = "grey70", shape = 21,
             stroke = 2.5) + 
  scale_fill_viridis(option = "plasma", name = "IDW Pop Scalar") +
  scale_size_continuous(range = c(1, 12), guide = "none") +
  xlab("Longitude") +
  ylab("Latitude") +
  facet_wrap(~month) +
  theme(legend.key.height = unit(1, "in"),
        legend.key.width = unit(0.65, "in"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24),
        panel.background = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 20))

