library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(ggstar)
library(cowplot)
library(ggrepel)
library(viridis)

communities <- st_read("../raw_data/Flathead sampling locations.kml")

tsidw_pop <- read.csv("../cleaned_data/averaged_temporally_scaled_inverse_distance_weighted_population_metrics.csv")

ppcp <- read.csv("../cleaned_data/ppcp.csv") %>%
  filter(month != "may") %>%
  group_by(tourist_season, stt, site, peri_sampling, sampling_event, ppcp) %>%
  summarize(mean_conc = mean(concentration, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(tourist_season, stt, site, peri_sampling, sampling_event) %>%
  summarize(total_concentration = sum(mean_conc, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(month = ifelse(peri_sampling == "peri_1", "june", NA),
         month = ifelse(peri_sampling == "peri_2", "july", month),
         month = ifelse(peri_sampling == "peri_3", "august", month),
         month = ifelse(peri_sampling == "peri_4", "september", month)) %>%
  filter(!is.na(month)) %>%
  group_by(tourist_season, stt, site, month) %>%
  summarize(mean_conc = mean(total_concentration, na.rm = TRUE)) 

sample_points_sf <- communities %>%
  filter(Name %in% c("BB", "YB", "FLBS", "BD", "WB", "WF", "HO",
                     "LK", "WS", "DA", "DU", "FI", "SJ", "SL", "BO")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  inner_join(., tsidw_pop, by = c("Name" = "site")) %>%
  mutate(month = factor(month, levels = c("may", "june", "july", "august", "september"), 
                        labels = c("May", "June", "July", "August", "September")))

# sample_points_sf <- communities %>%
#   filter(Name %in% c("BB", "YB", "FLBS", "BD", "WB", "WF", "HO",
#                      "LK", "WS", "DA", "DU", "FI", "SJ", "SL", "BO")) %>%
#   st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
#   inner_join(., ppcp, by = c("Name" = "site")) %>%
#   mutate(month = factor(month, levels = c("june", "july", "august", "september"), 
#                         labels = c("June", "July", "August", "September"))) 


sample_points <- sample_points_sf %>%
  st_coordinates() %>%
  data.frame() %>%
  cbind(., as.character(sample_points_sf$Name)) %>%
  rename("site" = "as.character(sample_points_sf$Name)",
         "lat" = "Y",
         "lon" = "X")

state_parks <- communities %>%
  filter(grepl(pattern = "State_Park", x = Name))

developments <- communities %>%
  filter(Name %in% c("Bigfork_shapefile", "Polson_shapefile", "Dayton_shapefile",
                     "Woods_Bay_shapefile", "Lakeside_shapefile"))

main_map <- ggplot() +
  annotation_map_tile(
    type = 'https://stamen-tiles.a.ssl.fastly.net/terrain-background/${z}/${x}/${y}.png',
    zoom = 11,
    cachedir = "data/map_tiles/") +
  geom_sf(data = developments, color = "black", alpha = 0.5, fill = viridis(10)[2]) +
  geom_sf(data = sample_points_sf, color = "black", fill = viridis(20)[10], pch = 21,
          size = 3, alpha = 0.7) +
  geom_sf(data = state_parks, color = "black", fill = viridis(20)[18], pch = 23,
          size = 3, alpha = 0.7) +
  geom_text_repel(data = sample_points,
                  aes(x = lon, y = lat, label = site),
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  nudge_x = c(.05, .07, .06, .05, .08, .07, 0, -.08, -.05, -.05,
                              -.05, -.05, -.02, .02, .02),
                  nudge_y = c(-.02, -.01, .02, -.01, 0, 0, .03, 0, 0, 0, 0, .03,
                              -.04, .03, -.03)) +
  annotation_north_arrow(which_north = "true", location = "tl",
                         width = unit(1.1, "cm"), height = unit(1.1, "cm")) +
  annotation_scale(location = "bl") +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = c(-114.45, -114.25, -114.05, -113.85)) +
  scale_y_continuous(breaks = c(47.70, 47.82, 47.94, 48.06)) +
  coord_sf(xlim = c(-114.5268, -113.7768), ylim = c(47.65, 48.15), crs = 4326) +
  theme_bw()

# Set up polygon to show boundaries of sampling map
bounds <- data.frame(
  x = c(-114.5268, -113.7768, -113.7768, -114.5268),
  y = c(48.15, 48.15, 47.65, 47.65)) %>%
  summarize(x = mean(x),
            y = mean(y))

inset_map <- ne_states(country = 'United States of America', returnclass = "sf") %>%
  filter(!(name %in% c("Alaska", "Hawaii"))) %>%
  ggplot() +
  geom_sf(data = ne_countries(continent = "North America", returnclass = "sf"),
          fill = "gray70", size = 0.25) +
  geom_sf(fill = "gray99") +
  geom_star(data = bounds, aes(x = x, y = y),
            color = "black", fill = "#d24644", size = 3) +
  theme_bw() +
  coord_sf(xlim = c(-130, -60), ylim = c(22, 50)) +
  xlab("") +
  ylab("") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank())

flathead_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(plot = inset_map,
            x = 0.065,
            y = 0.105, 
            width = .45, 
            height = .27, 
            scale = .66) 

ggsave(filename = "../figures_tables/flathead_map.png", plot = flathead_map,
       width = 5, height = 5, units = "in", device = "png")



############ TESTING



map_month <- ggplot() +
  annotation_map_tile(
    type = 'https://stamen-tiles.a.ssl.fastly.net/terrain-background/${z}/${x}/${y}.png',
    zoom = 11,
    cachedir = "data/map_tiles/") +
  #geom_sf(data = developments, color = "black", alpha = 0.5, fill = viridis(10)[2]) +
  geom_sf(data = sample_points_sf, 
          aes(fill = (atsidw_pop),
              size = (atsidw_pop)), 
          color = "black", alpha = 0.7, pch = 21) +
  scale_fill_viridis_c(option = "plasma", name = "TSIDW Pop", 
                       trans = "log",
                       limits = c(200, 27000),
                       breaks = c(300, 1000, 3000, 10000, 25000), 
                       labels = c(300, 1000, 3000, 10000, 25000)) +
  scale_size_continuous(range = c(5,12), guide = NULL, trans = "log") +
  facet_wrap(~ month) +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = c(-114.45, -114.25, -114.05, -113.85)) +
  scale_y_continuous(breaks = c(47.70, 47.82, 47.94, 48.06)) +
  coord_sf(xlim = c(-114.5268, -113.7768), ylim = c(47.65, 48.15), crs = 4326) +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 24),
        strip.text = element_text(size = 24),
        legend.key.height = unit(x = 0.75, units = "in"),
        legend.key.width = unit(x = 0.5, units = "in"))
  
ggsave(filename = "../figures_tables/flathead_tsidw_map_monthly.png", plot = map_month,
       width = 14, height = 10, units = "in", device = "png")







