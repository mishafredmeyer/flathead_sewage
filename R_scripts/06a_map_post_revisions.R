## This script is associated with an analysis of spatially and 
## and temporally heterogeneous nutrient loading into Flathead
## Lake's nearshore zone. This script aggregates data necessary to make
## maps that detail the human population around Flathead Lake over space and 
## time. The scripts also make maps that highlight where samples were collected,
## and where human lakeside developments are located. 
##Questions about this script should be directed to 
## Michael F. Meyer (mfmeyer2@wisc.edu). 

## This script is divided into the following sections
## 1. Load Packages
## 2. Load and aggregate data
## 3. Create map of sampling locations by sewage treatment method
## 4. Create map of lakeside developments
## 5. Create map of littoral PPCP concentrations

# 1. Load Packages --------------------------------------------------------

library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(ggstar)
library(cowplot)
library(ggrepel)
library(viridis)


# 2. Load and aggregate data ----------------------------------------------

communities <- st_read("../raw_data/Flathead sampling locations.kml")

tsidw_pop <- read.csv("../cleaned_data/averaged_temporally_scaled_inverse_distance_weighted_population_metrics.csv")

ppcp <- read.csv("../cleaned_data/ppcp.csv") %>%
  #filter(month != "may") %>%
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


# 3. Create map of sampling locations by sewage treatment method ----------

site_map <- ggplot() +
  annotation_map_tile(
    type = 'https://stamen-tiles.a.ssl.fastly.net/terrain-background/${z}/${x}/${y}.png',
    zoom = 11,
    cachedir = "data/map_tiles/") +
  ##geom_sf(data = developments, color = "black", alpha = 0.5, fill = viridis(10)[2]) +
  geom_sf(data = sample_points_sf %>%
            mutate(stt = ifelse(Name %in% c("HO", "DU"), "none", stt)) %>%
            st_zm(drop = TRUE) %>%
            dplyr::select(Name, stt, geometry) %>%
            mutate(stt = stringr::str_to_sentence(stt)) %>%
            unique(), 
          aes(fill = stt), 
          color = "black", pch = 21,
          size = 6, alpha = 0.7) +
  scale_fill_manual(values = c(plasma(20)[c(19, 10, 2)], "grey35"), name = element_blank()) +
  ##geom_sf(data = state_parks, color = "black", fill = viridis(20)[18], pch = 23,
  ##        size = 3, alpha = 0.7) +
  geom_text_repel(data = unique(sample_points)[, c("lon", "lat", "site")],
                  aes(x = lon, y = lat, label = site),
                  min.segment.length = 0, size = 6,
                  max.overlaps = Inf,
                  nudge_x = c(.05, .07, .06, .05, .08, .07, 0, -.08, -.05, -.05,
                              -.05, -.05, -.02, .02, .02),
                  nudge_y = c(-.02, -.01, .02, -.01, 0, 0, .03, 0, 0, 0, 0, .03,
                              -.04, .03, -.03)) +
  annotation_north_arrow(which_north = "true", location = "tl",
                         width = unit(1.1, "cm"), height = unit(1.1, "cm")) +
  annotation_scale(location = "bl") +
  ggtitle("A. Sampling Locations") +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = c(-114.45, -114.25, -114.05, -113.85)) +
  scale_y_continuous(breaks = c(47.70, 47.82, 47.94, 48.06)) +
  coord_sf(xlim = c(-114.35, -113.85), ylim = c(47.65, 48.15), crs = 4326) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.91),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 14))


# 4. Create map of lakeside developments ----------------------------------

development_map <- ggplot() +
  annotation_map_tile(
    type = 'https://stamen-tiles.a.ssl.fastly.net/terrain-background/${z}/${x}/${y}.png',
    zoom = 11,
    cachedir = "data/map_tiles/") +
  geom_sf(data = developments %>%
            mutate(Development = "Municipality") %>%
            bind_rows(., 
                     state_parks %>%
                       mutate(Development = "State Park")) %>%
          filter(Development != "State Park"), 
          aes(fill = Development), alpha = 0.7, show.legend = "point") +
  geom_sf(data = developments %>%
            mutate(Development = "Municipality") %>%
            bind_rows(., 
                      state_parks %>%
                        mutate(Development = "State Park")) %>%
          filter(Development == "State Park"), 
          aes(fill = Development), shape = 21, alpha = 0.7, size = 5) +
  geom_text_repel(data = data.frame(site = c("Swan R.", "Flathead R."),
                                         lon = c(-114.08086273131164, -114.12845000247533),
                                         lat = c(48.05991806713545, 48.061358255662874)),
                  aes(x = lon, y = lat, label = site),
                  min.segment.length = 0, size = 6, 
                  max.overlaps = Inf,
                  nudge_x = c(.03, -.03),
                  nudge_y = c( .03, -.02)) +
  scale_fill_manual(values = c(viridis(10)[2], viridis(20)[18]), name = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  ggtitle("B. Human Development and Hydrologic Context") +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = c(-114.45, -114.25, -114.05, -113.85)) +
  scale_y_continuous(breaks = c(47.70, 47.82, 47.94, 48.06)) +
  coord_sf(xlim = c(-114.35, -113.85), ylim = c(47.65, 48.15), crs = 4326) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.91),
        legend.key = element_blank(),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 14))

arranged_plots <- ggarrange(plotlist = list(site_map, development_map), ncol = 2)

ggsave(filename = "../figures_tables/combined_maps_revised.png", plot = arranged_plots, 
       width = 14, height = 10, units = "in", bg = "white")


# 5. Create map of littoral PPCP concentrations ---------------------------

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







