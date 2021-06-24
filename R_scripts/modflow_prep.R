library(raster)
library(sf)
library(tidyverse)
library(data.table)
library(FedData)

hydrolakes_shp <- st_read("E:/GLOBAL_LAKES_AREA/HydroLAKES/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp")

flathead <- hydrolakes_shp %>% 
  filter(grepl(pattern = "flathead", x = Lake_name, ignore.case = TRUE))

## Flathead Hylak_id == 730

plot(flathead)

glcp <- fread("C:/Users/michael.f.meyer/Dropbox/flathead_sewage/cleaned_data/hplc_cleaned_output/Derived-products/glcp.csv",
              integer64 = "character")

glcp_flathead <- glcp %>%
  filter(Hylak_id == 730)

## HYBAS_ID == 7060348000, bsn_lvl = "lev06"

hydrobasins <- st_read("E:/GLOBAL_LAKES_AREA/HydroSat_copied_from_labou/HB_basic/HydroBASINS/hybas_na_lev01-12_v1c/hybas_na_lev06_v1c.shp")

basin_flathead <- hydrobasins %>%
  filter(HYBAS_ID == 7060348000)

ggplot() +
  geom_sf(data = basin_flathead, alpha = 0.4) +
  geom_sf(data = flathead, fill = "blue", alpha = 0.6) +
  geom_sf(data = basins_grid, alpha = 0.2) +
  theme_minimal()

basins_grid <- st_sf(st_make_grid(basin_flathead, what = "polygons", n = c(10,20)))

plot(st_intersection(basins_grid[36,1], basin_flathead))

st_write(obj = basin_flathead, dsn = "C:/Users/michael.f.meyer/Dropbox/flathead_sewage/cleaned_data/flathead_basin.shp",
         driver = "ESRI shapefile")

basin_grid_subset <- st_intersection(basins_grid[36,1], basin_flathead) %>%
  rename("GEOMETRY" = "st_make_grid.basin_flathead..what....polygons...n...c.10..20..")

st_write(obj = basin_grid_subset, dsn = "C:/Users/michael.f.meyer/Dropbox/flathead_sewage/cleaned_data/bains_grid_test.shp",
         driver = "ESRI shapefile")

##### playing with FedData package

test_download <- get_ssurgo(template = st_transform(basin_grid_subset,
                                                    crs = st_crs(basin_grid_subset)), label = "test", 
                            raw.dir = "../../../raw_data/raw_ssurgo_data/")
