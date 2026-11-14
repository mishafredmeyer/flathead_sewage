Last edited:  November 2026

# Food web implications of heterogeneous lakeside human disturbance on nearshore periphyton communities

Authors: Michael F. Meyer, Shawn P. Devlin, Aaron W. E. Galloway, Julie B. Schram, Stephanie E. Hampton

Points of Contact: Michael Meyer (mfmeyer2@wisc.edu)

Repository Type: R scripts supporting publication

## Respository Contents 

This repository contains scripts for the associated manuscript 
"Food web implications of heterogeneous lakeside human disturbance on nearshore periphyton communities".
The goal of the scripts contained within this repository is to support the main analyses 
presented in the manuscript. 

The repository is meant to contain all scripts, data, tables, and figures
relating to the project. The authors structured the project with a few
components that they hope will allow for others to be able to work most
efficiently, but they recognize that this may not be realized. If you
have suggestions for improvements, please feel free to send them over
email to Michael (mfmeeyr2@wisc.edu).

The directories for this analysis need to be structured according to the tree:
```
├───cleaned_data
├───cleaned_disaggregated_data
├───figures_tables
└───R_scripts
```
where data are read from `cleaned_disgaggreated_data` and written to `cleaned_data`. 

The ultimate goal of the analysis is to understand how spatially and temporally heterogeneous 
nutrient addition from seasonal tourism may alter periphyton community composition and, thus, 
the available nutrition to the food web. 

## Repository Organization

- `R_scripts/00_clean_disaggregated_data.R`: This script harmonizes the assortment of datasets into an a flexible format suitable for publication but not an analysis-ready format. 
    - `Inputs`: Pre-aggregated forms of the input datasets described in the manuscript. 
    - `Outputs`: One CSV per observable unit.
- `R_scripts/00a_landsat_temperature_data_aggregation.R`: This script harmonizes in situ temperature data and corrects remote sensed temperature data for future analyses. 
    - `Inputs`: Pre-aggregated forms of the input datasets described in the manuscript as well as aggregated outputs from Google Earth Engine. 
    - `Outputs`: One CSV per observable unit. 
- `R_scripts/01_data_aggregation.R`: This script harmonizes the assortment of datasets into an an analysis-ready format. 
    - `Inputs`: The final outputs of `00_clean_disaggregated_data.R` and `00a_landsat_temperature_data_aggregation.R`.
    - `Outputs`: One CSV per observable unit.   
- `R_scripts/02_analysis.R`: This script visualizes and performs univariate analyses of nutrients and various sewage indicators. 
    - `Inputs`: The final output of `01_data_aggregation.R`.
    - `Outputs`: Several pngs and summary CSVs that report results from models.   
- `R_scripts/03_community_composition_analysis.R`: This script visualizes and performs univariate and multivariate analyses of periphyton community compositions. 
    - `Inputs`: The final output of `01_data_aggregation.R`.
    - `Outputs`: Several pngs and models that evaluate spatio-temporal variation in the periphyton communities.  
- `R_scripts/04_fatty_acid_analysis.R`: This script visualizes and performs univariate and multivariate analyses of periphyton fatty acid compositions. 
    - `Inputs`: The final output of `01_data_aggregation.R`. 
    - `Ouputs`: Several pngs and models that evaluate spatio-temporal variation in the periphyton fatty acids.
- `R_scripts/05_stoichiometry_analysis.R`: This script visualizes and performs univariate analyses of periphyton carbon, nitrogen, and phosphorus stoichiometries. 
    - `Inputs`: The final output of `01_data_aggregation.R`. 
    - `Ouputs`: Several pngs and models that evaluate spatio-temporal variation in the periphyton stoichiometry.
- `R_scripts/06_map.R`: This script generates publication quality maps that present sampling locations and lakeside developments in one figure. 
    - `Inputs`: The final output of `01_data_aggregation.R`. 
    - `Ouputs`: Several maps of the study region.
- `R_scripts/06a_map_post_revisions.R`: This script generates publication quality maps that present sampling locations and lakeside developments. This script differs from `06_map.R` in that these aesthetic differences were made per suggestions from the manuscript's reviewers. 
    - `Inputs`: The final output of `01_data_aggregation.R`. 
    - `Ouputs`: Several maps of the study region. 
- `.gitignore`: A file specifying which files should not be tracked in the git repository. 

## Setup and Package Installation

The scripts contained herein were run using R version 4.4.2.

Necessary Packages to run this codebase include:
- [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html)
- [sf](https://cran.r-project.org/web/packages/sf/index.html)
- [stringi](https://cran.r-project.org/web/packages/stringi/index.html)
- [janitor](https://cran.r-project.org/web/packages/janitor/index.html)
- [car](https://cran.r-project.org/web/packages/car/index.html)
- [OpenStreetMap](https://cran.r-project.org/web/packages/OpenStreetMap/index.html)
- [viridis](https://cran.r-project.org/web/packages/viridis/index.html)
- [readxl](https://cran.r-project.org/web/packages/readxl/index.html)
- [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
- [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
- [MVN](https://cran.r-project.org/web/packages/MVN/index.html)
- [ape](https://cran.r-project.org/web/packages/ape/index.html)
- [ggrepe](https://cran.r-project.org/web/packages/ggrepel/index.html)
- [broom](https://cran.r-project.org/web/packages/broom/index.html)
- [lubridate](https://cran.r-project.org/web/packages/lubridate/index.html)
- [vegan](https://cran.r-project.org/web/packages/vegan/index.html)
- [emmeans](https://cran.r-project.org/web/packages/emmeans/index.html)
- [rstatix](https://cran.r-project.org/web/packages/rstatix/index.html)


## Running the code

This code has been organized in piecemeal scripts. Providing that the end user has downloaded the data, they can be run in the order detailed above. 

## Built With

-   [R](https://www.r-project.org/about.html) - data formatting and processing

